use algebra::bls12_381::Fr;
use bitvec::prelude::*;
use ff_fft::domain::EvaluationDomain;
use ff_fft::polynomial::{DensePolynomial as Polynomial};
use num_traits::identities::{One, Zero};

use crate::utils;

pub fn compute_g(domain: &EvaluationDomain<Fr>, z: &Fr) -> Polynomial<Fr> {
    // get bits for z
    // consider only the first `n` bits
    let mut z_bits = utils::to_bits(&z);
    BitVec::truncate(&mut z_bits, domain.size());

    // push the first evaluation point, i.e. (n-1)th bit of z
    let mut evaluations: Vec<Fr> = Vec::with_capacity(domain.size() - 1);
    let z_n_minus_1 = Fr::from(*z_bits.last().unwrap() as u8);
    evaluations.push(z_n_minus_1);

    // for the rest of bits (n-2 .. 0)
    // g_i = 2* g_(i+1) + z_i
    let mut prev_eval = z_n_minus_1;
    for z_bit in z_bits.iter().rev().skip(1) {
        let eval: Fr = (Fr::from(2u8) * prev_eval) + Fr::from(*z_bit as u8);
        evaluations.push(eval);

        prev_eval = eval;
    }

    // reverse the evaluations
    let evaluations: Vec<Fr> = evaluations.iter().rev().cloned().collect();

    Polynomial::from_coefficients_vec(domain.ifft(&evaluations))
}

pub fn compute_w1_w2(
    domain: &EvaluationDomain<Fr>,
    g: &Polynomial<Fr>,
    z: &Fr
) -> (Polynomial<Fr>, Polynomial<Fr>) {
    let one = Fr::one();
    let w_n_minus_1 = domain.elements().last().unwrap();

    // polynomial: P(x) = x - w^(n-1)
    let x_minus_w_n_minus_1 =
        Polynomial::<Fr>::from_coefficients_vec(vec![-w_n_minus_1, one]);

    // polynomial: P(x) = x^n - 1
    // domain.vanishing_polynomial() outputs a sparse polynomial
    let mut vanishing_coefficients = vec![Fr::zero(); domain.size() + 1];
    vanishing_coefficients[0] = -one;
    vanishing_coefficients[domain.size()] = one;
    let x_n_minus_1: Polynomial::<Fr> =
        Polynomial::<Fr>::from_coefficients_vec(vanishing_coefficients);

    // polynomial: P(x) = x - 1
    let x_minus_1 = Polynomial::<Fr>::from_coefficients_vec(vec![-one, one]);

    // f is the constant polynomial f = z
    let f = Polynomial::<Fr>::from_coefficients_vec(vec![z.clone()]);
    let g_minus_f = g - &f;
    let w1: Polynomial::<Fr> = &(&g_minus_f * &x_n_minus_1) / &x_minus_1;

    // polynomial: P(x) = 1
    let poly_one: Polynomial::<Fr> = Polynomial::<Fr>::from_coefficients_vec(vec![one]);
    let one_minus_g = &poly_one - g;
    let w2: Polynomial::<Fr> = &(&(g * &one_minus_g) * &x_n_minus_1) / &x_minus_w_n_minus_1;

    (w1, w2)
}

pub fn compute_w3(
    domain: &EvaluationDomain<Fr>,
    g: &Polynomial<Fr>
) -> Polynomial<Fr> {
    // w3: [g(X) - 2g(Xw)] * [1 - g(X) + 2g(Xw)] * [X - w^(n-1)]
    // degree of g = n - 1
    // degree of w3 = (n - 1) + (n - 1) + 1 = 2n - 1
    // the new domain can be of size 2n
    let domain_2n = EvaluationDomain::<Fr>::new(2 * domain.size()).unwrap();

    // find evaluations of g in the new domain
    let mut g_evals = domain_2n.fft(&g);

    // since we have doubled the domain size
    // the roots of unity of the new domain will also occur
    // in between the roots of unity of the original domain.
    // hence, if g(X) <- g_evals[i]
    // then g(Xw) <- g_evals[i+2]
    g_evals.push(g_evals[0]);
    g_evals.push(g_evals[1]);

    // calculate evaluations of w3
    let w_n_minus_1 = domain.elements().last().unwrap();
    let two = Fr::from(2u8);
    let w3_evals: Vec<Fr> = (0..domain_2n.size())
        .zip(domain_2n.elements())
        .into_iter()
        .map(|(i, x_i)| {
            let part_a = g_evals[i] - (two * g_evals[i + 2]);
            let part_b = Fr::one() - g_evals[i] + (two * g_evals[i + 2]);
            let part_c = x_i - w_n_minus_1;
            part_a * part_b * part_c
        })
        .collect();

    Polynomial::<Fr>::from_coefficients_vec(domain_2n.ifft(&w3_evals))
}

#[cfg(test)]
mod test {
    use super::*;

    use algebra::{to_bytes, ToBytes};
    use algebra_core::fields::Field;
    use num_traits::identities::One;
    use rand::Rng;

    #[test]
    fn test_compute_g() {
        // n = 8, 2^n = 256, 0 <= z < 2^n
        // degree of polynomial should be (n - 1)
        // it should also evaluate to `z` at x = 1
        let n = 8usize;
        let domain: EvaluationDomain<Fr> = EvaluationDomain::<Fr>::new(n).unwrap();
        let z = Fr::from(100u8);

        let g = compute_g(&domain, &z);
        assert_eq!(g.degree(), n - 1);
        assert_eq!(g.evaluate(Fr::one()), z);

        // n2 = 4, 2^n2 = 16, 0 <= z < 2^n2
        // degree of polynomial should be (n2 - 1)
        // it should also evaluate to `z2` at x = 1
        let n2 = 4usize;
        let domain2: EvaluationDomain<Fr> = EvaluationDomain::<Fr>::new(n2).unwrap();
        let z2 = Fr::from(13u8);

        let g2 = compute_g(&domain2, &z2);
        assert_eq!(g2.degree(), n2 - 1);
        assert_eq!(g2.evaluate(Fr::one()), z2);
    }

    #[test]
    fn test_compute_w1_w2() {
        let n = 8usize;
        let domain: EvaluationDomain<Fr> = EvaluationDomain::<Fr>::new(n).unwrap();
        let z = Fr::from(92u8);
        let g = compute_g(&domain, &z);

        let (w1, w2) = compute_w1_w2(&domain, &g, &z);

        // both w1 and w2 should evaluate to 0 at x = 1
        assert!(w1.evaluate(Fr::one()).is_zero());
        assert!(w2.evaluate(Fr::one()).is_zero());

        // both w1 and w2 should evaluate to 0 at all roots of unity
        for root in domain.elements() {
            assert!(w1.evaluate(root).is_zero());
            assert!(w2.evaluate(root).is_zero());
        }

        // evaluate w1 at a random field element
        let n_as_ref = utils::to_raw_bytes(
            &utils::slice_to_array32(to_bytes![Fr::from(n as u8)].unwrap().as_slice())
        );
        let mut rng = rand::thread_rng();
        let r = Fr::from(rng.gen::<u64>());
        let part_a = g.evaluate(r);
        let part_b = z.clone();
        let part_c = (r.pow(&n_as_ref) - Fr::one()) / (r - Fr::one());
        let w1_expected = (part_a - part_b) * part_c;
        assert_eq!(w1.evaluate(r), w1_expected);

        // evaluate w2 at a random field element
        let r = Fr::from(rng.gen::<u64>());
        let part_a = g.evaluate(r);
        let part_b = Fr::one() - part_a;
        let part_c = {
            let w_n_minus_1 = domain.elements().last().unwrap();
            (r.pow(&n_as_ref) - Fr::one()) / (r - w_n_minus_1)
        };
        let w2_expected = part_a * part_b * part_c;
        assert_eq!(w2.evaluate(r), w2_expected);
    }

    #[test]
    fn test_compute_w3() {
        let n = 8usize;
        let domain: EvaluationDomain<Fr> = EvaluationDomain::<Fr>::new(n).unwrap();
        let z = Fr::from(83u8);
        let g = compute_g(&domain, &z);

        let w3 = compute_w3(&domain, &g);

        // w3 should evaluate to 0 at all roots of unity for original domain
        for root in domain.elements() {
            assert!(w3.evaluate(root).is_zero());
        }

        // w3 degree should be 2n - 1
        assert_eq!(w3.degree(), 2 * domain.size() - 1);

        // evaluate w3 at a random field element
        let w_n_minus_1 = domain.elements().last().unwrap();
        let mut rng = rand::thread_rng();
        let r = Fr::from(rng.gen::<u64>());
        let part_a = g.evaluate(r) - (Fr::from(2u8) * g.evaluate(r * domain.group_gen));
        let part_b = Fr::one() - g.evaluate(r) + (Fr::from(2u8) * g.evaluate(r * domain.group_gen));
        let part_c = r - w_n_minus_1;
        let w3_expected = part_a * part_b * part_c;
        assert_eq!(w3.evaluate(r), w3_expected);

        // evaluate w3 at another random field element
        let mut rng = rand::thread_rng();
        let r = Fr::from(rng.gen::<u64>());
        let part_a = g.evaluate(r) - (Fr::from(2u8) * g.evaluate(r * domain.group_gen));
        let part_b = Fr::one() - g.evaluate(r) + (Fr::from(2u8) * g.evaluate(r * domain.group_gen));
        let part_c = r - w_n_minus_1;
        let w3_expected = part_a * part_b * part_c;
        assert_eq!(w3.evaluate(r), w3_expected);
    }
}
