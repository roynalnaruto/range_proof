use crate::utils;
use ark_bls12_381::Fr;
use ark_ff::Field;
use ark_poly::EvaluationDomain;
use ark_poly::{
    univariate::DensePolynomial as Polynomial, DenseUVPolynomial, GeneralEvaluationDomain,
};
use bitvec::prelude::*;
use num_traits::identities::{One, Zero};

pub fn compute_f(domain: &GeneralEvaluationDomain<Fr>, z: &Fr, r: &Fr) -> Polynomial<Fr> {
    // f is a linear polynomial: f(1) = z
    Polynomial::<Fr>::from_coefficients_vec(domain.ifft(&vec![z.clone(), r.clone()]))
}

pub fn compute_g(
    domain: &GeneralEvaluationDomain<Fr>,
    z: &Fr,
    alpha: &Fr,
    beta: &Fr,
) -> Polynomial<Fr> {
    // get bits for z
    // consider only the first `n` bits
    let mut z_bits = utils::to_bits(&z);
    BitVec::truncate(&mut z_bits, domain.size());

    // push the first evaluation point, i.e. (n-1)th bit of z
    let mut evaluations: Vec<Fr> = Vec::with_capacity(domain.size());
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

    // compute g
    let g_poly = Polynomial::from_coefficients_vec(domain.ifft(&evaluations));

    // extended domain
    let domain_2n: GeneralEvaluationDomain<Fr> =
        GeneralEvaluationDomain::<Fr>::new(domain.size() + 1).unwrap();

    // map the original g_poly to domain(n+1)
    // add random values alpha and beta as evaluations of g
    // at all even indices, g_evals[2k] matches
    // the evaluation at some original root of unity
    // Hence only update two odd indices with alpha and beta
    // this makes g evaluate to the expected evaluations at all
    // roots of unity of domain size `n`, but makes is a different polynomial
    let mut g_evals = domain_2n.fft(&g_poly);
    g_evals[1] = *alpha;
    g_evals[3] = *beta;

    Polynomial::from_coefficients_vec(domain_2n.ifft(&g_evals))
}

pub fn compute_w1_w2(
    domain: &GeneralEvaluationDomain<Fr>,
    g: &Polynomial<Fr>,
    f: &Polynomial<Fr>,
) -> (Polynomial<Fr>, Polynomial<Fr>) {
    let one = Fr::one();
    let w_n_minus_1 = domain.elements().last().unwrap();

    // polynomial: P(x) = x - w^(n-1)
    let x_minus_w_n_minus_1 = Polynomial::<Fr>::from_coefficients_vec(vec![-w_n_minus_1, one]);

    // polynomial: P(x) = x^n - 1
    // domain.vanishing_polynomial() outputs a sparse polynomial
    let mut vanishing_coefficients = vec![Fr::zero(); domain.size() + 1];
    vanishing_coefficients[0] = -one;
    vanishing_coefficients[domain.size()] = one;
    let x_n_minus_1: Polynomial<Fr> =
        Polynomial::<Fr>::from_coefficients_vec(vanishing_coefficients);

    // polynomial: P(x) = x - 1
    let x_minus_1 = Polynomial::<Fr>::from_coefficients_vec(vec![-one, one]);

    let g_minus_f = g - f;
    let w1: Polynomial<Fr> = &(&g_minus_f * &x_n_minus_1) / &x_minus_1;

    // polynomial: P(x) = 1
    let poly_one: Polynomial<Fr> = Polynomial::<Fr>::from_coefficients_vec(vec![one]);
    let one_minus_g = &poly_one - g;
    let w2: Polynomial<Fr> = &(&(g * &one_minus_g) * &x_n_minus_1) / &x_minus_w_n_minus_1;

    (w1, w2)
}

pub fn compute_w3(
    domain: &GeneralEvaluationDomain<Fr>,
    domain_2n: &GeneralEvaluationDomain<Fr>,
    g: &Polynomial<Fr>,
) -> Polynomial<Fr> {
    // w3: [g(X) - 2g(Xw)] * [1 - g(X) + 2g(Xw)] * [X - w^(n-1)]
    // degree of g = n - 1
    // degree of w3 = (2n - 1) + (2n - 1) + 1 = 4n - 1
    // the new domain can be of size 4n
    let domain_4n = GeneralEvaluationDomain::<Fr>::new(2 * domain_2n.size()).unwrap();

    // find evaluations of g in the new domain
    let mut g_evals = domain_4n.fft(&g);

    // since we have doubled the domain size
    // the roots of unity of the new domain will also occur
    // in between the roots of unity of the original domain.
    // hence, if g(X) <- g_evals[i]
    // then g(Xw) <- g_evals[i+4]
    g_evals.push(g_evals[0]);
    g_evals.push(g_evals[1]);
    g_evals.push(g_evals[2]);
    g_evals.push(g_evals[3]);

    // calculate evaluations of w3
    let w_n_minus_1 = domain.elements().last().unwrap();
    let two = Fr::from(2u8);
    let w3_evals: Vec<Fr> = (0..domain_4n.size())
        .zip(domain_4n.elements())
        .into_iter()
        .map(|(i, x_i)| {
            let part_a = g_evals[i] - (two * g_evals[i + 4]);
            let part_b = Fr::one() - g_evals[i] + (two * g_evals[i + 4]);
            let part_c = x_i - w_n_minus_1;
            part_a * part_b * part_c
        })
        .collect();

    Polynomial::<Fr>::from_coefficients_vec(domain_4n.ifft(&w3_evals))
}

pub fn compute_q(
    domain: &GeneralEvaluationDomain<Fr>,
    w1: &Polynomial<Fr>,
    w2: &Polynomial<Fr>,
    w3: &Polynomial<Fr>,
    tau: &Fr,
) -> (Polynomial<Fr>, Polynomial<Fr>) {
    // find constant polynomials for tau and tau^2
    let poly_tau: Polynomial<Fr> = Polynomial::<Fr>::from_coefficients_vec(vec![tau.clone()]);
    let poly_tau_2: Polynomial<Fr> = Polynomial::<Fr>::from_coefficients_vec(vec![tau.square()]);

    // find linear combination of w1, w2, w3
    let lc = &(w1 + &(w2 * &poly_tau)) + &(w3 * &poly_tau_2);

    lc.divide_by_vanishing_poly(*domain).unwrap()
}

pub fn compute_w_cap(
    domain: &GeneralEvaluationDomain<Fr>,
    f: &Polynomial<Fr>,
    q: &Polynomial<Fr>,
    rho: &Fr,
) -> Polynomial<Fr> {
    let n_as_ref = utils::as_ref(&Fr::from(domain.size() as u8));
    let one = Fr::one();
    let rho_n_minus_1 = rho.pow(&n_as_ref) - one;
    let rho_n_minus_1_by_rho_minus_1 = rho_n_minus_1 / ((*rho) - one);

    let rho_poly_1: Polynomial<Fr> =
        Polynomial::<Fr>::from_coefficients_vec(vec![rho_n_minus_1_by_rho_minus_1]);
    let rho_poly_2: Polynomial<Fr> = Polynomial::<Fr>::from_coefficients_vec(vec![rho_n_minus_1]);

    &(f * &rho_poly_1) + &(q * &rho_poly_2)
}

pub fn compute_w(
    domain: &GeneralEvaluationDomain<Fr>,
    w1: &Polynomial<Fr>,
    w2: &Polynomial<Fr>,
    w3: &Polynomial<Fr>,
    q: &Polynomial<Fr>,
    tau: &Fr,
) -> Polynomial<Fr> {
    let poly_tau: Polynomial<Fr> = Polynomial::<Fr>::from_coefficients_vec(vec![*tau]);
    let poly_tau_2: Polynomial<Fr> = Polynomial::<Fr>::from_coefficients_vec(vec![tau.square()]);

    &(&(w1 + &(w2 * &poly_tau)) + &(w3 * &poly_tau_2)) - &q.mul_by_vanishing_poly(*domain)
}

pub fn compute_w2_w3_parts(
    w2: &Polynomial<Fr>,
    w3: &Polynomial<Fr>,
    tau: &Fr,
) -> (Polynomial<Fr>, Polynomial<Fr>) {
    let poly_tau = Polynomial::<Fr>::from_coefficients_vec(vec![tau.clone()]);
    let poly_tau_2 = Polynomial::<Fr>::from_coefficients_vec(vec![tau.square()]);

    (w2 * &poly_tau, w3 * &poly_tau_2)
}

pub fn compute_w1_part(domain: &GeneralEvaluationDomain<Fr>, g: &Polynomial<Fr>) -> Polynomial<Fr> {
    let one = Fr::one();
    let divisor = Polynomial::<Fr>::from_coefficients_vec(vec![-one, one]);

    &g.mul_by_vanishing_poly(*domain) / &divisor
}

#[cfg(test)]
mod test {
    use super::*;
    use ark_serialize::CanonicalSerialize;

    use crate::{
        commitment_scheme::{commit, trusted_setup},
        range_proof::verification::compute_w_cap_commitment,
    };
    use ark_bls12_381::Fr;
    use ark_ff::Field;
    use ark_poly::Polynomial;
    use ark_serialize::Compress;
    use num_traits::identities::One;
    use rand::Rng;

    #[test]
    fn test_compute_f() {
        let n = 8usize;
        let domain: GeneralEvaluationDomain<Fr> = GeneralEvaluationDomain::<Fr>::new(n).unwrap();
        let z = Fr::from(2u8);
        let r = Fr::from(4u8);
        let f = compute_f(&domain, &z, &r);

        let mut rng = rand::thread_rng();
        let rho = Fr::from(rng.gen::<u64>());

        assert_eq!(f.evaluate(&Fr::one()), z);
        assert_eq!(f.evaluate(&domain.group_gen()), r);
        assert!(f.evaluate(&rho).ne(&z));
        assert!(f.evaluate(&rho).ne(&r));
    }

    #[test]
    fn test_compute_g() {
        // n = 8, 2^n = 256, 0 <= z < 2^n
        // degree of polynomial should be (n - 1)
        // it should also evaluate to `z` at x = 1
        let n = 8usize;
        let domain: GeneralEvaluationDomain<Fr> = GeneralEvaluationDomain::<Fr>::new(n).unwrap();
        let z = Fr::from(100u8);

        let mut rng = rand::thread_rng();
        let alpha = Fr::from(rng.gen::<u64>());
        let beta = Fr::from(rng.gen::<u64>());
        let g = compute_g(&domain, &z, &alpha, &beta);
        assert_eq!(g.degree(), 2usize * n - 1);
        assert_eq!(g.evaluate(&Fr::one()), z);

        // n2 = 4, 2^n2 = 16, 0 <= z < 2^n2
        // degree of polynomial should be (n2 - 1)
        // it should also evaluate to `z2` at x = 1
        let n2 = 4usize;
        let domain2: GeneralEvaluationDomain<Fr> = GeneralEvaluationDomain::<Fr>::new(n2).unwrap();
        let z2 = Fr::from(13u8);

        let alpha = Fr::from(rng.gen::<u64>());
        let beta = Fr::from(rng.gen::<u64>());
        let g2 = compute_g(&domain2, &z2, &alpha, &beta);
        assert_eq!(g2.degree(), 2usize * n2 - 1);
        assert_eq!(g2.evaluate(&Fr::one()), z2);
    }

    #[test]
    fn test_compute_w1_w2() {
        let n = 8usize;
        let domain: GeneralEvaluationDomain<Fr> = GeneralEvaluationDomain::<Fr>::new(n).unwrap();

        let mut rng = rand::thread_rng();
        let r = Fr::from(rng.gen::<u64>());
        let alpha = Fr::from(rng.gen::<u64>());
        let beta = Fr::from(rng.gen::<u64>());
        let z = Fr::from(92u8);
        let f = compute_f(&domain, &z, &r);
        let g = compute_g(&domain, &z, &alpha, &beta);

        let (w1, w2) = compute_w1_w2(&domain, &g, &f);

        // both w1 and w2 should evaluate to 0 at x = 1
        assert!(w1.evaluate(&Fr::one()).is_zero());
        assert!(w2.evaluate(&Fr::one()).is_zero());

        // both w1 and w2 should evaluate to 0 at all roots of unity
        for root in domain.elements() {
            assert!(w1.evaluate(&root).is_zero());
            assert!(w2.evaluate(&root).is_zero());
        }

        let mut bytes = Vec::new();
        let z = Fr::from(n as u8);
        z.serialize_with_mode(&mut bytes, Compress::No).unwrap();
        // evaluate w1 at a random field element
        let n_as_ref = utils::to_raw_bytes(&utils::slice_to_array32(bytes.as_slice()));
        let r = Fr::from(rng.gen::<u64>());
        let part_a = g.evaluate(&r);
        let part_b = f.evaluate(&r);
        let part_c = (r.pow(&n_as_ref) - Fr::one()) / (r - Fr::one());
        let w1_expected = (part_a - part_b) * part_c;
        assert_eq!(w1.evaluate(&r), w1_expected);

        // evaluate w2 at a random field element
        let r = Fr::from(rng.gen::<u64>());
        let part_a = g.evaluate(&r);
        let part_b = Fr::one() - part_a;
        let part_c = {
            let w_n_minus_1 = domain.elements().last().unwrap();
            (r.pow(&n_as_ref) - Fr::one()) / (r - w_n_minus_1)
        };
        let w2_expected = part_a * part_b * part_c;
        assert_eq!(w2.evaluate(&r), w2_expected);
    }

    #[test]
    fn test_compute_w3() {
        let n = 8usize;
        let domain: GeneralEvaluationDomain<Fr> = GeneralEvaluationDomain::<Fr>::new(n).unwrap();
        let domain_2n: GeneralEvaluationDomain<Fr> =
            GeneralEvaluationDomain::<Fr>::new(2usize * n).unwrap();

        let z = Fr::from(83u8);
        let mut rng = rand::thread_rng();
        let alpha = Fr::from(rng.gen::<u64>());
        let beta = Fr::from(rng.gen::<u64>());
        let g = compute_g(&domain, &z, &alpha, &beta);

        let w3 = compute_w3(&domain, &domain_2n, &g);

        // w3 should evaluate to 0 at all roots of unity for original domain
        for root in domain.elements() {
            assert!(w3.evaluate(&root).is_zero());
        }

        // w3 degree should be 4n - 1
        assert_eq!(w3.degree(), 4 * domain.size() - 1);

        // evaluate w3 at a random field element
        let w_n_minus_1 = domain.elements().last().unwrap();
        let r = Fr::from(rng.gen::<u64>());
        let part_a = g.evaluate(&r) - (Fr::from(2u8) * g.evaluate(&(r * domain.group_gen())));
        let part_b =
            Fr::one() - g.evaluate(&r) + (Fr::from(2u8) * g.evaluate(&(r * domain.group_gen())));
        let part_c = r - w_n_minus_1;
        let w3_expected = part_a * part_b * part_c;
        assert_eq!(w3.evaluate(&r), w3_expected);

        // evaluate w3 at another random field element
        let mut rng = rand::thread_rng();
        let r = Fr::from(rng.gen::<u64>());
        let part_a = g.evaluate(&r) - (Fr::from(2u8) * g.evaluate(&(r * domain.group_gen())));
        let part_b =
            Fr::one() - g.evaluate(&r) + (Fr::from(2u8) * g.evaluate(&(r * domain.group_gen())));
        let part_c = r - w_n_minus_1;
        let w3_expected = part_a * part_b * part_c;
        assert_eq!(w3.evaluate(&r), w3_expected);
    }

    #[test]
    fn test_compute_q() {
        let n = 8usize;
        let domain: GeneralEvaluationDomain<Fr> = GeneralEvaluationDomain::<Fr>::new(n).unwrap();
        let domain_2n: GeneralEvaluationDomain<Fr> =
            GeneralEvaluationDomain::<Fr>::new(2usize * n).unwrap();
        let mut rng = rand::thread_rng();

        let z = Fr::from(68u8);
        let r = Fr::from(rng.gen::<u64>());
        let alpha = Fr::from(rng.gen::<u64>());
        let beta = Fr::from(rng.gen::<u64>());
        let f = compute_f(&domain, &z, &r);
        let g = compute_g(&domain, &z, &alpha, &beta);
        let (w1, w2) = compute_w1_w2(&domain, &g, &f);
        let w3 = compute_w3(&domain, &domain_2n, &g);

        let mut rng = rand::thread_rng();
        let tau = Fr::from(rng.gen::<u64>());

        let (_q, q_rem) = compute_q(&domain, &w1, &w2, &w3, &tau);

        // since the linear combination should also
        // satisfy all roots of unity,
        // q_rem should be a zero polynomial
        assert!(q_rem.is_zero());
    }

    #[test]
    fn test_compute_w_cap() {
        // initial setup
        let n = 8usize;
        let (pk, _) = trusted_setup(4usize * n).unwrap();
        let domain: GeneralEvaluationDomain<Fr> = GeneralEvaluationDomain::<Fr>::new(n).unwrap();
        let domain_2n: GeneralEvaluationDomain<Fr> =
            GeneralEvaluationDomain::<Fr>::new(2usize * n).unwrap();

        // random numbers
        let mut rng = rand::thread_rng();
        let r = Fr::from(rng.gen::<u64>());
        let tau = Fr::from(rng.gen::<u64>());
        let rho = Fr::from(rng.gen::<u64>());
        let alpha = Fr::from(rng.gen::<u64>());
        let beta = Fr::from(rng.gen::<u64>());

        // compute polynomials
        let z = Fr::from(68u8);
        let f = compute_f(&domain, &z, &r);
        let g = compute_g(&domain, &z, &alpha, &beta);
        let (w1, w2) = compute_w1_w2(&domain, &g, &f);
        let w3 = compute_w3(&domain, &domain_2n, &g);
        let (q, _) = compute_q(&domain, &w1, &w2, &w3, &tau);
        let w_cap = compute_w_cap(&domain, &f, &q, &rho);

        // compute commitments
        let f_commitment = commit(&pk, &f);
        let q_commitment = commit(&pk, &q);
        let w_cap_commitment_expected = commit(&pk, &w_cap);

        // calculate w_cap commitment
        // fact that commitment scheme is additively homomorphic
        let w_cap_commitment_calculated =
            compute_w_cap_commitment(&domain, f_commitment, q_commitment, &rho);

        assert_eq!(w_cap_commitment_expected, w_cap_commitment_calculated);
    }

    #[test]
    fn test_compute_w1_part() {
        let n = 8usize;
        let domain: GeneralEvaluationDomain<Fr> = GeneralEvaluationDomain::<Fr>::new(n).unwrap();
        let z = Fr::from(92u8);
        let mut rng = rand::thread_rng();
        let alpha = Fr::from(rng.gen::<u64>());
        let beta = Fr::from(rng.gen::<u64>());
        let g = compute_g(&domain, &z, &alpha, &beta);

        let rho = Fr::from(rng.gen::<u64>());
        let g_eval = g.evaluate(&rho);

        let n_as_ref = utils::as_ref(&Fr::from(domain.size() as u8));
        let one = Fr::one();
        let rho_n_minus_1 = rho.pow(&n_as_ref) - one;

        let w1_part_poly = compute_w1_part(&domain, &g);

        assert_eq!(
            w1_part_poly.evaluate(&rho),
            g_eval * rho_n_minus_1 / (rho - one)
        )
    }

    #[test]
    fn test_compute_w2_w3_part() {
        let n = 8usize;
        let domain: GeneralEvaluationDomain<Fr> = GeneralEvaluationDomain::<Fr>::new(n).unwrap();
        let domain_2n: GeneralEvaluationDomain<Fr> =
            GeneralEvaluationDomain::<Fr>::new(2usize * n).unwrap();
        let mut rng = rand::thread_rng();

        let z = Fr::from(92u8);
        let r = Fr::from(rng.gen::<u64>());
        let alpha = Fr::from(rng.gen::<u64>());
        let beta = Fr::from(rng.gen::<u64>());
        let f = compute_f(&domain, &z, &r);
        let g = compute_g(&domain, &z, &alpha, &beta);
        let (_, w2) = compute_w1_w2(&domain, &g, &f);
        let w3 = compute_w3(&domain, &domain_2n, &g);

        let tau = Fr::from(rng.gen::<u64>());
        let rho = Fr::from(rng.gen::<u64>());
        let g_eval = g.evaluate(&rho);
        let g_omega_eval = g.evaluate(&(rho * domain.group_gen()));

        let n_as_ref = utils::as_ref(&Fr::from(domain.size() as u8));
        let one = Fr::one();
        let two = Fr::from(2u8);
        let rho_n_minus_1 = rho.pow(&n_as_ref) - one;
        let w_n_minus_1 = domain.elements().last().unwrap();

        let (w2_part_poly, w3_part_poly) = compute_w2_w3_parts(&w2, &w3, &tau);

        assert_eq!(
            w2_part_poly.evaluate(&rho),
            tau * g_eval * (one - g_eval) * (rho_n_minus_1) / (rho - w_n_minus_1)
        );

        assert_eq!(w3_part_poly.evaluate(&rho), {
            let part_a = g_eval - (two * g_omega_eval);
            let part_b = one - part_a;
            let part_c = rho - w_n_minus_1;
            tau * tau * part_a * part_b * part_c
        });
    }
}
