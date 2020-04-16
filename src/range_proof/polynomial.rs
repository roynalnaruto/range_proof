use algebra::bls12_381::Fr;
use bitvec::prelude::*;
use ff_fft::domain::EvaluationDomain;
use ff_fft::polynomial::{DensePolynomial as Polynomial};

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

#[cfg(test)]
mod test {
    use super::*;

    use num_traits::identities::One;

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
}
