use bls12_381::Scalar;

use num_bigint::BigUint;
use rand;
use rand_core::{CryptoRng, RngCore};

use crate::errors::Error;
use crate::utils;

/// q - 1 = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000000
const MODULUS_MINUS_1: Scalar = Scalar::from_raw([
    0xffffffff00000000,
    0x53bda402fffe5bfe,
    0x3339d80809a1d805,
    0x73eda753299d7d48,
]);

/// Convenience function to generate random scalar
pub fn random<R: CryptoRng + RngCore>(rng: &mut R) -> Scalar {
    let mut scalar_bytes = [0u8; 64];
    rng.fill_bytes(&mut scalar_bytes);

    Scalar::from_bytes_wide(&scalar_bytes)
}

/// BLS12_381's scalar belongs to a finite field of
/// prime order q
pub fn primitive_nth_root_of_unity(n: usize) -> Result<Scalar, Error> {
    // n should not be zero
    if n == 0usize {
        return Err(Error::ZeroethRootOfUnity)
    }

    // proceed only if n divides (q-1)
    let q_minus_1 = BigUint::from_bytes_le(&MODULUS_MINUS_1.to_bytes());
    let remainder = q_minus_1.clone() % n;
    if remainder.ne(&BigUint::from(0u32)) {
        return Err(Error::TrivialRootOfUnity)
    }

    // get a random scalar
    let mut rng = rand::thread_rng();
    let x = random(&mut rng);

    // calculate the exponent and exponentiate the
    // random scalar
    let exponent = q_minus_1 / n;
    let exponent = exponent.to_bytes_le();
    let exponent = exponent.as_slice();
    let exponent = Scalar::from_bytes(&utils::slice_to_array32(&exponent)).unwrap();
    let w = x.pow(&utils::to_raw_bytes(&exponent));

    // if w exponentiated by n/2 (mod q) is 1
    // it is not the primitive root.
    // we will need to recalculate a random x in this case
    let n_by_2 = Scalar::from((n / 2) as u64);
    let w_pow_n_by_2 = w.pow(&utils::to_raw_bytes(&n_by_2));
    match w_pow_n_by_2.eq(&Scalar::one()) {
        true => primitive_nth_root_of_unity(n),
        _    => Ok(w)
    }
}

/// Calculate the nth roots of unity from the
/// primitive nth root of unity
pub fn nth_roots_of_unity(n: usize, w: &Scalar) -> Vec<Scalar> {
    let mut roots: Vec<Scalar> = Vec::with_capacity(n);
    roots.push(w.clone());

    for i in 2..=n {
        let e = utils::to_raw_bytes(&Scalar::from(i as u64));
        let w_pow_e = w.pow(&e);
        roots.push(w_pow_e);
    }

    roots
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_primitive_nth_root_of_unity() {
        let one = Scalar::one();

        // 6 divides (q - 1)
        // w^6 will be roll back to unity
        // w^7 will be w again
        let sixth_root = primitive_nth_root_of_unity(6usize).unwrap();
        let scalar_6 = utils::to_raw_bytes(&Scalar::from(6u64));
        let scalar_7 = utils::to_raw_bytes(&Scalar::from(7u64));
        let w_pow_6 = sixth_root.pow(&scalar_6);
        let w_pow_7 = sixth_root.pow(&scalar_7);
        assert_eq!(w_pow_6.eq(&one), true);
        assert_eq!(w_pow_7.eq(&one), false);
        assert_eq!(w_pow_7.eq(&sixth_root), true);

        // 32 divides (q - 1)
        // w^32 will be roll back to unity
        let root = primitive_nth_root_of_unity(32usize).unwrap();
        let scalar_32 = utils::to_raw_bytes(&Scalar::from(32u64));
        let scalar_33 = utils::to_raw_bytes(&Scalar::from(33u64));
        let w_pow_32 = root.pow(&scalar_32);
        let w_pow_33 = root.pow(&scalar_33);
        assert_eq!(w_pow_32.eq(&one), true);
        assert_eq!(w_pow_33.eq(&one), false);
        assert_eq!(w_pow_33.eq(&root), true);

        // 5 does not divide (q-1)
        let fifth_root = primitive_nth_root_of_unity(5usize);
        assert!(fifth_root.is_err());

        // zeroeth root of unity does not exist
        let zeroth = primitive_nth_root_of_unity(0usize);
        assert!(zeroth.is_err());
    }

    #[test]
    fn test_nth_roots_of_unity() {
        let one = Scalar::one();

        // n = 16
        let primitive_16 = primitive_nth_root_of_unity(16usize).unwrap();
        let nth_roots_16 = nth_roots_of_unity(16, &primitive_16);
        let sixteen = utils::to_raw_bytes(&Scalar::from(16u64));
        for root in nth_roots_16.iter() {
            let root_pow = root.pow(&sixteen);
            assert_eq!(root_pow.eq(&one), true);
        }

        // n = 6
        let primitive_6 = primitive_nth_root_of_unity(6usize).unwrap();
        let nth_roots_6 = nth_roots_of_unity(6, &primitive_6);
        let six = utils::to_raw_bytes(&Scalar::from(6u64));
        for root in nth_roots_6.iter() {
            let root_pow = root.pow(&six);
            assert_eq!(root_pow.eq(&one), true);
        }

        // n = 64
        let primitive_64 = primitive_nth_root_of_unity(64usize).unwrap();
        let nth_roots_64 = nth_roots_of_unity(64, &primitive_64);
        let sixtyfour = utils::to_raw_bytes(&Scalar::from(64u64));
        for root in nth_roots_64.iter() {
            let root_pow = root.pow(&sixtyfour);
            assert_eq!(root_pow.eq(&one), true);
        }
    }
}
