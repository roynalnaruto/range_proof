use std::borrow::Cow;

use algebra::{bls12_381::Fr, Bls12_381};
use ff_fft::polynomial::{DensePolynomial as Polynomial};
use poly_commit::kzg10::{Commitment, KZG10, Powers, UniversalParams, VerifierKey};

use crate::errors::Error;

type Kzg10Bls12_381 = KZG10<Bls12_381>;

pub fn trusted_setup<'a>(
    max_degree: usize
) -> Result<(Powers<'a, Bls12_381>, VerifierKey<Bls12_381>), Error> {
    let mut rng = rand::thread_rng();
    let universal_params =
        Kzg10Bls12_381::setup(max_degree, false, &mut rng).unwrap();

    trim(&universal_params, max_degree)
}

pub fn commit(
    powers: &Powers<Bls12_381>,
    p: &Polynomial<Fr>
) -> Commitment<Bls12_381> {
    let (commitment, _) = KZG10::commit(&powers, &p, None, None).unwrap();

    commitment
}

// Copy from https://github.com/scipr-lab/poly-commit/blob/master/src/kzg10/mod.rs#L473
fn trim<'a>(
    pp: &UniversalParams<Bls12_381>,
    mut supported_degree: usize,
) -> Result<(Powers<'a, Bls12_381>, VerifierKey<Bls12_381>), Error> {
    if supported_degree == 1 {
        supported_degree += 1;
    }
    let powers_of_g = pp.powers_of_g[..=supported_degree].to_vec();
    let powers_of_gamma_g = pp.powers_of_gamma_g[..=supported_degree].to_vec();

    let powers = Powers {
        powers_of_g: Cow::Owned(powers_of_g),
        powers_of_gamma_g: Cow::Owned(powers_of_gamma_g),
    };
    let vk = VerifierKey {
        g: pp.powers_of_g[0],
        gamma_g: pp.powers_of_gamma_g[0],
        h: pp.h,
        beta_h: pp.beta_h,
        prepared_h: pp.prepared_h.clone(),
        prepared_beta_h: pp.prepared_beta_h.clone(),
    };
    Ok((powers, vk))
}
