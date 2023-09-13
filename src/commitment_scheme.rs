use std::{borrow::Cow, ops::Mul};

use ark_bls12_381::{g1::G1Projective, Bls12_381, Fr};
use ark_poly::{univariate::DensePolynomial as Polynomial, DenseUVPolynomial};
use ark_poly_commit::kzg10::{Commitment, Powers, Proof, UniversalParams, VerifierKey, KZG10};
use ark_std::test_rng;
use num_traits::identities::{One, Zero};

use crate::errors::Error;

type Kzg10Bls12_381 = KZG10<Bls12_381, Polynomial<Fr>>;

pub fn trusted_setup<'a>(
    max_degree: usize,
) -> Result<(Powers<'a, Bls12_381>, VerifierKey<Bls12_381>), Error> {
    let universal_params = Kzg10Bls12_381::setup(max_degree, false, &mut test_rng()).unwrap();

    trim(&universal_params, max_degree)
}

pub fn commit(powers: &Powers<Bls12_381>, p: &Polynomial<Fr>) -> Commitment<Bls12_381> {
    let (commitment, _) = KZG10::commit(&powers, p, None, None).unwrap();

    commitment
}

pub fn create_witness(polynomial: &Polynomial<Fr>, point: &Fr) -> Polynomial<Fr> {
    let divisor: Polynomial<Fr> =
        Polynomial::<Fr>::from_coefficients_vec(vec![-point.clone(), Fr::one()]);

    polynomial / &divisor
}

pub fn create_aggregate_witness(
    polynomials: Vec<Polynomial<Fr>>,
    point: &Fr,
    challenge: &Fr,
) -> Polynomial<Fr> {
    let mut power = Fr::one();
    let mut result = Polynomial::<Fr>::zero();

    for polynomial in polynomials {
        let tmp_polynomial = &polynomial * &Polynomial::<Fr>::from_coefficients_vec(vec![power]);
        result += &tmp_polynomial;
        power *= challenge;
    }

    let divisor: Polynomial<Fr> =
        Polynomial::<Fr>::from_coefficients_vec(vec![-point.clone(), Fr::one()]);

    &result / &divisor
}

pub fn aggregate_commitments(
    commitments: Vec<&Commitment<Bls12_381>>,
    aggregation_challenge: Fr,
) -> Commitment<Bls12_381> {
    let mut powers = Fr::one();
    let mut result = G1Projective::zero();

    for commitment in commitments {
        let intermediate_comm = commitment.0.mul(powers);
        result += &intermediate_comm;
        powers = powers * aggregation_challenge;
    }

    Commitment(result.into())
}

pub fn aggregate_values(values: Vec<&Fr>, aggregation_challenge: Fr) -> Fr {
    let mut powers = Fr::one();
    let mut result = Fr::zero();

    for value in values {
        let intermediate_value = *value * powers;
        result += &intermediate_value;
        powers = powers * aggregation_challenge;
    }

    result
}

pub fn check(
    vk: &VerifierKey<Bls12_381>,
    poly_commitment: &Commitment<Bls12_381>,
    witness_commitment: &Commitment<Bls12_381>,
    point: Fr,
    value: Fr,
) -> bool {
    let proof = Proof {
        w: witness_commitment.0,
        random_v: Some(Fr::from(0u8)),
    };

    Kzg10Bls12_381::check(vk, poly_commitment, point, value, &proof).unwrap()
}

pub fn batch_check(
    vk: &VerifierKey<Bls12_381>,
    poly_commitments: Vec<Commitment<Bls12_381>>,
    witness_commitments: Vec<Commitment<Bls12_381>>,
    points: Vec<Fr>,
    values: Vec<Fr>,
) -> bool {
    let mut proofs: Vec<Proof<Bls12_381>> = Vec::new();
    for witness in witness_commitments {
        let proof = Proof {
            w: witness.0,
            random_v: Some(Fr::zero()),
        };
        proofs.push(proof);
    }

    Kzg10Bls12_381::batch_check(
        vk,
        poly_commitments.as_slice(),
        &points,
        &values,
        &proofs,
        &mut test_rng(),
    )
    .unwrap()
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
    let powers_of_gamma_g = (0..=supported_degree)
        .map(|i| pp.powers_of_gamma_g[&i])
        .collect();

    let powers = Powers {
        powers_of_g: Cow::Owned(powers_of_g),
        powers_of_gamma_g: Cow::Owned(powers_of_gamma_g),
    };
    let vk = VerifierKey {
        g: pp.powers_of_g[0],
        gamma_g: pp.powers_of_gamma_g[&0],
        h: pp.h,
        beta_h: pp.beta_h,
        prepared_h: pp.prepared_h.clone(),
        prepared_beta_h: pp.prepared_beta_h.clone(),
    };
    Ok((powers, vk))
}
