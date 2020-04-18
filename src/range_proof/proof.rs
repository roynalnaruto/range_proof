use algebra::{bls12_381::Fr, Bls12_381};
use algebra_core::fields::Field;
use ff_fft::domain::EvaluationDomain;
use num_traits::identities::One;
use poly_commit::kzg10::{Commitment, Powers};

use crate::{
    commitment_scheme::{commit},
    range_proof::polynomial,
    transcript::TranscriptProtocol,
    utils
};

pub fn prove(
    pk: &Powers<Bls12_381>,
    domain: &EvaluationDomain<Fr>,
    z: &Fr,
    transcript: &mut dyn TranscriptProtocol
) -> (
    Fr, Fr, Fr,
    Commitment<Bls12_381>,
    Commitment<Bls12_381>,
    Commitment<Bls12_381>
) {
    // compute all polynomials
    let f_poly = polynomial::compute_f(&z);
    let g_poly = polynomial::compute_g(&domain, &z);
    let (w1_poly, w2_poly) = polynomial::compute_w1_w2(&domain, &g_poly, &z);
    let w3_poly = polynomial::compute_w3(&domain, &g_poly);

    // aggregate w1, w2 and w3 to compute quotient polynomial
    // `tau` is the random scalar for aggregation
    let tau = transcript.challenge_scalar(b"tau");
    let (q_poly, _) =
        polynomial::compute_q(&domain, &w1_poly, &w2_poly, &w3_poly, &tau);

    // compute commitments to polynomials
    let f_commitment = commit(&pk, &f_poly);
    let g_commitment = commit(&pk, &g_poly);
    let q_commitment = commit(&pk, &q_poly);

    // `rho` is the random evaluation point
    let rho = transcript.challenge_scalar(b"rho");
    let g_eval = g_poly.evaluate(rho);
    let g_omega_eval = g_poly.evaluate(rho * domain.group_gen);

    // compute evaluation of w_cap at `rho`
    let n_as_ref = utils::as_ref(&Fr::from(domain.size() as u8));
    let one = Fr::one();
    let rho_n_minus_1 = rho.pow(&n_as_ref) - Fr::one();
    let part_a = z.clone() * rho_n_minus_1 / (rho - one);
    let part_b = q_poly.evaluate(rho) * rho_n_minus_1;
    let w_cap_eval = part_a + part_b;

    (
        g_eval,
        g_omega_eval,
        w_cap_eval,
        f_commitment,
        g_commitment,
        q_commitment
    )
}
