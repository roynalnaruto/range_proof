use algebra::{bls12_381::Fr, Bls12_381};
use ff_fft::domain::EvaluationDomain;
use poly_commit::kzg10::{Commitment, Powers};
use rand::Rng;

use crate::{
    commitment_scheme::{
        commit,
        create_witness,
        create_aggregate_witness
    },
    range_proof::polynomial,
    transcript::TranscriptProtocol
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
    Commitment<Bls12_381>,
    Commitment<Bls12_381>,
    Commitment<Bls12_381>
) {
    // compute all polynomials
    let mut rng = rand::thread_rng();
    let r = Fr::from(rng.gen::<u64>());
    let f_poly = polynomial::compute_f(&domain, &z, &r);

    let g_poly = polynomial::compute_g(&domain, &z);
    let (w1_poly, w2_poly) = polynomial::compute_w1_w2(&domain, &g_poly, &f_poly);
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

    // evaluate g at `rho * omega`
    let rho_omega = rho * domain.group_gen;
    let g_omega_eval = g_poly.evaluate(rho_omega);

    // compute evaluation of w_cap at ρ
    let w_cap_poly = polynomial::compute_w_cap(&domain, &f_poly, &q_poly, &rho);
    let w_cap_eval = w_cap_poly.evaluate(rho);

    // compute witness for g(X) at ρw
    let shifted_witness_poly = create_witness(&g_poly, &rho_omega);
    let shifted_witness_commitment = commit(&pk, &shifted_witness_poly);

    // compute aggregate witness for
    // g(X) at ρ, f(X) at ρ, w_cap(X) at ρ
    let aggregation_challenge = transcript.challenge_scalar(b"aggregation_challenge");
    let aggregate_witness_poly = create_aggregate_witness(
        vec![g_poly, w_cap_poly],
        &rho, &aggregation_challenge
    );
    let aggregate_witness_commitment = commit(&pk, &aggregate_witness_poly);

    (
        g_eval,
        g_omega_eval,
        w_cap_eval,
        f_commitment,
        g_commitment,
        q_commitment,
        aggregate_witness_commitment,
        shifted_witness_commitment
    )
}
