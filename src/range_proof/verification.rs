use algebra::{
    bls12_381::{Fr, g1::{G1Affine, G1Projective}},
    Bls12_381
};
use algebra_core::fields::Field;
use ff_fft::domain::EvaluationDomain;
use num_traits::identities::{One, Zero};
use poly_commit::kzg10::{Commitment, VerifierKey};

use crate::{
    commitment_scheme,
    errors::Error,
    range_proof::{Commitments, Evaluations},
    transcript::TranscriptProtocol,
    utils
};

pub fn verify(
    vk: &VerifierKey<Bls12_381>,
    domain: &EvaluationDomain<Fr>,
    evaluations: &Evaluations,
    commitments: &Commitments,
    aggregate_witness_commitment: Commitment<Bls12_381>,
    shifted_witness_commitment: Commitment<Bls12_381>,
    transcript: &mut dyn TranscriptProtocol
) -> Result<(), Error> {
    // get the evaluation point (rho)
    // and random aggregation challenge for q(X) (tau)
    let tau = transcript.challenge_scalar(b"tau");
    let rho = transcript.challenge_scalar(b"rho");

    // calculate w_cap_commitment
    let w_cap_commitment: Commitment<Bls12_381> =
        compute_w_cap_commitment(&domain, commitments.f, commitments.q, &rho);

    // calculate w2(ρ) and w3(ρ)
    let (w1_part, w2_part, w3_part) = compute_w1_w2_w3_evals(
        &domain, &evaluations.g,
        &evaluations.g_omega,
        &rho, &tau
    );

    // calculate w(ρ)
    // that should zero since w(X) is after all a zero polynomial
    let w_at_rho = w1_part + w2_part + w3_part - evaluations.w_cap;
    if !w_at_rho.is_zero() {
        return Err(Error::ZeroPolynomialCheckFailure);
    }

    // check aggregate witness commitment
    let aggregation_challenge = transcript.challenge_scalar(b"aggregation_challenge");
    let aggregate_poly_commitment =
        commitment_scheme::aggregate_commitments(vec![
            &commitments.g,
            &w_cap_commitment
        ], aggregation_challenge);
    let aggregate_value =
        commitment_scheme::aggregate_values(vec![
            &evaluations.g,
            &evaluations.w_cap
        ], aggregation_challenge);
    if commitment_scheme::check(
        &vk,
        &aggregate_poly_commitment,
        &aggregate_witness_commitment,
        rho,
        aggregate_value
    ) == false {
        return Err(Error::AggregateWitnessCheckFailure)
    }

    // check shifted witness commitment
    let rho_omega = rho * domain.group_gen;
    if commitment_scheme::check(
        &vk,
        &commitments.g,
        &shifted_witness_commitment,
        rho_omega,
        evaluations.g_omega
    ) == false {
        return Err(Error::ShiftedWitnessCheckFailure);
    }

    Ok(())
}

pub fn compute_w1_w2_w3_evals(
    domain: &EvaluationDomain<Fr>,
    g_eval: &Fr,
    g_omega_eval: &Fr,
    rho: &Fr,
    tau: &Fr
) -> (Fr, Fr, Fr) {
    let n_as_ref = utils::as_ref(&Fr::from(domain.size() as u8));
    let one = Fr::one();
    let two = Fr::from(2u8);
    let rho_n_minus_1 = rho.pow(&n_as_ref) - one;
    let w_n_minus_1 = domain.elements().last().unwrap();

    // w1_part
    let w1_eval = (*g_eval) * rho_n_minus_1 / ((*rho) - one);

    // w2
    let w2_eval = (*g_eval) * (one - (*g_eval)) * (rho_n_minus_1) /
        ((*rho) - w_n_minus_1);

    // w3
    let w3_eval = {
        let part_a = (*g_eval) - (two * (*g_omega_eval));
        let part_b = one - part_a;
        let part_c = (*rho) - w_n_minus_1;
        part_a * part_b * part_c
    };

    (w1_eval, (*tau) * w2_eval, (*tau) * (*tau) * w3_eval)
}

pub fn compute_w_cap_commitment(
    domain: &EvaluationDomain<Fr>,
    f_commitment: Commitment<Bls12_381>,
    q_commitment: Commitment<Bls12_381>,
    rho: &Fr
) -> Commitment<Bls12_381> {
    let mut f_commit: G1Projective = f_commitment.0.into();
    let mut q_commit: G1Projective = q_commitment.0.into();
    let (rho_relation_1, rho_relation_2) = compute_rho_relations(&domain, &rho);
    f_commit *= rho_relation_1;
    q_commit *= rho_relation_2;
    let w_cap_commitment = f_commit + q_commit;
    let w_cap_commitment: G1Affine = w_cap_commitment.into();

    Commitment(w_cap_commitment)
}

fn compute_rho_relations(
    domain: &EvaluationDomain<Fr>,
    rho: &Fr
) -> (Fr, Fr) {
    let n_as_ref = utils::as_ref(&Fr::from(domain.size() as u8));
    let one = Fr::one();
    let rho_n_minus_1 = rho.pow(&n_as_ref) - one;
    let rho_n_minus_1_by_rho_minus_1 = rho_n_minus_1 / (rho.clone() - one);

    (rho_n_minus_1_by_rho_minus_1, rho_n_minus_1)
}
