use ark_bls12_381::{Bls12_381, Fr};
use ark_poly::EvaluationDomain;
use ark_poly::GeneralEvaluationDomain;
use ark_poly_commit::kzg10::{Commitment, Powers, VerifierKey};
pub mod polynomial;
pub mod proof;
pub mod verification;

use crate::{errors::Error, transcript::TranscriptProtocol};

#[derive(Debug)]
pub struct Evaluations {
    pub g: Fr,
    pub g_omega: Fr,
    pub w_cap: Fr,
}

#[derive(Debug)]
pub struct Commitments {
    pub f: Commitment<Bls12_381>,
    pub g: Commitment<Bls12_381>,
    pub q: Commitment<Bls12_381>,
}

#[derive(Debug)]
pub struct RangeProof {
    pub evaluations: Evaluations,
    pub commitments: Commitments,
    pub aggregate_witness_commitment: Commitment<Bls12_381>,
    pub shifted_witness_commitment: Commitment<Bls12_381>,
}

impl RangeProof {
    // convenience function to prove: 0 <= z < 2^n
    pub fn prove(
        pk: &Powers<Bls12_381>,
        n: usize,
        z: &Fr,
        transcript: &mut dyn TranscriptProtocol,
    ) -> RangeProof {
        let domain: GeneralEvaluationDomain<Fr> = GeneralEvaluationDomain::<Fr>::new(n).unwrap();

        let (
            g_eval,
            g_omega_eval,
            w_cap_eval,
            f_commitment,
            g_commitment,
            q_commitment,
            aggregate_witness_commitment,
            shifted_witness_commitment,
        ) = proof::prove(&pk, &domain, &z, transcript);

        RangeProof {
            evaluations: Evaluations {
                g: g_eval,
                g_omega: g_omega_eval,
                w_cap: w_cap_eval,
            },
            commitments: Commitments {
                f: f_commitment,
                g: g_commitment,
                q: q_commitment,
            },
            aggregate_witness_commitment: aggregate_witness_commitment,
            shifted_witness_commitment: shifted_witness_commitment,
        }
    }

    // convenience function to verify a range proof
    pub fn verify(
        &self,
        vk: &VerifierKey<Bls12_381>,
        n: usize,
        transcript: &mut dyn TranscriptProtocol,
    ) -> Result<(), Error> {
        let domain: GeneralEvaluationDomain<Fr> = GeneralEvaluationDomain::<Fr>::new(n).unwrap();

        verification::verify(
            &vk,
            &domain,
            &self.evaluations,
            &self.commitments,
            self.aggregate_witness_commitment,
            self.shifted_witness_commitment,
            transcript,
        )
    }
}

#[cfg(test)]
mod test {
    use super::*;

    use merlin::Transcript;

    use crate::commitment_scheme::trusted_setup;

    #[test]
    fn test_prove() {
        let n = 8usize;
        let (pk, _vk) = trusted_setup(4usize * n).unwrap();

        let mut transcript = Transcript::new(b"range_proof");

        // prove that 0 <= z < 2^n
        let z = Fr::from(100u8);

        // get proof
        let _range_proof = RangeProof::prove(&pk, n, &z, &mut transcript);
    }

    #[test]
    fn test_true_verify() {
        let n = 8usize;
        let (pk, vk) = trusted_setup(4usize * n).unwrap();

        let mut proof_transcript = Transcript::new(b"range_proof");
        let mut verification_transcript = Transcript::new(b"range_proof");

        // prove that 0 <= z < 2^n
        let z = Fr::from(100u8);

        // get proof
        let range_proof = RangeProof::prove(&pk, n, &z, &mut proof_transcript);

        let result = RangeProof::verify(&range_proof, &vk, n, &mut verification_transcript);

        assert!(result.is_ok());
    }

    #[test]
    fn test_false_verify() {
        let n = 8usize;
        let (pk, vk) = trusted_setup(4usize * n).unwrap();

        let mut proof_transcript = Transcript::new(b"range_proof");
        let mut verification_transcript = Transcript::new(b"range_proof");

        // prove that 0 <= z < 2^n
        // z is not in the range
        let z = Fr::from(500u32);

        // get proof
        let range_proof = RangeProof::prove(&pk, n, &z, &mut proof_transcript);

        let result = RangeProof::verify(&range_proof, &vk, n, &mut verification_transcript);

        assert!(result.is_err());
    }
}
