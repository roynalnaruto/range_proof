use algebra::{bls12_381::Fr, Bls12_381};
use ff_fft::EvaluationDomain;
use poly_commit::kzg10::{Commitment, Powers};

pub mod polynomial;
pub mod proof;
pub mod verification;

use crate::{
    errors::Error,
    transcript::TranscriptProtocol
};

pub struct Evaluations {
    pub g: Fr,
    pub g_omega: Fr,
    pub w_cap: Fr,
}

pub struct Commitments {
    pub f: Commitment<Bls12_381>,
    pub g: Commitment<Bls12_381>,
    pub q: Commitment<Bls12_381>,
}

pub struct RangeProof {
    pub evaluations: Evaluations,
    pub commitments: Commitments,
}

impl RangeProof {
    // convenience function to prove: 0 <= z < 2^n
    pub fn prove(
        pk: &Powers<Bls12_381>,
        n: usize,
        z: &Fr,
        transcript: &mut dyn TranscriptProtocol
    ) -> RangeProof {
        let domain: EvaluationDomain<Fr> = EvaluationDomain::<Fr>::new(n).unwrap();

        let (
            g_eval, g_omega_eval, w_cap_eval,
            f_commitment, g_commitment, q_commitment
        ) = proof::prove(&pk, &domain, &z, transcript);

        RangeProof {
            evaluations: Evaluations {
                g: g_eval,
                g_omega: g_omega_eval,
                w_cap: w_cap_eval
            },
            commitments: Commitments {
                f: f_commitment,
                g: g_commitment,
                q: q_commitment
            }
        }
    }

    // convenience function to verify a range proof
    pub fn verify(&self) -> Result<(), Error> {
        // TODO
        verification::verify()
    }
}
