use algebra::{bls12_381::Fr, Bls12_381};
use poly_commit::kzg10::Commitment;

pub mod polynomial;

struct Evaluations {
    f: Fr,
    g: Fr,
    g_omega: Fr,
    w_cap: Fr,
}

struct Commitments {
    f: Commitment<Bls12_381>,
    g: Commitment<Bls12_381>,
    q: Commitment<Bls12_381>,
}

struct RangeProof {
    evaluations: Evaluations,
    commitments: Commitments,
}
