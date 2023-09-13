use ark_bls12_381::{Bls12_381, Fr};
use ark_poly_commit::kzg10::Commitment;
use ark_serialize::{CanonicalSerialize, Compress};
use merlin::Transcript;

// Copied from https://github.com/kevaundray/plookup
// Implements functions to interact with Merlin Transcript
pub trait TranscriptProtocol {
    /// Append a `commitment` with the given `label`.
    fn append_commitment(&mut self, label: &'static [u8], comm: &Commitment<Bls12_381>);

    /// Compute a `label`ed challenge variable.
    fn challenge_scalar(&mut self, label: &'static [u8]) -> Fr;
}

impl TranscriptProtocol for Transcript {
    fn append_commitment(&mut self, label: &'static [u8], comm: &Commitment<Bls12_381>) {
        let mut bytes = Vec::new();
        comm.0
            .serialize_with_mode(&mut bytes, Compress::No)
            .unwrap();
        self.append_message(label, &bytes);
    }

    fn challenge_scalar(&mut self, label: &'static [u8]) -> Fr {
        use ark_ff::UniformRand;

        let mut buf = [0u8; 32];
        self.challenge_bytes(label, &mut buf);

        let mut rng = ark_std::test_rng();
        Fr::rand(&mut rng)
    }
}
