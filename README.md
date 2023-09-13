# Range Proof
The code in this repository is a implementation of the range proof protocol described [here](https://hackmd.io/@dabo/B1U4kx8XI).

This work is only proof of concept, it is not audited or does not come with any claims.

# Example
* To prove that `z = 100` lies in the range `[0, 2^8)`
* Trusted setup
```rust
use commitment_scheme;
// range: [0, 2^n)
let n = 8usize;
// trusted setup
let (pk, vk) = commitment_scheme::trusted_setup(4usize * n).unwrap();
```
* [Prover] Create range proof
```rust
use ark_bls12_381::Fr;
use merlin::Transcript;
// number in the above range
let z = Fr::from(100u8);
// merlin transcript will be used to
// transform an interactive protocol
// into a non-interactive protocol
let mut proof_transcript = Transcript::new(b"range_proof");
// generate range proof
let proof = RangeProof::prove(&pk, n, &z, &mut proof_transcript);
```
* [Verifier] Verify range proof
```rust
// verification transcript
let mut verification_transcript = Transcript::new(b"range_proof");
// verify the above range proof
let result = RangeProof::verify(&proof, &vk, n, &mut verification_transcript);
// assert that the result is ok
assert!(result.is_ok());
```

# License
[In detail here](https://github.com/roynalnaruto/range_proof/blob/master/LICENSE.md)
