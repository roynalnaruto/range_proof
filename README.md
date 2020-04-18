# Range Proof
The code in this repository is a **WIP** implementation of the range proof protocol described [here](https://hackmd.io/@dabo/B1U4kx8XI).

This work is only proof of concept, and its purpose is for me to learn.

# Example
* Wish to prove that `z = 100` lies in the range `[0, 2^8)`
* Trusted setup
```rust
let n = 8usize;
let (pk, vk) = trusted_setup(n).unwrap();
```
* [Prover] Create range proof
```rust
let z = Fr::from(100u8);
let mut transcript = Transcript::new(b"range_proof");
let range_proof = RangeProof::prove(&pk, n, &z, &transcript);
```
