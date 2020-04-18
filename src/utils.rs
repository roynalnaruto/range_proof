use std::convert::TryInto;

use algebra::{bls12_381::Fr, to_bytes, ToBytes};
use bitvec::prelude::*;

macro_rules! slice_to_array {
    ($name: ident, $sty: ty) => {
        pub fn $name(slice: &[u8]) -> $sty {
            slice.try_into().expect("slice has incorrect length")
        }
    }
}

slice_to_array!(slice_to_array8, &[u8; 8]);
slice_to_array!(slice_to_array32, &[u8; 32]);

pub fn to_raw_bytes(in_bytes: &[u8; 32]) -> [u64; 4] {
    let chunks = in_bytes.chunks(8);

    let mut out = [0u64; 4];

    for (chunk, i) in chunks.into_iter().zip(0..4) {
        out[i] = u64::from_le_bytes(*slice_to_array8(chunk));
    }

    out
}

pub fn to_bits(z: &Fr) -> BitVec<Lsb0, u8> {
    let bytes = to_bytes![z].unwrap();

    BitVec::<Lsb0, u8>::from_slice(&bytes)
}

pub fn as_ref(z: &Fr) -> [u64; 4] {
    let z_bytes = to_bytes![z].unwrap();
    let slice_32 = slice_to_array32(&z_bytes.as_slice());

    to_raw_bytes(&slice_32)
}

#[cfg(test)]
mod test {
    use super::*;

    use bls12_381::Scalar;

    extern crate typename;
    use typename::TypeName;

    #[test]
    fn test_to_raw_bytes() {
        let three = Scalar::from(3 as u64);
        let three_raw_bytes = to_raw_bytes(&three.to_bytes());
        assert_eq!(three_raw_bytes, [3, 0, 0, 0]);

        let expected_raw_bytes = [
            0xffffffff00000000,
            0x53bda402fffe5bfe,
            0x3339d80809a1d805,
            0x73eda753299d7d48,
        ];
        let scalar = Scalar::from_raw(expected_raw_bytes);
        let raw_bytes = to_raw_bytes(&scalar.to_bytes());
        assert_eq!(raw_bytes, expected_raw_bytes);
    }

    #[test]
    fn test_slice_to_array() {
        let vector = vec![0u8, 1u8, 0u8, 1u8, 0u8, 1u8, 0u8, 1u8];
        let slice = vector.as_slice();
        let slice_8 = slice_to_array8(&slice);
        assert_eq!((*slice_8).type_name_of(), "[u8; 8]");

        let vector = vec![1u8; 32];
        let slice = vector.as_slice();
        let slice_32 = slice_to_array32(&slice);
        assert_eq!((*slice_32).type_name_of(), "[u8; 32]");
    }

    #[test]
    fn test_to_bits() {
        let bits = to_bits(&Fr::from(5u64));
        assert_eq!(bits[0], true);
        assert_eq!(bits[1], false);
        assert_eq!(bits[2], true);
        for bit in bits.iter().skip(3) {
            assert_eq!((*bit), false);
        }

        let bits = to_bits(&Fr::from(12u64));
        assert_eq!(bits[0], false);
        assert_eq!(bits[1], false);
        assert_eq!(bits[2], true);
        assert_eq!(bits[3], true);
        assert_eq!(bits[4], false);
        for bit in bits.iter().skip(5) {
            assert_eq!((*bit), false);
        }
    }
}
