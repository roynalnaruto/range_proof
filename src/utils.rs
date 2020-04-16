use bls12_381::Scalar;

use std::convert::TryInto;

pub fn to_raw_bytes(scalar: &Scalar) -> [u64; 4] {
    let s_bytes = scalar.to_bytes();
    let s_chunks = s_bytes.chunks(8);

    let mut out = [0u64; 4];

    for (s_chunk, i) in s_chunks.into_iter().zip(0..4) {
        out[i] = u64::from_le_bytes(*slice_to_array8(s_chunk));
    }

    out
}

pub fn slice_to_array8(slice: &[u8]) -> &[u8; 8] {
    slice.try_into().expect("slice with incorrect length")
}

pub fn slice_to_array32(slice: &[u8]) -> &[u8; 32] {
    slice.try_into().expect("slice with incorrect length")
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_to_raw_bytes() {
        let three = Scalar::from(3 as u64);
        let three_raw_bytes = to_raw_bytes(&three);
        assert_eq!(three_raw_bytes, [3, 0, 0, 0]);

        let expected_raw_bytes = [
            0xffffffff00000000,
            0x53bda402fffe5bfe,
            0x3339d80809a1d805,
            0x73eda753299d7d48,
        ];
        let scalar = Scalar::from_raw(expected_raw_bytes);
        let raw_bytes = to_raw_bytes(&scalar);
        assert_eq!(raw_bytes, expected_raw_bytes);
    }
}
