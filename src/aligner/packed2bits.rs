fn validate_wfa2_sequence_length(name: &str, len: usize) -> i32 {
    assert!(
        len <= i32::MAX as usize,
        "{name} logical length must fit in i32"
    );
    len as i32
}

fn packed2bits_min_bytes(logical_len: usize) -> usize {
    logical_len.div_ceil(4)
}

pub(crate) fn validate_packed2bits_sequence(
    name: &str,
    sequence: &[u8],
    logical_len: usize,
) -> i32 {
    let logical_len = validate_wfa2_sequence_length(name, logical_len);
    let min_bytes = packed2bits_min_bytes(logical_len as usize);
    assert!(
        sequence.len() >= min_bytes,
        "{name} packed 2-bit buffer is too short for logical length"
    );
    logical_len
}

/// Pack ASCII DNA into WFA2's 2-bit A/C/G/T layout.
///
/// Bases are encoded as A=0, C=1, G=2, and T=3. Base `i` is stored at bit
/// position `2 * (i % 4)` of byte `i / 4`, so earlier bases occupy lower-order
/// bits. Uppercase and lowercase A/C/G/T are accepted; any other byte panics.
pub fn pack_dna_2bits(sequence: &[u8]) -> Vec<u8> {
    let mut packed = vec![0; packed2bits_min_bytes(sequence.len())];
    for (i, &base) in sequence.iter().enumerate() {
        let encoded = match base {
            b'A' | b'a' => 0,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' => 3,
            _ => panic!("invalid DNA base for 2-bit packing: 0x{base:02X}"),
        };
        packed[i / 4] |= encoded << (2 * (i % 4));
    }
    packed
}
