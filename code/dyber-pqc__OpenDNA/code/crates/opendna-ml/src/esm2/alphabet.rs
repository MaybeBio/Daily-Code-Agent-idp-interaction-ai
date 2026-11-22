//! The fixed 33-token ESM-2 alphabet.
//!
//! ESM-2 does not use a learned tokenizer: the vocabulary is a fixed list of
//! amino-acid letters plus special tokens, so it is written out here rather
//! than loaded from a `tokenizer.json`. The order matches
//! `facebook/esm2_*_UR50D`'s vocabulary exactly and must not be reordered —
//! the indices are the embedding rows.

/// Vocabulary in index order.
pub const TOKENS: [&str; 33] = [
    "<cls>", "<pad>", "<eos>", "<unk>", "L", "A", "G", "V", "S", "E", "R", "T", "I", "D", "P", "K", "Q", "N",
    "F", "Y", "M", "H", "W", "C", "X", "B", "U", "Z", "O", ".", "-", "<null_1>", "<mask>",
];

/// Start-of-sequence token prepended to every sequence.
pub const CLS: u32 = 0;
/// Padding token.
pub const PAD: u32 = 1;
/// End-of-sequence token appended to every sequence.
pub const EOS: u32 = 2;
/// Token for residues outside the alphabet.
pub const UNK: u32 = 3;
/// Mask token. Used by masked-token scoring, not by plain embedding.
pub const MASK: u32 = 32;

/// Token id for a one-letter residue code (case-insensitive).
///
/// Unrecognised characters map to [`UNK`]. `X` has its own id and is *not*
/// the same as `<unk>`: it means "unknown residue" in the sequence rather
/// than "unknown token".
pub fn token_id(residue: char) -> u32 {
    let upper = residue.to_ascii_uppercase();
    let mut buf = [0u8; 4];
    let s = upper.encode_utf8(&mut buf);
    // Only single-character entries can match a residue; the specials are
    // longer, so this cannot accidentally match "<cls>" and friends.
    TOKENS.iter().position(|t| *t == s).map_or(UNK, |i| i as u32)
}

/// Encode a sequence as `<cls> … <eos>`.
pub fn encode(sequence: &str) -> Vec<u32> {
    let mut out = Vec::with_capacity(sequence.len() + 2);
    out.push(CLS);
    out.extend(sequence.chars().filter(|c| !c.is_whitespace()).map(token_id));
    out.push(EOS);
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn known_residues_match_the_reference_vocabulary() {
        assert_eq!(token_id('L'), 4);
        assert_eq!(token_id('A'), 5);
        assert_eq!(token_id('M'), 20);
        assert_eq!(token_id('C'), 23);
        assert_eq!(token_id('X'), 24);
        assert_eq!(token_id('l'), 4, "lower case must fold to the same id");
    }

    #[test]
    fn unknown_characters_become_unk() {
        assert_eq!(token_id('@'), UNK);
        assert_eq!(token_id('J'), UNK);
    }

    #[test]
    fn encoding_brackets_the_sequence() {
        // Matches the reference tokenisation of "M": [0, 20, 2].
        assert_eq!(encode("M"), vec![CLS, 20, EOS]);
        assert_eq!(encode("MQ IF"), vec![CLS, 20, 16, 12, 18, EOS]);
    }
}
