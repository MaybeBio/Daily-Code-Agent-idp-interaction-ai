//! Small deterministic PRNG used for split assignment.
//!
//! Splits must be reproducible across machines and across dependency upgrades.
//! General-purpose RNG crates do not promise that a given seed produces the
//! same shuffle across major versions, so the pipeline uses SplitMix64
//! (Steele, Lea and Flood, 2014), which is tiny, fully specified, and has
//! published reference outputs.

/// SplitMix64 generator.
#[derive(Debug, Clone)]
pub struct SplitMix64 {
    state: u64,
}

impl SplitMix64 {
    pub fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    pub fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Uniform integer in `0..bound` without modulo bias (rejection sampling).
    ///
    /// # Panics
    /// Panics if `bound == 0`.
    pub fn below(&mut self, bound: u64) -> u64 {
        assert!(bound > 0, "bound must be positive");
        // Largest multiple of `bound` that fits in u64; values at or above it are rejected.
        let zone = u64::MAX - (u64::MAX % bound);
        loop {
            let x = self.next_u64();
            if x < zone {
                return x % bound;
            }
        }
    }

    /// In-place Fisher-Yates shuffle.
    pub fn shuffle<T>(&mut self, items: &mut [T]) {
        for i in (1..items.len()).rev() {
            let j = self.below(i as u64 + 1) as usize;
            items.swap(i, j);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matches_reference_outputs() {
        // Reference values for seed 0 from the SplitMix64 reference implementation.
        let mut r = SplitMix64::new(0);
        assert_eq!(r.next_u64(), 0xE220_A839_7B1D_CDAF);
        assert_eq!(r.next_u64(), 0x6E78_9E6A_A1B9_65F4);
        assert_eq!(r.next_u64(), 0x06C4_5D18_8009_454F);
    }

    #[test]
    fn shuffle_is_deterministic_permutation() {
        let mut a: Vec<u32> = (0..100).collect();
        let mut b = a.clone();
        SplitMix64::new(42).shuffle(&mut a);
        SplitMix64::new(42).shuffle(&mut b);
        assert_eq!(a, b);
        let mut sorted = a.clone();
        sorted.sort();
        assert_eq!(sorted, (0..100).collect::<Vec<_>>());
        assert_ne!(a, (0..100).collect::<Vec<_>>());
    }

    #[test]
    fn below_stays_in_range() {
        let mut r = SplitMix64::new(7);
        for bound in [1u64, 2, 3, 10, 1 << 40] {
            for _ in 0..1000 {
                assert!(r.below(bound) < bound);
            }
        }
    }
}
