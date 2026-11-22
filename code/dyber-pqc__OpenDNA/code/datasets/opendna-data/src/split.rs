//! Group-aware, seeded split assignment.
//!
//! Records are never split individually. Every record belongs to a group
//! (an MMseqs2 cluster, or a connected component of clusters and CATH
//! topologies), and whole groups are assigned to one split. This is what keeps
//! homologous sequences from appearing on both sides of a train/test boundary.

use std::collections::BTreeMap;
use std::fmt;

use anyhow::{bail, Result};
use serde::Serialize;

use crate::rng::SplitMix64;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum Split {
    Train,
    Validation,
    Test,
}

impl Split {
    pub const ALL: [Split; 3] = [Split::Train, Split::Validation, Split::Test];

    pub fn as_str(self) -> &'static str {
        match self {
            Split::Train => "train",
            Split::Validation => "validation",
            Split::Test => "test",
        }
    }
}

impl fmt::Display for Split {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.pad(self.as_str())
    }
}

/// Target fraction of the total group weight for each split.
#[derive(Debug, Clone, Copy, Serialize)]
pub struct Fractions {
    pub train: f64,
    pub validation: f64,
    pub test: f64,
}

impl Fractions {
    pub fn new(train: f64, validation: f64, test: f64) -> Result<Self> {
        let f = Self { train, validation, test };
        for (name, v) in [("train", train), ("validation", validation), ("test", test)] {
            if !(0.0..=1.0).contains(&v) || !v.is_finite() {
                bail!("split fraction {name}={v} must be in [0, 1]");
            }
        }
        if ((train + validation + test) - 1.0).abs() > 1e-9 {
            bail!("split fractions must sum to 1, got {}", train + validation + test);
        }
        Ok(f)
    }

    fn get(&self, s: Split) -> f64 {
        match s {
            Split::Train => self.train,
            Split::Validation => self.validation,
            Split::Test => self.test,
        }
    }
}

/// Assigns each group to a split.
///
/// * `weights` maps group id to its weight (for example the number of records
///   or residues in the group). Groups with zero weight are still assigned.
/// * `fixed` pins specific groups to a split regardless of the fractions; this
///   is used to hold out external benchmark sets. Fixed groups do not count
///   toward the fractions, which apply to the remaining (free) weight only.
///
/// Free groups are visited in a seeded random order (sorted by id first, so the
/// result does not depend on map iteration order) and each goes to the split
/// whose current weight is furthest below its target. With many small groups
/// the realised fractions converge to the targets; a single very large group
/// can make them deviate, which is reported by [`SplitSummary`].
pub fn assign_groups(
    weights: &BTreeMap<String, u64>,
    fixed: &BTreeMap<String, Split>,
    fractions: Fractions,
    seed: u64,
) -> Result<BTreeMap<String, Split>> {
    for g in fixed.keys() {
        if !weights.contains_key(g) {
            bail!("fixed group {g} has no weight entry");
        }
    }
    let mut free: Vec<(&String, u64)> =
        weights.iter().filter(|(g, _)| !fixed.contains_key(*g)).map(|(g, w)| (g, *w)).collect();
    let total: u64 = free.iter().map(|(_, w)| *w).sum();
    SplitMix64::new(seed).shuffle(&mut free);

    let mut assigned: BTreeMap<String, Split> = fixed.clone();
    let mut current = [0u64; 3];
    for (group, w) in free {
        let mut best = None;
        let mut best_deficit = f64::NEG_INFINITY;
        for (i, s) in Split::ALL.iter().enumerate() {
            let target = fractions.get(*s);
            if target <= 0.0 {
                continue;
            }
            let deficit = target * total as f64 - current[i] as f64;
            // Strictly greater: ties go to the earlier split (train before validation before test).
            if deficit > best_deficit {
                best_deficit = deficit;
                best = Some(i);
            }
        }
        let i = best.expect("fractions sum to 1, so at least one is positive");
        current[i] += w;
        assigned.insert(group.clone(), Split::ALL[i]);
    }
    Ok(assigned)
}

/// Counts per split, used in dataset cards and manifests.
#[derive(Debug, Clone, Default, Serialize)]
pub struct SplitSummary {
    pub train: SplitCounts,
    pub validation: SplitCounts,
    pub test: SplitCounts,
}

#[derive(Debug, Clone, Copy, Default, Serialize)]
pub struct SplitCounts {
    pub rows: u64,
    pub groups: u64,
}

impl SplitSummary {
    pub fn get_mut(&mut self, s: Split) -> &mut SplitCounts {
        match s {
            Split::Train => &mut self.train,
            Split::Validation => &mut self.validation,
            Split::Test => &mut self.test,
        }
    }

    pub fn get(&self, s: Split) -> SplitCounts {
        match s {
            Split::Train => self.train,
            Split::Validation => self.validation,
            Split::Test => self.test,
        }
    }

    pub fn total_rows(&self) -> u64 {
        self.train.rows + self.validation.rows + self.test.rows
    }

    /// Builds a summary from a per-row group id iterator and the assignment.
    pub fn from_rows<'a>(
        row_groups: impl IntoIterator<Item = &'a str>,
        assignment: &BTreeMap<String, Split>,
    ) -> Self {
        let mut s = Self::default();
        let mut seen: BTreeMap<&str, ()> = BTreeMap::new();
        for g in row_groups {
            let split = assignment[g];
            s.get_mut(split).rows += 1;
            if seen.insert(g, ()).is_none() {
                s.get_mut(split).groups += 1;
            }
        }
        s
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn weights(n: usize) -> BTreeMap<String, u64> {
        (0..n).map(|i| (format!("g{i:04}"), 1 + (i as u64 % 5))).collect()
    }

    #[test]
    fn deterministic_for_seed() {
        let w = weights(500);
        let f = Fractions::new(0.8, 0.1, 0.1).unwrap();
        let a = assign_groups(&w, &BTreeMap::new(), f, 1).unwrap();
        let b = assign_groups(&w, &BTreeMap::new(), f, 1).unwrap();
        let c = assign_groups(&w, &BTreeMap::new(), f, 2).unwrap();
        assert_eq!(a, b);
        assert_ne!(a, c);
    }

    #[test]
    fn fractions_are_approximately_met() {
        let w = weights(2000);
        let total: u64 = w.values().sum();
        let f = Fractions::new(0.8, 0.1, 0.1).unwrap();
        let a = assign_groups(&w, &BTreeMap::new(), f, 9).unwrap();
        let mut per = [0u64; 3];
        for (g, s) in &a {
            per[*s as usize] += w[g];
        }
        let frac: Vec<f64> = per.iter().map(|x| *x as f64 / total as f64).collect();
        assert!((frac[0] - 0.8).abs() < 0.01, "{frac:?}");
        assert!((frac[1] - 0.1).abs() < 0.01, "{frac:?}");
        assert!((frac[2] - 0.1).abs() < 0.01, "{frac:?}");
    }

    #[test]
    fn every_group_assigned_once_and_fixed_respected() {
        let w = weights(100);
        let mut fixed = BTreeMap::new();
        fixed.insert("g0003".to_string(), Split::Test);
        fixed.insert("g0004".to_string(), Split::Test);
        let f = Fractions::new(0.9, 0.1, 0.0).unwrap();
        let a = assign_groups(&w, &fixed, f, 3).unwrap();
        assert_eq!(a.len(), w.len());
        assert_eq!(a["g0003"], Split::Test);
        assert_eq!(a["g0004"], Split::Test);
        // With a zero test fraction no free group may land in test.
        assert_eq!(a.values().filter(|s| **s == Split::Test).count(), 2);
    }

    #[test]
    fn rejects_bad_fractions_and_unknown_fixed() {
        assert!(Fractions::new(0.8, 0.1, 0.2).is_err());
        assert!(Fractions::new(-0.1, 0.6, 0.5).is_err());
        let mut fixed = BTreeMap::new();
        fixed.insert("missing".to_string(), Split::Test);
        let f = Fractions::new(0.8, 0.1, 0.1).unwrap();
        assert!(assign_groups(&weights(3), &fixed, f, 0).is_err());
    }

    #[test]
    fn summary_counts_rows_and_groups() {
        let mut a = BTreeMap::new();
        a.insert("x".to_string(), Split::Train);
        a.insert("y".to_string(), Split::Test);
        let s = SplitSummary::from_rows(["x", "x", "y"], &a);
        assert_eq!(s.train.rows, 2);
        assert_eq!(s.train.groups, 1);
        assert_eq!(s.test.rows, 1);
        assert_eq!(s.total_rows(), 3);
    }
}
