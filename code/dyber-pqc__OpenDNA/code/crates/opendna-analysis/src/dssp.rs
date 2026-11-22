//! Secondary-structure assignment following Kabsch & Sander (1983), "DSSP".
//!
//! Implements the electrostatic H-bond model, n-turns, α/3₁₀/π helices,
//! β-bridges and ladders, and bends. Not implemented: β-bulge linking of
//! ladders (residues inside a bulge are reported as loop rather than `E`)
//! and polyproline helices. Results otherwise follow the classic DSSP
//! precedence `H > B > E > G > I > T > S`.

use crate::backbone::{extract, BackboneResidue};
use opendna_core::{Structure, Vec3};
use serde::Serialize;

/// DSSP secondary-structure classes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize)]
pub enum SecondaryStructure {
    /// α-helix (`H`)
    AlphaHelix,
    /// Isolated β-bridge (`B`)
    Bridge,
    /// Extended strand in a ladder (`E`)
    Strand,
    /// 3₁₀ helix (`G`)
    Helix310,
    /// π helix (`I`)
    PiHelix,
    /// H-bonded turn (`T`)
    Turn,
    /// Bend (`S`)
    Bend,
    /// None of the above (`-`)
    Loop,
}

impl SecondaryStructure {
    /// The DSSP one-letter code.
    pub fn code(self) -> char {
        match self {
            Self::AlphaHelix => 'H',
            Self::Bridge => 'B',
            Self::Strand => 'E',
            Self::Helix310 => 'G',
            Self::PiHelix => 'I',
            Self::Turn => 'T',
            Self::Bend => 'S',
            Self::Loop => '-',
        }
    }
}

/// Per-residue assignment, aligned with [`crate::backbone::extract`].
#[derive(Debug, Clone, Serialize)]
pub struct Assignment {
    /// Residues that were assigned (those with complete backbones).
    pub residues: Vec<BackboneResidue>,
    /// Secondary structure for each residue.
    pub states: Vec<SecondaryStructure>,
}

impl Assignment {
    /// DSSP string, one character per residue.
    pub fn to_code_string(&self) -> String {
        self.states.iter().map(|s| s.code()).collect()
    }
}

/// Coupling constant q1·q2·f in kcal/mol·Å (0.42e · 0.20e · 332).
const COUPLING: f64 = 0.084 * 332.0;
/// Energy threshold for an H-bond, kcal/mol.
const HBOND_MAX: f64 = -0.5;
/// Lowest energy DSSP reports, kcal/mol.
const ENERGY_FLOOR: f64 = -9.9;
/// Residue pairs with Cα farther apart than this cannot H-bond.
const CA_CUTOFF: f64 = 9.0;

/// Assign secondary structure to every amino-acid residue with a complete backbone.
pub fn assign(structure: &Structure) -> Assignment {
    let residues = extract(structure);
    let states = assign_backbone(&residues);
    Assignment { residues, states }
}

fn assign_backbone(res: &[BackboneResidue]) -> Vec<SecondaryStructure> {
    let n = res.len();
    let hbonds = HBonds::compute(res);
    // `no_break(a, b)`: residues a..=b are one continuous chain segment.
    let mut segment = vec![0usize; n];
    for i in 1..n {
        segment[i] = if res[i].linked_to_previous { segment[i - 1] } else { segment[i - 1] + 1 };
    }
    let no_break = |a: usize, b: usize| a <= b && b < n && segment[a] == segment[b];

    // n-turns: turn[k][i] means an H-bond from CO(i) to NH(i+n), n = 3, 4, 5.
    let mut turn = [vec![false; n], vec![false; n], vec![false; n]];
    for (row, len) in turn.iter_mut().zip([3usize, 4, 5]) {
        for (i, t) in row.iter_mut().enumerate().take(n.saturating_sub(len)) {
            *t = no_break(i, i + len) && hbonds.is_bond(i, i + len);
        }
    }

    // Helices: two consecutive n-turns at i-1 and i make residues i..i+n-1 helical.
    let helix = |k: usize, len: usize| {
        let mut h = vec![false; n];
        for i in 1..n.saturating_sub(len) {
            if turn[k][i - 1] && turn[k][i] {
                h[i..i + len].iter_mut().for_each(|x| *x = true);
            }
        }
        h
    };
    let alpha = helix(1, 4);
    let h310 = helix(0, 3);
    let pi = helix(2, 5);

    // β-bridges.
    let hb = |a: usize, b: usize| hbonds.is_bond(a, b);
    let mut bridge_partners: Vec<Vec<(usize, bool)>> = vec![Vec::new(); n]; // (partner, parallel)
    for i in 1..n.saturating_sub(1) {
        for j in (i + 3)..n.saturating_sub(1) {
            if !(no_break(i - 1, i + 1) && no_break(j - 1, j + 1)) {
                continue;
            }
            let parallel = (hb(i - 1, j) && hb(j, i + 1)) || (hb(j - 1, i) && hb(i, j + 1));
            let antiparallel = (hb(i, j) && hb(j, i)) || (hb(i - 1, j + 1) && hb(j - 1, i + 1));
            if parallel || antiparallel {
                bridge_partners[i].push((j, parallel));
                bridge_partners[j].push((i, parallel));
            }
        }
    }
    // A bridge is part of a ladder if an adjacent residue bridges the adjacent partner.
    let in_ladder = |i: usize| {
        bridge_partners[i].iter().any(|&(j, parallel)| {
            let step: isize = if parallel { 1 } else { -1 };
            [-1isize, 1].iter().any(|&d| {
                let (ni, nj) = (i as isize + d, j as isize + d * step);
                ni >= 0
                    && nj >= 0
                    && (ni as usize) < n
                    && bridge_partners[ni as usize].contains(&(nj as usize, parallel))
            })
        })
    };

    // Bend: Cα(i-2)→Cα(i) vs Cα(i)→Cα(i+2) angle above 70°.
    let bend = |i: usize| {
        i >= 2 && no_break(i - 2, i + 2) && {
            let u = res[i].ca - res[i - 2].ca;
            let v = res[i + 2].ca - res[i].ca;
            Vec3::angle(u, Vec3::default(), v) > 70.0
        }
    };
    let in_turn = |i: usize| {
        [3usize, 4, 5].iter().enumerate().any(|(k, &len)| (1..len).any(|off| i >= off && turn[k][i - off]))
    };

    (0..n)
        .map(|i| {
            if alpha[i] {
                SecondaryStructure::AlphaHelix
            } else if !bridge_partners[i].is_empty() && !in_ladder(i) {
                SecondaryStructure::Bridge
            } else if !bridge_partners[i].is_empty() {
                SecondaryStructure::Strand
            } else if h310[i] {
                SecondaryStructure::Helix310
            } else if pi[i] {
                SecondaryStructure::PiHelix
            } else if in_turn(i) {
                SecondaryStructure::Turn
            } else if bend(i) {
                SecondaryStructure::Bend
            } else {
                SecondaryStructure::Loop
            }
        })
        .collect()
}

/// For each donor N-H, the two strongest acceptor C=O partners.
struct HBonds {
    acceptors: Vec<[(usize, f64); 2]>,
}

impl HBonds {
    fn compute(res: &[BackboneResidue]) -> Self {
        let n = res.len();
        // Amide H placed 1 Å from N, opposite the previous carbonyl (DSSP convention).
        let hydrogens: Vec<Option<Vec3>> = (0..n)
            .map(|i| {
                let r = &res[i];
                let is_pro = r.amino_acid == opendna_core::AminoAcid::Pro;
                (i > 0 && r.linked_to_previous && !is_pro)
                    .then(|| r.n + (res[i - 1].c - res[i - 1].o).normalized())
            })
            .collect();

        let mut acceptors = vec![[(usize::MAX, 0.0), (usize::MAX, 0.0)]; n];
        let c2 = CA_CUTOFF * CA_CUTOFF;
        for (donor, h) in hydrogens.iter().enumerate() {
            let Some(h) = *h else { continue };
            let nd = res[donor].n;
            for (acc, ra) in res.iter().enumerate() {
                // An NH cannot bond to its own residue's or the preceding residue's
                // carbonyl; DSSP skips both pairs.
                if acc == donor || acc + 1 == donor || res[donor].ca.distance_squared(ra.ca) >= c2 {
                    continue;
                }
                let e = COUPLING
                    * (1.0 / ra.o.distance(nd) + 1.0 / ra.c.distance(h)
                        - 1.0 / ra.o.distance(h)
                        - 1.0 / ra.c.distance(nd));
                let e = e.max(ENERGY_FLOOR);
                let best = &mut acceptors[donor];
                if e < best[0].1 {
                    best[1] = best[0];
                    best[0] = (acc, e);
                } else if e < best[1].1 {
                    best[1] = (acc, e);
                }
            }
        }
        Self { acceptors }
    }

    /// `true` if C=O of `acceptor` is H-bonded to N-H of `donor`.
    fn is_bond(&self, acceptor: usize, donor: usize) -> bool {
        self.acceptors
            .get(donor)
            .is_some_and(|best| best.iter().any(|&(a, e)| a == acceptor && e < HBOND_MAX))
    }
}
