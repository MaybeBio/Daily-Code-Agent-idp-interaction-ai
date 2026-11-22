//! PDB format (fixed-column, version 3.3).
//!
//! Reads `ATOM`/`HETATM` records from the first `MODEL` only. Header records
//! are ignored. Columns follow the wwPDB 3.3 specification; short lines are
//! accepted as long as they reach the end of the z coordinate (column 54).

use crate::builder::{infer_element, AtomRecord, StructureBuilder};
use crate::IoError;
use opendna_core::{Atom, ResidueId, Structure, Vec3};
use std::fmt::Write;

/// Parse PDB text.
///
/// # Errors
/// [`IoError::Malformed`] for an unparsable coordinate or serial, and
/// [`IoError::Empty`] if there are no atom records.
pub fn read(content: &str) -> Result<Structure, IoError> {
    let mut builder = StructureBuilder::default();
    let mut seen_model = false;

    for (idx, line) in content.lines().enumerate() {
        let lineno = idx + 1;
        let record = line.get(..6).unwrap_or(line);
        match record.trim_end() {
            "MODEL" => {
                if seen_model {
                    break;
                }
                seen_model = true;
            }
            "ENDMDL" => break,
            "ATOM" | "HETATM" => builder.push(parse_atom(line, lineno)?),
            _ => {}
        }
    }

    if builder.is_empty() {
        return Err(IoError::Empty("ATOM/HETATM records"));
    }
    Ok(builder.finish())
}

/// Columns are 1-based and inclusive in the spec; this takes them as such.
fn col(line: &str, from: usize, to: usize) -> &str {
    line.get(from - 1..to.min(line.len())).unwrap_or("")
}

fn parse_atom(line: &str, lineno: usize) -> Result<AtomRecord, IoError> {
    if !line.is_ascii() {
        return Err(IoError::malformed(lineno, "non-ASCII characters in atom record"));
    }
    if line.len() < 54 {
        return Err(IoError::malformed(lineno, "atom record shorter than 54 columns"));
    }
    let coord = |from, to, axis| {
        col(line, from, to)
            .trim()
            .parse::<f64>()
            .map_err(|_| IoError::malformed(lineno, format!("invalid {axis} coordinate")))
    };
    let hetero = line.starts_with("HETATM");
    let name = col(line, 13, 16).trim().to_string();
    let element = match col(line, 77, 78).trim() {
        "" => infer_element(&name, hetero),
        e => e.to_ascii_uppercase(),
    };
    let number = col(line, 23, 26)
        .trim()
        .parse::<i32>()
        .map_err(|_| IoError::malformed(lineno, "invalid residue number"))?;

    Ok(AtomRecord {
        hetero,
        alt_loc: col(line, 17, 17).chars().next().filter(|c| *c != ' '),
        chain: col(line, 22, 22).trim().to_string(),
        residue_name: col(line, 18, 20).trim().to_string(),
        residue_id: ResidueId { number, insertion: col(line, 27, 27).chars().next().filter(|c| *c != ' ') },
        atom: Atom {
            // Serial numbers overflow to hybrid-36 or "*****" in very large files;
            // they are not used for identity, so fall back to 0.
            serial: col(line, 7, 11).trim().parse().unwrap_or(0),
            name,
            element,
            position: Vec3::new(coord(31, 38, "x")?, coord(39, 46, "y")?, coord(47, 54, "z")?),
            occupancy: col(line, 55, 60).trim().parse().unwrap_or(1.0),
            b_factor: col(line, 61, 66).trim().parse().unwrap_or(0.0),
        },
    })
}

/// Serialise a structure as PDB text. Atoms are renumbered from 1 and a
/// `TER` record is written after each polymer chain.
pub fn write(structure: &Structure) -> String {
    let mut out = String::new();
    let mut serial = 1u32;
    for chain in &structure.chains {
        let chain_char = chain.id.chars().next().unwrap_or(' ');
        let mut last_polymer = None;
        for residue in &chain.residues {
            let record = if residue.hetero { "HETATM" } else { "ATOM  " };
            for atom in &residue.atoms {
                // Names shorter than 4 characters start in column 14 unless the
                // element symbol is two letters (e.g. " CA " vs "FE  ").
                let name = if atom.name.len() >= 4 || atom.element.len() == 2 {
                    format!("{:<4}", atom.name)
                } else {
                    format!(" {:<3}", atom.name)
                };
                let _ = writeln!(
                    out,
                    "{record}{serial:>5} {name} {res:>3} {chain_char}{num:>4}{ins}   {x:>8.3}{y:>8.3}{z:>8.3}{occ:>6.2}{b:>6.2}          {el:>2}",
                    res = residue.name,
                    num = residue.id.number,
                    ins = residue.id.insertion.unwrap_or(' '),
                    x = atom.position.x,
                    y = atom.position.y,
                    z = atom.position.z,
                    occ = atom.occupancy,
                    b = atom.b_factor,
                    el = atom.element,
                );
                serial += 1;
            }
            if !residue.hetero {
                last_polymer = Some(residue);
            }
        }
        if let Some(r) = last_polymer {
            let _ = writeln!(out, "TER   {serial:>5}      {:>3} {chain_char}{:>4}", r.name, r.id.number);
            serial += 1;
        }
    }
    out.push_str("END\n");
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    const TWO_RESIDUES: &str = "\
ATOM      1  N   ALA A   1      11.104   6.134  -6.504  1.00  0.00           N
ATOM      2  CA  ALA A   1      11.639   6.071  -5.147  1.00  0.00           C
ATOM      3  N   GLY A   2      12.000   7.000  -4.000  1.00 10.00           N
ATOM      4  CA AGLY A   2      12.500   7.500  -3.500  0.60 10.00           C
ATOM      5  CA BGLY A   2      12.600   7.600  -3.600  0.40 10.00           C
HETATM    6 ZN    ZN A 101       0.000   0.000   0.000  1.00 20.00
";

    #[test]
    fn reads_hierarchy_and_drops_second_altloc() {
        let s = read(TWO_RESIDUES).unwrap();
        assert_eq!(s.chains.len(), 1);
        let chain = &s.chains[0];
        assert_eq!(chain.residues.len(), 3);
        assert_eq!(chain.residues[1].atoms.len(), 2);
        assert!((chain.residues[1].atoms[1].position.x - 12.5).abs() < 1e-9);
        assert_eq!(chain.sequence().to_string(), "AG");
    }

    #[test]
    fn infers_element_for_hetero_ion() {
        let s = read(TWO_RESIDUES).unwrap();
        let zn = &s.chains[0].residues[2];
        assert!(zn.hetero);
        assert_eq!(zn.atoms[0].element, "ZN");
    }

    #[test]
    fn stops_after_first_model() {
        let text = format!("MODEL        1\n{TWO_RESIDUES}ENDMDL\nMODEL        2\n{TWO_RESIDUES}ENDMDL\n");
        assert_eq!(read(&text).unwrap().atom_count(), 5);
    }

    #[test]
    fn bad_coordinate_reports_line() {
        let text = "ATOM      1  N   ALA A   1      11.104   abcde  -6.504  1.00  0.00           N\n";
        let err = read(text).unwrap_err();
        assert!(matches!(err, IoError::Malformed { at, .. } if at.line == 1), "{err}");
    }

    #[test]
    fn empty_input_is_an_error() {
        assert_eq!(read("HEADER    NOTHING\n"), Err(IoError::Empty("ATOM/HETATM records")));
    }

    #[test]
    fn write_then_read_round_trips() {
        let s = read(TWO_RESIDUES).unwrap();
        let again = read(&write(&s)).unwrap();
        assert_eq!(again.atom_count(), s.atom_count());
        for (a, b) in s.atoms().zip(again.atoms()) {
            assert_eq!(a.name, b.name);
            assert_eq!(a.element, b.element);
            assert!(a.position.distance(b.position) < 1e-3);
        }
    }
}
