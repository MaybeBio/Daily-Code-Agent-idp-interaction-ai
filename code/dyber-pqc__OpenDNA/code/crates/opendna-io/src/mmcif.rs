//! mmCIF (`PDBx`) coordinate reader.
//!
//! Only the `_atom_site` loop is read. Author chain ids and residue numbers
//! (`auth_*`) are preferred over `label_*` so numbering matches the PDB
//! format and the literature. The first model is kept.

use crate::builder::{AtomRecord, StructureBuilder};
use crate::IoError;
use opendna_core::{Atom, ResidueId, Structure, Vec3};

/// Parse mmCIF text.
///
/// # Errors
/// [`IoError::Empty`] if there is no `_atom_site` loop, or
/// [`IoError::Malformed`] for missing required columns or bad numbers.
pub fn read(content: &str) -> Result<Structure, IoError> {
    let (header_line, columns, values) = find_atom_site_loop(content)?;
    let cols = Columns::resolve(&columns, header_line)?;
    let width = columns.len();
    if values.len() % width != 0 {
        return Err(IoError::malformed(header_line, "atom_site loop has a partial row"));
    }

    let mut builder = StructureBuilder::default();
    let mut first_model: Option<&str> = None;
    for (row_idx, row) in values.chunks(width).enumerate() {
        let line = row[0].1;
        let get = |i: usize| row[i].0;
        if let Some(m) = cols.model {
            let model = get(m);
            match first_model {
                None => first_model = Some(model),
                Some(first) if first != model => break,
                Some(_) => {}
            }
        }
        let num = |i: usize, what: &str| {
            get(i).parse::<f64>().map_err(|_| IoError::malformed(line, format!("invalid {what}")))
        };
        let optional = |i: Option<usize>| i.map(get).filter(|v| !is_null(v));

        let residue_number = optional(cols.auth_seq)
            .or_else(|| optional(Some(cols.label_seq)))
            // Waters and ligands have no label_seq_id; give them a stable row-based number.
            .map_or(Ok(i32::try_from(row_idx).unwrap_or(i32::MAX)), |v| {
                v.parse::<i32>().map_err(|_| IoError::malformed(line, "invalid residue number"))
            })?;

        builder.push(AtomRecord {
            hetero: get(cols.group) == "HETATM",
            alt_loc: optional(cols.alt).and_then(|v| v.chars().next()),
            chain: optional(cols.auth_asym).unwrap_or_else(|| get(cols.label_asym)).to_string(),
            residue_name: get(cols.comp).to_string(),
            residue_id: ResidueId {
                number: residue_number,
                insertion: optional(cols.ins).and_then(|v| v.chars().next()),
            },
            atom: Atom {
                serial: get(cols.id).parse().unwrap_or(0),
                name: get(cols.atom_name).to_string(),
                element: get(cols.element).to_ascii_uppercase(),
                position: Vec3::new(
                    num(cols.x, "Cartn_x")?,
                    num(cols.y, "Cartn_y")?,
                    num(cols.z, "Cartn_z")?,
                ),
                occupancy: optional(cols.occupancy).and_then(|v| v.parse().ok()).unwrap_or(1.0),
                b_factor: optional(cols.b_factor).and_then(|v| v.parse().ok()).unwrap_or(0.0),
            },
        });
    }

    if builder.is_empty() {
        return Err(IoError::Empty("atom_site rows"));
    }
    Ok(builder.finish())
}

fn is_null(v: &str) -> bool {
    v == "." || v == "?"
}

struct Columns {
    group: usize,
    id: usize,
    element: usize,
    atom_name: usize,
    alt: Option<usize>,
    comp: usize,
    label_asym: usize,
    auth_asym: Option<usize>,
    label_seq: usize,
    auth_seq: Option<usize>,
    ins: Option<usize>,
    x: usize,
    y: usize,
    z: usize,
    occupancy: Option<usize>,
    b_factor: Option<usize>,
    model: Option<usize>,
}

impl Columns {
    fn resolve(names: &[&str], line: usize) -> Result<Self, IoError> {
        let find = |n: &str| names.iter().position(|c| *c == n);
        let need = |n: &str| {
            find(n).ok_or_else(|| IoError::malformed(line, format!("missing column _atom_site.{n}")))
        };
        Ok(Self {
            group: need("group_PDB")?,
            id: need("id")?,
            element: need("type_symbol")?,
            atom_name: need("label_atom_id")?,
            alt: find("label_alt_id"),
            comp: need("label_comp_id")?,
            label_asym: need("label_asym_id")?,
            auth_asym: find("auth_asym_id"),
            label_seq: need("label_seq_id")?,
            auth_seq: find("auth_seq_id"),
            ins: find("pdbx_PDB_ins_code"),
            x: need("Cartn_x")?,
            y: need("Cartn_y")?,
            z: need("Cartn_z")?,
            occupancy: find("occupancy"),
            b_factor: find("B_iso_or_equiv"),
            model: find("pdbx_PDB_model_num"),
        })
    }
}

type Token<'a> = (&'a str, usize);

/// Returns the line of the loop header, the column names (without the
/// `_atom_site.` prefix) and every value token with its line number.
fn find_atom_site_loop(content: &str) -> Result<(usize, Vec<&str>, Vec<Token<'_>>), IoError> {
    let lines: Vec<&str> = content.lines().collect();
    let mut i = 0;
    while i < lines.len() {
        if lines[i].trim() == "loop_"
            && lines.get(i + 1).is_some_and(|l| l.trim_start().starts_with("_atom_site."))
        {
            let header_line = i + 1;
            i += 1;
            let mut columns = Vec::new();
            while let Some(name) = lines.get(i).and_then(|l| l.trim().strip_prefix("_atom_site.")) {
                columns.push(name.trim());
                i += 1;
            }
            let mut values = Vec::new();
            while i < lines.len() {
                let line = lines[i];
                let t = line.trim_start();
                if t.starts_with('_')
                    || t.starts_with("loop_")
                    || t.starts_with("data_")
                    || t.starts_with('#')
                {
                    break;
                }
                if let Some(first) = line.strip_prefix(';') {
                    // Multi-line text field, terminated by a line that is just ';'.
                    let start = i + 1;
                    let mut text_end = i + 1;
                    while text_end < lines.len() && !lines[text_end].starts_with(';') {
                        text_end += 1;
                    }
                    // The value text is not needed for coordinates; keep the first line as a token.
                    values.push((first.trim(), start));
                    i = text_end + 1;
                    continue;
                }
                tokenize(line, i + 1, &mut values).map_err(|m| IoError::malformed(i + 1, m))?;
                i += 1;
            }
            return Ok((header_line, columns, values));
        }
        i += 1;
    }
    Err(IoError::Empty("_atom_site loop"))
}

/// Split one line into CIF value tokens. Quoted values end at a matching quote
/// that is followed by whitespace or end of line, so `"O5'"` is one token.
fn tokenize<'a>(line: &'a str, lineno: usize, out: &mut Vec<Token<'a>>) -> Result<(), String> {
    let bytes = line.as_bytes();
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i].is_ascii_whitespace() {
            i += 1;
            continue;
        }
        let quote = bytes[i];
        if quote == b'\'' || quote == b'"' {
            let start = i + 1;
            let mut j = start;
            loop {
                if j >= bytes.len() {
                    return Err("unterminated quoted value".into());
                }
                if bytes[j] == quote && (j + 1 == bytes.len() || bytes[j + 1].is_ascii_whitespace()) {
                    break;
                }
                j += 1;
            }
            out.push((&line[start..j], lineno));
            i = j + 1;
        } else {
            let start = i;
            while i < bytes.len() && !bytes[i].is_ascii_whitespace() {
                i += 1;
            }
            out.push((&line[start..i], lineno));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE: &str = r#"data_TEST
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_asym_id
_atom_site.pdbx_PDB_model_num
ATOM   1 N N   . MET A 1 1 ? 27.340 24.430 2.614 1.00 9.67 1 A 1
ATOM   2 C CA  . MET A 1 1 ? 26.266 25.413 2.842 1.00 10.38 1 A 1
ATOM   3 C CA  . GLN A 1 2 ? 26.913 26.639 -3.198 1.00 9.62 2 A 1
HETATM 4 O "O5'" . HOH C 2 . ? 1.0 2.0 3.0 1.00 30.0 201 A 1
ATOM   5 N N   . MET A 1 1 ? 0.0 0.0 0.0 1.00 9.67 1 A 2
#
"#;

    #[test]
    fn reads_first_model_with_auth_numbering() {
        let s = read(SAMPLE).unwrap();
        assert_eq!(s.atom_count(), 4);
        let chain = s.chain("A").unwrap();
        assert_eq!(chain.sequence().to_string(), "MQ");
        let water = chain.residues.last().unwrap();
        assert!(water.is_water());
        assert_eq!(water.id.number, 201);
        assert_eq!(water.atoms[0].name, "O5'");
    }

    #[test]
    fn missing_loop_is_empty() {
        assert_eq!(read("data_X\n_cell.length_a 10\n"), Err(IoError::Empty("_atom_site loop")));
    }

    #[test]
    fn tokenizer_handles_embedded_quotes() {
        let mut out = Vec::new();
        tokenize(r#"a 'it's' "x y" b"#, 1, &mut out).unwrap();
        let vals: Vec<_> = out.iter().map(|t| t.0).collect();
        assert_eq!(vals, ["a", "it's", "x y", "b"]);
    }
}
