//! Shared logic for turning a flat stream of atom records into the
//! chain → residue → atom hierarchy, applying model/altloc filtering.

use opendna_core::{Atom, Chain, Residue, ResidueId, Structure};

/// One atom record as read from any coordinate format.
pub(crate) struct AtomRecord {
    pub hetero: bool,
    pub alt_loc: Option<char>,
    pub chain: String,
    pub residue_name: String,
    pub residue_id: ResidueId,
    pub atom: Atom,
}

#[derive(Default)]
pub(crate) struct StructureBuilder {
    chains: Vec<Chain>,
    /// The altloc chosen for the current residue, so all of its atoms agree.
    chosen_alt: Option<char>,
}

impl StructureBuilder {
    pub fn push(&mut self, rec: AtomRecord) {
        let same_residue = self.chains.last().is_some_and(|c| {
            c.id == rec.chain
                && c.residues.last().is_some_and(|r| r.id == rec.residue_id && r.name == rec.residue_name)
        });

        if !same_residue {
            self.chosen_alt = None;
        }
        // Keep the first alternate location seen for each residue and drop the rest.
        if let Some(alt) = rec.alt_loc {
            match self.chosen_alt {
                None => self.chosen_alt = Some(alt),
                Some(chosen) if chosen != alt => return,
                Some(_) => {}
            }
        }

        if self.chains.last().is_none_or(|c| c.id != rec.chain) {
            self.chains.push(Chain { id: rec.chain.clone(), residues: Vec::new() });
        }
        let chain = self.chains.last_mut().expect("chain pushed above");

        if same_residue {
            let residue = chain.residues.last_mut().expect("checked by same_residue");
            // A duplicate atom name inside one residue means a second altloc
            // written without an altloc flag; keep the first.
            if residue.atom(&rec.atom.name).is_none() {
                residue.atoms.push(rec.atom);
            }
        } else {
            chain.residues.push(Residue {
                name: rec.residue_name,
                id: rec.residue_id,
                hetero: rec.hetero,
                atoms: vec![rec.atom],
            });
        }
    }

    pub fn finish(self) -> Structure {
        Structure { chains: self.chains }
    }

    pub fn is_empty(&self) -> bool {
        self.chains.is_empty()
    }
}

/// Infer an element symbol from an atom name when the element column is
/// blank (common in files written by older tools).
///
/// For polymer atoms the first letter is the element (`CA` is C-alpha,
/// `OG1` is oxygen). For hetero groups a name that is exactly a common ion or
/// halogen symbol is taken as that element.
pub(crate) fn infer_element(atom_name: &str, hetero: bool) -> String {
    const IONS: [&str; 12] = ["CL", "BR", "FE", "ZN", "MG", "MN", "CA", "NA", "CU", "CO", "NI", "SE"];
    let letters: String =
        atom_name.chars().filter(char::is_ascii_alphabetic).collect::<String>().to_ascii_uppercase();
    if hetero && IONS.contains(&letters.as_str()) {
        return letters;
    }
    letters.chars().next().map(String::from).unwrap_or_default()
}
