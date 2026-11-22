//! Convert library analysis results into protocol reports.

use crate::commands::fail;
use opendna_analysis::{backbone, dssp, interactions, properties, sasa, shape};
use opendna_app_protocol::{CommandError, Contact, ResidueAnnotation, SequenceReport, StructureReport};
use opendna_core::{ResidueId, Sequence};
use opendna_store::StructureFile;
use std::collections::HashMap;

pub fn sequence(seq: &Sequence) -> Result<SequenceReport, CommandError> {
    let p = properties::compute(seq).map_err(fail)?;
    let window = 9.min(seq.len());
    Ok(SequenceReport {
        length: p.length,
        molecular_weight: p.molecular_weight,
        isoelectric_point: p.isoelectric_point,
        charge_at_ph7: p.charge_at_ph7,
        gravy: p.gravy,
        aromaticity: p.aromaticity,
        aliphatic_index: p.aliphatic_index,
        extinction_reduced: p.extinction_reduced,
        extinction_cystines: p.extinction_cystines,
        hydropathy: properties::hydropathy_profile(seq, window).map_err(fail)?,
    })
}

fn contact(p: &interactions::ResiduePair) -> Contact {
    let label = |r: &interactions::ResidueRef| format!("{}:{}{}", r.chain, r.amino_acid.three_letter(), r.id);
    Contact { a: label(&p.a), b: label(&p.b), distance: p.distance }
}

pub fn structure(file: &StructureFile) -> Result<StructureReport, CommandError> {
    let structure = opendna_io::read_structure(&file.content).map_err(fail)?;
    let assignment = dssp::assign(&structure);
    let torsions = backbone::torsions(&assignment.residues);
    let area = sasa::compute(&structure, sasa::SasaParams::default()).map_err(fail)?;
    let relative: HashMap<(&str, ResidueId), Option<f64>> =
        area.residues.iter().map(|r| ((r.chain.as_str(), r.id), r.relative)).collect();

    let mut residues = Vec::new();
    let mut ca = Vec::new();
    for ((r, state), t) in assignment.residues.iter().zip(&assignment.states).zip(&torsions) {
        if r.chain != file.chain {
            continue;
        }
        ca.push(r.ca);
        residues.push(ResidueAnnotation {
            label: r.id.to_string(),
            number: r.id.number,
            amino_acid: r.amino_acid.one_letter(),
            dssp: state.code(),
            relative_sasa: relative.get(&(r.chain.as_str(), r.id)).copied().flatten(),
            phi: t.phi,
            psi: t.psi,
        });
    }

    Ok(StructureReport {
        chain: file.chain.clone(),
        radius_of_gyration: shape::radius_of_gyration(&ca).map_err(fail)?,
        total_sasa: area.total,
        disulfides: interactions::disulfides(&structure).iter().map(contact).collect(),
        salt_bridges: interactions::salt_bridges(&structure).iter().map(contact).collect(),
        residues,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use opendna_store::StructureFormat;

    fn crambin() -> StructureFile {
        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../crates/opendna-analysis/tests/data/1crn.pdb");
        StructureFile {
            chain: "A".into(),
            format: StructureFormat::Pdb,
            content: std::fs::read_to_string(path).unwrap(),
        }
    }

    #[test]
    fn structure_report_covers_chain() {
        let r = structure(&crambin()).unwrap();
        assert_eq!(r.residues.len(), 46);
        assert_eq!(r.disulfides.len(), 3);
        assert_eq!(r.disulfides[0].a, "A:CYS3");
        assert!(r.residues.iter().filter(|x| x.dssp == 'H').count() > 15);
        assert!(r.residues.iter().all(|x| x.relative_sasa.is_some()));
    }

    #[test]
    fn short_sequences_get_a_shorter_window() {
        let r = sequence(&Sequence::parse("MKT").unwrap()).unwrap();
        assert_eq!(r.hydropathy.len(), 1);
    }
}
