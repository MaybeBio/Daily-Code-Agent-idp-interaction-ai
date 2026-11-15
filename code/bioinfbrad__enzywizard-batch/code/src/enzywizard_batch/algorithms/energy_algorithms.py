from __future__ import annotations

from typing import Dict, Any

from openmm import Platform,Context,VerletIntegrator,LocalEnergyMinimizer,unit
from openmm.app import PDBFile, PDBxFile,Modeller,ForceField,NoCutoff,HBonds
from ..utils.logging_utils import Logger


def compute_energy_terms(struct: PDBFile | PDBxFile, logger: Logger, minimize_energy: bool = True, minimization_iteration:int = 100,force_field_file="charmm36.xml") -> Dict[str, float] | None:


    terms = [
        "total_potential_energy",
        "harmonic_bond_force",
        "harmonic_angle_force",
        "custom_bond_force",
        "custom_torsion_force",
        "custom_nonbonded_force",
        "nonbonded_force",
        "periodic_torsion_force",
        "cmap_torsion_force",
    ]

    mapping = {
        "total_potential_energy": "TotalPotentialEnergy",
        "harmonic_bond_force": "HarmonicBondForce",
        "harmonic_angle_force": "HarmonicAngleForce",
        "custom_bond_force": "CustomBondForce",
        "custom_torsion_force": "CustomTorsionForce",
        "custom_nonbonded_force": "CustomNonbondedForce",
        "nonbonded_force": "NonbondedForce",
        "periodic_torsion_force": "PeriodicTorsionForce",
        "cmap_torsion_force": "CMAPTorsionForce",
    }

    result: Dict[str, float] = {}

    try:
        ff = ForceField(force_field_file)
    except Exception as e:
        logger.print(f"[ERROR] Failed to load OpenMM force field: {e}")
        return None

    try:
        modeller = Modeller(struct.topology, struct.positions)
        system = ff.createSystem(
            modeller.topology,
            nonbondedMethod=NoCutoff,
            constraints=HBonds,
            rigidWater=True,
            removeCMMotion=True,
        )
    except Exception as e:
        logger.print(f"[ERROR] Failed to create OpenMM system: {e}")
        return None

    try:
        group_map: Dict[int, str] = {}
        for i, force in enumerate(system.getForces()):
            force.setForceGroup(i)
            group_map[i] = type(force).__name__
    except Exception as e:
        logger.print(f"[ERROR] Failed to assign OpenMM force groups: {e}")
        return None

    try:
        integrator = VerletIntegrator(1.0 * unit.femtoseconds)
        platform = Platform.getPlatformByName("CPU")
        context = Context(system, integrator, platform)
        context.setPositions(modeller.positions)
    except Exception as e:
        logger.print(f"[ERROR] Failed to initialize OpenMM context: {e}")
        return None

    if minimize_energy and minimization_iteration:
        try:
            LocalEnergyMinimizer.minimize(
                context,
                tolerance=10.0 * unit.kilojoule_per_mole / unit.nanometer,
                maxIterations=minimization_iteration,
            )
        except Exception as e:
            logger.print(f"[ERROR] Failed during energy minimization: {e}")
            return None

    try:
        state_total = context.getState(getEnergy=True)
        total_energy = state_total.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        result["total_potential_energy"] = float(total_energy)
    except Exception as e:
        logger.print(f"[ERROR] Failed to calculate total potential energy: {e}")
        return None

    for term in terms[1:]:
        force_name = mapping[term]

        matched_group = None
        for gid, name in group_map.items():
            if name == force_name:
                matched_group = gid
                break

        if matched_group is None:
            logger.print(f"[WARNING] OpenMM force {force_name} not found in system. Set to 0.0.")
            result[term] = 0.0
            continue

        try:
            state = context.getState(getEnergy=True, groups={matched_group})
            energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
            result[term] = float(energy)
        except Exception as e:
            logger.print(f"[ERROR] Failed to calculate energy term {term}: {e}")
            return None

    return result



def postprocess_energy_report_to_schema(raw_report: Dict[str, Any]) -> Dict[str, Any]:
    """
    Postprocess the raw EnzyWizard-Energy report into the schema-compliant report.
    """

    raw_energy_terms = raw_report.get("energy_terms", {})

    schema_report: Dict[str, Any] = {
        "report_type": raw_report.get("output_type"),
        "energy_terms": {
            "total_potential_energy": raw_energy_terms.get("total_potential_energy"),
            "harmonic_bond_potential_energy": raw_energy_terms.get("harmonic_bond_force"),
            "harmonic_angle_potential_energy": raw_energy_terms.get("harmonic_angle_force"),
            "custom_bond_potential_energy": raw_energy_terms.get("custom_bond_force"),
            "custom_torsion_potential_energy": raw_energy_terms.get("custom_torsion_force"),
            "custom_nonbonded_potential_energy": raw_energy_terms.get("custom_nonbonded_force"),
            "nonbonded_potential_energy": raw_energy_terms.get("nonbonded_force"),
            "periodic_torsion_potential_energy": raw_energy_terms.get("periodic_torsion_force"),
            "cmap_torsion_potential_energy": raw_energy_terms.get("cmap_torsion_force"),
        },
    }

    return schema_report

def generate_energy_report(energy_terms: Dict[str, float],logger: Logger) -> Dict[str, Any] | None:

    if energy_terms is None or len(energy_terms) == 0:
        logger.print("[ERROR] Energy terms input is empty.")
        return None

    if not isinstance(energy_terms, dict):
        logger.print("[ERROR] Energy terms must be a dictionary.")
        return None

    raw_report = {
        "output_type": "enzywizard_energy",
        "energy_terms": energy_terms,
    }

    return postprocess_energy_report_to_schema(
        raw_report=raw_report
    )
