"""Force-field layer"""

# rq-b45bedcc rq-2691f166 rq-d843e829 rq-ad30757b
from ._amber import read_amber_parameters
from ._charges import charge_template, load_charge_file
from ._core import parametrize
from ._dft import DftChargeError, DftCharges, export_dft_inputs, export_dft_molecule, load_dft_charges
from ._errors import (
    ChargeFileError,
    ForceFieldError,
    ForceFieldUnavailableError,
    MissingParameterError,
    ParameterFileError,
    ParametrizationError,
    TypingError,
    UnknownForceFieldError,
)
from ._gaff2 import Gaff2Backend
from ._gromacs import write_gromacs
from ._ions import IonGroup, ionic_groups
from ._model import AtomType, ChargeReport, ForceFieldBackend, FragmentParameters, ParameterSet, ResidueBlock, Topology
from ._oplsaa import GromacsForceField, OplsaaBackend, TypingRules, load_typing_rules, read_gromacs_forcefield
from ._params import ionic_parameter_template, load_parameter_file
from ._registry import available_backends, get_backend, register_backend

register_backend(Gaff2Backend())
register_backend(OplsaaBackend())

__all__ = [
    "AtomType",
    "ChargeFileError",
    "ChargeReport",
    "DftChargeError",
    "DftCharges",
    "ForceFieldBackend",
    "ForceFieldError",
    "ForceFieldUnavailableError",
    "FragmentParameters",
    "Gaff2Backend",
    "GromacsForceField",
    "IonGroup",
    "MissingParameterError",
    "OplsaaBackend",
    "ParameterFileError",
    "ParameterSet",
    "ParametrizationError",
    "ResidueBlock",
    "Topology",
    "TypingError",
    "TypingRules",
    "UnknownForceFieldError",
    "available_backends",
    "charge_template",
    "export_dft_inputs",
    "export_dft_molecule",
    "get_backend",
    "ionic_groups",
    "ionic_parameter_template",
    "load_charge_file",
    "load_dft_charges",
    "load_parameter_file",
    "load_typing_rules",
    "parametrize",
    "read_amber_parameters",
    "read_gromacs_forcefield",
    "register_backend",
    "write_gromacs",
]
