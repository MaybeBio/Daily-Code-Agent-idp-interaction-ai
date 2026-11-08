"""Linear chains of repeat units with initiator end groups (2D structure only)"""

from __future__ import annotations

import re
from collections.abc import Sequence

from rdkit import Chem, rdBase

from . import _chem
from ._errors import ChainError, EndGroupError
from ._library import Library, default_library
from ._model import Chain, Ends, Initiator, Zwitterion

_LABEL = re.compile(r"\[\*:([12])\]")
_HYDROGEN_CAP = "[H][*:1]"
_METHYL_CAP = "C[*:1]"
CAPS = ("CH3", "H")


def _relabel(smiles: str, inner: int, outer: int | None = None) -> str:
    """Rename [*:1] -> [*:inner] and [*:2] -> [*:outer] in one pass"""
    return _LABEL.sub(lambda m: f"[*:{inner if m.group(1) == '1' else outer}]", smiles)


def _resolve_ends(ends: Ends | str | None, initiator: Initiator | None) -> Ends:
    choice = None if ends is None else Ends(ends)
    if initiator is None:
        if choice not in (None, Ends.NONE):
            raise EndGroupError(choice.value, None, "end caps need an initiator")
        return Ends.NONE
    has_head, has_tail = initiator.head_smiles is not None, initiator.tail_smiles is not None
    if choice is None:
        if has_head and has_tail:
            return Ends.BOTH
        return Ends.HEAD if has_head else Ends.TAIL
    needs_head = choice in (Ends.BOTH, Ends.HEAD)
    needs_tail = choice in (Ends.BOTH, Ends.TAIL)
    if needs_head and not has_head:
        raise EndGroupError(choice.value, initiator.name, "the initiator has no head fragment")
    if needs_tail and not has_tail:
        raise EndGroupError(choice.value, initiator.name, "the initiator has no tail fragment")
    return choice


# rq-a2534bc5 rq-fbb9e221
def build_chain(
    units: Sequence[Zwitterion | str],
    *,
    initiator: Initiator | str | None = None,
    ends: Ends | str | None = None,
    library: Library | None = None,
    cap: str = "CH3",
) -> Chain:
    """Link repeat units head-to-tail and cap the chain ends with initiator fragments"""
    if cap not in CAPS:
        raise ValueError(f"cap must be 'CH3' or 'H', got {cap!r}")
    lib = library if library is not None else default_library()
    resolved = [lib.get(u) if isinstance(u, str) else u for u in units]
    if not resolved:
        raise ChainError("A chain needs at least one unit.")
    for unit in resolved:
        if not unit.polymerizable:
            raise ChainError(f"'{unit.name}' is not polymerizable and cannot be a chain unit.")
    init = lib.initiator(initiator) if isinstance(initiator, str) else initiator
    choice = _resolve_ends(ends, init)

    head = init.head_smiles if init is not None and choice in (Ends.BOTH, Ends.HEAD) else None
    tail = init.tail_smiles if init is not None and choice in (Ends.BOTH, Ends.TAIL) else None
    n = len(resolved)
    # labels: head <-1-> unit 0 <-2-> ... <-(n+1)-> tail
    pieces = [_relabel(head or _HYDROGEN_CAP, 1)]
    pieces += [_relabel(u.repeat_unit_smiles(), i + 1, i + 2) for i, u in enumerate(resolved)]
    # rq-8994396e: free ends: head H, tail cap
    pieces.append(_relabel(tail or (_METHYL_CAP if cap == "CH3" else _HYDROGEN_CAP), n + 1))

    params = Chem.SmilesParserParams()
    params.removeHs = False
    with rdBase.BlockLogs():
        mol = Chem.MolFromSmiles(".".join(pieces), params)
    if mol is None:
        raise ChainError("Could not assemble the chain SMILES.")
    mol = Chem.molzip(mol)
    if any(a.GetAtomicNum() == 0 for a in mol.GetAtoms()):
        raise ChainError("Unmatched attachment points remain in the chain.")
    mol = Chem.RemoveHs(mol)
    Chem.SanitizeMol(mol)
    formula, weight = _chem.formula_and_weight(mol)
    return Chain(
        smiles=Chem.MolToSmiles(mol),
        units=tuple(u.name for u in resolved),
        initiator=init.name if init is not None else None,
        ends=choice,
        head_group=head,
        tail_group=tail,
        net_charge=_chem.net_charge(mol),
        counterion=_counterion(_chem.net_charge(mol), init if head or tail else None, resolved),
        formula=formula,
        molecular_weight=weight,
        cap=cap,
    )


# rq-ddcf8658
def _counterion(net: int, init, units) -> str | None:
    """The ion that neutralizes the chain"""
    if net == 0:
        return None
    candidates = ([init.counterion] if init is not None else []) + [u.counterion for u in units]
    for ion in candidates:
        if ion and ion.endswith("-") == (net > 0):
            return ion
    return init.counterion if init is not None else None
