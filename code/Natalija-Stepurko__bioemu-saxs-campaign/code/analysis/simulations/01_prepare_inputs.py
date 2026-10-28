"""The input table: the 41 clean disordered profiles with a BioEmu-1 ensemble, their construct sequences, the
conditions from each SASBDB record by the design's rule (§4), the run settings and seeds.
Writes results/simulations/inputs.csv. Runs in the main environment or the simulation environment."""

import pandas as pd
from common import OUT, RES
from simlib import prepare_inputs

inp = prepare_inputs(
    pd.read_csv(RES / "resolvability.csv"),
    pd.read_csv(RES / "entries.csv"),
    pd.read_csv(RES / "coverage" / "sasbdb_metadata.csv"),
)
inp.to_csv(OUT / "inputs.csv", index=False)
print(
    f"{len(inp)} proteins, {inp['length'].min()}-{inp['length'].max()} residues; "
    f"ionic strength from the record for {(inp['ionic_source'] == 'parsed').sum()}, "
    f"denaturant in {inp['denaturant'].sum()}; {inp['ns'].sum():.0f} ns in total"
)
