# Skills-BIOKT

## gmx-prep
- [ ] Improve gmx-prep skill with force-field-specific mdp templates (AMBER ff99SBws-STQ, ff19SB, a99SBdisp, CHARMM36m)
- [ ] Replace H-REMD in `references/hremd.md` with REST2. Current template is broken: no atom has a distinct B state, so `vdw-lambdas` perturbs nothing and every replica runs the same Hamiltonian (verified on GROMACS 2026.3, 2026-10-03). Base the rewrite on the working protocol in `Colabs/Cyclic/Peptides-mal/Linear/prepare_runs.py` (`plumed partial_tempering`, geometric T ladder, `-hrex -replex 100`). Decide: put REST2 in `plumed-remd` (needs PLUMED-patched gmx) with a pointer from gmx-prep? 8 replicas / Tmax 400 K as defaults or per-system examples?
- [ ] Confirm water model for ff99SBws-STQ: skill uses TIP4P/2005, but the FF's `watermodels.dat` marks `tip4p2005s` (protein–water ×1.1) as recommended.
