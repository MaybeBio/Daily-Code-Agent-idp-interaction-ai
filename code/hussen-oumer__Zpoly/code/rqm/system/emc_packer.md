# Feature: EMC-Style Packer (`PACKER = emc`) <!-- rq-e6bbf42d -->

An Enhanced-Monte-Carlo-style packer (after EMC by P. J. in 't Veld and the Theodorou–Suter
amorphous-cell method) for large polymer systems: it starts from the chain built by the polymer
builder (the same molecule and coordinates as the `JOB = build` PDB), makes `CHAINS` copies with
different conformations, inserts them with a configurational bias, and removes the remaining
overlaps with Monte Carlo moves. Water and salt are added by GROMACS (`gmx solvate`, `gmx genion`),
as with the built-in packer.

## Method <!-- rq-f4bb437c -->

1. **Chain once, copies by pivots.** The chain is built once. Each copy gets its own conformation from
   `pivots_per_copy` random **pivot moves**: a rotation by a random angle about a rotatable backbone
   bond (single, not in a ring) of the atoms on one side of it. Bond lengths and angles stay exact; a
   pivot that makes the chain overlap itself is rejected.
2. **Cell lists.** Heavy atoms are binned into periodic cells of edge ≥ `min_distance`; an overlap test
   looks only at the 27 neighbouring cells. An **overlap** is a pair of heavy atoms (H excluded, as in
   the built-in packer) closer than `min_distance` (2 Å) that belong to different molecules, or to the
   same chain and are more than 3 bonds apart.
3. **Configurational-bias insertion.** For each copy, `insertion_trials` random orientations and
   positions (inside the slab region for `slab`) are tried; the trial with the fewest overlaps is kept
   (ties: the first), even if it still overlaps.
4. **MC relaxation.** Sweeps over the molecules that still overlap; each takes a random move from
   the move table (`translate`, `rotate`, `pivot`; `local_pivot`: a pivot about one of the 3 bonds
   nearest to an overlapping atom; `reinsert`: the best of `insertion_trials / 4` new random spots;
   weighted by `move_weights`; step sizes drawn between 10 % and 100 % of the maximum) and keeps it if its
   overlap count does not increase (and every atom stays inside the slab region). It stops when there
   are no overlaps; after `max_sweeps` it raises `SystemBuildError` with the remaining overlaps and the
   advice to enlarge the box or the slab.
5. **Inflation.** Insertion and a first relaxation use a small overlap radius (`start_radius`,
   0.8 Å); the radius then grows to `min_distance` in `inflation_steps` steps (6), relaxing after
   each (the Lubachevsky–Stillinger idea), so each step only has a few new overlaps to remove.
6. Ions (counterions) are packed after the chains in the same way (rigid, single atom).
7. **Range.** With 12-mers it packs to about 0.45 g/cm³ in seconds. Higher densities are reached
   with `DENSITY` (`rqm/system/bulk.md`, `rqm/system/slab.md`): the packer works at 0.35 g/cm³ and
   GROMACS shrinks the box (`deform`) to the requested density.

Random numbers come from `seed`, so a result is reproducible.

## Settings <!-- rq-aeb72655 -->

`EmcSettings` (frozen dataclass, all in one place): `min_distance` 2.0 Å, `cell` 3.0 Å,
`insertion_trials` 200, `pivots_per_copy` 20, `max_sweeps` 500 (per inflation step),
`inflation_steps` 6, `start_radius` 0.8 Å, `max_shift` 2.0 Å,
`max_angle` 30°, `move_weights` {translate: 1, rotate: 1, pivot: 1, local_pivot: 4, reinsert: 0.5}.

## Extension Points <!-- rq-481bb439 -->

- `MOVES`: name → function `(molecule, rng, settings, world) -> new coordinates | None` (`world.box`,
  `world.overlap(xyz)`); a new move is one function and one entry (and a weight).
- `EmcSettings` fields are the tuning knobs; `emc_pack(..., settings=EmcSettings(...))`.
- `PACKERS` in `system/_pack.py`: name → packer function, so another packer is one entry.

## Feature API <!-- rq-dff41c8f -->

- `emc_pack(pieces, box, seed=0, settings=None) -> list[np.ndarray]` (same contract as the other <!-- rq-5f40b4b3 -->
  packers: one coordinate array per piece, Å). `Piece.pivots` lists `(a, b, moving atom indices)` for
  the pivot bonds; `chain_piece` fills it.
- `bulk_system` / `slab_system(..., packer="emc")`; input key `PACKER = emc`. With `emc`, the chain is <!-- rq-8182c0fe -->
  built once and copied.

## Gherkin Scenarios <!-- rq-b0b88b25 -->

```gherkin
Feature: EMC-style packer

  @rq-fb226d8e
  Scenario: Pivot moves keep the covalent geometry
    Given a DMAPS 5-mer
    When 50 pivot moves are applied
    Then every bond length and bond angle is unchanged within 1e-6

  @rq-53277b50
  Scenario: Dense slab without overlaps
    Given the fake backend
    When a slab of 4 DMAPS 6-mers with KPS is packed with packer="emc" in a 3.4 x 3.4 x 10.2 nm box with a 4 nm slab
    Then no two heavy atoms of different molecules are closer than 2 Å and every chain atom is inside the slab region
    And run.sh solvates with gmx solvate as for the built-in packer

  @rq-758b974f
  Scenario: Faster than regrowing
    When a bulk system with packer="emc" is built
    Then the chain is built once (one build_polymer call for the chains)

  @rq-41f479b2
  Scenario: Reproducible
    When emc_pack is called twice with the same seed
    Then the coordinates are identical

  @rq-df4bcf3b
  Scenario: Too dense
    When the molecules cannot fit (box smaller than the polymer volume)
    Then SystemBuildError is raised naming the remaining overlaps

  @rq-fdfb3cb2
  Scenario: New move
    When a function is added to MOVES with a weight in EmcSettings
    Then emc_pack uses it
```
