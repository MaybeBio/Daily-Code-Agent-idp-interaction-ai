"""
The gene locus first drawn in scene 4 and reused, shifted, in later scenes:
compacted chromatin on both sides of an open stretch of DNA carrying one gene
and its promoter.

The storyboard copies the same drawing from scene to scene, so the geometry is
kept here once, in scene-4 slide coordinates; a scene places the locus with
`Locus(offset)`.
"""

import numpy as np

from .molecules import (DNA_BOND, RNA_BOND, TF_RADIUS, Polymer, circle_points, dna_obstacles,
                        draw_dna, draw_rna, draw_rnap, draw_tf, grow_rna)

# The DNA from its lower-left end to its upper end; transcription runs this way.
DNA_TRACE = [
    ("coil", [(1.71, 8.13), (1.93, 8.02), (2.1, 7.74), (2.21, 7.52), (2.43, 7.3), (2.71, 7.19),
              (2.93, 7.22), (3.15, 7.27), (3.28, 7.47), (3.28, 7.69), (3.23, 7.96), (3.2, 8.19),
              (3.26, 8.46), (3.43, 8.68), (3.65, 8.9), (3.81, 9.1), (3.98, 9.27), (4.14, 9.37),
              (4.37, 9.26)]),
    ("linker", [(4.37, 9.26), (4.51, 9.35), (4.59, 9.53), (4.81, 9.62), (5.14, 9.67),
                (5.27, 9.6), (5.38, 9.4), (5.49, 9.18), (5.56, 8.9), (5.6, 8.7), (5.71, 8.54),
                (5.91, 8.42), (6.13, 8.36), (6.3, 8.38)]),
    ("gene", [(6.3, 8.38), (6.69, 8.48), (7.13, 8.59), (7.4, 8.74), (7.51, 8.93), (7.54, 9.16),
              (7.5, 9.37), (7.48, 9.54)]),
    ("linker", [(7.48, 9.54), (7.46, 9.84), (7.51, 10.06), (7.62, 10.25), (7.66, 10.47),
                (7.61, 10.65), (7.51, 10.8)]),
    ("coil", [(7.51, 10.8), (7.4, 10.95), (7.59, 11.17), (7.76, 11.39), (7.87, 11.61),
              (7.92, 11.83), (8.07, 12.11), (8.12, 12.33), (8.09, 12.53), (7.9, 12.69),
              (7.73, 12.83), (7.51, 12.91), (7.24, 13.05), (7.02, 13.13), (6.69, 13.19)]),
]

# The transcript, from the polymerase (where it grows) to its free end; it is
# cut to the length of the gene (see Locus.full_transcript).
FULL_TRANSCRIPT_TRACE = [
    (7.46, 8.85), (7.43, 8.57), (7.26, 8.37), (7.19, 8.19), (6.93, 8.06), (6.65, 8.06),
    (6.4, 8.02), (6.44, 7.74), (6.33, 7.6), (6.08, 7.67), (5.8, 7.85), (5.58, 7.83),
    (5.54, 7.67), (5.75, 7.52), (6.08, 7.43), (6.38, 7.27), (6.57, 7.09), (6.72, 6.95),
    (6.65, 6.69), (6.41, 6.55), (6.08, 6.53), (5.8, 6.59), (5.58, 6.46)]

POLYMERASE_AT_END = (7.48, 9.01)  # where the sketch has it, near the end of the gene
POLYMERASE_START_OFFSET = 2.0     # nodes past the start of the gene, at initiation
PROMOTER_OFFSET = 4.0             # nodes before the start of the gene, where a TF binds
PROMOTER_SPAN = (11, 2)           # the promoter runs from 11 to 2 nodes before the gene
POLYMERASE_RADIUS = 0.4
BOUND_TF_OFFSET = 0.22            # distance from the DNA of a TF bound beside it

# Transcripts leave the back of the polymerase and trail behind it: they head
# backwards along the gene, tilted away from the DNA by TRAIL_TILT (0 = straight
# back, 1 = 45 degrees), wiggling little and pulled firmly back to that heading.
TRAIL_TILT = 0.7
TRAIL_WIGGLE = 0.3
TRAIL_PULL = 0.3


def shifted(points, offset):
    return [(x + offset[0], y + offset[1]) for x, y in points]


class Locus:
    # `trace` replaces DNA_TRACE, e.g. by DNA_TRACE continued further along the
    # chain (scene 20); it must keep DNA_TRACE up to the end of the gene.
    def __init__(self, offset=(0.0, 0.0), trace=DNA_TRACE):
        self.offset = np.asarray(offset, float)
        self.dna = Polymer.from_trace([(kind, shifted(pts, offset)) for kind, pts in trace],
                                      DNA_BOND)
        gene_start = self.dna.span("gene")[0]
        # Node indices: polymerase at initiation and termination, TF binding site.
        self.start = gene_start + POLYMERASE_START_OFFSET
        self.end = self.dna.nearest(self.offset + POLYMERASE_AT_END)
        self.promoter = gene_start - PROMOTER_OFFSET
        # Drawn as an orange stretch of the DNA just before the gene.
        first, last = (int(gene_start - n) for n in PROMOTER_SPAN)
        self.dna.kinds[first:last] = "promoter"
        gene_first, gene_last = self.dna.span("gene")
        self.gene_length = (gene_last - gene_first) * DNA_BOND

    def progress(self, fraction):
        """Node index of a polymerase that has transcribed `fraction` of the gene."""
        return self.start + fraction * (self.end - self.start)

    def full_transcript(self):
        """The finished mRNA: no longer than the gene it copies."""
        return Polymer.from_trace([("rna", shifted(FULL_TRANSCRIPT_TRACE, self.offset))],
                                  RNA_BOND).cut(self.gene_length)

    def transcribe(self, fractions, seed, blocked=()):
        """(node index, transcript or None) for polymerases on the gene, given the
        fraction of the gene each has transcribed (0 = initiating, 1 = terminating).

        A transcript's length is proportional to how far its polymerase has gone;
        the finished one is the traced scene-4 mRNA, cut to the length of the gene;
        shorter ones are grown at random, trailing behind their polymerase, clear of
        the DNA, of the other polymerases, of the `blocked` points (e.g. bound
        proteins) and of each other.
        """
        rng = np.random.default_rng(seed)
        full = self.full_transcript()
        positions = [self.progress(f) for f in fractions]
        blocked = [self.dna.curve()[0], *blocked]
        states = []
        for fraction, u in sorted(zip(fractions, positions), reverse=True):
            n_beads = round(fraction * len(full.nodes))
            if fraction == 1:
                rna = full
            elif n_beads >= 2:
                p, t = self.dna.at(u)
                away = np.array([t[1], -t[0]])  # right of the direction of transcription
                others = [circle_points(self.dna.at(v)[0], POLYMERASE_RADIUS)
                          for v in positions if v != u]
                rna = grow_rna(p, -t + TRAIL_TILT * away, n_beads, rng,
                               avoid=np.vstack(blocked + others), wiggle=TRAIL_WIGGLE,
                               pull=TRAIL_PULL)
            else:
                rna = None
            if rna is not None:
                blocked.append(rna.nodes)
            states.append((u, rna))
        return states

    def bound_tfs(self, sites):
        """Centres of TFs bound beside the DNA. Each site is (nodes from the start
        of the gene, negative = upstream; side, +1 = left of the direction of
        transcription, -1 = right)."""
        gene_start = self.dna.span("gene")[0]
        centers = []
        for nodes, side in sites:
            p, t = self.dna.at(gene_start + nodes)
            centers.append(p + side * BOUND_TF_OFFSET * np.array([-t[1], t[0]]))
        return centers

    def obstacles(self, states=(), tfs=()):
        """What scattered molecules must keep clear of, as (points, radius) pairs:
        the DNA, polymerases and transcripts (from transcribe()), bound TFs."""
        avoid = dna_obstacles(self.dna)
        for u, rna in states:
            avoid.append((self.dna.at(u)[0], POLYMERASE_RADIUS))
            if rna is not None:
                avoid.append((rna.nodes, 0.06))
        avoid += [(c, TF_RADIUS) for c in tfs]
        return avoid

    def draw(self, ax, states=(), tfs=(), tf_angle=20):
        """The DNA, polymerases with their transcripts, and bound TFs."""
        draw_dna(ax, self.dna, z=1)
        for u, rna in states:
            if rna is not None:
                draw_rna(ax, rna, z=3)
            draw_rnap(ax, self.dna, u, z=4)
        for i, c in enumerate(tfs):
            draw_tf(ax, c, angle=tf_angle + 70 * i, z=5)
