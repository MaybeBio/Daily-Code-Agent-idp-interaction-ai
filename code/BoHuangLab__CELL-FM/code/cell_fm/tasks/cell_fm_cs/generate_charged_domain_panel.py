import numpy as np
import random
import pandas as pd
import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.extend([".", ".."])

from cell_fm.tasks.cell_fm_cs.generate_mutation_scd import make_blocky_pattern, make_alternating_pattern


CHARGED_SET = set("KRDE")
POS_SET = set("KR")
NEG_SET = set("DE")
SPACER_SET = set("GSQNT")

NUM_POS = 15
NUM_NEG = 15
NUM_SPACER = 36
NUM_GENERATED = 100
SEED = 42

save_path = "cell_fm/tasks/cell_fm_cs/data/random_KRDE/"

# randomly generate a sequence with specified number of positive and negative charges
def generate_sequence(num_pos=NUM_POS, num_neg=NUM_NEG, spacer_len=NUM_SPACER):
    positives = random.choices(list(POS_SET), k=num_pos)
    negatives = random.choices(list(NEG_SET), k=num_neg)
    spacers = random.choices(list(SPACER_SET), k=spacer_len)
    all_residues = positives + negatives + spacers
    random.shuffle(all_residues)
    return "".join(all_residues)

seq_list = []
alt_list = []
blocky_list = []

random.seed(SEED)
while len(seq_list) < NUM_GENERATED:
    seq = generate_sequence()

    alt_seq = make_alternating_pattern(seq)
    blocky_seq = make_blocky_pattern(seq)

    if (seq not in seq_list) and (alt_seq not in alt_list) and (blocky_seq not in blocky_list):
        seq_list.append(seq)
        alt_list.append(alt_seq)
        blocky_list.append(blocky_seq)
    else:
        continue

# save to CSV
os.makedirs(save_path, exist_ok=True)

df = pd.DataFrame({"protein_seq": seq_list})
output_path = os.path.join(save_path, "panel.csv")
df.to_csv(output_path, index=True)
print(f"Saved {len(seq_list)} sequences to {output_path}")

