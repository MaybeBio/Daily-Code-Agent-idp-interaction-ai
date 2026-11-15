# 15-Minute Online Meeting Guide — Prof. Daisuke Tsugama (UTokyo)

## 1. 3-Minute Live Screen-Share Walkthrough (`https://soumyajitmandal06-gif.github.io/unippi-owl/`)
1. **Show Instant 6-Complex Cross-Mechanism Leaderboard (60 sec):**
   - Point out Rank #1 (`1_VIP1_Homodimer_Positive`, UIS = 0.777, dG = -10.35 kcal/mol, Pattern A) and Rank #3 (`2_VIP1_bZIP61_GroupE_Positive`, UIS = 0.732, dG = -10.10 kcal/mol) vs. Rank #5 (`5_bZIP34_Proline_Homodimer_Negative`, UIS = 0.120, dG = +5.00 kcal/mol, 6 Helical Proline Kinks, Pattern C).
2. **Show IDR-Diluted Rescue (`Pattern B`) & CPK21 Phosphorylation Switch (60 sec):**
   - Click Row #4 (`4_VIP1_GRF6_14-3-3_IDR_Diluted_Pos`: ipTM = 0.52, LIS = 0.783, dG = -7.61 kcal/mol).
   - Switch **PTM Gating** dropdown to `CPK21-Phosphorylated (pSer-14-3-3 Lock)`: show dG jump to `-9.81 kcal/mol` and UIS rise to `0.724` while Row #6 (`IDR_Decoy`) stays locked at `+2.10 kcal/mol`.
3. **Show Double-Mutant Epistatic Rescue & Real Crystal Structure `1BRS.cif` (60 sec):**
   - Click Row #5 (`5_bZIP34`), click **Auto-Load Rescue Preset** in the Double-Mutant simulator (`A:PRO5->LEU + B:PRO5->LEU`), and show how removing the Proline kink restores binding from `+5.00` to `-0.10 kcal/mol`.
   - Drag and drop `demo_structures_exp/1BRS.cif` into the browser dropzone to show live 6-chain crystal parsing (`UIS = 0.882, dG = -12.33 kcal/mol`).

## 2. Three Technical Questions to Ask Prof. Tsugama
1. **CHUMA Output Format:** *"What is the exact column/JSON output format of your laboratory's CHUMA algorithm so I can add direct one-click CHUMA file import into UniPPI-OWL?"*
2. **Target Protein Library for EMF3 / VIP1:** *"Which specific candidate protein family or AF3 output batch from your students should we run through `unippi_cli.py` first?"*
3. **Software Manuscript Co-Authorship:** *"Once we test UniPPI-OWL on your laboratory's empirical AF3 models, would you be open to co-authoring a short software application note (e.g., in Bioinformatics or JOSS)?"*
