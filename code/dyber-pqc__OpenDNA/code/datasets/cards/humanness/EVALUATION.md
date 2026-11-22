# Humanness table: does it separate the cohorts it should?

Scored with a table of 141 donors and 1096088 peptides against the four antibody cohorts, whose donors are disjoint from the table's. Each number is a **per-sequence identity**: the fraction of a sequence's scorable 9-mers occurring in at least that share of the table's donors. Nothing is tuned and no threshold was chosen after seeing the numbers.

What the right answer is, fixed before measuring:

- `human` -- high: donors the table has never seen
- `mouse` -- lower: mouse germline, scored against human donors
- `diagnostic` -- high: mice with human immunoglobulin loci, so the molecule is human
- `out_of_distribution` -- lower, and not a ranking between the species

## loose (>=1% of donors)

### By cohort and organism

| Cohort | group | sequences | p5 | p25 | median | p75 | p95 |
|---|---|---:|---:|---:|---:|---:|---:|
| **human** | all | 94923 | 53.8% | 70.3% | **80.8%** | 88.6% | 93.8% |
|  | human | 94923 | 53.8% | 70.3% | **80.8%** | 88.6% | 93.8% |
| **mouse** | all | 81620 | 29.6% | 45.8% | **58.9%** | 70.0% | 89.1% |
|  | mouse_BALB/c | 30629 | 25.2% | 39.7% | **53.1%** | 63.6% | 88.0% |
|  | mouse_C57BL/6 | 45991 | 33.3% | 48.7% | **60.9%** | 71.7% | 88.2% |
|  | mouse_RAG2-GFP/129Sve | 5000 | 46.9% | 64.7% | **73.0%** | 83.3% | 92.7% |
| **diagnostic** | all | 30000 | 82.9% | 90.3% | **93.0%** | 94.8% | 98.2% |
|  | HIS-Mouse | 10000 | 77.8% | 85.6% | **89.7%** | 92.9% | 96.2% |
|  | Kymouse | 20000 | 89.8% | 92.3% | **93.9%** | 95.5% | 98.3% |
| **out_of_distribution** | all | 40000 | 25.9% | 36.1% | **51.3%** | 72.3% | 83.2% |
|  | bactrian_camel | 15000 | 21.8% | 30.8% | **37.2%** | 43.8% | 54.3% |
|  | rabbit | 5000 | 25.2% | 30.7% | **34.5%** | 39.1% | 51.3% |
|  | rat | 5000 | 39.8% | 53.2% | **64.9%** | 71.8% | 80.2% |
|  | rhesus | 15000 | 55.5% | 67.8% | **74.1%** | 79.6% | 86.6% |

### By isotype

| Cohort | group | sequences | p5 | p25 | median | p75 | p95 |
|---|---|---:|---:|---:|---:|---:|---:|
| **human** | all | 94923 | 53.8% | 70.3% | **80.8%** | 88.6% | 93.8% |
|  | Bulk | 10000 | 49.1% | 68.6% | **83.5%** | 90.6% | 94.6% |
|  | IGHA | 19923 | 48.7% | 65.5% | **73.4%** | 80.3% | 88.6% |
|  | IGHG | 40000 | 52.8% | 67.4% | **76.5%** | 84.4% | 91.3% |
|  | IGHM | 25000 | 77.3% | 86.2% | **89.4%** | 92.0% | 95.1% |
| **mouse** | all | 81620 | 29.6% | 45.8% | **58.9%** | 70.0% | 89.1% |
|  | Bulk | 28793 | 23.8% | 35.6% | **50.5%** | 65.7% | 90.4% |
|  | IGHA | 5000 | 38.0% | 49.1% | **59.5%** | 71.1% | 86.2% |
|  | IGHG | 25000 | 38.7% | 50.9% | **60.0%** | 68.5% | 83.3% |
|  | IGHM | 22827 | 38.0% | 52.3% | **64.4%** | 75.2% | 90.7% |
| **diagnostic** | all | 30000 | 82.9% | 90.3% | **93.0%** | 94.8% | 98.2% |
|  | IGHA | 5000 | 78.6% | 86.0% | **89.7%** | 92.7% | 95.3% |
|  | IGHG | 5000 | 76.9% | 85.5% | **89.6%** | 92.9% | 97.3% |
|  | IGHM | 20000 | 89.8% | 92.3% | **93.9%** | 95.5% | 98.3% |
| **out_of_distribution** | all | 40000 | 25.9% | 36.1% | **51.3%** | 72.3% | 83.2% |
|  | Bulk | 20000 | 22.5% | 30.8% | **36.4%** | 42.7% | 53.6% |
|  | IGHA | 5000 | 39.8% | 53.2% | **64.9%** | 71.8% | 80.2% |
|  | IGHG | 5000 | 61.1% | 69.4% | **74.1%** | 78.8% | 84.9% |
|  | IGHM | 10000 | 52.9% | 66.7% | **74.1%** | 80.2% | 87.4% |


## medium (>=10% of donors)

### By cohort and organism

| Cohort | group | sequences | p5 | p25 | median | p75 | p95 |
|---|---|---:|---:|---:|---:|---:|---:|
| **human** | all | 94923 | 38.3% | 56.6% | **70.7%** | 83.5% | 89.9% |
|  | human | 94923 | 38.3% | 56.6% | **70.7%** | 83.5% | 89.9% |
| **mouse** | all | 81620 | 11.4% | 19.2% | **25.9%** | 36.2% | 68.4% |
|  | mouse_BALB/c | 30629 | 9.6% | 17.3% | **23.6%** | 34.3% | 81.6% |
|  | mouse_C57BL/6 | 45991 | 12.6% | 19.8% | **26.6%** | 35.8% | 62.2% |
|  | mouse_RAG2-GFP/129Sve | 5000 | 20.7% | 27.1% | **34.5%** | 40.8% | 64.3% |
| **diagnostic** | all | 30000 | 76.1% | 87.2% | **90.0%** | 91.8% | 94.1% |
|  | HIS-Mouse | 10000 | 69.0% | 79.3% | **84.4%** | 88.9% | 92.9% |
|  | Kymouse | 20000 | 87.3% | 89.4% | **90.9%** | 92.2% | 94.4% |
| **out_of_distribution** | all | 40000 | 14.8% | 23.5% | **36.4%** | 58.8% | 72.8% |
|  | bactrian_camel | 15000 | 11.8% | 18.3% | **23.7%** | 29.7% | 40.7% |
|  | rabbit | 5000 | 15.5% | 20.7% | **24.3%** | 28.6% | 37.9% |
|  | rat | 5000 | 22.7% | 35.1% | **45.5%** | 52.7% | 61.1% |
|  | rhesus | 15000 | 41.1% | 55.0% | **62.5%** | 69.1% | 77.0% |

### By isotype

| Cohort | group | sequences | p5 | p25 | median | p75 | p95 |
|---|---|---:|---:|---:|---:|---:|---:|
| **human** | all | 94923 | 38.3% | 56.6% | **70.7%** | 83.5% | 89.9% |
|  | Bulk | 10000 | 32.8% | 54.1% | **74.1%** | 86.7% | 91.5% |
|  | IGHA | 19923 | 33.0% | 50.6% | **60.7%** | 69.9% | 81.2% |
|  | IGHG | 40000 | 37.6% | 53.1% | **64.9%** | 75.8% | 87.0% |
|  | IGHM | 25000 | 65.8% | 80.2% | **85.3%** | 88.0% | 91.7% |
| **mouse** | all | 81620 | 11.4% | 19.2% | **25.9%** | 36.2% | 68.4% |
|  | Bulk | 28793 | 8.7% | 15.2% | **22.1%** | 39.0% | 86.0% |
|  | IGHA | 5000 | 14.7% | 21.8% | **27.6%** | 33.3% | 64.0% |
|  | IGHG | 25000 | 14.8% | 20.2% | **25.0%** | 32.1% | 51.9% |
|  | IGHM | 22827 | 14.4% | 23.0% | **30.2%** | 39.6% | 58.4% |
| **diagnostic** | all | 30000 | 76.1% | 87.2% | **90.0%** | 91.8% | 94.1% |
|  | IGHA | 5000 | 70.9% | 80.3% | **84.6%** | 88.9% | 92.7% |
|  | IGHG | 5000 | 67.5% | 77.9% | **83.9%** | 89.2% | 93.6% |
|  | IGHM | 20000 | 87.3% | 89.4% | **90.9%** | 92.2% | 94.4% |
| **out_of_distribution** | all | 40000 | 14.8% | 23.5% | **36.4%** | 58.8% | 72.8% |
|  | Bulk | 20000 | 12.4% | 19.0% | **23.9%** | 29.4% | 40.0% |
|  | IGHA | 5000 | 22.7% | 35.1% | **45.5%** | 52.7% | 61.1% |
|  | IGHG | 5000 | 46.6% | 56.8% | **62.7%** | 68.5% | 75.5% |
|  | IGHM | 10000 | 38.5% | 53.7% | **62.5%** | 69.6% | 77.8% |


## strict (>=50% of donors)

### By cohort and organism

| Cohort | group | sequences | p5 | p25 | median | p75 | p95 |
|---|---|---:|---:|---:|---:|---:|---:|
| **human** | all | 94923 | 20.7% | 37.7% | **53.0%** | 70.6% | 85.0% |
|  | human | 94923 | 20.7% | 37.7% | **53.0%** | 70.6% | 85.0% |
| **mouse** | all | 81620 | 1.9% | 6.4% | **10.0%** | 15.6% | 41.0% |
|  | mouse_BALB/c | 30629 | 1.8% | 6.7% | **9.9%** | 15.9% | 63.8% |
|  | mouse_C57BL/6 | 45991 | 1.9% | 6.3% | **9.8%** | 15.1% | 36.1% |
|  | mouse_RAG2-GFP/129Sve | 5000 | 5.8% | 9.3% | **12.5%** | 18.0% | 38.1% |
| **diagnostic** | all | 30000 | 58.1% | 73.6% | **77.7%** | 80.7% | 88.5% |
|  | HIS-Mouse | 10000 | 50.8% | 63.2% | **70.5%** | 76.5% | 84.0% |
|  | Kymouse | 20000 | 73.0% | 77.0% | **78.8%** | 83.0% | 89.2% |
| **out_of_distribution** | all | 40000 | 5.7% | 10.2% | **17.9%** | 40.4% | 56.0% |
|  | bactrian_camel | 15000 | 3.7% | 8.0% | **11.5%** | 15.7% | 24.1% |
|  | rabbit | 5000 | 6.2% | 8.0% | **8.9%** | 10.7% | 15.1% |
|  | rat | 5000 | 8.3% | 14.5% | **18.9%** | 23.9% | 32.4% |
|  | rhesus | 15000 | 23.0% | 37.1% | **45.2%** | 52.2% | 60.0% |

### By isotype

| Cohort | group | sequences | p5 | p25 | median | p75 | p95 |
|---|---|---:|---:|---:|---:|---:|---:|
| **human** | all | 94923 | 20.7% | 37.7% | **53.0%** | 70.6% | 85.0% |
|  | Bulk | 10000 | 13.6% | 31.9% | **52.8%** | 73.5% | 84.4% |
|  | IGHA | 19923 | 16.1% | 31.0% | **41.1%** | 51.7% | 65.8% |
|  | IGHG | 40000 | 21.6% | 35.5% | **47.5%** | 59.6% | 76.8% |
|  | IGHM | 25000 | 45.7% | 65.8% | **75.8%** | 83.0% | 87.8% |
| **mouse** | all | 81620 | 1.9% | 6.4% | **10.0%** | 15.6% | 41.0% |
|  | Bulk | 28793 | 1.0% | 5.8% | **8.9%** | 18.1% | 71.4% |
|  | IGHA | 5000 | 2.8% | 6.9% | **9.8%** | 13.2% | 36.1% |
|  | IGHG | 25000 | 3.7% | 6.3% | **9.2%** | 13.0% | 28.7% |
|  | IGHM | 22827 | 1.8% | 8.2% | **11.8%** | 17.5% | 33.3% |
| **diagnostic** | all | 30000 | 58.1% | 73.6% | **77.7%** | 80.7% | 88.5% |
|  | IGHA | 5000 | 54.7% | 65.8% | **71.8%** | 76.9% | 84.4% |
|  | IGHG | 5000 | 48.7% | 60.6% | **69.0%** | 75.5% | 83.6% |
|  | IGHM | 20000 | 73.0% | 77.0% | **78.8%** | 83.0% | 89.2% |
| **out_of_distribution** | all | 40000 | 5.7% | 10.2% | **17.9%** | 40.4% | 56.0% |
|  | Bulk | 20000 | 4.3% | 8.0% | **10.5%** | 14.5% | 22.5% |
|  | IGHA | 5000 | 8.3% | 14.5% | **18.9%** | 23.9% | 32.4% |
|  | IGHG | 5000 | 28.6% | 39.8% | **46.4%** | 52.7% | 60.4% |
|  | IGHM | 10000 | 20.7% | 35.5% | **44.3%** | 51.8% | 59.6% |


## The confound, and what is left after it

The diagnostic cohort is mostly IGHM, the human cohort spans isotypes, and a class-switched sequence has been mutated away from germline -- which is where the peptides shared by many donors are. So the diagnostic cohort scoring above the human one could be a difference in somatic mutation rather than in the species of the loci. Matched isotype by isotype at the strict (>=50% of donors) threshold, median identity:

| Isotype | human | diagnostic | difference |
|---|---:|---:|---:|
| IGHA | 41.1% | 71.8% | +30.7 pp |
| IGHG | 47.5% | 69.0% | +21.5 pp |
| IGHM | 75.8% | 78.8% | +3.0 pp |

The diagnostic cohort scores at least as high as the human cohort on 3 of 3 shared isotypes, so the result is not the isotype difference. These are mice, and the measure calls their antibodies human: it tracks the immunoglobulin rather than the animal, which is the behaviour a humanness measure needs and the reason this cohort carries no species label. A classifier trained on donor species would have been trained to call them mouse.

What it does not show is that the measure beats a germline lookup here. These sequences carry human V genes, so a V-gene table also answers "human". The peptide measure earns its place on resolution instead: it scores every window, including the junction and CDR3, where there is no germline gene to look up and where a humanised antibody actually departs from anything a human has made.

That the diagnostic cohort scores *higher* is not a claim that transgenic mice make better human antibodies. Their repertoires are less affinity-matured, so they sit closer to germline, and germline is where prevalence is highest.

The same effect sets a limit on reading any single score. Within the human cohort alone, median identity runs from 41.1% (IGHA) to 75.8% (IGHM) -- a spread of 35 percentage points between isotypes of the same species. An absolute identity therefore cannot be read as "how human" without knowing how mutated the sequence is, which is another reason the per-window output is the useful one.

## Rows that could not be scored

| Cohort | rows | no scorable window |
|---|---:|---:|
| human | 94923 | 0 |
| mouse | 81620 | 0 |
| diagnostic | 30000 | 0 |
| out_of_distribution | 40000 | 0 |

A sequence whose every 9-mer holds a non-standard residue has no measurement, and is left out of both numerator and denominator rather than scored zero.

## What this does and does not establish

It establishes that the table behaves as the method claims on donors it has never seen, and that it answers the humanisation question a species classifier cannot: the diagnostic cohort is mice, and a measure of humanness should call their antibodies human.

It establishes nothing about immunogenicity. These are prevalences in sampled repertoires, and whether an antibody provokes a response in a person depends on MHC presentation, dose, route and the recipient. It also establishes nothing comparable with published OASis figures: this table's donor units are shallow, so every prevalence here is biased down. The table's own card gives the size of that bias.
