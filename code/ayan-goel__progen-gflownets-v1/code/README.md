# Property-guided IDR generation with ProGen2 and GFlowNets

Train a GFlowNet to guide ProGen2 toward specified properties of intrinsically
disordered regions (IDRs). The working method jointly adapts ProGen's final
transformer block and a guidance head. The original
ProGen remains fixed as the reward reference.

**The matched supervised baseline has higher success on our current charge task.**
With the ninefold success bonus and 1,600 updates, GFN achieves 38.87% joint
charge/predicted-disorder success versus 44.43% for weighted supervision.
Excluding exact training copies gives 38.87% versus 38.83%, with 1,194 versus
1,193 distinct successful noncopy tails. GFN preserves more predicted disorder
and overall unique outputs; a substantial GFN advantage remains unestablished.
See the [completed results](https://github.com/ayan-goel/progen-gflownets-v1/blob/HEAD/docs/experiments/idr_success_bonus_2026-11-19.md#completed-results--19-november-2026).
New-prompt/property generalization and experimental biology remain untested.

Completed stability experiments established accurate local reward-target sampling
and showed how reward errors limit biological utility. They remain the foundation
for this pivot; a general cross-protein stability predictor is no longer the
immediate objective.

- [Current status](https://github.com/ayan-goel/progen-gflownets-v1/blob/HEAD/docs/status.md)
- [Scientific specification](https://github.com/ayan-goel/progen-gflownets-v1/blob/HEAD/spec.md)
- [Pivot rationale](https://github.com/ayan-goel/progen-gflownets-v1/blob/HEAD/docs/decisions.md#2026-11-12--pivot-from-stability-to-idpidr-property-guidance)
- [Research documentation](https://github.com/ayan-goel/progen-gflownets-v1/blob/HEAD/docs/README.md)
- [Completed experiment records](https://github.com/ayan-goel/progen-gflownets-v1/blob/HEAD/docs/experiments/README.md)
- [Existing commands](https://github.com/ayan-goel/progen-gflownets-v1/blob/HEAD/docs/runbook.md)

## Completed experiment: stronger success reward

Three fresh paired repeats gave both methods a ninefold multiplier for meeting
charge and predicted-disorder requirements. All three miss the declared GFN
advantage screen. The [protocol and results](https://github.com/ayan-goel/progen-gflownets-v1/blob/HEAD/docs/experiments/idr_success_bonus_2026-11-19.md)
preserve the design, outcome and next decision. No further GPU experiment is
selected. The [runbook](https://github.com/ayan-goel/progen-gflownets-v1/blob/HEAD/docs/runbook.md#matched-training-with-a-success-bonus)
provides optional reproduction and scoring recovery.

## Completed comparison: GFlowNet versus supervised training

The experiment compared GFlowNet against reward-weighted supervised learning, with the
same adaptable ProGen block/head, initialization and scored training examples.
Two fresh repeats trained each method for 1,600 updates and evaluated fresh IDRs.
See the [protocol](https://github.com/ayan-goel/progen-gflownets-v1/blob/HEAD/docs/experiments/idr_objective_2026-11-18.md) for the comparison
and its limits.

Keep 1,600 as the working budget and weighted supervision as a comparator.
The [runbook](https://github.com/ayan-goel/progen-gflownets-v1/blob/HEAD/docs/runbook.md#matched-gfn-versus-supervised-training)
provides optional reproduction and saved-generation scoring recovery.

Real experiments and all GPU execution are user initiated. Development checks
are CPU-only:

~~~bash
scripts/test_cpu.sh
~~~
