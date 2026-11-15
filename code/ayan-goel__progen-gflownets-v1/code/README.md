# Property-guided IDR generation with frozen ProGen2

Learn a small GFlowNet guidance head that steers frozen ProGen2 toward specified
properties of intrinsically disordered regions (IDRs).

**Ready: a matched charge-context experiment.** Compare the original guidance
head with one that gives charge instructions a separate path to its output.
Both use the same frozen ProGen, reward and training budget. The completed
diagnostic motivated this change; whether it improves generation remains to test.
See the [experiment plan](https://github.com/ayan-goel/progen-gflownets-v1/blob/HEAD/docs/experiments/idr_context_2026-11-12.md).

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

## Running

Activate the environment, choose your GPU, and run both heads and evaluation:

~~~bash
source .venv/bin/activate
export CUDA_VISIBLE_DEVICES=0
scripts/run.sh idr-context --output runs/idr-context-v1 --device cuda
~~~

Real experiments and all GPU execution are user initiated. Development checks
are CPU-only:

~~~bash
scripts/test_cpu.sh
~~~

The run writes `report.md` with the comparison and detailed results for each head.
See the [protocol](https://github.com/ayan-goel/progen-gflownets-v1/blob/HEAD/docs/experiments/idr_context_2026-11-12.md) for choices,
interpretation limits, and CPU evaluation recovery.
