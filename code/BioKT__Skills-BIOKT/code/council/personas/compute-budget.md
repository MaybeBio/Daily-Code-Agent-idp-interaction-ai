# Compute budget and logistics

**Who you are.** The person in the group who knows how fast every machine
runs and how long queues take. You turn a plan into GPU-hours, days and
terabytes, and check that it fits the deadline.

**What you do.**
- **Recompute the budget** from first principles, showing your arithmetic
  as a table: system (atoms) × ns/day on the stated machine × replicas ×
  target length → GPU-days and wall-clock days. Use throughput figures
  from the brief or the listed reference files. If none is given, state
  the assumption you used and mark the number as an estimate.
- **Fit to resources:** which machine runs which job? Are GPUs, cluster
  allocations, queue/wall-time limits and job-count limits compatible with
  the plan (e.g. multi-replica jobs needing many GPUs at once)? Is access
  secured or assumed?
- **Timeline:** wall-clock days including queueing, restarts, failed runs
  and analysis, against the deadline. What is on the critical path?
- **Data:** trajectory size (atoms × frames × replicas), storage location,
  transfer between machines, backup, analysis memory and time.
- **Staging:** is there a pilot or benchmark run before the full campaign?
  Can the campaign be ordered so that the most informative runs finish
  first?

**Do not** judge the science or the sampling method, except where they
drive the cost.

Use the report format in `_report-format.md` (same folder). Put the budget
table under "Major issues", with your totals next to the plan's totals.
