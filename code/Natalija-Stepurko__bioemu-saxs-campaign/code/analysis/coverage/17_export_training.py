"""Copy the small training-data outputs into results/coverage/: the per-construct nearest-neighbour table
(all 439 entries), the training-set composition and overlap summary, and the versions file (which now also
records the training-data sources)."""
import json
import shutil

from common import EXPORT, TAB, VERSIONS

EXPORT.mkdir(parents=True, exist_ok=True)
shutil.copyfile(TAB / "training_nearest.csv", EXPORT / "training_nearest.csv")
s = json.loads((TAB / "training_summary.json").read_text())
(EXPORT / "training_summary.json").write_text(json.dumps(s, indent=1) + "\n")
(EXPORT / "versions.json").write_text(json.dumps(json.loads(VERSIONS.read_text()), indent=2) + "\n")
print("exported to", EXPORT, sorted(p.name for p in EXPORT.iterdir()))
