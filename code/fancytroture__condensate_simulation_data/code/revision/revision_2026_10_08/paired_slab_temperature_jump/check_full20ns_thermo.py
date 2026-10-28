"""Verify thermostat and fixed periodic volume without touching trajectories."""
import csv
import json
import numpy as np
from extract_full20ns import BASE, CASES, OUT, sha, write_json


def main():
    rows = []
    for case, (_, low, high) in CASES.items():
        for snapshot in range(1, 9):
            for temp in (low, high):
                root = BASE / f"{case}_slab_temperature_jump" / "runs" / f"snapshot{snapshot:02d}_{temp}K"
                identity = json.loads((root / "identity.json").read_text())
                assert identity["snapshot"] == snapshot and identity["temperature_K"] == temp
                values = []
                for part in sorted(root.glob("part[0-9][0-9][0-9][0-9]")):
                    meta = json.loads((part / "complete.json").read_text())
                    assert meta["identity"] == identity
                    assert sha(part / "thermo.csv") == meta["files"]["thermo.csv"]
                    values.append(np.loadtxt(part / "thermo.csv", delimiter=",", comments="#", usecols=range(7)))
                data = np.concatenate(values)
                assert data.shape == (20000, 7) and np.isfinite(data).all()
                np.testing.assert_allclose(data[:, 0], np.arange(1, 20001)*100)
                np.testing.assert_allclose(data[:, 1], np.arange(1, 20001), atol=1e-4)
                np.testing.assert_allclose(data[:, 5], 112000.)
                np.testing.assert_allclose(data[:, 6], data[0, 6], atol=1e-12)
                selected = data[data[:, 1] > 100+1e-5]
                mean_temp = selected[:, 4].mean()
                assert abs(mean_temp-temp) < 2
                rows.append(dict(case=case, snapshot=snapshot, target_temperature_K=temp,
                    mean_temperature_K=mean_temp, temperature_sd_K=selected[:, 4].std(),
                    volume_nm3=data[0, 5], whole_box_density_g_ml=data[0, 6],
                    density_note="Whole-box density is not dense-phase core density",
                    n_thermo_records=len(data), all_hashes_verified=True))
        print(f"Thermo verified: {case}", flush=True)
    with (OUT / "thermostat_and_volume_checks.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    write_json(OUT / "thermostat_validation.json", dict(status="complete", branches=len(rows),
        maximum_absolute_mean_temperature_error_K=max(abs(r["mean_temperature_K"]-r["target_temperature_K"]) for r in rows),
        all_volumes_fixed_nm3=112000., checked_thermo_files=1920, script_sha256=sha(__file__)))


if __name__ == "__main__":
    main()
