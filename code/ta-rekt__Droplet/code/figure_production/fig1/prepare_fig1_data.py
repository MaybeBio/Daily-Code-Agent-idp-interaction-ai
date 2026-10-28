#!/usr/bin/env python3
"""Reduce the inputs of Fig. 1 to small files in processed_data/.

    python figure_production/fig1/prepare_fig1_data.py --data-root ~/Archive/Work/droplet_data

Reads
  <data-root>/20260423/proximal_promoter_transc_d_50_c_<c>/output_cpu<k>.json
      the proximal-promoter lattice runs (as Fig. 2), for panel E:
        E  trace,  c = 0.006, output_cpu3   (viz.ipynb cell 155, frames 0-40000; its `run = 12`
           does not name the published trace -- cpu3 is the one that matches it: the
           smoothed activity peaks at 0.58 near t = 1.71e5 inside the long distance dip)
  --svg  the paper canvas fig_sources/fig1.svg, for panel F only: F is a hand-drawn
         sketch with no code behind it, so its step traces, shaded intervals and
         interval markers are read back from the drawing's path data.
(B and C, the 3D configurations -- c = 0.003 run 13 and c = 0.017 run 2 -- are
pasted from the source canvas by extract_art.py, so they need no data here.)
Writes
  processed_data/fig1_trace.csv     smoothed transcription, distance, one row per frame;
                                    time = frame * mc_samples_per_frame (fig1_scalars.json)
  processed_data/fig1_sketch.json   panel F, times as fractions of the trace span
  processed_data/fig1_scalars.json  provenance and constants
"""
import argparse, json, os, re, sys
import xml.etree.ElementTree as ET
import numpy as np
from scipy.ndimage import uniform_filter1d

EXP = "20260423/proximal_promoter_transc_d_50_c_{:.3f}/output_cpu{}.json"
TRACE = (0.006, 3)
START, WIDTH = 0, 40000          # cell 155
BOX, RNAP = 62, 0.02             # cell 110: total volume fraction (c + 0.02) * 62^3
SMOOTH_MCS = 500                 # cell 155: window = 500 * int(sampling multiplier) frames


def resolve_repo(start):
    d = os.path.abspath(start)
    while d != os.path.dirname(d):
        if os.path.isfile(os.path.join(d, "Project.toml")) and os.path.isdir(os.path.join(d, "spec")):
            return d
        d = os.path.dirname(d)
    sys.exit("repo root not found above " + start)


def run(root, c, k):
    return json.load(open(os.path.join(root, EXP.format(c, k))))[str(k)]["1"]


# ---- panel F: path data of the sketch -----------------------------------------
NUM = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"


def path_points(d):
    """Vertices of an SVG path made of M/m/L/l/H/h/V/v commands (all F uses)."""
    pts, x, y = [], 0.0, 0.0
    for cmd, args in re.findall(r"([MmLlHhVv])([^MmLlHhVvZz]*)", d):
        v = [float(a) for a in re.findall(NUM, args)]
        if cmd in "Hh":
            for a in v:
                x = a if cmd == "H" else x + a; pts.append((x, y))
        elif cmd in "Vv":
            for a in v:
                y = a if cmd == "V" else y + a; pts.append((x, y))
        else:
            for i in range(0, len(v), 2):
                if cmd in "ML":
                    x, y = v[i], v[i + 1]
                else:
                    x, y = x + v[i], y + v[i + 1]
                pts.append((x, y))
    return np.array(pts)


def sketch(svg):
    """Panel F of the paper canvas, in units of the trace span (0..1 in time,
    0 = Off, 1 = On)."""
    els = {e.get("id"): e for e in ET.parse(svg).getroot().iter() if e.get("id")}
    out = {}
    for tag, trace, on_y, off_y in (("polymer", "path268-6-4-1", 278.79, 291.66),
                                    ("free", "path268-6-2", 317.43, 330.22)):
        p = path_points(els[trace].get("d"))
        out[tag] = {"x": p[:, 0], "level": np.round((off_y - p[:, 1]) / (off_y - on_y)).astype(int)}
    # common time frame: the polymer trace's span
    x0 = min(out["polymer"]["x"].min(), out["free"]["x"].min())
    x1 = max(out["polymer"]["x"].max(), out["free"]["x"].max())
    t = lambda x: (np.asarray(x, float) - x0) / (x1 - x0)
    for tag in out:
        out[tag] = {"t": t(out[tag]["x"]).round(5).tolist(), "level": out[tag]["level"].tolist()}

    def rect(i):
        e = els[i]
        return float(e.get("x")), float(e.get("x")) + float(e.get("width"))
    spans = []
    for i, state in (("rect1370-60-1", "open"), ("rect1370-1", "closed"), ("rect1370-60-4-9", "open"),
                     ("rect1370-6-7", "closed"), ("rect1370-60-3-9", "open")):
        a, b = rect(i)
        spans.append({"state": state, "t0": round(float(t(a)), 5), "t1": round(float(t(b)), 5)})
    first_x = lambda i: path_points(els[i].get("d"))[0, 0]
    out["spans"] = spans
    # interval markers (the faint vertical ticks) and the brackets between them
    out["polymer_marks"] = [round(float(t(first_x(i))), 5) for i in
                            ("path270-0-3-6-5", "path270-0-3-6-3-0", "path270-0-3-6-3-6-5")]
    out["free_marks"] = [round(float(t(first_x(i))), 5) for i in
                         ("path270-0-3-6-8-6", "path270-0-3-6-8-3-7",
                          "path270-0-3-6-3-8-5", "path270-0-3-6-3-6-8-3")]
    out["brackets"] = {"polymer": [["T^x", 0, 1], ["T^o", 1, 2]],
                       "free": [["T^-", 0, 1], ["T^+", 2, 3]]}
    return out


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    repo = resolve_repo(here)
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", default=os.path.expanduser("~/Archive/Work/droplet_data"))
    ap.add_argument("--svg", default=os.path.join(repo, "..", "Tarek_Droplets", "fig_sources", "fig1.svg"))
    ap.add_argument("--out", default=os.path.join(os.path.dirname(here), "processed_data"))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    scal = {"data_root_layout": EXP}

    # ---- E: transcription and probe-enhancer distance -------------------------------
    c, k = TRACE
    r = run(a.data_root, c, k)
    mult = 50000 / ((c + RNAP) * BOX**3)          # MC samples per frame (cell 110)
    act = np.asarray(r["transcription_activity"], float)[START:START + WIDTH]
    dist = np.asarray(r["pattern_distances"], float)[START:START + WIDTH]
    t = np.arange(START, START + WIDTH) * mult / 1e5
    window = SMOOTH_MCS * int(mult)
    smooth = uniform_filter1d(act, size=window)
    with open(os.path.join(a.out, "fig1_trace.csv"), "w") as f:
        f.write("activity_smoothed,distance\n")
        for row in zip(smooth, dist):
            f.write("%.4g,%.4g\n" % row)      # 4 digits: far below a printed line width
    scal["trace"] = {"c": c, "run": k, "start": START, "width": WIDTH,
                     "mc_samples_per_frame": mult, "smoothing_window_frames": window,
                     "smoothing_window_mc_samples": window * mult}
    print("trace  c = %.3f run %2d: %d frames, t to %.3f e5, smoothing %d frames (%.0f MC samples)"
          % (c, k, len(t), t[-1], window, window * mult))

    # ---- F: the sketch -----------------------------------------------------------
    json.dump(sketch(a.svg), open(os.path.join(a.out, "fig1_sketch.json"), "w"), indent=1)
    scal["sketch_source"] = "Tarek_Droplets/fig_sources/fig1.svg (path data of panel F)"
    json.dump(scal, open(os.path.join(a.out, "fig1_scalars.json"), "w"), indent=1)
    print("wrote fig1_* in", os.path.relpath(a.out, repo))


if __name__ == "__main__":
    main()
