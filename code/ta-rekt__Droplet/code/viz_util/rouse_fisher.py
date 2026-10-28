"""
Fisher information carried by the escape and return sojourns of the Rouse chain.

Used by viz-rouse.ipynb: "Polymer Fisher rates" (F_poly vs concentration) and
"Fisher rate versus genomic separation" (F_poly vs s).

IDEA
The polymer alternates between two kinds of sojourns:
  * escape : from contact until the chain end first reaches the dividing surface x_div
  * return : from x_div back to contact
One escape followed by one return is a "cycle".  How much one cycle tells us about the
bulk concentration phi is the Fisher information of the sojourn-time distribution:

      I(phi) = average of score(t)^2,      score(t) = d log P(t | phi) / d phi

The Rouse runs are done at discrete concentrations, so the derivative is a finite
difference between two of them, phi_base < phi_plus, and the average runs over the
sojourns simulated at phi_base:

      score(t) = [ log P(t | phi_plus) - log P(t | phi_base) ] / (phi_plus - phi_base)

For each channel (escape, return) we report
      information per cycle     J = phi_base^2 * I        (information about log phi)
      Fisher information rate   F = J / cycle_time,       cycle_time = <escape> + <return>

HOW EACH PIECE IS ESTIMATED
  1. P(t | phi) is a plain histogram -- no smoothing -- on bins that each contain the same
     number of pooled (phi_base + phi_plus) sojourns.  The result hardly depends on the
     number of bins (J changes by <= 2% between 10 and 160 bins at 20261013 high_c).
  2. Cross-fitting.  The sojourns at phi_base are split into folds.  A sojourn's score is
     computed from a histogram that does NOT contain its own fold, so no sojourn is scored
     by a histogram it helped to build.
  3. Noise floor.  Counting noise in the histograms makes score^2 positive even when P
     does not depend on phi at all.  We measure this floor by shuffling the phi labels of
     the sojourns (which removes any real phi dependence) and subtract it.
  4. Error bars.  Block bootstrap of the whole calculation: blocks of consecutive cycles of
     one run are redrawn with replacement and everything is recomputed.  As a check, each
     independent run is also analysed on its own; the spread between runs should match the
     bootstrap error (their ratio is printed).
  5. Trend along the sweep.  A smoothing spline through log10 F, weighted by the error
     bars, with the smoothest curve that still agrees with the points (reduced chi^2 <= 1).
     Its band is the spread of the same spline refitted to every bootstrap replicate.
  6. Quality cuts.  A point with fewer than `min_cycles_per_phi` cycles at either
     concentration is skipped, and a channel is reported only where J is at least
     `min_signal_to_noise` standard errors above zero.

The same finite-difference information computed from the 1D FPT solver's exact PDFs is
reported next to J for reference only; it does not enter the estimate.

USAGE
    settings = FisherSettings()                      # or override any field
    pairs = [dict(runs_base=..., runs_plus=..., solver_base=..., solver_plus=...,
                  phi_base=..., phi_plus=...), ...]  # one per point of the sweep (or None)
    results = analyse_sweep(x_values, pairs, settings, x_name="s")

`runs_base` / `runs_plus` are lists of Rouse output records (one dict per independent
run, with "escape_times" and "return_times" in seconds); `solver_base` / `solver_plus`
are the matching per-phi records of fisher_rate_from_fpt_solver.json.  The contents of
`results` are documented in analyse_sweep.
"""

from dataclasses import dataclass, field

import numpy as np

CHANNELS = ("escape", "return")


@dataclass
class FisherSettings:
    """Every tunable number of the estimator, with the values used in viz-rouse.ipynb."""
    time_rescale: float = 50          # simulated seconds -> plotting time units (= `coef`, cell 97)
    n_histogram_bins: int = 40        # equal-count bins of the sojourn-time histograms
    min_cycles_per_phi: int = 300     # skip a point with fewer cycles at either concentration
    n_folds: int = 5                  # cross-fitting folds
    block_length: int = 50            # consecutive cycles per bootstrap block
    n_bootstrap: int = 200            # bootstrap replicates
    n_permutations: int = 50          # label shuffles used to measure the noise floor
    min_signal_to_noise: float = 3    # report a channel only where J >= this many standard errors
    spline_lambdas: np.ndarray = field(default_factory=lambda: np.logspace(-10, 1, 23))
    spline_subdivision: int = 6       # spline grid points per interval between neighbouring points
    seed: int = 20261013              # seed of the random numbers (shuffles, bootstrap)
    rng: np.random.Generator = field(init=False, repr=False)

    def __post_init__(self):
        self.rng = np.random.default_rng(self.seed)


# -----------------------------------------------------------------------------------------
# 1. Sojourn samples
# -----------------------------------------------------------------------------------------
def cycles_per_run(runs):
    """Escape and return times of every run, trimmed to whole cycles.

    A run usually stops in the middle of a cycle, so its escape and return lists can differ
    by one entry.  Trimming both to the shorter length keeps escape i and return i in the
    same cycle.  Returns {"escape": [array per run], "return": [array per run]}.
    """
    times = {"escape": [], "return": []}
    for run in runs:
        n_cycles = min(len(run["escape_times"]), len(run["return_times"]))
        times["escape"].append(np.asarray(run["escape_times"][:n_cycles], float))
        times["return"].append(np.asarray(run["return_times"][:n_cycles], float))
    return times


def block_labels(n_cycles_per_run, settings):
    """Block number of every cycle of the concatenated runs.

    A block is `block_length` consecutive cycles of ONE run; blocks never span two runs.
    Bootstrapping whole blocks keeps any correlation between successive cycles intact.
    Block b also defines the cross-fitting fold of its cycles: fold = b % n_folds.
    """
    labels, next_block = [], 0
    for n_cycles in n_cycles_per_run:
        block_within_run = np.arange(n_cycles) // settings.block_length
        labels.append(block_within_run + next_block)
        if n_cycles > 0:
            next_block += block_within_run[-1] + 1
    return np.concatenate(labels)


# -----------------------------------------------------------------------------------------
# 2. Fisher information from histograms
# -----------------------------------------------------------------------------------------
def equal_count_bin_edges(times_base, times_plus, settings):
    """Bin edges such that every bin holds the same number of pooled sojourns.

    The outermost edges are -inf and +inf so that every sojourn, also in a bootstrap
    resample, falls into a bin.
    """
    pooled = np.concatenate([times_base, times_plus])
    quantiles = np.linspace(0, 1, settings.n_histogram_bins + 1)[1:-1]
    return np.concatenate([[-np.inf], np.quantile(pooled, quantiles), [np.inf]])


def histogram_probabilities(times, edges, settings):
    """Fraction of the sojourns `times` that falls into each bin."""
    bin_of_sojourn = np.searchsorted(edges, times) - 1
    return np.bincount(bin_of_sojourn, minlength=settings.n_histogram_bins) / times.size


def cross_fitted_fisher_info(times_base, fold_of_sojourn, times_plus, edges, delta_phi, settings):
    """Estimate of I: the squared score averaged over the sojourns at phi_base.

    For a sojourn t in fold k,
        score(t) = [log P_plus(bin of t) - log P_base(bin of t)] / delta_phi,
    where P_base is the histogram of the phi_base sojourns OUTSIDE fold k.  A sojourn in a
    bin that is empty in either histogram has no finite score and contributes zero.
    """
    bin_of_sojourn = np.searchsorted(edges, times_base) - 1
    prob_plus = histogram_probabilities(times_plus, edges, settings)
    squared_score = np.zeros(times_base.size)
    for fold in range(settings.n_folds):
        in_fold = fold_of_sojourn == fold
        prob_base = histogram_probabilities(times_base[~in_fold], edges, settings)   # fold left out
        p_base = prob_base[bin_of_sojourn[in_fold]]
        p_plus = prob_plus[bin_of_sojourn[in_fold]]
        both_nonzero = (p_base > 0) & (p_plus > 0)
        log_ratio = np.log(np.where(both_nonzero, p_plus, 1.0) / np.where(both_nonzero, p_base, 1.0))
        squared_score[in_fold] = (log_ratio / delta_phi) ** 2
    return squared_score.mean()


def permutation_noise_floor(times_base, times_plus, edges, delta_phi, settings):
    """Average estimate of I after shuffling which sojourn belongs to which concentration.

    Shuffling removes any real dependence on phi, so what remains is the contribution of
    the histograms' counting noise to the squared score.  It is subtracted from the real
    estimate.
    """
    pooled = np.concatenate([times_base, times_plus])
    n_base = times_base.size
    floor_estimates = []
    for _ in range(settings.n_permutations):
        shuffled = pooled[settings.rng.permutation(pooled.size)]
        random_folds = settings.rng.integers(settings.n_folds, size=n_base)
        floor_estimates.append(cross_fitted_fisher_info(shuffled[:n_base], random_folds,
                                                        shuffled[n_base:], edges, delta_phi, settings))
    return np.mean(floor_estimates)


def exact_1d_fisher_info(solver_base, solver_plus, channel, delta_phi):
    """The same finite-difference information from the 1D FPT solver's exact PDFs.

    phi_plus's PDF is interpolated (in log-log) onto phi_base's time grid; the numerically
    empty tails are ignored.  Reference only.
    """
    t = np.asarray(solver_base[f"{channel}_t"])
    pdf_base = np.asarray(solver_base[f"{channel}_pdf"])
    pdf_plus = np.exp(np.interp(np.log(t), np.log(solver_plus[f"{channel}_t"]),
                                np.log(np.maximum(solver_plus[f"{channel}_pdf"], 1e-300))))
    tiny = pdf_base.max() * 1e-10
    usable = (pdf_base > tiny) & (pdf_plus > tiny)
    log_ratio = np.log(np.maximum(pdf_plus, tiny)) - np.log(np.maximum(pdf_base, tiny))
    return np.trapezoid(np.where(usable, pdf_base * (log_ratio / delta_phi) ** 2, 0.0), t)


# -----------------------------------------------------------------------------------------
# 3. Bootstrap
# -----------------------------------------------------------------------------------------
def resample_blocks(block_of_cycle, settings):
    """Cycle indices of one block-bootstrap resample.

    Blocks are drawn with replacement separately within each cross-fitting fold, so every
    drawn block keeps its fold and every fold keeps its number of blocks.
    """
    n_blocks = block_of_cycle[-1] + 1
    first_cycle = np.searchsorted(block_of_cycle, np.arange(n_blocks))    # blocks are contiguous
    end_cycle = np.append(first_cycle[1:], block_of_cycle.size)
    drawn_blocks = []
    for fold in range(min(settings.n_folds, n_blocks)):
        blocks_in_fold = np.arange(fold, n_blocks, settings.n_folds)
        drawn_blocks.append(settings.rng.choice(blocks_in_fold, size=blocks_in_fold.size))
    return np.concatenate([np.arange(first_cycle[b], end_cycle[b]) for b in np.concatenate(drawn_blocks)])


# -----------------------------------------------------------------------------------------
# 4. Smooth trend along the sweep
# -----------------------------------------------------------------------------------------
def smoothing_spline(x, y, weights, smoothing, settings):
    """Weighted cubic smoothing spline through the points (x, y), with x evenly spaced.

    Minimises   sum_i w_i (y_i - f(x_i))^2 / sum_i w_i  +  smoothing * integral f''^2,
    with f represented on a fine grid (`spline_subdivision` points per interval, x mapped
    to [0, 1]).  Points with zero weight are ignored.  Returns the fitted values at the data
    points and the effective number of parameters (n_points when the curve goes through
    every point, 2 for a straight line).  Same algorithm as fit_smooth in cell 110.
    """
    n_points = len(x)
    n_grid = (n_points - 1) * settings.spline_subdivision + 1
    spacing = 1.0 / (n_grid - 1)
    data_index = np.arange(n_points) * settings.spline_subdivision    # grid point of each data point

    pick_data = np.zeros((n_points, n_grid))                          # grid values -> data points
    pick_data[np.arange(n_points), data_index] = 1.0
    second_difference = (np.eye(n_grid) - 2 * np.eye(n_grid, k=1) + np.eye(n_grid, k=2))[:-2]
    norm_weights = weights / weights.sum()

    lhs = (pick_data.T @ (norm_weights[:, None] * pick_data)
           + smoothing * spacing * (second_difference.T @ second_difference) / spacing**4)
    lhs += 1e-12 * np.trace(lhs) / n_grid * np.eye(n_grid)           # numerical safety
    hat = np.linalg.solve(lhs, pick_data.T * norm_weights)           # grid values = hat @ y
    fitted_on_grid = hat @ np.where(weights > 0, y, 0.0)
    return fitted_on_grid[data_index], float(np.trace(hat[data_index]))


def trend_along(x_values, rate, rate_se, rate_bootstrap, settings):
    """Smooth trend of a positive rate along an evenly spaced sweep, and its +-1 sigma band.

    The spline is fitted to log10(rate) with weights 1 / SE^2 (SE in log10 units).  Of the
    smoothing strengths in `spline_lambdas` (increasing), the largest whose fit still has a
    reduced chi^2 <= 1 is used: the curve follows the points as closely as their error bars
    require and no closer.  The band is the spread of the same spline, with the same
    smoothing, refitted to every bootstrap replicate of the rates (n_points x n_bootstrap).
    Returns dict(fit, band (2 x n_points), fit_edf, fit_chi2, fit_replicates); NaN where
    there is no data.  fit_replicates (n_bootstrap x n_points) holds the refitted curve of
    every replicate, so that anything built on the trend can be propagated replicate by
    replicate.
    """
    n_points = len(x_values)
    trend = dict(fit=np.full(n_points, np.nan), band=np.full((2, n_points), np.nan),
                 fit_edf=np.nan, fit_chi2=np.nan,
                 fit_replicates=np.full((rate_bootstrap.shape[1], n_points), np.nan))
    se_log10 = rate_se / (rate * np.log(10))
    weights = np.where(np.isfinite(se_log10) & (rate > 0), 1 / se_log10**2, 0.0)
    n_used = np.sum(weights > 0)
    if n_used < 3:
        return trend                                                  # too few points for a trend
    log_rate = np.log10(np.where(weights > 0, rate, 1.0))

    def reduced_chi2(fitted, edf, min_dof):
        return np.sum(weights * (log_rate - fitted) ** 2) / max(n_used - edf, min_dof)

    fits = [smoothing_spline(x_values, log_rate, weights, lam, settings) for lam in settings.spline_lambdas]
    acceptable = [k for k, (fitted, edf) in enumerate(fits) if reduced_chi2(fitted, edf, 1e-9) <= 1]
    chosen = max(acceptable, default=0)                               # the smoothest acceptable fit
    fitted, edf = fits[chosen]

    replicate_fits = []
    for replicate in rate_bootstrap.T:
        replicate_weights = weights * (replicate > 0)
        log_replicate = np.log10(np.where(replicate_weights > 0, replicate, 1.0))
        replicate_fits.append(smoothing_spline(x_values, log_replicate, replicate_weights,
                                               settings.spline_lambdas[chosen], settings)[0])
    spread = np.std(replicate_fits, axis=0, ddof=1)

    trend["fit"] = np.where(weights > 0, 10**fitted, np.nan)
    trend["band"] = np.where(weights > 0, [10**(fitted - spread), 10**(fitted + spread)], np.nan)
    trend["fit_edf"], trend["fit_chi2"] = edf, reduced_chi2(fitted, edf, 1)
    trend["fit_replicates"] = np.where(weights > 0, 10**np.array(replicate_fits), np.nan)
    return trend


# -----------------------------------------------------------------------------------------
# 5. Analysis of one concentration pair
# -----------------------------------------------------------------------------------------
def rate_of_single_run(times_base, times_plus, escape_base, return_base, phi_base, delta_phi, settings):
    """F from one run alone, with its own bins, folds and noise floor.

    Only used to check that the bootstrap error matches the spread between runs.
    """
    edges = equal_count_bin_edges(times_base, times_plus, settings)
    folds = (np.arange(times_base.size) // settings.block_length) % settings.n_folds
    info = phi_base**2 * (cross_fitted_fisher_info(times_base, folds, times_plus, edges, delta_phi, settings)
                          - permutation_noise_floor(times_base, times_plus, edges, delta_phi, settings))
    return info / ((escape_base.mean() + return_base.mean()) * settings.time_rescale)


def analyse_pair(runs_base, runs_plus, solver_base, solver_plus, phi_base, phi_plus, settings):
    """Every estimate for one pair of concentrations, or None if there are too few cycles.

    runs_base, runs_plus      : the Rouse runs at phi_base and phi_plus (one dict per run)
    solver_base, solver_plus  : the 1D solver's per-phi records at the two concentrations
    """
    delta_phi = phi_plus - phi_base

    # ---- sojourn samples: all runs pooled ----
    per_run_base, per_run_plus = cycles_per_run(runs_base), cycles_per_run(runs_plus)
    times_base = {ch: np.concatenate(per_run_base[ch]) for ch in CHANNELS}
    times_plus = {ch: np.concatenate(per_run_plus[ch]) for ch in CHANNELS}
    n_cycles_base, n_cycles_plus = times_base["escape"].size, times_plus["escape"].size
    if min(n_cycles_base, n_cycles_plus) < settings.min_cycles_per_phi:
        return None
    block_base = block_labels([run.size for run in per_run_base["escape"]], settings)
    block_plus = block_labels([run.size for run in per_run_plus["escape"]], settings)
    fold_base = block_base % settings.n_folds

    # ---- cycle time and closed probability at phi_base ----
    mean_escape = times_base["escape"].mean() * settings.time_rescale
    mean_return = times_base["return"].mean() * settings.time_rescale
    cycle_time = mean_escape + mean_return

    # ---- information per cycle: raw estimate minus the noise floor ----
    edges = {ch: equal_count_bin_edges(times_base[ch], times_plus[ch], settings) for ch in CHANNELS}
    raw_info = {ch: cross_fitted_fisher_info(times_base[ch], fold_base, times_plus[ch], edges[ch],
                                             delta_phi, settings) for ch in CHANNELS}
    noise_floor = {ch: permutation_noise_floor(times_base[ch], times_plus[ch], edges[ch], delta_phi, settings)
                   for ch in CHANNELS}
    info = {ch: phi_base**2 * (raw_info[ch] - noise_floor[ch]) for ch in CHANNELS}

    # ---- block bootstrap of all of the above (the noise floor is kept fixed) ----
    boot_info = {ch: np.empty(settings.n_bootstrap) for ch in CHANNELS}
    boot_escape, boot_return = np.empty(settings.n_bootstrap), np.empty(settings.n_bootstrap)
    for b in range(settings.n_bootstrap):
        pick_base, pick_plus = resample_blocks(block_base, settings), resample_blocks(block_plus, settings)
        boot_escape[b] = times_base["escape"][pick_base].mean() * settings.time_rescale
        boot_return[b] = times_base["return"][pick_base].mean() * settings.time_rescale
        for ch in CHANNELS:
            raw = cross_fitted_fisher_info(times_base[ch][pick_base], fold_base[pick_base],
                                           times_plus[ch][pick_plus], edges[ch], delta_phi, settings)
            boot_info[ch][b] = phi_base**2 * (raw - noise_floor[ch])
    boot_cycle_time = boot_escape + boot_return

    result = dict(phi_base=phi_base, phi_plus=phi_plus,
                  n_cycles_base=n_cycles_base, n_cycles_plus=n_cycles_plus,
                  mean_escape_time=mean_escape, mean_return_time=mean_return,
                  cycle_time=cycle_time, cycle_time_se=boot_cycle_time.std(ddof=1),
                  p_closed=mean_escape / cycle_time, p_closed_bootstrap=boot_escape / boot_cycle_time)
    result["p_closed_se"] = result["p_closed_bootstrap"].std(ddof=1)

    # ---- per channel: keep it only if J is clearly above zero ----
    for ch in CHANNELS:
        info_se = boot_info[ch].std(ddof=1)
        channel = dict(info_1d=phi_base**2 * exact_1d_fisher_info(solver_base, solver_plus, ch, delta_phi),
                       noise_floor_fraction=noise_floor[ch] / raw_info[ch])
        if info[ch] >= settings.min_signal_to_noise * info_se:
            channel["info"], channel["info_se"] = info[ch], info_se
            channel["rate"] = info[ch] / cycle_time
            channel["rate_bootstrap"] = boot_info[ch] / boot_cycle_time
            channel["rate_se"] = channel["rate_bootstrap"].std(ddof=1)
            channel["rate_per_run"] = [
                rate_of_single_run(base, plus, escape, ret, phi_base, delta_phi, settings)
                for base, plus, escape, ret in zip(per_run_base[ch], per_run_plus[ch],
                                                   per_run_base["escape"], per_run_base["return"])]
        result[ch] = channel
    return result


# -----------------------------------------------------------------------------------------
# 6. A whole sweep (over s or over phi), with its trend and a summary table
# -----------------------------------------------------------------------------------------
POINT_KEYS = ("phi_base", "phi_plus", "n_cycles_base", "n_cycles_plus", "mean_escape_time",
              "mean_return_time", "cycle_time", "cycle_time_se", "p_closed", "p_closed_se")
CHANNEL_KEYS = ("info", "info_se", "info_1d", "noise_floor_fraction", "rate", "rate_se")


def analyse_sweep(x_values, pairs, settings, x_name="x"):
    """Analyse every point of a sweep and fit the trend of F along it.

    x_values : the swept variable (s or phi_base), evenly spaced
    pairs    : one dict per point with runs_base, runs_plus, solver_base, solver_plus,
               phi_base, phi_plus (see analyse_pair), or None for a point without data
    Returns a dict with, for every point:
        x                             the swept values
        phi_base, phi_plus            the two concentrations of the pair
        n_cycles_base, n_cycles_plus  cycles available at each concentration
        mean_escape_time              mean escape (= closed) time at phi_base
        mean_return_time              mean return (= open) time at phi_base
        cycle_time, cycle_time_se     their sum and its bootstrap error
        p_closed, p_closed_se         mean_escape_time / cycle_time and its error
        p_closed_bootstrap            its bootstrap replicates (n_points x n_bootstrap)
        "escape" and "return", each a dict of
            info, info_se             J and its error
            info_1d                   J of the 1D solver (reference)
            noise_floor_fraction      subtracted noise floor / raw estimate
            rate, rate_se             F and its error
            rate_bootstrap            bootstrap replicates of F (n_points x n_bootstrap); they
                                      share their draws with p_closed_bootstrap
            rate_per_run              F of each run analysed alone (n_points x n_runs)
            fit, band, fit_edf, fit_chi2   smooth trend of F along x (see trend_along)
            fit_replicates            the trend refitted to every bootstrap replicate
    Times are in plotting units (simulated seconds x time_rescale).  NaN marks a skipped
    point or an unresolved channel.
    """
    n_points = len(x_values)
    n_runs = max(len(pair["runs_base"]) for pair in pairs if pair is not None)

    results = dict(x=np.asarray(x_values), p_closed_bootstrap=np.full((n_points, settings.n_bootstrap), np.nan))
    for key in POINT_KEYS:
        results[key] = np.full(n_points, np.nan)
    for ch in CHANNELS:
        results[ch] = {key: np.full(n_points, np.nan) for key in CHANNEL_KEYS}
        results[ch]["rate_bootstrap"] = np.full((n_points, settings.n_bootstrap), np.nan)
        results[ch]["rate_per_run"] = np.full((n_points, n_runs), np.nan)

    for i, (x, pair) in enumerate(zip(x_values, pairs)):
        point = None if pair is None else analyse_pair(settings=settings, **pair)
        if point is None:
            print(f"{x_name} = {x:g}: fewer than {settings.min_cycles_per_phi} cycles at one of the "
                  f"concentrations -- skipped")
            continue
        for key in POINT_KEYS + ("p_closed_bootstrap",):
            results[key][i] = point[key]
        for ch in CHANNELS:
            for key, value in point[ch].items():
                if key == "rate_per_run":
                    results[ch][key][i, :len(value)] = value
                else:
                    results[ch][key][i] = value

    for ch in CHANNELS:
        results[ch].update(trend_along(results["x"], results[ch]["rate"], results[ch]["rate_se"],
                                       results[ch]["rate_bootstrap"], settings))

    print_summary(results, x_name, n_runs, settings)
    return results


def print_summary(results, x_name, n_runs, settings):
    """One line per point with the cycle time and, per channel, J and F with relative errors."""
    print(f"{n_runs} run(s) pooled; {settings.n_histogram_bins} equal-count bins\n")
    channel_header = f"{'floor':>5} {'J':>7} {'J_1D':>6} {'SE_J':>5} {'F':>7} {'SE_F':>5}"
    print(f"{x_name:>10} {'cycles':>6} {'':>6} {'cycle time':>15} | {'escape':^47} | {'return':^47}")
    print(f"{'':>10} {'base':>6} {'plus':>6} {'mean':>8} {'SE':>6} | {channel_header} | {channel_header}")
    for i, x in enumerate(results["x"]):
        if not np.isfinite(results["cycle_time"][i]):
            continue
        line = (f"{x:>10.4g} {results['n_cycles_base'][i]:6.0f} {results['n_cycles_plus'][i]:6.0f} "
                f"{results['cycle_time'][i]:8.4f} {results['cycle_time_se'][i] / results['cycle_time'][i]:6.3f}")
        for ch in CHANNELS:
            c = results[ch]
            line += (f" | {c['noise_floor_fraction'][i]:5.2f} {c['info'][i]:7.3f} {c['info_1d'][i]:6.2f} "
                     f"{c['info_se'][i] / c['info'][i]:5.3f} {c['rate'][i]:7.2f} {c['rate_se'][i] / c['rate'][i]:5.3f}")
        print(line)
    print("\nfloor = subtracted noise floor / raw estimate; J_1D = 1D-solver value (reference);"
          " SE columns are relative errors")

    for ch in CHANNELS:
        c = results[ch]
        if np.isfinite(c["fit_edf"]):
            message = f"{ch}: trend has {c['fit_edf']:.1f} effective parameters, reduced chi^2 {c['fit_chi2']:.2f}"
        else:
            message = f"{ch}: no trend (fewer than 3 resolved points)"
        resolved = np.isfinite(c["rate"])
        if n_runs > 1 and resolved.any():
            run_spread = c["rate_per_run"][resolved].std(axis=1, ddof=1) / np.sqrt(n_runs)
            message += (f";  spread between runs / sqrt({n_runs}) = "
                        f"{np.nanmean(run_spread / c['rate_se'][resolved]):.2f} x bootstrap error")
        unresolved = [f"{x:g}" for i, x in enumerate(results["x"])
                      if np.isfinite(results["cycle_time"][i]) and not resolved[i]]
        if unresolved:
            message += (f";  not resolved (J < {settings.min_signal_to_noise:g} SE) at "
                        f"{x_name} = {', '.join(unresolved)}")
        print(message)
