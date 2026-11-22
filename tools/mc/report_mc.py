#!/usr/bin/env python3
"""
Summarize a Monte Carlo campaign written by run_mc.py.

    python3 tools/mc/report_mc.py mc_runs.csv [--p95-limit 10] [--handoff-limit-orbits 3]

Prints: pass rates against the requirements below, percentiles of each metric, which inputs the
outcomes are most sensitive to (Spearman rank correlation), and the worst runs (to re-run in
Simulink: python3 tools/mc/run_mc.py --start <run> --n 1 ... reproduces one exactly).
Standard library only.
"""
import argparse
import csv
import math

METRICS = ["handoff_s", "acquire_s", "point_frac_settled", "nadir_rms", "nadir_p95", "nadir_max", "est_rms",
           "est_p95", "rate_rms_degps", "resets", "rej_gated"]
NOT_INPUTS = set(METRICS) | {"run", "sim_seed", "returncode", "error", "point_frac", "uplinks", "tle_resid_m",
                             "magcal_ready", "magcal_uT", "coil_sign",
                             "coil_on", "rej_gyro_stuck", "rej_gyro_spike", "rej_mag_stuck", "rej_mag_norm", "rej_mag_nan",
                             "rej_pd", "rej_pd_diodes", "tle_accepted", "tle_rejected", "tle_last_reject",
                             "orbit_fail_streak", "orbit_fail_reason"}


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def pct(xs, q):
    xs = sorted(xs)
    if not xs:
        return float("nan")
    k = (len(xs) - 1) * q
    f, c = math.floor(k), math.ceil(k)
    return xs[f] if f == c else xs[f] + (xs[c] - xs[f]) * (k - f)


def ranks(xs):
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    r = [0.0] * len(xs)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and xs[order[j + 1]] == xs[order[i]]:
            j += 1
        for k in range(i, j + 1):
            r[order[k]] = (i + j) / 2.0
        i = j + 1
    return r


def spearman(x, y):
    pairs = [(a, b) for a, b in zip(x, y) if not math.isnan(a) and not math.isnan(b)]
    if len(pairs) < 10:
        return float("nan")
    rx, ry = ranks([p[0] for p in pairs]), ranks([p[1] for p in pairs])
    n = len(pairs)
    mx, my = sum(rx) / n, sum(ry) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    sxx = sum((a - mx) ** 2 for a in rx)
    syy = sum((b - my) ** 2 for b in ry)
    return sxy / math.sqrt(sxx * syy) if sxx > 0 and syy > 0 else float("nan")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("--p95-limit", type=float, default=10.0, help="settled nadir p95 requirement (deg)")
    ap.add_argument("--handoff-limit-orbits", type=float, default=3.0)
    ap.add_argument("--orbit-s", type=float, default=5700.0)
    ap.add_argument("--worst", type=int, default=10)
    a = ap.parse_args()
    rows = list(csv.DictReader(open(a.csv)))
    ok = [r for r in rows if r.get("returncode") == "0" and r.get("nadir_rms", "") != ""]
    print(f"{len(rows)} runs, {len(rows) - len(ok)} failed to run")
    if not ok:
        return

    # ---- requirements ----
    lim_h = a.handoff_limit_orbits * a.orbit_s
    handed = [r for r in ok if fnum(r["handoff_s"]) <= lim_h]
    settled_ok = [r for r in ok if fnum(r["nadir_p95"]) <= a.p95_limit]
    no_reset = [r for r in ok if fnum(r["resets"]) == 0]
    n = len(ok)
    print("\nRequirement pass rates")
    print(f"  hand-off to pointing within {a.handoff_limit_orbits:g} orbits : {len(handed)}/{n} ({100 * len(handed) / n:.1f}%)")
    print(f"  settled nadir p95 <= {a.p95_limit:g} deg              : {len(settled_ok)}/{n} ({100 * len(settled_ok) / n:.1f}%)")
    print(f"  no filter resets                         : {len(no_reset)}/{n} ({100 * len(no_reset) / n:.1f}%)")
    never = [r for r in ok if math.isnan(fnum(r["handoff_s"]))]
    print(f"  never handed off                         : {len(never)}/{n}")

    # ---- percentiles ----
    print("\nMetric percentiles (over runs; NaN = never happened, excluded)")
    print(f"  {'metric':20s} {'n':>5s} {'P50':>10s} {'P90':>10s} {'P99':>10s} {'max':>10s}")
    for m in METRICS:
        xs = [fnum(r.get(m)) for r in ok]
        xs = [x for x in xs if not math.isnan(x)]
        if xs:
            print(f"  {m:20s} {len(xs):5d} {pct(xs, .5):10.3g} {pct(xs, .9):10.3g} {pct(xs, .99):10.3g} {max(xs):10.3g}")

    # ---- sensitivity ----
    inputs = [k for k in ok[0].keys() if k not in NOT_INPUTS]
    print("\nSensitivity: Spearman rank correlation (|rho| > 0.1 shown, top 8 per outcome)")
    for m in ["handoff_s", "acquire_s", "nadir_p95", "est_rms"]:
        y = [fnum(r.get(m)) for r in ok]
        if m in ("handoff_s", "acquire_s"):  # never = worst
            y = [1e9 if math.isnan(v) else v for v in y]
        rho = []
        for k in inputs:
            x = [fnum(r.get(k)) for r in ok]
            if all(math.isnan(v) for v in x):
                continue
            # magnitude matters for signed vectors (bias, residual dipole, ...)
            if any(k.startswith(s) for s in ("mag_bias", "gyro_bias0", "m_res", "sf_err", "w0_vec")):
                x = [abs(v) for v in x]
            c = spearman(x, y)
            if not math.isnan(c) and abs(c) > 0.1:
                rho.append((abs(c), c, k))
        rho.sort(reverse=True)
        print(f"  {m}: " + (", ".join(f"{k} {c:+.2f}" for _, c, k in rho[:8]) or "nothing above 0.1"))

    # ---- worst runs ----
    def badness(r):
        h = fnum(r["handoff_s"])
        return (0 if not math.isnan(h) else 1, fnum(r["nadir_p95"]) if not math.isnan(fnum(r["nadir_p95"])) else 1e9)
    worst = sorted(ok, key=badness, reverse=True)[: a.worst]
    print(f"\nWorst {len(worst)} runs (rerun one: run_mc.py --start RUN --n 1 --seed <campaign seed>)")
    print(f"  {'run':>5s} {'handoff_s':>9s} {'acq_s':>8s} {'nadir_p95':>9s} {'est_rms':>7s} {'resets':>6s} "
          f"{'|w0|':>5s} {'|gyro_b|':>8s} {'alt':>5s} {'inc':>5s}")
    for r in worst:
        w0 = math.sqrt(sum(fnum(r[f"w0_vec_{j}"]) ** 2 for j in range(3)))
        gb = math.sqrt(sum(fnum(r[f"gyro_bias0_{j}"]) ** 2 for j in range(3)))
        print(f"  {r['run']:>5s} {fnum(r['handoff_s']):9.0f} {fnum(r['acquire_s']):8.0f} {fnum(r['nadir_p95']):9.2f} "
              f"{fnum(r['est_rms']):7.2f} {r['resets']:>6s} {w0:5.1f} {gb:8.4f} {fnum(r['alt_km']):5.0f} "
              f"{fnum(r['inc_deg']):5.0f}")


if __name__ == "__main__":
    main()
