#!/usr/bin/env python3
"""
Monte Carlo driver for the host closed loop (build/adcs-plant-sim: plant.c truth + body()).

Each run draws every uncertain truth parameter from the ranges in PARAMS below, runs the sim
for --orbits orbits, and records the inputs plus the sim's SUMMARY line in one CSV row.

    python3 tools/mc/run_mc.py --n 1000 --jobs 10 --out mc_runs.csv
    python3 tools/mc/report_mc.py mc_runs.csv

Reproducible: run i always uses the parameters drawn from (--seed, i), and sim seed = i + 1.
Standard library only.

THE RANGES BELOW ARE ASSUMPTIONS. Each one says where it comes from; replace them with measured /
spec values (and say so in the comment) before quoting results.
"""
import argparse
import csv
import math
import os
import random
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
DEFAULT_SIM = os.path.join(REPO, "build", "adcs-plant-sim")

# Nominal inertia (kg m^2): Simulink Spacecraft Dynamics block (symmetric part) == plant.c J_SC
J_NOM = [[0.029730, 24.4275e-6, 3672.061e-6],
         [24.4275e-6, 36414.034e-6, -285.496e-6],
         [3672.061e-6, -285.496e-6, 12429.899e-6]]


# ---------------- small linear algebra (no numpy) ----------------
def jacobi_eig(A, sweeps=50):
    """Eigen-decomposition of a symmetric 3x3: returns (eigenvalues, V) with A = V diag(w) V^T."""
    a = [row[:] for row in A]
    V = [[1.0 if i == j else 0.0 for j in range(3)] for i in range(3)]
    for _ in range(sweeps):
        off = sum(a[i][j] ** 2 for i in range(3) for j in range(3) if i != j)
        if off < 1e-30:
            break
        for p in range(3):
            for q in range(p + 1, 3):
                if abs(a[p][q]) < 1e-300:
                    continue
                th = 0.5 * math.atan2(2 * a[p][q], a[q][q] - a[p][p])
                c, s = math.cos(th), math.sin(th)
                for k in range(3):
                    akp, akq = a[k][p], a[k][q]
                    a[k][p], a[k][q] = c * akp - s * akq, s * akp + c * akq
                for k in range(3):
                    apk, aqk = a[p][k], a[q][k]
                    a[p][k], a[q][k] = c * apk - s * aqk, s * apk + c * aqk
                for k in range(3):
                    vkp, vkq = V[k][p], V[k][q]
                    V[k][p], V[k][q] = c * vkp - s * vkq, s * vkp + c * vkq
    return [a[i][i] for i in range(3)], V


def matmul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(3)) for j in range(3)] for i in range(3)]


def transpose(A):
    return [[A[j][i] for j in range(3)] for i in range(3)]


def axis_angle(axis, ang):
    x, y, z = axis
    c, s, C = math.cos(ang), math.sin(ang), 1 - math.cos(ang)
    return [[c + x * x * C, x * y * C - z * s, x * z * C + y * s],
            [y * x * C + z * s, c + y * y * C, y * z * C - x * s],
            [z * x * C - y * s, z * y * C + x * s, c + z * z * C]]


def unit_vector(rng):
    z = rng.uniform(-1, 1)
    t = rng.uniform(0, 2 * math.pi)
    r = math.sqrt(1 - z * z)
    return [r * math.cos(t), r * math.sin(t), z]


def random_quaternion(rng):
    """Uniform on SO(3) (Shoemake), scalar first"""
    u1, u2, u3 = rng.random(), rng.random(), rng.random()
    a, b = math.sqrt(1 - u1), math.sqrt(u1)
    return [b * math.cos(2 * math.pi * u3), a * math.sin(2 * math.pi * u2),
            a * math.cos(2 * math.pi * u2), b * math.sin(2 * math.pi * u3)]


# ---------------- parameter draws ----------------
def draw_inertia(rng):
    """Principal moments each scaled U(0.8, 1.2) and principal axes rotated by up to 10 deg about a
    random axis. ASSUMPTION: +-20% / 10 deg covers CAD-vs-flight (harness, deployables, battery).
    Re-drawn until the triangle inequality holds (physically realizable)."""
    w, V = jacobi_eig(J_NOM)
    while True:
        ws = [wi * rng.uniform(0.8, 1.2) for wi in w]
        if all(ws[i] + ws[j] > ws[k] for i, j, k in ((0, 1, 2), (0, 2, 1), (1, 2, 0))):
            break
    R = axis_angle(unit_vector(rng), math.radians(rng.uniform(0, 10)))
    VR = matmul(R, V)
    D = [[ws[i] if i == j else 0.0 for j in range(3)] for i in range(3)]
    J = matmul(matmul(VR, D), transpose(VR))
    return J


def u(rng, lo, hi):
    return rng.uniform(lo, hi)


def draw(i, base_seed):
    """All parameters for run i. Returns (args for adcs-plant-sim, dict of drawn values for the CSV)."""
    rng = random.Random(base_seed * 1_000_003 + i)
    p = {}
    J = draw_inertia(rng)
    p["J"] = [J[0][0], J[1][1], J[2][2], J[0][1], J[0][2], J[1][2]]

    # orbit. ASSUMPTION: rideshare SSO/ISS-like LEO; replace with the launch provider's range.
    p["alt_km"] = u(rng, 450, 650)
    p["inc_deg"] = u(rng, 45, 98)
    p["raan_deg"] = u(rng, 0, 360)
    p["argp_deg"] = u(rng, 0, 360)
    p["ma_deg"] = u(rng, 0, 360)
    p["epoch_offset_days"] = u(rng, 0, 365)   # sun direction / eclipse season

    # deployment. ASSUMPTION: tip-off up to 30 deg/s, any direction, any attitude.
    w0 = u(rng, 0, 30)
    p["w0_vec"] = [w0 * c for c in unit_vector(rng)]
    p["q0"] = random_quaternion(rng)

    # sensors: each parameter independently between the nominal and worst-case values of the
    # Simulink "Sensors (noise)" block (plant.c plant_error_model).
    p["sigma_mag"] = u(rng, 30e-9, 50e-9)
    p["mag_bias"] = [u(rng, -500e-9, 500e-9) for _ in range(3)]          # residual hard iron
    p["gyro_bias0"] = [u(rng, -5e-3, 5e-3) for _ in range(3)]            # NOTE: hand-off threshold issue (#2)
    gscale = u(rng, 1.0, 3.0)
    p["sigma_gyro"] = [3.2e-5 * gscale, 3.2e-5 * gscale, 4.8e-5 * gscale]  # LPF1 13 Hz nominal .. 3x
    p["sigma_rrw"] = u(rng, 1.4e-7, 5e-7)
    p["sigma_pd"] = u(rng, 0.01, 0.02)
    p["pd_spread"] = u(rng, 0.0, 0.05)
    p["albedo"] = u(rng, 0.0, 0.4)

    # actuators (ADCS_Actuators block nominal .. worst)
    p["sf_err"] = [u(rng, -0.2, 0.2) for _ in range(3)]
    p["m_res"] = [u(rng, -15e-3, 15e-3) for _ in range(3)]
    p["sigma_D"] = u(rng, 1e-3, 5.5e-3)
    p["coupling"] = u(rng, 0.0, 2.0)          # uT/A coil -> magnetometer. ASSUMPTION: hardware said "small"

    # ground segment
    p["uplink_days"] = u(rng, 1.0, 3.0)
    p["sigma_at"] = u(rng, 0.07, 0.27)
    return p


def to_args(p, sim_seed, orbits):
    args = [f"orbits={orbits}", f"seed={sim_seed}", "csv=none", "quiet=1", "worst=0"]
    for k, v in p.items():
        if isinstance(v, list):
            args.append(f"{k}=" + ",".join(f"{x:.9g}" for x in v))
        else:
            args.append(f"{k}={v:.9g}")
    return args


def flatten(p):
    out = {}
    for k, v in p.items():
        if isinstance(v, list):
            for j, x in enumerate(v):
                out[f"{k}_{j}"] = x
        else:
            out[k] = v
    return out


def run_one(sim, i, base_seed, orbits, extra):
    p = draw(i, base_seed)
    args = [sim] + to_args(p, i + 1, orbits) + extra
    res = subprocess.run(args, capture_output=True, text=True)
    summary = {}
    for line in res.stdout.splitlines():
        if line.startswith("SUMMARY"):
            for tok in line.split()[1:]:
                k, _, v = tok.partition("=")
                summary[k] = v
    row = {"run": i, "sim_seed": i + 1, "returncode": res.returncode}
    row.update(flatten(p))
    row.update(summary)
    if res.returncode != 0 or not summary:
        row["error"] = (res.stderr or "no SUMMARY").strip().replace("\n", " ")[:200]
    return row


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n", type=int, default=200)
    ap.add_argument("--start", type=int, default=0, help="first run index (to extend a campaign)")
    ap.add_argument("--seed", type=int, default=2026, help="campaign seed (parameter draws)")
    ap.add_argument("--orbits", type=float, default=6.0)
    ap.add_argument("--jobs", type=int, default=os.cpu_count() or 4)
    ap.add_argument("--sim", default=DEFAULT_SIM)
    ap.add_argument("--out", default="mc_runs.csv")
    ap.add_argument("extra", nargs="*", help="extra key=value passed to every run (e.g. a fault)")
    a = ap.parse_args()
    if not os.path.exists(a.sim):
        sys.exit(f"sim not found: {a.sim} (build it: cmake --build build)")

    rows = []
    with ThreadPoolExecutor(max_workers=a.jobs) as ex:
        futs = [ex.submit(run_one, a.sim, i, a.seed, a.orbits, a.extra) for i in range(a.start, a.start + a.n)]
        for k, f in enumerate(as_completed(futs), 1):
            rows.append(f.result())
            if k % max(1, a.n // 20) == 0 or k == a.n:
                print(f"  {k}/{a.n} runs done", flush=True)
    rows.sort(key=lambda r: r["run"])
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(a.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {len(rows)} runs to {a.out}")


if __name__ == "__main__":
    main()
