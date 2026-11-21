# Monte Carlo (host closed loop)

`run_mc.py` runs `build/adcs-plant-sim` (plant.c truth + the real `body()`) many times with every
uncertain truth parameter drawn at random, and writes one CSV row per run (inputs + the sim's
`SUMMARY` metrics). `report_mc.py` turns that CSV into pass rates, percentiles, sensitivities and
the worst runs. Standard-library Python only.

```bash
cmake --build build
python3 tools/mc/run_mc.py --n 1000 --jobs 10 --orbits 6 --out mc_runs.csv
python3 tools/mc/report_mc.py mc_runs.csv --p95-limit 10 --handoff-limit-orbits 3
# add a fault to every run, e.g. a reversed x coil from launch:
python3 tools/mc/run_mc.py --n 200 --out mc_flip.csv coil_fault_t=0 coil_gain=-1,1,1
# re-run run 137 exactly (same campaign --seed):
python3 tools/mc/run_mc.py --start 137 --n 1 --out one.csv
#   ...and also save its 10 s time history (extra key=value args override the defaults):
python3 tools/mc/run_mc.py --start 137 --n 1 --out one.csv csv=run137.csv
```

A 6-orbit run takes ~10 s, so 1000 runs on 10 cores is ~15-20 min. Simulink is ~20x real time,
so use it to cross-check the worst cases, not for the sweep.

## Parameter ranges (ASSUMPTIONS - replace with spec / measured values)

| Parameter | Range | Source |
|---|---|---|
| Inertia | principal moments x U(0.8, 1.2), principal axes rotated up to 10 deg | guess: CAD vs flight |
| Altitude | U(450, 650) km | guess: rideshare LEO |
| Inclination | U(45, 98) deg | guess |
| RAAN, arg. perigee, mean anomaly | U(0, 360) deg | - |
| Epoch (sun geometry, eclipses) | U(0, 365) days after 2025-01-01 | - |
| Tip-off rate | \|w0\| U(0, 30) deg/s, random direction | guess |
| Initial attitude | uniform random | - |
| Magnetometer noise | U(30, 50) nT | Simulink nominal .. worst |
| Magnetometer hard iron | U(-500, 500) nT per axis | Simulink worst |
| Gyro turn-on bias | U(-5, 5) mrad/s per axis | Simulink worst (see hand-off threshold) |
| Gyro white noise | 1x..3x [3.2, 3.2, 4.8]e-5 rad/s | SCH16T LPF1 13 Hz, 10 Hz sampling |
| Gyro rate random walk | U(1.4e-7, 5e-7) rad/s/sqrt(s) | Simulink nominal .. worst |
| Photodiode noise / gain spread / albedo | U(0.01, 0.02) / U(0, 5%) / U(0, 0.4) | Simulink nominal .. worst |
| Coil scale-factor error | U(-20, 20)% per axis | Simulink worst |
| Residual dipole | U(-15, 15) mA m^2 per axis | Simulink worst |
| Coil current noise | U(1, 5.5) mA m^2 | Simulink nominal .. worst |
| Coil -> magnetometer coupling | U(0, 2) uT/A | guess ("small" in hardware test) |
| TLE uplink period / OD timing error | U(1, 3) days / U(0.07, 0.27) s | Simulink nominal .. worst |

## Metrics (SUMMARY line of adcs-plant-sim)

- `handoff_s`: first switch from detumbling to pointing (NaN = never).
- `acquire_s`: time from that hand-off until nadir error stays < 10 deg for 600 s.
- Settled window = from one orbit after the hand-off to the end: `point_frac_settled`, `nadir_rms`,
  `nadir_p95`, `nadir_max` (truth +z vs nadir), `est_rms`, `est_p95` (filter vs truth, while pointing).
- `resets`, the sensor-check reject counters, `rej_gated` (innovation gate), `magcal_*`, `coil_*`.
