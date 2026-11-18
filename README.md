# ADCS Unix Testing

## Build 
```bash
git clone https://github.com/BrownSpaceEngineering/adcs-unix-testing.git
cd adcs-unix-testing
git submodule update --init --recursive
cmake -S . -B build
cmake --build build -j4
```

## Run the tests
```bash
./build/adcs-test
```
Every check prints `[PASS]` or `[FAIL]` (numeric checks also print the measured value and the
limit), and a summary at the end. A failing check never stops the run; the exit status is 1 if
anything failed. Tests live in `src/test.c`.

### Known failing tests
12 checks currently fail (357/369 pass). All 12 already failed before the sensor-frame change, on
`c516e7f`, with identical values, so they are not caused by it. Not yet fixed:

- `photodiode sun vector: NaN reading -> rejected`: a NaN photodiode reading is not rejected.
  This may be a real bug in `photodiode_determination.c` and needs a look.
- 11 checks in `ADCS step (body.c)`. These look like stale tests rather than an algorithm problem:
  `body()` now waits `BODY_READY_TO_POINT_TICKS` (50) consecutive ready ticks before handing off
  from detumble to pointing, but the tests expect the hand-off after a single step.
  - Hand-off after one step: `low rates + sun visible -> switches to pointing`,
    `attitude is initialized`, `QUEST initial attitude error` (91 deg, measured against an attitude
    QUEST never set), `transition step commands zero moment`.
  - Tick bookkeeping that assumes pointing started at tick 0: `ticks n through K command zero
    magnetic moment`, `predict-only ticks leave gyro bias unchanged`, `only measurement ticks make
    a meaningful gyro-bias correction`.
  - Recovery expected within a few steps: `valid magnetometer after a persistent fault resumes
    pointing`, and the three `gyro glitch:` checks (no filter is running to reset, and 5 steps is
    fewer than the 50-tick wait).

  The pointing-accuracy checks still pass: sunlit error 0.37 deg, eclipse max 1.0 deg, gyro bias
  within 0.002 deg/s. The plant sim also behaves normally.

On Linux (glibc), strict C99 hides `M_PI` / `M_SQRT1_2`, so `test.c` and `plant_sim.c` don't
compile as-is. For now, configure with `cmake -S . -B build -DCMAKE_C_FLAGS=-D_DEFAULT_SOURCE`.

Some expected values come from `tools/reference_values.py`, a pure-Python (no numpy)
double-precision reference that re-derives them independently (e.g. the WMM field as
`-grad(V)` of the spherical-harmonic potential). Re-run it if you change a model:
```bash
python3 tools/reference_values.py
```

## Conventions
- Quaternions are Hamilton, scalar-first `[w, x, y, z]`, active (`quat_apply(q, v) = q v q^-1`).
  The filter estimates body->ECI.
- `body()` takes **raw sensor-axis** readings and maps them to the body frame itself
  (`include/sensor_frames.h`): the magnetometer as `[U8, U9, U10]` in **Tesla** (same units as the
  WMM model), the gyro in its own chip axes in rad/s. Output currents are body X/Y/Z.
  Orbital elements are `[a (m), e, i, RAAN, argp, nu]` with angles in **degrees**.
  `orbital_to_eci*` take angles in radians.

  | Body axis | Magnetometer | Gyro | Magnetorquer |
  |---|---|---|---|
  | X | U8 | −Y | C2 |
  | Y | U9 | −X | C1 |
  | Z | −U10 | −Z | Air |

  The S-band antenna is on +Z (the pointing axis in `newpointing.c`); photodiode normals are
  already in body axes.
- `iterate()` (the UKF) takes 0 (predict only), 1 (magnetometer only, eclipse) or 2 (magnetometer
  and sun) measurement vectors, which should be unit vectors.

## C vs MATLAB parity tests
The `MATLAB parity` suite in `src/test.c` replays inputs through the C code and compares against
outputs recorded from the team's MATLAB (`../adcs/MATLAB/Algorithms`) in
`include/matlab_reference.h`. The header is committed, so `adcs-test` never needs MATLAB. To
regenerate it after changing the MATLAB (runs MATLAB headless, then rebuilds and runs the tests):
```bash
tools/matlab/run_matlab_parity.sh            # or: tools/matlab/run_matlab_parity.sh /path/to/adcs
```
- Generation is deterministic (fixed seed), so the header only changes when the MATLAB does.
- Known MATLAB bugs are printed as `[NOTE]` lines, not failures. For those, the C code is compared
  with the MATLAB after one exact, named text fix applied to a temporary copy
  (`generate_matlab_reference.m`); nothing in `../adcs` is modified. If upstream MATLAB changes so
  a fix no longer applies, the header records it and the notes say so.
- `pointing_error.m` needs `quatmultiply`/`quatinv`/`rotm2quat` (Aerospace/Robotics toolboxes).
  If they're missing, the plain implementations in `tools/matlab/shims` are used.
- The UKF comparison uses the local functions of `simulink.m` directly (not `fcn`, which needs the
  Statistics Toolbox for `normrnd`).
