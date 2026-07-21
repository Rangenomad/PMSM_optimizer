# MTPA Calibration & Efficiency MAP — Design Spec

> Date: 2026-07-21
> Branch: `subflow-update`
> Scope: Refactor Sub-flows A/D/E to use ψd/ψq flux-linkage table + MTPA calibration
> Sub-flow B: not modified in this branch (stable)
> Sub-flow C: removed — rated torque is superseded by Phase 2 MTPA calibration

---

## Architecture Overview

```
Step 0: θr Calibration (test_theta_r.py)
           │
           ▼
Phase 1: ψd/ψq 2D Flux Table (MS, Id×Iq grid)
           │  psi_dq_table.npz
           ▼
Phase 2: MTPA Calibration
           │  Input: (n, Is)   Output: (Thet_opt, Id_opt, Iq_opt, T_mtpa, Vs, M)
           │  mtpa_table.npz
           ├────► Phase 3-E: External Characteristic
           │         Input: Vdc, Imax (user-configurable)
           │         Output: external_characteristic.csv (T-n curve)
           │
           └────► Phase 3-D: Efficiency MAP
                     Input: mtpa_table → Thet at each (n, Is)
                     FEA: TR solver on speed×Is grid
                     Post: resample to speed×torque grid
                     Output: 6 MAP CSVs
```

**Key insight:** Phase 2 produces an (n,Is)→Thet table once. Phase 3-E reads feasibility boundaries for any Vdc/Imax. Phase 3-D runs FEA only with pre-computed MTPA angles — no voltage divergence.

---

## Step 0: θr Calibration

**Script:** `scripts/test_theta_r.py` (standalone, runs before all other steps)

**Purpose:** Find the permanent-magnet d-axis offset angle θr relative to phase-A winding axis.

**Method:**
1. Open `4_Partial_motor_MS2` (Magnetostatic)
2. Set all 6 coil currents to 0 A (zero-current excitation)
3. Solve — only PM flux present
4. `export_matrix()` → phase flux linkages ψa, ψb, ψc
5. Sweep θ ∈ [−90°, 90°], step 0.5°:
   - Park transform (ψa, ψb, ψc, θ) → ψd(θ), ψq(θ)
6. Find θr = argmin |ψq(θ)| (ideal: ψq=0, ψd=Φ)

**Validation:**
- ψd(θr) ≈ 0.099 Wb (PM flux linkage)
- ψq(θr) ≈ 0

**Output:** prints θr (deg), ψd(θr), ψq(θr) to stdout.

**Fallback:** If validation fails, user adjusts template so d-axis aligns with A-axis (θr=0).

---

## Phase 1: ψd/ψq 2D Flux Linkage Table

**Script:** integrated into `scripts/subflow_a_mtpa_cal.py`

**Solver:** MagnetostaticXY (`4_Partial_motor_MS2`)

**Verification grid (10×10):**

| Parameter | Range | Points | Step |
|-----------|-------|--------|------|
| Id        | −400 ~ +400 A | 10 | ~89 A |
| Iq        |  0 ~ +400 A   | 10 | ~44 A |
| Total     |               | 100  | ~5-8 min |

**Per-point procedure:**
```
1. Ia, Ib, Ic = id_iq_to_abc(Id, Iq, θr)    # use calibrated θr
2. _set_coil_currents(m2d, boundaries, Ia, Ib, Ic)
3. m2d.analyze('Setup1')
4. ψa, ψb, ψc = _extract_flux_linkages(m2d)  # via export_matrix()
5. ψd, ψq = park_abc_to_dq(ψa, ψb, ψc, θr)
6. Store in psi_dq[i_id, i_iq]
```

**Output file: `psi_dq_table.npz`**
```python
np.savez(path,
    id_grid  = np.array([-400, -311, ..., 400]),  # (n_id,)
    iq_grid  = np.array([0, 44, ..., 400]),        # (n_iq,)
    psi_d    = np.array([...]),                     # (n_id, n_iq) Wb
    psi_q    = np.array([...]),                     # (n_id, n_iq) Wb
    theta_r  = θr)                                  # float (rad)
```

**Key differences from old subflow_a:**
- No Ld/Lq calculation (no secant compensation, no ψd(0,Iq) interpolation)
- Raw ψd/ψq stored directly — cross-saturation naturally captured in 2D table
- θr injected from Step 0, not self-measured

---

## Phase 2: MTPA Calibration

**Script:** integrated into `scripts/subflow_a_mtpa_cal.py`

**Method:** Pure mathematical search on the ψd/ψq table (no FEA).

**Input:** (n [rpm], Is [A])  
**Output:** (Thet_opt, Id_opt, Iq_opt, T_mtpa, Vs, M, feasible)

### Core Equations

For any (Id, Iq):
```
ψd, ψq  = bilinear_interpolate(psi_dq_table, Id, Iq)
T       = 1.5 × P × (ψd·Iq − ψq·Id)
ω       = n · π/30 · P
Vd      = Rs·Id − ω·ψq
Vq      = Rs·Iq + ω·ψd
Vs      = √(Vd² + Vq²)
Vmax    = Vdc / √3
M       = Vs / Vmax
```

### Search Algorithm

For each (n, Is):
```
Given Is, sweep Thet ∈ [−90°, 90°] (fine grid, e.g. 200 points)
  Id = Is·cos(Thet), Iq = Is·sin(Thet)
  Compute T, Vs, M

Select Thet that maximizes T subject to:
  (1) Vs ≤ Vmax (voltage constraint)
  (2) Is ≤ Imax (current constraint, implicit)

If no Thet satisfies voltage constraint → feasible = False
If no Thet produces T > 0 → feasible = False
```

This is **maximum torque at given Is** constrained by voltage — the classic MTPA definition.

### Output Grid

| Parameter | Range | Points |
|-----------|-------|--------|
| n (speed) | 500 ~ 8000 rpm | 11 |
| Is (current) | 0 ~ Imax (e.g. 350 A) | 21 |

### Output File: `mtpa_table.npz`

```python
np.savez(path,
    speed_grid  = np.array([500, 1250, 2000, ..., 8000]),  # (n_spd,) rpm
    Is_grid     = np.array([0, 17.5, 35, ..., 350]),       # (n_Is,) A
    Thet_opt    = np.array([...]),  # (n_spd, n_Is) deg
    T_mtpa      = np.array([...]),  # (n_spd, n_Is) Nm — torque from flux table
    Vs          = np.array([...]),  # (n_spd, n_Is) V
    M           = np.array([...]),  # (n_spd, n_Is) dimensionless
    feasible    = np.array([...]),  # (n_spd, n_Is) bool
)
```

---

## Phase 3-E: External Characteristic (T-n Curve)

**Script:** refactored `scripts/subflow_e_external.py`

**Input:**
- `mtpa_table.npz` (Phase 2 output)
- User-configurable: `Vdc`, `Imax`, `Rs`

**Method:** Pure table lookup (no FEA).

For each speed n:
```
Vmax = Vdc / √3
Rank all feasible points at speed n by T_mtpa descending
Select the one with max T_mtpa satisfying:
  Is ≤ Imax (user's new current limit)
  Vs ≤ Vmax (user's new voltage limit)
```

The same `mtpa_table.npz` works for any Vdc/Imax combination — no re-calibration needed.

**Output: `external_characteristic.csv`**
```csv
n_rpm, T_max_Nm, Id_A, Iq_A, Thet_deg, Is_A, Vs_V, M, limit_type
500,   248,      -30,  248,  97,      250,   45,  0.17,  current
4000,  180,      -120, 215,  119,     246,   260, 1.00,  voltage
...
```

`limit_type`: `current` = Imax-constrained, `voltage` = Vdc-constrained.

---

## Phase 3-D: Efficiency MAP

**Script:** refactored `scripts/subflow_d_efficiency_map.py`

**Solver:** TransientXY (`5_Partial_motor_TR`)

### Step 1: FEA Scan on Speed×Is Grid

For each (n, Is) where `mtpa_table.feasible[n, Is] = True`:
```
1. Look up Thet_opt from mtpa_table
2. Set FEA variables: Speed_rpm, Imax=Is, Thet_deg=Thet_opt
3. TR solve (1 elec period × 50 steps)
4. Extract Moving1.Torque → T_actual
5. Store: (n, Is, Thet_opt, T_actual)
```

Feasible=False points are skipped entirely (solver would diverge).

### Step 2: Resample to Speed×Torque Uniform Grid

```
For each speed n:
  Collect (T_actual, Is) data points from FEA
  Filter: Is > 1A and T_actual > 0.5 Nm (exclude failed solves)
  Sort by T_actual ascending
  Interpolate to uniform T_grid
  Compute analytical quantities via ψd/ψq table interpolation:
    - Vs, M, PF, P_cu, efficiency
```

### Output: 6 MAP CSVs (speed×torque coordinates)

```
efficiency_map.csv    # η (%)
torque_map.csv        # T_actual (Nm, FEA measured)
current_map.csv       # Is (A)
modulation_map.csv    # M (dimensionless)
pf_map.csv            # Power Factor
voltage_map.csv       # Vs (V)
```

Format: `row=n_rpm, col=T_Nm, comma-separated values`
(Matches existing output format for backward compatibility.)

---

## File Structure

```
scripts/
├── test_theta_r.py              (NEW — Step 0, standalone)
├── subflow_a_mtpa_cal.py        (RENAMED from subflow_a_ldlq.py)
│   ├── Phase 1: psi_dq_scan()
│   └── Phase 2: mtpa_calibrate()
├── subflow_b_bemf.py            (UNCHANGED)
├── subflow_c_torque.py          (REMOVED — superseded by Phase 2 MTPA)
├── subflow_d_efficiency_map.py  (REFACTORED — uses mtpa_table)
└── subflow_e_external.py        (REFACTORED — uses mtpa_table)

pmsm_projects/<date>_mtpa_cal/
├── Prius_2D_Practice.aedt       (template copy)
├── psi_dq_table.npz             (Phase 1 output)
├── mtpa_table.npz               (Phase 2 output)
├── external_characteristic.csv  (Phase 3-E output)
└── efficiency_map/              (Phase 3-D output)
    ├── efficiency_map.csv
    ├── torque_map.csv
    ├── current_map.csv
    ├── modulation_map.csv
    ├── pf_map.csv
    └── voltage_map.csv
```

---

## Error Handling

- **gRPC retry:** 3 retries with exponential backoff (2^attempt seconds) for all `m2d.analyze()` calls
- **Checkpoint resume:** Save after each FEA point in Phase 3-D; resume from last completed point on restart
- **Feasibility guard:** Phase 3-D only submits FEA jobs for feasible (n, Is) points — never points where Vs > Vmax
- **Interpolation boundary:** Bilinear interpolation in psi_dq_table falls back to nearest-neighbor at grid edges
- **θr validation:** Step 0 fails fast if ψq(θr) / ψd(θr) > 0.05 (d-axis alignment suspect)

---

## Verification Criteria

| Step | Pass Condition |
|------|---------------|
| θr test | ψq(θr) ≈ 0, ψd(θr) ≈ 0.099 Wb |
| Phase 1 (sparse) | ψd/ψq table shows expected trends (ψq grows with Iq, ψd decreases with −Id) |
| Phase 2 (MTPA) | T_mtpa values are physically reasonable; feasible boundary matches expected speed range |
| Phase 3-E | T-n curve shape matches motor physics: flat torque at low speed, hyperbolic decay at high speed |
| Phase 3-D | FEA T_actual matches T_mtpa within ~10%; efficiency MAP shows peak > 95% |
