# MTPA Calibration & Efficiency MAP — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace Ld/Lq indirect model with direct ψd/ψq flux-linkage table → MTPA calibration → refactored efficiency MAP and external characteristic.

**Architecture:** Three-stage pipeline: (1) MS magnetostatic Id×Iq scan → ψd/ψq 2D table, (2) pure-math MTPA search on the table → (n,Is)→Thet lookup, (3) downstream E (T-n curve) and D (FEA TR efficiency MAP) consume the MTPA table.

**Tech Stack:** PyAEDT 0.25.1 (gRPC), Ansys Maxwell 2D 2023.2 (MagnetostaticXY + TransientXY), Python 3.13, NumPy, SciPy.

## Global Constraints

- All scripts use `sys.stdout.reconfigure(encoding='utf-8')` as first line
- All scripts open with `non_graphical=False, new_desktop=False, close_on_exit=False`
- Template: `references/Prius_2D_Practice.aedt` → copied to project dir
- Winding formula: `Imax*sin(Omega_rad*time + Thet)`, 9× turns factor per coil
- Thet_deg assigned as `str(angle)` without `°` suffix (gRPC second-assignment silent failure)
- gRPC: `boundary.update()` for MS current setting; `m2d[...]=value` for design variables
- Project path interface: `project_path=None` defaults to references template
- Phases 1 & 2 are combined in `subflow_a_mtpa_cal.py`; Phases 3-E & 3-D are separate scripts
- Sub-flow C (`subflow_c_torque.py`) removed — MTPA table supersedes it
- Sub-flow B unchanged

---

### Task 1: Write θr Calibration Test Script

**Files:**
- Create: `scripts/test_theta_r.py`

**Interfaces:**
- Consumes: `references/Prius_2D_Practice.aedt` (design `4_Partial_motor_MS2`)
- Produces: prints θr (deg), ψd(θr), ψq(θr) to stdout

- [ ] **Step 1: Create `scripts/test_theta_r.py`**

```python
"""θr Calibration — find PM d-axis offset angle.

Zero-current Magnetostatic solve → extract ψa/ψb/ψc → sweep Park angle
to find θ where ψq=0 (d-axis aligned with rotor PM field).

Usage:
    python scripts/test_theta_r.py
    python scripts/test_theta_r.py project_path=pmsm_projects/2026-07-19_efficiency_map
"""

import sys
sys.stdout.reconfigure(encoding='utf-8')
import os, sys, tempfile
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
_DEFAULT_TEMPLATE = str(ROOT / 'references' / 'Prius_2D_Practice.aedt')

# -- reused from subflow_a_ldlq.py (copied inline for standalone use) --

def park_abc_to_dq(psi_a, psi_b, psi_c, theta=0):
    """Park transform: abc → dq. q-axis LAGS d-axis by 90°."""
    c = np.cos(theta)
    s = np.sin(theta)
    psi_d = 2/3 * (psi_a * c + psi_b * np.cos(theta - 2*np.pi/3) + psi_c * np.cos(theta + 2*np.pi/3))
    psi_q = 2/3 * (psi_a * s + psi_b * np.sin(theta - 2*np.pi/3) + psi_c * np.sin(theta + 2*np.pi/3))
    return psi_d, psi_q


def main(project_path=None):
    from ansys.aedt.core import Maxwell2d

    if project_path:
        template = str(Path(project_path) / 'Prius_2D_Practice.aedt')
    else:
        template = _DEFAULT_TEMPLATE

    print(f'\n{"=" * 60}')
    print(f'  θr Calibration Test')
    print(f'{"=" * 60}')
    print(f'  Template: {template}')
    print(f'  Design:   4_Partial_motor_MS2 (Magnetostatic)')

    m2d = Maxwell2d(
        project=template, design='4_Partial_motor_MS2',
        solution_type='MagnetostaticXY',
        non_graphical=False, new_desktop=False, close_on_exit=False
    )

    # Phase coil boundaries (6 coils)
    coil_names = ['PhaseA1', 'PhaseA2', 'PhaseB1', 'PhaseB2', 'PhaseC1', 'PhaseC2']
    boundaries = []
    for name in coil_names:
        b = [b for b in m2d.boundaries if b.name == name][0]
        boundaries.append(b)

    # Set all coil currents to 0
    TURNS = 9
    for bnd in boundaries:
        bnd.props['Current'] = '0A'
        bnd.update()

    # Solve
    print('  Solving zero-current magnetostatic...')
    m2d.analyze('Setup1')

    # Export matrix → extract flux linkages
    out_file = os.path.join(tempfile.gettempdir(), 'pmsm_theta_r_test.txt')
    m2d.export_matrix('Matrix1', out_file)
    with open(out_file, 'r') as f:
        text = f.read()

    # Parse flux linkages per coil
    result = {}
    in_flux = False
    for line in text.splitlines():
        if line.strip().startswith('Flux Linkage'):
            in_flux = True
            continue
        if in_flux and line.strip() and not line.strip().startswith('Flux'):
            parts = line.strip().split()
            if len(parts) >= 2 and parts[0] in coil_names:
                result[parts[0]] = float(parts[1])

    # Total phase flux = 9 × (coil1 + coil2)
    psi_a = TURNS * (result['PhaseA1'] + result['PhaseA2'])
    psi_b = TURNS * (result['PhaseB1'] + result['PhaseB2'])
    psi_c = TURNS * (result['PhaseC1'] + result['PhaseC2'])

    print(f'  ψa = {psi_a:.6f} Wb')
    print(f'  ψb = {psi_b:.6f} Wb')
    print(f'  ψc = {psi_c:.6f} Wb')

    # Sweep θ to find d-axis (ψq minimized)
    theta_deg = np.linspace(-90, 90, 361)  # 0.5° steps
    psi_q_vals = []
    psi_d_vals = []
    for th in theta_deg:
        pd, pq = park_abc_to_dq(psi_a, psi_b, psi_c, theta=np.radians(th))
        psi_d_vals.append(pd)
        psi_q_vals.append(pq)

    psi_q_vals = np.array(psi_q_vals)
    psi_d_vals = np.array(psi_d_vals)

    idx_min = np.argmin(np.abs(psi_q_vals))
    theta_r_deg = theta_deg[idx_min]
    theta_r_rad = np.radians(theta_r_deg)
    psi_d_at_theta_r = psi_d_vals[idx_min]
    psi_q_at_theta_r = psi_q_vals[idx_min]

    print(f'\n  {"=" * 40}')
    print(f'  RESULTS')
    print(f'  {"=" * 40}')
    print(f'  θr      = {theta_r_deg:.1f}° ({theta_r_rad:.4f} rad)')
    print(f'  ψd(θr)  = {psi_d_at_theta_r:.6f} Wb  (expect ≈ 0.099 Wb)')
    print(f'  ψq(θr)  = {psi_q_at_theta_r:.6f} Wb  (expect ≈ 0)')

    # Validation
    qd_ratio = abs(psi_q_at_theta_r) / max(abs(psi_d_at_theta_r), 1e-12)
    if qd_ratio < 0.05:
        print(f'\n  ✓ PASS: |ψq/ψd| = {qd_ratio:.4f} < 0.05')
    else:
        print(f'\n  ✗ FAIL: |ψq/ψd| = {qd_ratio:.4f} ≥ 0.05')
        print(f'    d-axis alignment suspect. Consider adjusting template.')

    m2d.close_project()
    return theta_r_rad, psi_d_at_theta_r, psi_q_at_theta_r


if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument('project_path', nargs='?', default=None)
    args = p.parse_args()
    main(project_path=args.project_path)
```

- [ ] **Step 2: Run test to verify**

Run: `python scripts/test_theta_r.py`
Expected: prints θr ≈ 0° (if d-axis aligned) or some offset, ψq≈0

- [ ] **Step 3: Commit**

```bash
git add scripts/test_theta_r.py
git commit -m "feat: add θr calibration test script"
```

---

### Task 2: Rename + Refactor Sub-flow A — Phase 1 ψd/ψq Table

**Files:**
- Rename: `scripts/subflow_a_ldlq.py` → `scripts/subflow_a_mtpa_cal.py`
- Modify: `scripts/param_guard.py` (add 'A' entry for new variables)

**Interfaces:**
- Produces: `phase1_psi_dq_table(id_range, iq_range, n_id, n_iq, theta_r, project_path)` → (grids, psi_d, psi_q) + saves `psi_dq_table.npz`

- [ ] **Step 1: Copy existing file under new name**

```powershell
Copy-Item scripts/subflow_a_ldlq.py scripts/subflow_a_mtpa_cal.py
```

- [ ] **Step 2: Remove old Ld/Lq logic, replace with ψd/ψq table builder**

Edit `scripts/subflow_a_mtpa_cal.py` — replace entire content after the helper functions (keep `id_iq_to_abc`, `park_abc_to_dq`, `_set_coil_currents`, `_parse_matrix_flux_linkage`, `_extract_flux_linkages`, constants `TURNS`, `COIL_CONFIG`, `COIL_NAMES`) with new Phase 1 + Phase 2 code:

```python
"""Sub-flow A: MTPA Calibration
Phase 1: ψd/ψq 2D flux-linkage table (Magnetostatic, Id×Iq grid)
Phase 2: MTPA calibration → (n, Is) → (Thet_opt, Id_opt, T_mtpa, Vs, M)

Usage:
    python -c "from scripts.subflow_a_mtpa_cal import run; run(n_id=10, n_iq=10, project_path='pmsm_projects/xxx')"
"""

import sys, os, re, tempfile
sys.stdout.reconfigure(encoding='utf-8')
from pathlib import Path
import numpy as np
from ansys.aedt.core import Maxwell2d

ROOT = Path(__file__).resolve().parent.parent
_DEFAULT_TEMPLATE = str(ROOT / 'references' / 'Prius_2D_Practice.aedt')

TURNS = 9
COIL_CONFIG = [
    ('PhaseA1', True), ('PhaseA2', True),
    ('PhaseB1', True), ('PhaseB2', True),
    ('PhaseC1', True), ('PhaseC2', True),
]
COIL_NAMES = [c[0] for c in COIL_CONFIG]

# ── Coordinate transforms (reused from old subflow_a) ──

def id_iq_to_abc(Id, Iq, theta=0):
    """Inverse Park + Clarke: dq → abc."""
    Ialpha = Id * np.cos(theta) + Iq * np.sin(theta)
    Ibeta  = Id * np.sin(theta) - Iq * np.cos(theta)
    Ia = Ialpha
    Ib = -0.5 * Ialpha + np.sqrt(3) / 2 * Ibeta
    Ic = -0.5 * Ialpha - np.sqrt(3) / 2 * Ibeta
    return Ia, Ib, Ic


def park_abc_to_dq(psi_a, psi_b, psi_c, theta=0):
    """Park transform: abc → dq. q-axis LAGS d-axis by 90°."""
    c = np.cos(theta)
    s = np.sin(theta)
    psi_d = 2/3 * (psi_a*c + psi_b*np.cos(theta-2*np.pi/3) + psi_c*np.cos(theta+2*np.pi/3))
    psi_q = 2/3 * (psi_a*s + psi_b*np.sin(theta-2*np.pi/3) + psi_c*np.sin(theta+2*np.pi/3))
    return psi_d, psi_q


def _set_coil_currents(m2d, boundaries, Ia, Ib, Ic):
    """Set coil currents via boundary.update() (含 9 匝因子)."""
    currents = [TURNS * Ia, TURNS * Ia,
                TURNS * Ib, TURNS * Ib,
                TURNS * Ic, TURNS * Ic]
    for bnd, (name, is_pos), val in zip(boundaries, COIL_CONFIG, currents):
        bnd.props['Current'] = f'{val}A'
        bnd.props['IsPositive'] = is_pos
        bnd.update()


def _parse_matrix_flux_linkage(text):
    """Parse matrix export .txt → dict {coil_name: flux_Wb}."""
    in_flux = False
    result = {}
    for line in text.splitlines():
        if line.strip().startswith('Flux Linkage'):
            in_flux = True
            continue
        if in_flux and line.strip() and not line.strip().startswith('Flux'):
            parts = line.strip().split()
            if len(parts) >= 2 and parts[0] in COIL_NAMES:
                result[parts[0]] = float(parts[1])
    return result


def _extract_flux_linkages(m2d):
    """Solve → export matrix → parse flux linkage per coil."""
    m2d.analyze('Setup1')
    out_file = os.path.join(tempfile.gettempdir(), 'pmsm_matrix_export.txt')
    m2d.export_matrix('Matrix1', out_file)
    with open(out_file, 'r') as f:
        text = f.read()
    return _parse_matrix_flux_linkage(text)


# ═══════════════════════════════════════════
# Phase 1: ψd/ψq 2D Flux-Linkage Table
# ═══════════════════════════════════════════

def phase1_psi_dq_table(id_range=(-400, 400), iq_range=(0, 400),
                        n_id=10, n_iq=10, theta_r=0.0,
                        project_path=None, confirm=True):
    """Scan Id×Iq grid in Magnetostatic, build ψd/ψq 2D table.

    Parameters
    ----------
    id_range : tuple
        (Id_min, Id_max) in A
    iq_range : tuple
        (Iq_min, Iq_max) in A
    n_id : int
        Number of Id steps
    n_iq : int
        Number of Iq steps
    theta_r : float
        Rotor d-axis angle in radians (from θr calibration)
    project_path : str, optional
        Project directory path
    confirm : bool
        If True, prompt user before running

    Returns
    -------
    id_grid : np.ndarray (n_id,)
    iq_grid : np.ndarray (n_iq,)
    psi_d : np.ndarray (n_id, n_iq)
    psi_q : np.ndarray (n_id, n_iq)
    """
    from scripts.project_utils import get_template_path

    if project_path:
        template = get_template_path(Path(project_path))
    else:
        template = _DEFAULT_TEMPLATE

    id_grid = np.linspace(id_range[0], id_range[1], n_id)
    iq_grid = np.linspace(iq_range[0], iq_range[1], n_iq)
    total = n_id * n_iq

    print(f'\n{"=" * 60}')
    print(f'  Phase 1: ψd/ψq 2D Flux-Linkage Table')
    print(f'  Solver: MagnetostaticXY (4_Partial_motor_MS2)')
    print(f'{"=" * 60}')
    print(f'  Id grid: {id_range[0]} ~ {id_range[1]} A, {n_id} points')
    print(f'  Iq grid: {iq_range[0]} ~ {iq_range[1]} A, {n_iq} points')
    print(f'  Total:   {total} points, ~{total*3//60} min')
    print(f'  θr:      {np.degrees(theta_r):.1f}°')
    print(f'{"=" * 60}')

    if confirm:
        resp = input('\n  继续执行? [Y/n] ').strip().lower()
        if resp and resp != 'y':
            print('  已取消。')
            return None, None, None, None

    m2d = Maxwell2d(
        project=template, design='4_Partial_motor_MS2',
        solution_type='MagnetostaticXY',
        non_graphical=False, new_desktop=False, close_on_exit=False
    )

    boundaries = []
    for name, _ in COIL_CONFIG:
        b = [b for b in m2d.boundaries if b.name == name][0]
        boundaries.append(b)

    psi_d = np.full((n_id, n_iq), np.nan)
    psi_q = np.full((n_id, n_iq), np.nan)

    for i, Id_val in enumerate(id_grid):
        for j, Iq_val in enumerate(iq_grid):
            point = i * n_iq + j + 1
            Ia, Ib, Ic = id_iq_to_abc(Id_val, Iq_val, theta=theta_r)
            print(f'  [P1] {point}/{total} Id={Id_val:+7.1f}A Iq={Iq_val:+7.1f}A '
                  f'Ia={Ia:.1f}A ...', end='', flush=True)

            _set_coil_currents(m2d, boundaries, Ia, Ib, Ic)
            psi_coil = _extract_flux_linkages(m2d)

            psi_a = TURNS * (psi_coil['PhaseA1'] + psi_coil['PhaseA2'])
            psi_b = TURNS * (psi_coil['PhaseB1'] + psi_coil['PhaseB2'])
            psi_c = TURNS * (psi_coil['PhaseC1'] + psi_coil['PhaseC2'])

            pd, pq = park_abc_to_dq(psi_a, psi_b, psi_c, theta=theta_r)
            psi_d[i, j] = pd
            psi_q[i, j] = pq
            print(f' ψd={pd:.6f} ψq={pq:.6f}')

    m2d.close_project()

    # Save
    out_dir = Path(project_path) if project_path else Path.cwd()
    out_path = out_dir / 'psi_dq_table.npz'
    np.savez(out_path,
             id_grid=id_grid, iq_grid=iq_grid,
             psi_d=psi_d, psi_q=psi_q,
             theta_r=theta_r)
    print(f'\n  Saved: {out_path}')

    return id_grid, iq_grid, psi_d, psi_q


# ═══════════════════════════════════════════
# Phase 2: MTPA Calibration
# ═══════════════════════════════════════════

def _bilinear_interp(x_grid, y_grid, z_table, x, y):
    """Bilinear interpolation on a 2D regular grid.

    Parameters
    ----------
    x_grid : np.ndarray (nx,) — column axis (Id grid)
    y_grid : np.ndarray (ny,) — row axis (Iq grid)
    z_table : np.ndarray (nx, ny)
    x, y : float

    Returns
    -------
    float — interpolated value
    """
    from scipy.interpolate import RegularGridInterpolator
    interp = RegularGridInterpolator((x_grid, y_grid), z_table,
                                     bounds_error=False, fill_value=None)
    val = interp([[x, y]])[0]
    if np.isnan(val):
        # Fallback: nearest neighbor at edges
        ix = max(0, min(len(x_grid)-1, np.searchsorted(x_grid, x)))
        iy = max(0, min(len(y_grid)-1, np.searchsorted(y_grid, y)))
        return float(z_table[ix, iy])
    return float(val)


def phase2_mtpa_calibrate(psi_dq_path, vdc=450, imax=350,
                          speed_min=500, speed_max=8000, n_speed=11, n_Is=21,
                          Rs=0.05, pole_pairs=4):
    """MTPA calibration: search (n, Is) → optimal Thet/Id/Iq using ψd/ψq table.

    For each (speed, Is), sweeps current angle β ∈ [-90°, 90°],
    selects β that maximizes torque subject to voltage constraint Vs ≤ Vmax.

    Parameters
    ----------
    psi_dq_path : str or Path
        Path to psi_dq_table.npz (Phase 1 output)
    vdc : float
        DC bus voltage (V)
    imax : float
        Maximum phase current (A)
    speed_min, speed_max : float
        Speed range (rpm)
    n_speed : int
        Number of speed points
    n_Is : int
        Number of current points
    Rs : float
        Phase resistance (Ω)
    pole_pairs : int

    Returns
    -------
    dict with keys: speed_grid, Is_grid, Thet_opt, T_mtpa, Vs, M, feasible
    """
    import scipy
    Vmax = vdc / np.sqrt(3)

    # Load ψd/ψq table
    data = np.load(psi_dq_path)
    id_grid = data['id_grid']
    iq_grid = data['iq_grid']
    psi_d_tab = data['psi_d']
    psi_q_tab = data['psi_q']

    from scipy.interpolate import RegularGridInterpolator
    psi_d_interp = RegularGridInterpolator(
        (id_grid, iq_grid), psi_d_tab, bounds_error=False, fill_value=None)
    psi_q_interp = RegularGridInterpolator(
        (id_grid, iq_grid), psi_q_tab, bounds_error=False, fill_value=None)

    def interp_psi(Id, Iq):
        """Interpolate ψd, ψq at (Id, Iq). Falls back to NN at edges."""
        try:
            pd = float(psi_d_interp([[Id, Iq]])[0])
            pq = float(psi_q_interp([[Id, Iq]])[0])
        except Exception:
            ix = np.clip(np.searchsorted(id_grid, Id), 0, len(id_grid)-1)
            iy = np.clip(np.searchsorted(iq_grid, Iq), 0, len(iq_grid)-1)
            pd = float(psi_d_tab[ix, iy])
            pq = float(psi_q_tab[ix, iy])
        return pd, pq

    speed_grid = np.linspace(speed_min, speed_max, n_speed)
    Is_grid = np.linspace(0, imax, n_Is)

    Thet_opt = np.full((n_speed, n_Is), np.nan)
    T_mtpa   = np.full((n_speed, n_Is), np.nan)
    Vs       = np.full((n_speed, n_Is), np.nan)
    M        = np.full((n_speed, n_Is), np.nan)
    feasible = np.zeros((n_speed, n_Is), dtype=bool)

    N_BETA = 200  # angle resolution
    beta_grid = np.linspace(-np.pi/2, np.pi/2, N_BETA)

    print(f'\n{"=" * 60}')
    print(f'  Phase 2: MTPA Calibration')
    print(f'{"=" * 60}')
    print(f'  Vdc={vdc:.0f}V, Vmax={Vmax:.1f}V, Imax={imax:.0f}A')
    print(f'  Speed: {speed_min}~{speed_max} rpm, {n_speed} points')
    print(f'  Is:    0~{imax} A, {n_Is} points')
    print(f'  ψd/ψq table: {psi_dq_path}')
    print(f'{"=" * 60}')

    for i, n in enumerate(speed_grid):
        omega_e = n * np.pi / 30 * pole_pairs
        for j, Is_val in enumerate(Is_grid):
            if Is_val < 1e-6:
                Thet_opt[i, j] = 0.0
                T_mtpa[i, j] = 0.0
                Vs[i, j] = 0.0
                M[i, j] = 0.0
                feasible[i, j] = True
                continue

            T_best = -1e12
            best = None  # (beta, Id, Iq, Vs)

            for beta in beta_grid:
                Id = Is_val * np.sin(beta)
                Iq = Is_val * np.cos(beta)

                psi_d_val, psi_q_val = interp_psi(Id, Iq)

                T = 1.5 * pole_pairs * (psi_d_val * Iq - psi_q_val * Id)
                if T <= 0:
                    continue

                Vd = Rs * Id - omega_e * psi_q_val
                Vq = Rs * Iq + omega_e * psi_d_val
                Vs_val = np.sqrt(Vd**2 + Vq**2)

                if Vs_val <= Vmax and T > T_best:
                    T_best = T
                    best = (beta, Id, Iq, Vs_val)

            if best is not None:
                beta_opt, Id_opt, Iq_opt, Vs_val = best
                Thet_opt[i, j] = np.degrees(beta_opt)
                T_mtpa[i, j] = T_best
                Vs[i, j] = Vs_val
                M[i, j] = Vs_val / Vmax
                feasible[i, j] = True
            # else: stays NaN, feasible=False

        n_feas = np.sum(feasible[i, :])
        print(f'  n={n:.0f} rpm: {n_feas}/{n_Is} feasible')

    # Save
    out_dir = Path(psi_dq_path).parent
    out_path = out_dir / 'mtpa_table.npz'
    np.savez(out_path,
             speed_grid=speed_grid, Is_grid=Is_grid,
             Thet_opt=Thet_opt, T_mtpa=T_mtpa,
             Vs=Vs, M=M, feasible=feasible,
             vdc=vdc, imax=imax, Rs=Rs, pole_pairs=pole_pairs)
    print(f'\n  Saved: {out_path}')

    return {'speed_grid': speed_grid, 'Is_grid': Is_grid,
            'Thet_opt': Thet_opt, 'T_mtpa': T_mtpa,
            'Vs': Vs, 'M': M, 'feasible': feasible}


# ═══════════════════════════════════════════
# Combined run(): Phase 1 + Phase 2
# ═══════════════════════════════════════════

def run(id_range=(-400, 400), iq_range=(0, 400), n_id=10, n_iq=10,
        theta_r=0.0, vdc=450, imax=350,
        speed_min=500, speed_max=8000, n_speed=11, n_Is=21,
        Rs=0.05, pole_pairs=4, project_path=None):
    """Run full MTPA calibration pipeline: Phase 1 + Phase 2.

    Phase 1: MS ψd/ψq flux-linkage table
    Phase 2: MTPA search → (n,Is)→Thet lookup table
    """
    out_dir = Path(project_path) if project_path else Path.cwd()
    psi_dq_path = out_dir / 'psi_dq_table.npz'

    # Phase 1
    phase1_psi_dq_table(
        id_range=id_range, iq_range=iq_range,
        n_id=n_id, n_iq=n_iq, theta_r=theta_r,
        project_path=project_path, confirm=True
    )

    # Phase 2
    phase2_mtpa_calibrate(
        psi_dq_path=str(psi_dq_path), vdc=vdc, imax=imax,
        speed_min=speed_min, speed_max=speed_max,
        n_speed=n_speed, n_Is=n_Is,
        Rs=Rs, pole_pairs=pole_pairs
    )


if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument('--project-path', default=None)
    p.add_argument('--n-id', type=int, default=10)
    p.add_argument('--n-iq', type=int, default=10)
    p.add_argument('--theta-r', type=float, default=0.0)
    p.add_argument('--vdc', type=float, default=450)
    p.add_argument('--imax', type=float, default=350)
    args = p.parse_args()
    run(n_id=args.n_id, n_iq=args.n_iq, theta_r=args.theta_r,
        vdc=args.vdc, imax=args.imax, project_path=args.project_path)
```

- [ ] **Step 3: Update `scripts/param_guard.py` — add MTPA Calibration to allowed vars**

```python
# Replace the ALLOWED_VARS dict:
ALLOWED_VARS: dict[str, set[str]] = {
    'A': set(),  # MTPA Calibration — uses boundary.update() for MS, no m2d[...] access
    'B': {'Speed_rpm', 'Imax'},
    'D': {'Speed_rpm', 'Imax', 'Thet_deg'},
}
# Remove 'C' entry entirely
```

- [ ] **Step 4: Commit**

```bash
git add scripts/subflow_a_ldlq.py scripts/subflow_a_mtpa_cal.py scripts/param_guard.py
git commit -m "refactor: rename subflow_a → MTPA calibration (Phase 1+2)"
```

---

### Task 3: Remove Sub-flow C

**Files:**
- Delete: `scripts/subflow_c_torque.py`
- Modify: `scripts/param_guard.py` (already done in Task 2)
- Modify: `SKILL.md` (remove Sub-flow C references)

- [ ] **Step 1: Delete the file**

```powershell
Remove-Item scripts/subflow_c_torque.py
```

- [ ] **Step 2: Update SKILL.md — remove C references**

Read SKILL.md, find and remove:
- `| C (额定扭矩) | ... | Imax, Speed_rpm, Thet_deg |` from the allowed params table
- `Sub-flow C: 额定点扭矩` from the sub-flow mappings
- References to `subflow_c_torque`

- [ ] **Step 3: Commit**

```bash
git add scripts/subflow_c_torque.py SKILL.md
git commit -m "refactor: remove Sub-flow C (superseded by MTPA calibration)"
```

---

### Task 4: Refactor Sub-flow E — External Characteristic via MTPA Table

**Files:**
- Modify: `scripts/subflow_e_external.py`

**Interfaces:**
- Consumes: `mtpa_table.npz` (Phase 2 output)
- Produces: `external_characteristic.csv`

- [ ] **Step 1: Rewrite `scripts/subflow_e_external.py`**

```python
"""Sub-flow E: 外特性 T-n 曲线 (based on MTPA table)

Method: Pure table lookup (no FEA).
         Reads mtpa_table.npz (Phase 2), applies user Vdc/Imax constraints,
         outputs T-n external characteristic with boundary classification.

Usage:
    python -c "from scripts.subflow_e_external import run; run(vdc=450, imax=350)"
    python -c "from scripts.subflow_e_external import run; run(vdc=300, imax=250, mtpa_table_path='pmsm_projects/xxx/mtpa_table.npz')"
"""

import sys
sys.stdout.reconfigure(encoding='utf-8')
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent


def run(vdc=450, imax=350, speed_max=8000, speed_points_n=40,
        mtpa_table_path=None, pole_pairs=4, project_path=None):
    """Calculate external T-n characteristic from MTPA table.

    Parameters
    ----------
    vdc : float — DC bus voltage (V)
    imax : float — max phase current (A)
    speed_max : float — max speed (rpm)
    speed_points_n : int
    mtpa_table_path : str, optional — path to mtpa_table.npz
    pole_pairs : int
    project_path : str, optional — project dir for default table lookup
    """
    Vmax = vdc / np.sqrt(3)

    # Resolve mtpa_table path
    if mtpa_table_path:
        table_path = Path(mtpa_table_path)
    elif project_path:
        table_path = Path(project_path) / 'mtpa_table.npz'
    else:
        # Find latest project
        proj_dir = ROOT / 'pmsm_projects'
        candidates = sorted(proj_dir.glob('*/mtpa_table.npz'), reverse=True)
        if not candidates:
            raise FileNotFoundError(
                'No mtpa_table.npz found. Run Phase 2 MTPA calibration first.')
        table_path = candidates[0]

    data = np.load(table_path)
    speed_grid = data['speed_grid']
    Is_grid = data['Is_grid']
    T_mtpa = data['T_mtpa']
    Vs = data['Vs']
    feasible = data['feasible']

    P = pole_pairs

    print(f'\n{"=" * 70}')
    print(f'  Sub-flow E: 外特性 T-n 曲线 (MTPA table lookup)')
    print(f'{"=" * 70}')
    print(f'  MTPA table:  {table_path}')
    print(f'  Vdc={vdc:.0f}V  Vmax={Vmax:.1f}V  Imax={imax:.0f}A  P={P}')
    print(f'{"=" * 70}')

    # Generate output speeds — reuse table grid points where possible,
    # interpolate T_n curve elsewhere
    out_speeds = sorted(set(
        np.linspace(min(speed_grid), speed_max, speed_points_n)
    ))

    results = []
    prev_T = None

    for n in out_speeds:
        # For each speed, find feasible (Is, T) pairs satisfying new constraints
        # Use the nearest speed index in the table
        i_speed = np.argmin(np.abs(speed_grid - n))

        valid_Is = []
        valid_T = []
        valid_idx = []

        for j, Is_val in enumerate(Is_grid):
            if not feasible[i_speed, j]:
                continue
            if Is_val > imax:
                continue
            Vs_val = Vs[i_speed, j]
            Vmax_new = vdc / np.sqrt(3)  # Use user's Vdc for constraint
            if Vs_val > Vmax_new:
                continue
            valid_Is.append(Is_val)
            valid_T.append(T_mtpa[i_speed, j])
            valid_idx.append(j)

        if not valid_T:
            results.append({
                'speed': n, 'torque': 0.0, 'Id': 0.0, 'Iq': 0.0,
                'Is': 0.0, 'Thet_deg': 0.0, 'Vs': 0.0, 'M': 0.0,
                'limit_type': 'none'
            })
            continue

        # Pick max torque among feasible points
        best_j = valid_idx[np.argmax(valid_T)]
        T_max = float(T_mtpa[i_speed, best_j])
        Is_op = float(Is_grid[best_j])
        Thet_op = float(Thet_opt[i_speed, best_j]) if 'Thet_opt' in data else 0.0
        Vs_op = float(Vs[i_speed, best_j])
        Id_op = Is_op * np.sin(np.radians(Thet_op))
        Iq_op = Is_op * np.cos(np.radians(Thet_op))

        # Determine limit type
        if abs(Is_op - imax) < 1.0 or abs(Is_op - max(valid_Is)) < 1.0:
            limit_type = 'current'
        else:
            limit_type = 'voltage'

        results.append({
            'speed': n, 'torque': T_max,
            'Id': Id_op, 'Iq': Iq_op, 'Thet_deg': Thet_op,
            'Is': Is_op, 'Vs': Vs_op,
            'M': Vs_op / Vmax,
            'limit_type': limit_type,
        })

    # ── Output ──
    print(f'\n{"n_rpm":>8s} {"T_Nm":>8s} {"Id_A":>7s} {"Iq_A":>7s} {"Thet":>6s} {"Is_A":>7s} {"Vs_V":>7s} {"M":>6s} {"limit":>8s}')
    print('-' * 72)
    for r in results:
        print(f'{r["speed"]:8.0f} {r["torque"]:8.1f} {r["Id"]:7.1f} {r["Iq"]:7.1f} '
              f'{r["Thet_deg"]:6.1f} {r["Is"]:7.1f} {r["Vs"]:7.1f} {r["M"]:6.3f} {r["limit_type"]:>8s}')

    # Save CSV
    out_dir = Path(project_path) if project_path else table_path.parent
    out_path = out_dir / 'external_characteristic.csv'
    with open(out_path, 'w', encoding='utf-8') as f:
        f.write('n_rpm,T_max_Nm,Id_A,Iq_A,Thet_deg,Is_A,Vs_V,M,limit_type\n')
        for r in results:
            f.write(f'{r["speed"]:.1f},{r["torque"]:.2f},{r["Id"]:.2f},{r["Iq"]:.2f},'
                    f'{r["Thet_deg"]:.1f},{r["Is"]:.2f},{r["Vs"]:.1f},{r["M"]:.4f},{r["limit_type"]}\n')
    print(f'\n  Saved: {out_path}')

    # Stats
    valid = [r for r in results if r['torque'] > 0.1]
    if valid:
        peak = max(valid, key=lambda r: r['torque'])
        fw = [r for r in valid if r['limit_type'] == 'voltage']
        print(f'\n  峰值扭矩: {peak["torque"]:.1f} Nm @ {peak["speed"]:.0f} rpm')
        if fw:
            print(f'  基速/转折点: {fw[0]["speed"]:.0f} rpm')

    return results


if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument('--vdc', type=float, default=450)
    p.add_argument('--imax', type=float, default=350)
    p.add_argument('--mtpa-table', default=None)
    p.add_argument('--project-path', default=None)
    args = p.parse_args()
    run(vdc=args.vdc, imax=args.imax,
        mtpa_table_path=args.mtpa_table, project_path=args.project_path)
```

- [ ] **Step 2: Commit**

```bash
git add scripts/subflow_e_external.py
git commit -m "refactor: subflow_e uses MTPA table lookup for external characteristic"
```

---

### Task 5: Refactor Sub-flow D — Efficiency MAP via MTPA + FEA TR

**Files:**
- Modify: `scripts/subflow_d_efficiency_map.py`

**Interfaces:**
- Consumes: `mtpa_table.npz`, `psi_dq_table.npz`
- Produces: 6 MAP CSVs (speed×torque)

- [ ] **Step 1: Rewrite `scripts/subflow_d_efficiency_map.py`**

Keep `_sync_motion_angular_velocity()` and `_resample_to_torque_grid()` unchanged. Replace the rest:

```python
"""Sub-flow D: 全域效率 MAP
Solver: TransientXY (5_Partial_motor_TR)
Method: FEA scan on speed×Is grid with MTPA Thet, resample to speed×torque.

Usage:
    python -c "from scripts.subflow_d_efficiency_map import run; run(speed_steps=10, Is_steps=10, mtpa_table_path='pmsm_projects/xxx/mtpa_table.npz')"
"""

import sys, re, json, time as _time
sys.stdout.reconfigure(encoding='utf-8')
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
_DEFAULT_TEMPLATE = str(ROOT / 'references' / 'Prius_2D_Practice.aedt')

# ── Default motor parameters (used only for Rs in efficiency calc) ──
DEFAULT_Rs = 0.05  # Ω
DEFAULT_P = 4


def _sync_motion_angular_velocity(aedt_path, target_rpm):
    """Sync .aedt file Angular Velocity and Speed_rpm (unchanged from old version)."""
    file = Path(aedt_path)
    if not file.exists():
        return
    content = file.read_text(encoding='utf-8')
    new_content, n_av = re.subn(
        r"'Angular Velocity'='[^']*'",
        f"'Angular Velocity'='{target_rpm}rpm'", content)
    new_content, n_sp = re.subn(
        r"VariableProp\('Speed_rpm', 'UD', '', '[^']*'",
        f"VariableProp('Speed_rpm', 'UD', '', '{target_rpm}'", new_content)
    if n_av > 0 or n_sp > 0:
        file.write_text(new_content, encoding='utf-8')
        print(f'  [D] .aedt synced: AV={target_rpm}rpm ({n_av}), Speed_rpm={target_rpm} ({n_sp})')


def _resample_to_torque_grid(results, speed_points, torque_steps):
    """Resample (speed×Is) → (speed×torque) uniform grid.

    UNCHANGED from old version. Kept as-is.
    """
    n_speeds = len(speed_points)
    n_is = results.shape[1]

    all_tq = []
    for i in range(n_speeds):
        for j in range(n_is):
            r = results[i, j]
            if r and r.get('T_avg', 0) > 0.1:
                all_tq.append(r['T_avg'])
    if not all_tq:
        return None, None
    T_max = max(all_tq)

    T_points = np.linspace(0, T_max, torque_steps)

    interp_keys = ['Is', 'beta', 'eta', 'M', 'PF', 'Vs',
                   'P_cu', 'P_fe', 'P_mag', 'P_total']

    resampled = np.zeros((n_speeds, torque_steps), dtype=object)

    for i in range(n_speeds):
        tq_row, vals = [], {k: [] for k in interp_keys}
        for j in range(n_is):
            r = results[i, j]
            if r is None:
                continue
            if r.get('Is', 0) > 1 and r.get('T_avg', 0) < 0.5:
                continue
            tq_row.append(r.get('T_avg', 0))
            for k in interp_keys:
                vals[k].append(r.get(k, 0))

        if len(tq_row) < 2:
            continue

        order = np.argsort(tq_row)
        tq_sorted = np.array(tq_row)[order]

        for j, T_target in enumerate(T_points):
            if T_target < 1e-6:
                entry = {k: 0.0 for k in interp_keys}
                entry['T_avg'] = 0.0
                resampled[i, j] = entry
            elif T_target <= tq_sorted[-1]:
                entry = {'T_avg': T_target}
                for k in interp_keys:
                    v_sorted = np.array(vals[k])[order]
                    entry[k] = float(np.interp(T_target, tq_sorted, v_sorted))
                resampled[i, j] = entry
            else:
                resampled[i, j] = None

    return resampled, T_points


def run(speed_min=500, speed_max=8000, speed_steps=10, Is_steps=10,
        vdc=450, imax=350, mtpa_table_path=None, psi_dq_table_path=None,
        pole_pairs=4, Rs=0.05, project_path=None):
    """Efficiency MAP via FEA Transient with MTPA table lookup.

    Parameters
    ----------
    speed_steps : int — number of speed points
    Is_steps : int — number of current points
    vdc : float — DC bus voltage
    imax : float — max phase current
    mtpa_table_path : str — path to mtpa_table.npz
    psi_dq_table_path : str — path to psi_dq_table.npz (for analytical calc)
    pole_pairs : int
    Rs : float — phase resistance (Ω)
    project_path : str — project directory
    """
    from scripts.project_utils import get_template_path
    from scripts.param_guard import check_var
    from ansys.aedt.core import Maxwell2d

    # Resolve paths
    proj_dir = Path(project_path) if project_path else Path.cwd()
    if mtpa_table_path:
        mtpa_path = Path(mtpa_table_path)
    else:
        mtpa_path = proj_dir / 'mtpa_table.npz'
    if psi_dq_table_path:
        psi_dq_path = Path(psi_dq_table_path)
    else:
        psi_dq_path = proj_dir / 'psi_dq_table.npz'

    # Load MTPA table
    mtpa = np.load(mtpa_path)
    speed_grid = mtpa['speed_grid']
    Is_grid = mtpa['Is_grid']
    Thet_ref = mtpa['Thet_opt']
    feasible_ref = mtpa['feasible']

    Vmax = vdc / np.sqrt(3)

    # Build speed×Is scan grid
    speed_points = np.linspace(speed_min, speed_max, speed_steps)
    Is_points = np.linspace(0, imax, Is_steps)  # FEA grid (not necessarily same as mtpa)

    print(f'\n{"=" * 60}')
    print(f'  Sub-flow D: Efficiency MAP (MTPA + FEA TR)')
    print(f'{"=" * 60}')
    print(f'  Speed: {speed_min}~{speed_max} rpm, {speed_steps} points')
    print(f'  Is:    0~{imax} A, {Is_steps} points')
    print(f'  Total: {speed_steps * Is_steps} potential points')
    print(f'  MTPA table: {mtpa_path}')
    print(f'{"=" * 60}')

    # Template setup
    if project_path:
        template = get_template_path(Path(project_path))
    else:
        template = _DEFAULT_TEMPLATE
    _sync_motion_angular_velocity(template, int(speed_min))

    m2d = Maxwell2d(
        project=template, design='5_Partial_motor_TR',
        solution_type='TransientXY',
        non_graphical=False, new_desktop=False, close_on_exit=False
    )
    P = pole_pairs

    results = np.full((speed_steps, Is_steps), None, dtype=object)
    total = speed_steps * Is_steps
    count = 0

    # ── Checkpoint resume ──
    ckpt_file = proj_dir / '_subflow_d_checkpoint.json'
    start_i, start_j = 0, 0
    if ckpt_file.exists():
        ckpt = json.loads(ckpt_file.read_text(encoding='utf-8'))
        start_i = ckpt.get('i', 0)
        start_j = ckpt.get('j', 0)
        for entry in ckpt.get('results', []):
            results[entry['i'], entry['j']] = entry['data']
        done = sum(1 for r in results.flat if r is not None)
        print(f'  [Checkpoint] 恢复 {done} 个已计算点，从 ({start_i},{start_j}) 继续')

    for i in range(start_i, len(speed_points)):
        n = speed_points[i]
        # Update speed via boundary API
        try:
            setup = m2d.setups[0]
            motion = [b for b in setup.boundaries if b.name == 'MotionSetup1'][0]
            motion.props['AngularVelocity'] = f'{n}rpm'
            motion.update()
        except Exception:
            pass

        for j in range(start_j if i == start_i else 0, len(Is_points)):
            Is = Is_points[j]
            count += 1
            current = count
            print(f'\n[D] {current}/{total} n={n:.0f}rpm Is={Is:.1f}A ...', flush=True)

            # Look up MTPA Thet from table (nearest speed, nearest Is)
            i_spd = np.argmin(np.abs(speed_grid - n))
            i_Is = np.argmin(np.abs(Is_grid - Is))

            if not feasible_ref[i_spd, i_Is] or Is < 1e-6:
                print(f'  → 跳过 (不可行或零电流)')
                results[i, j] = {
                    'T_avg': 0.0, 'Is': Is, 'beta': 0.0,
                    'eta': 0.0, 'M': 0.0, 'PF': 0.0, 'Vs': 0.0,
                    'P_cu': 0.0, 'P_fe': 0.0, 'P_mag': 0.0, 'P_total': 0.0,
                }
                _save_checkpoint(ckpt_file, results, i, j, speed_points, Is_points)
                continue

            Thet = Thet_ref[i_spd, i_Is]
            if np.isnan(Thet):
                print(f'  → 跳过 (MTPA 无解)')
                results[i, j] = None
                _save_checkpoint(ckpt_file, results, i, j, speed_points, Is_points)
                continue

            # Set FEA parameters
            check_var('Thet_deg', 'D')
            m2d['Thet_deg'] = str(round(Thet, 1))
            check_var('Imax', 'D')
            m2d['Imax'] = f'{Is}A'

            # Time setup
            freq = n / 60 * P
            stop_time = 1 / freq
            time_step = 1 / (freq * 50)
            setup = m2d.setups[0]
            setup.props['StopTime'] = f'{stop_time}s'
            setup.props['TimeStep'] = f'{time_step}s'
            setup.update()

            m2d.save_project()

            # Solve with gRPC retry
            T_avg = _solve_one_point(m2d)
            print(f'  T_avg={T_avg:.2f} Nm')

            # Analytical calculations via ψd/ψq interpolation
            analytical = None
            if psi_dq_path.exists():
                analytical = _calc_from_psi_table(
                    psi_dq_path, P, n, Is, Thet, T_avg, Vmax, Rs)

            if analytical:
                analytical['T_avg'] = T_avg
                results[i, j] = analytical
            else:
                results[i, j] = {
                    'T_avg': T_avg, 'Is': Is, 'beta': Thet,
                    'eta': 0.0, 'M': 0.0, 'PF': 0.0, 'Vs': 0.0,
                    'P_cu': 0.0, 'P_fe': 0.0, 'P_mag': 0.0, 'P_total': 0.0,
                }

            _save_checkpoint(ckpt_file, results, i, j, speed_points, Is_points)

    m2d.close_project()

    # ── Resample to speed×torque grid ──
    resampled, T_points = _resample_to_torque_grid(results, speed_points, Is_steps)

    if resampled is None:
        print('\n  No valid FEA points — skipping MAP output.')
        return None, None

    # ── Save MAP CSVs ──
    out_dir = proj_dir / 'results'
    out_dir.mkdir(exist_ok=True)

    _save_maps(resampled, T_points, speed_points, out_dir)

    # Stats
    eta_vals = [r['eta'] for r in (resampled.flat if resampled is not None else [])
                if r and r.get('eta', 0) > 0.1]
    if eta_vals:
        print(f'\n  Peak η = {np.max(eta_vals):.1f}%')

    return resampled, T_points


def _solve_one_point(m2d, max_retries=3):
    """FEA solve with gRPC retry."""
    import time
    for attempt in range(max_retries):
        try:
            m2d.analyze('Setup1')
            data = m2d.post.get_solution_data_per_variation(
                expressions=['Moving1.Torque'])
            if data:
                tq = np.array(data.data_real('Moving1.Torque'))
                if len(tq) >= 1:
                    return float(np.mean(tq[-max(len(tq)//2, 1):]))
            return 0.0
        except Exception as e:
            err_msg = str(e)
            if 'GrpcApiError' in type(e).__name__ or 'gRPC' in err_msg:
                if attempt < max_retries - 1:
                    wait = 2 ** attempt
                    print(f' [gRPC retry {attempt+1}/{max_retries}, wait {wait}s]', end='', flush=True)
                    time.sleep(wait)
                    continue
            raise
    return 0.0


def _calc_from_psi_table(psi_dq_path, pole_pairs, speed, Is, Thet, T_avg, Vmax, Rs):
    """Compute Vs, M, PF, efficiency from ψd/ψq table interpolation."""
    import scipy
    from scipy.interpolate import RegularGridInterpolator

    data = np.load(psi_dq_path)
    id_grid = data['id_grid']
    iq_grid = data['iq_grid']

    psi_d_ip = RegularGridInterpolator(
        (id_grid, iq_grid), data['psi_d'], bounds_error=False, fill_value=None)
    psi_q_ip = RegularGridInterpolator(
        (id_grid, iq_grid), data['psi_q'], bounds_error=False, fill_value=None)

    beta = np.radians(Thet)
    Id = Is * np.sin(beta)
    Iq = Is * np.cos(beta)

    try:
        pd = float(psi_d_ip([[Id, Iq]])[0])
        pq = float(psi_q_ip([[Id, Iq]])[0])
    except Exception:
        ix = np.clip(np.searchsorted(id_grid, Id), 0, len(id_grid)-1)
        iy = np.clip(np.searchsorted(iq_grid, Iq), 0, len(iq_grid)-1)
        pd = float(data['psi_d'][ix, iy])
        pq = float(data['psi_q'][ix, iy])

    omega_e = speed * np.pi / 30 * pole_pairs
    Vd = Rs * Id - omega_e * pq
    Vq = Rs * Iq + omega_e * pd
    Vs = np.sqrt(Vd**2 + Vq**2)
    M = Vs / Vmax if Vmax > 0 else 0.0

    phi_v = np.arctan2(Vd, Vq)
    phi_i = np.arctan2(Id, Iq)
    PF = np.cos(phi_v - phi_i)

    P_out = T_avg * speed * 2 * np.pi / 60
    P_cu = 3 * (Is / np.sqrt(2))**2 * Rs
    P_loss = P_cu
    eta = P_out / (P_out + P_loss) * 100 if (P_out + P_loss) > 0 else 0.0

    return {
        'Is': Is, 'beta': Thet,
        'Vs': Vs, 'M': M, 'PF': PF,
        'P_cu': P_cu / 1000, 'P_fe': 0.0, 'P_mag': 0.0, 'P_total': P_loss / 1000,
        'eta': eta,
    }


def _save_checkpoint(ckpt_file, results, i, j, speed_points, Is_points):
    """Persist progress after each point."""
    entries = []
    for ii in range(len(speed_points)):
        for jj in range(len(Is_points)):
            r = results[ii, jj]
            if r is not None:
                entries.append({'i': int(ii), 'j': int(jj), 'data': r})
    ckpt = {
        'i': int(i), 'j': int(j),
        'speed_points': [float(s) for s in speed_points],
        'Is_points': [float(s) for s in Is_points],
        'results': entries,
    }
    ckpt_file.write_text(json.dumps(ckpt, ensure_ascii=False, indent=2),
                         encoding='utf-8')


def _save_maps(resampled, T_points, speed_points, out_dir):
    """Save 6 MAP CSVs in speed×torque format."""
    T_labels = [f'{t:.0f}' for t in T_points]
    col_labels = ',' + ','.join(T_labels)

    maps = {
        'efficiency_map': 'eta',
        'torque_map': 'T_avg',
        'current_map': 'Is',
        'modulation_map': 'M',
        'pf_map': 'PF',
        'voltage_map': 'Vs',
    }
    for fname, key in maps.items():
        rows = [f'# {fname} — 行:转速(rpm), 列:扭矩(Nm)']
        rows.append(col_labels)
        for i, n in enumerate(speed_points):
            vals = [str(resampled[i, j][key]) if resampled[i, j] else ''
                    for j in range(len(T_points))]
            rows.append(f'{n:.0f},{",".join(vals)}')
        (out_dir / f'{fname}.csv').write_text('\n'.join(rows), encoding='utf-8')
    print(f'\n  MAPs saved to {out_dir}')
```

- [ ] **Step 2: Commit**

```bash
git add scripts/subflow_d_efficiency_map.py
git commit -m "refactor: subflow_d uses MTPA table for FEA angle selection + ψd/ψq analytical calc"
```

---

### Task 6: Update SKILL.md

**Files:**
- Modify: `SKILL.md`

- [ ] **Step 1: Update sub-flow mappings, naming, and references**

Update SKILL.md to reflect:
- Sub-flow A renamed to "MTPA 标定 (Ld/Lq MAP + 主磁链 Φ)" → "MTPA 标定 (ψd/ψq 磁链表 + MTPA 搜索)"
- Sub-flow C removed
- Sub-flow D updated interface description
- Sub-flow E updated interface description
- Allowed params table updated (remove C)
- Add `test_theta_r.py` to the allowed scripts list

- [ ] **Step 2: Commit**

```bash
git add SKILL.md
git commit -m "docs: update SKILL.md for MTPA calibration refactor"
```

---

### Task 7: Final Validation — Step 0 + Phase 1 (Sparse) Dry Run

- [ ] **Step 1: Run θr calibration test**

```bash
python scripts/test_theta_r.py
```

Expected: passes validation (|ψq/ψd| < 0.05)

- [ ] **Step 2: Run Phase 1 + Phase 2 dry run (sparse 10×10)**

```bash
python -c "from scripts.subflow_a_mtpa_cal import run; run(n_id=10, n_iq=10)"
```

Check: `psi_dq_table.npz` and `mtpa_table.npz` are created, no NaNs in expected regions.

---

### Task 8: Final Git Push

```bash
git push -u origin subflow-update
```
