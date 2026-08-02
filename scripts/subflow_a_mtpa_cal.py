"""Sub-flow A: MTPA Calibration
Phase 1: ψd/ψq 2D flux-linkage table (Magnetostatic, Id×Iq grid)
         Default: FW-region densified (20 Id × 10 Iq), ~200 points, ~10 min
Phase 2: MTPA calibration → (n, Is) → (Thet_opt, Id_opt, T_mtpa, Vs, M)

Usage:
    # Full pipeline (Phase 1 + Phase 2) — default FW-densified 20×10 grid
    python -c "from scripts.subflow_a_mtpa_cal import run; run(project_path='pmsm_projects/xxx')"

    # Phase 1 with custom grid
    python -c "from scripts.subflow_a_mtpa_cal import phase1_psi_dq_table; phase1_psi_dq_table(n_id=10, n_iq=10, project_path='pmsm_projects/xxx')"

    # Phase 2 only (requires psi_dq_table.npz)
    python -c "from scripts.subflow_a_mtpa_cal import phase2_mtpa_calibrate; phase2_mtpa_calibrate('pmsm_projects/xxx/psi_dq_table.npz')"
"""

import sys, os, re, tempfile, json, time
sys.stdout.reconfigure(encoding='utf-8')
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
_DEFAULT_TEMPLATE = str(ROOT / 'references' / 'Prius_2D_Practice.aedt')

TURNS = 9
SYMMETRY_MULTIPLIER = 8  # 1/8 partial motor model → scale to full motor

COIL_CONFIG = [
    ('PhaseA1', True), ('PhaseA2', True),
    ('PhaseB1', True), ('PhaseB2', True),
    ('PhaseC1', True), ('PhaseC2', True),
]
COIL_NAMES = [c[0] for c in COIL_CONFIG]

# ── Coordinate transforms ──

def id_iq_to_abc(Id, Iq, theta=0):
    """Inverse Park + Clarke: dq → abc (theta in radians).

    Convention: q-axis LAGS d-axis by 90° (Maxwell motor convention).
    """
    Ialpha = Id * np.cos(theta) + Iq * np.sin(theta)
    Ibeta  = Id * np.sin(theta) - Iq * np.cos(theta)
    Ia = Ialpha
    Ib = -0.5 * Ialpha + np.sqrt(3) / 2 * Ibeta
    Ic = -0.5 * Ialpha - np.sqrt(3) / 2 * Ibeta
    return Ia, Ib, Ic


def park_abc_to_dq(psi_a, psi_b, psi_c, theta=0):
    """Park transform: abc → dq (theta in radians).

    Convention: q-axis LAGS d-axis by 90° (consistent with id_iq_to_abc).
    """
    c = np.cos(theta)
    s = np.sin(theta)
    psi_d = 2/3 * (psi_a * c + psi_b * np.cos(theta - 2*np.pi/3) + psi_c * np.cos(theta + 2*np.pi/3))
    psi_q = 2/3 * (psi_a * s + psi_b * np.sin(theta - 2*np.pi/3) + psi_c * np.sin(theta + 2*np.pi/3))
    return psi_d, psi_q


def _set_coil_currents(m2d, boundaries, Ia, Ib, Ic):
    """Set coil currents via boundary.update() (含 9 匝因子)."""
    currents = [TURNS * Ia, TURNS * Ia,
                TURNS * Ib, TURNS * Ib,
                TURNS * Ic, TURNS * Ic]
    for bnd, (name, is_pos), val in zip(boundaries, COIL_CONFIG, currents):
        bnd.props['Current'] = f'{val}A'
        bnd.update()


def calibrate_theta_r(project_path=None, stack_length=None):
    """Calibrate rotor electrical angle θ_r from PM flux at zero current.

    Method: zero-current solve → measure ψ_abc → solve for θ_r where
    ψ_q = park_abc_to_dq(ψ_a, ψ_b, ψ_c, θ_r) = 0 and ψ_d > 0.

    Analytical solution:
      θ_r = atan2( √3/2·(ψ_b - ψ_c),  ψ_a - ½·ψ_b - ½·ψ_c )

    Parameters
    ----------
    project_path : str, optional — project directory path
    stack_length : float, optional — motor stack length in mm.
                   Sets m2d.model_depth before solve. If None, uses file default.

    Returns
    -------
    theta_r : float — rotor electrical angle in radians [0, 2π)
    phi_0  : float — PM flux linkage ψ_d(0,0) in Wb (scaled to full motor)
    """
    from scripts.project_utils import get_template_path
    from ansys.aedt.core import Maxwell2d

    if project_path:
        template = get_template_path(Path(project_path))
    else:
        template = _DEFAULT_TEMPLATE

    print(f'\n{"=" * 60}')
    print(f'  θ_r Calibration: Zero-Current PM Flux Method')
    print(f'{"=" * 60}')

    m2d = Maxwell2d(
        project=template, design='4_Partial_motor_MS2',
        solution_type='MagnetostaticXY',
        non_graphical=False, new_desktop=False, close_on_exit=False
    )

    # Apply stack length before solve
    if stack_length is not None:
        m2d.model_depth = f'{stack_length}mm'
        print(f'  ModelDepth set to {stack_length} mm')

    boundaries = []
    for name, _ in COIL_CONFIG:
        b = [b for b in m2d.boundaries if b.name == name][0]
        boundaries.append(b)

    # Zero current
    print('  Setting Ia=Ib=Ic=0...')
    _set_coil_currents(m2d, boundaries, 0, 0, 0)

    # Solve
    print('  Solving...')
    m2d.analyze('Setup1')
    out_file = os.path.join(tempfile.gettempdir(), 'pmsm_matrix_export.txt')
    m2d.export_matrix('Matrix1', out_file)
    with open(out_file, 'r') as f:
        psi_coil = _parse_matrix_flux_linkage(f.read())

    psi_a = TURNS * (psi_coil['PhaseA1'] + psi_coil['PhaseA2'])
    psi_b = TURNS * (psi_coil['PhaseB1'] + psi_coil['PhaseB2'])
    psi_c = TURNS * (psi_coil['PhaseC1'] + psi_coil['PhaseC2'])

    # Scale from 1/8 partial model to full motor
    psi_a *= SYMMETRY_MULTIPLIER
    psi_b *= SYMMETRY_MULTIPLIER
    psi_c *= SYMMETRY_MULTIPLIER

    print(f'  ψ_a = {psi_a:.6f} Wb')
    print(f'  ψ_b = {psi_b:.6f} Wb')
    print(f'  ψ_c = {psi_c:.6f} Wb')

    # Analytical solution for θ_r where ψ_q=0
    num = np.sqrt(3) / 2 * (psi_b - psi_c)
    den = psi_a - 0.5 * psi_b - 0.5 * psi_c
    theta_r = np.arctan2(num, den)

    # Verify ψ_d at this angle (should be > 0 for PM flux)
    pd, pq = park_abc_to_dq(psi_a, psi_b, psi_c, theta=theta_r)
    if pd < 0:
        theta_r += np.pi  # flip 180° if d-axis points opposite

    phi_0 = abs(pd)

    print(f'  θ_r = {np.degrees(theta_r):.4f}° ({theta_r:.6f} rad)')
    print(f'  ψ_d = {pd:.6f} Wb, ψ_q = {pq:.6f} Wb')
    print(f'  Φ_0 (PM flux linkage) = {phi_0:.6f} Wb')
    print(f'{"=" * 60}')

    try:
        m2d.close_project()
    except Exception:
        pass

    return theta_r, phi_0


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


def _extract_flux_linkages(m2d, max_retries=2):
    """Solve → export matrix → parse flux linkage per coil.

    Retries on gRPC export failures (intermittent comm issues in pyaedt 0.25.1).
    """
    m2d.analyze('Setup1')
    out_file = os.path.join(tempfile.gettempdir(), 'pmsm_matrix_export.txt')

    last_err = None
    for attempt in range(max_retries + 1):
        try:
            m2d.export_matrix('Matrix1', out_file)
            with open(out_file, 'r') as f:
                text = f.read()
            return _parse_matrix_flux_linkage(text)
        except Exception as e:
            last_err = e
            if attempt < max_retries:
                print(f' [retry {attempt+1}]', end='', flush=True)
                time.sleep(2)
    raise last_err


# ═══════════════════════════════════════════
# Phase 1: ψd/ψq 2D Flux-Linkage Table
# ═══════════════════════════════════════════

def _default_fw_dense_id_grid():
    """Default densified Id grid with FW-region focus.

    Dense in [-400, 0] (FW + MTPA region), sparse in [0, 400].
    20 points total: 13 neg (33A step) + 7 pos (58A step).
    """
    id_neg = np.linspace(-400, 0, 13)
    id_pos = np.linspace(50, 400, 7)
    return np.concatenate([id_neg, id_pos])


def _default_fw_dense_iq_grid():
    """Default densified Iq grid with low-Iq focus for FW cross-saturation.

    Dense in [0, 200] (cross-sat sensitive region), sparse in [200, 400].
    15 points total: 11 low (20A step) + 5 high (50A step).
    """
    iq_low = np.linspace(0, 200, 11)
    iq_high = np.linspace(250, 400, 4)
    return np.concatenate([iq_low, iq_high])


def _default_iq_grid():
    """Default uniform Iq grid (legacy 10-point)."""
    return np.linspace(0, 400, 10)


def phase1_psi_dq_table(id_range=(-400, 400), iq_range=(0, 400),
                        n_id=None, n_iq=None,
                        id_points=None, iq_points=None,
                        theta_r='auto',
                        project_path=None, confirm=True,
                        stack_length=None):
    """Scan Id×Iq grid in Magnetostatic, build ψd/ψq 2D table.

    Grid can be specified in two ways:
      - Explicit points: pass id_points / iq_points arrays directly.
      - Uniform grid: pass id_range / n_id and iq_range / n_iq.
      - Default (no args): densified FW-region grid (20 Id × 10 Iq).

    Parameters
    ----------
    id_range : tuple — (Id_min, Id_max) in A (ignored if id_points given)
    iq_range : tuple — (Iq_min, Iq_max) in A (ignored if iq_points given)
    n_id : int — number of Id steps (ignored if id_points given)
    n_iq : int — number of Iq steps (ignored if iq_points given)
    id_points : np.ndarray, optional — explicit Id grid points
    iq_points : np.ndarray, optional — explicit Iq grid points
    theta_r : float or 'auto' — rotor d-axis angle in radians.
              'auto' → auto-calibrate via zero-current PM flux method
    project_path : str, optional — project directory path
    confirm : bool — if True, prompt user before running
    stack_length : float, optional — motor stack length in mm.
                   Sets m2d.model_depth before FEA. If None, uses file default.

    Returns
    -------
    id_grid : np.ndarray (n_id,)
    iq_grid : np.ndarray (n_iq,)
    psi_d : np.ndarray (n_id, n_iq) — scaled to full motor (×8 symmetry)
    psi_q : np.ndarray (n_id, n_iq) — scaled to full motor (×8 symmetry)
    """
    from scripts.project_utils import get_template_path
    from ansys.aedt.core import Maxwell2d

    # Auto-calibrate θ_r if requested
    if theta_r == 'auto':
        theta_r, phi_0 = calibrate_theta_r(project_path=project_path,
                                            stack_length=stack_length)

    if project_path:
        template = get_template_path(Path(project_path))
    else:
        template = _DEFAULT_TEMPLATE

    # Build Id/Iq grids
    if id_points is not None:
        id_grid = np.asarray(id_points, dtype=float)
    elif n_id is not None:
        id_grid = np.linspace(id_range[0], id_range[1], n_id)
    else:
        id_grid = _default_fw_dense_id_grid()

    if iq_points is not None:
        iq_grid = np.asarray(iq_points, dtype=float)
    elif n_iq is not None:
        iq_grid = np.linspace(iq_range[0], iq_range[1], n_iq)
    else:
        iq_grid = _default_fw_dense_iq_grid()

    n_id, n_iq = len(id_grid), len(iq_grid)
    total = n_id * n_iq

    print(f'\n{"=" * 60}')
    print(f'  Phase 1: ψd/ψq 2D Flux-Linkage Table')
    print(f'  Solver: MagnetostaticXY (4_Partial_motor_MS2)')
    print(f'{"=" * 60}')
    print(f'  Id grid: {id_grid[0]:.0f} ~ {id_grid[-1]:.0f} A, {n_id} points')
    if n_id <= 30:
        print(f'           {np.array2string(id_grid, precision=0, max_line_width=120)}')
    print(f'  Iq grid: {iq_grid[0]:.0f} ~ {iq_grid[-1]:.0f} A, {n_iq} points')
    print(f'  Total:   {total} points, ~{total * 3 // 60} min')
    print(f'  θr:      {np.degrees(theta_r):.1f}°')
    print(f'{"=" * 60}')

    # Confirm mode
    config_path = ROOT / 'config.json'
    if confirm and config_path.exists():
        with open(config_path) as f:
            mode = json.load(f).get('execution_mode', 'confirm')
    else:
        mode = 'auto' if not confirm else 'confirm'

    if mode == 'confirm':
        resp = input('\n  继续执行? [Y/n] ').strip().lower()
        if resp and resp != 'y':
            print('  已取消。')
            return None, None, None, None

    m2d = Maxwell2d(
        project=template, design='4_Partial_motor_MS2',
        solution_type='MagnetostaticXY',
        non_graphical=False, new_desktop=False, close_on_exit=False
    )

    # Apply stack length before FEA scan
    if stack_length is not None:
        m2d.model_depth = f'{stack_length}mm'
        print(f'  ModelDepth set to {stack_length} mm')

    boundaries = []
    for name, _ in COIL_CONFIG:
        b = [b for b in m2d.boundaries if b.name == name][0]
        boundaries.append(b)

    psi_d = np.full((n_id, n_iq), np.nan)
    psi_q = np.full((n_id, n_iq), np.nan)

    out_dir = Path(project_path) if project_path else Path.cwd()
    out_path = out_dir / 'psi_dq_table.npz'

    for i, Id_val in enumerate(id_grid):
        for j, Iq_val in enumerate(iq_grid):
            point = i * n_iq + j + 1
            Ia, Ib, Ic = id_iq_to_abc(Id_val, Iq_val, theta=theta_r)
            print(f'  [P1] {point}/{total} Id={Id_val:+7.1f}A Iq={Iq_val:+7.1f}A '
                  f'Ia={Ia:.1f}A ...', end='', flush=True)

            try:
                _set_coil_currents(m2d, boundaries, Ia, Ib, Ic)
                psi_coil = _extract_flux_linkages(m2d)

                psi_a = TURNS * (psi_coil['PhaseA1'] + psi_coil['PhaseA2'])
                psi_b = TURNS * (psi_coil['PhaseB1'] + psi_coil['PhaseB2'])
                psi_c = TURNS * (psi_coil['PhaseC1'] + psi_coil['PhaseC2'])

                # Scale from 1/8 partial model to full motor
                psi_a *= SYMMETRY_MULTIPLIER
                psi_b *= SYMMETRY_MULTIPLIER
                psi_c *= SYMMETRY_MULTIPLIER

                pd, pq = park_abc_to_dq(psi_a, psi_b, psi_c, theta=theta_r)
                psi_d[i, j] = pd
                psi_q[i, j] = pq
                print(f' ψd={pd:.6f} ψq={pq:.6f}')
            except Exception as e:
                print(f' FAIL: {e}')
                psi_d[i, j] = np.nan
                psi_q[i, j] = np.nan

            # Checkpoint save every 10 points (guard against mid-run crashes)
            if point % 10 == 0:
                np.savez(out_path,
                         id_grid=id_grid, iq_grid=iq_grid,
                         psi_d=psi_d, psi_q=psi_q,
                         theta_r=theta_r,
                         symmetry_multiplier=SYMMETRY_MULTIPLIER,
                         stack_length=stack_length)
                print(f'  [P1] checkpoint saved ({point}/{total})')

    # Final save (must be BEFORE close_project — pyaedt 0.25.1 close bug)
    np.savez(out_path,
             id_grid=id_grid, iq_grid=iq_grid,
             psi_d=psi_d, psi_q=psi_q,
             theta_r=theta_r,
             symmetry_multiplier=SYMMETRY_MULTIPLIER,
             stack_length=stack_length)
    print(f'\n  Saved: {out_path}')

    try:
        m2d.close_project()
    except Exception as e:
        print(f'  [WARN] close_project failed (known pyaedt 0.25.1 bug): {e}')

    return id_grid, iq_grid, psi_d, psi_q


# ═══════════════════════════════════════════
# Phase 2: MTPA Calibration
# ═══════════════════════════════════════════

def _build_interpolator(x_grid, y_grid, z_table):
    """Build a RegularGridInterpolator with NN fallback."""
    from scipy.interpolate import RegularGridInterpolator
    interp = RegularGridInterpolator(
        (x_grid, y_grid), z_table, bounds_error=False, fill_value=None)

    def fn(x, y):
        try:
            val = float(interp([[x, y]])[0])
            if np.isnan(val):
                raise ValueError('NaN')
            return val
        except Exception:
            ix = np.clip(np.searchsorted(x_grid, x), 0, len(x_grid) - 1)
            iy = np.clip(np.searchsorted(y_grid, y), 0, len(y_grid) - 1)
            return float(z_table[ix, iy])
    return fn


def phase2_mtpa_calibrate(psi_dq_path, vdc=450, imax=350,
                          speed_min=500, speed_max=8000, n_speed=11, n_Is=21,
                          Rs=0.05, pole_pairs=4):
    """MTPA calibration: search (n, Is) → optimal Thet/Id/Iq using ψd/ψq table.

    For each (speed, Is), sweeps current angle β ∈ [-90°, 90°],
    selects β that maximizes torque subject to voltage constraint Vs ≤ Vmax.

    Parameters
    ----------
    psi_dq_path : str or Path — path to psi_dq_table.npz (Phase 1 output)
    vdc : float — DC bus voltage (V)
    imax : float — maximum phase current (A)
    speed_min, speed_max : float — speed range (rpm)
    n_speed : int — number of speed points
    n_Is : int — number of current points
    Rs : float — phase resistance (Ω)
    pole_pairs : int

    Returns
    -------
    dict with keys: speed_grid, Is_grid, Thet_opt, T_mtpa, Vs, M, feasible
    """
    Vmax = vdc / np.sqrt(3)

    # Load ψd/ψq table
    data = np.load(psi_dq_path)
    id_grid = data['id_grid']
    iq_grid = data['iq_grid']
    psi_d_tab = data['psi_d']
    psi_q_tab = data['psi_q']

    psi_d_fn = _build_interpolator(id_grid, iq_grid, psi_d_tab)
    psi_q_fn = _build_interpolator(id_grid, iq_grid, psi_q_tab)

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
            best = None  # (beta, Id, Iq, Vs_val)

            for beta in beta_grid:
                Id = Is_val * np.sin(beta)
                Iq = Is_val * np.cos(beta)

                psi_d_val = psi_d_fn(Id, Iq)
                psi_q_val = psi_q_fn(Id, Iq)

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
# Phase 3: Outer Characteristic (smooth T-n / P-n curve)
# ═══════════════════════════════════════════

def outer_characteristic(psi_dq_path, vdc=500, imax=200,
                         speed_min=500, speed_max=10000, n_speed=100,
                         Rs=0.05, pole_pairs=4):
    """Compute smooth outer characteristic (T-n, P-n curves).

    Vectorized search: for each speed, builds a 2D (Is × β) grid and
    evaluates all points at once via batch interpolation, then finds
    the torque peak with sub-step accuracy via quadratic interpolation.

    Parameters
    ----------
    psi_dq_path : str or Path — path to psi_dq_table.npz (Phase 1 output)
    vdc : float — DC bus voltage (V)
    imax : float — maximum phase current (A)
    speed_min, speed_max : float — speed range (rpm)
    n_speed : int — number of speed points
    Rs : float — phase resistance (Ω)
    pole_pairs : int

    Returns
    -------
    dict with keys: speed, T_max, P_max_kW, Is_opt, Id_opt, Iq_opt, beta_opt, Vs, region
    """
    from scipy.interpolate import RegularGridInterpolator

    Vmax = vdc / np.sqrt(3)

    # Load ψd/ψq table
    data = np.load(psi_dq_path)
    id_grid = data['id_grid']
    iq_grid = data['iq_grid']
    psi_d_tab = data['psi_d']
    psi_q_tab = data['psi_q']

    # Build vectorized interpolators
    psi_d_interp = RegularGridInterpolator(
        (id_grid, iq_grid), psi_d_tab, bounds_error=False, fill_value=None)
    psi_q_interp = RegularGridInterpolator(
        (id_grid, iq_grid), psi_q_tab, bounds_error=False, fill_value=None)

    # Fine search grids
    n_Is = 200        # 1A steps for 0-200A
    n_beta = 300      # angle resolution
    speed_grid = np.linspace(speed_min, speed_max, n_speed)
    Is_fine = np.linspace(0, imax, n_Is)
    beta_grid = np.linspace(-np.pi/2, np.pi/2, n_beta)

    # Meshgrid for vectorized eval: (n_Is, n_beta)
    Is_mesh, beta_mesh = np.meshgrid(Is_fine, beta_grid, indexing='ij')
    Id_mesh = Is_mesh * np.sin(beta_mesh)   # (n_Is, n_beta)
    Iq_mesh = Is_mesh * np.cos(beta_mesh)   # (n_Is, n_beta)

    # Flatten for batch interpolation
    Id_flat = Id_mesh.ravel()
    Iq_flat = Iq_mesh.ravel()
    pts = np.column_stack([Id_flat, Iq_flat])

    # Pre-evaluate ψd, ψq on the full (Id, Iq) grid (independent of speed!)
    print('  Evaluating ψd/ψq on fine grid...', end=' ', flush=True)
    psi_d_flat = psi_d_interp(pts).reshape(n_Is, n_beta)
    psi_q_flat = psi_q_interp(pts).reshape(n_Is, n_beta)
    print('done.')

    # NaN fallback: nearest-neighbor fill
    nan_mask_d = np.isnan(psi_d_flat)
    nan_mask_q = np.isnan(psi_q_flat)
    if nan_mask_d.any() or nan_mask_q.any():
        print(f'  Filling {nan_mask_d.sum()} NaN(s) in ψd, '
              f'{nan_mask_q.sum()} in ψq via NN...')
        # Simple NN: clip Id/Iq to grid bounds
        Id_clipped = np.clip(Id_mesh, id_grid[0], id_grid[-1])
        Iq_clipped = np.clip(Iq_mesh, iq_grid[0], iq_grid[-1])
        i_id = np.clip(np.searchsorted(id_grid, Id_clipped), 1, len(id_grid)-1)
        i_iq = np.clip(np.searchsorted(iq_grid, Iq_clipped), 1, len(iq_grid)-1)
        # Pick the nearest grid point
        for ii in range(n_Is):
            for jj in range(n_beta):
                if np.isnan(psi_d_flat[ii, jj]):
                    i0 = np.clip(i_id[ii, jj], 0, len(id_grid)-1)
                    j0 = np.clip(i_iq[ii, jj], 0, len(iq_grid)-1)
                    psi_d_flat[ii, jj] = psi_d_tab[i0, j0]
                    psi_q_flat[ii, jj] = psi_q_tab[i0, j0]
        print('  done.')

    # Torque on the grid (independent of speed)
    T_mesh = 1.5 * pole_pairs * (psi_d_flat * Iq_mesh - psi_q_flat * Id_mesh)
    T_mesh[T_mesh <= 0] = np.nan

    # Result arrays
    T_max_arr = np.full(n_speed, np.nan)
    Is_opt_arr = np.full(n_speed, np.nan)
    Id_opt_arr = np.full(n_speed, np.nan)
    Iq_opt_arr = np.full(n_speed, np.nan)
    beta_opt_arr = np.full(n_speed, np.nan)
    Vs_opt_arr = np.full(n_speed, np.nan)
    region_arr = np.full(n_speed, '', dtype=object)

    print(f'\n{"=" * 60}')
    print(f'  Outer Characteristic')
    print(f'  Vdc={vdc:.0f}V, Vmax={Vmax:.1f}V, Imax={imax:.0f}A')
    print(f'  Grid: {n_Is} Is × {n_beta} β × {n_speed} speed')
    print(f'{"=" * 60}')

    for i, n in enumerate(speed_grid):
        omega_e = n * np.pi / 30 * pole_pairs

        # Compute Vs on the full grid
        Vd = Rs * Id_mesh - omega_e * psi_q_flat
        Vq = Rs * Iq_mesh + omega_e * psi_d_flat
        Vs_mesh = np.sqrt(Vd**2 + Vq**2)

        # Mask: feasible (voltage constraint) + valid torque
        feasible_mask = (Vs_mesh <= Vmax) & (~np.isnan(T_mesh))
        feasible_mask[0, :] = True  # Is=0 always feasible

        if not feasible_mask.any():
            continue

        # For each Is, find max T among feasible β
        T_best_per_Is = np.full(n_Is, np.nan)
        beta_best_idx = np.full(n_Is, -1, dtype=int)

        for j in range(n_Is):
            feasible_betas = np.where(feasible_mask[j, :])[0]
            if len(feasible_betas) == 0:
                continue
            T_feas = T_mesh[j, feasible_betas]
            best_local = np.argmax(T_feas)
            T_best_per_Is[j] = T_feas[best_local]
            beta_best_idx[j] = feasible_betas[best_local]

        valid = ~np.isnan(T_best_per_Is)
        if valid.sum() == 0:
            continue

        peak_idx = np.argmax(T_best_per_Is[valid])
        T_peak_raw = T_best_per_Is[valid][peak_idx]

        # Quadratic interpolation of T(Is) around the peak → sub-step accuracy
        if peak_idx > 0 and peak_idx < valid.sum() - 1:
            Is_vals_valid = Is_fine[valid]
            T_vals_valid = T_best_per_Is[valid]
            x = Is_vals_valid[[peak_idx-1, peak_idx, peak_idx+1]]
            y = T_vals_valid[[peak_idx-1, peak_idx, peak_idx+1]]

            A = np.column_stack([x**2, x, np.ones(3)])
            try:
                a, b, c = np.linalg.solve(A, y)
                if a < 0:
                    Is_peak = -b / (2 * a)
                    T_peak = a * Is_peak**2 + b * Is_peak + c
                    if x[0] <= Is_peak <= x[2] and T_peak > T_peak_raw * 0.95:
                        Is_opt_arr[i] = Is_peak
                        T_max_arr[i] = T_peak
                        # Interpolate β, Id, Iq, Vs at Is_peak
                        beta_opt_arr[i] = np.interp(Is_peak, Is_fine, beta_grid[beta_best_idx])
                        Id_opt_arr[i] = Is_peak * np.sin(beta_opt_arr[i])
                        Iq_opt_arr[i] = Is_peak * np.cos(beta_opt_arr[i])
                        # Recompute Vs at interpolated point
                        pd = float(psi_d_interp([[Id_opt_arr[i], Iq_opt_arr[i]]])[0])
                        pq = float(psi_q_interp([[Id_opt_arr[i], Iq_opt_arr[i]]])[0])
                        vd = Rs * Id_opt_arr[i] - omega_e * pq
                        vq = Rs * Iq_opt_arr[i] + omega_e * pd
                        Vs_opt_arr[i] = np.sqrt(vd**2 + vq**2)
                        region_arr[i] = 'MTPA' if (Vs_opt_arr[i] < Vmax * 0.99
                                            and Is_peak >= imax * 0.98) else 'FW'
                        continue
            except np.linalg.LinAlgError:
                pass

        # Fallback: discrete peak
        abs_idx = np.where(valid)[0][peak_idx]
        Is_opt_arr[i] = Is_fine[abs_idx]
        T_max_arr[i] = T_peak_raw
        beta_opt_arr[i] = beta_grid[beta_best_idx[abs_idx]]
        Id_opt_arr[i] = Is_opt_arr[i] * np.sin(beta_opt_arr[i])
        Iq_opt_arr[i] = Is_opt_arr[i] * np.cos(beta_opt_arr[i])
        Vs_opt_arr[i] = Vs_mesh[abs_idx, beta_best_idx[abs_idx]]
        region_arr[i] = 'MTPA' if (Vs_opt_arr[i] < Vmax * 0.99
                                   and Is_opt_arr[i] >= imax * 0.98) else 'FW'

        if (i + 1) % 20 == 0:
            print(f'  [{i+1}/{n_speed}] n={n:.0f} rpm: '
                  f'Is*={Is_opt_arr[i]:.1f}A, β*={np.degrees(beta_opt_arr[i]):.1f}°, '
                  f'T*={T_max_arr[i]:.1f}Nm, {region_arr[i]}')

    P_max_kW = T_max_arr * speed_grid * np.pi / 30 / 1000

    # Save
    out_dir = Path(psi_dq_path).parent
    out_path = out_dir / 'outer_characteristic.npz'
    np.savez(out_path,
             speed=speed_grid, T_max=T_max_arr, P_max_kW=P_max_kW,
             Is_opt=Is_opt_arr, Id_opt=Id_opt_arr, Iq_opt=Iq_opt_arr,
             beta_opt=beta_opt_arr, Vs=Vs_opt_arr, region=region_arr,
             vdc=vdc, imax=imax, Rs=Rs, pole_pairs=pole_pairs)
    print(f'  Saved: {out_path}')

    return {'speed': speed_grid, 'T_max': T_max_arr, 'P_max_kW': P_max_kW,
            'Is_opt': Is_opt_arr, 'Id_opt': Id_opt_arr, 'Iq_opt': Iq_opt_arr,
            'beta_opt': beta_opt_arr, 'Vs': Vs_opt_arr, 'region': region_arr}


# ═══════════════════════════════════════════
# Combined run(): Phase 1 + Phase 2
# ═══════════════════════════════════════════

def run(id_range=(-400, 400), iq_range=(0, 400), n_id=None, n_iq=None,
        id_points=None, iq_points=None,
        theta_r='auto', vdc=450, imax=350,
        speed_min=500, speed_max=8000, n_speed=11, n_Is=21,
        Rs=0.05, pole_pairs=4, project_path=None, stack_length=None):
    """Run full MTPA calibration pipeline: Phase 1 + Phase 2.

    Phase 1: MS ψd/ψq flux-linkage table (auto-calibrates θ_r)
    Phase 2: MTPA search → (n,Is)→Thet lookup table

    Grid: pass id_points/iq_points for custom grids, or n_id/n_iq for
    uniform spacing. Default (no args) uses FW-region densified grid.

    Parameters
    ----------
    stack_length : float, optional — motor stack length in mm.
                   Sets m2d.model_depth before FEA. If None, uses file default.
    """
    out_dir = Path(project_path) if project_path else Path.cwd()
    psi_dq_path = out_dir / 'psi_dq_table.npz'

    # Phase 1
    result_p1 = phase1_psi_dq_table(
        id_range=id_range, iq_range=iq_range,
        n_id=n_id, n_iq=n_iq,
        id_points=id_points, iq_points=iq_points,
        theta_r=theta_r,
        project_path=project_path, confirm=True,
        stack_length=stack_length
    )

    if result_p1[0] is None:
        print('\n  Phase 1 cancelled. Skipping Phase 2.')
        return

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
    p.add_argument('--theta-r', default='auto', help='Rotor angle in rad, or "auto" for calibration')
    p.add_argument('--vdc', type=float, default=450)
    p.add_argument('--imax', type=float, default=350)
    p.add_argument('--phase', type=int, default=0, help='1=Phase1 only, 2=Phase2 only, 0=both')
    p.add_argument('--psi-dq-path', default=None, help='Path to psi_dq_table.npz (Phase 2 only)')
    args = p.parse_args()

    theta_r = args.theta_r if args.theta_r == 'auto' else float(args.theta_r)

    if args.phase == 1:
        phase1_psi_dq_table(
            n_id=args.n_id, n_iq=args.n_iq,
            theta_r=theta_r, project_path=args.project_path)
    elif args.phase == 2:
        psi_path = args.psi_dq_path
        if not psi_path:
            proj = Path(args.project_path) if args.project_path else Path.cwd()
            psi_path = str(proj / 'psi_dq_table.npz')
        phase2_mtpa_calibrate(
            psi_path, vdc=args.vdc, imax=args.imax,
            n_speed=11, n_Is=21)
    else:
        run(n_id=args.n_id, n_iq=args.n_iq, theta_r=theta_r,
            vdc=args.vdc, imax=args.imax, project_path=args.project_path)
