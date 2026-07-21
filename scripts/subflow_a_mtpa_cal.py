"""Sub-flow A: MTPA Calibration
Phase 1: ψd/ψq 2D flux-linkage table (Magnetostatic, Id×Iq grid)
Phase 2: MTPA calibration → (n, Is) → (Thet_opt, Id_opt, T_mtpa, Vs, M)

Usage:
    # Full pipeline (Phase 1 + Phase 2)
    python -c "from scripts.subflow_a_mtpa_cal import run; run(n_id=10, n_iq=10, project_path='pmsm_projects/xxx')"

    # Phase 1 only
    python -c "from scripts.subflow_a_mtpa_cal import phase1_psi_dq_table; phase1_psi_dq_table(n_id=10, n_iq=10, project_path='pmsm_projects/xxx')"

    # Phase 2 only (requires psi_dq_table.npz)
    python -c "from scripts.subflow_a_mtpa_cal import phase2_mtpa_calibrate; phase2_mtpa_calibrate('pmsm_projects/xxx/psi_dq_table.npz')"
"""

import sys, os, re, tempfile, json
sys.stdout.reconfigure(encoding='utf-8')
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
_DEFAULT_TEMPLATE = str(ROOT / 'references' / 'Prius_2D_Practice.aedt')

TURNS = 9
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
    id_range : tuple — (Id_min, Id_max) in A
    iq_range : tuple — (Iq_min, Iq_max) in A
    n_id : int — number of Id steps
    n_iq : int — number of Iq steps
    theta_r : float — rotor d-axis angle in radians (from θr calibration)
    project_path : str, optional — project directory path
    confirm : bool — if True, prompt user before running

    Returns
    -------
    id_grid : np.ndarray (n_id,)
    iq_grid : np.ndarray (n_iq,)
    psi_d : np.ndarray (n_id, n_iq)
    psi_q : np.ndarray (n_id, n_iq)
    """
    from scripts.project_utils import get_template_path
    from ansys.aedt.core import Maxwell2d

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
    result_p1 = phase1_psi_dq_table(
        id_range=id_range, iq_range=iq_range,
        n_id=n_id, n_iq=n_iq, theta_r=theta_r,
        project_path=project_path, confirm=True
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
    p.add_argument('--theta-r', type=float, default=0.0)
    p.add_argument('--vdc', type=float, default=450)
    p.add_argument('--imax', type=float, default=350)
    p.add_argument('--phase', type=int, default=0, help='1=Phase1 only, 2=Phase2 only, 0=both')
    p.add_argument('--psi-dq-path', default=None, help='Path to psi_dq_table.npz (Phase 2 only)')
    args = p.parse_args()

    if args.phase == 1:
        phase1_psi_dq_table(
            n_id=args.n_id, n_iq=args.n_iq,
            theta_r=args.theta_r, project_path=args.project_path)
    elif args.phase == 2:
        psi_path = args.psi_dq_path
        if not psi_path:
            proj = Path(args.project_path) if args.project_path else Path.cwd()
            psi_path = str(proj / 'psi_dq_table.npz')
        phase2_mtpa_calibrate(
            psi_path, vdc=args.vdc, imax=args.imax,
            n_speed=11, n_Is=21)
    else:
        run(n_id=args.n_id, n_iq=args.n_iq, theta_r=args.theta_r,
            vdc=args.vdc, imax=args.imax, project_path=args.project_path)
