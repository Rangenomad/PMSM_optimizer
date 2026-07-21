"""θr Calibration — find PM d-axis offset angle.

Zero-current Magnetostatic solve → extract ψa/ψb/ψc → sweep Park angle
to find θ where ψq=0 (d-axis aligned with rotor PM field).

Usage:
    python scripts/test_theta_r.py
    python scripts/test_theta_r.py pmsm_projects/2026-07-21_mtpa_test
"""

import sys
sys.stdout.reconfigure(encoding='utf-8')
import os
import tempfile
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
_DEFAULT_TEMPLATE = str(ROOT / 'references' / 'Prius_2D_Practice.aedt')

TURNS = 9
COIL_NAMES = ['PhaseA1', 'PhaseA2', 'PhaseB1', 'PhaseB2', 'PhaseC1', 'PhaseC2']


def park_abc_to_dq(psi_a, psi_b, psi_c, theta=0):
    """Park transform: abc → dq (theta in radians).

    Convention: q-axis LAGS d-axis by 90° (Maxwell motor convention).
    """
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
    boundaries = []
    for name in COIL_NAMES:
        b = [b for b in m2d.boundaries if b.name == name][0]
        boundaries.append(b)

    # Set all coil currents to 0 (preserve original IsPositive)
    for bnd in boundaries:
        bnd.props['Current'] = '0A'
        bnd.update()

    # Solve
    print('  Solving zero-current magnetostatic...')
    m2d.analyze('Setup1')

    # Export matrix → extract flux linkages
    out_file = os.path.join(tempfile.gettempdir(), 'pmsm_theta_r_test.txt')
    m2d.export_matrix('Matrix1', out_file)
    with open(out_file, 'r', encoding='utf-8') as f:
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
            if len(parts) >= 2 and parts[0] in COIL_NAMES:
                result[parts[0]] = float(parts[1])

    # Total phase flux = TURNS * (coil1 + coil2) — export_matrix is per-turn
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
