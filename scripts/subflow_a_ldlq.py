"""Sub-flow A: Ld/Lq MAP + 主磁链 Φ
求解器: MagnetostaticXY (4_Partial_motor_MS2)
方法: 电流角扫描法, 通过 boundary.update() 修改已有 Current 激励 (含 9 匝因子)
      通过 export_matrix() 提取每线圈 Flux Linkage

用法:
    python -c "from scripts.subflow_a_ldlq import run; run(rated_current=250, current_steps=3, angle_steps=3)"
"""

import sys, os, re, tempfile
sys.stdout.reconfigure(encoding='utf-8')
from pathlib import Path
import numpy as np
from ansys.aedt.core import Maxwell2d

ROOT = Path(__file__).resolve().parent.parent
TEMPLATE = str(ROOT / 'references' / 'Prius_2D_Practice.aedt')

TURNS = 9  # Matrix NumberOfTurns per coil object

# 线圈配置: (boundary_name, IsPositive)
COIL_CONFIG = [
    ('PhaseA1', True), ('PhaseA2', True),
    ('PhaseB1', True), ('PhaseB2', True),
    ('PhaseC1', True), ('PhaseC2', True),
]
COIL_NAMES = [c[0] for c in COIL_CONFIG]


def id_iq_to_abc(Id, Iq, theta=0):
    """Inverse Park + Clarke: dq → abc (theta in radians)"""
    # θ=0: d-axis aligned with Phase A
    Ialpha = Id * np.cos(theta) - Iq * np.sin(theta)
    Ibeta  = Id * np.sin(theta) + Iq * np.cos(theta)
    Ia = Ialpha
    Ib = -0.5 * Ialpha + np.sqrt(3) / 2 * Ibeta
    Ic = -0.5 * Ialpha - np.sqrt(3) / 2 * Ibeta
    return Ia, Ib, Ic


def park_abc_to_dq(psi_a, psi_b, psi_c, theta=0):
    """Park transform: abc → dq (theta in radians)"""
    c = np.cos(theta)
    s = np.sin(theta)
    psi_d = 2 / 3 * (psi_a * c + psi_b * np.cos(theta - 2 * np.pi / 3) + psi_c * np.cos(theta + 2 * np.pi / 3))
    psi_q = -2 / 3 * (psi_a * s + psi_b * np.sin(theta - 2 * np.pi / 3) + psi_c * np.sin(theta + 2 * np.pi / 3))
    return psi_d, psi_q


def _set_coil_currents(m2d, boundaries, Ia, Ib, Ic):
    """Set coil currents via boundary.update() (含 9 匝因子)"""
    currents = [TURNS * Ia, TURNS * Ia,
                TURNS * Ib, TURNS * Ib,
                TURNS * Ic, TURNS * Ic]
    for bnd, (name, is_pos), val in zip(boundaries, COIL_CONFIG, currents):
        bnd.props['Current'] = f'{val}A'
        bnd.props['IsPositive'] = is_pos
        bnd.update()


def _parse_matrix_flux_linkage(text):
    """Parse matrix export .txt → dict {coil_name: flux_Wb}"""
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
    """Solve → export matrix → parse flux linkage per coil"""
    m2d.analyze('Setup1')
    out_file = os.path.join(tempfile.gettempdir(), 'pmsm_matrix_export.txt')
    m2d.export_matrix('Matrix1', out_file)
    with open(out_file, 'r') as f:
        text = f.read()
    return _parse_matrix_flux_linkage(text)


def run(rated_current=250, max_current=3.0, current_steps=6, angle_steps=7):
    I_rated = rated_current
    I_max = max_current * I_rated
    currents = np.linspace(I_max / current_steps, I_max, current_steps)
    angles_deg = np.linspace(0, 90, angle_steps)

    m2d = Maxwell2d(
        project=TEMPLATE, design='4_Partial_motor_MS2',
        solution_type='MagnetostaticXY',
        non_graphical=False, new_desktop=False, close_on_exit=False
    )

    # Get references to all 6 coil Current boundaries (once, before loop)
    boundaries = []
    for name, _ in COIL_CONFIG:
        b = [b for b in m2d.boundaries if b.name == name][0]
        boundaries.append(b)

    results = []
    total_points = len(currents) * len(angles_deg)
    point_num = 0
    for i, Ia_mag in enumerate(currents):
        for theta_deg in angles_deg:
            point_num += 1
            theta = np.radians(theta_deg)
            Id = Ia_mag * np.sin(theta)
            Iq = Ia_mag * np.cos(theta)

            # Transform Id/Iq → three-phase currents
            Ia, Ib, Ic = id_iq_to_abc(Id, Iq, theta=0)

            print(f'  [A] {point_num}/{total_points} ({100*point_num//total_points}%) '
                  f'求解 Id={Id:+7.2f}A, Iq={Iq:+7.2f}A, Ia={Ia:.2f}A ...')

            # Set coil currents via boundary.update() (works with gRPC)
            _set_coil_currents(m2d, boundaries, Ia, Ib, Ic)

            # Solve + export matrix → extract per-coil flux linkages
            psi_coil = _extract_flux_linkages(m2d)

            # Total phase flux linkage = turns × sum of coil flux linkages
            psi_a = TURNS * (psi_coil['PhaseA1'] + psi_coil['PhaseA2'])
            psi_b = TURNS * (psi_coil['PhaseB1'] + psi_coil['PhaseB2'])
            psi_c = TURNS * (psi_coil['PhaseC1'] + psi_coil['PhaseC2'])

            # Park transform → dq flux linkages
            psi_d, psi_q = park_abc_to_dq(psi_a, psi_b, psi_c, theta=0)

            # Inductances (avoid division by zero)
            Ld = psi_d / Id if abs(Id) > 1e-6 else 0
            Lq = psi_q / Iq if abs(Iq) > 1e-6 else 0

            results.append({
                'Id': Id, 'Iq': Iq,
                'Psi_d': psi_d, 'Psi_q': psi_q,
                'Ld': Ld, 'Lq': Lq,
            })

            print(f'  Id={Id:+7.2f}  Iq={Iq:+7.2f}  Ld={Ld:.5f}  Lq={Lq:.5f}')

    # 主磁链 Φ: Id=0, Iq≈0 时的 Ψq（永磁体贡献）
    near_zero = [r for r in results if abs(r['Id']) < 1e-6 and abs(r['Iq']) < 1e-6]
    phi = near_zero[0]['Psi_q'] if near_zero else results[0]['Psi_q']

    # 输出表格
    print(f'\n{"=" * 60}')
    print(f'主磁链 Φ = {phi:.6f} Wb')
    print(f'{"=" * 60}')
    print(f'{"Id\\Iq":>8s}', end='')
    for a in angles_deg:
        print(f'  {a:7.0f}°    ', end='')
    print()
    for i, Ia_mag in enumerate(currents):
        print(f'{Ia_mag:8.2f}', end='')
        for j in range(len(angles_deg)):
            idx = i * len(angles_deg) + j
            print(f'  Ld={results[idx]["Ld"]:.4f}', end='')
        print()
        print(f'{"":>8s}', end='')
        for j in range(len(angles_deg)):
            idx = i * len(angles_deg) + j
            print(f'  Lq={results[idx]["Lq"]:.4f}', end='')
        print()

    m2d.close_project()
    return results, phi


if __name__ == '__main__':
    run()
