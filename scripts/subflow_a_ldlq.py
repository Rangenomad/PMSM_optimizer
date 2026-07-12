"""Sub-flow A: Ld/Lq MAP + 主磁链 Φ
求解器: MagnetostaticXY (4_Partial_motor_MS2)
方法: 电流角扫描法, assign_current() 设置激励 (含 9 匝因子)

用法:
    python -c "from scripts.subflow_a_ldlq import run; run(rated_current=10, max_current=3, current_steps=6, angle_steps=7)"
"""

import sys
sys.stdout.reconfigure(encoding='utf-8')
from pathlib import Path
import numpy as np
from ansys.aedt.core import Maxwell2d

ROOT = Path(__file__).resolve().parent.parent
TEMPLATE = str(ROOT / 'references' / 'Prius_2D_Practice.aedt')

TURNS = 9  # Matrix NumberOfTurns per coil object


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


def run(rated_current=10, max_current=3.0, current_steps=6, angle_steps=7):
    I_rated = rated_current
    I_max = max_current * I_rated
    currents = np.linspace(I_max / current_steps, I_max, current_steps)
    angles_deg = np.linspace(0, 90, angle_steps)

    m2d = Maxwell2d(
        project=TEMPLATE, design='4_Partial_motor_MS2',
        solution_type='MagnetostaticXY',
        non_graphical=False, new_desktop=False, close_on_exit=False
    )

    results = []
    for Ia_mag in currents:
        for theta_deg in angles_deg:
            theta = np.radians(theta_deg)
            Id = Ia_mag * np.sin(theta)
            Iq = Ia_mag * np.cos(theta)

            # Transform Id/Iq → three-phase currents
            Ia, Ib, Ic = id_iq_to_abc(Id, Iq, theta=0)

            # Assign currents with turns factor (9 turns per coil)
            m2d.assign_current(assignment=['PhaseA1', 'PhaseA2'], current=TURNS * Ia, units='A')
            m2d.assign_current(assignment=['PhaseB1', 'PhaseB2'], current=TURNS * Ib, units='A')
            m2d.assign_current(assignment=['PhaseC1', 'PhaseC2'], current=TURNS * Ic, units='A')

            # Solve
            m2d.analyze('Setup1')

            # Extract flux linkages
            data = m2d.post.get_solution_data(
                expressions=[
                    'FluxLinkage(PhaseA1)', 'FluxLinkage(PhaseA2)',
                    'FluxLinkage(PhaseB1)', 'FluxLinkage(PhaseB2)',
                    'FluxLinkage(PhaseC1)', 'FluxLinkage(PhaseC2)',
                ],
                variations=m2d.post.get_solution_data_variation()
            )

            # Total phase flux linkage = turns × sum of coil flux linkages
            psi_a = TURNS * (data.data('FluxLinkage(PhaseA1)')[0] + data.data('FluxLinkage(PhaseA2)')[0])
            psi_b = TURNS * (data.data('FluxLinkage(PhaseB1)')[0] + data.data('FluxLinkage(PhaseB2)')[0])
            psi_c = TURNS * (data.data('FluxLinkage(PhaseC1)')[0] + data.data('FluxLinkage(PhaseC2)')[0])

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
