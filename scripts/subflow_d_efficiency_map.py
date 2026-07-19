"""Sub-flow D: 全域工作特性效率 MAP
求解器: TransientXY (5_Partial_motor_TR)
方法: (转速 × 扭矩) 二维网格扫描
      每个工作点的 Id/Iq 由 **MTPA 算法** 从目标扭矩生成,
      而非经验系数估算.
      复用同一 m2d 对象提升求解速度.
      FEA 提取 Moving1.Torque 单值平均值,
      Vs / M / PF / 损耗 全部基于 MTPA 的 Id/Iq 解析计算.
      (PyAEDT gRPC 在加载 Transient 下只返回单值时点)

用法:
    python -c "from scripts.subflow_d_efficiency_map import run; r=run(speed_steps=3, torque_steps=3)"
    python -c "from scripts.subflow_d_efficiency_map import run; r=run(speed_steps=8, torque_steps=8)"
"""

import sys
sys.stdout.reconfigure(encoding='utf-8')
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
_DEFAULT_TEMPLATE = str(ROOT / 'references' / 'Prius_2D_Practice.aedt')

# --- 电机参数 (默认值, 从 Sub-flow A/B/C 结果校准) ---
DEFAULT_Ld = 0.00035    # H
DEFAULT_Lq = 0.00080    # H
DEFAULT_Phi = 0.07      # Wb
DEFAULT_Rs = 0.0        # Ω


def _mtpa_for_torque(T_target, pole_pairs, Ld, Lq, Phi, imax):
    """MTPA 算法: 给定目标扭矩, 返回 (Id, Iq, Is, β) 使电流最小.

    IPM 扭矩方程:  T = 1.5·P·(Φ·Iq + (Ld-Lq)·Id·Iq)
    MTPA 轨迹:     Id = A - √(A² + Iq²)   其中 A = Φ / (2·(Lq-Ld))
    方法: 沿 MTPA 轨迹对 Iq 二分搜索, 找到扭矩 = T_target 的工作点.
    """
    if T_target <= 0:
        return 0.0, 0.0, 0.0, 0.0

    delta_L = Lq - Ld                          # > 0 for IPM
    A = Phi / (2 * delta_L) if delta_L > 1e-12 else 1e12

    def _torque_at_Iq(Iq):
        """沿 MTPA 轨迹的扭矩 (给定 Iq, Id = A - √(A²+Iq²))"""
        Id = A - np.sqrt(A * A + Iq * Iq) if delta_L > 1e-12 else 0.0
        return 1.5 * pole_pairs * (Phi * Iq - delta_L * Id * Iq)

    # MTPA 轨迹上可达到的最大扭矩 (Iq=Imax 时)
    T_max = _torque_at_Iq(imax)

    if T_target >= T_max:
        # 饱和: 工作在 Imax 下的 MTPA 最优角
        Iq = imax
        Id = A - np.sqrt(A * A + Iq * Iq) if delta_L > 1e-12 else 0.0
    else:
        # 二分搜索: 在 MTPA 轨迹上找 Iq 使扭矩 = T_target
        lo, hi = 0.0, imax
        for _ in range(60):
            mid = (lo + hi) / 2.0
            if _torque_at_Iq(mid) > T_target:
                hi = mid
            else:
                lo = mid
            if abs(_torque_at_Iq(mid) - T_target) < 1e-6:
                break
        Iq = (lo + hi) / 2.0
        Id = A - np.sqrt(A * A + Iq * Iq) if delta_L > 1e-12 else 0.0

    Is = np.sqrt(Id * Id + Iq * Iq)
    beta = np.arctan2(-Id, Iq)                 # MTPA 电流角 (正值 = 弱磁)
    return Id, Iq, Is, beta


def _calc_analytical(pole_pairs, speed, Id, Iq, T_avg,
                     vdc=300, Ld=DEFAULT_Ld, Lq=DEFAULT_Lq,
                     Phi=DEFAULT_Phi, Rs=DEFAULT_Rs):
    """解析电压方程: 基于 (Id, Iq) 计算 Vs / M / PF / 损耗 / 效率"""
    Imax = np.sqrt(Id**2 + Iq**2)
    omega_e = speed * np.pi / 30 * pole_pairs   # 电角速度 (rad/s)

    # 电压方程 (稳态, 忽略电阻压降)
    Vd = -omega_e * Lq * Iq
    Vq = omega_e * (Ld * Id + Phi)
    Vs = np.sqrt(Vd**2 + Vq**2)

    Vmax = vdc / np.sqrt(3)
    M = Vs / Vmax if Vmax > 0 else 0.0

    # 功率因数角 = 电压角 - 电流角
    phi_v = np.arctan2(Vd, Vq)
    phi_i = np.arctan2(Id, Iq)
    PF = np.cos(phi_v - phi_i)

    # 功率 / 损耗
    P_out = T_avg * speed * 2 * np.pi / 60        # W
    P_cu = 3 * (Imax / np.sqrt(2))**2 * Rs        # W
    P_loss = P_cu
    eta = P_out / (P_out + P_loss) * 100 if (P_out + P_loss) > 0 else 0.0

    return {
        'Id': Id, 'Iq': Iq, 'Is': Imax,
        'beta': np.degrees(np.arctan2(-Id, Iq)),  # MTPA 电流角(°)
        'Vd': Vd, 'Vq': Vq,
        'Vs': Vs / 1000 if Vs > 1000 else Vs,
        'PF': PF, 'M': M,
        'P_out': P_out / 1000,
        'P_cu': P_cu / 1000,
        'P_loss': P_loss / 1000,
        'eta': eta,
    }


def run(speed_min=1000, speed_max=8000, speed_steps=5, torque_steps=5,
        vdc=300, imax=250,
        Ld=DEFAULT_Ld, Lq=DEFAULT_Lq, Phi=DEFAULT_Phi, Rs=DEFAULT_Rs,
        project_path=None):
    from scripts.project_utils import get_template_path
    from scripts.param_guard import check_var
    from ansys.aedt.core import Maxwell2d

    if project_path:
        template = get_template_path(Path(project_path))
    else:
        template = _DEFAULT_TEMPLATE

    m2d = Maxwell2d(
        project=template, design='5_Partial_motor_TR',
        solution_type='TransientXY',
        non_graphical=False, new_desktop=False, close_on_exit=False
    )
    pole_pairs = float(m2d['Poles']) / 2

    speed_points = np.linspace(speed_min, speed_max, speed_steps)
    T_max_analytical = 1.5 * pole_pairs * (
        Phi * imax
        + (Lq - Ld) * (Phi/(2*(Lq-Ld)) - np.sqrt((Phi/(2*(Lq-Ld)))**2 + imax**2)) * imax
    ) if Lq > Ld else 1.5 * pole_pairs * Phi * imax
    T_max_analytical = max(100, T_max_analytical)

    torque_points = np.linspace(0, T_max_analytical, torque_steps)

    results = np.zeros((speed_steps, torque_steps), dtype=object)
    total = speed_steps * torque_steps
    count = 0

    print(f'效率 MAP 扫描: {speed_steps}×{torque_steps} = {total} 点')
    print(f'转速范围: {speed_min}~{speed_max} rpm')
    print(f'Vdc={vdc}V, Imax={imax}A')
    print(f'Ld={Ld*1000:.4f}mH, Lq={Lq*1000:.4f}mH, Phi={Phi:.4f}Wb, Rs={Rs:.4f}Ω')
    print(f'MTPA 最大扭矩: {T_max_analytical:.0f} Nm')
    print(f'{"=" * 70}')

    for i, n in enumerate(speed_points):
        for j in range(torque_steps):
            count += 1
            T_target = torque_points[j]

            # ── MTPA: 从目标扭矩生成 (Id, Iq, β) ──
            Id_mtpa, Iq_mtpa, Is_mtpa, beta_mtpa = _mtpa_for_torque(
                T_target, pole_pairs, Ld, Lq, Phi, imax
            )
            angle_deg = int(np.degrees(beta_mtpa))

            # 设置 FEA 参数
            check_var('Speed_rpm', 'D')
            m2d['Speed_rpm'] = f'{n}rpm'
            check_var('Imax', 'D')
            m2d['Imax'] = f'{Is_mtpa}A'
            check_var('Thet_deg', 'D')
            m2d['Thet_deg'] = str(-angle_deg)

            freq = n / 60 * pole_pairs
            setup = m2d.setups[0]
            setup.props['StopTime'] = f'{1/freq}s'
            setup.props['TimeStep'] = f'{1/(freq*50)}s'
            setup.update()

            print(f'[D] {count:3d}/{total} ({100*count//total:3d}%) '
                  f'n={n:5.0f} T_target={T_target:5.1f} '
                  f'Id={Id_mtpa:6.1f} Iq={Iq_mtpa:6.1f} '
                  f'β={np.degrees(beta_mtpa):5.1f}°',
                  end='', flush=True)

            m2d.analyze('Setup1')

            # ── FEA 扭矩 ──
            T_avg = 0.0
            data = m2d.post.get_solution_data_per_variation(
                expressions=['Moving1.Torque']
            )
            if data:
                tq = np.array(data.data_real('Moving1.Torque'))
                if len(tq) >= 1:
                    T_avg = float(tq[0])

            # ── 解析电参数 (基于 MTPA 的 Id/Iq) ──
            pt = _calc_analytical(
                pole_pairs, n, Id_mtpa, Iq_mtpa, T_avg,
                vdc=vdc, Ld=Ld, Lq=Lq, Phi=Phi, Rs=Rs
            )
            pt['T_avg'] = T_avg
            pt['T_target'] = T_target
            pt['beta'] = np.degrees(beta_mtpa)
            results[i, j] = pt

            tau_ok = abs(T_avg / T_target - 1) < 0.1 if T_target > 0 else True
            tau_mark = '✓' if tau_ok else 'Δ'
            print(f' → T={T_avg:6.1f} ({T_avg/T_target:.2f}){tau_mark} '
                  f'Vs={pt["Vs"]:.0f}V M={pt["M"]:.2f} '
                  f'PF={pt["PF"]:.3f} η={pt["eta"]:.1f}%')

    m2d.close_project()

    # ============ 输出 MAP ============
    print(f'\n{"=" * 70}')
    print(f'=== 效率 MAP (Vdc={vdc}V, Imax={imax}A) ===')
    print(f'{"=" * 70}')

    T_labels = [f'{t:.0f}' for t in torque_points]

    def _print_map(title, key, fmt='8.2f'):
        print(f'\n--- {title} ---')
        print(f'{"n\\T":>6s}', end='')
        for tl in T_labels:
            print(f'{tl:>8s}', end='')
        print()
        for i, n in enumerate(speed_points):
            print(f'{n:6.0f}', end='')
            for j in range(torque_steps):
                v = results[i, j][key] if results[i, j] else 0
                print(f'{v:{fmt}}', end='')
            print()

    _print_map('FEA 扭矩 T_avg (Nm)', 'T_avg', '8.1f')
    _print_map('效率 η (%)', 'eta', '8.1f')
    _print_map('调制比 M', 'M', '8.3f')
    _print_map('功率因数 PF', 'PF', '8.3f')
    _print_map('相电压 Vs (V)', 'Vs', '8.0f')
    _print_map('电流 Is (A)', 'Is', '8.1f')
    _print_map('电流角 β (°)', 'beta', '8.1f')

    eta_vals = np.array([r['eta'] for r in results.flat if r])
    if len(eta_vals):
        print(f'\n峰值效率: {np.max(eta_vals):.1f}%')

    # 保存 CSV(带行列标签)
    out_dir = ROOT / 'results'
    out_dir.mkdir(exist_ok=True)
    col_labels = ',' + ','.join(T_labels)

    for fname, (hdr, key) in {
        'torque_map': ('T_avg (Nm)', 'T_avg'),
        'efficiency_map': ('eta (%)', 'eta'),
        'modulation_map': ('M', 'M'),
        'pf_map': ('PF', 'PF'),
        'voltage_map': ('Vs (V)', 'Vs'),
        'current_map': ('Is (A)', 'Is'),
        'beta_map': ('beta (deg)', 'beta'),
    }.items():
        rows = [f'# {hdr}  ─ 行:转速(rpm), 列:T_target(Nm)']
        rows.append(col_labels)
        for i, n in enumerate(speed_points):
            vals = [str(results[i, j][key]) if results[i, j] else '0'
                    for j in range(torque_steps)]
            rows.append(f'{n:.0f},{",".join(vals)}')
        (out_dir / f'{fname}.csv').write_text('\n'.join(rows), encoding='utf-8')

    print(f'\n结果已保存到 results/*.csv')
    return results


if __name__ == '__main__':
    run()
