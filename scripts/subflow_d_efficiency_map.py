"""Sub-flow D: 全域工作特性效率 MAP
求解器: TransientXY (5_Partial_motor_TR)
方法: (转速 × 扭矩) 二维网格扫描
      复用同一 m2d 对象提升求解速度 (避免每次开/关工程重创建网格)
      每个工作点: FEA 提取 Moving1.Torque 单值平均值
      电压 Vs / 调制比 M / 功率因数 PF 采用解析电压方程
      (PyAEDT gRPC 在加载 Transient 下只返回单值时点)

用法:
    # 快速测试
    python -c "from scripts.subflow_d_efficiency_map import run; r=run(speed_steps=3, torque_steps=3)"
    # 完整扫描
    python -c "from scripts.subflow_d_efficiency_map import run; r=run(speed_steps=8, torque_steps=8)"
"""

import sys
sys.stdout.reconfigure(encoding='utf-8')
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
TEMPLATE = str(ROOT / 'references' / 'Prius_2D_Practice.aedt')

# --- 电机参数 (默认值, 可从 Sub-flow A/B/C 结果校准) ---
DEFAULT_Ld = 0.00035    # H
DEFAULT_Lq = 0.00080    # H
DEFAULT_Phi = 0.07      # Wb
DEFAULT_Rs = 0.0        # Ω


def _calc_analytical(speed, imax, angle_deg, pole_pairs, T_avg,
                     vdc=300, Ld=DEFAULT_Ld, Lq=DEFAULT_Lq,
                     Phi=DEFAULT_Phi, Rs=DEFAULT_Rs):
    """基于解析电压方程计算电参数, 无需 FEA 波形数据"""
    angle_rad = np.radians(angle_deg)
    Id = -imax * np.sin(angle_rad)     # 电机惯例: d 轴去磁
    Iq = imax * np.cos(angle_rad)

    omega_e = speed * np.pi / 30 * pole_pairs   # 电角速度 (rad/s)

    # 电压方程 (稳态, 忽略电阻压降)
    Vd = -omega_e * Lq * Iq
    Vq = omega_e * (Ld * Id + Phi)
    Vs = np.sqrt(Vd**2 + Vq**2)

    Vmax = vdc / np.sqrt(3)
    M = Vs / Vmax if Vmax > 0 else 0.0

    # 功率因数
    phi_v = np.arctan2(Vd, Vq)
    phi_i = np.arctan2(Id, Iq)
    PF = np.cos(phi_v - phi_i)

    # 功率 / 损耗
    P_out = T_avg * speed * 2 * np.pi / 60        # W
    P_cu = 3 * (imax / np.sqrt(2))**2 * Rs        # W
    P_loss = P_cu   # 可扩展: + P_fe + P_mech
    eta = P_out / (P_out + P_loss) * 100 if (P_out + P_loss) > 0 else 0.0

    return {
        'Id': Id, 'Iq': Iq, 'Vd': Vd, 'Vq': Vq,
        'Vs': Vs / 1000 if Vs > 1000 else Vs,
        'PF': PF, 'M': M,
        'P_out': P_out / 1000, 'P_cu': P_cu / 1000,
        'P_loss': P_loss / 1000, 'eta': eta,
    }


def run(speed_min=1000, speed_max=8000, speed_steps=5, torque_steps=5,
        vdc=300, imax=250):
    from ansys.aedt.core import Maxwell2d

    m2d = Maxwell2d(
        project=TEMPLATE, design='5_Partial_motor_TR',
        solution_type='TransientXY',
        non_graphical=False, new_desktop=False, close_on_exit=False
    )
    pole_pairs = float(m2d['Poles']) / 2

    speed_points = np.linspace(speed_min, speed_max, speed_steps)
    torque_max = 250  # Nm (估算上限)

    results = np.zeros((speed_steps, torque_steps), dtype=object)
    total_points = speed_steps * torque_steps
    count = 0

    print(f'效率 MAP 扫描: {speed_steps}×{torque_steps} = {total_points} 点')
    print(f'转速范围: {speed_min}~{speed_max} rpm')
    print(f'Vdc={vdc}V, Imax={imax}A')
    print(f'Ld={DEFAULT_Ld*1000:.4f}mH, Lq={DEFAULT_Lq*1000:.4f}mH, Phi={DEFAULT_Phi:.4f}Wb')
    print(f'{"=" * 70}')

    for i, n in enumerate(speed_points):
        for j in range(torque_steps):
            count += 1
            T_target = torque_max * j / (torque_steps - 1) if torque_steps > 1 else 0

            # 从目标扭矩估算 Id/Iq (MTPA 近似)
            T_ratio = T_target / torque_max if torque_max > 0 else 0
            Id_est = -0.5 * T_ratio * imax
            Iq_est = 0.8 * (1 - 0.2 * T_ratio) * imax
            Is_est = np.sqrt(Id_est**2 + Iq_est**2)
            if Is_est > imax:
                Is_est = imax
                ratio = imax / Is_est
                Id_est *= ratio
                Iq_est *= ratio
            angle_deg = int(np.degrees(np.arctan2(Id_est, Iq_est)))

            # 设置参数
            m2d['Speed_rpm'] = f'{n}rpm'
            m2d['Imax'] = f'{Is_est}A'
            m2d['Thet_deg'] = str(angle_deg)

            freq = n / 60 * pole_pairs
            setup = m2d.setups[0]
            setup.props['StopTime'] = f'{1/freq}s'
            setup.props['TimeStep'] = f'{1/(freq*200)}s'
            setup.update()

            print(f'[D] {count:3d}/{total_points} ({100*count//total_points:3d}%) '
                  f'n={n:5.0f} T={T_target:5.1f} Is={Is_est:5.1f} θ={angle_deg:3d}°',
                  end='', flush=True)

            m2d.analyze('Setup1')

            # --- FEA 扭矩 ---
            T_avg = 0.0
            data = m2d.post.get_solution_data_per_variation(
                expressions=['Moving1.Torque']
            )
            if data:
                tq = np.array(data.data_real('Moving1.Torque'))
                if len(tq) >= 1:
                    T_avg = float(tq[0])

            # --- 解析电参数 ---
            pt = _calc_analytical(n, Is_est, angle_deg, pole_pairs, T_avg,
                                 vdc=vdc, Ld=DEFAULT_Ld, Lq=DEFAULT_Lq, Phi=DEFAULT_Phi)
            pt['T_avg'] = T_avg
            pt['Is'] = Is_est
            results[i, j] = pt

            tau_print = '?' if T_target == 0 else f'{T_avg/T_target:.2f}'
            print(f' → T={T_avg:6.1f} ({tau_print}) Vs={pt["Vs"]:.0f}V '
                  f'M={pt["M"]:.2f} PF={pt["PF"]:.3f} η={pt["eta"]:.1f}%')

    m2d.close_project()

    # ============ 输出 MAP ============
    print(f'\n{"=" * 70}')
    print(f'=== 效率 MAP (Vdc={vdc}V, Imax={imax}A) ===')
    print(f'{"=" * 70}')

    T_labels = [f'{torque_max*j/(torque_steps-1):.0f}'
                for j in range(torque_steps)]

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

    _print_map('扭矩 T_avg (Nm)', 'T_avg', '8.1f')
    _print_map('效率 η (%)', 'eta', '8.1f')
    _print_map('调制比 M', 'M', '8.3f')
    _print_map('功率因数 PF', 'PF', '8.3f')
    _print_map('相电压 Vs (V)', 'Vs', '8.0f')
    _print_map('电流 Is (A)', 'Is', '8.1f')

    # 峰值效率
    eta_vals = np.array([r['eta'] for r in results.flat if r])
    if len(eta_vals):
        print(f'\n峰值效率: {np.max(eta_vals):.1f}%')

    # 保存 CSV
    out_dir = ROOT / 'results'
    out_dir.mkdir(exist_ok=True)
    for fname, (hdr, key) in {
        'torque_map': ('T_avg (Nm)', 'T_avg'),
        'efficiency_map': ('η (%)', 'eta'),
        'modulation_map': ('M', 'M'),
        'pf_map': ('PF', 'PF'),
        'voltage_map': ('Vs (V)', 'Vs'),
        'current_map': ('Is (A)', 'Is'),
    }.items():
        mat = np.array([[r[key] for r in row] for row in results])
        np.savetxt(out_dir / f'{fname}.csv', mat, delimiter=',',
                   header=hdr, fmt='%.3f')

    print(f'\n结果已保存到 results/*.csv')
    return results


if __name__ == '__main__':
    run()
