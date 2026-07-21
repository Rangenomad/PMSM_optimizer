"""Sub-flow E: 外特性计算 (T-n 曲线)
方法: 纯数学计算 (不跑 FEA)
原理: 电压圆 Vs ≤ Vdc/√3, 电流圆 Is ≤ Imax 约束下, 基于 Ld/Lq + Φ 模型计算

支持三种饱和模型:
  - constant:   Ld = Ld0 (固定值)
  - linear:     Ld(Id) = Ld0 + k1*|Id|
  - quadratic:  Ld(Id) = Ld0 + k1*|Id| + k2*|Id|²  (默认, FEA 标定)

用法:
    # 默认二次饱和模型 (推荐)
    python -c "from scripts.subflow_e_external import run; run(vdc=300, imax=153)"
    # 350A 全范围
    python -c "from scripts.subflow_e_external import run; run(vdc=300, imax=350, speed_max=12000)"
"""

import sys
sys.stdout.reconfigure(encoding='utf-8')
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent

# ── FEA 标定的默认参数 (二次饱和模型, @350A 反凸极电机) ──
_DEFAULT_Ld0 = 0.002056    # H — 二次截距 (Ld at Id=0)
_DEFAULT_Lq  = 0.001550    # H — MS 锚点
_DEFAULT_Phi = 0.098857    # Wb — FEA 扭矩标定 (T@0deg=207.6 Nm)
_DEFAULT_Rs  = 0.0         # Ω
_DEFAULT_k1  = -2.09e-6    # H/A — 线性饱和系数
_DEFAULT_k2  = 1.03e-8     # H/A² — 二次饱和系数
_DEFAULT_P   = 4           # 极对数


def run(vdc=300, imax=250, speed_max=12000, rs=None, speed_points_n=40,
        pole_pairs=None, phi=None, ld_lq_data=None, project_path=None,
        ld0=None, lq=None, k1=None, k2=None, model='quadratic',
        mtpa_table_path=None):
    """计算外特性 T-n 曲线。

    两种模式:
      1. MTPA table (推荐): 传入 mtpa_table_path, 基于 ψd/ψq 磁链表查找
      2. Ld/Lq 模型 (向后兼容): 传入 ld0/lq/k1/k2 或使用默认值

    Parameters
    ----------
    vdc : float — 直流母线电压 (V)
    imax : float — 最大相电流幅值 (A)
    speed_max : float — 最高机械转速 (rpm)
    rs : float, optional — 相电阻 (Ω)
    speed_points_n : int — 转速分点数
    pole_pairs : int, optional — 极对数, 默认 4
    mtpa_table_path : str, optional — MTPA 表路径 (Phase 2 输出)
    project_path : str, optional — 项目目录
    (以下为 Ld/Lq 模型参数, mtpa_table_path 传入时忽略)
    phi, ld0, lq, k1, k2, model, ld_lq_data
    """
    # ── MTPA table mode ──
    if mtpa_table_path:
        return _run_mtpa(vdc=vdc, imax=imax, speed_max=speed_max,
                         speed_points_n=speed_points_n, pole_pairs=pole_pairs,
                         rs=rs, mtpa_table_path=mtpa_table_path,
                         project_path=project_path)

    # ── 参数解析 (Ld/Lq 模型) ──
    P = pole_pairs if pole_pairs is not None else _DEFAULT_P
    Rs = rs if rs is not None else _DEFAULT_Rs
    Lq = lq if lq is not None else _DEFAULT_Lq
    Ld0 = ld0 if ld0 is not None else _DEFAULT_Ld0
    Phi = phi if phi is not None else _DEFAULT_Phi
    K1 = k1 if k1 is not None else _DEFAULT_k1
    K2 = k2 if k2 is not None else _DEFAULT_k2

    # 从 Sub-flow A 数据覆盖 (向后兼容, 使用平均值)
    if ld_lq_data is not None:
        ld_vals = np.array([d['Ld'] for d in ld_lq_data])
        lq_vals = np.array([d['Lq'] for d in ld_lq_data])
        Ld0 = np.mean(ld_vals[ld_vals > 0]) if np.any(ld_vals > 0) else Ld0
        Lq = np.mean(lq_vals[lq_vals > 0]) if np.any(lq_vals > 0) else Lq

    Vmax = vdc / np.sqrt(3)  # SVPWM 线性调制区上限

    # ── 饱和模型 ──
    if model == 'quadratic':
        def Ld_at(Id_val):
            return Ld0 + K1 * abs(Id_val) + K2 * Id_val**2
    elif model == 'linear':
        def Ld_at(Id_val):
            return Ld0 + K1 * abs(Id_val)
    else:  # constant
        def Ld_at(Id_val):
            return Ld0

    # ── 打印参数 ──
    dL_sign = 'Ld > Lq (反凸极)' if (Ld0 - Lq) > 0 else 'Ld < Lq (传统 IPM)'
    print(f'\n{"=" * 70}')
    print(f'  Sub-flow E: 外特性 T-n 曲线')
    print(f'{"=" * 70}')
    print(f'  Vdc    = {vdc:.0f} V        Vmax(SVPWM) = {Vmax:.1f} V')
    print(f'  Imax   = {imax:.0f} A        Rs = {Rs:.4f} Ω')
    print(f'  P      = {P}           转速范围: 1 ~ {speed_max:.0f} rpm')
    print(f'  Phi    = {Phi*1000:.4f} mWb    Lq = {Lq*1000:.4f} mH')
    print(f'  Ld0    = {Ld0*1000:.4f} mH     ({dL_sign})')
    if model != 'constant':
        print(f'  k1     = {K1:.2e} H/A')
    if model == 'quadratic':
        print(f'  k2     = {K2:.2e} H/A²')
    print(f'  饱和模型: {model}')
    print(f'{"=" * 70}')

    # ── 转速范围 ──
    speeds = np.geomspace(max(1, speed_max / speed_points_n), speed_max, speed_points_n)

    # ── 2D 搜索网格: (Id, Iq) 覆盖整个电流圆内部 ──
    N_BETA = 360     # 角度分辨率
    N_IS = 200       # 电流幅值分辨率

    # 预计算网格 (一次, 所有转速共用)
    beta_grid = np.linspace(-np.pi/2, np.pi/2, N_BETA)  # rad
    Is_grid = np.linspace(0, imax, N_IS)                 # A

    # beta × Is 网格 → Id, Iq 矩阵
    Beta_mesh, Is_mesh = np.meshgrid(beta_grid, Is_grid, indexing='ij')
    Id_all = Is_mesh * np.sin(Beta_mesh)   # (N_BETA, N_IS)
    Iq_all = Is_mesh * np.cos(Beta_mesh)   # (N_BETA, N_IS)

    # 预计算 Ld 矩阵 (依赖 Id)
    if model == 'quadratic':
        Ld_all = Ld0 + K1 * np.abs(Id_all) + K2 * Id_all**2
    elif model == 'linear':
        Ld_all = Ld0 + K1 * np.abs(Id_all)
    else:
        Ld_all = np.full_like(Id_all, Ld0)

    dL_all = Ld_all - Lq  # (N_BETA, N_IS)

    # 扭矩系数 (不含 Iq 的部分): T = 1.5*P * (Phi*Iq + dL*Id*Iq)
    # 预计算 torque_base = 1.5*P*(Phi + dL*Id) * Iq

    results = []
    for n in speeds:
        omega_e = 2 * np.pi * n / 60 * P  # 电角速度 (rad/s)

        # ── 电压约束: 全网格计算 ──
        Vd = -omega_e * Lq * Iq_all + Rs * Id_all
        Vq = omega_e * (Ld_all * Id_all + Phi) + Rs * Iq_all
        Vs = np.sqrt(Vd**2 + Vq**2)

        # 可行 mask
        feasible = Vs <= Vmax

        if not np.any(feasible):
            # 完全不可行 — 返回零扭矩
            results.append({
                'speed': n, 'torque': 0.0, 'power': 0.0,
                'Id': 0.0, 'Iq': 0.0, 'beta_deg': 0.0,
                'modulation': 0.0, 'Is': 0.0, 'region': 'FW',
            })
            continue

        # ── 扭矩: 仅计算可行点 ──
        T_all = np.full_like(Id_all, -1e12)
        T_feasible = 1.5 * P * (Phi * Iq_all + dL_all * Id_all * Iq_all)
        T_all[feasible] = T_feasible[feasible]

        # 找最大扭矩点
        idx_flat = np.argmax(T_all)
        idx_beta, idx_Is = np.unravel_index(idx_flat, Id_all.shape)

        best_T = T_all[idx_beta, idx_Is]
        best_Id = Id_all[idx_beta, idx_Is]
        best_Iq = Iq_all[idx_beta, idx_Is]
        best_Is = Is_grid[idx_Is]
        best_beta = beta_grid[idx_beta]
        Ld_op = Ld_all[idx_beta, idx_Is]

        # 区域判定
        region = 'MTPA' if abs(best_Is - imax) < 1e-3 else 'FW'

        # 输出功率
        P_out = best_T * n * 2 * np.pi / 60 / 1000  # kW

        # 调制比
        Ld_op = Ld_at(best_Id)
        Vd_op = -omega_e * Lq * best_Iq + Rs * best_Id
        Vq_op = omega_e * (Ld_op * best_Id + Phi) + Rs * best_Iq
        Vs_op = np.sqrt(Vd_op**2 + Vq_op**2)
        M = Vs_op / Vmax

        results.append({
            'speed': n,
            'torque': best_T,
            'power': P_out,
            'Id': best_Id,
            'Iq': best_Iq,
            'beta_deg': np.degrees(best_beta),
            'modulation': M,
            'Is': best_Is,
            'region': region,
        })

    # ── 输出表格 ──
    print(f'\n{"转速":>8s} {"扭矩":>8s} {"功率":>7s} {"Id":>7s} {"Iq":>7s} {"Is":>7s} {"β":>6s} {"M":>7s} {"区域":>5s}')
    print(f'{"rpm":>8s} {"Nm":>8s} {"kW":>7s} {"A":>7s} {"A":>7s} {"A":>7s} {"°":>6s} {"":>7s} {"":>5s}')
    print('-' * 72)

    turn_point_reported = False
    for r in results:
        print(f'{r["speed"]:8.0f} {r["torque"]:8.1f} {r["power"]:7.2f} '
              f'{r["Id"]:7.2f} {r["Iq"]:7.2f} {r["Is"]:7.1f} {r["beta_deg"]:6.1f} '
              f'{r["modulation"]:7.3f} {r["region"]:>5s}')

        if not turn_point_reported and (r['region'] == 'FW' or r['modulation'] > 0.95):
            print(f'\n  ▲ 基速/转折点: {r["speed"]:.0f} rpm (进入弱磁区)')
            turn_point_reported = True

    # ── 统计 ──
    peak_power = max(results, key=lambda r: r['power'])
    max_torque = max(results, key=lambda r: r['torque'])
    fw_entries = [r for r in results if r['region'] == 'FW']
    fw_start = fw_entries[0]['speed'] if fw_entries else None

    print(f'\n{"=" * 70}')
    print(f'  峰值扭矩: {max_torque["torque"]:.1f} Nm @ {max_torque["speed"]:.0f} rpm')
    print(f'  峰值功率: {peak_power["power"]:.2f} kW @ {peak_power["speed"]:.0f} rpm')
    if fw_start:
        print(f'  基速 (转折点): {fw_start:.0f} rpm')
        last = results[-1]
        print(f'  最高转速扭矩: {last["torque"]:.1f} Nm @ {last["speed"]:.0f} rpm')
    print(f'{"=" * 70}')

    return results


def _run_mtpa(vdc=450, imax=350, speed_max=8000, speed_points_n=40,
              pole_pairs=None, rs=None, mtpa_table_path=None, project_path=None):
    """MTPA-table-based external characteristic (T-n curve).

    Pure table lookup — reads mtpa_table.npz, applies user Vdc/Imax,
    outputs T-n curve with boundary classification.
    """
    P = pole_pairs if pole_pairs is not None else _DEFAULT_P
    Rs = rs if rs is not None else _DEFAULT_Rs
    Vmax = vdc / np.sqrt(3)

    # Resolve mtpa_table path
    if mtpa_table_path:
        table_path = Path(mtpa_table_path)
    elif project_path:
        table_path = Path(project_path) / 'mtpa_table.npz'
    else:
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

    print(f'\n{"=" * 70}')
    print(f'  Sub-flow E: 外特性 T-n 曲线 (MTPA table lookup)')
    print(f'{"=" * 70}')
    print(f'  MTPA table:  {table_path}')
    print(f'  Vdc={vdc:.0f}V  Vmax={Vmax:.1f}V  Imax={imax:.0f}A  P={P}')
    print(f'  Table: {len(speed_grid)} speeds × {len(Is_grid)} Is')
    print(f'{"=" * 70}')

    # Generate output speeds
    out_speeds = sorted(set(
        np.linspace(min(speed_grid), speed_max, speed_points_n)
    ))

    results = []
    for n in out_speeds:
        i_speed = np.argmin(np.abs(speed_grid - n))

        valid_T = []
        valid_idx = []

        for j, Is_val in enumerate(Is_grid):
            if not feasible[i_speed, j]:
                continue
            if Is_val > imax:
                continue
            if Vs[i_speed, j] > Vmax:
                continue
            valid_T.append(T_mtpa[i_speed, j])
            valid_idx.append(j)

        if not valid_T:
            results.append({
                'speed': n, 'torque': 0.0, 'power': 0.0,
                'Id': 0.0, 'Iq': 0.0, 'beta_deg': 0.0,
                'modulation': 0.0, 'Is': 0.0, 'region': 'none',
            })
            continue

        best_j = valid_idx[np.argmax(valid_T)]
        T_max = float(T_mtpa[i_speed, best_j])
        Is_op = float(Is_grid[best_j])
        Thet_op = float(data['Thet_opt'][i_speed, best_j])
        Vs_op = float(Vs[i_speed, best_j])
        Id_op = Is_op * np.sin(np.radians(Thet_op))
        Iq_op = Is_op * np.cos(np.radians(Thet_op))

        # Determine limit type
        if abs(Is_op - imax) < 1.0 or abs(Is_op - max(Is_grid[valid_idx])) < 1.0:
            limit_type = 'current'
        else:
            limit_type = 'voltage'

        P_out = T_max * n * 2 * np.pi / 60 / 1000

        results.append({
            'speed': n, 'torque': T_max, 'power': P_out,
            'Id': Id_op, 'Iq': Iq_op, 'beta_deg': Thet_op,
            'modulation': Vs_op / Vmax, 'Is': Is_op,
            'region': limit_type,
        })

    # ── Output ──
    print(f'\n{"n_rpm":>8s} {"T_Nm":>8s} {"kW":>7s} {"Id_A":>7s} {"Iq_A":>7s} '
          f'{"Is_A":>7s} {"Thet":>6s} {"M":>6s} {"limit":>8s}')
    print('-' * 76)
    for r in results:
        print(f'{r["speed"]:8.0f} {r["torque"]:8.1f} {r["power"]:7.2f} '
              f'{r["Id"]:7.1f} {r["Iq"]:7.1f} {r["Is"]:7.1f} '
              f'{r["beta_deg"]:6.1f} {r["modulation"]:6.3f} {r["region"]:>8s}')

    # Save CSV
    out_dir = Path(project_path) if project_path else table_path.parent
    out_path = out_dir / 'external_characteristic.csv'
    with open(out_path, 'w', encoding='utf-8') as f:
        f.write('n_rpm,T_max_Nm,P_kW,Id_A,Iq_A,Thet_deg,Is_A,Vs_V,M,limit_type\n')
        for r in results:
            f.write(f'{r["speed"]:.1f},{r["torque"]:.2f},{r["power"]:.2f},'
                    f'{r["Id"]:.2f},{r["Iq"]:.2f},{r["beta_deg"]:.1f},'
                    f'{r["Is"]:.2f},{r["modulation"]:.4f},{r["region"]}\n')
    print(f'\n  Saved: {out_path}')

    # Stats
    valid = [r for r in results if r['torque'] > 0.1]
    if valid:
        peak = max(valid, key=lambda r: r['torque'])
        fw = [r for r in valid if r['region'] == 'voltage']
        print(f'\n  峰值扭矩: {peak["torque"]:.1f} Nm @ {peak["speed"]:.0f} rpm')
        peak_pwr = max(valid, key=lambda r: r['power'])
        print(f'  峰值功率: {peak_pwr["power"]:.2f} kW @ {peak_pwr["speed"]:.0f} rpm')
        if fw:
            print(f'  基速/转折点: {fw[0]["speed"]:.0f} rpm')

    return results


if __name__ == '__main__':
    run()
