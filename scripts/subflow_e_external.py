"""Sub-flow E: 外特性计算 (T-n 曲线)
方法: 纯解析 (不跑 FEA)
原理: 基于 Sub-flow A 输出的 ψd/ψq 磁链表, 在电压圆 Vs≤Vdc/√3 和
      电流圆 Is≤Imax 约束下, MTPA + FW 控制求取各转速的最大扭矩.

主要模式:
  - ψd/ψq 磁链表 (推荐): 直接读取 psi_dq_table.npz, 对任意 Vdc/Imax 重算
  - Ld/Lq 参数模型 (向后兼容): 固定参数 + 饱和模型

用法:
    # 主管线: 基于磁链表 (推荐)
    python -c "from scripts.subflow_e_external import run; run(
        project_path='pmsm_projects/xxx', vdc=500, imax=300)"

    # 向后兼容: Ld/Lq 参数模型
    python -c "from scripts.subflow_e_external import run; run(vdc=300, imax=153)"
"""

import sys
sys.stdout.reconfigure(encoding='utf-8')
from pathlib import Path
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import minimize
import warnings
warnings.filterwarnings('ignore')

ROOT = Path(__file__).resolve().parent.parent

# ── Ld/Lq 参数模型默认值 (向后兼容) ──
_DEFAULT_Ld0 = 0.002056
_DEFAULT_Lq  = 0.001550
_DEFAULT_Phi = 0.098857
_DEFAULT_Rs  = 0.0
_DEFAULT_k1  = -2.09e-6
_DEFAULT_k2  = 1.03e-8
_DEFAULT_P   = 4


# ═══════════════════════════════════════════════════════════════
#  ψd/ψq 磁链表模式 (主管线推荐)
# ═══════════════════════════════════════════════════════════════

def _run_psidq(psi_dq_npz, vdc=500, imax=300, speed_max=10000,
               speed_points_n=150, pole_pairs=4, rs=0.05,
               out_dir=None, plot=True):
    """基于 ψd/ψq 磁链表的外特性计算.

    对每个转速, 在精细 (Id,Iq) 网格上搜索满足电压/电流约束的最大扭矩点.
    纯解析, ~1 秒完成 150 转速点.
    """
    # ── 加载磁链表 ──
    data = np.load(psi_dq_npz)
    id_grid = data['id_grid']
    iq_grid = data['iq_grid']
    psi_d_tbl = data['psi_d']
    psi_q_tbl = data['psi_q']

    psi_d_fn = RegularGridInterpolator((id_grid, iq_grid), psi_d_tbl,
                                        bounds_error=False, fill_value=None)
    psi_q_fn = RegularGridInterpolator((id_grid, iq_grid), psi_q_tbl,
                                        bounds_error=False, fill_value=None)

    Vmax = vdc / np.sqrt(3)
    print(f'\n{"=" * 70}')
    print(f'  Sub-flow E: 外特性 T-n 曲线 (基于 ψd/ψq 磁链表)')
    print(f'{"=" * 70}')
    print(f'  磁链表:  {psi_dq_npz}')
    print(f'  Vdc={vdc}V  Vmax={Vmax:.1f}V  Imax={imax}A  P={pole_pairs}')
    print(f'  转速: 1 ~ {speed_max} rpm, {speed_points_n} 点')
    print(f'{"=" * 70}')

    # ── 精细 (Id,Iq) 网格 ──
    N_fine = 200
    id_f = np.linspace(id_grid[0], id_grid[-1], N_fine)
    iq_f = np.linspace(iq_grid[0], iq_grid[-1], N_fine)
    ID, IQ = np.meshgrid(id_f, iq_f, indexing='ij')
    pts = np.column_stack([ID.ravel(), IQ.ravel()])
    psi_d = psi_d_fn(pts).reshape(N_fine, N_fine)
    psi_q = psi_q_fn(pts).reshape(N_fine, N_fine)
    T_grid = 1.5 * pole_pairs * (psi_d * IQ - psi_q * ID)
    Is_grid = np.sqrt(ID**2 + IQ**2)

    # 电流约束 mask
    current_ok = Is_grid <= imax

    # ── 转速扫描 ──
    speeds = np.geomspace(1, speed_max, speed_points_n)
    results = []
    base_speed = None

    for n in speeds:
        omega_e = n * np.pi / 30 * pole_pairs
        Vs_grid = omega_e * np.sqrt(psi_d**2 + psi_q**2)
        feasible = current_ok & (Vs_grid <= Vmax)

        if not np.any(feasible):
            results.append({
                'speed': n, 'T_max': 0.0, 'P_max_kW': 0.0,
                'Id_opt': 0.0, 'Iq_opt': 0.0, 'Is_opt': 0.0,
                'beta_opt': 0.0, 'Vs': 0.0, 'M': 0.0, 'region': 'none',
            })
            continue

        # ── Step 1: grid search for initial guess ──
        idx = np.argmax(np.where(feasible, T_grid, -np.inf))
        i0, j0 = np.unravel_index(idx, (N_fine, N_fine))
        Id0, Iq0 = float(ID[i0, j0]), float(IQ[i0, j0])

        # ── Step 2: SLSQP refinement along voltage boundary ──
        # Only refine if we are near the voltage limit (FW region)
        Vs0 = float(omega_e * np.sqrt(psi_d[i0,j0]**2 + psi_q[i0,j0]**2))
        M0 = Vs0 / Vmax

        def _torque_objective(x):
            """Negative torque (for minimization)."""
            Id_, Iq_ = x[0], x[1]
            pd = float(psi_d_fn([(Id_, Iq_)])[0])
            pq = float(psi_q_fn([(Id_, Iq_)])[0])
            T = 1.5 * pole_pairs * (pd * Iq_ - pq * Id_)
            return -T

        def _voltage_constraint(x):
            """Vs ≤ Vmax → Vmax - Vs ≥ 0."""
            Id_, Iq_ = x[0], x[1]
            pd = float(psi_d_fn([(Id_, Iq_)])[0])
            pq = float(psi_q_fn([(Id_, Iq_)])[0])
            Vs = omega_e * np.sqrt(pd**2 + pq**2)
            return Vmax - Vs

        def _current_constraint(x):
            """Is ≤ Imax → Imax - Is ≥ 0."""
            return imax - np.sqrt(x[0]**2 + x[1]**2)

        constraints = [
            {'type': 'ineq', 'fun': _voltage_constraint},
            {'type': 'ineq', 'fun': _current_constraint},
        ]
        bounds = [(id_grid[0], id_grid[-1]), (iq_grid[0], iq_grid[-1])]

        try:
            res = minimize(
                _torque_objective, x0=[Id0, Iq0],
                method='SLSQP', bounds=bounds, constraints=constraints,
                options={'ftol': 1e-8, 'maxiter': 100, 'disp': False}
            )
            if res.success:
                Id_opt, Iq_opt = float(res.x[0]), float(res.x[1])
            else:
                Id_opt, Iq_opt = Id0, Iq0
        except Exception:
            Id_opt, Iq_opt = Id0, Iq0

        # ── Evaluate final operating point ──
        pd_opt = float(psi_d_fn([(Id_opt, Iq_opt)])[0])
        pq_opt = float(psi_q_fn([(Id_opt, Iq_opt)])[0])
        Is_opt = np.sqrt(Id_opt**2 + Iq_opt**2)
        Vs_opt = float(omega_e * np.sqrt(pd_opt**2 + pq_opt**2))
        T_max = 1.5 * pole_pairs * (pd_opt * Iq_opt - pq_opt * Id_opt)
        M = Vs_opt / Vmax
        beta = np.degrees(np.arctan2(Id_opt, Iq_opt))
        P_max = T_max * n * 2 * np.pi / 60 / 1000

        # Region classification
        if abs(Is_opt - imax) < imax * 0.01:
            region = 'MTPA'
        else:
            region = 'FW'
            if base_speed is None:
                base_speed = n

        results.append({
            'speed': n, 'T_max': T_max, 'P_max_kW': P_max,
            'Id_opt': Id_opt, 'Iq_opt': Iq_opt, 'Is_opt': Is_opt,
            'beta_opt': beta, 'Vs': Vs_opt, 'M': M, 'region': region,
        })

    # ── 输出表格 ──
    print(f'\n{"n_rpm":>8s} {"T_Nm":>8s} {"kW":>7s} {"Id_A":>7s} {"Iq_A":>7s} '
          f'{"Is_A":>7s} {"Beta":>6s} {"M":>6s} {"区域":>5s}')
    print('-' * 76)
    for r in results:
        if r['T_max'] > 0.1:
            print(f'{r["speed"]:8.0f} {r["T_max"]:8.1f} {r["P_max_kW"]:7.2f} '
                  f'{r["Id_opt"]:7.1f} {r["Iq_opt"]:7.1f} {r["Is_opt"]:7.1f} '
                  f'{r["beta_opt"]:6.1f} {r["M"]:6.3f} {r["region"]:>5s}')

    # ── 统计 ──
    valid = [r for r in results if r['T_max'] > 0.1]
    peak_T = max(valid, key=lambda r: r['T_max'])
    peak_P = max(valid, key=lambda r: r['P_max_kW'])
    fw_pts = [r for r in valid if r['region'] == 'FW']
    fw_start = fw_pts[0]['speed'] if fw_pts else None
    last = valid[-1]
    taper = last['T_max'] / peak_T['T_max'] * 100 if peak_T['T_max'] > 0 else 0

    print(f'\n{"=" * 70}')
    print(f'  峰值扭矩:  {peak_T["T_max"]:.1f} Nm @ {peak_T["speed"]:.0f} rpm')
    print(f'  峰值功率:  {peak_P["P_max_kW"]:.2f} kW @ {peak_P["speed"]:.0f} rpm')
    if fw_start:
        print(f'  弱磁拐点:  {fw_start:.0f} rpm')
    print(f'  万转扭矩:  {last["T_max"]:.1f} Nm (下垂率 {taper:.1f}%)')
    print(f'  MTPA 区间:  {peak_T["speed"]:.0f} ~ {fw_start:.0f} rpm' if fw_start else '')
    print(f'{"=" * 70}')

    # ── 保存 NPZ ──
    out = Path(out_dir) if out_dir else Path(psi_dq_npz).parent
    out.mkdir(parents=True, exist_ok=True)
    npz_path = out / 'outer_characteristic.npz'
    np.savez_compressed(
        npz_path,
        speed=np.array([r['speed'] for r in results]),
        T_max=np.array([r['T_max'] for r in results]),
        P_max_kW=np.array([r['P_max_kW'] for r in results]),
        Is_opt=np.array([r['Is_opt'] for r in results]),
        Id_opt=np.array([r['Id_opt'] for r in results]),
        Iq_opt=np.array([r['Iq_opt'] for r in results]),
        beta_opt=np.array([r['beta_opt'] for r in results]),
        Vs=np.array([r['Vs'] for r in results]),
        M=np.array([r['M'] for r in results]),
        region=np.array([r['region'] for r in results], dtype=object),
        vdc=vdc, imax=imax, pole_pairs=pole_pairs, rs=rs,
    )
    print(f'\n  已保存: {npz_path}')

    # ── 绘图 ──
    if plot:
        png_path = out / f'outer_char_{vdc:.0f}V_{imax:.0f}A.png'
        _plot_outer_char(results, vdc, imax, str(png_path))

    return results


def _plot_outer_char(results, vdc, imax, out_path):
    """生成 4-panel 外特性图."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'DejaVu Sans']
    plt.rcParams['axes.unicode_minus'] = False

    speeds = np.array([r['speed'] for r in results if r['T_max'] > 0.1])
    T = np.array([r['T_max'] for r in results if r['T_max'] > 0.1])
    P = np.array([r['P_max_kW'] for r in results if r['T_max'] > 0.1])
    Id = np.array([r['Id_opt'] for r in results if r['T_max'] > 0.1])
    Iq = np.array([r['Iq_opt'] for r in results if r['T_max'] > 0.1])
    Is = np.array([r['Is_opt'] for r in results if r['T_max'] > 0.1])
    Vs = np.array([r['Vs'] for r in results if r['T_max'] > 0.1])
    regions = [r['region'] for r in results if r['T_max'] > 0.1]

    # FW corner
    fw_idx = next((i for i, r in enumerate(regions) if r == 'FW'), len(speeds)-1)
    mtpa_mask = np.array([r == 'MTPA' for r in regions])
    fw_mask = np.array([r == 'FW' for r in regions])

    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle(f'PMSM 外特性  (Vdc={vdc}V, Imax={imax}A)',
                 fontsize=15, fontweight='bold')

    # Panel 1: T-n + P-n
    ax = axes[0, 0]
    ax_twin = ax.twinx()
    ax.plot(speeds, T, 'b-', lw=2, label='Torque')
    ax_twin.plot(speeds, P, 'r-', lw=2, label='Power')
    ax.axvline(speeds[fw_idx], color='gray', ls='--', alpha=0.7,
               label=f'FW corner ≈ {speeds[fw_idx]:.0f} rpm')
    ax.set(xlabel='Speed (rpm)', ylabel='Torque (Nm)', title='T-n / P-n')
    ax.legend(loc='center left', fontsize=8)
    ax_twin.legend(loc='upper right', fontsize=8)
    ax_twin.set_ylabel('Power (kW)')
    ax.grid(True, alpha=0.3)

    # Panel 2: Is + Vs
    ax = axes[0, 1]
    ax_twin = ax.twinx()
    ax.plot(speeds, Is, 'g-', lw=2, label='Is (A)')
    ax_twin.plot(speeds, Vs, 'm-', lw=2, label='Vs (V)')
    ax.axhline(imax, color='g', ls='--', alpha=0.5, label=f'Imax={imax}A')
    ax.axvline(speeds[fw_idx], color='gray', ls='--', alpha=0.7)
    ax.set(xlabel='Speed (rpm)', ylabel='Is (A)', title='Is & Vs')
    ax.legend(loc='center left', fontsize=8)
    ax_twin.legend(loc='upper right', fontsize=8)
    ax_twin.set_ylabel('Vs (V)')
    ax.grid(True, alpha=0.3)

    # Panel 3: Id/Iq trajectory
    ax = axes[1, 0]
    ax.plot(speeds[mtpa_mask], Id[mtpa_mask], 'r-', lw=2, label='Id (MTPA)')
    ax.plot(speeds[fw_mask], Id[fw_mask], 'r--', lw=2, label='Id (FW)')
    ax.plot(speeds[mtpa_mask], Iq[mtpa_mask], 'b-', lw=2, label='Iq (MTPA)')
    ax.plot(speeds[fw_mask], Iq[fw_mask], 'b--', lw=2, label='Iq (FW)')
    ax.axvline(speeds[fw_idx], color='gray', ls='--', alpha=0.7)
    ax.set(xlabel='Speed (rpm)', ylabel='Current (A)', title='Id / Iq trajectory')
    ax.legend(loc='best', fontsize=8)
    ax.grid(True, alpha=0.3)

    # Panel 4: dq-plane
    ax = axes[1, 1]
    theta_c = np.linspace(0, 2*np.pi, 200)
    ax.plot(imax * np.sin(theta_c), imax * np.cos(theta_c), 'gray', ls='--', alpha=0.5)
    ax.text(imax*0.7, -imax*0.15, f'{imax}A', fontsize=8, color='gray')
    sc = ax.scatter(Id, Iq, c=speeds, cmap='coolwarm', s=8, alpha=0.8)
    ax.plot(Id[0], Iq[0], 'go', ms=8, label='Start (MTPA)')
    ax.plot(Id[-1], Iq[-1], 'ro', ms=8, label='End (FW)')
    ax.set(xlabel='Id (A)', ylabel='Iq (A)', title='dq-plane trajectory')
    ax.axhline(0, color='k', alpha=0.2)
    ax.axvline(0, color='k', alpha=0.2)
    ax.legend(loc='best', fontsize=8)
    ax.grid(True, alpha=0.3)
    cb = fig.colorbar(sc, ax=ax, label='Speed (rpm)', shrink=0.8)

    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'  图表已保存: {out_path}')


# ═══════════════════════════════════════════════════════════════
#  Ld/Lq 参数模型 (向后兼容)
# ═══════════════════════════════════════════════════════════════

def _run_ldlq_model(vdc, imax, speed_max, speed_points_n, pole_pairs, rs,
                    phi, ld0, lq, k1, k2, model, ld_lq_data):
    """旧的 Ld/Lq 参数模型外特性计算 (向后兼容)."""
    P = pole_pairs if pole_pairs is not None else _DEFAULT_P
    Rs = rs if rs is not None else _DEFAULT_Rs
    Lq = lq if lq is not None else _DEFAULT_Lq
    Ld0 = ld0 if ld0 is not None else _DEFAULT_Ld0
    Phi = phi if phi is not None else _DEFAULT_Phi
    K1 = k1 if k1 is not None else _DEFAULT_k1
    K2 = k2 if k2 is not None else _DEFAULT_k2

    if ld_lq_data is not None:
        ld_vals = np.array([d['Ld'] for d in ld_lq_data])
        lq_vals = np.array([d['Lq'] for d in ld_lq_data])
        Ld0 = np.mean(ld_vals[ld_vals > 0]) if np.any(ld_vals > 0) else Ld0
        Lq = np.mean(lq_vals[lq_vals > 0]) if np.any(lq_vals > 0) else Lq

    Vmax = vdc / np.sqrt(3)

    if model == 'quadratic':
        def Ld_at(Id_val):
            return Ld0 + K1 * abs(Id_val) + K2 * Id_val**2
    elif model == 'linear':
        def Ld_at(Id_val):
            return Ld0 + K1 * abs(Id_val)
    else:
        def Ld_at(Id_val):
            return Ld0

    dL_sign = 'Ld > Lq (反凸极)' if (Ld0 - Lq) > 0 else 'Ld < Lq (传统 IPM)'
    print(f'\n{"=" * 70}')
    print(f'  Sub-flow E: 外特性 T-n 曲线 (Ld/Lq 参数模型)')
    print(f'{"=" * 70}')
    print(f'  Vdc={vdc:.0f}V  Vmax={Vmax:.1f}V  Imax={imax}A  P={P}')
    print(f'  Phi={Phi*1000:.4f}mWb  Lq={Lq*1000:.4f}mH  Ld0={Ld0*1000:.4f}mH ({dL_sign})')
    print(f'  模型: {model}')
    print(f'{"=" * 70}')

    speeds = np.geomspace(max(1, speed_max / speed_points_n), speed_max, speed_points_n)
    N_BETA, N_IS = 360, 200
    beta_grid = np.linspace(-np.pi/2, np.pi/2, N_BETA)
    Is_grid = np.linspace(0, imax, N_IS)
    Beta_mesh, Is_mesh = np.meshgrid(beta_grid, Is_grid, indexing='ij')
    Id_all = Is_mesh * np.sin(Beta_mesh)
    Iq_all = Is_mesh * np.cos(Beta_mesh)

    if model == 'quadratic':
        Ld_all = Ld0 + K1 * np.abs(Id_all) + K2 * Id_all**2
    elif model == 'linear':
        Ld_all = Ld0 + K1 * np.abs(Id_all)
    else:
        Ld_all = np.full_like(Id_all, Ld0)
    dL_all = Ld_all - Lq

    results = []
    for n in speeds:
        omega_e = 2 * np.pi * n / 60 * P
        Vd = -omega_e * Lq * Iq_all + Rs * Id_all
        Vq = omega_e * (Ld_all * Id_all + Phi) + Rs * Iq_all
        Vs = np.sqrt(Vd**2 + Vq**2)
        feasible = Vs <= Vmax

        if not np.any(feasible):
            results.append({
                'speed': n, 'torque': 0.0, 'power': 0.0,
                'Id': 0.0, 'Iq': 0.0, 'beta_deg': 0.0,
                'modulation': 0.0, 'Is': 0.0, 'region': 'FW',
            })
            continue

        T_all = np.full_like(Id_all, -1e12)
        T_feasible = 1.5 * P * (Phi * Iq_all + dL_all * Id_all * Iq_all)
        T_all[feasible] = T_feasible[feasible]
        idx_flat = np.argmax(T_all)
        idx_beta, idx_Is = np.unravel_index(idx_flat, Id_all.shape)
        best_T = T_all[idx_beta, idx_Is]
        best_Id = Id_all[idx_beta, idx_Is]
        best_Iq = Iq_all[idx_beta, idx_Is]
        best_Is = Is_grid[idx_Is]
        best_beta = beta_grid[idx_beta]
        region = 'MTPA' if abs(best_Is - imax) < 1e-3 else 'FW'
        P_out = best_T * n * 2 * np.pi / 60 / 1000
        Ld_op = Ld_at(best_Id)
        Vd_op = -omega_e * Lq * best_Iq + Rs * best_Id
        Vq_op = omega_e * (Ld_op * best_Id + Phi) + Rs * best_Iq
        Vs_op = np.sqrt(Vd_op**2 + Vq_op**2)
        M = Vs_op / Vmax

        results.append({
            'speed': n, 'torque': best_T, 'power': P_out,
            'Id': best_Id, 'Iq': best_Iq, 'beta_deg': np.degrees(best_beta),
            'modulation': M, 'Is': best_Is, 'region': region,
        })

    # 输出 (与 _run_psidq 格式兼容)
    print(f'\n{"n_rpm":>8s} {"T_Nm":>8s} {"kW":>7s} {"Id_A":>7s} {"Iq_A":>7s} '
          f'{"Is_A":>7s} {"Beta":>6s} {"M":>6s} {"区域":>5s}')
    print('-' * 76)
    for r in results:
        print(f'{r["speed"]:8.0f} {r["torque"]:8.1f} {r["power"]:7.2f} '
              f'{r["Id"]:7.2f} {r["Iq"]:7.2f} {r["Is"]:7.1f} '
              f'{r["beta_deg"]:6.1f} {r["modulation"]:7.3f} {r["region"]:>5s}')

    return results


# ═══════════════════════════════════════════════════════════════
#  主入口
# ═══════════════════════════════════════════════════════════════

def run(vdc=300, imax=250, speed_max=12000, speed_points_n=150,
        pole_pairs=None, rs=None, project_path=None,
        # Ld/Lq 模型参数 (向后兼容)
        phi=None, ld0=None, lq=None, k1=None, k2=None, model='quadratic',
        ld_lq_data=None,
        # 磁链表模式参数
        psi_dq_npz=None):
    """计算外特性 T-n 曲线.

    两种模式 (自动选择):
      1. ψd/ψq 磁链表 (推荐): project_path 或 psi_dq_npz 指定时自动使用
      2. Ld/Lq 参数模型 (向后兼容): 无磁链表时使用

    Parameters
    ----------
    vdc : float — 直流母线电压 (V)
    imax : float — 最大相电流 (A)
    speed_max : float — 最高转速 (rpm)
    speed_points_n : int — 转速点数
    pole_pairs : int — 极对数, 默认 4
    rs : float — 相电阻 (Ω)
    project_path : str — 项目目录 (自动查找 psi_dq_table.npz)
    psi_dq_npz : str — 磁链表路径 (优先于 project_path)
    (以下为 Ld/Lq 模型参数)
    phi, ld0, lq, k1, k2, model, ld_lq_data

    Returns
    -------
    list of dict — 每个转速点的工作点信息
    """
    # ── 尝试 ψd/ψq 磁链表模式 ──
    if psi_dq_npz is None and project_path is not None:
        candidate = Path(project_path) / 'psi_dq_table.npz'
        if candidate.exists():
            psi_dq_npz = str(candidate)

    if psi_dq_npz is not None:
        P = pole_pairs if pole_pairs is not None else _DEFAULT_P
        Rs = rs if rs is not None else 0.05
        return _run_psidq(
            psi_dq_npz, vdc=vdc, imax=imax, speed_max=speed_max,
            speed_points_n=speed_points_n, pole_pairs=P, rs=Rs,
            out_dir=project_path)

    # ── 回退到 Ld/Lq 参数模型 ──
    return _run_ldlq_model(
        vdc=vdc, imax=imax, speed_max=speed_max,
        speed_points_n=speed_points_n, pole_pairs=pole_pairs, rs=rs,
        phi=phi, ld0=ld0, lq=lq, k1=k1, k2=k2, model=model,
        ld_lq_data=ld_lq_data)


if __name__ == '__main__':
    run()
