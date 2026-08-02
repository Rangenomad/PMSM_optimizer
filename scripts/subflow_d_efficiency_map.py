"""Sub-flow D: 全域工作特性效率 MAP

基于 Sub-flow A 输出的 ψd/ψq 磁链表，通过 MTPA + 弱磁控制确定每个
(转速, 扭矩) 工作点的最优 (Id, Iq)，解析计算效率、损耗、功率因数等。

方法:
  - 默认: 纯解析，基于 ψd(Id,Iq) / ψq(Id,Iq) RegularGridInterpolator
  - 可选: FEA 验证 — 在解析最优 (Id,Iq) 点运行 Transient FEA（待实现）

依赖:
  - Sub-flow A → psi_dq_table.npz（ψd/ψq 磁链表）
  - scipy ≥ 1.7 (RegularGridInterpolator)

用法:
    # 从项目目录自动查找 psi_dq_table.npz
    python -c "from scripts.subflow_d_efficiency_map import run; run(
        project_path='pmsm_projects/2026-08-02_mtpa_cal_dense_v2',
        vdc=500, imax=300)"

    # 直接指定 psi_dq_table.npz 路径
    python -c "from scripts.subflow_d_efficiency_map import run; run(
        psi_dq_npz='pmsm_projects/xxx/psi_dq_table.npz',
        vdc=500, imax=300, n_speed=30, n_torque=30)"
"""

import sys
sys.stdout.reconfigure(encoding='utf-8')

from pathlib import Path
import numpy as np
from scipy.interpolate import RegularGridInterpolator
import warnings
warnings.filterwarnings('ignore')

ROOT = Path(__file__).resolve().parent.parent


# ═══════════════════════════════════════════════════════════════
# 铁耗模型参数
# ═══════════════════════════════════════════════════════════════
# P_fe = k_fe * (f_e / f0)^alpha * (psi / psi_pm)^beta
#   f_e   = speed * pole_pairs / 60  (Hz)
#   psi   = sqrt(psi_d^2 + psi_q^2)  — 当前工作点磁链
#   psi_pm = PM 磁链 (Id=0, Iq=0)    — 基准磁链
# 典型值 (50-70kW PMSM): ~2% 额定功率 @ 3000rpm 为铁耗
DEFAULT_k_fe = 300.0      # W — 基准铁耗 @ (f0, psi_pm)
DEFAULT_f0 = 200.0        # Hz — 基准电频率 (3000rpm × 4poles / 60)
DEFAULT_alpha = 1.6       # — 频率指数 (1.5-2.0)
DEFAULT_beta = 1.8        # — 磁链指数 (1.5-2.0)


# ═══════════════════════════════════════════════════════════════
# 核心计算
# ═══════════════════════════════════════════════════════════════

class EfficiencyMapComputer:
    """基于 ψd/ψq 磁链表 + MTPA/FW 控制的效率 MAP 计算器."""

    def __init__(self, psi_dq_npz):
        """加载磁链表并构建插值器。

        Parameters
        ----------
        psi_dq_npz : str or Path
            Sub-flow A 输出的 psi_dq_table.npz
        """
        data = np.load(psi_dq_npz)
        self.id_grid = data['id_grid']       # (N_id,)
        self.iq_grid = data['iq_grid']       # (N_iq,)
        self.psi_d_tbl = data['psi_d']       # (N_id, N_iq)
        self.psi_q_tbl = data['psi_q']       # (N_id, N_iq)

        self.psi_d_fn = RegularGridInterpolator(
            (self.id_grid, self.iq_grid), self.psi_d_tbl,
            bounds_error=False, fill_value=None)
        self.psi_q_fn = RegularGridInterpolator(
            (self.id_grid, self.iq_grid), self.psi_q_tbl,
            bounds_error=False, fill_value=None)

        # PM 磁链 (Id=0, Iq=0)
        self.psi_pm = float(self.psi_d_fn([[0.0, 0.0]])[0])

        # 缓存: 精细网格预计算数据
        self._cache = None

    def _ensure_cache(self, pole_pairs, imax, n_fine=250):
        """预计算精细 (Id,Iq) 网格上的所有标量场."""
        if self._cache is not None and self._cache['n_fine'] == n_fine:
            return

        id_f = np.linspace(self.id_grid[0], self.id_grid[-1], n_fine)
        iq_f = np.linspace(self.iq_grid[0], self.iq_grid[-1], n_fine)
        ID, IQ = np.meshgrid(id_f, iq_f, indexing='ij')
        pts = np.column_stack([ID.ravel(), IQ.ravel()])

        psi_d = self.psi_d_fn(pts).reshape(n_fine, n_fine)
        psi_q = self.psi_q_fn(pts).reshape(n_fine, n_fine)
        T = 1.5 * pole_pairs * (psi_d * IQ - psi_q * ID)
        Is = np.sqrt(ID**2 + IQ**2)
        Vs_per_rpm = (np.pi / 30 * pole_pairs) * np.sqrt(psi_d**2 + psi_q**2)
        psi_mag = np.sqrt(psi_d**2 + psi_q**2)

        self._cache = {
            'n_fine': n_fine,
            'id_f': id_f, 'iq_f': iq_f,
            'ID': ID, 'IQ': IQ,
            'psi_d': psi_d, 'psi_q': psi_q,
            'T': T, 'Is': Is,
            'Vs_per_rpm': Vs_per_rpm,
            'psi_mag': psi_mag,
        }

    def _build_mtpa_lut(self, imax, pole_pairs, n_levels=500):
        """构建 MTPA 查找表: T → (Id,Iq,Is,Vs_per_rpm,psi_mag)."""
        self._ensure_cache(pole_pairs, imax)
        c = self._cache
        T_f, Is_f = c['T'], c['Is']

        T_max = np.max(T_f[Is_f <= imax])
        T_levels = np.linspace(0, T_max, n_levels)

        lut = np.zeros(n_levels, dtype=[
            ('T', 'f8'), ('Id', 'f8'), ('Iq', 'f8'),
            ('Is', 'f8'), ('Vs_per_rpm', 'f8'), ('psi_mag', 'f8'),
            ('found', 'bool')
        ])
        lut['T'] = T_levels

        for k, T_tgt in enumerate(T_levels):
            if T_tgt < 0.5:
                lut['found'][k] = True
                continue

            tol = max(0.5, T_tgt * 0.003)
            mask = (np.abs(T_f - T_tgt) < tol) & (Is_f <= imax)

            if not np.any(mask):
                # 最近可行点
                feasible_flat = np.flatnonzero(Is_f <= imax)
                if len(feasible_flat) == 0:
                    continue
                idx = feasible_flat[np.argmin(np.abs(T_f.ravel()[feasible_flat] - T_tgt))]
                i_idx, j_idx = np.unravel_index(idx, (c['n_fine'], c['n_fine']))
            else:
                idx = np.argmin(Is_f[mask])
                flat_idx = np.flatnonzero(mask)[idx]
                i_idx, j_idx = np.unravel_index(flat_idx, (c['n_fine'], c['n_fine']))

            lut['Id'][k] = c['ID'][i_idx, j_idx]
            lut['Iq'][k] = c['IQ'][i_idx, j_idx]
            lut['Is'][k] = c['Is'][i_idx, j_idx]
            lut['Vs_per_rpm'][k] = c['Vs_per_rpm'][i_idx, j_idx]
            lut['psi_mag'][k] = c['psi_mag'][i_idx, j_idx]
            lut['found'][k] = True

        self._mtpa_lut = lut
        self._T_max_global = T_max
        return lut, T_max

    def _solve_fw(self, T_target, omega_e, Vmax, imax, pole_pairs,
                  id_start, iq_start, n_scan=400):
        """弱磁区搜索: 沿恒扭矩曲线找 Vs=Vmax 的工作点.

        从 MTPA 点的 Id 开始, 向更负方向扫描, 对每个 Id 找 Iq 使 T=T_target,
        然后检查 Vs ≤ Vmax 且 Is ≤ imax.

        Returns
        -------
        (Id, Iq, Is, Vs, psi_mag, success) : tuple
        """
        id_scan = np.linspace(id_start, self.id_grid[0], n_scan)
        iq_try_arr = np.linspace(0, self.iq_grid[-1], n_scan)

        for id_try in id_scan:
            # 计算 T(Id=id_try, Iq) 沿 Iq 的分布
            psi_d_iq = self.psi_d_fn(
                np.column_stack([np.full_like(iq_try_arr, id_try), iq_try_arr]))
            psi_q_iq = self.psi_q_fn(
                np.column_stack([np.full_like(iq_try_arr, id_try), iq_try_arr]))
            T_iq = 1.5 * pole_pairs * (psi_d_iq * iq_try_arr - psi_q_iq * id_try)

            if np.max(T_iq) < T_target:
                continue

            iq_for_t = np.interp(T_target, T_iq, iq_try_arr)
            Is_try = np.sqrt(id_try**2 + iq_for_t**2)
            if Is_try > imax:
                continue

            psi_d_val = float(self.psi_d_fn([[id_try, iq_for_t]])[0])
            psi_q_val = float(self.psi_q_fn([[id_try, iq_for_t]])[0])
            Vs_try = omega_e * np.sqrt(psi_d_val**2 + psi_q_val**2)

            if Vs_try <= Vmax * 1.001:
                psi_mag = np.sqrt(psi_d_val**2 + psi_q_val**2)
                return (float(id_try), float(iq_for_t), float(Is_try),
                        float(Vs_try), float(psi_mag), True)

        return (float(id_start), float(iq_start), 0.0, 0.0, 0.0, False)

    def compute(self, vdc=500, imax=300, speed_min=500, speed_max=8000,
                n_speed=20, n_torque=20, Rs=0.05, pole_pairs=4,
                k_fe=DEFAULT_k_fe, f0=DEFAULT_f0,
                alpha=DEFAULT_alpha, beta=DEFAULT_beta):
        """计算效率 MAP.

        Returns
        -------
        dict with keys:
            speed_grid, T_grid, eta_map, pfe_map, pcu_map, ploss_map,
            pout_map, is_map, id_map, iq_map, vs_map, m_map, pf_map,
            psi_map, feasible_map, region_map, points, T_max_global, ...
        """
        Vmax = vdc / np.sqrt(3)
        print(f'[D] Vdc={vdc}V, Imax={imax}A, Vmax={Vmax:.1f}V')
        print(f'[D] 转速: {speed_min}~{speed_max} rpm, {n_speed} 点')
        print(f'[D] 扭矩: {n_torque} 点, Rs={Rs}Ω, P={pole_pairs}')
        print(f'[D] 铁耗: k_fe={k_fe}W @ {f0}Hz, α={alpha}, β={beta}')

        # 构建 MTPA LUT
        lut, T_max = self._build_mtpa_lut(imax, pole_pairs)
        T_levels = lut['T']
        print(f'[D] 全局最大扭矩: {T_max:.1f} Nm')

        # 网格
        speed_grid = np.linspace(speed_min, speed_max, n_speed)
        T_grid = np.linspace(0, T_max, n_torque)

        # 预分配
        shape = (n_speed, n_torque)
        eta_map = np.full(shape, np.nan)
        pfe_map = np.full(shape, np.nan)
        pcu_map = np.full(shape, np.nan)
        ploss_map = np.full(shape, np.nan)
        pout_map = np.full(shape, np.nan)
        is_map = np.full(shape, np.nan)
        id_map = np.full(shape, np.nan)
        iq_map = np.full(shape, np.nan)
        vs_map = np.full(shape, np.nan)
        m_map = np.full(shape, np.nan)
        pf_map = np.full(shape, np.nan)
        psi_map = np.full(shape, np.nan)
        feasible_map = np.full(shape, False)
        region_map = np.full(shape, '', dtype=object)

        n_mtpa = n_fw = n_infeas = 0
        points = []

        for i_s, speed in enumerate(speed_grid):
            omega_e = speed * np.pi / 30 * pole_pairs

            for j_t, T_target in enumerate(T_grid):
                pt = self._eval_one_point(
                    T_target, speed, omega_e, Vmax, imax, pole_pairs,
                    Rs, k_fe, f0, alpha, beta, lut, T_levels)
                points.append(pt)

                feasible_map[i_s, j_t] = pt['feasible']
                if pt['feasible']:
                    eta_map[i_s, j_t] = pt['eta']
                    pfe_map[i_s, j_t] = pt['P_fe']
                    pcu_map[i_s, j_t] = pt['P_cu']
                    ploss_map[i_s, j_t] = pt['P_loss']
                    pout_map[i_s, j_t] = pt['P_out']
                    is_map[i_s, j_t] = pt['Is']
                    id_map[i_s, j_t] = pt['Id']
                    iq_map[i_s, j_t] = pt['Iq']
                    vs_map[i_s, j_t] = pt['Vs']
                    m_map[i_s, j_t] = pt['M']
                    pf_map[i_s, j_t] = pt['PF']
                    psi_map[i_s, j_t] = pt['psi_mag']
                    region_map[i_s, j_t] = pt['region']
                    if pt['region'] == 'MTPA':
                        n_mtpa += 1
                    elif pt['region'] == 'FW':
                        n_fw += 1
                else:
                    n_infeas += 1

        total = n_speed * n_torque
        print(f'[D] {n_speed}×{n_torque}={total} 点: '
              f'MTPA={n_mtpa}, FW={n_fw}, 不可行={n_infeas}')

        return {
            'speed_grid': speed_grid,
            'T_grid': T_grid,
            'eta_map': eta_map, 'pfe_map': pfe_map, 'pcu_map': pcu_map,
            'ploss_map': ploss_map, 'pout_map': pout_map,
            'is_map': is_map, 'id_map': id_map, 'iq_map': iq_map,
            'vs_map': vs_map, 'm_map': m_map, 'pf_map': pf_map,
            'psi_map': psi_map, 'feasible_map': feasible_map,
            'region_map': region_map,
            'points': points,
            'T_max_global': T_max, 'psi_pm': self.psi_pm,
            'vdc': vdc, 'imax': imax, 'Vmax': Vmax,
            'Rs': Rs, 'pole_pairs': pole_pairs,
            'k_fe': k_fe, 'f0': f0, 'alpha': alpha, 'beta': beta,
        }

    def _eval_one_point(self, T_target, speed, omega_e, Vmax, imax, pole_pairs,
                        Rs, k_fe, f0, alpha, beta, lut, T_levels):
        """计算单个工作点."""
        base = {
            'speed': speed, 'T_target': T_target,
            'Id': 0.0, 'Iq': 0.0, 'Is': 0.0,
            'psi_d': float(self.psi_d_fn([[0.0, 0.0]])[0]),
            'psi_q': 0.0, 'psi_mag': float(self.psi_d_fn([[0.0, 0.0]])[0]),
            'Vd': 0.0, 'Vq': 0.0, 'Vs': 0.0, 'M': 0.0, 'PF': 0.0,
            'P_out': 0.0, 'P_cu': 0.0, 'P_fe': 0.0, 'P_loss': 0.0,
            'eta': 0.0, 'feasible': False, 'region': 'idle',
        }

        if T_target < 0.5:
            # 零扭矩 / 空载点
            psi_pm = self.psi_pm
            base['Vq'] = omega_e * psi_pm
            base['Vs'] = omega_e * psi_pm
            base['M'] = base['Vs'] / Vmax if Vmax > 0 else 0.0
            base['feasible'] = True
            base['region'] = 'idle'
            return base

        # ── 1. 查 MTPA LUT ──
        k = np.argmin(np.abs(T_levels - T_target))
        Id_mtpa = float(lut['Id'][k])
        Iq_mtpa = float(lut['Iq'][k])
        Is_mtpa = float(lut['Is'][k])

        if not lut['found'][k]:
            base['Id'] = Id_mtpa
            base['Iq'] = Iq_mtpa
            base['Is'] = Is_mtpa
            return base  # infeasible

        # ── 2. 计算 ψd, ψq ──
        psi_d_val = float(self.psi_d_fn([[Id_mtpa, Iq_mtpa]])[0])
        psi_q_val = float(self.psi_q_fn([[Id_mtpa, Iq_mtpa]])[0])
        psi_mag = np.sqrt(psi_d_val**2 + psi_q_val**2)

        Vs_mtpa = omega_e * psi_mag

        # ── 3. 判断 MTPA vs FW ──
        if Vs_mtpa <= Vmax:
            Id, Iq, Is = Id_mtpa, Iq_mtpa, Is_mtpa
            Vs = Vs_mtpa
            region = 'MTPA'
        else:
            Id_fw, Iq_fw, Is_fw, Vs_fw, psi_fw, ok = self._solve_fw(
                T_target, omega_e, Vmax, imax, pole_pairs,
                Id_mtpa, Iq_mtpa)
            if ok:
                Id, Iq, Is = Id_fw, Iq_fw, Is_fw
                Vs = Vs_fw
                psi_mag = psi_fw
                psi_d_val = float(self.psi_d_fn([[Id, Iq]])[0])
                psi_q_val = float(self.psi_q_fn([[Id, Iq]])[0])
                region = 'FW'
            else:
                base['Id'] = Id_mtpa
                base['Iq'] = Iq_mtpa
                base['Is'] = Is_mtpa
                base['psi_d'] = psi_d_val
                base['psi_q'] = psi_q_val
                base['psi_mag'] = psi_mag
                base['Vs'] = Vs_mtpa
                base['M'] = Vs_mtpa / Vmax
                base['region'] = 'infeasible'
                return base

        # ── 4. 电气量 ──
        Vd = -omega_e * psi_q_val
        Vq = omega_e * psi_d_val
        M = Vs / Vmax
        phi_v = np.arctan2(Vd, Vq)
        phi_i = np.arctan2(Id, Iq)
        PF = np.cos(phi_v - phi_i)

        # ── 5. 损耗 ──
        P_out = T_target * speed * 2 * np.pi / 60
        P_cu = 3.0 * (Is / np.sqrt(2))**2 * Rs
        f_e = speed * pole_pairs / 60
        P_fe = k_fe * (f_e / f0)**alpha * (psi_mag / self.psi_pm)**beta
        P_loss = P_cu + P_fe
        eta = P_out / (P_out + P_loss) * 100.0 if (P_out + P_loss) > 0 else 0.0

        return {
            'speed': speed, 'T_target': T_target,
            'Id': float(Id), 'Iq': float(Iq), 'Is': float(Is),
            'psi_d': psi_d_val, 'psi_q': psi_q_val, 'psi_mag': psi_mag,
            'Vd': Vd, 'Vq': Vq, 'Vs': Vs, 'M': M, 'PF': PF,
            'P_out': P_out, 'P_cu': P_cu, 'P_fe': P_fe,
            'P_loss': P_loss, 'eta': eta,
            'feasible': True, 'region': region,
        }


# ═══════════════════════════════════════════════════════════════
# 绘图
# ═══════════════════════════════════════════════════════════════

def plot_efficiency_map(results, out_path):
    """生成 4-panel 效率 MAP 图.

    Panel 1: 效率 η (%)         Panel 2: 总损耗 P_loss (kW)
    Panel 3: 相电流 Is (A) + M  Panel 4: 功率因数 PF
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei',
                                        'DejaVu Sans']
    plt.rcParams['axes.unicode_minus'] = False

    speed = results['speed_grid']
    T_pts = results['T_grid']
    eta_map = results['eta_map']
    ploss_map = results['ploss_map']
    is_map = results['is_map']
    m_map = results['m_map']
    pf_map = results['pf_map']
    feasible = results['feasible_map']
    vdc = results['vdc']
    imax = results['imax']
    T_max = results['T_max_global']

    # Mask infeasible
    eta_m = np.where(feasible, eta_map, np.nan)
    ploss_m = np.where(feasible, ploss_map / 1000, np.nan)
    is_m = np.where(feasible, is_map, np.nan)
    m_m = np.where(feasible, m_map, np.nan)
    pf_m = np.where(feasible, pf_map, np.nan)

    S, T = np.meshgrid(speed, T_pts, indexing='ij')

    # 外特性包络
    T_env = np.array([np.max(T_pts[feasible[i, :]]) if np.any(feasible[i, :])
                       else 0.0 for i in range(len(speed))])

    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle(f'PMSM 效率 MAP  (Vdc={vdc}V, Imax={imax}A, '
                 f'T$_{{\\rm max}}$={T_max:.0f} Nm)',
                 fontsize=15, fontweight='bold')

    # ── Panel 1: Efficiency ──
    ax = axes[0, 0]
    lv = np.arange(60, 101, 2)
    cf = ax.contourf(S, T, eta_m, levels=lv, cmap='RdYlGn', extend='min')
    ct = ax.contour(S, T, eta_m, levels=lv[::3], colors='black',
                    linewidths=0.3)
    ax.clabel(ct, inline=True, fontsize=7, fmt='%.0f%%')
    ax.plot(speed, T_env, 'k--', lw=1.5, label=u'外特性包络')
    ax.set(xlabel='转速 (rpm)', ylabel='扭矩 (Nm)', title=u'效率 η (%)',
           xlim=(speed[0], speed[-1]), ylim=(0, T_max * 1.05))
    ax.legend(loc='lower left', fontsize=8)
    ax.grid(True, alpha=0.3)
    fig.colorbar(cf, ax=ax, label=u'效率 η (%)', shrink=0.8)

    # ── Panel 2: Total Loss ──
    ax = axes[0, 1]
    lv2 = np.linspace(0, np.nanmax(ploss_m) * 1.1, 15)
    cf = ax.contourf(S, T, ploss_m, levels=lv2, cmap='YlOrRd', extend='max')
    ct = ax.contour(S, T, ploss_m, levels=lv2[::3], colors='black',
                    linewidths=0.3)
    ax.clabel(ct, inline=True, fontsize=7, fmt='%.1f kW')
    ax.plot(speed, T_env, 'k--', lw=1.5, label=u'外特性包络')
    ax.set(xlabel='转速 (rpm)', ylabel='扭矩 (Nm)', title=u'总损耗 P$_{loss}$ (kW)',
           xlim=(speed[0], speed[-1]), ylim=(0, T_max * 1.05))
    ax.legend(loc='upper left', fontsize=8)
    ax.grid(True, alpha=0.3)
    fig.colorbar(cf, ax=ax, label='总损耗 (kW)', shrink=0.8)

    # ── Panel 3: Is + M ──
    ax = axes[1, 0]
    cf = ax.contourf(S, T, is_m, levels=15, cmap='Blues', extend='max')
    ct_i = ax.contour(S, T, is_m, levels=np.arange(0, imax + 1, 30),
                      colors='navy', linewidths=0.8)
    ax.clabel(ct_i, inline=True, fontsize=7, fmt='%.0f A')
    ct_m = ax.contour(S, T, m_m, levels=[0.7, 0.8, 0.9, 1.0],
                      colors='red', linewidths=0.8, linestyles='--')
    ax.clabel(ct_m, inline=True, fontsize=7, fmt='M=%.1f')
    ax.plot(speed, T_env, 'k--', lw=1.5, label=u'外特性包络')
    ax.set(xlabel='转速 (rpm)', ylabel='扭矩 (Nm)',
           title=u'相电流 Is (A) + 调制比 M (红线)',
           xlim=(speed[0], speed[-1]), ylim=(0, T_max * 1.05))
    ax.legend(loc='lower left', fontsize=8)
    ax.grid(True, alpha=0.3)
    fig.colorbar(cf, ax=ax, label='相电流 Is (A)', shrink=0.8)

    # ── Panel 4: Power Factor ──
    ax = axes[1, 1]
    lv4 = np.arange(0.0, 1.05, 0.05)
    cf = ax.contourf(S, T, pf_m, levels=lv4, cmap='PuOr', extend='min')
    ct = ax.contour(S, T, pf_m, levels=lv4[::4], colors='black',
                    linewidths=0.3)
    ax.clabel(ct, inline=True, fontsize=7, fmt='%.2f')
    ax.plot(speed, T_env, 'k--', lw=1.5, label=u'外特性包络')
    ax.set(xlabel='转速 (rpm)', ylabel='扭矩 (Nm)', title=u'功率因数 PF',
           xlim=(speed[0], speed[-1]), ylim=(0, T_max * 1.05))
    ax.legend(loc='lower left', fontsize=8)
    ax.grid(True, alpha=0.3)
    fig.colorbar(cf, ax=ax, label='PF', shrink=0.8)

    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'[D] 图表已保存: {out_path}')


# ═══════════════════════════════════════════════════════════════
# 报告
# ═══════════════════════════════════════════════════════════════

def print_summary(results):
    """打印效率 MAP 文字汇总."""
    speed = results['speed_grid']
    T_pts = results['T_grid']
    eta = results['eta_map']
    feas = results['feasible_map']

    # 峰值效率
    eta_valid = eta[feas]
    if len(eta_valid) == 0:
        print('\n⚠ 无可行工作点')
        return

    max_flat = np.nanargmax(eta)
    max_i, max_j = np.unravel_index(max_flat, eta.shape)
    print(f'\n{"=" * 70}')
    print(f'=== 效率 MAP 汇总 ===')
    print(f'峰值效率: {eta_valid.max():.1f}% '
          f'@ {speed[max_i]:.0f} rpm, {T_pts[max_j]:.0f} Nm')
    print(f'>98%: {np.sum(eta_valid>98)/len(eta_valid)*100:.1f}%  '
          f'>97%: {np.sum(eta_valid>97)/len(eta_valid)*100:.1f}%  '
          f'>95%: {np.sum(eta_valid>95)/len(eta_valid)*100:.1f}%  '
          f'>90%: {np.sum(eta_valid>90)/len(eta_valid)*100:.1f}%')

    # 效率表
    print(f'\n--- 效率 η (%) ---')
    hdr = f'{"n\\T":>7s}'
    for t in T_pts:
        hdr += f'{t:>8.0f}'
    print(hdr)
    for i, s in enumerate(speed):
        row = f'{s:7.0f}'
        for j in range(len(T_pts)):
            row += f'{eta[i,j]:8.1f}' if feas[i, j] else f'{"  --":>8s}'
        print(row)


def save_results(results, out_dir, prefix='efficiency_map'):
    """保存 NPZ + 报告 + 图表."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    vdc = results['vdc']
    imax = results['imax']

    # NPZ
    npz_path = out_dir / f'{prefix}.npz'
    np.savez_compressed(
        npz_path,
        speed_grid=results['speed_grid'],
        T_grid=results['T_grid'],
        eta_map=results['eta_map'],
        pfe_map=results['pfe_map'],
        pcu_map=results['pcu_map'],
        ploss_map=results['ploss_map'],
        pout_map=results['pout_map'],
        is_map=results['is_map'],
        id_map=results['id_map'],
        iq_map=results['iq_map'],
        vs_map=results['vs_map'],
        m_map=results['m_map'],
        pf_map=results['pf_map'],
        psi_map=results['psi_map'],
        feasible_map=results['feasible_map'],
        region_map=results['region_map'],
        T_max_global=results['T_max_global'],
        vdc=vdc, imax=imax, Vmax=results['Vmax'],
        Rs=results['Rs'], pole_pairs=results['pole_pairs'],
    )
    print(f'[D] 数据已保存: {npz_path}')

    # 图表
    png_path = out_dir / f'{prefix}_{vdc:.0f}V_{imax:.0f}A.png'
    plot_efficiency_map(results, str(png_path))

    # 报告
    eta = results['eta_map']
    feas = results['feasible_map']
    eta_valid = eta[feas]
    speed = results['speed_grid']
    T_pts = results['T_grid']
    max_flat = np.nanargmax(eta)
    max_i, max_j = np.unravel_index(max_flat, eta.shape)

    lines = [
        f'# PMSM 效率 MAP 报告',
        f'Vdc={vdc}V, Imax={imax}A, T_max={results["T_max_global"]:.1f} Nm',
        f'',
        f'## 峰值效率',
        f'- η = {eta_valid.max():.1f}% @ {speed[max_i]:.0f} rpm, '
        f'{T_pts[max_j]:.0f} Nm',
        f'- P_cu = {results["pcu_map"][max_i,max_j]:.0f} W, '
        f'P_fe = {results["pfe_map"][max_i,max_j]:.0f} W',
        f'',
        f'## 高效区占比（占可行域）',
    ]
    for th in [98, 97, 95, 93, 90]:
        pct = np.sum(eta_valid > th) / len(eta_valid) * 100
        lines.append(f'- >{th}%: {pct:.1f}%')

    report_path = out_dir / f'{prefix}_{vdc:.0f}V_{imax:.0f}A.md'
    report_path.write_text('\n'.join(lines), encoding='utf-8')
    print(f'[D] 报告已保存: {report_path}')

    return npz_path, png_path


# ═══════════════════════════════════════════════════════════════
# 主入口
# ═══════════════════════════════════════════════════════════════

def run(psi_dq_npz=None, project_path=None,
        speed_min=500, speed_max=8000, n_speed=20, n_torque=20,
        vdc=500, imax=300, Rs=0.05, pole_pairs=4,
        k_fe=DEFAULT_k_fe, f0=DEFAULT_f0,
        alpha=DEFAULT_alpha, beta=DEFAULT_beta,
        out_dir=None):
    """计算 PMSM 全域效率 MAP.

    依赖 Sub-flow A 输出的 psi_dq_table.npz。

    Parameters
    ----------
    psi_dq_npz : str or Path, optional
        Sub-flow A 输出的磁链表路径。
        如果为 None，从 project_path 下自动查找 psi_dq_table.npz。
    project_path : str or Path, optional
        项目目录（包含 psi_dq_table.npz）。
    speed_min, speed_max : float
        转速范围 (rpm)
    n_speed, n_torque : int
        转速/扭矩网格点数
    vdc : float
        直流母线电压 (V)
    imax : float
        峰值相电流限制 (A)
    Rs : float
        相电阻 (Ω)
    pole_pairs : int
        极对数
    k_fe, f0, alpha, beta : float
        铁耗模型参数:
        P_fe = k_fe * (f_e/f0)^alpha * (psi/psi_pm)^beta
    out_dir : str or Path, optional
        输出目录（默认与 psi_dq_npz 同目录）

    Returns
    -------
    dict — 完整的效率 MAP 结果
    """
    # ── 解析输入路径 ──
    if psi_dq_npz is None:
        if project_path is None:
            raise ValueError('必须提供 psi_dq_npz 或 project_path 之一')
        psi_dq_npz = Path(project_path) / 'psi_dq_table.npz'
    psi_path = Path(psi_dq_npz)
    if not psi_path.exists():
        raise FileNotFoundError(f'ψd/ψq 磁链表不存在: {psi_path}')

    if out_dir is None:
        out_dir = psi_path.parent
    else:
        out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f'╔══════════════════════════════════════════════════════════╗')
    print(f'║  PMSM 效率 MAP — Sub-flow D                            ║')
    print(f'╠══════════════════════════════════════════════════════════╣')
    print(f'║  输入: {psi_path.name}')
    print(f'║  Vdc={vdc}V, Imax={imax}A, Rs={Rs}Ω, P={pole_pairs}')
    print(f'║  转速: {speed_min}~{speed_max} rpm ({n_speed} 点)')
    print(f'║  扭矩: {n_torque} 点')
    print(f'╚══════════════════════════════════════════════════════════╝')

    # ── 计算 ──
    computer = EfficiencyMapComputer(str(psi_path))
    results = computer.compute(
        vdc=vdc, imax=imax,
        speed_min=speed_min, speed_max=speed_max,
        n_speed=n_speed, n_torque=n_torque,
        Rs=Rs, pole_pairs=pole_pairs,
        k_fe=k_fe, f0=f0, alpha=alpha, beta=beta)

    # ── 输出 ──
    print_summary(results)
    save_results(results, out_dir)

    return results


# ═══════════════════════════════════════════════════════════════
# 命令行入口
# ═══════════════════════════════════════════════════════════════

if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser(
        description='PMSM 效率 MAP (Sub-flow D) — 基于 ψd/ψq 磁链表')
    ap.add_argument('psi_dq_npz', nargs='?', default=None,
                    help='psi_dq_table.npz 路径')
    ap.add_argument('--project', '-p', default=None,
                    help='项目目录 (自动查找 psi_dq_table.npz)')
    ap.add_argument('--vdc', type=float, default=500)
    ap.add_argument('--imax', type=float, default=300)
    ap.add_argument('--speed-min', type=float, default=500)
    ap.add_argument('--speed-max', type=float, default=8000)
    ap.add_argument('--n-speed', type=int, default=20)
    ap.add_argument('--n-torque', type=int, default=20)
    ap.add_argument('--Rs', type=float, default=0.05)
    ap.add_argument('--k-fe', type=float, default=DEFAULT_k_fe)
    ap.add_argument('--out-dir', default=None)
    args = ap.parse_args()

    run(psi_dq_npz=args.psi_dq_npz, project_path=args.project,
        speed_min=args.speed_min, speed_max=args.speed_max,
        n_speed=args.n_speed, n_torque=args.n_torque,
        vdc=args.vdc, imax=args.imax, Rs=args.Rs,
        k_fe=args.k_fe, out_dir=args.out_dir)
