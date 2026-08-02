"""Sub-flow F: Report Generation

Generates:
  1. External characteristic CSV table (speed, torque, power)
  2. Loss table CSV (speed, torque, loss, Udc, Iac_rms, M, PF)
  3. Simulation report (markdown) with embedded charts and table links

Usage:
    python -c "from scripts.subflow_f_report import run; run(project_path='pmsm_projects/xxx')"
"""

import sys
sys.stdout.reconfigure(encoding='utf-8')
from pathlib import Path
from datetime import datetime
import numpy as np

ROOT = Path(__file__).resolve().parent.parent


# ═══════════════════════════════════════════════════════════════
#  Table 1: External Characteristic
# ═══════════════════════════════════════════════════════════════

def generate_external_char_table(project_path):
    """Generate external characteristic CSV table from outer_characteristic.npz.

    Columns: Speed (rpm), Torque (Nm), Power (kW)
    """
    proj = Path(project_path)
    npz_path = proj / 'outer_characteristic.npz'
    if not npz_path.exists():
        raise FileNotFoundError(f'outer_characteristic.npz not found in {project_path}')

    data = np.load(npz_path, allow_pickle=True)
    speeds = data['speed']
    torques = data['T_max']
    powers = data['P_max_kW']

    # Filter valid points (torque > 0.1 Nm)
    valid = torques > 0.1

    csv_path = proj / 'outer_characteristic_table.csv'
    with open(csv_path, 'w', encoding='utf-8-sig') as f:
        f.write('Speed (rpm),Torque (Nm),Power (kW)\n')
        for s, t, p in zip(speeds[valid], torques[valid], powers[valid]):
            f.write(f'{s:.0f},{t:.1f},{p:.2f}\n')

    n_pts = valid.sum()
    print(f'[F] Table 1 (External Characteristic): {n_pts} points → {csv_path.name}')
    return csv_path


# ═══════════════════════════════════════════════════════════════
#  Table 2: Loss Table
# ═══════════════════════════════════════════════════════════════

def generate_loss_table(project_path):
    """Generate loss table CSV from efficiency_map.npz.

    Columns: Speed (rpm), Torque (Nm), Loss (W), Udc (V), Iac_rms (A), M, PF
    Covers all feasible operating points in the efficiency map grid.
    """
    proj = Path(project_path)
    npz_path = proj / 'efficiency_map.npz'
    if not npz_path.exists():
        raise FileNotFoundError(f'efficiency_map.npz not found in {project_path}')

    data = np.load(npz_path, allow_pickle=True)
    speed_grid = data['speed_grid']
    T_grid = data['T_grid']
    ploss_map = data['ploss_map']
    is_map = data['is_map']
    m_map = data['m_map']
    pf_map = data['pf_map']
    feasible = data['feasible_map']
    vdc = float(data['vdc'])

    csv_path = proj / 'loss_table.csv'
    with open(csv_path, 'w', encoding='utf-8-sig') as f:
        f.write('Speed (rpm),Torque (Nm),Loss (W),Udc (V),Iac_rms (A),M,PF\n')
        for i, speed in enumerate(speed_grid):
            for j, torque in enumerate(T_grid):
                if feasible[i, j]:
                    loss = ploss_map[i, j]
                    Is = is_map[i, j]
                    Iac_rms = Is / np.sqrt(2)
                    M = m_map[i, j]
                    PF = pf_map[i, j]
                    f.write(f'{speed:.0f},{torque:.1f},{loss:.1f},{vdc:.0f},'
                            f'{Iac_rms:.1f},{M:.3f},{PF:.3f}\n')

    n_pts = feasible.sum()
    print(f'[F] Table 2 (Loss Table): {n_pts} points → {csv_path.name}')
    return csv_path


# ═══════════════════════════════════════════════════════════════
#  Simulation Report (Markdown)
# ═══════════════════════════════════════════════════════════════

def _extract_summary(project_path):
    """Extract key metrics from project data files."""
    proj = Path(project_path)
    summary = {}

    # ── External characteristic ──
    oc_path = proj / 'outer_characteristic.npz'
    if oc_path.exists():
        oc = np.load(oc_path, allow_pickle=True)
        summary['vdc'] = float(oc['vdc'])
        summary['imax'] = float(oc['imax'])
        speeds = oc['speed']
        torques = oc['T_max']
        powers = oc['P_max_kW']
        regions = oc['region'] if 'region' in oc else None

        valid = torques > 0.1
        if valid.any():
            peak_idx = np.argmax(torques[valid])
            summary['peak_T'] = float(torques[valid][peak_idx])
            summary['peak_P'] = float(np.max(powers[valid]))

            # FW corner: first point where region == 'FW'
            if regions is not None:
                for k, r in enumerate(regions):
                    if r == 'FW':
                        summary['fw_speed'] = float(speeds[k])
                        break
                if 'fw_speed' not in summary:
                    summary['fw_speed'] = None

            # Torque taper at max speed
            last_T = torques[valid][-1]
            summary['taper'] = last_T / summary['peak_T'] * 100 if summary['peak_T'] > 0 else 0
            summary['last_speed'] = float(speeds[valid][-1])
            summary['last_T'] = float(last_T)

    # ── Efficiency map ──
    em_path = proj / 'efficiency_map.npz'
    if em_path.exists():
        em = np.load(em_path, allow_pickle=True)
        eta_map = em['eta_map']
        feas = em['feasible_map']
        eta_valid = eta_map[feas]

        if len(eta_valid) > 0:
            summary['peak_eta'] = float(np.nanmax(eta_valid))
            max_flat = np.nanargmax(eta_map)
            max_i, max_j = np.unravel_index(max_flat, eta_map.shape)
            summary['eta_speed'] = float(em['speed_grid'][max_i])
            summary['eta_torque'] = float(em['T_grid'][max_j])
            summary['eta_gt97'] = np.sum(eta_valid > 97) / len(eta_valid) * 100
            summary['eta_gt95'] = np.sum(eta_valid > 95) / len(eta_valid) * 100
            summary['eta_gt90'] = np.sum(eta_valid > 90) / len(eta_valid) * 100
            summary['n_feasible'] = len(eta_valid)
            summary['global_T_max'] = float(em['T_max_global'])

    # ── Stack length from psi_dq_table ──
    psi_path = proj / 'psi_dq_table.npz'
    if psi_path.exists():
        psi = np.load(psi_path, allow_pickle=True)
        if 'stack_length' in psi:
            sl = psi['stack_length']
            summary['stack_length'] = f'{sl:.2f}' if sl is not None else '83.82 (default)'
        else:
            summary['stack_length'] = '83.82 (default)'
        if 'symmetry_multiplier' in psi:
            summary['sym_mult'] = int(psi['symmetry_multiplier'])
        else:
            summary['sym_mult'] = 8

    return summary


def generate_report(project_path):
    """Generate simulation report in markdown format."""
    proj = Path(project_path)
    s = _extract_summary(project_path)

    # ── Find chart images ──
    oc_png = None
    em_png = None
    for f in sorted(proj.glob('*.png')):
        name = f.name.lower()
        if 'outer_char' in name and oc_png is None:
            oc_png = f
        elif ('efficiency_map' in name or 'efficiency' in name) and em_png is None:
            em_png = f

    L = []  # report lines

    def w(line=''):
        L.append(line)

    w('# PMSM 仿真报告')
    w()
    w(f'**生成时间**: {datetime.now().strftime("%Y-%m-%d %H:%M")}')
    w()

    # ── 1. 仿真目的 ──
    w('## 1. 仿真目的')
    w()
    w('基于 Prius 2D 模板，对永磁同步电机进行全工况电磁性能评估：')
    w()
    w('- **Sub-flow A**: MTPA 标定 — FEA 磁链表 (Id×Iq 扫描) + MTPA 优化')
    w('- **Sub-flow E**: 外特性 T-n 曲线 — 纯解析，基于磁链表')
    w('- **Sub-flow D**: 全域效率 MAP — 纯解析，含 Steinmetz 铁耗模型')
    w()

    # ── 2. 仿真设置 ──
    w('## 2. 仿真设置')
    w()
    w('| 参数 | 数值 |')
    w('|------|------|')
    w(f'| 直流母线电压 (Vdc) | {s.get("vdc", "—"):.0f} V |' if 'vdc' in s else '| 直流母线电压 | — |')
    w(f'| 峰值相电流 (Imax) | {s.get("imax", "—"):.0f} A |' if 'imax' in s else '| 峰值相电流 | — |')
    w(f'| 叠片长度 (StackLength) | {s.get("stack_length", "83.82")} mm |')
    w(f'| 极对数 | 4 |')
    w(f'| 相电阻 (Rs) | 0.05 Ω |')
    w(f'| 铁耗模型 | Steinmetz: k_fe=300W @ 200Hz, α=1.6, β=1.8 |')
    w(f'| FEA 求解器 | MagnetostaticXY |')
    w(f'| FEA 网格 | 20 Id × 15 Iq = 300 点 |')
    w(f'| 对称倍数 | {s.get("sym_mult", 8)}× (1/8 局部模型 → 完整电机) |')
    w()

    # ── 3. 仿真结果 ──
    w('## 3. 仿真结果')
    w()

    # 3.1 外特性
    w('### 3.1 外特性 (T-n 曲线)')
    w()
    if 'peak_T' in s:
        w('| 指标 | 数值 |')
        w('|------|------|')
        w(f'| 峰值扭矩 | {s["peak_T"]:.1f} Nm |')
        w(f'| 峰值功率 | {s["peak_P"]:.2f} kW |')
        if s.get('fw_speed'):
            w(f'| 弱磁拐点 | {s["fw_speed"]:.0f} rpm |')
        w(f'| 最高转速扭矩 | {s.get("last_T", 0):.1f} Nm @ {s.get("last_speed", 0):.0f} rpm ({s.get("taper", 0):.1f}% 下垂) |')
        w()

    if oc_png:
        w(f'![外特性]({oc_png.name})')
        w()

    # 3.2 效率 MAP
    w('### 3.2 效率 MAP')
    w()
    if 'peak_eta' in s:
        w('| 指标 | 数值 |')
        w('|------|------|')
        w(f'| 峰值效率 | {s["peak_eta"]:.1f}% @ {s["eta_speed"]:.0f} rpm, {s["eta_torque"]:.0f} Nm |')
        w(f'| >97% 高效区占比 | {s.get("eta_gt97", 0):.1f}% |')
        w(f'| >95% 高效区占比 | {s.get("eta_gt95", 0):.1f}% |')
        w(f'| >90% 高效区占比 | {s.get("eta_gt90", 0):.1f}% |')
        w(f'| 全局最大扭矩 | {s.get("global_T_max", 0):.1f} Nm |')
        w(f'| 可行工作点数 | {s.get("n_feasible", 0)} |')
        w()

    if em_png:
        w(f'![效率MAP]({em_png.name})')
        w()

    # ── 4. 数据表 ──
    w('## 4. 数据表')
    w()
    w(f'- **[外特性表](outer_characteristic_table.csv)** — 各转速下的最大扭矩和功率')
    w(f'- **[损耗表](loss_table.csv)** — 全工况网格的损耗、电流、调制比、功率因数')
    w()

    # ── 5. 备注 ──
    w('## 5. 备注')
    w()
    w('- 电机模型为 1/8 局部模型 (6 coils)，采用 Master/Slave 反周期边界条件。')
    w('- 磁链数据自动缩放 8× 等效为完整电机值。')
    w('- 铁耗为 Steinmetz 解析估计，不含 PWM 谐波损耗和磁钢涡流损耗。')
    w('- 效率绝对值可能偏高 1-3%（因未建模的损耗机制）。')
    w(f'- 主管线: ① Sub-flow A (FEA) → ② Sub-flow E (外特性) → ③ Sub-flow D (效率MAP) → ④ Sub-flow F (报告)')
    w()

    report_path = proj / 'simulation_report.md'
    report_path.write_text('\n'.join(L), encoding='utf-8')
    print(f'[F] Report: {report_path.name}')
    return report_path


# ═══════════════════════════════════════════════════════════════
#  Main Entry Point
# ═══════════════════════════════════════════════════════════════

def run(project_path):
    """Generate all reports and tables for a project.

    Requires outer_characteristic.npz and efficiency_map.npz in the project directory
    (outputs of Sub-flow E and Sub-flow D respectively).

    Parameters
    ----------
    project_path : str or Path
        Project directory containing the .npz result files.

    Returns
    -------
    dict with keys: table_ext, table_loss, report
    """
    proj = Path(project_path)
    if not proj.exists():
        raise FileNotFoundError(f'Project directory not found: {project_path}')

    print(f'╔══════════════════════════════════════════════════════════╗')
    print(f'║  PMSM Report Generation — Sub-flow F                    ║')
    print(f'╠══════════════════════════════════════════════════════════╣')
    print(f'║  Project: {proj.name}')
    print(f'╚══════════════════════════════════════════════════════════╝')
    print()

    # 1. External characteristic table
    csv1 = str(generate_external_char_table(project_path))

    # 2. Loss table
    csv2 = str(generate_loss_table(project_path))

    # 3. Simulation report
    report = str(generate_report(project_path))

    print(f'\n[F] Sub-flow F complete:')
    print(f'    Table 1 (外特性): {Path(csv1).name}')
    print(f'    Table 2 (损耗):   {Path(csv2).name}')
    print(f'    Report:           {Path(report).name}')

    return {'table_ext': csv1, 'table_loss': csv2, 'report': report}


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser(description='PMSM Report Generation (Sub-flow F)')
    ap.add_argument('project_path', help='Project directory path')
    args = ap.parse_args()
    run(args.project_path)
