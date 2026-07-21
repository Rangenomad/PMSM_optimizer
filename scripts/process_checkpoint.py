"""从 checkpoint 加载 Sub-flow D 结果，重采样到 speed×torque 网格并输出 MAP"""
import sys
sys.stdout.reconfigure(encoding='utf-8')
import json, numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

# 加载 checkpoint
ckpt_file = ROOT / 'pmsm_projects/2026-07-19_efficiency_map/_subflow_d_checkpoint.json'
ckpt = json.loads(ckpt_file.read_text(encoding='utf-8'))

speed_points = np.array(ckpt['speed_points'])
Imax_points = np.array(ckpt['Imax_points'])
n_speeds = len(speed_points)
n_is = len(Imax_points)

# 重建 results 数组
results = np.full((n_speeds, n_is), None, dtype=object)
for entry in ckpt['results']:
    results[entry['i'], entry['j']] = entry['data']

# 统计有效点
valid = sum(1 for r in results.flat if r is not None)
zero_tq = sum(1 for r in results.flat if r and r.get('T_avg', 0) < 0.1 and r.get('Is', 0) > 1)
print(f'加载 {valid} 个已计算点')
print(f'其中 Is>0 且 T≈0: {zero_tq} 个 (求解失败)')

# 打印原始数据摘要
print(f'\n{"="*70}')
print('原始 FEA 数据 (speed × Is)')
print(f'{"="*70}')
print(f'{"n\\Is":>6s}', end='')
for isv in Imax_points:
    print(f'{isv:>8.0f}', end='')
print(f'  {"":>6s}')
for i, n in enumerate(speed_points):
    print(f'{n:6.0f}', end='')
    for j, isv in enumerate(Imax_points):
        r = results[i, j]
        if r and r.get('T_avg', 0) > 0.1:
            print(f'{r["T_avg"]:8.1f}', end='')
        elif r and r.get('Is', 0) < 1:
            print(f'{"0":>8s}', end='')
        else:
            print(f'{"-":>8s}', end='')
    print()

print(f'\n{"n\\Is":>6s}', end='')
for isv in Imax_points:
    print(f'{isv:>8.0f}', end='')
print(f'  {"η%":>6s}')
for i, n in enumerate(speed_points):
    print(f'{n:6.0f}', end='')
    for j, isv in enumerate(Imax_points):
        r = results[i, j]
        if r and r.get('eta', 0) > 0.1:
            print(f'{r["eta"]:8.1f}', end='')
        elif r and r.get('Is', 0) < 1:
            print(f'{"0":>8s}', end='')
        else:
            print(f'{"-":>8s}', end='')
    print()

# ── 重采样到 speed × torque ──
sys.path.insert(0, str(ROOT))
from scripts.subflow_d_efficiency_map import _resample_to_torque_grid

resampled, T_points = _resample_to_torque_grid(results, speed_points, len(Imax_points))

if resampled is not None:
    T_labels = [f'{t:.0f}' for t in T_points]

    print(f'\n{"="*70}')
    print(f'=== 重采样 MAP (speed × 扭矩) ===')
    print(f'行: {len(speed_points)} 转速, 列: {len(T_points)} 扭矩')
    print(f'扭矩范围: 0~{T_points[-1]:.0f} Nm')
    print(f'{"="*70}')

    def _print_map(title, key, fmt='8.2f'):
        print(f'\n--- {title} ---')
        print(f'{"n\\T":>6s}', end='')
        for tl in T_labels:
            print(f'{tl:>8s}', end='')
        print()
        for i, n in enumerate(speed_points):
            print(f'{n:6.0f}', end='')
            for j in range(len(T_points)):
                r = resampled[i, j]
                v = r[key] if r else float('nan')
                if isinstance(v, float) and np.isnan(v):
                    print(f'{"-":>8s}', end='')
                else:
                    print(f'{v:{fmt}}', end='')
            print()

    _print_map('效率 η (%)', 'eta', '8.1f')
    _print_map('FEA 扭矩 T_avg (Nm)', 'T_avg', '8.1f')
    _print_map('电流 Is (A)', 'Is', '8.1f')
    _print_map('调制比 M', 'M', '8.3f')
    _print_map('功率因数 PF', 'PF', '8.3f')
    _print_map('相电压 Vs (V)', 'Vs', '8.0f')

    eta_vals = [r['eta'] for r in resampled.flat if r and r.get('eta', 0) > 0.1]
    if eta_vals:
        print(f'\n峰值效率: {np.max(eta_vals):.1f}%')

    # ── 保存 CSV ──
    out_dir = ROOT / 'results'
    out_dir.mkdir(exist_ok=True)
    col_labels = ',' + ','.join(T_labels)

    csv_maps = {
        'efficiency_map': ('eta (%)', 'eta'),
        'torque_map': ('T_avg (Nm)', 'T_avg'),
        'current_map': ('Is (A)', 'Is'),
        'modulation_map': ('M', 'M'),
        'pf_map': ('PF', 'PF'),
        'voltage_map': ('Vs (V)', 'Vs'),
    }
    for fname, (hdr, key) in csv_maps.items():
        rows = [f'# {hdr}  ─ 行:转速(rpm), 列:扭矩(Nm), Thet=45°']
        rows.append(col_labels)
        for i, n in enumerate(speed_points):
            vals = []
            for j in range(len(T_points)):
                r = resampled[i, j]
                vals.append(str(r[key]) if r else '')
            rows.append(f'{n:.0f},{",".join(vals)}')
        (out_dir / f'{fname}.csv').write_text('\n'.join(rows), encoding='utf-8')

    print(f'\n结果已保存到 results/*.csv')
