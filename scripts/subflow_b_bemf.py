"""Sub-flow B: 空载反电势 (Back EMF)
求解器: TransientXY (5_Partial_motor_TR)
方法: Imax=0 空载, 额定转速旋转, 提取三相感应电压波形 + FFT

用法:
    python -c "from scripts.subflow_b_bemf import run; run(rated_speed=3000)"
"""

import sys
sys.stdout.reconfigure(encoding='utf-8')
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
_DEFAULT_TEMPLATE = str(ROOT / 'references' / 'Prius_2D_Practice.aedt')


def run(rated_speed=3000, elec_periods=2, time_steps_per_cycle=50, project_path=None):
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

    print(f'[B] 步骤 1/4: 设置空载工况 Imax=0, Speed={rated_speed}rpm')
    # 设置空载 (Imax=0) 和转速
    m2d['Imax'] = '0A'
    check_var('Speed_rpm', 'B')
    m2d['Speed_rpm'] = f'{rated_speed}rpm'

    # 计算仿真时间 (PolePairs = Poles/2)
    pole_pairs = float(m2d['Poles']) / 2
    freq = rated_speed / 60 * pole_pairs  # 电频率 (Hz)
    stop_time = elec_periods / freq
    time_step = 1 / (freq * time_steps_per_cycle)

    # 更新求解设置
    print(f'[B] 步骤 2/4: 更新求解设置 (StopTime={stop_time:.4f}s, TimeStep={time_step:.6f}s)')
    setup = m2d.setups[0]
    setup.props['StopTime'] = f'{stop_time}s'
    setup.props['TimeStep'] = f'{time_step}s'
    setup.update()

    print(f'[B] 步骤 3/4: 正在求解... (电频率 {freq:.1f} Hz, {elec_periods} 个电周期)')
    m2d.analyze('Setup1')

    print(f'[B] 步骤 4/4: 提取结果 + FFT 分析')
    # 提取三相感应电压 (InducedVoltage)
    data = m2d.post.get_solution_data_per_variation(
        expressions=['InducedVoltage(Phase_A)', 'InducedVoltage(Phase_B)', 'InducedVoltage(Phase_C)']
    )

    # Transient 主扫描轴为 Time (单位: ns)
    time_ns = np.array(data.primary_sweep_values, dtype=float)
    time_vals = time_ns * 1e-9  # 转换为秒
    va = np.array(data.data_real('InducedVoltage(Phase_A)')) / 1000  # mV → V
    vb = np.array(data.data_real('InducedVoltage(Phase_B)')) / 1000
    vc = np.array(data.data_real('InducedVoltage(Phase_C)')) / 1000

    # FFT 分析 (一次性计算, 避免重复 FFT)
    results = {}
    for name, v in [('A', va), ('B', vb), ('C', vc)]:
        fft_vals = np.fft.rfft(v)
        fft_mag = np.abs(fft_vals)
        freqs = np.fft.rfftfreq(len(v), d=time_step)

        # 基波幅值 (找峰值, 排除 DC)
        fundamental_idx = np.argmax(fft_mag[1:]) + 1
        fundamental_raw = fft_mag[fundamental_idx]
        fundamental_mag = fundamental_raw * 2 / len(v)

        # THD: 排除 DC(索引0) 和基波, 其余均为谐波
        mask = np.ones(len(fft_mag), dtype=bool)
        mask[0] = False       # 排除 DC
        mask[fundamental_idx] = False  # 排除基波
        harmonics = np.sqrt(np.sum(fft_mag[mask] ** 2))
        thd = harmonics / fundamental_raw * 100

        results[name] = {
            'fundamental': fundamental_mag,
            'thd': thd,
            'freq': freqs[fundamental_idx],
        }

    # 线反电势常数 Ke (V/krpm)
    ke_line = results['A']['fundamental'] / (rated_speed / 1000) * np.sqrt(3)

    print(f'\n{"=" * 60}')
    print(f'=== 空载反电势结果 ===')
    print(f'额定转速: {rated_speed} rpm')
    print(f'频率: {results["A"]["freq"]:.1f} Hz')
    print(f'{"=" * 60}')
    for name in ['A', 'B', 'C']:
        r = results[name]
        print(f'{name} 相: 基波幅值={r["fundamental"]:.2f} V, THD={r["thd"]:.2f}%')
    print(f'线反电势常数 Ke = {ke_line:.4f} V/(krpm)')
    print(f'{"=" * 60}')

    m2d.close_project()
    return results


if __name__ == '__main__':
    run()
