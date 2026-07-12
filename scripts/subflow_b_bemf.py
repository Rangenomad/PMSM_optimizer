"""Sub-flow B: 空载反电势 (Back EMF)
求解器: TransientXY (5_Partial_motor_TR)
方法: Imax=0 空载, 额定转速旋转, 提取三相电压波形 + FFT

用法:
    python -c "from scripts.subflow_b_bemf import run; run(rated_speed=3000)"
"""

import sys
sys.stdout.reconfigure(encoding='utf-8')
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
TEMPLATE = str(ROOT / 'references' / 'Prius_2D_Practice.aedt')


def compute_thd(signal):
    """计算 THD (Total Harmonic Distortion)"""
    fft_vals = np.fft.rfft(signal)
    fft_mag = np.abs(fft_vals)
    fundamental = fft_mag[1] if len(fft_mag) > 1 else fft_mag[0]
    if fundamental == 0:
        return 0
    harmonics = np.sqrt(np.sum(fft_mag[2:] ** 2))
    return harmonics / fundamental * 100


def run(rated_speed=3000, elec_periods=2, time_steps_per_cycle=200):
    from ansys.aedt.core import Maxwell2d

    m2d = Maxwell2d(
        project=TEMPLATE, design='5_Partial_motor_TR',
        solution_type='TransientXY',
        non_graphical=False, new_desktop=False, close_on_exit=False
    )

    # 设置空载 (Imax=0) 和转速
    m2d['Imax'] = '0A'
    m2d['Speed_rpm'] = f'{rated_speed}rpm'

    # 计算仿真时间
    pole_pairs = float(m2d['PolePairs'])
    freq = rated_speed / 60 * pole_pairs  # 电频率 (Hz)
    stop_time = elec_periods / freq
    time_step = 1 / (freq * time_steps_per_cycle)

    # 更新求解设置
    setup = m2d.setups[0]
    setup.props['StopTime'] = f'{stop_time}s'
    setup.props['TimeStep'] = f'{time_step}s'
    setup.update()

    m2d.analyze('Setup1')

    # 提取三相电压波形
    data = m2d.post.get_solution_data(
        expressions=['Voltage(Phase_A)', 'Voltage(Phase_B)', 'Voltage(Phase_C)'],
        variations=m2d.post.get_solution_data_variation()
    )

    time_vals = np.array(data.data('Time'))
    va = np.array(data.data('Voltage(Phase_A)'))
    vb = np.array(data.data('Voltage(Phase_B)'))
    vc = np.array(data.data('Voltage(Phase_C)'))

    # FFT 分析
    fs = 1.0 / time_step
    results = {}
    for name, v in [('A', va), ('B', vb), ('C', vc)]:
        fft_vals = np.fft.rfft(v)
        fft_mag = np.abs(fft_vals)
        n = len(fft_vals)
        freqs = np.fft.rfftfreq(len(v), d=time_step)

        # 基波幅值 (第一个非直流分量)
        fundamental_idx = np.argmax(fft_mag[1:]) + 1
        fundamental_mag = fft_mag[fundamental_idx] * 2 / len(v)
        thd = compute_thd(v)

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
