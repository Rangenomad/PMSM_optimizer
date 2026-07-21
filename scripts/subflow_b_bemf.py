"""Sub-flow B: 空载反电势 (Back EMF)
求解器: TransientXY (5_Partial_motor_TR)
方法: Imax=0 空载, 额定转速旋转, 提取三相感应电压波形 + FFT

用法:
    python -c "from scripts.subflow_b_bemf import run; run(rated_speed=3000)"
"""

import sys
sys.stdout.reconfigure(encoding='utf-8')
import re
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
_DEFAULT_TEMPLATE = str(ROOT / 'references' / 'Prius_2D_Practice.aedt')


def _sync_motion_angular_velocity(aedt_path, target_rpm):
    """开项目前同步 .aedt 文件中的 Angular Velocity 和 Speed_rpm 默认值。

    gRPC 下边界属性修改不持久化，且 m2d['Speed_rpm']='2000rpm' 带单位赋值
    会破坏绕组公式 Omega=360*speed_rpm*PolePairs/60 的表达式计算。

    此函数直接编辑 .aedt 文本文件，确保：
    1. Angular Velocity = '{target_rpm}rpm'（MotionSetup 机械角速度）
    2. Speed_rpm 默认值 = '{target_rpm}'（无单位，绕组公式依赖）
    """
    file = Path(aedt_path)
    if not file.exists():
        return
    content = file.read_text(encoding='utf-8')

    # 1. 同步 Angular Velocity（可能是硬编码值或变量引用）
    new_content, n_av = re.subn(
        r"'Angular Velocity'='[^']*'",
        f"'Angular Velocity'='{target_rpm}rpm'",
        content
    )

    # 2. 同步 Speed_rpm 默认值（强制无单位，避免绕组公式计算异常）
    new_content, n_sp = re.subn(
        r"VariableProp\('Speed_rpm', 'UD', '', '[^']*'",
        f"VariableProp('Speed_rpm', 'UD', '', '{target_rpm}'",
        new_content
    )

    if n_av > 0 or n_sp > 0:
        file.write_text(new_content, encoding='utf-8')
        print(f'  [B] .aedt 文件已同步: AV={target_rpm}rpm ({n_av}处), Speed_rpm={target_rpm} ({n_sp}处)')


def _plot_bemf(time_s, va, vb, vc, rated_speed, save_dir=None):
    """Plot three-phase BEMF waveform (time-domain only, no FFT).

    Parameters
    ----------
    time_s : np.ndarray
        Time vector (seconds).
    va, vb, vc : np.ndarray
        Phase A/B/C voltage (V).
    rated_speed : int
        Speed label for plot title and filename.
    save_dir : Path or str, optional
        Directory to save the plot. Defaults to cwd.
    """
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print('[B] 警告: matplotlib 未安装，跳过波形图生成')
        return

    fig, ax = plt.subplots(figsize=(10, 4.5))
    t_ms = time_s * 1000
    ax.plot(t_ms, va, label='Phase A', color='#e74c3c', lw=0.8)
    ax.plot(t_ms, vb, label='Phase B', color='#2ecc71', lw=0.8)
    ax.plot(t_ms, vc, label='Phase C', color='#3498db', lw=0.8)
    ax.set_xlabel('Time (ms)')
    ax.set_ylabel('Voltage (V)')
    ax.set_title(f'Back EMF @ {rated_speed} rpm')
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)

    ymin = min(va.min(), vb.min(), vc.min())
    ymax = max(va.max(), vb.max(), vc.max())
    margin = (ymax - ymin) * 0.1
    ax.set_ylim(ymin - margin, ymax + margin)

    plt.tight_layout()

    out_dir = Path(save_dir) if save_dir else Path.cwd()
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = str(out_dir / f'bemf_{rated_speed}rpm.png')
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'[B] 波形图已保存: {out_path}')


def run(rated_speed=3000, elec_periods=2, time_steps_per_cycle=50, project_path=None):
    from scripts.project_utils import get_template_path
    from scripts.param_guard import check_var
    from ansys.aedt.core import Maxwell2d

    if project_path:
        template = get_template_path(Path(project_path))
    else:
        template = _DEFAULT_TEMPLATE

    # 打开项目前同步 MotionSetup Angular Velocity 和 Speed_rpm 默认值
    # gRPC 无法持久化边界属性修改，且 m2d['Speed_rpm']='2000rpm' 带单位赋值会
    # 破坏绕组公式 Omega=360*speed_rpm*PolePairs/60 的表达式计算
    _sync_motion_angular_velocity(template, rated_speed)

    m2d = Maxwell2d(
        project=template, design='5_Partial_motor_TR',
        solution_type='TransientXY',
        non_graphical=False, new_desktop=False, close_on_exit=False
    )

    print(f'[B] 步骤 1/4: 设置空载工况 Imax=0, Speed={rated_speed}rpm')
    # 设置空载 (Imax=0)
    # Speed_rpm 已通过 .aedt 文件预设，不通过 PyAEDT 设置（避免单位破坏绕组公式）
    check_var('Imax', 'B')
    m2d['Imax'] = '0A'

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

    # 保存波形图
    _plot_bemf(time_vals, va, vb, vc, rated_speed,
               save_dir=Path(project_path) if project_path else None)

    m2d.close_project()
    return results


if __name__ == '__main__':
    run()
