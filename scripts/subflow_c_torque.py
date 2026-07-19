"""Sub-flow C: 额定点扭矩
求解器: TransientXY (5_Partial_motor_TR)
方法: 额定负载, 额定转速, 提取 Moving1.Torque 波形
      电流角未指定时自动 MTPA 搜索

用法:
    # 指定电流角
    python -c "from scripts.subflow_c_torque import run; run(rated_current=250, current_angle=-30, rated_speed=3000)"
    # 自动 MTPA
    python -c "from scripts.subflow_c_torque import run; run(rated_current=250, rated_speed=3000)"
"""

import sys
sys.stdout.reconfigure(encoding='utf-8')
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
_DEFAULT_TEMPLATE = str(ROOT / 'references' / 'Prius_2D_Practice.aedt')


def _mtpa_scan(m2d, rated_speed, pole_pairs):
    """MTPA 角度扫描: 候选角各跑 1 个电周期, 选扭矩最大者"""
    candidate_angles = [0, -15, -25, -35, -45, -60]
    freq = rated_speed / 60 * pole_pairs

    # 粗网格: 100 步/周期
    setup = m2d.setups[0]
    setup.props['StopTime'] = f'{1/freq}s'
    setup.props['TimeStep'] = f'{1/(freq*50)}s'
    setup.update()

    best_T = -1e9
    best_angle = candidate_angles[0]

    for idx, angle in enumerate(candidate_angles):
        print(f'  [C-MTPA] {idx+1}/{len(candidate_angles)} 候选角 θ={angle}°, 正在求解...')
        m2d['Thet_deg'] = str(angle)  # 不用 ° 后缀，gRPC 赋值第二次后静默失败
        m2d.save_project()
        m2d.analyze('Setup1')

        data = m2d.post.get_solution_data_per_variation(
            expressions=['Moving1.Torque']
        )
        torque = np.array(data.data_real('Moving1.Torque'))
        time_vals = np.array(data.primary_sweep_values)
        dt_actual = np.mean(np.diff(time_vals))
        steps_per_period_actual = int(round((1/freq) / dt_actual))
        # 最后半周期稳态平均
        half_period = max(steps_per_period_actual // 2, 1)
        T_avg = np.mean(torque[-half_period:])
        print(f'    θ={angle:3d}° → T={T_avg:.2f} Nm')

        if T_avg > best_T:
            best_T = T_avg
            best_angle = angle

    print(f'  MTPA 最优角: θ={best_angle}° (T={best_T:.2f} Nm)')
    return best_angle


def run(rated_current=250, current_angle=None, rated_speed=3000, elec_periods=3, project_path=None):
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

    # 设置负载和转速
    check_var('Imax', 'C')
    m2d['Imax'] = f'{rated_current}A'
    check_var('Speed_rpm', 'C')
    m2d['Speed_rpm'] = f'{rated_speed}rpm'
    pole_pairs = float(m2d['Poles']) / 2

    # MTPA 自动搜索（未指定电流角时）
    auto_mtpa = current_angle is None
    if auto_mtpa:
        print('  [C] 未指定电流角, 正在搜索 MTPA 最优角...')
        current_angle = _mtpa_scan(m2d, rated_speed, pole_pairs)

    check_var('Thet_deg', 'C')
    m2d['Thet_deg'] = str(current_angle)  # 整数, 不用 ° 后缀

    # 完整仿真
    freq = rated_speed / 60 * pole_pairs
    stop_time = elec_periods / freq
    time_step = 1 / (freq * 50)  # 50 steps/period

    print(f'[C] 步骤 1/3: 设置工况 Imax={rated_current}A, θ={current_angle}°, Speed={rated_speed}rpm')
    print(f'[C] 步骤 2/3: 正在求解... (电频率 {freq:.1f}Hz, {elec_periods} 个电周期)')
    setup = m2d.setups[0]
    setup.props['StopTime'] = f'{stop_time}s'
    setup.props['TimeStep'] = f'{time_step}s'
    setup.update()

    m2d.analyze('Setup1')

    print(f'[C] 步骤 3/3: 提取扭矩波形 + 纹波分析')

    # 提取扭矩波形
    data = m2d.post.get_solution_data_per_variation(
        expressions=['Moving1.Torque']
    )
    torque = np.array(data.data_real('Moving1.Torque'))
    time_vals = np.array(data.primary_sweep_values)  # 实际时间轴

    # 根据实际时间步计算每周期步数（更可靠，不依赖预设值）
    dt_actual = np.mean(np.diff(time_vals))  # 实际时间步长
    T_period = 1 / freq                     # 电周期
    steps_per_period_actual = int(round(T_period / dt_actual))

    # 取最后一个完整周期的稳态数据
    if len(torque) >= steps_per_period_actual:
        steady_torque = torque[-steps_per_period_actual:]
    else:
        steady_torque = torque  # 数据不足时全部使用

    avg_torque = np.mean(steady_torque)
    max_torque = np.max(steady_torque)
    min_torque = np.min(steady_torque)
    ripple_pp = (max_torque - min_torque) / avg_torque * 100

    print(f'\n{"=" * 60}')
    print(f'=== 额定点扭矩结果 ===')
    print(f'平均扭矩: {avg_torque:.2f} Nm')
    print(f'扭矩脉动: {ripple_pp:.2f}% (峰峰值)')
    print(f'最大扭矩: {max_torque:.2f} Nm')
    print(f'最小扭矩: {min_torque:.2f} Nm')
    print(f'电流角: {current_angle}° ({"自动 MTPA" if auto_mtpa else "用户指定"})')
    print(f'转速: {rated_speed} rpm')
    print(f'电流: {rated_current} A')
    print(f'{"=" * 60}')

    m2d.close_project()
    return {
        'current_angle': current_angle,
        'avg_torque': avg_torque,
        'ripple_pp': ripple_pp,
        'max_torque': max_torque,
        'min_torque': min_torque,
    }


if __name__ == '__main__':
    run()
