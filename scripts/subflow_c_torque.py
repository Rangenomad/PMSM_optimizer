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
        print(f'  [C] .aedt 文件已同步: AV={target_rpm}rpm ({n_av}处), Speed_rpm={target_rpm} ({n_sp}处)')


def _mtpa_scan(m2d, rated_speed, pole_pairs):
    """MTPA 角度扫描: 候选角各跑 1 个电周期, 选扭矩最大者.

    候选角覆盖正负区间:
    - 负角 → 弱磁方向 (Lq > Ld 的传统 IPM 最优)
    - 正角 → 增磁方向 (Ld > Lq 的反凸极电机最优)
    """
    from scripts.param_guard import check_var
    candidate_angles = [-20, 0, 15, 30, 45, 55]
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
        check_var('Thet_deg', 'C')
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

    # 打开项目前同步 MotionSetup Angular Velocity 和 Speed_rpm 默认值
    # gRPC 无法持久化边界属性修改，且 m2d['Speed_rpm']='2000rpm' 带单位赋值会
    # 破坏绕组公式 Omega=360*speed_rpm*PolePairs/60 的表达式计算
    _sync_motion_angular_velocity(template, rated_speed)

    m2d = Maxwell2d(
        project=template, design='5_Partial_motor_TR',
        solution_type='TransientXY',
        non_graphical=False, new_desktop=False, close_on_exit=False
    )

    # 设置负载电流
    # Speed_rpm 已通过 .aedt 文件预设，不通过 PyAEDT 设置（避免单位破坏绕组公式）
    check_var('Imax', 'C')
    m2d['Imax'] = f'{rated_current}A'

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
