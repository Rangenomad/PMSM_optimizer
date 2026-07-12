"""Sub-flow D: 全域工作特性 MAP
求解器: 批量 TransientXY (5_Partial_motor_TR)
方法: (转速 × 扭矩) 二维网格扫描

用法:
    python -c "from scripts.subflow_d_efficiency_map import run; run()"
"""

import sys
sys.stdout.reconfigure(encoding='utf-8')
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
TEMPLATE = str(ROOT / 'references' / 'Prius_2D_Practice.aedt')


def run(speed_min=500, speed_max=10000, speed_steps=10, torque_steps=10, vdc=300, imax=250):
    from ansys.aedt.core import Maxwell2d

    speed_points = np.linspace(speed_min, speed_max, speed_steps)
    # 扭矩网格从 0 到估算最大扭矩（通过电流和永磁体估算）
    # 实际应基于 Sub-flow A/C 结果校准; 这里用占位比例
    torque_max = 300  # Nm (approximate for Prius motor)
    torque_points = np.linspace(0, torque_max, torque_steps)

    # results[n_speed, n_torque] = [P_loss, Is, PF, M]
    results = np.zeros((len(speed_points), len(torque_points), 4))

    m2d = Maxwell2d(
        project=TEMPLATE, design='5_Partial_motor_TR',
        solution_type='TransientXY',
        non_graphical=False, new_desktop=False, close_on_exit=False
    )

    Rs = 0.0  # 相电阻 (Ω), 默认忽略, 高速时可忽略

    for i, n in enumerate(speed_points):
        for j, T_req in enumerate(torque_points):
            print(f'  转速={n:.0f} rpm, 扭矩={T_req:.0f} Nm...')

            # 设置转速
            m2d['Speed_rpm'] = f'{n}rpm'

            # 估算工作点 (简化: 线性比例估算 Is 和电流角)
            # 如需精确 Id/Iq → 集成 Sub-flow A 的 Ld/Lq MAP 反算
            Id_approx = -0.5 * T_req / torque_max * imax
            Iq_approx = 0.8 * T_req / torque_max * imax
            Is = np.sqrt(Id_approx**2 + Iq_approx**2)

            m2d['Imax'] = f'{Is}A'
            m2d['Thet_deg'] = f'{np.degrees(np.arctan2(Id_approx, Iq_approx)):.1f}°'

            # 仿真一个电周期
            pole_pairs = float(m2d['PolePairs'])
            freq = n / 60 * pole_pairs
            setup = m2d.setups[0]
            setup.props['StopTime'] = f'{1/freq}s'
            setup.props['TimeStep'] = f'{1/(freq*200)}s'
            setup.update()

            m2d.analyze('Setup1')

            # 提取结果
            data = m2d.post.get_solution_data(
                expressions=['Moving1.Torque', 'InputVoltage(Phase_A)',
                             'InputCurrent(Phase_A)'],
                variations=m2d.post.get_solution_data_variation()
            )

            torque_wave = np.array(data.data('Moving1.Torque'))
            T_avg = np.mean(torque_wave[-200:])

            # 铜损
            Pcu = 3 * Is**2 * Rs

            # 电压幅值 (从 InputVoltage 波形提取)
            v_wave = np.array(data.data('InputVoltage(Phase_A)'))
            Vs = np.max(np.abs(v_wave))

            # 功率因数 (简化: 从电压电流相位差计算)
            i_wave = np.array(data.data('InputCurrent(Phase_A)'))
            # 零交叉法估算功率因数角
            v_sign = np.sign(v_wave - np.mean(v_wave))
            i_sign = np.sign(i_wave - np.mean(i_wave))
            PF = np.mean(v_sign * i_sign)  # 简化功率因数

            # 调制比
            M = Vs / (vdc / np.sqrt(3))

            results[i, j] = [Pcu, Is, PF, M]

    # 输出 4 张 MAP
    print('\n' + '=' * 60)
    print('=== 损耗 MAP P_loss (W) ===')
    print(f'{"n\\T":>8s}', end='')
    for T_v in torque_points:
        print(f'{T_v:8.0f}', end='')
    print()
    for i, n in enumerate(speed_points):
        print(f'{n:8.0f}', end='')
        for j in range(len(torque_points)):
            print(f'{results[i,j,0]:8.1f}', end='')
        print()

    print('\n=== Is MAP (A) ===')
    print(f'{"n\\T":>8s}', end='')
    for T_v in torque_points:
        print(f'{T_v:8.0f}', end='')
    print()
    for i, n in enumerate(speed_points):
        print(f'{n:8.0f}', end='')
        for j in range(len(torque_points)):
            print(f'{results[i,j,1]:8.1f}', end='')
        print()

    print('\n=== 功率因数 MAP ===')
    for i, n in enumerate(speed_points):
        print(f'{n:8.0f}', end='')
        for j in range(len(torque_points)):
            print(f'{results[i,j,2]:8.3f}', end='')
        print()

    print('\n=== 调制比 MAP ===')
    for i, n in enumerate(speed_points):
        print(f'{n:8.0f}', end='')
        for j in range(len(torque_points)):
            print(f'{results[i,j,3]:8.3f}', end='')
        print()

    # 保存 CSV
    np.savetxt(ROOT / 'results' / 'loss_map.csv', results[:, :, 0], delimiter=',',
               header='Loss MAP (W)', fmt='%.1f')
    np.savetxt(ROOT / 'results' / 'is_map.csv', results[:, :, 1], delimiter=',',
               header='Is MAP (A)', fmt='%.1f')
    np.savetxt(ROOT / 'results' / 'pf_map.csv', results[:, :, 2], delimiter=',',
               header='Power Factor MAP', fmt='%.3f')
    np.savetxt(ROOT / 'results' / 'modulation_map.csv', results[:, :, 3], delimiter=',',
               header='Modulation MAP', fmt='%.3f')
    print('\n结果已保存到 results/')

    m2d.close_project()
    return results


if __name__ == '__main__':
    run()
