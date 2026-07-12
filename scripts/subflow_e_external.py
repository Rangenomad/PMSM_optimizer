"""Sub-flow E: 外特性计算 (T-n 曲线)
方法: 纯数学计算 (不跑 FEA)
原理: 电压圆 Vs ≤ Vdc/√3, 电流圆 Is ≤ Imax 约束下, 基于 Ld/Lq MAP + Φ 计算

用法:
    python -c "from scripts.subflow_e_external import run; run(vdc=300, imax=250, speed_max=12000)"
"""

import sys
sys.stdout.reconfigure(encoding='utf-8')
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent


def run(vdc=300, imax=250, speed_max=12000, rs=0.0, speed_points_n=20, pole_pairs=4,
        ld_lq_data=None, phi=None):
    """
    Parameters
    ----------
    vdc : float
        直流母线电压 (V)
    imax : float
        最大相电流幅值 (A)
    speed_max : float
        最高机械转速 (rpm) — 外特性曲线的截止转速
    rs : float
        相电阻 (Ω)
    speed_points_n : int
        转速分点数
    pole_pairs : int
        极对数
    ld_lq_data : list of dict, optional
        Sub-flow A 结果: [{'Id':..., 'Iq':..., 'Ld':..., 'Lq':...}, ...]
    phi : float, optional
        主磁链 (Wb)
    """
    # 电压限幅 (SVPWM 线性调制区)
    Vmax = vdc / np.sqrt(3)  # M_max = 1.0

    # 如果没有 Sub-flow A 数据, 使用近似值
    if ld_lq_data is None:
        Ld = 0.00035  # H
        Lq = 0.0008   # H
        Phi = phi if phi else 0.12  # Wb
    else:
        # 从 Sub-flow A 结果取平均 Ld/Lq (不饱和区)
        # 进阶: 可构建 2D 插值模型 (RegularGridInterpolator) 获取饱和特性
        ld_vals = np.array([d['Ld'] for d in ld_lq_data])
        lq_vals = np.array([d['Lq'] for d in ld_lq_data])
        Ld = np.mean(ld_vals[ld_vals > 0]) if np.any(ld_vals > 0) else 0.00035
        Lq = np.mean(lq_vals[lq_vals > 0]) if np.any(lq_vals > 0) else 0.0008
        Phi = phi if phi else 0.12

    print(f'Vdc={vdc}V, Imax={imax}A, Rs={rs}Ω')
    print(f'Ld={Ld*1000:.4f}mH, Lq={Lq*1000:.4f}mH, Phi={Phi:.4f}Wb')

    # 转速范围 (从 1 rpm 避免除零, 几何分布: 低速密高速疏)
    speeds = np.geomspace(max(1, speed_max / speed_points_n), speed_max, speed_points_n)

    results = []
    for n in speeds:
        omega_e = 2 * np.pi * n / 60 * pole_pairs  # 电角速度 (rad/s)

        # MTPA 工作点搜索 (恒扭矩区)
        # 电压方程: Vd = -omega_e * Lq * Iq
        #           Vq = omega_e * (Ld * Id + Phi)
        #           Vs = sqrt(Vd^2 + Vq^2)
        # 电流约束: Id^2 + Iq^2 <= Imax^2
        # 电压约束: Vs <= Vmax

        # 简化 MTPA: 在电流圆上找最大化 T = 1.5 * pole_pairs * (Phi*Iq + (Ld-Lq)*Id*Iq) 的点
        n_search = 500
        best_T = 0
        best_Id = 0
        best_Iq = 0
        region = 'MTPA'

        for angle in np.linspace(np.pi/2, 0, n_search):  # Id 从 -Imax 到 0
            Id_candidate = -imax * np.sin(angle)
            Iq_candidate = imax * np.cos(angle)

            # 电压方程 (忽略电阻压降)
            Vd = -omega_e * Lq * Iq_candidate
            Vq = omega_e * (Ld * Id_candidate + Phi)
            Vs = np.sqrt(Vd**2 + Vq**2)
            M = Vs / Vmax  # 调制比 (M > 1.0 表示进入过调制/需弱磁)

            if M > 1.0:
                continue  # 超出电压限, 需要弱磁

            # 扭矩方程 (IPM: 永磁转矩 + 磁阻转矩)
            T = 1.5 * pole_pairs * (Phi * Iq_candidate + (Ld - Lq) * Id_candidate * Iq_candidate)

            if T > best_T:
                best_T = T
                best_Id = Id_candidate
                best_Iq = Iq_candidate

        # 弱磁区: 在电压圆和电流圆交点上找最大扭矩
        if best_T == 0:
            region = '弱磁'
            for Id_candidate in np.linspace(-imax, 0, n_search):
                Iq_candidate = np.sqrt(imax**2 - Id_candidate**2)

                Vd = -omega_e * Lq * Iq_candidate
                Vq = omega_e * (Ld * Id_candidate + Phi)
                Vs = np.sqrt(Vd**2 + Vq**2)

                if Vs <= Vmax:
                    T = 1.5 * pole_pairs * (Phi * Iq_candidate + (Ld - Lq) * Id_candidate * Iq_candidate)
                    if T > best_T:
                        best_T = T
                        best_Id = Id_candidate
                        best_Iq = Iq_candidate

        # 输出功率
        P = best_T * n * 2 * np.pi / 60 / 1000  # kW

        # 调制比
        Vd_op = -omega_e * Lq * best_Iq
        Vq_op = omega_e * (Ld * best_Id + Phi)
        Vs_op = np.sqrt(Vd_op**2 + Vq_op**2)
        M = Vs_op / Vmax  # M_max = 1.0 (SVPWM 线性调制)

        results.append({
            'speed': n,
            'torque': best_T,
            'power': P,
            'Id': best_Id,
            'Iq': best_Iq,
            'modulation': M,
            'region': region,
        })

    # 输出
    print(f'\n{"=" * 70}')
    print(f'=== 外特性 (T-n 曲线) ===')
    print(f'Vdc = {vdc} V, Imax = {imax} A, Rs = {rs} Ω')
    print(f'{"=" * 70}')
    print(f'{"转速(rpm)":>10s} {"扭矩(Nm)":>8s} {"功率(kW)":>8s} {"Id(A)":>8s} {"Iq(A)":>8s} {"调制比":>8s} {"区域":>8s}')
    for r in results:
        print(f'{r["speed"]:10.0f} {r["torque"]:8.1f} {r["power"]:8.2f} {r["Id"]:8.2f} {r["Iq"]:8.2f} {r["modulation"]:8.3f} {r["region"]:>8s}')

    # 转折点
    for i, r in enumerate(results):
        if r['region'] == '弱磁' or r['modulation'] > 0.95:
            print(f'\n转折点转速: {r["speed"]:.0f} rpm')
            break

    # 峰值功率
    peak = max(results, key=lambda r: r['power'])
    print(f'峰值功率: {peak["power"]:.2f} kW @ {peak["speed"]:.0f} rpm')

    return results


if __name__ == '__main__':
    run()
