# PMSM Optimizer Sub-flows Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Create standalone sub-flow script modules (`scripts/`) implementing all 5 PMSM sub-flows, update `references/pmsm-methods.md` and `SKILL.md` to reference them, and verify each against the AEDT template.

**Architecture:** Template-formula approach — one-time GUI winding setup, then fully automated via design variable changes. Sub-flow A (Magnetostatic) uses `assign_current()` with 9× turns factor. Sub-flows B/C/D (Transient) modify design variables `Imax`, `Speed_rpm`, `Thet_deg`. Sub-flow E (pure math) uses interpolation on Sub-flow A results.

**Tech Stack:** PyAEDT 0.25.1 (gRPC), Ansys Maxwell 2D 2023.2, Python 3.13, Windows 11.

## Global Constraints

- All scripts use `sys.stdout.reconfigure(encoding='utf-8')` as first line
- All scripts open with `non_graphical=True, new_desktop=True, close_on_exit=True` unless specified otherwise
- Template file: `references/Prius_2D_Practice.aedt` (absolute path relative to skill directory)
- Design variables: `Imax` (current amplitude), `Speed_rpm` (speed), `Thet_deg` (current angle), `Omega_rad` (electrical angular velocity), `PolePairs` (pole pairs)
- Winding formula (Transient): `Imax*sin(Omega_rad*time + Thet)` (Phase A), `-2*pi/3` (Phase B), `+2*pi/3` (Phase C)
- Magnetostatic assign_current must multiply by turns factor: **9 × Ia** per coil object
- PyAEDT gRPC limitations: `ChangeProperty`, `SetPropertyValue`, `DeleteBoundary`, `GetPropertyValue` all fail — use `m2d[...]=value` for variables, `assign_current()` for current excitations

---

### Task 1: Project structure — add sub-flow script modules

**Files:**
- Create: `C:\Users\Wesley\.claude\skills\pmsm-optimizer\scripts\__init__.py`
- Create: `C:\Users\Wesley\.claude\skills\pmsm-optimizer\scripts\subflow_a_ldlq.py`
- Create: `C:\Users\Wesley\.claude\skills\pmsm-optimizer\scripts\subflow_b_bemf.py`
- Create: `C:\Users\Wesley\.claude\skills\pmsm-optimizer\scripts\subflow_c_torque.py`
- Create: `C:\Users\Wesley\.claude\skills\pmsm-optimizer\scripts\subflow_d_efficiency_map.py`
- Create: `C:\Users\Wesley\.claude\skills\pmsm-optimizer\scripts\subflow_e_external.py`
- Overwrite: `C:\Users\Wesley\.claude\skills\pmsm-optimizer\references\pmsm-methods.md`
- Modify: `C:\Users\Wesley\.claude\skills\pmsm-optimizer\SKILL.md`

**Interfaces:**
- Consumes: Template file `references/Prius_2D_Practice.aedt` with 2 designs (4_MS2, 5_TR)
- Produces: Runnable scripts for each sub-flow

- [ ] **Step 1: Create `scripts/__init__.py`**

```python
# PMSM Optimizer - Sub-flow scripts package
```

- [ ] **Step 2: Create `scripts/subflow_a_ldlq.py`**

```python
"""Sub-flow A: Ld/Lq MAP + 主磁链 Φ
求解器: MagnetostaticXY (4_Partial_motor_MS2)
方法: 电流角扫描法, assign_current() 设置激励 (含 9 匝因子)

用法:
    python -c "from scripts.subflow_a_ldlq import run; run(rated_current=10, max_current=3, current_steps=6, angle_steps=7)"
"""

import sys
sys.stdout.reconfigure(encoding='utf-8')
from pathlib import Path
import numpy as np
from ansys.aedt.core import Maxwell2d

ROOT = Path(__file__).resolve().parent.parent
TEMPLATE = str(ROOT / 'references' / 'Prius_2D_Practice.aedt')

TURNS = 9  # Matrix NumberOfTurns per coil object


def id_iq_to_abc(Id, Iq, theta=0):
    """Inverse Park + Clarke: dq → abc (theta in radians)"""
    # θ=0: d-axis aligned with Phase A
    Ialpha = Id * np.cos(theta) - Iq * np.sin(theta)
    Ibeta  = Id * np.sin(theta) + Iq * np.cos(theta)
    Ia = Ialpha
    Ib = -0.5 * Ialpha + np.sqrt(3) / 2 * Ibeta
    Ic = -0.5 * Ialpha - np.sqrt(3) / 2 * Ibeta
    return Ia, Ib, Ic


def park_abc_to_dq(psi_a, psi_b, psi_c, theta=0):
    """Park transform: abc → dq (theta in radians)"""
    c = np.cos(theta)
    s = np.sin(theta)
    psi_d = 2 / 3 * (psi_a * c + psi_b * np.cos(theta - 2 * np.pi / 3) + psi_c * np.cos(theta + 2 * np.pi / 3))
    psi_q = -2 / 3 * (psi_a * s + psi_b * np.sin(theta - 2 * np.pi / 3) + psi_c * np.sin(theta + 2 * np.pi / 3))
    return psi_d, psi_q


def run(rated_current=10, max_current=3.0, current_steps=6, angle_steps=7):
    I_rated = rated_current
    I_max = max_current * I_rated
    currents = np.linspace(I_max / current_steps, I_max, current_steps)
    angles_deg = np.linspace(0, 90, angle_steps)

    m2d = Maxwell2d(
        project=TEMPLATE, design='4_Partial_motor_MS2',
        solution_type='MagnetostaticXY',
        non_graphical=False, new_desktop=False, close_on_exit=False
    )

    results = []
    for Ia_mag in currents:
        for theta_deg in angles_deg:
            theta = np.radians(theta_deg)
            Id = Ia_mag * np.sin(theta)
            Iq = Ia_mag * np.cos(theta)

            # Transform Id/Iq → three-phase currents
            Ia, Ib, Ic = id_iq_to_abc(Id, Iq, theta=0)

            # Assign currents with turns factor (9 turns per coil)
            m2d.assign_current(assignment=['PhaseA1', 'PhaseA2'], current=TURNS * Ia, units='A')
            m2d.assign_current(assignment=['PhaseB1', 'PhaseB2'], current=TURNS * Ib, units='A')
            m2d.assign_current(assignment=['PhaseC1', 'PhaseC2'], current=TURNS * Ic, units='A')

            # Solve
            m2d.analyze('Setup1')

            # Extract flux linkages
            data = m2d.post.get_solution_data(
                expressions=[
                    'FluxLinkage(PhaseA1)', 'FluxLinkage(PhaseA2)',
                    'FluxLinkage(PhaseB1)', 'FluxLinkage(PhaseB2)',
                    'FluxLinkage(PhaseC1)', 'FluxLinkage(PhaseC2)',
                ],
                variations=m2d.post.get_solution_data_variation()
            )

            # Total phase flux linkage = turns × sum of coil flux linkages
            psi_a = TURNS * (data.data('FluxLinkage(PhaseA1)')[0] + data.data('FluxLinkage(PhaseA2)')[0])
            psi_b = TURNS * (data.data('FluxLinkage(PhaseB1)')[0] + data.data('FluxLinkage(PhaseB2)')[0])
            psi_c = TURNS * (data.data('FluxLinkage(PhaseC1)')[0] + data.data('FluxLinkage(PhaseC2)')[0])

            # Park transform → dq flux linkages
            psi_d, psi_q = park_abc_to_dq(psi_a, psi_b, psi_c, theta=0)

            # Inductances (avoid division by zero)
            Ld = psi_d / Id if abs(Id) > 1e-6 else 0
            Lq = psi_q / Iq if abs(Iq) > 1e-6 else 0

            results.append({
                'Id': Id, 'Iq': Iq,
                'Psi_d': psi_d, 'Psi_q': psi_q,
                'Ld': Ld, 'Lq': Lq,
            })

            print(f'  Id={Id:+7.2f}  Iq={Iq:+7.2f}  Ld={Ld:.5f}  Lq={Lq:.5f}')

    # 主磁链 Φ: Id=0, Iq≈0 时的 Ψq（永磁体贡献）
    near_zero = [r for r in results if abs(r['Id']) < 1e-6 and abs(r['Iq']) < 1e-6]
    phi = near_zero[0]['Psi_q'] if near_zero else results[0]['Psi_q']

    # 输出表格
    print(f'\n{"=" * 60}')
    print(f'主磁链 Φ = {phi:.6f} Wb')
    print(f'{"=" * 60}')
    print(f'{"Id\\Iq":>8s}', end='')
    for a in angles_deg:
        print(f'  {a:7.0f}°    ', end='')
    print()
    for i, Ia_mag in enumerate(currents):
        print(f'{Ia_mag:8.2f}', end='')
        for j in range(len(angles_deg)):
            idx = i * len(angles_deg) + j
            print(f'  Ld={results[idx]["Ld"]:.4f}', end='')
        print()
        print(f'{"":>8s}', end='')
        for j in range(len(angles_deg)):
            idx = i * len(angles_deg) + j
            print(f'  Lq={results[idx]["Lq"]:.4f}', end='')
        print()

    m2d.close_project()
    return results, phi


if __name__ == '__main__':
    run()
```

- [ ] **Step 3: Create `scripts/subflow_b_bemf.py`**

```python
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
```

- [ ] **Step 4: Create `scripts/subflow_c_torque.py`**

```python
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
TEMPLATE = str(ROOT / 'references' / 'Prius_2D_Practice.aedt')


def _mtpa_scan(m2d, rated_speed, pole_pairs):
    """MTPA 角度扫描: 候选角各跑 1 个电周期, 选扭矩最大者"""
    candidate_angles = [0, -15, -25, -35, -45, -60]
    freq = rated_speed / 60 * pole_pairs

    # 粗网格: 100 步/周期
    setup = m2d.setups[0]
    setup.props['StopTime'] = f'{1/freq}s'
    setup.props['TimeStep'] = f'{1/(freq*100)}s'
    setup.update()

    best_T = -1e9
    best_angle = candidate_angles[0]

    for angle in candidate_angles:
        m2d['Thet_deg'] = f'{angle}°'
        m2d.analyze('Setup1')

        data = m2d.post.get_solution_data(
            expressions=['Moving1.Torque'],
            variations=m2d.post.get_solution_data_variation()
        )
        torque = np.array(data.data('Moving1.Torque'))
        T_avg = np.mean(torque[-50:])  # 最后半周期平均
        print(f'    θ={angle:3d}° → T={T_avg:.2f} Nm')

        if T_avg > best_T:
            best_T = T_avg
            best_angle = angle

    print(f'  MTPA 最优角: θ={best_angle}° (T={best_T:.2f} Nm)')
    return best_angle


def run(rated_current=250, current_angle=None, rated_speed=3000, elec_periods=3):
    from ansys.aedt.core import Maxwell2d

    m2d = Maxwell2d(
        project=TEMPLATE, design='5_Partial_motor_TR',
        solution_type='TransientXY',
        non_graphical=False, new_desktop=False, close_on_exit=False
    )

    # 设置负载和转速
    m2d['Imax'] = f'{rated_current}A'
    m2d['Speed_rpm'] = f'{rated_speed}rpm'
    pole_pairs = float(m2d['PolePairs'])

    # MTPA 自动搜索（未指定电流角时）
    auto_mtpa = current_angle is None
    if auto_mtpa:
        print('  未指定电流角, 正在搜索 MTPA 最优角...')
        current_angle = _mtpa_scan(m2d, rated_speed, pole_pairs)

    m2d['Thet_deg'] = f'{current_angle}°'

    # 完整仿真
    freq = rated_speed / 60 * pole_pairs
    stop_time = elec_periods / freq
    time_step = 1 / (freq * 200)  # 200 steps/period

    setup = m2d.setups[0]
    setup.props['StopTime'] = f'{stop_time}s'
    setup.props['TimeStep'] = f'{time_step}s'
    setup.update()

    m2d.analyze('Setup1')

    # 提取扭矩波形
    data = m2d.post.get_solution_data(
        expressions=['Moving1.Torque'],
        variations=m2d.post.get_solution_data_variation()
    )
    torque = np.array(data.data('Moving1.Torque'))

    # 取最后一个周期的稳态数据
    steps_per_period = 200
    steady_torque = torque[-steps_per_period:]

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
```

- [ ] **Step 5: Create `scripts/subflow_d_efficiency_map.py`**

```python
"""Sub-flow D: 全域工作特性效率 MAP
求解器: TransientXY (5_Partial_motor_TR)
方法: (转速 × 扭矩) 二维网格扫描
      每个工作点的 Id/Iq 由 **MTPA 算法** 从目标扭矩生成,
      而非经验系数估算.

核心逻辑:
  MTPA 算法 (_mtpa_for_torque):
    - 沿 MTPA 轨迹 Id = A - √(A²+Iq²) 对 Iq 二分搜索匹配 T_target
    - 其中 A = Φ / (2·(Lq-Ld))
  FEA: 以 MTPA 给出的 (Imax, β) 求解, 取 Moving1.Torque 单值
  解析: 基于同一 (Id, Iq) 计算 Vs/M/PF/损耗

gRPC 限制: 加载 Transient 只返回单值时点, 故电参数用解析而非波形.
"""

import sys
sys.stdout.reconfigure(encoding='utf-8')
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
TEMPLATE = str(ROOT / 'references' / 'Prius_2D_Practice.aedt')

DEFAULT_Ld = 0.00035   # H
DEFAULT_Lq = 0.00080   # H
DEFAULT_Phi = 0.07     # Wb
DEFAULT_Rs = 0.0       # Ω


def _mtpa_for_torque(T_target, pole_pairs, Ld, Lq, Phi, imax):
    """MTPA: T_target → (Id, Iq, Is, β) 使电流最小"""
    if T_target <= 0:
        return 0.0, 0.0, 0.0, 0.0
    delta_L = Lq - Ld
    A = Phi / (2 * delta_L) if delta_L > 1e-12 else 1e12

    def torque_at_Iq(Iq):
        Id = A - np.sqrt(A * A + Iq * Iq) if delta_L > 1e-12 else 0.0
        return 1.5 * pole_pairs * (Phi * Iq - delta_L * Id * Iq)

    if T_target >= torque_at_Iq(imax):
        Iq = imax
    else:
        lo, hi = 0.0, imax
        for _ in range(60):
            mid = (lo + hi) / 2.0
            if torque_at_Iq(mid) > T_target:
                hi = mid
            else:
                lo = mid
        Iq = (lo + hi) / 2.0
    Id = A - np.sqrt(A * A + Iq * Iq) if delta_L > 1e-12 else 0.0
    Is = np.sqrt(Id * Id + Iq * Iq)
    beta = np.arctan2(-Id, Iq)
    return Id, Iq, Is, beta


def run(speed_min=1000, speed_max=8000, speed_steps=5, torque_steps=5,
        vdc=300, imax=250,
        Ld=DEFAULT_Ld, Lq=DEFAULT_Lq, Phi=DEFAULT_Phi, Rs=DEFAULT_Rs):
    # ... (创建 m2d, 网格循环, 每点调 _mtpa_for_torque + FEA + 解析计算 + CSV 输出)
```

- [ ] **Step 6: Create `scripts/subflow_e_external.py`**

```python
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
```

- [ ] **Step 7: Overwrite `references/pmsm-methods.md` with formulas + module refs**

Replace the entire file (keep the header). Each sub-flow gets: method summary, key formulas, and module reference.

```markdown
# PMSM Optimizer — 子流程计算模板

本文档包含 5 个子流程的方法说明和关键公式。完整代码见 `scripts/` 目录。

---

## Sub-flow A：Ld/Lq MAP + 主磁链 Φ

**求解器**: MagnetostaticXY (`4_Partial_motor_MS2`)

**方法**: 电流角扫描法。在每个 (Id, Iq) 工作点用 `assign_current()` 设置三相电流（含 9 匝因子），求解后提取磁链做 Park 变换。

**关键公式**:
```
Id = Ia_mag * sin(θ),  Iq = Ia_mag * cos(θ)          θ = 0~90°
Ia = Id,  Ib = -0.5×Id + 0.866×Iq,  Ic = -0.5×Id - 0.866×Iq
assign_current 电流值 = 9 × Ia/Ib/Ic                 (9 匝/线圈)
Ψa = 9 × (Ψ_A1 + Ψ_A2),  Ψb = 9 × (Ψ_B1 + Ψ_B2),  Ψc = 9 × (Ψ_C1 + Ψ_C2)
Ψd = Park(Ψa,Ψb,Ψc),  Ψq = Park(Ψa,Ψb,Ψc)           (θ=0°)
Ld = Ψd/Id,  Lq = Ψq/Iq
```

**模块**: `scripts/subflow_a_ldlq.py` — `run(rated_current, max_current, current_steps, angle_steps)`

---

## Sub-flow B：空载反电势

**求解器**: TransientXY (`5_Partial_motor_TR`)

**方法**: Imax=0 空载，额定转速旋转，提取三相电压波形做 FFT 谐波分析。

**关键公式**:
```
freq = Speed_rpm / 60 × PolePairs
StopTime = elec_periods / freq
Ke = V_fundamental / (Speed_rpm/1000)
Ke_line = Ke × √3
THD = √(Σ|Hk|²) / |H1|          k ≥ 2
```

**模块**: `scripts/subflow_b_bemf.py` — `run(rated_speed, elec_periods, time_steps_per_cycle)`

---

## Sub-flow C：额定点扭矩

**求解器**: TransientXY (`5_Partial_motor_TR`)

**方法**: 额定电流、额定转速。电流角自动 MTPA（候选 [0,-15,-25,-35,-45,-60]° 各 1 周期扫描）或用户指定。

**关键公式**:
```
Id = Imax × sin(Thet),  Iq = Imax × cos(Thet)
T = 1.5 × PolePairs × (Φ×Iq + (Ld-Lq)×Id×Iq)
T_ripple = (T_max - T_min) / T_avg × 100%
```

**模块**: `scripts/subflow_c_torque.py` — `run(rated_current, current_angle, rated_speed, elec_periods)`

---

## Sub-flow D：全域工作特性 MAP

**求解器**: 批量 TransientXY (`5_Partial_motor_TR`)

**方法**: (转速 × 扭矩) 二维网格扫描，每点一个电周期。估算 Id/Iq 工作点提取 4 张 MAP。

**关键公式**:
```
P_cu = 3 × Is² × Rs
PF ≈ sign_correlation(v, i)
M = Vs / (Vdc/√3)
```

**模块**: `scripts/subflow_d_efficiency_map.py` — `run(speed_min, speed_max, speed_steps, torque_steps, vdc, imax)`

---

## Sub-flow E：外特性计算

**方法**: 纯数学计算（不跑 FEA）。电压圆 Vs ≤ Vdc/√3 (M_max=1.0), 电流圆 Is ≤ Imax。

**关键公式**:
```
Vd = -ω × Lq × Iq,  Vq = ω × (Ld × Id + Φ)
Vs = √(Vd² + Vq²) ≤ Vdc/√3
Is = √(Id² + Iq²) ≤ Imax
T = 1.5 × PolePairs × (Φ×Iq + (Ld-Lq)×Id×Iq)
```

**模块**: `scripts/subflow_e_external.py` — `run(vdc, imax, speed_max, rs, speed_points_n, pole_pairs, ld_lq_data, phi)`
```

- [ ] **Step 8: Update `SKILL.md`**

Three targeted edits:

1. Line ~47: Change template name `pmsm_template` → `Prius_2D_Practice`
2. Line ~88-100: Replace parameter table with actual template variables
3. Line ~146-150: Update sub-flow code reference to `scripts/` modules

```python
# Edit 1: template name (in-memory replacement)
with open('SKILL.md', 'r', encoding='utf-8') as f:
    content = f.read()
content = content.replace('pmsm_template', 'Prius_2D_Practice')
with open('SKILL.md', 'w', encoding='utf-8') as f:
    f.write(content)
```

Edit 2: Replace the parameter table in "第二步" with:
```markdown
| 参数 | 变量名 | 修改方式 |
|------|--------|---------|
| 相电流幅值 | `Imax` | `m2d['Imax'] = '250A'` |
| 电流角 | `Thet_deg` | `m2d['Thet_deg'] = '-30°'` |
| 机械转速 | `Speed_rpm` | `m2d['Speed_rpm'] = '3000rpm'` |
| 极对数 | `PolePairs` | 固定值 (Poles/2=4) |
```

Edit 3: In "子流程执行规范" section, replace:
```
所有子流程的代码模板位于 `references/pmsm-methods.md`。
```
with:
```
所有子流程的代码位于 `scripts/` 目录，方法说明见 `references/pmsm-methods.md`。
AI 根据用户指令选择模块调用。详见设计文档 `docs/superpowers/specs/2026-07-12-pmsm-optimizer-design.md`。
```

- [ ] **Step 9: Create results directory**

```python
python -c "import os; os.makedirs('results', exist_ok=True)"
```

- [ ] **Step 10: Commit**

```bash
git add scripts/ SKILL.md references/pmsm-methods.md results/
git commit -m "feat: implement all 5 PMSM sub-flows
- Create scripts/ modules with complete run() functions
- Sub-flow A: Ld/Lq MAP via assign_current() with 9x turns
- Sub-flow B: Back EMF with FFT + THD + Ke
- Sub-flow C: Rated torque with auto MTPA scan
- Sub-flow D: Efficiency MAP grid scan (4 CSVs)
- Sub-flow E: External characteristic (math-based)
- Update SKILL.md with template variable table
- Update pmsm-methods.md with formulas and module refs"
```

---

### Verification

After all steps, run each sub-flow and validate output:

1. **Sub-flow A**: `python -c "from scripts.subflow_a_ldlq import run; run(rated_current=250, current_steps=3, angle_steps=3)"` — Verify Ld < Lq (IPM typical), saturation with increasing current.

2. **Sub-flow B**: `python -c "from scripts.subflow_b_bemf import run; run(rated_speed=3000)"` — Verify THD < 15%, Ke ~50-80 V/krpm for Prius-class motor.

3. **Sub-flow C**: `python -c "from scripts.subflow_c_torque import run; run(rated_current=250, rated_speed=3000)"` — Verify torque > 0, auto MTPA finds optimal angle.

4. **Sub-flow D**: `python -c "from scripts.subflow_d_efficiency_map import run; run(speed_steps=3, torque_steps=3)"` — Verify 4 CSVs saved with smooth variation.

5. **Sub-flow E**: `python -c "from scripts.subflow_e_external import run; run(vdc=300, imax=250, speed_max=12000)"` — Verify torque decreases at high speed (flux weakening).
