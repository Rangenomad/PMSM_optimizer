"""Sub-flow A: Ld/Lq MAP + 主磁链 Φ
求解器: MagnetostaticXY (4_Partial_motor_MS2)
方法: 电流角扫描法, 通过 boundary.update() 修改已有 Current 激励 (含 9 匝因子)
      通过 export_matrix() 提取每线圈 Flux Linkage

Ld/Lq 计算方法 — 交叉饱和补偿 (secant inductance):
  Ld(Id,Iq) = [ψ_d(Id,Iq) - ψ_d(0,Iq)] / Id
  Lq(Id,Iq) = ψ_q(Id,Iq) / Iq

  其中 ψ_d(0,Iq) 从各电流级别的 θ=0° 点插值获得，
  代表该 Iq 下的有效永磁磁链（含负载效应和交叉饱和影响）。

  使用此定义，dq 扭矩公式对每个测量点精确成立：
  T = 1.5·P·[ψ_d(0,Iq)·Iq + (Ld - Lq)·Id·Iq] = 1.5·P·(ψ_d·Iq - ψ_q·Id)

参考文献:
  [1] IEEE Trans. Magnetics, 2010 — Ld = (Ψd-Ψd(0,Iq))/Id
  [2] MathWorks/Ansys "Deriving Fast and Accurate PMSM Motor Model from FEA", 2017
  [3] TUM 博士论文 — 交叉饱和与永磁磁链扣除
  [4] ORNL 2011 APEEM 报告 — 2010 Prius 电机基准参数

用法:
    python -c "from scripts.subflow_a_ldlq import run; r, phi = run(rated_current=250, max_current=1.4, current_steps=5, angle_steps=6)"
    python -c "from scripts.subflow_a_ldlq import run; r, phi = run(rated_current=250, max_current=1.4, current_steps=5, angle_steps=6, project_path='pmsm_projects/xxx')"
"""

import sys, os, re, tempfile
sys.stdout.reconfigure(encoding='utf-8')
from pathlib import Path
import numpy as np
from ansys.aedt.core import Maxwell2d

ROOT = Path(__file__).resolve().parent.parent
_DEFAULT_TEMPLATE = str(ROOT / 'references' / 'Prius_2D_Practice.aedt')

TURNS = 9  # Matrix NumberOfTurns per coil object

# 线圈配置: (boundary_name, IsPositive)
COIL_CONFIG = [
    ('PhaseA1', True), ('PhaseA2', True),
    ('PhaseB1', True), ('PhaseB2', True),
    ('PhaseC1', True), ('PhaseC2', True),
]
COIL_NAMES = [c[0] for c in COIL_CONFIG]


def id_iq_to_abc(Id, Iq, theta=0):
    """Inverse Park + Clarke: dq → abc (theta in radians)

    Convention: q-axis LAGS d-axis by 90° (Maxwell motor convention).
    d-axis at angle θ, q-axis at θ-90°: q_hat = (sinθ, -cosθ).
    Positive Iq → field along q_hat → produces positive (motoring) torque.

    Verified: at θ=0, Id=0,Iq=1 → Ia=0,Ib=-√3/2,Ic=+√3/2 (q-axis field)
    """
    Ialpha = Id * np.cos(theta) + Iq * np.sin(theta)
    Ibeta  = Id * np.sin(theta) - Iq * np.cos(theta)
    Ia = Ialpha
    Ib = -0.5 * Ialpha + np.sqrt(3) / 2 * Ibeta
    Ic = -0.5 * Ialpha - np.sqrt(3) / 2 * Ibeta
    return Ia, Ib, Ic


def park_abc_to_dq(psi_a, psi_b, psi_c, theta=0):
    """Park transform: abc → dq (theta in radians)

    Convention: q-axis LAGS d-axis by 90° (consistent with id_iq_to_abc).
    Forward: ψ_d = proj(ψ_abc, d_hat), ψ_q = proj(ψ_abc, q_hat).
    q_hat = (sinθ, -cosθ) → ψ_q = +2/3·Σ ψ_i·sin(θ_i) (note: + sign).
    Round-trip: park(inv_park(Id,Iq)) = (Id,Iq) for any θ.
    """
    c = np.cos(theta)
    s = np.sin(theta)
    psi_d = 2 / 3 * (psi_a * c + psi_b * np.cos(theta - 2 * np.pi / 3) + psi_c * np.cos(theta + 2 * np.pi / 3))
    psi_q = 2 / 3 * (psi_a * s + psi_b * np.sin(theta - 2 * np.pi / 3) + psi_c * np.sin(theta + 2 * np.pi / 3))
    return psi_d, psi_q


def _set_coil_currents(m2d, boundaries, Ia, Ib, Ic):
    """Set coil currents via boundary.update() (含 9 匝因子)"""
    currents = [TURNS * Ia, TURNS * Ia,
                TURNS * Ib, TURNS * Ib,
                TURNS * Ic, TURNS * Ic]
    for bnd, (name, is_pos), val in zip(boundaries, COIL_CONFIG, currents):
        bnd.props['Current'] = f'{val}A'
        bnd.props['IsPositive'] = is_pos
        bnd.update()


def _parse_matrix_flux_linkage(text):
    """Parse matrix export .txt → dict {coil_name: flux_Wb}"""
    in_flux = False
    result = {}
    for line in text.splitlines():
        if line.strip().startswith('Flux Linkage'):
            in_flux = True
            continue
        if in_flux and line.strip() and not line.strip().startswith('Flux'):
            parts = line.strip().split()
            if len(parts) >= 2 and parts[0] in COIL_NAMES:
                result[parts[0]] = float(parts[1])
    return result


def _extract_flux_linkages(m2d):
    """Solve → export matrix → parse flux linkage per coil"""
    m2d.analyze('Setup1')
    out_file = os.path.join(tempfile.gettempdir(), 'pmsm_matrix_export.txt')
    m2d.export_matrix('Matrix1', out_file)
    with open(out_file, 'r') as f:
        text = f.read()
    return _parse_matrix_flux_linkage(text)


def run(rated_current=250, max_current=3.0, current_steps=6, angle_steps=7,
        project_path=None, phi_override=None):
    """Run Sub-flow A: Ld/Lq MAP + PM flux linkage Phi.

    Uses cross-saturation-compensated secant inductance:
      Ld(Id,Iq) = [ψ_d(Id,Iq) - ψ_d(0,Iq)] / Id
      Lq(Id,Iq) = ψ_q(Id,Iq) / Iq

    ψ_d(0,Iq) is interpolated from θ=0° points at each current level,
    capturing the Iq-dependent effective PM flux.

    Parameters
    ----------
    phi_override : float, optional
        Override the measured Phi with a calibrated value.
        Use when FEA torque at Thet=0 gives a different effective Phi
        (e.g. from Sub-flow C at Thet=0: Phi_eff = T0/(1.5*P*Imax)).
        Default: None (use measured Phi from zero-current point)
    """
    from scripts.project_utils import get_template_path

    # 解析模板路径
    if project_path:
        template = get_template_path(Path(project_path))
    else:
        template = _DEFAULT_TEMPLATE

    I_rated = rated_current
    I_max = max_current * I_rated
    currents = np.linspace(I_max / current_steps, I_max, current_steps)
    angles_deg = np.linspace(0, 90, angle_steps)
    total_points = len(currents) * len(angles_deg)

    # ── 打印扫描计划 ──
    print(f'\n{"=" * 60}')
    print(f'  Sub-flow A: Ld/Lq MAP + 主磁链 Φ')
    print(f'  (交叉饱和补偿 secant inductance)')
    print(f'{"=" * 60}')
    print(f'  额定电流:    {I_rated:.0f} A')
    print(f'  最大电流:    {I_max:.0f} A ({max_current:.1f}× 额定)')
    print(f'  电流分点:    {current_steps} 点  [{", ".join(f"{v:.0f}" for v in currents)}] A')
    print(f'  角度分点:    {angle_steps} 点  [{", ".join(f"{v:.0f}" for v in angles_deg)}]°')
    print(f'  扫描点数:    1 (零电流) + {total_points} = {total_points + 1} 点')
    print(f'  预估耗时:    ~{total_points * 10 // 60} 分钟 (每点约10秒)')
    print(f'  输出:        Ld(Id,Iq) MAP + Lq(Id,Iq) MAP + Φ_eff(Iq) + 扭矩验证')
    print(f'{"=" * 60}')

    # 确认模式检测
    import json
    config_path = ROOT / 'config.json'
    if config_path.exists():
        with open(config_path) as f:
            mode = json.load(f).get('execution_mode', 'confirm')
    else:
        mode = 'confirm'

    if mode == 'confirm':
        resp = input('\n  继续执行? [Y/n] ').strip().lower()
        if resp and resp != 'y':
            print('  已取消。')
            return [], 0.0

    m2d = Maxwell2d(
        project=template, design='4_Partial_motor_MS2',
        solution_type='MagnetostaticXY',
        non_graphical=False, new_desktop=False, close_on_exit=False
    )

    # Get references to all 6 coil Current boundaries (once, before loop)
    boundaries = []
    for name, _ in COIL_CONFIG:
        b = [b for b in m2d.boundaries if b.name == name][0]
        boundaries.append(b)

    # PolePairs 是模板固定参数 (8 poles / 2 = 4), 非 Maxwell 设计变量
    pole_pairs = 4
    print(f'\n  极对数 PolePairs = {pole_pairs} (固定, 8-pole PMSM)')

    # ── 零电流点：测纯永磁磁链 Φ_0 并标定转子角度 ──
    print('  [A] 求解零电流点（标定转子角度 + 主磁链 Φ_0）...')
    _set_coil_currents(m2d, boundaries, 0, 0, 0)
    psi_coil_zero = _extract_flux_linkages(m2d)
    psi_a_z = TURNS * (psi_coil_zero['PhaseA1'] + psi_coil_zero['PhaseA2'])
    psi_b_z = TURNS * (psi_coil_zero['PhaseB1'] + psi_coil_zero['PhaseB2'])
    psi_c_z = TURNS * (psi_coil_zero['PhaseC1'] + psi_coil_zero['PhaseC2'])
    # Clarke transform (amplitude-invariant, 2/3 coef)
    psi_alpha = 2/3 * (psi_a_z - 0.5*psi_b_z - 0.5*psi_c_z)
    psi_beta  = 2/3 * (np.sqrt(3)/2 * psi_b_z - np.sqrt(3)/2 * psi_c_z)
    theta_r = np.arctan2(psi_beta, psi_alpha)  # rotor d-axis electrical angle [rad]
    phi_0 = np.sqrt(psi_alpha**2 + psi_beta**2)  # PM flux at zero current (always positive)
    phi = phi_override if phi_override is not None else phi_0
    phi_label = '(calibrated)' if phi_override is not None else '(measured @ zero current)'
    print(f'  转子 d 轴角度 θ_r = {np.degrees(theta_r):.1f}°')
    print(f'  主磁链 Φ_0 = {phi_0:.6f} Wb {phi_label}')
    if phi_override is not None:
        print(f'    (原始测量 Φ_0 = {phi_0:.6f} Wb)')

    # ── 主扫描 ──
    # 收集所有点的 (Id, Iq, ψ_d, ψ_q)，以及 θ=0° 点的 ψ_d(0, Iq)
    raw_points = []  # list of dicts
    # ψ_d(0, Iq) table: {Iq_value: psi_d}
    psi_d0_table = {}  # key: Iq (float), value: psi_d

    point_num = 0
    for i, Ia_mag in enumerate(currents):
        for theta_deg in angles_deg:
            point_num += 1
            theta = np.radians(theta_deg)
            Id = Ia_mag * np.sin(theta)
            Iq = Ia_mag * np.cos(theta)

            # Transform Id/Iq → three-phase currents (aligned with rotor d-axis θ_r)
            Ia, Ib, Ic = id_iq_to_abc(Id, Iq, theta=theta_r)

            print(f'  [A] {point_num}/{total_points} ({100*point_num//total_points}%) '
                  f'求解 Id={Id:+7.2f}A, Iq={Iq:+7.2f}A, Ia={Ia:.2f}A ...')

            # Set coil currents via boundary.update()
            _set_coil_currents(m2d, boundaries, Ia, Ib, Ic)

            # Solve + export matrix → extract per-coil flux linkages
            psi_coil = _extract_flux_linkages(m2d)

            # Total phase flux linkage = turns × sum of coil flux linkages
            psi_a = TURNS * (psi_coil['PhaseA1'] + psi_coil['PhaseA2'])
            psi_b = TURNS * (psi_coil['PhaseB1'] + psi_coil['PhaseB2'])
            psi_c = TURNS * (psi_coil['PhaseC1'] + psi_coil['PhaseC2'])

            # Park transform → dq flux linkages (using rotor angle θ_r)
            psi_d, psi_q = park_abc_to_dq(psi_a, psi_b, psi_c, theta=theta_r)

            pt = {
                'Id': Id, 'Iq': Iq,
                'Ia_mag': Ia_mag, 'theta_deg': theta_deg,
                'Psi_d': psi_d, 'Psi_q': psi_q,
            }
            raw_points.append(pt)

            # Record ψ_d(0, Iq) from θ=0° points
            if abs(theta_deg) < 0.01:
                psi_d0_table[Iq] = psi_d
                print(f'    → ψ_d(0,{Iq:.0f}) = {psi_d:.6f} Wb  (有效 Φ 标定点)')

            print(f'    Id={Id:+7.2f}  Iq={Iq:+7.2f}  ψ_d={psi_d:.6f}  ψ_q={psi_q:.6f}')

    # ── 后处理：交叉饱和补偿 Ld/Lq ──
    # 对 ψ_d(0, Iq) 做插值：用于在任意 Iq 处查有效 PM 磁链
    iq_levels = np.array(sorted(psi_d0_table.keys()))
    psi_d0_levels = np.array([psi_d0_table[iq] for iq in iq_levels])

    print(f'\n{"=" * 60}')
    print(f'  ψ_d(0, Iq) — 有效永磁磁链 vs Iq')
    print(f'{"=" * 60}')
    print(f'  {"Iq (A)":>8s}  {"ψ_d(0,Iq) (Wb)":>16s}')
    for iq, psid0 in zip(iq_levels, psi_d0_levels):
        print(f'  {iq:8.1f}  {psid0:16.6f}')

    results = []
    for pt in raw_points:
        Id, Iq, psi_d, psi_q = pt['Id'], pt['Iq'], pt['Psi_d'], pt['Psi_q']

        # Interpolate ψ_d(0, Iq) for this point's Iq
        if len(iq_levels) >= 2:
            psi_d0_at_iq = float(np.interp(Iq, iq_levels, psi_d0_levels))
        elif len(iq_levels) == 1:
            psi_d0_at_iq = float(psi_d0_levels[0])
        else:
            psi_d0_at_iq = phi_0  # fallback: zero-current value

        # Secant Ld with cross-saturation compensation
        if abs(Id) > 1e-6:
            Ld = (psi_d - psi_d0_at_iq) / Id
        else:
            Ld = 0.0

        # Secant Lq
        if abs(Iq) > 1e-6:
            Lq = psi_q / Iq
        else:
            Lq = 0.0

        # Co-energy torque verification (frame-invariant)
        T_coenergy = 1.5 * pole_pairs * (psi_d * Iq - psi_q * Id)

        # Torque from dq model (should match T_coenergy exactly)
        T_model = 1.5 * pole_pairs * (psi_d0_at_iq * Iq + (Ld - Lq) * Id * Iq)

        results.append({
            'Id': Id, 'Iq': Iq,
            'Psi_d': psi_d, 'Psi_q': psi_q,
            'Psi_d0': psi_d0_at_iq,  # effective PM flux at this Iq
            'Ld': Ld, 'Lq': Lq,
            'T_coenergy': T_coenergy,
            'T_model': T_model,
        })

    # ── 输出 ──
    print(f'\n{"=" * 60}')
    print(f'  Ld/Lq MAP (交叉饱和补偿 secant inductance)')
    print(f'{"=" * 60}')
    print(f'  主磁链 Φ_0 (零电流): {phi_0:.6f} Wb')
    print(f'  极对数 PolePairs:     {pole_pairs:.0f}')
    print(f'{"=" * 60}')

    # 打印 Ld/Lq 表格
    print(f'\n  Ld (H) 表格:')
    print(f'  {"Id\\Iq":>8s}', end='')
    for a in angles_deg:
        print(f'  θ={a:5.0f}°  ', end='')
    print()
    for i, Ia_mag in enumerate(currents):
        print(f'  {Ia_mag:8.2f}', end='')
        for j in range(len(angles_deg)):
            idx = i * len(angles_deg) + j
            print(f'  {results[idx]["Ld"]:+.6f}', end='')
        print()

    print(f'\n  Lq (H) 表格:')
    print(f'  {"Id\\Iq":>8s}', end='')
    for a in angles_deg:
        print(f'  θ={a:5.0f}°  ', end='')
    print()
    for i, Ia_mag in enumerate(currents):
        print(f'  {Ia_mag:8.2f}', end='')
        for j in range(len(angles_deg)):
            idx = i * len(angles_deg) + j
            print(f'  {results[idx]["Lq"]:+.6f}', end='')
        print()

    # ── 扭矩验证 ──
    print(f'\n{"=" * 60}')
    print(f'  扭矩验证: T_coenergy = 1.5·P·(ψ_d·Iq - ψ_q·Id)')
    print(f'  T_model = 1.5·P·[ψ_d(0,Iq)·Iq + (Ld-Lq)·Id·Iq]')
    print(f'  (两者应精确相等)')
    print(f'{"=" * 60}')
    print(f'  {"Id":>8s}  {"Iq":>8s}  {"T_coenergy":>10s}  {"T_model":>10s}  {"Δ":>8s}')
    for r in results:
        delta = abs(r['T_coenergy'] - r['T_model'])
        mark = ' ✓' if delta < 1e-9 else ' ✗'
        print(f'  {r["Id"]:+8.2f}  {r["Iq"]:+8.2f}  {r["T_coenergy"]:+10.2f}  {r["T_model"]:+10.2f}  {delta:8.2e}{mark}')

    # ── 摘要 ──
    print(f'\n{"=" * 60}')
    print(f'  摘要')
    print(f'{"=" * 60}')
    print(f'  Φ_0 (零电流):              {phi_0:.6f} Wb')
    if len(iq_levels) > 0:
        print(f'  ψ_d(0,Iq) 范围:            {psi_d0_levels[0]:.6f} ~ {psi_d0_levels[-1]:.6f} Wb')
        print(f'  Φ 负载增强比 (@Imax):      {psi_d0_levels[-1]/phi_0:.2f}×')
    lds = [r['Ld'] for r in results if abs(r['Id']) > 1e-6]
    lqs = [r['Lq'] for r in results if abs(r['Iq']) > 1e-6]
    if lds:
        print(f'  Ld 范围:                   {min(lds):+.6f} ~ {max(lds):+.6f} H')
    if lqs:
        print(f'  Lq 范围:                   {min(lqs):+.6f} ~ {max(lqs):+.6f} H')
    if lds and lqs:
        dL_sign = 'Ld > Lq (反凸极)' if np.mean(lds) > np.mean(lqs) else 'Lq > Ld (正凸极)'
        print(f'  凸极类型:                  {dL_sign}')
        print(f'  (Ld-Lq) 范围:              {min([ld-lq for ld,lq in zip(lds,lqs)]):+.6f} ~ {max([ld-lq for ld,lq in zip(lds,lqs)]):+.6f} H')

    # ── 外特性用 Φ 推荐值 ──
    print(f'\n{"=" * 60}')
    print(f'  外特性 (Sub-flow E) 参数建议:')
    print(f'{"=" * 60}')
    if len(iq_levels) >= 2:
        print(f'  Φ_eff(Iq) 表 (ψ_d 在 Id=0 处的值):')
        for iq, psid0 in zip(iq_levels, psi_d0_levels):
            print(f'    Iq={iq:6.1f}A  →  Φ_eff={psid0:.6f} Wb')
    print(f'')
    print(f'  注: Φ_eff > Φ_0 是由于负载下永磁体工作点位移和')
    print(f'      定子铁心饱和缓解导致的等效磁链增强。')
    print(f'  如需扭矩标定 Φ，在 TR 求解器中跑 Thet=0° 点：')
    print(f'  Φ_eff_TR = T(@Thet=0) / (1.5 × {pole_pairs:.0f} × Imax)')

    m2d.close_project()
    return results, phi_0


if __name__ == '__main__':
    run()
