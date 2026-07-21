"""Sub-flow D: 全域工作特性效率 MAP
求解器: TransientXY (5_Partial_motor_TR)
方法: (转速 × 扭矩) 二维网格扫描
      每个工作点的 Id/Iq 由 **MTPA 算法** 从目标扭矩生成,
      而非经验系数估算.
      复用同一 m2d 对象提升求解速度.
      FEA 提取 Moving1.Torque 单值平均值,
      Vs / M / PF / 损耗 全部基于 MTPA 的 Id/Iq 解析计算.
      (PyAEDT gRPC 在加载 Transient 下只返回单值时点)

用法:
    python -c "from scripts.subflow_d_efficiency_map import run; r=run(speed_steps=3, torque_steps=3)"
    python -c "from scripts.subflow_d_efficiency_map import run; r=run(speed_steps=8, torque_steps=8)"
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
        print(f'  [D] .aedt 文件已同步: AV={target_rpm}rpm ({n_av}处), Speed_rpm={target_rpm} ({n_sp}处)')

# --- 电机参数 (默认值, 从 Sub-flow A/B/C 结果校准) ---
# 校准后的默认参数 (来自 FEA @350A 验证数据, 2026-07-19)
# MS 测量: Φ=0.0815 Wb, Lq≈0.00155 H. FEA 扭矩校准: Φ_eff=0.0989 Wb.
# 二次饱和模型: Ld(Id)=Ld0 + k1*|Id| + k2*|Id|^2
#   拟合自 FEA 7 点扭矩数据, 总误差 45.2 Nm (线性模型 118.2 Nm)
#   MTPA 角误差 1.5° vs FEA 50° (线性模型 3.5°)
DEFAULT_Ld = 0.002056    # H — 二次饱和模型截距 (Ld at Id=0)
DEFAULT_Lq = 0.00155     # H — MS 测量值 (Id=0, Iq=350)
DEFAULT_Phi = 0.098857   # Wb — FEA 扭矩校准值 (T@0deg=207.6 Nm)
DEFAULT_Rs = 0.05        # Ω — 相电阻
DEFAULT_Ld_slope = -2.09e-6   # H/A — 线性饱和系数 (k1)
DEFAULT_Ld_slope2 = 1.03e-8   # H/A² — 二次饱和系数 (k2)


def _resample_to_torque_grid(results, speed_points, torque_steps):
    """将 (speed × Is) FEA 网格重采样到 (speed × 扭矩) 均匀网格.

    对每个转速行, 利用该行 T_avg(Is) 实测关系, 反插值到均匀扭矩点上.
    超出该转速最大扭矩的点返回 None.
    """
    n_speeds = len(speed_points)
    n_is = results.shape[1]

    # 找全局最大扭矩
    all_tq = []
    for i in range(n_speeds):
        for j in range(n_is):
            r = results[i, j]
            if r and r.get('T_avg', 0) > 0.1:
                all_tq.append(r['T_avg'])
    if not all_tq:
        return None, None
    T_max = max(all_tq)

    T_points = np.linspace(0, T_max, torque_steps)

    # 需要插值的量 (在 results dict 中的 key)
    interp_keys = ['Is', 'beta', 'eta', 'M', 'PF', 'Vs', 'Id', 'Iq',
                   'P_cu', 'P_fe', 'P_mag', 'P_total']

    resampled = np.zeros((n_speeds, torque_steps), dtype=object)

    for i in range(n_speeds):
        # 收集该转速行所有有效点的 T_avg 和各量
        tq_row, vals = [], {k: [] for k in interp_keys}
        for j in range(n_is):
            r = results[i, j]
            if r is None:
                continue
            # 过滤求解失败点: Is>0 但 T≈0
            if r.get('Is', 0) > 1 and r.get('T_avg', 0) < 0.5:
                continue
            tq_row.append(r.get('T_avg', 0))
            for k in interp_keys:
                vals[k].append(r.get(k, 0))

        if len(tq_row) < 2:
            continue

        # 按扭矩单调排序 (确保插值有效)
        order = np.argsort(tq_row)
        tq_sorted = np.array(tq_row)[order]

        for j, T_target in enumerate(T_points):
            if T_target < 1e-6:
                # 零扭矩点
                entry = {k: 0.0 for k in interp_keys}
                entry['T_avg'] = 0.0
                resampled[i, j] = entry
            elif T_target <= tq_sorted[-1]:
                entry = {'T_avg': T_target}
                for k in interp_keys:
                    v_sorted = np.array(vals[k])[order]
                    entry[k] = float(np.interp(T_target, tq_sorted, v_sorted))
                resampled[i, j] = entry
            else:
                # 超出该转速最大扭矩 → NaN
                resampled[i, j] = None

    return resampled, T_points


def _mtpa_for_torque(T_target, pole_pairs, Ld, Lq, Phi, imax, Ld_slope=0.0, Ld_slope2=0.0):
    """MTPA 算法: 给定目标扭矩, 返回 (Id, Iq, Is, β) 使电流最小.

    支持两种凸极类型:
      - Lq > Ld (传统 IPM): MTPA 角 < 0 (弱磁, 利用磁阻扭矩)
      - Ld > Lq (反凸极):   MTPA 角 > 0 (增磁, 利用磁阻扭矩)

    饱和模型: Ld(Id) = Ld + Ld_slope*|Id| + Ld_slope2*|Id|^2
    - Ld_slope2=0, Ld_slope=0 → 恒定 Ld/Lq, 解析二次求解
    - Ld_slope2=0, Ld_slope>0 → 线性饱和, 解析三次求解
    - Ld_slope2>0            → 二次饱和, 数值优化 (黄金分割搜索)
    """
    if T_target <= 0:
        return 0.0, 0.0, 0.0, 0.0

    dL = Ld - Lq  # Ld-Lq: >0 反凸极, <0 正常 IPM
    GOLDEN = (np.sqrt(5) - 1) / 2  # 黄金分割比

    def _Ld_at(Id):
        return Ld + Ld_slope * abs(Id) + Ld_slope2 * Id * Id

    def _torque(Id, Iq):
        Ld_eff = _Ld_at(Id)
        return 1.5 * pole_pairs * (Phi * Iq + (Ld_eff - Lq) * Id * Iq)

    def _torque_at_Is_beta(Is_val, beta_rad):
        """给定 Is 和电流角 β, 返回扭矩 + (Id, Iq)."""
        Id_val = Is_val * np.sin(beta_rad)
        Iq_val = Is_val * np.cos(beta_rad)
        return _torque(Id_val, Iq_val), Id_val, Iq_val

    if Ld_slope2 > 1e-15:
        # ── 二次饱和模式: Ld(Id)=Ld0 + k1*|Id| + k2*|Id|^2 ──
        # 黄金分割搜索找最优 β (0~80° for 反凸极, -60~0° for 传统 IPM)
        # 外层: 二分搜索 Is 匹配目标扭矩.
        def _torque_mtpa_at_Is(Is_val):
            """给定 Is, 黄金分割搜索最优 β, 返回 (T_max, Id_opt, Iq_opt)."""
            if dL > 0:
                lo, hi = 0.0, np.radians(80.0)
            else:
                lo, hi = np.radians(-60.0), 0.0

            # 黄金分割搜索: 30 次迭代足够收敛到 1e-6 rad
            b = hi - GOLDEN * (hi - lo)
            a = lo + GOLDEN * (hi - lo)
            T_b, Id_b, Iq_b = _torque_at_Is_beta(Is_val, b)
            T_a, Id_a, Iq_a = _torque_at_Is_beta(Is_val, a)
            for _ in range(35):
                if T_b > T_a:
                    hi = a; a = b; T_a = T_b; Id_a, Iq_a = Id_b, Iq_b
                    b = hi - GOLDEN * (hi - lo)
                    T_b, Id_b, Iq_b = _torque_at_Is_beta(Is_val, b)
                else:
                    lo = b; b = a; T_b = T_a; Id_b, Iq_b = Id_a, Iq_a
                    a = lo + GOLDEN * (hi - lo)
                    T_a, Id_a, Iq_a = _torque_at_Is_beta(Is_val, a)
            beta_mid = (lo + hi) / 2.0
            return _torque_at_Is_beta(Is_val, beta_mid)

        T_max, Id_max, Iq_max = _torque_mtpa_at_Is(imax)
        if T_target >= T_max:
            Is, Id, Iq = imax, Id_max, Iq_max
        else:
            lo, hi = 0.0, imax
            for _ in range(60):
                mid = (lo + hi) / 2.0
                T_mid, _, _ = _torque_mtpa_at_Is(mid)
                if T_mid > T_target:
                    hi = mid
                else:
                    lo = mid
                if abs(T_mid - T_target) < 1e-6:
                    break
            Is = (lo + hi) / 2.0
            _, Id, Iq = _torque_mtpa_at_Is(Is)
        beta = np.degrees(np.arctan2(Id, Iq))

    elif Ld_slope > 1e-12:
        # ── 线性饱和模式: 基于 dL(Id)=dL0+k*|Id| 的 MTPA 条件 ──
        # dT/dβ=0 → 三次方程: 3k*Is²*s³ + 2dL0*Is*s² + (Φ-2k*Is²)*s - dL0*Is = 0
        # 其中 s=sinβ, dL0=Ld-Lq, k=Ld_slope.
        #
        # 方法: 外层二分搜索 Is, 内层解三次方程求 β_opt(Is).
        def _torque_mtpa_at_Is(Is_val):
            """给定 Is, 返回 MTPA 最优角下的扭矩."""
            dL0 = Ld - Lq
            k = Ld_slope
            a3 = 3 * k * Is_val * Is_val
            a2 = 2 * dL0 * Is_val
            a1 = Phi - 2 * k * Is_val * Is_val
            a0 = -dL0 * Is_val
            if abs(a3) < 1e-15:
                # Fall through to constant-dL case
                if abs(dL0) > 1e-12:
                    disc = Phi*Phi + 8*dL0*dL0*Is_val*Is_val
                    s = (-Phi + np.sqrt(disc)) / (4*dL0*Is_val) if dL0 > 0 else \
                        (Phi - np.sqrt(disc)) / (4*dL0*Is_val)
                    s = max(-1.0, min(1.0, s))
                else:
                    s = 0.0
            else:
                roots = np.roots([a3, a2, a1, a0])
                s = 0.0
                for r in roots:
                    if abs(r.imag) < 1e-10 and -1.0 <= r.real <= 1.0:
                        # Choose the root giving max torque
                        beta_try = np.arcsin(r.real)
                        Id_try = Is_val * np.sin(beta_try)
                        Iq_try = Is_val * np.cos(beta_try)
                        T_try = _torque(Id_try, Iq_try)
                        cur_max = _torque(Is_val*np.sin(np.arcsin(s)), Is_val*np.cos(np.arcsin(s))) if abs(s)>1e-12 else 0
                        if T_try > cur_max:
                            s = r.real
            beta_opt = np.arcsin(s)
            Id_opt = Is_val * np.sin(beta_opt)
            Iq_opt = Is_val * np.cos(beta_opt)
            return _torque(Id_opt, Iq_opt), Id_opt, Iq_opt

        T_max, Id_max, Iq_max = _torque_mtpa_at_Is(imax)
        if T_target >= T_max:
            Is, Id, Iq = imax, Id_max, Iq_max
        else:
            lo, hi = 0.0, imax
            for _ in range(60):
                mid = (lo + hi) / 2.0
                T_mid, _, _ = _torque_mtpa_at_Is(mid)
                if T_mid > T_target:
                    hi = mid
                else:
                    lo = mid
                if abs(T_mid - T_target) < 1e-6:
                    break
            Is = (lo + hi) / 2.0
            _, Id, Iq = _torque_mtpa_at_Is(Is)
        beta = np.degrees(np.arctan2(Id, Iq))

    elif abs(dL) > 1e-12:
        # ── 标准 MTPA 轨迹法 (恒定 Ld/Lq) ──
        # 统一公式: sinβ = [-Φ + √(Φ²+8dL²Is²)] / (4·dL·Is)
        # 对给定 Is 求 MTPA 角 → 扭矩 → 二分搜索 Is 匹配 T_target.
        def _torque_mtpa_const_at_Is(Is_val):
            if abs(dL) < 1e-12:
                return Phi * Is_val, 0.0, Is_val
            disc = Phi*Phi + 8*dL*dL*Is_val*Is_val
            # Pick the root giving max torque: +sqrt for dL>0, -sqrt for dL<0
            if dL > 0:
                s = (-Phi + np.sqrt(disc)) / (4 * dL * Is_val)
            else:
                s = (Phi - np.sqrt(disc)) / (4 * dL * Is_val)
            s = max(-1.0, min(1.0, s))
            beta_opt = np.arcsin(s)
            Id_opt = Is_val * np.sin(beta_opt)
            Iq_opt = Is_val * np.cos(beta_opt)
            T_opt = 1.5 * pole_pairs * (Phi * Iq_opt + dL * Id_opt * Iq_opt)
            return T_opt, Id_opt, Iq_opt

        T_max, Id_max, Iq_max = _torque_mtpa_const_at_Is(imax)
        if T_target >= T_max:
            Is, Id, Iq = imax, Id_max, Iq_max
        else:
            lo, hi = 0.0, imax
            for _ in range(60):
                mid = (lo + hi) / 2.0
                T_mid, _, _ = _torque_mtpa_const_at_Is(mid)
                if T_mid > T_target:
                    hi = mid
                else:
                    lo = mid
                if abs(T_mid - T_target) < 1e-6:
                    break
            Is = (lo + hi) / 2.0
            _, Id, Iq = _torque_mtpa_const_at_Is(Is)
        beta = np.degrees(np.arctan2(Id, Iq))
    else:
        # ── 非凸极 (Ld ≈ Lq) ──
        # T = 1.5*P*Φ*Iq, Id=0
        Iq = T_target / (1.5 * pole_pairs * Phi)
        if Iq > imax:
            Iq = imax
        Id = 0.0
        Is = Iq
        beta = 0.0

    return Id, Iq, Is, beta


def _calc_analytical(pole_pairs, speed, Id, Iq, T_avg,
                     vdc=300, Ld=DEFAULT_Ld, Lq=DEFAULT_Lq,
                     Phi=DEFAULT_Phi, Rs=DEFAULT_Rs):
    """解析电压方程: 基于 (Id, Iq) 计算 Vs / M / PF / 损耗 / 效率"""
    Imax = np.sqrt(Id**2 + Iq**2)
    omega_e = speed * np.pi / 30 * pole_pairs   # 电角速度 (rad/s)

    # 电压方程 (稳态, 忽略电阻压降)
    Vd = -omega_e * Lq * Iq
    Vq = omega_e * (Ld * Id + Phi)
    Vs = np.sqrt(Vd**2 + Vq**2)

    Vmax = vdc / np.sqrt(3)
    M = Vs / Vmax if Vmax > 0 else 0.0

    # 功率因数角 = 电压角 - 电流角
    phi_v = np.arctan2(Vd, Vq)
    phi_i = np.arctan2(Id, Iq)
    PF = np.cos(phi_v - phi_i)

    # 功率 / 损耗
    P_out = T_avg * speed * 2 * np.pi / 60        # W
    P_cu = 3 * (Imax / np.sqrt(2))**2 * Rs        # W
    P_loss = P_cu
    eta = P_out / (P_out + P_loss) * 100 if (P_out + P_loss) > 0 else 0.0

    return {
        'Id': Id, 'Iq': Iq, 'Is': Imax,
        'beta': np.degrees(np.arctan2(Id, Iq)),  # MTPA 电流角(°)
        'Vd': Vd, 'Vq': Vq,
        'Vs': Vs / 1000 if Vs > 1000 else Vs,
        'PF': PF, 'M': M,
        'P_out': P_out / 1000,
        'P_cu': P_cu / 1000,
        'P_loss': P_loss / 1000,
        'eta': eta,
    }


def run(speed_min=1000, speed_max=8000, speed_steps=5, torque_steps=5,
        vdc=300, imax=250,
        Ld=DEFAULT_Ld, Lq=DEFAULT_Lq, Phi=DEFAULT_Phi, Rs=DEFAULT_Rs,
        Ld_slope=DEFAULT_Ld_slope, Ld_slope2=DEFAULT_Ld_slope2,
        project_path=None):
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
    # D 流程多转速扫描，先以最低转速初始化文件，循环中再通过 boundary API 逐点同步 AV
    _sync_motion_angular_velocity(template, int(speed_min))

    m2d = Maxwell2d(
        project=template, design='5_Partial_motor_TR',
        solution_type='TransientXY',
        non_graphical=False, new_desktop=False, close_on_exit=False
    )
    # PolePairs = 4 (固定, 8-pole PMSM, 模板参数)
    pole_pairs = 4
    theta_45 = np.radians(45)  # 固定电流角

    speed_points = np.linspace(speed_min, speed_max, speed_steps)
    # Imax 线性扫描 (替代 MTPA 扭矩目标)
    Imax_points = np.linspace(0, imax, torque_steps)

    results = np.full((speed_steps, torque_steps), None, dtype=object)
    total = speed_steps * torque_steps
    count = 0

    # ── 断点恢复: 检查 checkpoint 文件 ──
    import json, time as _time
    ckpt_dir = Path(project_path if project_path else '.').resolve()
    ckpt_file = ckpt_dir / '_subflow_d_checkpoint.json'
    start_i, start_j = 0, 0
    if ckpt_file.exists():
        try:
            ckpt = json.loads(ckpt_file.read_text(encoding='utf-8'))
            if (ckpt.get('speed_points') == list(speed_points)
                    and ckpt.get('Imax_points') == list(Imax_points)):
                for entry in ckpt.get('results', []):
                    i, j, d = entry['i'], entry['j'], entry['data']
                    results[i, j] = d
                    if d: count += 1
                start_i = ckpt.get('next_i', 0)
                start_j = ckpt.get('next_j', 0)
                print(f'[D] 从断点恢复: ({start_i},{start_j}), 已有 {count}/{total} 点')
        except Exception:
            pass

    if count == 0:
        print(f'效率 MAP 扫描: {speed_steps}×{torque_steps} = {total} 点')
    print(f'转速范围: {speed_min}~{speed_max} rpm')
    print(f'Vdc={vdc}V, Imax={imax}A, Thet=45° (固定电流角)')
    print(f'Ld={Ld*1000:.4f}mH, Lq={Lq*1000:.4f}mH, Phi={Phi:.4f}Wb, Rs={Rs:.4f}Ω')
    print(f'Imax 步进: {", ".join(f"{v:.0f}" for v in Imax_points)} A')
    print(f'{"=" * 70}')

    def _save_checkpoint():
        entries = []
        for i2 in range(speed_steps):
            for j2 in range(torque_steps):
                if results[i2, j2] is not None:
                    entries.append({'i': i2, 'j': j2, 'data': results[i2, j2]})
        ckpt = {
            'speed_points': list(speed_points),
            'Imax_points': list(Imax_points),
            'next_i': i, 'next_j': j + 1 if j + 1 < torque_steps else 0,
            'results': entries,
        }
        ckpt_file.write_text(json.dumps(ckpt, ensure_ascii=False), encoding='utf-8')

    def _solve_one_point(max_retries=3):
        """带 gRPC 重试的单点 FEA 求解"""
        for attempt in range(max_retries):
            try:
                m2d.analyze('Setup1')
                data = m2d.post.get_solution_data_per_variation(
                    expressions=['Moving1.Torque']
                )
                if data:
                    tq = np.array(data.data_real('Moving1.Torque'))
                    if len(tq) >= 1:
                        return float(np.mean(tq))
                return 0.0
            except Exception as e:
                err_msg = str(e)
                if 'GrpcApiError' in type(e).__name__ or 'gRPC' in err_msg:
                    if attempt < max_retries - 1:
                        wait = 2 ** attempt
                        print(f' [gRPC重试{attempt+1}/{max_retries}, 等{wait}s]',
                              end='', flush=True)
                        _time.sleep(wait)
                        continue
                raise
        return 0.0

    for i in range(start_i, speed_steps):
        n = speed_points[i]
        j_start = start_j if i == start_i else 0
        for j in range(j_start, torque_steps):
            Is_val = Imax_points[j]
            count = sum(1 for r in results.flat if r is not None)
            current = count + 1  # 当前正在做的点编号
            Id_45 = Is_val * np.sin(theta_45)
            Iq_45 = Is_val * np.cos(theta_45)

            # 设置 FEA 参数 — 固定 Thet=45°, Imax 线性扫描
            check_var('Speed_rpm', 'D')
            m2d['Speed_rpm'] = str(int(n))
            b = [b for b in m2d.boundaries if b.name == 'MotionSetup1'][0]
            b.props['Angular Velocity'] = f'{int(n)}rpm'
            b.update()
            check_var('Imax', 'D')
            m2d['Imax'] = f'{Is_val}A'
            check_var('Thet_deg', 'D')
            m2d['Thet_deg'] = '45'

            freq = n / 60 * pole_pairs
            setup = m2d.setups[0]
            setup.props['StopTime'] = f'{1/freq}s'
            setup.props['TimeStep'] = f'{1/(freq*50)}s'
            setup.update()

            print(f'[D] {current:3d}/{total} ({100*current//total:3d}%) '
                  f'n={n:5.0f} Is={Is_val:5.0f}A '
                  f'Id={Id_45:+6.1f} Iq={Iq_45:+6.1f} Thet=45°',
                  end='', flush=True)

            if Is_val < 1e-6:
                T_avg = 0.0
            else:
                T_avg = _solve_one_point()

            # ── 解析电参数 (固定 Thet=45°) ──
            pt = _calc_analytical(
                pole_pairs, n, Id_45, Iq_45, T_avg,
                vdc=vdc, Ld=Ld, Lq=Lq, Phi=Phi, Rs=Rs
            )
            pt['T_avg'] = T_avg
            pt['Is'] = Is_val
            pt['beta'] = 45.0
            results[i, j] = pt

            print(f' → T={T_avg:6.1f} Nm '
                  f'Vs={pt["Vs"]:.0f}V M={pt["M"]:.2f} '
                  f'PF={pt["PF"]:.3f} η={pt["eta"]:.1f}%')

            # 每个点完成后保存断点
            _save_checkpoint()

    # 全部完成, 清理断点文件
    if ckpt_file.exists():
        ckpt_file.unlink()

    m2d.close_project()

    # ============ 重采样: (speed × Is) → (speed × 扭矩) ============
    resampled, T_points = _resample_to_torque_grid(results, speed_points, torque_steps)

    out_dir = ROOT / 'results'
    out_dir.mkdir(exist_ok=True)

    # ── 先保存原始 FEA 数据 (speed×Is) ──
    Is_labels = [f'{t:.0f}' for t in Imax_points]
    raw_keys = {
        'raw_torque': ('T_avg (Nm)', 'T_avg'),
        'raw_efficiency': ('eta (%)', 'eta'),
    }
    for fname, (hdr, key) in raw_keys.items():
        rows = [f'# {hdr}  ─ 行:转速(rpm), 列:Is(A), Thet=45° (FEA原始)']
        rows.append(',' + ','.join(Is_labels))
        for i, n in enumerate(speed_points):
            vals = [str(results[i, j][key]) if results[i, j] else ''
                    for j in range(torque_steps)]
            rows.append(f'{n:.0f},{",".join(vals)}')
        (out_dir / f'{fname}.csv').write_text('\n'.join(rows), encoding='utf-8')

    # ── 输出重采样 MAP (speed × 扭矩) ──
    if resampled is None:
        print('⚠ 无有效 FEA 数据, 无法生成 MAP')
        return results

    T_labels = [f'{t:.0f}' for t in T_points]

    print(f'\n{"=" * 70}')
    print(f'=== 效率 MAP  重采样到 (转速 × 扭矩) 网格 ===')
    print(f'Vdc={vdc}V, Imax={imax}A, Thet=45°')
    print(f'行: {len(speed_points)} 转速 ({speed_points[0]:.0f}~{speed_points[-1]:.0f} rpm)')
    print(f'列: {len(T_points)} 扭矩 (0~{T_points[-1]:.0f} Nm)')
    print(f'{"=" * 70}')

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
                s = f'{v:{fmt}}' if not (isinstance(v, float) and np.isnan(v)) else '     nan'
                print(f'{s:>8s}', end='')
            print()

    _print_map('效率 η (%)', 'eta', '8.1f')
    _print_map('FEA 扭矩 T_avg (Nm)', 'T_avg', '8.1f')
    _print_map('电流 Is (A)', 'Is', '8.1f')
    _print_map('调制比 M', 'M', '8.3f')
    _print_map('功率因数 PF', 'PF', '8.3f')
    _print_map('相电压 Vs (V)', 'Vs', '8.0f')

    # 峰值效率
    eta_vals = [r['eta'] for r in resampled.flat if r]
    if eta_vals:
        print(f'\n峰值效率: {np.max(eta_vals):.1f}%')

    # ── 保存 CSV (speed × 扭矩) ──
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
    print(f'  原始数据: raw_torque.csv, raw_efficiency.csv (speed×Is)')
    print(f'  重采样:   efficiency_map.csv, torque_map.csv, current_map.csv,')
    print(f'            modulation_map.csv, pf_map.csv, voltage_map.csv (speed×T)')
    return resampled


if __name__ == '__main__':
    run()
