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
