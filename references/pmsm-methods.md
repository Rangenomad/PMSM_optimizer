# PMSM Optimizer — 子流程计算模板

本文档包含 5 个子流程的方法说明和关键公式。完整代码见 `scripts/` 目录。

---

## Transient 求解器公共配置

所有基于 TransientXY (`5_Partial_motor_TR`) 的子流程共用以下时间参数配置规则。

### 电频率

```
freq = Speed_rpm / 60 × PolePairs      (Hz)
```

### 停止时间

每个工况点的停止时间 = N 个完整电周期：

```
StopTime = elec_periods / freq           (s)
```

### 时间步长

每电周期固定步数，步长随转速自动缩放：

```
TimeStep = StopTime / (elec_periods × steps_per_period)    (s)
```

或等价于：

```
TimeStep = 1 / (freq × steps_per_period)                   (s)
```

### 各子流程参数

| 子流程 | elec_periods | steps_per_period | 说明 |
|--------|-------------|-----------------|------|
| B (空载反电势) | 2-3 | 50 | 需多周期 FFT 分析 |
| C (额定扭矩 — MTPA 扫描) | 1 | 50 | 快速筛选 |
| C (额定扭矩 — 精确结算) | 3 | 50 | 稳态波形 + 扭矩脉动 |
| D (效率 MAP) | 1 | 50 | 单周期取平均扭矩 |

### 物理含义

- **每点固定 N 个电周期**，不随转速变化。转速高则 StopTime 短、TimeStep 短；转速低则 StopTime 长、TimeStep 长。**每工况点的计算量基本恒定**。
- 永磁同步电机的电磁转矩波动频率为 6× 电频率，每个电周期有 6 个纹波峰谷。50 步/周期下每个纹波周期约有 8 步，可分辨平均值和粗略脉动率。
- Sub-flow D 效率 MAP 只用平均扭矩，不关心波形，故取 elec_periods=1 最小化求解时间。

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

**时间配置**: Transient 公共配置 (elec_periods=2~3, steps_per_period=50)

**关键公式**:
```
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

**方法**: (转速 × 扭矩) 二维网格扫描。**电流 (Id, Iq) 由 MTPA 算法从目标扭矩生成**，而非经验系数估算。FEA 仅取 Moving1.Torque 单值，电参数用解析电压方程计算（因 PyAEDT gRPC 在加载 Transient 下只返回单值时点）。

**MTPA 算法** (`_mtpa_for_torque`):
```
IPM 扭矩方程:     T = 1.5·P·(Φ·Iq + (Ld-Lq)·Id·Iq)
MTPA 轨迹:       Id = A - √(A² + Iq²)    其中 A = Φ / (2·(Lq-Ld))
方法:           沿 MTPA 轨迹二分搜索 Iq, 使扭矩 = T_target
```

**时间配置**: Transient 公共配置 (elec_periods=1, steps_per_period=50)

**关键公式** (解析电参数):
```
Vd = -ω·Lq·Iq,   Vq = ω·(Ld·Id+Φ),   Vs = √(Vd²+Vq²)
PF = cos(arctan2(Vd,Vq) - arctan2(Id,Iq))
M = Vs / (Vdc/√3)
P_cu = 3 × (Is/√2)² × Rs
```

**输出**: 7 张 MAP（T_avg / η / M / PF / Vs / Is / β），每张 CSV 带行列标签

**模块**: `scripts/subflow_d_efficiency_map.py` — `run(speed_min, speed_max, speed_steps, torque_steps, vdc, imax, Ld, Lq, Phi, Rs)`

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
