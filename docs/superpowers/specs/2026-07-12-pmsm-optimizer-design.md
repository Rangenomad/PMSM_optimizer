# PMSM Optimizer — 设计文档

## 概述

基于 Ansys Maxwell 2D (PyAEDT) 的永磁同步电机电磁参数提取与优化工具。用户提供一个预先配置好的 PMSM 模板模型（.aedt），通过修改设计变量和运行预设子流程，自动计算 Ld/Lq MAP、反电势、额定扭矩、效率 MAP、外特性等。

---

## 架构策略

### 核心决策：模板公式法

由于 PyAEDT 通过 gRPC 协议与 AEDT 通信时，大量 COM 遗留 API（`ChangeProperty`、`SetPropertyValue`、`DeleteBoundary`、`GetPropertyValue`）不可用，采用**一次 GUI 设置 + 自动化变量修改**的策略：

```
一次 GUI 设置（用户手动）       后续自动化（Python 脚本）
─────────────────────────      ─────────────────────────
1. 绕组激励公式写入              1. 修改设计变量（电流、转速等）
   引用设计变量                   2. 求解 → 后处理提取结果
2. 材料属性设置                  3. 批量扫描（Sub-flow A/D）
3. 边界条件                     4. 数学计算（Sub-flow E）
```

**用户只需在 GUI 中设置一次绕组公式，后续所有操作通过 Python 修改变量完成。**

### 为什么不用 assign_current()？

`assign_current()` 需要用户知道拓扑细节（每相线圈数、每个槽的导体数、并联支路数），而用户只提供相电流。Winding Group 方式允许在 GUI 中一次性配置好绕组拓扑关系，后续只需通过变量改变电流幅值和角度。

---

## 模板结构

模板文件 `Prius_2D_Practice.aedt` 包含 3 个设计，从原模板清理后保留：

| 设计名 | 求解器类型 | 用途 | 边界条件 | 绕组 |
|--------|-----------|------|---------|------|
| `4_Partial_motor_MS2` | MagnetostaticXY | Sub-flow A: Ld/Lq MAP | Balloon (Region 外边界) | Current 激励 ×6, Matrix(9匝/线圈) |
| `5_Partial_motor_TR` | TransientXY | Sub-flow B/C/D: 反电势/额定扭矩/效率 MAP | Vector Potential=0 (外弧边), Master/Slave (径向边) | Winding Group + Coil(9匝/线圈) |
| `6_Partial_motor_CT` | TransientXY | Sub-flow B/C/D 备用/对照 | Vector Potential=0 (外弧边), Master/Slave (径向边) | (同上) |

### 关键设计变量

| 变量名 | 公式 | 当前值 | 用途 |
|--------|------|--------|------|
| `Poles` | 8 | 8 | 极数 (固定) |
| `PolePairs` | Poles/2 | 4 | 极对数 (固定) |
| `Speed_rpm` | 3000 | 3000 | 机械转速 (rpm) |
| `Omega` | 360*speed_rpm*PolePairs/60 | 72000 | 电角速度 (°/s) |
| `Omega_rad` | Omega*pi/180 | 1256.63 | 电角速度 (rad/s) |
| `Thet_deg` | 20 | 20 | 电流角/初始角度 (°) |
| `Thet` | thet_deg*pi/180 | 0.349 | 电流角/初始角度 (rad) |
| `Imax` | 250 | 250 | 相电流幅值 (A) |

### 绕组激励公式

已在 GUI 中设置的 Winding Group 公式（用户一次性配置）：

- **Phase_A**: `Imax*sin(Omega_rad*time + Thet)`
- **Phase_B**: `Imax*sin(Omega_rad*time + Thet - 2*pi/3)`
- **Phase_C**: `Imax*sin(Omega_rad*time + Thet + 2*pi/3)`

> 自动化时只需修改 `Imax`、`Thet_deg`、`Speed_rpm` 三个变量即可改变激励。

> **匝数**：每个线圈对象代表 9 匝（Coil 的 Conductor number / Matrix 的 NumberOfTurns = 9）。assign_current() 设的电流是总 MMF，需乘以 9。Winding Group 自动处理该乘法。

### 材料

| 材料 | 用途 | 关键属性 |
|------|------|---------|
| `n36z_20` | 永磁体 (NdFeB N36Z) | Br=1.17-1.21 T, Hc≈-860 kA/m |
| 硅钢片 | 定转子铁芯 | 模板预设 |
| 铜 | 绕组 | 模板预设 |

---

## 子流程设计

### Sub-flow A：Ld/Lq MAP + 主磁链 Φ

**求解器**: MagnetostaticXY (`4_Partial_motor_MS2`)
**原理**: 电流角扫描法。在每个 (Id, Iq) 工作点求解静磁场，从磁链中提取 d/q 轴分量。

**4_MS2 绕组结构**（已通过 Matrix 和 Current 激励预设）：

```
PhaseA: PhaseA1(9匝, 串联) + PhaseA2(9匝, 串联), 1并联支路
PhaseB: PhaseB1(9匝, 串联) + PhaseB2(9匝, 串联), 1并联支路
PhaseC: PhaseC1(9匝, 串联) + PhaseC2(9匝, 串联), 1并联支路
```

每相有 2 个 coil 对象，每个 coil 代表 **9 匝**（Conductor number / NumberOfTurns = 9）。

**实现要点**:

1. Magnetostatic 求解器不支持 Winding Group 激励 → 需要用 `assign_current()` 直接设置电流

2. **电流变换（Id/Iq → Ia/Ib/Ic → 匝数放大）**：
   ```
   d 轴与 A 相对齐（θ=0°）：
   
   Ia =  Id
   Ib = -0.5*Id + 0.866*Iq  
   Ic = -0.5*Id - 0.866*Iq
   
   # assign_current 设置的是总 MMF，每线圈代表 9 匝
   assign_current('PhaseA1', current=9*Ia)
   assign_current('PhaseA2', current=9*Ia)
   assign_current('PhaseB1', current=9*Ib)
   assign_current('PhaseB2', current=9*Ib)
   assign_current('PhaseC1', current=9*Ic)
   assign_current('PhaseC2', current=9*Ic)
   ```

3. 求解后提取三相磁链，需乘匝数得到总磁链：
   ```python
   psi_a = 9 * (FluxLinkage(PhaseA1) + FluxLinkage(PhaseA2))
   psi_b = 9 * (FluxLinkage(PhaseB1) + FluxLinkage(PhaseB2))
   psi_c = 9 * (FluxLinkage(PhaseC1) + FluxLinkage(PhaseC2))
   ```

4. Park 变换将三相磁链转为 d/q 轴分量（θ=0°）：
   ```python
   psi_d =  2/3 * (psi_a*cos(0) + psi_b*cos(-2π/3) + psi_c*cos(+2π/3))  =  2/3*(psi_a - 0.5*psi_b - 0.5*psi_c)
   psi_q = -2/3 * (psi_a*sin(0) + psi_b*sin(-2π/3) + psi_c*sin(+2π/3))  =  2/3*(0.866*psi_b - 0.866*psi_c)
   ```

5. Ld = Ψd/Id, Lq = Ψq/Iq（Id, Iq ≠ 0 时）

**输入参数**:
```
rated_current: 额定电流幅值 (A) — 每根导体的电流，非匝数放大值
max_current:   最大电流倍率 (×额定)
current_steps: 电流幅值分点数
angle_steps:   电流角分点数 (0-90°)
```

**输出**: Ld(Id,Iq) MAP 表格 + Lq(Id,Iq) MAP 表格 + 主磁链 Φ (Wb)

---

### Sub-flow B：空载反电势

**求解器**: TransientXY (`5_Partial_motor_TR`)
**原理**: 电枢绕组开路（Imax=0），永磁体激励，额定转速旋转。

**实现要点**:
1. 设置 `Imax=0`（通过变量），绕组不贡献磁场
2. 设置 `Speed_rpm` 到用户指定转速
3. 求解 2-3 个电周期
4. 提取三相感应电压 `Voltage(A)`, `Voltage(B)`, `Voltage(C)`
5. FFT 分析基波幅值和 THD

**输入参数**:
```
rated_speed:     额定转速 (rpm)
elec_periods:    仿真电周期数 (默认 2)
time_steps:      每周期步数 (默认 200)
```

**输出**: 三相 BEMF 波形 + 基波幅值 + THD + 线反电势常数 Ke (V/krpm)

---

### Sub-flow C：额定点扭矩

**求解器**: TransientXY (`5_Partial_motor_TR`)
**原理**: 额定负载、额定转速，仿真一个完整电周期。

**实现要点**:
1. 设置 `Imax` 到额定电流
2. 如果用户没指定电流角，用 `Thet_deg=0`（或 MTPA 查找）
3. 设置 `Speed_rpm` 到额定转速
4. 求解 3 个电周期（取最后一个周期的稳态数据）
5. 提取 `Moving1.Torque` 波形

**输入参数**:
```
rated_current:  额定电流 (A)
current_angle:  电流角 (°)
rated_speed:    额定转速 (rpm)
elec_periods:   仿真电周期数 (默认 3)
```

**输出**: 平均扭矩 + 扭矩脉动率 + 最大/最小扭矩

---

### Sub-flow D：全域工作特性 MAP

**求解器**: 批量 TransientXY (`5_Partial_motor_TR`)
**原理**: (转速 × 扭矩) 二维网格扫描，每点跑一个电周期。

**实现要点**:
1. 基于 Sub-flow A 的 Ld/Lq MAP 反算每个 (n, T) 工作点的 Id/Iq
2. 批量提交仿真（逐点或并行）
3. 从每个工作点提取：相电流幅值 Is、铜损 Pcu=3×Is²×Rs、铁损（如果模型配置了铁损计算）、电压幅值 Vs、功率因数
4. 计算调制比 M = Vs / (Vdc/√3)

**输入参数**:
```
speed_min/max/steps:  转速网格
torque_steps:         扭矩分点数
vdc:                  直流母线电压 (V)
imax:                 最大相电流 (A)
```

**输出**: 4 张 MAP（损耗 / Is / 功率因数 / 调制比），每张 CSV + 终端预览

---

### Sub-flow E：外特性计算

**求解器**: 纯数学计算（不跑 FEA）
**原理**: 在电压圆 Vs ≤ Vdc/√3 和电流圆 Is ≤ Imax 约束下，基于 Sub-flow A 的 Ld/Lq MAP 计算各转速最大输出扭矩。

**实现要点**:
1. 从 Sub-flow A 获取 Ld(Id,Iq), Lq(Id,Iq), Φ 数据（插值模型）
2. 对每个转速：
   a. 计算电频率 ω = 2π × n × PolePairs / 60
   b. 求解 MTPA 工作点（恒扭矩区）：在电流圆上找扭矩最大点
   c. 求解弱磁区工作点：超出电压限时沿电压圆边界计算
3. 输出 T-n 特性曲线、工作点轨迹、转折点转速、峰值功率

**输入参数**:
```
Ld/Lq/Φ:    来自 Sub-flow A（必选前置）
Vdc:         直流母线电压 (V)
Imax:        最大相电流 (A)
Rs:          相电阻 (Ω)
speed_points: 转速分点数 (默认 20)
```

**输出**: T-n 曲线表 + 工作点轨迹 (Id, Iq) 表 + 转折点 + 峰值功率

---

## 执行流程

```
用户指令 → Step 0: 确认模式 (auto/confirm)
         → Step 1: 连接模板项目
         → Step 2: 修改参数 (可选)
         → Step 3: 判定需求 → 映射子流程
         → 执行子流程 (A/B/C/D/E)
         → 输出结果
```

### 子流程依赖关系

```
Sub-flow A (Ld/Lq MAP)     ← 独立，但 D/E 依赖其结果
Sub-flow B (反电势)        ← 独立
Sub-flow C (额定扭矩)      ← 独立
Sub-flow D (工作特性 MAP)  ← 可选依赖 A 的结果
Sub-flow E (外特性)        ← 必须依赖 A 的输出
```

---

## 模板初始化脚本 (`setup_template.py`)

已在之前会话中完成并验证。功能：

1. 从备份恢复模板文件
2. 清理多余设计（保留 3 个）
3. 4_MS2 → 设置 Balloon 边界
4. 5_TR → 设置 MotionSetup (Speed_rpm 变量引用) + Transient Setup + Vector Potential
5. 6_CT → 同上

**用户仍需在 GUI 中手动完成**：
- Winding Group 电流公式输入（一次性的）
- 磁钢材料属性确认（Br/Hc 是否正确）

---

## 已知的 PyAEDT gRPC API 限制

| API | 状态 | 替代方案 |
|-----|------|---------|
| `ChangeProperty()` | ❌ 不可用 | 通过 `m2d[...] = value` 修改变量 |
| `SetPropertyValue()` | ❌ 不可用 | 同上 |
| `DeleteBoundary()` | ❌ 不可用 | 重建模板 |
| `GetPropertyValue()` | ❌ 不可用 | 通过 `m2d.variable_manager` 读取 |
| `assign_current()` | ✅ 可用 | 用于 Magnetostatic 激励设置 |
| `assign_winding()` | ✅ 可用 | 用于 Transient Winding Group |
| `assign_coil()` | ✅ 可用 | 用于指定线圈对象 |
| `m2d[...]=value` | ✅ 可用 | 设计变量读写 |
| `assign_balloon()` | ✅ 可用 | Magnetostatic 边界 |
| `assign_vector_potential()` | ✅ 可用 | Transient 边界 |
| `analyze()` | ✅ 可用 | 求解控制 |
| `post.get_solution_data()` | ✅ 可用 | 后处理提取数据 |

---

## 关于 Skill.md 和 pmsm-methods.md 的修改计划

### SKILL.md
- 第零步（执行模式）→ 无需修改
- 第一步（连接模板）→ 更新为使用 `Prius_2D_Practice` 而非 `pmsm_template`
- 第二步（修改参数）→ 补充当前模板的实际变量名列表
- 第三步（判定需求）→ 无需修改
- 子流程执行规范 → 需要补充：
  - Sub-flow A 中 Magnetostatic 需用 `assign_current()` 而非 Winding Group 的说明
  - Sub-flow B/C/D 中通过变量修改激励的具体方式
  - 各子流程的实际参数名映射

### pmsm-methods.md
- **Sub-flow A** → 改完善 TODO 代码：
  - 补充 `assign_current()` 设置三相电流
  - 补充 Park 变换磁链提取
  - 补充主磁链 Φ 提取
- **Sub-flow B** → 改完善 TODO 代码：
  - 补充通过变量设置 Imax=0
  - 补充瞬态求解时间设置
  - 补充电压波形提取和 FFT
- **Sub-flow C** → 补充完整代码（目前只有输出格式）
- **Sub-flow D** → 补充完整代码框架
- **Sub-flow E** → 补充完整数学计算代码

---

## 验证状态

| 项目 | 状态 |
|------|------|
| 模板清理 | ✅ 已验证 |
| Balloon 边界 (4_MS2) | ✅ 已验证 |
| Vector Potential (5_TR/6_CT) | ✅ 已验证 |
| MotionSetup + Setup1 | ✅ 已验证 |
| 绕组公式 GUI 设置 | ✅ 用户已完成 |
| 磁钢材料属性 | ✅ 用户已手动设置 Br/Hc |
| Sub-flow A 代码 | ❌ 待完成 |
| Sub-flow B 代码 | ❌ 待完成 |
| Sub-flow C 代码 | ❌ 待完成 |
| Sub-flow D 代码 | ❌ 待完成 |
| Sub-flow E 代码 | ❌ 待完成 |
