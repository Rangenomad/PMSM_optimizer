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

模板文件 `Prius_2D_Practice.aedt` 包含 2 个设计，从原模板清理后保留：

| 设计名 | 求解器类型 | 用途 | 边界条件 | 绕组 |
|--------|-----------|------|---------|------|
| `4_Partial_motor_MS2` | MagnetostaticXY | Sub-flow A: Ld/Lq MAP | Balloon (Region 外边界) | Current 激励 ×6, Matrix(9匝/线圈) |
| `5_Partial_motor_TR` | TransientXY | Sub-flow B/C/D: 反电势/额定扭矩/效率 MAP | Vector Potential=0 (外弧边), Master/Slave (径向边) | Winding Group + Coil(9匝/线圈) |

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

## Transient 求解器公共配置

所有基于 `5_Partial_motor_TR` (TransientXY) 的子流程使用统一的时间参数配置规则。

### 公式

```
freq        = Speed_rpm / 60 × PolePairs      (Hz, 电频率)
StopTime    = elec_periods / freq              (s, N 个完整电周期)
TimeStep    = 1 / (freq × steps_per_period)    (s, 每电周期固定步数)
```

### 各子流程参数

| 子流程 | elec_periods | steps_per_period | 说明 |
|--------|-------------|-----------------|------|
| B (空载反电势) | 2~3 | 50 | 多周期 FFT 分析 |
| C MTPA 候选角扫描 | 1 | 50 | 快速筛选 |
| C 精确结算 | 3 | 50 | 稳态波形 + 脉动率 |
| D (效率 MAP) | 1 | 50 | 单周期取平均扭矩 |

### 物理含义

- **每工况点固定计算 N 个电周期**，不因转速变化而增加步数。转速高则 StopTime 和 TimeStep 同比例缩短，计算量恒定。
- 电磁转矩纹波频率为 6× 电频率，50 步/周期下每个纹波周期约 8 步，可分辨平均值和粗略脉动率。
- Sub-flow D 只取平均扭矩，elec_periods=1 最小化求解时间。

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

1. Magnetostatic 求解器不支持 Winding Group 激励 → 需要用 Current 边界设定电流

2. **gRPC 兼容**：`assign_current()` 在 gRPC 下不稳定，改用 `boundary.update()` 修改已有 Current 边界对象的 `Current` / `IsPositive` 属性值。

3. **转子位置标定 + 主磁链 Φ 测量（零电流点）**：
   在扫描开始前先求解一个零电流点，从三相磁链中反算转子实际 d 轴位置：
   ```
   ψa, ψb, ψc ← 零电流 Matrix 导出
   ψα = 2/3·(ψa - 0.5·ψb - 0.5·ψc)              (Clarke 变换)
   ψβ = 2/3·(√3/2·ψb - √3/2·ψc)
   θr = atan2(ψβ, ψα)                              (转子 d 轴电角度)
   Φ  = |ψαβ|                                       (永磁磁链幅值)
   ```
   **原因**：转子初始位置不一定是 0°，使用实测 θr 做 Park 变换可避免 d/q 轴分量串扰。

4. **电流变换（Id/Iq → Ia/Ib/Ic → 匝数放大）**（使用实测转子角度 θr）：
   ```
   Iα = Id·cos(θr) - Iq·sin(θr)
   Iβ = Id·sin(θr) + Iq·cos(θr)
   Ia = Iα
   Ib = -0.5·Iα + 0.866·Iβ
   Ic = -0.5·Iα - 0.866·Iβ
   
   每个 coil 对象代表 9 匝 → 设定电流 = 9 × Ia/Ib/Ic
   ```

5. **磁链提取**：因 gRPC 下 `get_solution_data()` 的 `FluxLinkage()` 表达式不稳定，改用 `export_matrix('Matrix1', file)` 导出矩阵结果，从 .txt 文件中解析每线圈的 Flux Linkage 值。

6. 磁链汇总（乘匝数）：
   ```python
   psi_a = 9 * (psi_coil['PhaseA1'] + psi_coil['PhaseA2'])
   psi_b = 9 * (psi_coil['PhaseB1'] + psi_coil['PhaseB2'])
   psi_c = 9 * (psi_coil['PhaseC1'] + psi_coil['PhaseC2'])
   ```

7. Park 变换将三相磁链转为 d/q 轴分量（使用实测转子角度 θr）：
   ```python
   psi_d = 2/3·(ψa·cos(θr) + ψb·cos(θr-120°) + ψc·cos(θr+120°))
   psi_q = -2/3·(ψa·sin(θr) + ψb·sin(θr-120°) + ψc·sin(θr+120°))
   ```

8. **扣除永磁磁链计算 Ld/Lq**（Id, Iq ≠ 0 时）：
   ```
   Ld = (Ψd - Φ) / Id     — 扣除永磁磁链，仅剩电枢反应贡献
   Lq = Ψq / Iq            — PM 在 q 轴无贡献，无需扣除
   ```
   > 文献 [1] 标准做法是 Ld = (Ψd(Id,Iq) - Ψd(0,Iq))/Id 以计入交叉饱和，当前使用固定 Φ = Ψd(0,0) 简化处理 [3]。

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
3. **时间配置**: Transient 公共配置 (elec_periods=2~3, steps_per_period=50)
4. **求解 & 数据提取**（gRPC 兼容 API）：
   - 使用 `get_solution_data_per_variation()`（非 `get_solution_data()`）
   - 表达式：`'InducedVoltage(Phase_A)'`（三相）
   - 数据读取：`data.data_real('InducedVoltage(Phase_A)')` → 返回 **mV**，需 `/1000` 转 V
   - 时间轴：`data.primary_sweep_values` → 返回 **ns**，需 `×1e-9` 转 s
5. **FFT 分析**（Python numpy）：
   - `np.fft.rfft()` 计算幅值谱，`2/N` 幅值校正
   - 基波幅值：DC（索引 0）后最大峰值
   - THD：`√(Σ|Hk|²) / |H1|`，k ≥ 2（排除 DC 和基波）
6. 线反电势常数：`Ke_line = V_fund × √3 / (Speed_rpm/1000)` V/(krpm)

**输入参数**:
```
rated_speed:           额定转速 (rpm)
elec_periods:          仿真电周期数 (默认 2)
time_steps_per_cycle:  每周期步数 (默认 50)
```

**输出**: 三相 BEMF 波形 + 基波幅值 + THD + 线反电势常数 Ke (V/krpm)

---

### Sub-flow C：额定点扭矩

**求解器**: TransientXY (`5_Partial_motor_TR`)
**原理**: 额定负载、额定转速，仿真一个完整电周期。
**时间配置**: Transient 公共配置 (elec_periods=1~3, steps_per_period=50)

**电流角说明**：
```
Thet 是 A 相电流超前 A 相反电势的相位角（电角度）

Id = Imax * sin(Thet)     Thet < 0 → 弱磁方向（Id 负）
Iq = Imax * cos(Thet)     Thet > 0 → 增磁方向（Id 正）
IPM 电机 MTPA 通常在 -20° ~ -45°（负 Id，利用磁阻转矩）
```

**实现要点**:
1. 设置 `Imax` 到额定电流，`Speed_rpm` 到额定转速
2. **MTPA 自动搜索**（用户未指定电流角时）：
   - 以粗网格（50 步/周期）跑 6 个候选角的短仿真（1 个电周期）
   - 候选角：`[0, -15, -25, -35, -45, -60]` 度
   - 每角设置 `m2d['Thet_deg'] = str(angle)`（**无 `°` 后缀**，gRPC 下带 `°` 后缀第二次赋值后静默失败）
   - 提取 Moving1.Torque，取**最后半周期平均**（`torque[-half_period:]`），选扭矩最大的角度
   - 如果 Sub-flow A 的 Ld/Lq/Φ 数据可用，用 MTPA 公式直接算最优角（跳过 FEA 扫描）
3. 以 MTPA 角（或用户指定角）跑完整仿真（3 个电周期，50 步/周期）
4. 取最后一个完整电周期的 `Moving1.Torque` 稳态数据，计算平均值和峰峰值脉动率

**输入参数**:
```
rated_current:  额定电流 (A)
current_angle:  电流角 (°), 可选 — 不指定则自动 MTPA
rated_speed:    额定转速 (rpm)
elec_periods:   仿真电周期数 (默认 3)
```

**输出**: MTPA 最优角 + 平均扭矩 + 扭矩脉动率 + 最大/最小扭矩

---

### Sub-flow D：全域工作特性 MAP

**求解器**: 批量 TransientXY (`5_Partial_motor_TR`)
**原理**: (转速 × 扭矩) 二维网格扫描，每点跑一个电周期。
**时间配置**: Transient 公共配置 (elec_periods=1, steps_per_period=50)
        电流 (Id, Iq) 由 **MTPA 算法** 从目标扭矩生成，而非经验系数估算。

**MTPA 算法** (`_mtpa_for_torque()`):
```
IPM 扭矩方程:  T = 1.5·P·(Φ·Iq + (Ld-Lq)·Id·Iq)
MTPA 轨迹:     Id = A - √(A² + Iq²)   其中 A = Φ / (2·(Lq-Ld))
```
沿 MTPA 轨迹对 Iq 二分搜索匹配目标扭矩，确保每点工作在最优电流角。

**实现要点**:
1. **MTPA 生成 Id/Iq**: `T_target → (Id, Iq, Is, β)` 沿 MTPA 轨迹二分搜索
2. **FEA 扭矩验证**: 以 (Imax=Is, Thet_deg=-β) 设置 FEA 参数并求解，取 `Moving1.Torque` 单值平均值
3. **解析电参数** (基于 MTPA 的 Id/Iq，无需 FEA 波形):
   - Vs = √(Vd²+Vq²), Vd = -ω·Lq·Iq, Vq = ω·(Ld·Id+Φ)
   - 调制比 M = Vs / (Vdc/√3)
   - 功率因数 PF = cos(φv - φi)
   - 铜损 Pcu = 3 × Is²/2 × Rs
4. **gRPC 限制**: PyAEDT 0.25.1 gRPC 在加载 Transient (Imax>0) 下只返回单值时点 → 解析计算替代波形提取
5. 复用同一 m2d 对象加速求解（避免每次开/关工程重创建网格）

**输入参数**:
```
speed_min/max/steps:  转速网格
torque_steps:         扭矩分点数
vdc:                  直流母线电压 (V)
imax:                 最大相电流 (A)
Ld, Lq, Phi, Rs:      电机参数 (可传值覆盖默认)
```

**输出**: 7 张 MAP（T_avg / η / M / PF / Vs / Is / β），每张 CSV 带行列标签 + 终端预览

CSV 文件格式（`results/` 目录）：
```
torque_map.csv       # 列: 扭矩目标值, 行: 转速, 值: FEA 扭矩 (Nm)
efficiency_map.csv   # 列: 扭矩目标值, 行: 转速, 值: 效率 (%)
modulation_map.csv   # 列: 扭矩目标值, 行: 转速, 值: 调制比 M
pf_map.csv           # 列: 扭矩目标值, 行: 转速, 值: 功率因数 PF
voltage_map.csv      # 列: 扭矩目标值, 行: 转速, 值: 相电压 Vs (V)
current_map.csv      # 列: 扭矩目标值, 行: 转速, 值: 电流 Is (A)
beta_map.csv         # 列: 扭矩目标值, 行: 转速, 值: 电流角 β (°)
```
每 CSV 首行为 `#` 注释头，次行为列标签，后续每行为转速行 + CSV 数值。

---

### Sub-flow E：外特性计算

**求解器**: 纯数学计算（不跑 FEA）
**原理**: 在电压圆 Vs ≤ Vdc/√3（调制比 M ≤ 1.0，SVPWM 线性区）和电流圆 Is ≤ Imax 约束下，基于 Sub-flow A 的 Ld/Lq MAP 计算各转速最大输出扭矩。

**调制比说明**：
```
调制比 M = Vs / (Vdc/√3)

线性调制区 (SVPWM):  M_max = 1.0  ← 外特性计算使用此限
过调制(方波):        M_max ≈ 1.10 (不用于外特性计算)
```

**实现要点**:
1. 从 Sub-flow A 获取 Ld(Id,Iq), Lq(Id,Iq), Φ 数据（插值模型）
2. 转速从 0 到用户指定的 n_max，几何分布（低速密高速疏）
3. 对每个转速：
   a. 计算电频率 ω = 2π × n × PolePairs / 60
   b. **MTPA 区**（恒扭矩）：在电流圆上搜索扭矩最大点，检查 Vs ≤ Vdc/√3
   c. **弱磁区**：超出电压限时沿电压圆与电流圆的交点计算
4. 输出 T-n 曲线、转折点、峰值功率

**输入参数**:
```
Ld/Lq/Φ:    来自 Sub-flow A（推荐前置 — 未提供时使用默认近似值）
Vdc:         直流母线电压 (V)
Imax:        最大相电流 (A)
speed_max:   最高转速 (rpm)     ← 外特性截止转速
Rs:          相电阻 (Ω), 可选 — 默认 0（高速时可忽略）
speed_points_n: 转速分点数 (默认 20)
```

**输出**: T-n 曲线表 + 工作点轨迹 (Id, Iq) 表 + 转折点 + 峰值功率

---

## 执行流程

```
用户指令 → [新增] Step -1: 需求澄清
         → [新增] Step 0: 创建项目文件夹，拷贝模板到项目目录
         → [新增] Step 0a: 确认模式 (auto/confirm)
         → Step 1: 打开项目模板（项目文件夹中的副本）
         → Step 2: 修改参数 (仅限白名单内参数)
         → Step 3: 判定需求 → 映射子流程
         → 执行子流程 (A/B/C/D/E)
         → 输出结果
```

### 项目隔离机制（Step -1 & Step 0）

为防止模板文件被污染，加入项目文件夹机制：

- **Step -1 需求澄清**：与用户确认计算目标和参数范围
- **Step 0 创建项目**：
  1. 在 `pmsm_projects/` 下以 `YYYY-MM-DD_<描述>` 格式创建项目文件夹
  2. 使用 `shutil.copy2()` 将 `references/Prius_2D_Practice.aedt` 复制到项目目录
  3. 所有后续操作基于副本进行，原始模板不被修改
  4. 项目文件夹内生成 `project.json`（创建时间、参数配置）

### 参数权限管控（Step 2 约束）

AI 只能调节与当前仿真任务相关的必要参数，禁止修改全局参数。

**各子流程允许参数：**

| 子流程 | 允许修改 | 说明 |
|--------|---------|------|
| A (Ld/Lq MAP) | `Imax` | Thet_deg 由脚本内部 Id/Iq→abc 控制 |
| B (反电势) | `Speed_rpm` | Imax 脚本自动设为 0 |
| C (额定扭矩) | `Imax`, `Speed_rpm`, `Thet_deg` | |
| D (效率 MAP) | `Speed_rpm`, `Imax`, `Thet_deg` | 脚本内部 MTPA 控制 Id/Iq |
| E (外特性) | 无（纯数学计算） | 全部通过函数参数传入 |

**全局禁止修改（所有子流程）：**
`Poles`, `PolePairs`, 几何尺寸, 材料属性, MotionSetup 初始位置, Master/Slave 边界

**实现机制：**
1. SKILL.md 中定义参数白名单，AI 在生成脚本前即可知悉
2. 脚本中增加 `set_var()` 守卫函数，运行时验证赋值目标
3. 用户坚持修改禁止参数时 → 建议在 GUI 中手动操作

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
2. 清理多余设计（保留 2 个核心设计 + 备份）
3. 4_MS2 → 设置 Balloon 边界 + Current 激励 + Matrix
4. 5_TR → 设置 MotionSetup (Speed_rpm 变量引用) + Transient Setup + Vector Potential

> 6_Partial_motor_CT 已移除 — 所有 Transient 子流程 (B/C/D) 统一使用 5_TR。

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
| `assign_current()` | ⚠️ gRPC 不稳定 | 改用 `boundary.update()` 修改已有 Current 边界（Sub-flow A） |
| `assign_winding()` | ✅ 可用 | 用于 Transient Winding Group |
| `assign_coil()` | ✅ 可用 | 用于指定线圈对象 |
| `m2d[...]=value` | ✅ 可用 | 设计变量读写 |
| `assign_balloon()` | ✅ 可用 | Magnetostatic 边界 |
| `assign_vector_potential()` | ✅ 可用 | Transient 边界 |
| `analyze()` | ✅ 可用 | 求解控制 |
| `post.get_solution_data()` | ⚠️ gRPC 单值限制 | 推荐 `get_solution_data_per_variation()` + `data.data_real()` + `data.primary_sweep_values` |
| `export_matrix()` | ✅ 可用 | 替代 FluxLinkage 表达式（Sub-flow A 磁链提取） |

---

## Skill.md 和 pmsm-methods.md

两个文件已完成更新：

### SKILL.md
- 模板文件名 → `Prius_2D_Practice.aedt`
- 参数表 → 实际模板变量（Imax/Thet_deg/Speed_rpm/PolePairs）
- 子流程引用 → `scripts/` 目录

### pmsm-methods.md
包含所有子流程的关键公式和模块引用：
- **Sub-flow A** → boundary.update() + 9 匝因子 + 转子位置标定 + Park 变换
- **Sub-flow B** → Imax=0 + InducedVoltage + FFT + THD
- **Sub-flow C** → MTPA 自动搜索 + 扭矩波形提取
- **Sub-flow D** → MTPA 算法 + 解析电参数 + 7 MAP CSV
- **Sub-flow E** → 纯数学外特性计算（Vdc/Imax/speed_max）

---

## 验证状态

| 项目 | 状态 |
|------|------|
| 模板清理 | ✅ 已验证 |
| Balloon 边界 (4_MS2) | ✅ 已验证 |
| Vector Potential (5_TR) | ✅ 已验证 |
| MotionSetup + Setup1 | ✅ 已验证 |
| 绕组公式 GUI 设置 | ✅ 用户已完成 |
| 磁钢材料属性 | ✅ 用户已手动设置 Br/Hc |
| Sub-flow A 代码 | ✅ 已完成 (scripts/subflow_a_ldlq.py) |
| Sub-flow B 代码 | ✅ 已完成 (scripts/subflow_b_bemf.py) |
| Sub-flow C 代码 | ✅ 已完成 (scripts/subflow_c_torque.py) |
| Sub-flow D 代码 | ✅ 已完成 (scripts/subflow_d_efficiency_map.py) |
| Sub-flow E 代码 | ✅ 已完成 (scripts/subflow_e_external.py) |
| Spec ↔ 代码同步 | ✅ 已完成 (2026-07-19) |
| Plan ↔ 代码同步 | ✅ 已完成 (2026-07-19) |
