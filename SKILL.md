---
name: pmsm-optimizer
description: 永磁同步电机(PMSM)电磁设计参数提取与优化。基于固定的模板模型，通过修改参数和运行预设子流程，计算Ld/Lq MAP、反电势、扭矩、全域工作特性MAP、外特性等。当用户提到PMSM参数计算、Ld/Lq、交直轴电感、反电势、扭矩分析、损耗MAP、效率MAP、外特性、T-n曲线、MTPA、弱磁控制、永磁同步电机优化时触发。
allowed-tools: Bash(python -c *) Bash(python << *) Read Write Glob
---

# PMSM Optimizer Skill

基于 Ansys Maxwell2D (PyAEDT) 的永磁同步电机电磁参数提取与优化工具。

## 与 maxwell2d-controller 的关系

本 skill 继承其设计原则（分步执行、auto/confirm 模式、跨步骤会话保持、代码规范），但专注于 PMSM 电磁设计这一个场景，不做通用电磁仿真。

---

## 第零步：需求澄清 + 创建项目

### 0a: 需求澄清

与用户确认计算目标和参数范围（算哪个子流程、电流/转速范围等），澄清后继续。

### 0b: 创建项目文件夹

在 `pmsm_projects/` 下创建独立项目文件夹（如 `pmsm_projects/2026-07-19_ldlq_scan/`），将模板文件 `references/Prius_2D_Practice.aedt` 复制到项目目录中。**所有后续操作基于项目目录中的副本，不修改原始模板文件。**

### 0c: 确认执行模式

读取 `config.json`：

- 文件存在且含 `"execution_mode"` → 直接使用，不再询问
- 文件不存在 → 向用户展示 A/C 选项，保存后继续

**模式说明：**
- **auto 模式**：计划展示后直接执行，不等待用户确认
- **confirm 模式**：每步展示代码，等待用户确认后执行

模式可随时通过"切换为确认模式"/"切换为自动模式"切换。

---

## 第一步：连接项目模板

打开项目目录中的模板副本（非原始模板）。

### 关键参数

| 参数 | 值 | 说明 |
|---|---|---|
| 模板文件 | `references/Prius_2D_Practice.aedt` | 用户预先配置好的 PMSM 模型 |
| 工作模式 | `non_graphical=False` | 图形模式，用户可实时查看 |
| 会话模式 | `new_desktop=False, close_on_exit=False` | 复用已有 AEDT 会话 |
| 求解器选择 | 根据子流程自动选择 | Sub-flow A → `4_Partial_motor_MS2` (Magnetostatic) |
|  |  | Sub-flow B/C/D → `5_Partial_motor_TR` (Transient) |

### 执行代码

```python
import sys
sys.stdout.reconfigure(encoding='utf-8')
from pathlib import Path
from ansys.aedt.core import Maxwell2d

# 项目目录（由 Step 0b 创建，包含模板副本）
PROJECT_DIR = Path('pmsm_projects/2026-07-19_ldlq_scan')
TEMPLATE = str(PROJECT_DIR / 'Prius_2D_Practice.aedt')

# Sub-flow A: Magnetostatic
m2d = Maxwell2d(
    project=TEMPLATE,
    design="4_Partial_motor_MS2",
    solution_type="MagnetostaticXY",
    non_graphical=False,
    new_desktop=False,
    close_on_exit=False
)

# 或 Sub-flow B/C/D: Transient
# m2d = Maxwell2d(
#     project=TEMPLATE,
#     design="5_Partial_motor_TR",
#     solution_type="TransientXY",
#     ...
# )

print(f"Project: {m2d.project_name}")
print(f"Design:  {m2d.design_name}")
print("模板项目加载完成")
```

---

## 第二步：按指令修改参数（可选）

用户要求修改参数时执行此步。以下为可修改参数及其对应操作：

| 参数 | 变量名 | 修改方式 |
|------|--------|---------|
| 相电流幅值 | `Imax` | `m2d['Imax'] = '250A'` |
| 电流角 | `Thet_deg` | `m2d['Thet_deg'] = '-30°'` |
| 机械转速 | `Speed_rpm` | `m2d['Speed_rpm'] = '3000rpm'` |
| 极对数 | `PolePairs` | 固定值 (Poles/2=4) |

**提示**：完整设计变量列表见设计文档 `docs/superpowers/specs/2026-07-12-pmsm-optimizer-design.md`。修改前先用 `m2d['VariableName']` 验证值。

### 参数权限管控

AI **只能**调节与当前仿真任务相关的必要参数。修改参数前先确认当前子流程，然后只操作白名单内的参数。

**各子流程允许修改的参数：**

| 子流程 | 允许修改 | 说明 |
|--------|---------|------|
| A (Ld/Lq MAP) | `Imax` | Thet_deg 由脚本内部 Id/Iq→abc 自动控制 |
| B (反电势) | `Speed_rpm` | Imax 脚本自动设为 0 |
| C (额定扭矩) | `Imax`, `Speed_rpm`, `Thet_deg` | |
| D (效率 MAP) | `Speed_rpm`, `Imax`, `Thet_deg` | 脚本内部 MTPA 控制 Id/Iq |
| E (外特性) | 无（纯数学计算） | 全部通过函数参数传入 |

**全局禁止修改（所有子流程）：** `Poles`, `PolePairs`, 几何尺寸, 材料属性, MotionSetup 初始位置, Master/Slave 边界

> 如果用户要求修改禁止参数，说明原因并建议用户在 GUI 中手动修改。不可绕过此限制。

---

## 第三步：判定用户需求 → 生成子流程计划

根据用户的指令判定需要执行的子流程，生成检视列表。

### 用户指令 → 子流程映射

| 用户说 | 执行计划 |
|---|---|
| "算 Ld/Lq" | Step1 → Sub-flow A |
| "算反电势" | Step1 → Sub-flow B |
| "算额定扭矩" | Step1 → Sub-flow C |
| "跑个损耗 MAP + 工作特性" | Step1 → Sub-flow D |
| "算外特性，母线 300V，电流 200A" | Step1 → Sub-flow A → Sub-flow E |
| "改长度到 80mm，再算反电势和扭矩" | Step1 → Step2 → Sub-flow B → Sub-flow C |
| "全部算一遍" | Step1 → Sub-flow A → B → C → D → E |
| "算外特性，需要 Sub-flow A 的数据" | Step1 → Sub-flow A → Sub-flow E |

### 子流程间依赖关系

```
Sub-flow A (Ld/Lq MAP)     ← 独立，但 Sub-flow D/E 依赖其结果
Sub-flow B (反电势)        ← 独立
Sub-flow C (额定扭矩)      ← 独立
Sub-flow D (工作特性 MAP)  ← 可重用 A 的结果优化工作点反算
Sub-flow E (外特性)        ← 依赖 A 的输出
```

### 计划展示格式（auto 模式示例）

```
PMSM 计算计划：
[x] Step 1: 打开模板项目
[>] Sub-flow A: Ld/Lq MAP + 主磁链 Φ
[ ] Sub-flow C: 额定点扭矩
[ ] Sub-flow E: 外特性计算（Vdc=300V, Imax=200A）
```

---

## 子流程执行规范

所有子流程的代码位于 `scripts/` 目录，方法说明见 `references/pmsm-methods.md`。
AI 根据用户指令选择模块调用。详见设计文档 `docs/superpowers/specs/2026-07-12-pmsm-optimizer-design.md`。

### 通用规范

1. **每个脚本独立运行**，通过 `new_desktop=False` 保持会话连接
2. **每脚本首行**：`sys.stdout.reconfigure(encoding='utf-8')`
3. **heredoc 分隔符**：`<< 'EOF'`（带引号，防止 shell 展开 `$`）
4. **不常用 API 前先查签名**：`python -c "from ansys.aedt.core import Maxwell2d; help(Maxwell2d.<method>)"`
5. **每步验证**：输出关键结果确认合理后再继续

### Sub-flow A：Ld/Lq MAP + 主磁链 Φ

- 求解器：`MagnetostaticXY`
- 方法：Id/Iq 矩阵扫描 + 磁链法
- 输出：Ld(Id,Iq) / Lq(Id,Iq) MAP + 主磁链 Φ
- 模板参考：`references/pmsm-methods.md` §A

### Sub-flow B：空载反电势

- 求解器：`TransientXY`
- 方法：空载额定转速旋转
- 输出：三相 BEMF 波形 + FFT 谐波幅值 + THD
- 模板参考：`references/pmsm-methods.md` §B

### Sub-flow C：额定点扭矩

- 求解器：`TransientXY`
- 方法：额定负载、额定转速
- 输出：扭矩波形 + 平均值 + 脉动率
- 模板参考：`references/pmsm-methods.md` §C

### Sub-flow D：全域工作特性 MAP

- 求解器：`TransientXY`（批量）
- 方法：T-n 网格扫描，每点一个电周期
- 输出：损耗 MAP / Is MAP / 功率因数 MAP / 调制比 MAP
- 模板参考：`references/pmsm-methods.md` §D

### Sub-flow E：外特性计算

- 方法：纯数学计算（不跑 FEA）
- 输入：Sub-flow A 结果 + 用户输入（Vdc/Imax/Rs）
- 约束：电压圆 Vs ≤ Vdc/√3，电流圆 Is ≤ Imax
- 区域：MTPA 恒扭矩区 + 弱磁区
- 输出：T-n 曲线 + 工作点轨迹
- 模板参考：`references/pmsm-methods.md` §E

---

## 已知问题（继承 + PMSM 补充）

### 继承自 maxwell2d-controller

1. **GrpcApiError**：脚本中断后对象状态不同步 → `m2d.modeler.refresh_all_ids()` 或断开重连
2. **API 参数名版本差异**：不常用 API 前先用 `help()` 确认签名
3. **批量复制后材料设置失败**：使用 `assign_material()` 替代直接属性赋值
4. **模块导入路径错误**：使用正确的 `ansys.aedt.core.modeler.modeler_2d` 等路径

### PMSM 特有

5. **冻结磁导率 API**：不同 PyAEDT 版本中冻结磁导率的 API 差异较大，使用前用 `help()` 确认
6. **Band 运动带**：Band 必须是闭合环面，几何不正确时 `assign_rotate_motion` 报错
7. **批量仿真资源管理**：Sub-flow D 可能跑几百个工况，建议：
   - 先跑小网格验证（3×3 点），再跑完整网格
   - 控制同时打开的 AEDT 实例数
   - 每点仿真后检查收敛状态
8. **电流角扫描的周期边界**：转子旋转后需确保 Master/Slave 边界对齐

---

## 执行模式切换

用户可随时切换执行模式：

- "切换为确认模式" / "confirm mode" → config.json 中设为 `"confirm"`，后续每步等待确认
- "切换为自动模式" / "auto mode" → config.json 中设为 `"auto"`，后续自动执行
- "当前是什么模式" → 读取 config.json 告知用户

---

## 自然语言示例

以下输入应被理解并正确映射到子流程：

```
"帮我算 PMSM 的 Ld/Lq MAP，额定电流 12A，最大电流 36A"
→ Step1 → Sub-flow A (rated_current=12, max_current=3x)

"算空载反电势，转速 4000rpm"
→ Step1 → Sub-flow B (rated_speed=4000)

"轴向长度改成 85mm，算额定扭矩和反电势"
→ Step1 → Step2(改长度) → Sub-flow C → Sub-flow B

"跑个外特性，母线 250V，峰值电流 300A，电阻 0.008Ω"
→ Step1 → Sub-flow A → Sub-flow E (Vdc=250, Imax=300, Rs=0.008)

"算损耗 MAP，转速 500-8000，10 个点"
→ Step1 → Sub-flow D (speed_min=500, speed_max=8000, speed_steps=10)
```
