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

## 第零步：确认执行模式

完全沿用 maxwell2d-controller 的机制。读取 `~/.claude/skills/pmsm-optimizer/config.json`：

- 文件存在且含 `"execution_mode"` → 直接使用，不再询问
- 文件不存在 → 向用户展示 A/C 选项，保存后继续

**模式说明：**
- **auto 模式**：计划展示后直接执行，不等待用户确认
- **confirm 模式**：每步展示代码，等待用户确认后执行

模式可随时通过"切换为确认模式"/"切换为自动模式"切换。

---

## 第一步：连接模板项目

将模板文件从 skill 目录复制到当前工作目录并打开。

### 关键参数

| 参数 | 值 | 说明 |
|---|---|---|
| 模板文件 | `references/pmsm_template.aedt` | 用户预先配置好的 PMSM 模型 |
| 工作模式 | `non_graphical=False` | 图形模式，用户可实时查看 |
| 会话模式 | `new_desktop=False, close_on_exit=False` | 复用已有 AEDT 会话 |
| 求解器选择 | 根据子流程自动选择 | Magnetostatic / Transient |

### 执行代码

```bash
# 复制模板到工作目录
cp references/pmsm_template.aedt ./pmsm_working.aedt

python << 'EOF'
import sys
sys.stdout.reconfigure(encoding='utf-8')
from ansys.aedt.core import Maxwell2d
import subprocess

m2d = Maxwell2d(
    project="pmsm_working",
    design="Maxwell2DDesign1",
    solution_type="MagnetostaticXY",   # 先以 Magnetostatic 打开，后续切换
    non_graphical=False,
    new_desktop=False,
    close_on_exit=False
)

# 将 AEDT 窗口置前
try:
    subprocess.run(
        ['powershell', '-Command',
         '$wshell = New-Object -ComObject wscript.shell; $wshell.AppActivate("Ansys Electronics Desktop")'],
        capture_output=True, timeout=5
    )
except Exception:
    pass

print(f"Project: {m2d.project_name}")
print(f"Design:  {m2d.design_name}")
print("模板项目加载完成")
EOF
```

---

## 第二步：按指令修改参数（可选）

用户要求修改参数时执行此步。以下为可修改参数及其对应操作：

| 参数 | 用户指令示例 | 操作代码 |
|---|---|---|
| 轴向长度 | "改到 80mm" | `m2d.model_depth = "80mm"` |
| 每槽匝数 | "每槽 25 匝" | `m2d["turns_per_slot"] = "25"` |
| 磁钢牌号 | "换成 N42SH" | 修改材料库中的磁钢属性 |
| 硅钢牌号 | "定子换 B35AV1900" | 修改材料库中的硅钢属性 |
| 额定电流 | "额定电流 15A" | 修改电流激励幅值 |
| 额定转速 | "转速 3000rpm" | `m2d.set_initial_angle()` + MotionSetup 转速修改 |
| 气隙长度 | "气隙改到 0.8mm" | 移动转子几何 |
| 磁钢厚度 | "磁钢 5mm" | 修改磁钢几何尺寸 |
| 绕组连接方式 | "Y 接" / "Δ 接" | 修改绕组设置 |
| 铜线直径 | "线径 1.2mm" | 修改绕组属性 |

**提示**：模板中的参数名（设计变量名）因人而异，AI 在执行修改前应先用 `m2d.variable_manager.design_variable_names` 列出可用变量让用户确认。

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

所有子流程的代码模板位于 `references/pmsm-methods.md`。AI 对应选取代码块，替换 `{{placeholder}}` 后执行。

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
