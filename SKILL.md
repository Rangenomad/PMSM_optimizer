---
name: pmsm-optimizer
description: 永磁同步电机(PMSM)电磁设计参数提取与优化。基于固定的Prius 2D模板，用户指定直流母线电压、逆变器电流、电机长度等外部条件后，按管线进行MTPA标定→外特性→效率MAP计算。也支持Ld/Lq、反电势、额定扭矩等单点参数提取。触发词：PMSM、Ld/Lq、反电势、扭矩、效率MAP、外特性、MTPA、弱磁控制、永磁同步电机优化。
allowed-tools: Bash(python -c *) Bash(python << *) Read Write Glob
---

# PMSM Optimizer Skill

基于 Ansys Maxwell2D (PyAEDT) 的永磁同步电机电磁性能评估工具。
**模板固定为 Prius 2D (Prius_2D_Practice.aedt)，不做几何修改。**

---

## 设计哲学

```                         
用户使用条件                    主管线 (A → E → D → F)              输出
─────────────                  ─────────────────────────               ──────
Vdc (母线电压)      ──→  ① Sub-flow A               ──→  ψd/ψq 磁链表
Imax (电流限制)            MTPA 标定                        MTPA 标定表
叠片长度 (可选)            Magnetostatic FEA

                         ② Sub-flow E               ──→  T-n 曲线
                           外特性计算                      工作点轨迹
                           纯解析 (基于磁链表)

                         ③ Sub-flow D               ──→  效率 MAP
                           效率 MAP                       损耗 MAP
                           纯解析 (基于磁链表)             PF/M/Is MAP

                         ④ Sub-flow F               ──→  仿真报告
                           报告生成                        外特性表 CSV
                           纯数据整合                      损耗表 CSV
                                                          Markdown 报告


次要管线 (按需)
─────────────
Sub-flow B: 空载反电势 (Transient FEA)
Sub-flow C: 额定点扭矩 (Transient FEA)
Ld/Lq/Φ 参数提取 (Sub-flow A 的中间产物)
```

**核心思路**：用户只需指定使用条件（Vdc, Imax），主管线自动完成从 FEA 标定到性能评估的全流程。FEA 只跑一次（Sub-flow A），后续外特性和效率 MAP 均为解析计算，秒级出结果。

---

## 用户使用条件

| 参数 | 含义 | 典型值 | 影响 |
|------|------|--------|------|
| **Vdc** | 直流母线电压 (V) | 300, 500, 650 | 决定弱磁拐点、电压极限圆 |
| **Imax** | 逆变器峰值电流 (A) | 200, 300, 350 | 决定最大扭矩、电流极限圆 |
| 叠片长度 | 电机轴向长度 (mm) | 83.82 (可调) | 等比例缩比扭矩/功率 |
| Rs | 相电阻 (Ω) | 0.05 | 影响铜耗、效率绝对值 |

> **注意**：
> - 叠片长度（StackLength）可通过函数参数传入，会自动更新 FEA 模型的 ModelDepth。默认值为 83.82 mm（Prius 基准）。
> - **局部模型**：模板使用 1/8 局部电机模型（6 coils，Master/Slave 反周期边界）。`subflow_a_mtpa_cal.py` 自动将磁链缩放 8 倍（`SYMMETRY_MULTIPLIER=8`）得到完整电机等效值，后续管线无需额外处理。
> - Vdc 和 Imax 也通过函数参数传入，不手动修改 FEA 模型。

---

## 环境要求

| 组件 | 版本 | 说明 |
|------|------|------|
| Python | ≥ 3.10 (推荐 3.13) | |
| AEDT | **2023.2** | 需预先安装 |
| **pyaedt** | **0.25.1** | ⚠️ 严格锁定，见下方 |
| numpy | ≥ 1.20 | |
| scipy | ≥ 1.7 | RegularGridInterpolator |
| matplotlib | ≥ 3.5 | 图表生成 |

### ⚠️ pyaedt 版本

| pyaedt | Python 3.13 | AEDT 2023.2 | 状态 |
|--------|------------|-------------|------|
| **0.25.1** | ✅ | ✅ | **唯一可用** |
| 1.3.0 | ✅ | ❌ gRPC 不兼容 | 不可用 |
| 0.7~0.8 | ❌ | ✅ | 不可用 |

```bash
pip install pyaedt==0.25.1
```

---

## 主管线

### 第零步：创建项目

```python
from scripts.project_utils import create_project
project_dir = create_project('my_design')  # → pmsm_projects/YYYY-MM-DD_my_design/
```

模板 `references/Prius_2D_Practice.aedt` 被复制到项目目录。所有后续操作基于副本。

### ① Sub-flow A：MTPA 标定

**这是唯一需要跑 FEA 的步骤。**

```
求解器: MagnetostaticXY (4_Partial_motor_MS2)
方法:   Id×Iq 网格扫描 → ψd(Id,Iq), ψq(Id,Iq) 磁链表
        基于磁链表做 MTPA 优化 → (n,Is) → (Id_opt, Iq_opt, T_mtpa)
默认网格: 20 Id × 15 Iq = 300 FEA 点, ~15 分钟
扭矩精度: 0.4-1.7% (FW 区 3× 优于 10×10 均匀网格)
```

**用法**：

```python
# 完整 Phase 1 + Phase 2 (默认 FW 加密 20×15 网格)
from scripts.subflow_a_mtpa_cal import run
results = run(project_path='pmsm_projects/my_design', vdc=500, imax=300)

# 仅 Phase 1: ψd/ψq 磁链表 (自定义网格)
from scripts.subflow_a_mtpa_cal import phase1_psi_dq_table
phase1_psi_dq_table(n_id=10, n_iq=10, project_path='pmsm_projects/my_design')

# 仅 Phase 2: MTPA 标定 (基于已有磁链表)
from scripts.subflow_a_mtpa_cal import phase2_mtpa_calibrate
phase2_mtpa_calibrate('pmsm_projects/my_design/psi_dq_table.npz', vdc=500, imax=300)
```

**输出**：
| 文件 | 内容 |
|------|------|
| `psi_dq_table.npz` | id_grid, iq_grid, psi_d(20×15), psi_q(20×15) |
| `mtpa_table.npz` | speed×Is → Thet_opt, T_mtpa, Vs, M, feasible |

### ② Sub-flow E：外特性

```
方法:   纯解析 (不跑 FEA)
输入:   psi_dq_table.npz (来自 ①) + Vdc/Imax
约束:   Vs ≤ Vdc/√3 (电压圆), Is ≤ Imax (电流圆)
区域:   MTPA 恒扭矩区 + FW 弱磁区
输出:   T-n 曲线, Id/Iq/Is/Vs 轨迹, 4-panel 图表
```

**用法**：

```python
from scripts.subflow_e_external import run
results = run(project_path='pmsm_projects/my_design', vdc=500, imax=300)
# → outer_characteristic.npz + outer_char_500V_300A.png
```

**典型输出 (Vdc=500V, Imax=300A)**：

| 指标 | 值 |
|------|-----|
| 峰值扭矩 | 495.5 Nm @ 500 rpm |
| 弱磁拐点 | ~1393 rpm |
| 峰值功率 | 71.8 kW @ 1456 rpm |
| 万转下垂率 | 8.9% |

### ③ Sub-flow D：效率 MAP

```
方法:   纯解析 (不跑 FEA)
输入:   psi_dq_table.npz (来自 ①) + Vdc/Imax
原理:   磁链表插值 + MTPA/FW 控制 → 每 (n,T) 点的最优 (Id,Iq)
损耗:   铜耗 3×(Is/√2)²×Rs + 铁耗 k×(f/f0)^α×(ψ/ψpm)^β
输出:   η/P_loss/Is/M/PF MAP, 4-panel 图表
```

**用法**：

```python
from scripts.subflow_d_efficiency_map import run
results = run(project_path='pmsm_projects/my_design', vdc=500, imax=300,
              n_speed=30, n_torque=30)
# → efficiency_map.npz + efficiency_map_500V_300A.png
```

**典型输出 (Vdc=500V, Imax=300A, 30×30 网格)**：

| 指标 | 值 |
|------|-----|
| 峰值效率 | 98.0% @ 2828 rpm, 85 Nm |
| >97% 高效区 | 44.9% |
| >95% 高效区 | 69.7% |
| 全局最大扭矩 | 495.3 Nm |

### ④ Sub-flow F：报告生成

```
方法:   纯数据整合 (不跑 FEA)
输入:   outer_characteristic.npz (来自 ②) + efficiency_map.npz (来自 ③)
输出:   外特性表 CSV + 损耗表 CSV + Markdown 仿真报告
```

**用法**：

```python
from scripts.subflow_f_report import run
run(project_path='pmsm_projects/my_design')
# → outer_characteristic_table.csv + loss_table.csv + simulation_report.md
```

**输出文件**：

| 文件 | 内容 |
|------|------|
| `outer_characteristic_table.csv` | 转速, 扭矩, 功率 |
| `loss_table.csv` | 转速, 扭矩, 损耗, Udc, Iac_rms, M, PF |
| `simulation_report.md` | 仿真目的、设置、结果汇总、贴图、表链接 |

### 一键运行主管线

```python
# 等价于 ①→②→③→④
from scripts.subflow_a_mtpa_cal import run as run_a
from scripts.subflow_e_external import run as run_e
from scripts.subflow_d_efficiency_map import run as run_d
from scripts.subflow_f_report import run as run_f

project = 'pmsm_projects/my_design'
vdc, imax = 500, 300

run_a(project_path=project, vdc=vdc, imax=imax)        # ① FEA ~15min
run_e(project_path=project, vdc=vdc, imax=imax)        # ② 解析 <1s
run_d(project_path=project, vdc=vdc, imax=imax)        # ③ 解析 ~5s
run_f(project_path=project)                             # ④ 报告 <1s
```

---

## 次要管线

### Sub-flow B：空载反电势

```
求解器: TransientXY (5_Partial_motor_TR)
方法:   Imax=0, 额定转速旋转一个电周期
输出:   三相 BEMF 波形图 (bemf_{speed}rpm.png)
```

```python
from scripts.subflow_b_bemf import run
run(speed=3000, project_path='pmsm_projects/my_design')
```

### Sub-flow C：额定点扭矩

```
求解器: TransientXY (5_Partial_motor_TR)
方法:   固定 (Imax, Thet_deg, Speed_rpm) 求解一个电周期
输出:   扭矩波形 + 平均值 + 脉动率
```

```python
from scripts.subflow_c_torque import run
run(imax=250, thet_deg=45, speed=3000, project_path='pmsm_projects/my_design')
```

### Ld/Lq 参数提取 (旧版 Sub-flow A)

```
求解器: MagnetostaticXY
方法:   电流角扫描法 (与 Id/Iq 网格法不同)
输出:   Ld(Id,Iq), Lq(Id,Iq), 主磁链 Φ
```

```python
from scripts.subflow_a_ldlq import run
r, phi = run(rated_current=250, max_current=1.4, current_steps=5, angle_steps=6,
             project_path='pmsm_projects/my_design')
```

> 这是旧版 Ld/Lq 计算方法（电流角扫描），与主管线使用的 Id/Iq 网格法不同。通常推荐使用主管线的 Sub-flow A 获得 ψd/ψq 磁链表。

---

## 参数权限

各子流程允许修改的 FEA 设计变量：

| 子流程 | 允许修改 | 说明 |
|--------|---------|------|
| A (MTPA 标定) | `Imax`, `StackLength` | Id/Iq 由脚本内部自动控制；StackLength 设置 ModelDepth |
| B (反电势) | `Speed_rpm`, `StackLength` | Imax 自动设为 0 |
| C (额定扭矩) | `Imax`, `Speed_rpm`, `Thet_deg`, `StackLength` | |
| D (效率 MAP) | `StackLength` | 纯解析，StackLength 需与 FEA 表一致（默认 83.82） |
| E (外特性) | `StackLength` | 纯解析，StackLength 需与 FEA 表一致（默认 83.82） |

**全局禁止修改**：`Poles`, `PolePairs`, 其他几何尺寸 (`AirGap`, `MagnetThickness` 等), 材料属性, MotionSetup 初始位置, Master/Slave 边界。

---

## 脚本索引

```
scripts/
├── subflow_a_mtpa_cal.py       # ① 主管线: MTPA 标定 (FEA + 解析)
├── subflow_e_external.py       # ② 主管线: 外特性 (解析)
├── subflow_d_efficiency_map.py # ③ 主管线: 效率 MAP (解析)
├── subflow_f_report.py         # ④ 主管线: 报告生成 (数据整合)
├── subflow_a_ldlq.py           # 次要: Ld/Lq 参数 (FEA, 旧版电流角扫描)
├── subflow_b_bemf.py           # 次要: 反电势 (FEA)
├── subflow_c_torque.py         # 次要: 额定扭矩 (FEA)
├── project_utils.py            # 项目创建/模板拷贝
├── param_guard.py              # 参数权限守卫
└── com_extract.py              # COM 后备提取 (实验性)
```

---

## 自然语言示例

| 用户说 | 执行计划 |
|--------|---------|
| "母线 500V，电流 300A，算全部性能" | ① A → ② E → ③ D → ④ F |
| "生成报告" | ④ F (需 E/D 已完成) |
| "算外特性，母线 300V，电流 200A" | ① A → ② E (Vdc=300, Imax=200) |
| "跑个效率 MAP，500V 300A" | ① A → ③ D (如果 A 已有则跳过) |
| "算空载反电势，3000rpm" | B (speed=3000) |
| "额定点扭矩，250A, 45°, 3000rpm" | C (imax=250, thet=45, speed=3000) |
| "算 Ld/Lq" | A (仅 Phase 1, 或旧版 subflow_a_ldlq) |
| "全部算一遍" | ① A → ② E → ③ D → ④ F → B → C |

---

## 已知问题

1. **GrpcApiError**：FEA 脚本中断后对象状态不同步 → `m2d.modeler.refresh_all_ids()` 或重连
2. **gRPC 重试**：Sub-flow A FEA 含自动重试（最多 3 次，指数退避）
3. **铁耗模型**：Sub-flow D 铁耗为 Steinmetz 解析估计（f 指数 + 磁链指数），不含 PWM 谐波、磁钢涡流损耗，效率绝对值偏高约 1-3%
4. **效率 MAP 无 FEA 验证**：当前效率为纯解析计算，电磁量基于 FEA 验证过的 ψd/ψq 表（扭矩误差 0.4-1.7%），但损耗未经过 FEA 直接验证
5. **COM 后备方案**：`com_extract.py` 为实验性，CreateReport 数组序列化问题待解决
6. **Band 运动带**：Band 必须是闭合环面，几何不正确时 `assign_rotate_motion` 报错
7. **1/8 局部模型**：模板为 48 槽 8 极电机的 1/8 模型（Master/Slave 反周期边界），磁链自动缩放 8×。若更换模板需调整 `SYMMETRY_MULTIPLIER`。

---

## 执行模式

读取 `config.json` 决定 auto/confirm 模式：

- **auto**：展示计划后直接执行
- **confirm**：每步等待用户确认

切换："切换为确认模式" / "切换为自动模式"
