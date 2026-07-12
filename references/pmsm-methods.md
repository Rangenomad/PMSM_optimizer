# PMSM Optimizer — 子流程计算模板

本文档包含 5 个子流程的完整可运行 Python 代码模板。每个模板用 `{{placeholder}}` 标记可变参数，AI 根据用户指令替换后执行。

---

## Sub-flow A：Ld/Lq MAP + 主磁链 Φ

### 方法说明

- **求解器**：MagnetostaticXY
- **方法**：电流角扫描法。在每个 (Id, Iq) 工作点求解静磁场，从磁链中提取 d/q 轴分量，计算 Ld = Ψd/Id，Lq = Ψq/Iq。
- **适用场景**：需要电感饱和特性时（效率 MAP、外特性计算的前置步骤）

### 可变参数

| 参数 | 说明 | 默认值 |
|---|---|---|
| `{{rated_current}}` | 额定电流幅值 (A) | 10 |
| `{{max_current}}` | 最大电流倍率 (× 额定) | 3.0 |
| `{{current_steps}}` | 电流幅值分点数 | 6 |
| `{{angle_steps}}` | 电流角分点数 (0-90°) | 7 |
| `{{project_name}}` | 打开的 AEDT 项目名 | pmsm_template |
| `{{design_name}}` | 设计名称 | Maxwell2DDesign1 |

### 输出格式

```
=== Ld/Lq MAP ===
Id\Iq │ -30A   -20A   -10A     0A    10A    20A    30A
──────┼────────────────────────────────────────────────
-30A  │ Ld=xx  Ld=xx  Ld=xx  Ld=xx  Ld=xx  Ld=xx  Ld=xx
      │ Lq=xx  Lq=xx  Lq=xx  Lq=xx  Lq=xx  Lq=xx  Lq=xx
...
主磁链 Φ = x.xxx Wb
```

### 代码框架

```python
import sys
sys.stdout.reconfigure(encoding='utf-8')
from ansys.aedt.core import Maxwell2d
import numpy as np

m2d = Maxwell2d(
    project="{{project_name}}",
    design="{{design_name}}",
    solution_type="MagnetostaticXY",
    non_graphical=False,
    new_desktop=False,
    close_on_exit=False
)

I_rated = {{rated_current}}          # 额定电流
I_max = {{max_current}} * I_rated    # 最大电流
n_I = {{current_steps}}              # 电流步数
n_theta = {{angle_steps}}            # 角度步数

currents = np.linspace(I_max / n_I, I_max, n_I)
angles = np.linspace(0, 90, n_theta)

# === 电流角扫描循环 ===
results = []
for Ia in currents:
    for theta_deg in angles:
        theta = np.radians(theta_deg)
        Id = Ia * np.sin(theta)
        Iq = Ia * np.cos(theta)

        # 设置 d/q 轴电流激励（具体 API 视模板中的激励设置而定）
        # TODO: 补充具体绕组电流设置代码
        # m2d[...] = Id  # 设置 d 轴电流分量
        # m2d[...] = Iq  # 设置 q 轴电流分量

        # 求解
        m2d.analyze("Setup1")

        # 提取 d/q 轴磁链
        # TODO: 补充磁链提取代码
        # psi_d = ...
        # psi_q = ...

        # 计算电感（Id ≠ 0 时）
        Ld = psi_d / Id if abs(Id) > 1e-6 else 0
        Lq = psi_q / Iq if abs(Iq) > 1e-6 else 0

        results.append({
            "Id": Id, "Iq": Iq,
            "Ld": Ld, "Lq": Lq,
            "Psi_d": psi_d, "Psi_q": psi_q
        })

# 主磁链：Id=0, Iq=0 时的 q 轴磁链（或永磁体单独贡献）
phi = ...  # TODO: 补充永磁体磁链提取

# 输出结果表格
print(f"主磁链 Φ = {phi:.6f} Wb")
print("=== Ld/Lq MAP ===")
# 格式化输出 Id×Iq 矩阵
...
```

---

## Sub-flow B：空载反电势

### 方法说明

- **求解器**：TransientXY
- **条件**：电枢绕组开路（电流激励设为 0），永磁体激励，额定转速机械旋转
- **输出**：三相 BEMF 波形 + FFT 谐波分析 + THD

### 可变参数

| 参数 | 说明 | 默认值 |
|---|---|---|
| `{{rated_speed}}` | 额定转速 (rpm) | 3000 |
| `{{elec_periods}}` | 仿真电周期数 | 2 |
| `{{time_steps_per_cycle}}` | 每电周期步数 | 200 |

### 输出格式

```
=== 空载反电势结果 ===
额定转速: 3000 rpm
频率: 200 Hz

A 相: 基波幅值=xxx V, THD=x.x%
B 相: 基波幅值=xxx V, THD=x.x%
C 相: 基波幅值=xxx V, THD=x.x%

线反电势常数 Ke = xxx V/(krpm)
```

### 代码框架

```python
import sys
sys.stdout.reconfigure(encoding='utf-8')
from ansys.aedt.core import Maxwell2d
import numpy as np

m2d = Maxwell2d(
    project="{{project_name}}",
    design="{{design_name}}",
    solution_type="TransientXY",
    non_graphical=False,
    new_desktop=False,
    close_on_exit=False
)

speed = {{rated_speed}}              # rpm
n_periods = {{elec_periods}}         # 电周期数
steps_per = {{time_steps_per_cycle}} # 每周期步数

# 计算仿真时间
# poles = m2d....  # 获取极对数 TODO
# freq = speed / 60 * poles
# stop_time = n_periods / freq
# time_step = 1 / (freq * steps_per)

# 设置空载（电流源=0）
# TODO: 设置三相电流为 0

# 设置转速
# TODO: 修改 MotionSetup 的转速

# 配置求解
# setup = m2d.create_setup()
# setup.props["StopTime"] = f"{stop_time}s"
# setup.props["TimeStep"] = f"{time_step}s"

m2d.analyze("Setup1")

# 提取三相电压
# TODO: 从后处理提取 Voltage(A), Voltage(B), Voltage(C) 波形

# FFT 谐波分析
# TODO: 对各相电压做 FFT，计算基波幅值和 THD
```

---

## Sub-flow C：额定点扭矩

### 方法说明

- **求解器**：TransientXY
- **条件**：额定负载（Id/Iq 根据 MTPA 或用户给定）、额定转速
- **输出**：扭矩波形 + 平均值 + 脉动率

### 可变参数

| 参数 | 说明 | 默认值 |
|---|---|---|
| `{{rated_current}}` | 额定电流 (A) | 10 |
| `{{current_angle}}` | 电流角 (°) | 0 (MTPA) |
| `{{rated_speed}}` | 额定转速 (rpm) | 3000 |
| `{{elec_periods}}` | 仿真电周期数 | 3 |
| `{{load_type}}` | 负载类型 | rated |

### 输出格式

```
=== 额定点扭矩结果 ===
平均扭矩: xx.x Nm
扭矩脉动: x.x% (峰峰值)
最大扭矩: xx.x Nm
最小扭矩: xx.x Nm
电流角: xx°
```

---

## Sub-flow D：全域工作特性 MAP

### 方法说明

- **求解器**：批量 TransientXY
- **网格**：(转速 × 扭矩) 二维网格，每点跑一个电周期
- **输出**：4 张 MAP（同一组计算数据的不同视角）

### 输入网格

| 参数 | 说明 | 默认值 |
|---|---|---|
| `{{speed_min}}` | 最低转速 (rpm) | 500 |
| `{{speed_max}}` | 最高转速 (rpm) | 10000 |
| `{{speed_steps}}` | 转速分点数 | 10 |
| `{{torque_steps}}` | 扭矩分点数 | 10 |
| `{{vdc}}` | 直流母线电压 (V) | 300 |
| `{{imax}}` | 最大相电流 (A) | 200 |

### 输出格式

每张 MAP 输出为一个 CSV 矩阵加一个终端表格预览：

```
[1/4] 损耗 MAP P_loss(n, T)        —  单位 W
[2/4] 电流幅值 MAP Is(n, T)         —  单位 A
[3/4] 功率因数 MAP cosφ(n, T)       —  无量纲
[4/4] 调制比 MAP M(n, T)            —  无量纲（>1 为过调制区）
```

### 代码框架

```python
import sys
sys.stdout.reconfigure(encoding='utf-8')
from ansys.aedt.core import Maxwell2d
import numpy as np

# 参数
speed_points = np.linspace({{speed_min}}, {{speed_max}}, {{speed_steps}})
torque_points = np.linspace(...)  # 根据额定扭矩推导
Vdc = {{vdc}}
Imax = {{imax}}

results = np.zeros((len(speed_points), len(torque_points), 4))
# [:, :, 0] = 损耗, [:, :, 1] = Is, [:, :, 2] = PF, [:, :, 3] = M

for i, n in enumerate(speed_points):
    for j, T in enumerate(torque_points):
        # 根据 Ld/Lq MAP (Sub-flow A) 反算 Id/Iq 工作点
        # 设置激励、转速 → 求解一个电周期
        # 提取：相电流幅值 Is、铜损 Pcu、铁损 Pfe、电压幅值 Vs、功率角

        P_loss = Pcu + Pfe
        PF = ...        # 功率因数
        M = Vs / (Vdc / np.sqrt(3))  # 调制比

        results[i, j] = [P_loss, Is, PF, M]

# 输出 4 张 MAP
print("=== 损耗 MAP (W) ===")
# ... 格式化矩阵 ...

print("=== Is MAP (A) ===")
# ...

print("=== 功率因数 MAP ===")
# ...

print("=== 调制比 MAP ===")
# ...
```

---

## Sub-flow E：外特性计算

### 方法说明

- **不跑 FEA**，纯数学计算
- **原理**：在电压圆（Vs ≤ Vdc/√3）和电流圆（Is ≤ Imax）约束下，基于 Ld/Lq MAP + Φ 计算各转速的最大输出扭矩
- **包含**：恒扭矩区（MTPA）+ 弱磁区
- 详见 `references/demagnetization.md`（如有）

### 输入

| 参数 | 说明 | 来源 |
|---|---|---|
| Ld(Id,Iq), Lq(Id,Iq), Φ | 电感与主磁链 | Sub-flow A 输出 |
| Vdc | 直流母线电压 (V) | 用户输入 |
| Imax | 最大相电流 (A) | 用户输入 |
| Rs | 相电阻 (Ω) | 模板预设或用户输入 |
| `{{speed_points_n}}` | 外特性转速分点数 | 20 |

### 输出格式

```
=== 外特性 (T-n 曲线) ===
Vdc = 300 V, Imax = 200 A, Rs = 0.01 Ω

转速(rpm) │ 扭矩(Nm) │ 功率(kW) │ Id(A) │ Iq(A) │ 调制比 │ 区域
──────────┼──────────┼──────────┼───────┼───────┼────────┼──────
    500   │    xx    │   xx     │  xx   │  xx   │  0.xx  │ MTPA
   1000   │    xx    │   xx     │  xx   │  xx   │  0.xx  │ MTPA
   ...
   6000   │    xx    │   xx     │  xx   │  xx   │  0.xx  │ 弱磁
   ...

额定转速: xxxx rpm (转折点)
最高转速: xxxx rpm
峰值功率: xx.x kW @ xxxx rpm
```
