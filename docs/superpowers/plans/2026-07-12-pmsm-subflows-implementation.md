# PMSM Optimizer Sub-flows Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Create standalone sub-flow script modules (`scripts/`) implementing all 5 PMSM sub-flows, update `references/pmsm-methods.md` and `SKILL.md` to reference them, and verify each against the AEDT template.

> **Updated (2026-07-19):** Plan has been synchronized with actual implementation. Key differences from original plan: step counts 200→50 (all sub-flows); Sub-flow A uses `boundary.update()` (gRPC compatible) and `export_matrix()` for flux linkage; Sub-flows B/C/D use `get_solution_data_per_variation()` + `data.data_real()` + `primary_sweep_values`; Thet_deg set without `°` suffix; pole_pairs from `float(m2d['Poles']) / 2`; Sub-flow D expanded to 7 MAP CSVs with full MTPA+analytical engine. See actual files in `scripts/` for current version.

**Architecture:** Template-formula approach — one-time GUI winding setup, then fully automated via design variable changes. Sub-flow A (Magnetostatic) uses `boundary.update()` on existing Current boundaries with 9× turns factor (gRPC compatible). Sub-flows B/C/D (Transient) modify design variables `Imax`, `Speed_rpm`, `Thet_deg`. Sub-flow E (pure math) uses interpolation on Sub-flow A results.

**Tech Stack:** PyAEDT 0.25.1 (gRPC), Ansys Maxwell 2D 2023.2, Python 3.13, Windows 11.

## Global Constraints

- All scripts use `sys.stdout.reconfigure(encoding='utf-8')` as first line
- All scripts open with `non_graphical=False, new_desktop=False, close_on_exit=False` (graphical mode, reuse AEDT session)
- Template file: `references/Prius_2D_Practice.aedt` (absolute path relative to skill directory)
- Design variables: `Imax` (current amplitude), `Speed_rpm` (speed), `Thet_deg` (current angle), `Omega_rad` (electrical angular velocity), `PolePairs`/`Poles` (pole pairs — use `float(m2d['Poles']) / 2` to read)
- Winding formula (Transient): `Imax*sin(Omega_rad*time + Thet)` (Phase A), `-2*pi/3` (Phase B), `+2*pi/3` (Phase C)
- Sub-flow A (Magnetostatic): uses `boundary.update()` (not `assign_current()`) on existing Current boundaries, multiplies by turns factor: **9 × Ia** per coil object. Flux linkage via `export_matrix('Matrix1', file)` parsing.
- Sub-flows B/C/D (Transient): use `get_solution_data_per_variation()` + `data.data_real()` + `data.primary_sweep_values` (not `get_solution_data()` + `data.data()` — gRPC requirement)
- Sub-flow B `InducedVoltage` data returned in mV, divide by 1000; time in ns, multiply by 1e-9
- `Thet_deg` assigned as `str(angle)` without `°` suffix — gRPC fails silently on second `°`-suffixed assignment
- PyAEDT gRPC limitations: `ChangeProperty`, `SetPropertyValue`, `DeleteBoundary`, `GetPropertyValue` all fail — use `m2d[...]=value` for variables

---

### Task 1: Project structure — add sub-flow script modules

**Files:**
- Create: `C:\Users\Wesley\.claude\skills\pmsm-optimizer\scripts\__init__.py`
- Create: `C:\Users\Wesley\.claude\skills\pmsm-optimizer\scripts\subflow_a_ldlq.py`
- Create: `C:\Users\Wesley\.claude\skills\pmsm-optimizer\scripts\subflow_b_bemf.py`
- Create: `C:\Users\Wesley\.claude\skills\pmsm-optimizer\scripts\subflow_c_torque.py`
- Create: `C:\Users\Wesley\.claude\skills\pmsm-optimizer\scripts\subflow_d_efficiency_map.py`
- Create: `C:\Users\Wesley\.claude\skills\pmsm-optimizer\scripts\subflow_e_external.py`
- Overwrite: `C:\Users\Wesley\.claude\skills\pmsm-optimizer\references\pmsm-methods.md`
- Modify: `C:\Users\Wesley\.claude\skills\pmsm-optimizer\SKILL.md`

**Interfaces:**
- Consumes: Template file `references/Prius_2D_Practice.aedt` with 2 designs (4_MS2, 5_TR)
- Produces: Runnable scripts for each sub-flow

- [ ] **Step 1: Create `scripts/__init__.py`**

```python
# PMSM Optimizer - Sub-flow scripts package
```

- [x] **Step 2: Create `scripts/subflow_a_ldlq.py`**

See current source at `scripts/subflow_a_ldlq.py`. Key implementation details:

- **gRPC 兼容**: 用 `boundary.update()` 修改已有 Current 边界对象，而非 `m2d.assign_current()`
- **磁链提取**: 用 `export_matrix('Matrix1', file)` → 解析 .txt 文件，而非 `get_solution_data()` 的 `FluxLinkage()` 表达式
- **极对数**: `float(m2d['Poles']) / 2` — 模板变量名是 `Poles` 而非 `PolePairs`
- **转子位置标定**: 零电流点测量 → Clarke 变换 → `atan2(ψβ, ψα)` 获取转子实际 d 轴角度 θr
- **Park 变换**: 使用实测 θr（非固定 θ=0°），避免 d/q 轴分量串扰
- **Ld 计算**: `(Ψd - Φ)/Id`（扣除永磁磁链）而非 `Ψd/Id`

- [x] **Step 3: Create `scripts/subflow_b_bemf.py`**

See current source at `scripts/subflow_b_bemf.py`. Key implementation details:

- **步数**: 默认 50 步/周期
- **表达式**: `InducedVoltage(Phase_A)`（Maxwell Transient 下用 InducedVoltage 提取反电势）
- **API**: `get_solution_data_per_variation()` + `data.data_real()` + `data.primary_sweep_values`（gRPC 兼容）
- **数据转换**: InducedVoltage 返回 mV，除 1000→V；时间轴返回 ns，乘 1e-9→s
- **极对数**: `float(m2d['Poles']) / 2`
- **FFT**: `np.fft.rfft()`，`2/N` 幅值校正，THD = √(Σ|Hk|²)/|H1| (k≥2)
- **Ke**: `V_fund × √3 / (Speed_rpm/1000)` V/(krpm)

- [x] **Step 4: Create `scripts/subflow_c_torque.py`**

See current source at `scripts/subflow_c_torque.py`. Key implementation details:

- **步数**: MTPA 扫描 50 步/周期，精确结算 50 步/周期
- **MTPA 候选角**: [0, -15, -25, -35, -45, -60] 度，各跑 1 电周期
- **稳态平均**: 取最后半周期 `torque[-half_period:]` 平均（基于实际步数）
- **`Thet_deg` 赋值**: 用 `str(angle)` 不带 `°` 后缀 — gRPC 下 `°` 后缀在第二次赋值后静默失败
- **API**: `get_solution_data_per_variation()` + `data.data_real()`
- **极对数**: `float(m2d['Poles']) / 2`

- [x] **Step 5: Create `scripts/subflow_d_efficiency_map.py`**

See current source at `scripts/subflow_d_efficiency_map.py`. 完整实现：

- `_mtpa_for_torque()` — MTPA 算法从目标扭矩生成 (Id, Iq, β)，沿 MTPA 轨迹二分搜索
- `_calc_analytical()` — 基于 (Id, Iq) 解析计算 Vs / M / PF / 损耗/效率
- FEA 仅提取 Moving1.Torque 单值（gRPC 限制），用于验证扭矩匹配度
- 完整的网格扫描循环 — **7 张 MAP** 输出（T_avg / η / M / PF / Vs / Is / β）
- 每点输出 CSV 文件（`results/*.csv`），带行列标签
- 使用 `pole_pairs = float(m2d['Poles']) / 2`
- **50 步/周期**

- [x] **Step 6: Create `scripts/subflow_e_external.py`**

See current source at `scripts/subflow_e_external.py`. 核心逻辑：
- 纯数学计算（不跑 FEA），MTPA 区 + 弱磁区双区搜索
- 500 点网格搜索 MTPA 最优工作点
- 转速使用 `np.geomspace` 几何分布（低速密高速疏）
- 输出 T-n 曲线表 + 工作点轨迹 + 转折点 + 峰值功率
- `pole_pairs` 是函数参数（纯数学计算，不连接 PyAEDT）

- [x] **Step 7: Overwrite `references/pmsm-methods.md` with formulas + module refs**

See current file at `references/pmsm-methods.md`. 包含：

- **Transient 求解器公共配置**章节（公式、参数表、物理含义）
- **Sub-flow A 的转子位置标定 + 磁链法**完整描述
- **Sub-flow D 的 MTPA 算法**（关键公式 + 二分搜索流程）
- **解析电参数公式**（Vd = -ω·Lq·Iq, Vq = ω·(Ld·Id+Φ), PF, M 等）
- 所有步数统一为 **50 步/周期**

- [x] **Step 8: Update `SKILL.md`**

已完成三项编辑：
1. 模板名 `pmsm_template` → `Prius_2D_Practice`
2. 参数表替换为实际模板变量（Imax/Thet_deg/Speed_rpm/PolePairs）
3. 子流程引用从 `references/pmsm-methods.md` → `scripts/` 目录

- [x] **Step 9: Create results directory**

```python
python -c "import os; os.makedirs('results', exist_ok=True)"
```

- [x] **Step 10: Commit**

```bash
git add scripts/ SKILL.md references/pmsm-methods.md results/
git commit -m "feat: implement all 5 PMSM sub-flows
- Create scripts/ modules with complete run() functions
- Sub-flow A: Ld/Lq MAP via boundary.update() with 9x turns
- Sub-flow B: Back EMF with FFT + THD + Ke
- Sub-flow C: Rated torque with auto MTPA scan
- Sub-flow D: Efficiency MAP grid scan (7 CSVs)
- Sub-flow E: External characteristic (math-based)
- Update SKILL.md with template variable table
- Update pmsm-methods.md with formulas and module refs"
```

---

### Verification

After all steps, run each sub-flow and validate output:

1. **Sub-flow A**: `python -c "from scripts.subflow_a_ldlq import run; run(rated_current=250, current_steps=3, angle_steps=3)"` — Verify Ld < Lq (IPM typical), saturation with increasing current.

2. **Sub-flow B**: `python -c "from scripts.subflow_b_bemf import run; run(rated_speed=3000)"` — Verify THD < 15%, Ke ~50-80 V/krpm for Prius-class motor.

3. **Sub-flow C**: `python -c "from scripts.subflow_c_torque import run; run(rated_current=250, rated_speed=3000)"` — Verify torque > 0, auto MTPA finds optimal angle.

4. **Sub-flow D**: `python -c "from scripts.subflow_d_efficiency_map import run; run(speed_steps=3, torque_steps=3)"` — Verify 7 MAP CSVs saved with smooth variation.

5. **Sub-flow E**: `python -c "from scripts.subflow_e_external import run; run(vdc=300, imax=250, speed_max=12000)"` — Verify torque decreases at high speed (flux weakening).

---

## 下一阶段优化实施

### 优化 1：项目隔离机制

- [ ] **Step 1: 创建项目辅助函数** — 在 `scripts/` 下新增 `project_utils.py`：
  - `create_project(description: str) -> Path`: 创建 `pmsm_projects/YYYY-MM-DD_<描述>/`，拷贝模板
  - `load_project(project_name: str) -> Path`: 按名称找到已有项目
  - 在项目目录内生成 `project.json`（记录创建时间、参数、模板路径）
- [ ] **Step 2: 修改子流程脚本** — 所有 5 个脚本的 `run()` 函数增加 `project_path` 参数：
  - `TEMPLATE` 从固定路径改为 `project_path / 'Prius_2D_Practice.aedt'`
  - 向后兼容：未传 `project_path` 时使用原始模板路径
- [ ] **Step 3: 更新 SKILL.md** — 在执行流程中加入项目创建步骤

### 优化 2：参数权限管控

- [ ] **Step 4: 创建参数守卫模块** — 在 `scripts/` 下新增 `param_guard.py`：
  ```python
  ALLOWED_VARS = {
      'A': {'Imax'},
      'B': {'Speed_rpm'},
      'C': {'Imax', 'Speed_rpm', 'Thet_deg'},
      'D': {'Speed_rpm', 'Imax', 'Thet_deg'},
  }
  FORBIDDEN_VARS = {'Poles', 'PolePairs'}  # 全局禁止
  
  def check_var(name, subflow):
      if name in FORBIDDEN_VARS:
          raise PermissionError(f"全局禁止修改参数 '{name}'")
      if name not in ALLOWED_VARS.get(subflow, set()):
          raise PermissionError(f"[{subflow}] 不允许修改 '{name}'")
  ```
- [ ] **Step 5: 在 sub-flow 脚本中集成检查** — 每个脚本在 `m2d[...] = value` 前调用 `check_var()`
- [ ] **Step 6: 更新 SKILL.md** — 加入参数权限章节，提前知会 AI
- [ ] **Step 7: 测试与迭代** — 跑一遍验证流程，发现问题及时调整
