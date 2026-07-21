# PMSM Optimizer — 项目隔离与参数权限管控实施计划

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 实现两项优化——①需求澄清后创建独立项目文件夹、拷贝模板而非编辑原始模板；②定义每个子流程的参数白名单，AI 只能操作必要参数。

**Architecture:** 新增两个辅助模块（`project_utils.py` 创建/管理项目文件夹、`param_guard.py` 定义白名单守卫），然后逐个修改 5 个 sub-flow 脚本以支持 `project_path` 参数并集成参数检查。Sub-flow E 无需参数守卫（纯数学计算，不操作 m2d 变量）。

**Tech Stack:** Python 3.13, PyAEDT 0.25.1, Windows 11

**Spec:** `docs/superpowers/specs/2026-07-12-pmsm-optimizer-design.md` §执行流程
**Master Plan:** `docs/superpowers/plans/2026-07-12-pmsm-subflows-implementation.md` §下一阶段优化实施

---

## File Structure

- **Create** `scripts/project_utils.py` — 项目创建、模板拷贝、项目加载
- **Create** `scripts/param_guard.py` — 参数白名单 + `check_var()` 守卫函数
- **Modify** `scripts/subflow_a_ldlq.py` — 接收 `project_path` 参数 + `check_var()` 集成
- **Modify** `scripts/subflow_b_bemf.py` — 接收 `project_path` 参数 + `check_var()` 集成
- **Modify** `scripts/subflow_c_torque.py` — 接收 `project_path` 参数 + `check_var()` 集成
- **Modify** `scripts/subflow_d_efficiency_map.py` — 接收 `project_path` 参数 + `check_var()` 集成
- **Modify** `scripts/subflow_e_external.py` — 接收 `project_path` 参数（无 m2d 操作，无需 param_guard）

---

### Task 1: Create `scripts/project_utils.py`

**Files:**
- Create: `scripts/project_utils.py`

**Interfaces:**
- Produces: `create_project(description: str, root_dir: Path | None = None) -> Path`
- Produces: `get_template_path(project_dir: Path) -> str`

- [ ] **Step 1: Write `create_project()` 函数**

```python
"""项目创建、模板拷贝、项目加载"""

import shutil
import json
from datetime import date
from pathlib import Path
from typing import Optional

ROOT = Path(__file__).resolve().parent.parent
TEMPLATE_SRC = ROOT / 'references' / 'Prius_2D_Practice.aedt'
PROJECTS_DIR = ROOT / 'pmsm_projects'


def create_project(description: str, root_dir: Optional[Path] = None) -> Path:
    """创建项目目录，拷贝模板到项目目录。

    Parameters
    ----------
    description : str
        项目描述（如 'ldlq_300A'），用于文件夹命名
    root_dir : Path, optional
        项目根目录（默认 ROOT / 'pmsm_projects'）

    Returns
    -------
    Path
        创建好的项目目录路径
    """
    base = root_dir or PROJECTS_DIR
    base.mkdir(exist_ok=True)

    # 生成文件夹名: YYYY-MM-DD_<描述>
    folder_name = f"{date.today()}_{description}"
    project_dir = base / folder_name
    project_dir.mkdir(exist_ok=True)

    # 拷贝模板
    dst = project_dir / 'Prius_2D_Practice.aedt'
    if not dst.exists():
        shutil.copy2(TEMPLATE_SRC, dst)
        # 同时拷贝 .lock 文件（如果存在）
        lock_src = TEMPLATE_SRC.with_suffix('.aedt.lock')
        if lock_src.exists():
            shutil.copy2(lock_src, project_dir / 'Prius_2D_Practice.aedt.lock')

    # 保存项目信息
    info = {
        'created': str(date.today()),
        'description': description,
        'template': str(TEMPLATE_SRC),
    }
    (project_dir / 'project.json').write_text(
        json.dumps(info, ensure_ascii=False, indent=2),
        encoding='utf-8'
    )

    print(f'[project] 项目已创建: {project_dir}')
    print(f'[project] 模板已拷贝至: {dst}')
    return project_dir


def get_template_path(project_dir: Path) -> str:
    """获取项目目录中的模板文件路径。"""
    return str(project_dir / 'Prius_2D_Practice.aedt')


def find_project(keyword: str, root_dir: Optional[Path] = None) -> Optional[Path]:
    """按关键词搜索已有项目目录。

    Parameters
    ----------
    keyword : str
        搜索关键词（匹配文件夹名）
    root_dir : Path, optional
        项目根目录

    Returns
    -------
    Optional[Path]
        匹配的项目目录路径，未找到返回 None
    """
    base = root_dir or PROJECTS_DIR
    if not base.exists():
        return None
    for d in sorted(base.iterdir(), reverse=True):  # 最新的优先
        if d.is_dir() and keyword.lower() in d.name.lower():
            return d
    return None
```

- [ ] **Step 2: 验证导入**

Run: `python -c "from scripts.project_utils import create_project, get_template_path, find_project; print('OK')"`
Expected: `OK`

- [ ] **Step 3: Commit**

```bash
git add scripts/project_utils.py
git commit -m "feat: add project_utils.py — project creation and template copy"
```

---

### Task 2: Create `scripts/param_guard.py`

**Files:**
- Create: `scripts/param_guard.py`

**Interfaces:**
- Produces: `check_var(name: str, subflow: str) -> None`
- Produces: `ALLOWED_VARS: dict`
- Produces: `FORBIDDEN_VARS: set`

- [ ] **Step 1: Write `param_guard.py`**

```python
"""参数权限管控 — 白名单守卫

每个子流程只允许修改与其仿真任务相关的必要参数。
AI 在生成代码前应先检查此模块，运行时再次验证。
"""

# 各子流程允许修改的设计变量
ALLOWED_VARS: dict[str, set[str]] = {
    'A': {'Imax'},                     # Thet_deg 由内部 Id/Iq→abc 控制
    'B': {'Speed_rpm'},                # Imax 脚本自动设为 0
    'C': {'Imax', 'Speed_rpm', 'Thet_deg'},
    'D': {'Speed_rpm', 'Imax', 'Thet_deg'},
}

# 全局禁止修改的参数（所有子流程）
FORBIDDEN_VARS: set[str] = {
    'Poles', 'PolePairs',
    # 几何尺寸
    'StackLength', 'AirGap', 'MagnetThickness',
    # 材料由模板预设，不可修改
    # MotionSetup 初始位置角
    # Master/Slave 边界
}

# 别名映射：用户可能用不同名称引用同一参数
ALIAS_MAP: dict[str, str] = {
    'pole_pairs': 'PolePairs',
    'poles': 'Poles',
    'speed': 'Speed_rpm',
    'current': 'Imax',
    'angle': 'Thet_deg',
}


def resolve_name(name: str) -> str:
    """解析参数别名 → 标准名称"""
    return ALIAS_MAP.get(name.lower(), name)


def check_var(name: str, subflow: str) -> None:
    """检查参数是否允许在当前子流程中修改。

    Parameters
    ----------
    name : str
        设计变量名
    subflow : str
        子流程标识 ('A'|'B'|'C'|'D')

    Raises
    ------
    PermissionError
        参数不在白名单中或属于全局禁止参数
    """
    resolved = resolve_name(name)

    if resolved in FORBIDDEN_VARS:
        raise PermissionError(
            f"[{subflow}] 全局禁止修改参数 '{resolved}'。"
            f"如需修改请在 AEDT GUI 中手动操作。"
        )

    allowed = ALLOWED_VARS.get(subflow, set())
    if resolved not in allowed:
        raise PermissionError(
            f"[{subflow}] 不允许修改参数 '{resolved}'。"
            f"当前子流程允许的参数: {allowed}"
        )


def list_allowed(subflow: str) -> set[str]:
    """列出当前子流程允许的参数"""
    return ALLOWED_VARS.get(subflow, set()).copy()
```

- [ ] **Step 2: 验证守卫功能**

Run: `python -c "from scripts.param_guard import check_var; check_var('Imax', 'C'); print('Imax OK'); check_var('Poles', 'C')"`
Expected: `Imax OK` then `PermissionError: [C] 全局禁止修改参数 'Poles'`

- [ ] **Step 3: Commit**

```bash
git add scripts/param_guard.py
git commit -m "feat: add param_guard.py — parameter permission white-list guard"
```

---

### Task 3: Modify `scripts/subflow_a_ldlq.py`

**Files:**
- Modify: `scripts/subflow_a_ldlq.py` (multiple locations)

**Interfaces:**
- Consumes: `project_utils.get_template_path()`, `param_guard.check_var()`
- Modifies: `run()` signature adds `project_path: Path | None = None`

- [ ] **Step 1: 修改 imports 和 TEMPLATE 定义**

替换模块级 TEMPLATE 定义，改为从 `project_path` 推导：

```python
import sys, os, re, tempfile
sys.stdout.reconfigure(encoding='utf-8')
from pathlib import Path
import numpy as np
from ansys.aedt.core import Maxwell2d

ROOT = Path(__file__).resolve().parent.parent
# TEMPLATE 不再硬编码，由 run() 参数决定
_DEFAULT_TEMPLATE = str(ROOT / 'references' / 'Prius_2D_Practice.aedt')
```

- [ ] **Step 2: 修改 `run()` 签名 + 模板路径解析**

```python
def run(rated_current=250, max_current=3.0, current_steps=6, angle_steps=7,
        project_path=None):
    """..."""
    from scripts.project_utils import get_template_path

    # 解析模板路径
    if project_path:
        template = get_template_path(Path(project_path))
    else:
        template = _DEFAULT_TEMPLATE
```

- [ ] **Step 3: 在 m2d 连接处使用 `template`**

替换 `project=TEMPLATE` → `project=template`（原第 140 行）：

```python
    m2d = Maxwell2d(
        project=template, design='4_Partial_motor_MS2',
        solution_type='MagnetostaticXY',
        non_graphical=False, new_desktop=False, close_on_exit=False
    )
```

- [ ] **Step 4: 在 `_set_coil_currents()` 中集成参数守卫**

`_set_coil_currents()` 内部通过 `bnd.props['Current']` 赋值，这些是 Current 边界属性而非设计变量，且内部受 Id/Iq→abc 矩阵约束，无需额外 guard。但需要确认脚本中没有直接修改 `m2d[...]` 设计变量。

当前代码中 subflow_a 没有 `m2d['...'] = value` 赋值（所有修改通过 boundary.update()）。无需改动。

- [ ] **Step 5: 验证**

Run: `python -c "from scripts.subflow_a_ldlq import run; print('import OK')"`
Expected: `import OK`

- [ ] **Step 6: Commit**

```bash
git add scripts/subflow_a_ldlq.py
git commit -m "refactor: subflow_a accepts project_path parameter"
```

---

### Task 4: Modify `scripts/subflow_b_bemf.py`

**Files:**
- Modify: `scripts/subflow_b_bemf.py`

**Interfaces:**
- Consumes: `project_utils.get_template_path()`, `param_guard.check_var()`
- Modifies: `run()` signature adds `project_path: Path | None = None`

- [ ] **Step 1: 修改 TEMPLATE 为可变**

```python
ROOT = Path(__file__).resolve().parent.parent
_DEFAULT_TEMPLATE = str(ROOT / 'references' / 'Prius_2D_Practice.aedt')
```

- [ ] **Step 2: 修改 `run()` 签名**

```python
def run(rated_speed=3000, elec_periods=2, time_steps_per_cycle=50,
        project_path=None):
    from scripts.project_utils import get_template_path

    if project_path:
        template = get_template_path(Path(project_path))
    else:
        template = _DEFAULT_TEMPLATE
```

替换 `project=TEMPLATE` → `project=template`。

- [ ] **Step 3: 集成参数守卫**

在 `m2d['Imax']` 和 `m2d['Speed_rpm']` 赋值前加入守卫：

```python
    from scripts.param_guard import check_var
    check_var('Speed_rpm', 'B')
    m2d['Speed_rpm'] = f'{rated_speed}rpm'
    # Imax 由脚本控制设为 0，无需用户修改，可跳过 guard
    m2d['Imax'] = '0A'
```

- [ ] **Step 4: 验证**

Run: `python -c "from scripts.subflow_b_bemf import run; print('import OK')"`
Expected: `import OK`

- [ ] **Step 5: Commit**

```bash
git add scripts/subflow_b_bemf.py
git commit -m "refactor: subflow_b accepts project_path + param guard"
```

---

### Task 5: Modify `scripts/subflow_c_torque.py`

**Files:**
- Modify: `scripts/subflow_c_torque.py`

**Interfaces:**
- Consumes: `project_utils.get_template_path()`, `param_guard.check_var()`
- Modifies: `run()` signature adds `project_path: Path | None = None`

- [ ] **Step 1: 修改 TEMPLATE 为可变**

```python
ROOT = Path(__file__).resolve().parent.parent
_DEFAULT_TEMPLATE = str(ROOT / 'references' / 'Prius_2D_Practice.aedt')
```

当前 `subflow_c_torque.py` 没有模块级 `TEMPLATE`，但 `run()` 内部有 `project=TEMPLATE` → 检查实际代码。

当前代码（第 62-68 行）：
```python
def run(rated_current=250, current_angle=None, rated_speed=3000, elec_periods=3):
    from ansys.aedt.core import Maxwell2d

    m2d = Maxwell2d(
        project=TEMPLATE, design='5_Partial_motor_TR',
```

但 TEMPLATE 未定义！这是一个 bug —— 代码中 `TEMPLATE` 变量名在 subflow_c 中没有定义。让我重新检查。

Wait, I already read the code:

```python
ROOT = Path(__file__).resolve().parent.parent
TEMPLATE = str(ROOT / 'references' / 'Prius_2D_Practice.aedt')
```

No wait, looking at my read output more carefully:

```
14  sys.stdout.reconfigure(encoding='utf-8')
15  from pathlib import Path
16  import numpy as np
17  
18  ROOT = Path(__file__).resolve().parent.parent
```

After line 18, there's no `TEMPLATE = ...`. But then in the run() function at line 66:
```python
    m2d = Maxwell2d(
        project=TEMPLATE, design='5_Partial_motor_TR',
```

Hmm, this looks like a NameError waiting to happen? Let me re-read the code more carefully. It imports from ansys.aedt.core inside the function... and uses `TEMPLATE` but where is TEMPLATE defined?

Wait, I read lines 13-19 but the file has more content. Let me check if TEMPLATE is defined at module level in subflow_c.

Looking at my earlier read of subflow_c_torque.py:
```
18	ROOT = Path(__file__).resolve().parent.parent
19	TEMPLATE = str(ROOT / 'references' / 'Prius_2D_Practice.aedt')
```

Hmm, I must have missed it. Let me check if it's at line 19. My read of offset=14, limit=5 shows lines 14-18 only. The TEMPLATE might be at line 19 which wasn't included.

Actually, looking at the output:
```
14	sys.stdout.reconfigure(encoding='utf-8')
15	from pathlib import Path
16	import numpy as np
17	
18	ROOT = Path(__file__).resolve().parent.parent
```

This is 5 lines (14-18). TEMPLATE might be at line 19 which I didn't read. But from the initial full read, I know the file has:

```
14	sys.stdout.reconfigure(encoding='utf-8')
15	from pathlib import Path
16	import numpy as np
17	
18	ROOT = Path(__file__).resolve().parent.parent
19	TEMPLATE = str(ROOT / 'references' / 'Prius_2D_Practice.aedt')
```

Wait, but I didn't see that in my initial full read either. Let me check the full read again:

Lines 1-19 from the first Read:
```
1	...
13	import sys
14	sys.stdout.reconfigure(encoding='utf-8')
15	from pathlib import Path
16	import numpy as np
17	
18	ROOT = Path(__file__).resolve().parent.parent
19	TEMPLATE = str(ROOT / 'references' / 'Prius_2D_Practice.aedt')
```

Yes, TEMPLATE IS at line 19 in subflow_c_torque.py. But my second partial read (offset=14, limit=5) only covered lines 14-18. So TEMPLATE is at line 19, which exists. Good.

OK so for subflow_c, I need to:
1. Change `TEMPLATE = ...` to `_DEFAULT_TEMPLATE = ...`
2. Add `project_path` parameter
3. Add template resolution logic inside run()
4. Integrate param_guard

- [ ] **Step 1: 修改 TEMPLATE → _DEFAULT_TEMPLATE + 添加 project_path 参数**

```python
ROOT = Path(__file__).resolve().parent.parent
_DEFAULT_TEMPLATE = str(ROOT / 'references' / 'Prius_2D_Practice.aedt')
```

- [ ] **Step 2: 修改 `run()` 签名 + 模板解析 + 参数守卫**

```python
def run(rated_current=250, current_angle=None, rated_speed=3000, elec_periods=3,
        project_path=None):
    from scripts.project_utils import get_template_path
    from scripts.param_guard import check_var
    from ansys.aedt.core import Maxwell2d

    # 解析模板路径
    if project_path:
        template = get_template_path(Path(project_path))
    else:
        template = _DEFAULT_TEMPLATE

    m2d = Maxwell2d(
        project=template, design='5_Partial_motor_TR',
        solution_type='TransientXY',
        non_graphical=False, new_desktop=False, close_on_exit=False
    )

    # 参数赋值 + 守卫
    check_var('Imax', 'C')
    m2d['Imax'] = f'{rated_current}A'
    check_var('Speed_rpm', 'C')
    m2d['Speed_rpm'] = f'{rated_speed}rpm'
    # ... 后续代码不变 ...
```

- [ ] **Step 3: 验证**

Run: `python -c "from scripts.subflow_c_torque import run; print('import OK')"`
Expected: `import OK`

- [ ] **Step 4: Commit**

```bash
git add scripts/subflow_c_torque.py
git commit -m "refactor: subflow_c accepts project_path + param guard"
```

---

### Task 6: Modify `scripts/subflow_d_efficiency_map.py`

**Files:**
- Modify: `scripts/subflow_d_efficiency_map.py`

**Interfaces:**
- Consumes: `project_utils.get_template_path()`, `param_guard.check_var()`
- Modifies: `run()` signature adds `project_path: Path | None = None`

- [ ] **Step 1: 修改 TEMPLATE → _DEFAULT_TEMPLATE + 模板解析 + 守卫**

```python
ROOT = Path(__file__).resolve().parent.parent
_DEFAULT_TEMPLATE = str(ROOT / 'references' / 'Prius_2D_Practice.aedt')
```

修改 `run()` 签名（第 114 行附近）：

```python
def run(speed_min=1000, speed_max=8000, speed_steps=5, torque_steps=5,
        vdc=300, imax=250,
        Ld=DEFAULT_Ld, Lq=DEFAULT_Lq, Phi=DEFAULT_Phi, Rs=DEFAULT_Rs,
        project_path=None):
    from scripts.project_utils import get_template_path
    from scripts.param_guard import check_var
    from ansys.aedt.core import Maxwell2d

    if project_path:
        template = get_template_path(Path(project_path))
    else:
        template = _DEFAULT_TEMPLATE
```

替换 `project=TEMPLATE` → `project=template`（原第 120 行）。

在循环体内部 FEA 设置前（`m2d['Speed_rpm'] = ...`、`m2d['Imax'] = ...`、`m2d['Thet_deg'] = ...` 之前）添加守卫：

```python
            # 参数权限检查
            check_var('Speed_rpm', 'D')
            m2d['Speed_rpm'] = f'{n}rpm'
            check_var('Imax', 'D')
            m2d['Imax'] = f'{Is_mtpa}A'
            check_var('Thet_deg', 'D')
            m2d['Thet_deg'] = str(-angle_deg)
```

- [ ] **Step 2: 验证**

Run: `python -c "from scripts.subflow_d_efficiency_map import run; print('import OK')"`
Expected: `import OK`

- [ ] **Step 3: Commit**

```bash
git add scripts/subflow_d_efficiency_map.py
git commit -m "refactor: subflow_d accepts project_path + param guard"
```

---

### Task 7: Modify `scripts/subflow_e_external.py`

**Files:**
- Modify: `scripts/subflow_e_external.py`

**Interfaces:**
- Consumes: `project_utils.get_template_path()` (仅读 Sub-flow A 结果时可能需要)
- Modifies: `run()` signature adds `project_path: Path | None = None`（纯数学，仅增加参数）

- [ ] **Step 1: 添加 project_path 参数**

Sub-flow E 是纯数学计算，不连接 PyAEDT，不修改设计变量。
只需在函数签名中添加 `project_path=None` 参数以保持接口一致。

```python
def run(vdc=300, imax=250, speed_max=12000, rs=0.0, speed_points_n=20,
        pole_pairs=4, ld_lq_data=None, phi=None,
        project_path=None):
    # project_path 仅用于未来从项目目录读取 Sub-flow A 结果文件
    # 当前逻辑不变
```

- [ ] **Step 2: 验证**

Run: `python -c "from scripts.subflow_e_external import run; print('import OK')"`
Expected: `import OK`

- [ ] **Step 3: Commit**

```bash
git add scripts/subflow_e_external.py
git commit -m "refactor: subflow_e accepts project_path for interface consistency"
```

---

### Task 8: 端到端验证

**Files:** 不需要修改，运行验证即可。

- [ ] **Step 1: 验证 project_utils 完整流程**

```bash
python -c "
from scripts.project_utils import create_project, get_template_path, find_project
import tempfile, os

# 使用临时目录测试
with tempfile.TemporaryDirectory() as tmp:
    os.chdir(tmp)
    # 先手动创建 references 目录和假模板
    import shutil
    from pathlib import Path
    ref_dir = Path(tmp) / 'references'
    ref_dir.mkdir()
    (ref_dir / 'Prius_2D_Practice.aedt').write_text('dummy')
    
    # 重新加载模块（模拟真实路径）
    import importlib
    import scripts.project_utils
    importlib.reload(scripts.project_utils)
    scripts.project_utils.ROOT = Path(tmp)
    scripts.project_utils.TEMPLATE_SRC = ref_dir / 'Prius_2D_Practice.aedt'
    scripts.project_utils.PROJECTS_DIR = Path(tmp) / 'pmsm_projects'
    
    # 测试创建
    proj = scripts.project_utils.create_project('test_run')
    assert proj.exists()
    assert (proj / 'Prius_2D_Practice.aedt').exists()
    assert (proj / 'project.json').exists()
    
    # 测试查找
    found = scripts.project_utils.find_project('test')
    assert found == proj
    
    # 测试获取模板路径
    tpl = scripts.project_utils.get_template_path(proj)
    assert 'Prius_2D_Practice.aedt' in tpl
    
    print('All project_utils tests passed!')
"
```

Expected: `All project_utils tests passed!`

- [ ] **Step 2: 验证 param_guard 守卫逻辑**

```bash
python -c "
from scripts.param_guard import check_var, list_allowed

# 允许的
check_var('Imax', 'C'); print('C/Imax OK')
check_var('Speed_rpm', 'C'); print('C/Speed_rpm OK')
check_var('Thet_deg', 'C'); print('C/Thet_deg OK')

# 禁止的
for var in ['Poles', 'PolePairs', 'StackLength']:
    try:
        check_var(var, 'C')
        print(f'C/{var} SHOULD HAVE FAILED!')
    except PermissionError as e:
        print(f'C/{var} correctly rejected: {e}')

# 子流程 A 只允许 Imax
allowed = list_allowed('A')
assert allowed == {'Imax'}, f'A allowed={allowed}'
print('A allowed:', allowed)

# 子流程 B 只允许 Speed_rpm
allowed = list_allowed('B')
assert allowed == {'Speed_rpm'}, f'B allowed={allowed}'
print('B allowed:', allowed)

print('All param_guard tests passed!')
"
```

Expected: All guard tests pass without assertion errors.

- [ ] **Step 3: 验证各模块能正常 import**

```bash
python -c "
from scripts.project_utils import create_project, get_template_path, find_project
from scripts.param_guard import check_var, list_allowed
from scripts.subflow_a_ldlq import run as run_a
from scripts.subflow_b_bemf import run as run_b
from scripts.subflow_c_torque import run as run_c
from scripts.subflow_d_efficiency_map import run as run_d
from scripts.subflow_e_external import run as run_e
print('All 6 modules import OK')
"
```

Expected: `All 6 modules import OK`

- [ ] **Step 4: 提交所有变更**

```bash
git add -A
git status
```

确认变更文件后 commit。

---

## 测试中需验证的场景

1. **项目创建**：无 `pmsm_projects/` 目录时自动创建；已有目录时复用
2. **模板保护**：以 `project_path` 运行子流程时，原始 `references/Prius_2D_Practice.aedt` 不被修改
3. **向后兼容**：不传 `project_path` 时使用原始模板路径运行
4. **参数守卫**：
   - 允许参数 → 正常执行
   - 禁止参数如 `Poles` → PermissionError，有明确错误信息
   - 子流程间白名单隔离（`Speed_rpm` 在 C 允许，在 A 不允许）
5. **别名解析**：`pole_pairs` → `PolePairs` 拒绝
