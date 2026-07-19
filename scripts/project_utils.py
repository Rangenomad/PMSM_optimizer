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
