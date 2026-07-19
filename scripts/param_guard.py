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
