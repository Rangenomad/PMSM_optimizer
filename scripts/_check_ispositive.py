"""Quick check: original IsPositive values from template."""
import sys
sys.stdout.reconfigure(encoding='utf-8')
from pathlib import Path
from ansys.aedt.core import Maxwell2d

TEMPLATE = str(Path('pmsm_projects/2026-07-21_mtpa_test') / 'Prius_2D_Practice.aedt')
m2d = Maxwell2d(project=TEMPLATE, design='4_Partial_motor_MS2',
    solution_type='MagnetostaticXY', non_graphical=False, new_desktop=False, close_on_exit=False)

for name in ['PhaseA1','PhaseA2','PhaseB1','PhaseB2','PhaseC1','PhaseC2']:
    b = [b for b in m2d.boundaries if b.name == name][0]
    cur = b.props.get('Current', '?')
    isp = b.props.get('IsPositive', '?')
    print(f'{name}: Current={cur}, IsPositive={isp}')

m2d.close_project()
