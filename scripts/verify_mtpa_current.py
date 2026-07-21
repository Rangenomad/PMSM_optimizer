"""MTPA verification: scan Thet=0~60° at Imax=100A, 200A"""
import sys
sys.stdout.reconfigure(encoding='utf-8')
from pathlib import Path
import numpy as np, re
from ansys.aedt.core import Maxwell2d

template = "pmsm_projects/2026-07-19_efficiency_map/Prius_2D_Practice.aedt"

# Sync .aedt file
file = Path(template)
content = file.read_text(encoding="utf-8")
content, _ = re.subn(r"'Angular Velocity'='[^']*'", "'Angular Velocity'='2000rpm'", content)
content, _ = re.subn(r"VariableProp\('Speed_rpm', 'UD', '', '[^']*'",
                      "VariableProp('Speed_rpm', 'UD', '', '2000'", content)
file.write_text(content, encoding="utf-8")
print("[init] .aedt synced: AV=2000rpm")

m2d = Maxwell2d(
    project=template, design='5_Partial_motor_TR',
    solution_type='TransientXY',
    non_graphical=False, new_desktop=False, close_on_exit=False
)
pole_pairs = 4
speed = 2000.0
freq = speed / 60 * pole_pairs
stop_time = 1.0 / freq
time_step = 1.0 / (freq * 50)

setup = m2d.setups[0]
setup.props['StopTime'] = f'{stop_time}s'
setup.props['TimeStep'] = f'{time_step}s'
setup.update()

b = [b for b in m2d.boundaries if b.name == 'MotionSetup1'][0]
b.props['Angular Velocity'] = '2000rpm'
b.update()
m2d['Speed_rpm'] = '2000'

angles = [0, 10, 20, 30, 40, 50, 60]
currents = [100, 200]

print(f"\n{'='*60}")
print(f"MTPA verification: 100A, 200A @ 2000rpm")
print(f"PolePairs={pole_pairs}, freq={freq:.1f}Hz, stop={stop_time:.4f}s")
print(f"{'='*60}")

for Imax_val in currents:
    print(f"\n--- Imax = {Imax_val}A ---")
    print(f"  {'Thet':>5s}  {'T_avg (Nm)':>10s}")
    print(f"  {'-'*20}")
    for thet in angles:
        m2d['Imax'] = f'{Imax_val}A'
        m2d['Thet_deg'] = str(thet)
        m2d.analyze('Setup1')
        data = m2d.post.get_solution_data_per_variation(
            expressions=['Moving1.Torque']
        )
        if data:
            tq = np.array(data.data_real('Moving1.Torque'))
            T_avg = float(np.mean(tq)) if len(tq) > 0 else 0.0
        else:
            T_avg = 0.0
        print(f"  {thet:5.0f}  {T_avg:10.1f}")

m2d.close_project()
print("\nDone.")
