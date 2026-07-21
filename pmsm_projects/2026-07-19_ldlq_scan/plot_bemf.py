"""Sub-flow B: BEMF @ 4000rpm — 仿真 + 绘图"""
import sys, os
sys.stdout.reconfigure(encoding='utf-8')
os.environ['MPLBACKEND'] = 'Agg'
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from ansys.aedt.core import Maxwell2d

ROOT = Path(__file__).resolve().parent.parent.parent
PROJECT_DIR = ROOT / 'pmsm_projects' / '2026-07-19_ldlq_scan'
TEMPLATE = str(PROJECT_DIR / 'Prius_2D_Practice.aedt')

rated_speed = 4000
elec_periods = 3
steps_per_cycle = 50

m2d = Maxwell2d(
    project=TEMPLATE, design='5_Partial_motor_TR',
    solution_type='TransientXY',
    non_graphical=False, new_desktop=False, close_on_exit=False
)

m2d['Imax'] = '0A'
m2d['Speed_rpm'] = f'{rated_speed}rpm'

pole_pairs = float(m2d['Poles']) / 2
freq = rated_speed / 60 * pole_pairs
stop_time = elec_periods / freq
time_step = 1 / (freq * steps_per_cycle)

setup = m2d.setups[0]
setup.props['StopTime'] = f'{stop_time}s'
setup.props['TimeStep'] = f'{time_step}s'
setup.update()

print(f'求解中... ({rated_speed}rpm, {freq:.1f}Hz, {elec_periods}周期)')
m2d.analyze('Setup1')

data = m2d.post.get_solution_data_per_variation(
    expressions=['InducedVoltage(Phase_A)', 'InducedVoltage(Phase_B)', 'InducedVoltage(Phase_C)',
                 'Moving1.Torque']
)
time_ns = np.array(data.primary_sweep_values, dtype=float)
time_s = time_ns * 1e-9
va = np.array(data.data_real('InducedVoltage(Phase_A)')) / 1000
vb = np.array(data.data_real('InducedVoltage(Phase_B)')) / 1000
vc = np.array(data.data_real('InducedVoltage(Phase_C)')) / 1000
torque = np.array(data.data_real('Moving1.Torque'))

# ── 稳态段 —— 丢弃第1周期 ──
n_per_period = steps_per_cycle
steady_start = n_per_period
va_s = va[steady_start:]
n = len(va_s)

# FFT
fft_mag = np.abs(np.fft.rfft(va_s)) * 2 / n
fft_freq = np.fft.rfftfreq(n, d=time_step)
fund_idx = np.argmax(fft_mag[1:]) + 1
fund_freq = fft_freq[fund_idx]
fund_amp = fft_mag[fund_idx]

# THD
mask = np.ones(len(fft_mag), dtype=bool)
mask[0] = False
mask[fund_idx] = False
harmonics_rms = np.sqrt(np.sum(fft_mag[mask] ** 2))
thd = harmonics_rms / fft_mag[fund_idx] * 100
ke_line = fund_amp * np.sqrt(3) / (rated_speed / 1000)

# ── 打印结果 ──
print(f'\n{"="*55}')
print(f'  空载反电势 @ {rated_speed} rpm')
print(f'{"="*55}')
print(f'  电频率:             {freq:.1f} Hz')
print(f'  A相基波幅值(稳态):  {fund_amp:.2f} V')
print(f'  THD(稳态):          {thd:.1f}%')
print(f'  线反电势常数 Ke:    {ke_line:.2f} V/(krpm)')
print(f'{"="*55}\n')

# ── 绘图 ──
fig, axes = plt.subplots(3, 1, figsize=(12, 10))
fig.suptitle(f'Back EMF @ {rated_speed} rpm  (freq = {freq:.1f} Hz)', fontsize=14)

t_ms = time_s * 1000
t_steady_ms = steady_start * time_step * 1000

# 1) 三相 BEMF 波形
ax = axes[0]
ax.plot(t_ms, va, label='Phase A', color='#e74c3c')
ax.plot(t_ms, vb, label='Phase B', color='#2ecc71')
ax.plot(t_ms, vc, label='Phase C', color='#3498db')
ax.axvline(x=t_steady_ms, color='gray', ls='--', lw=0.8, label='steady start')
ax.set_ylabel('Voltage (V)')
ax.legend(loc='upper right', ncol=4)
ax.grid(True, alpha=0.3)

# 2) FFT 频谱
ax2 = axes[1]
n_h = min(30, len(fft_freq))
ax2.stem(fft_freq[:n_h], fft_mag[:n_h], linefmt='C0-', markerfmt='C0o', basefmt=' ')
ax2.set_xlim(0, fft_freq[n_h-1] if n_h < len(fft_freq) else fft_freq[-1])
ax2.set_xlabel('Frequency (Hz)')
ax2.set_ylabel('Amplitude (V)')
ax2.set_title('FFT Spectrum (steady-state only)')
ax2.grid(True, alpha=0.3)
ax2.annotate(f'Fundamental\n{fund_freq:.1f}Hz = {fund_amp:.1f}V',
             xy=(fund_freq, fund_amp), xytext=(fund_freq+80, fund_amp*0.6),
             arrowprops=dict(arrowstyle='->', color='red'), color='red', fontsize=9)

for order, c in [(3, 'orange'), (5, 'green'), (7, 'purple')]:
    hz = fund_freq * order
    idx = np.argmin(np.abs(fft_freq - hz))
    if idx < len(fft_freq) and fft_mag[idx] > 1:
        ax2.annotate(f'{order}× ({fft_mag[idx]:.1f}V)', xy=(fft_freq[idx], fft_mag[idx]),
                     fontsize=8, color=c)

# 3) 齿槽/定位力矩
ax3 = axes[2]
ax3.plot(t_ms, torque, color='#9b59b6')
ax3.axvline(x=t_steady_ms, color='gray', ls='--', lw=0.8, label='steady start')
ax3.set_xlabel('Time (ms)')
ax3.set_ylabel('Torque (Nm)')
ax3.set_title('Torque (unloaded — cogging + reluctance)')
ax3.grid(True, alpha=0.3)
ax3.legend(loc='upper right')

fig.text(0.5, 0.005,
    f'A: {fund_amp:.1f}V @ {fund_freq:.1f}Hz  |  Ke(line) = {ke_line:.2f}V/(krpm)  |  '
    f'THD = {thd:.1f}%  |  PolePairs = {int(pole_pairs)}',
    ha='center', fontsize=9, bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9))

plt.tight_layout(rect=[0, 0.06, 1, 1])
plt.subplots_adjust(hspace=0.3)

out_png = str(PROJECT_DIR / 'bemf_4000rpm.png')
fig.savefig(out_png, dpi=150, bbox_inches='tight')
print(f'已保存: {out_png}')

# 也保存 CSV 波形数据
out_csv = str(PROJECT_DIR / 'bemf_4000rpm.csv')
header = 'time_s,PhaseA_V,PhaseB_V,PhaseC_V,Torque_Nm'
np.savetxt(out_csv, np.column_stack([time_s, va, vb, vc, torque]),
           delimiter=',', header=header, comments='')
print(f'已保存: {out_csv}')

plt.close()
m2d.close_project()
