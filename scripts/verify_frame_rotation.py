"""Verify MS-to-TR frame rotation and compute corrected Ld/Lq.

Key insight: MS rotor d-axis at -120°, TR d-axis at 0° at t=0.
For the same physical current distribution:
  (Id_tr, Iq_tr) = R(-120°) * (Id_ms, Iq_ms)

To measure Ld/Lq for predicting TR torque at Thet=X°:
  MS_angle = X - 120°  →  Id_ms=Imax*sin(X-120), Iq_ms=Imax*cos(X-120)

These (Id_ms, Iq_ms) correspond to Id_tr = Imax*sin(X), Iq_tr = Imax*cos(X).
Then Ld, Lq computed in MS frame should predict torque via: T = 1.5*P*(Φ*Iq_tr + (Ld-Lq)*Id_tr*Iq_tr)
"""

import sys, os, re, tempfile
sys.stdout.reconfigure(encoding='utf-8')
from pathlib import Path
import numpy as np
from ansys.aedt.core import Maxwell2d

ROOT = Path(__file__).resolve().parent.parent
TURNS = 9

COIL_CONFIG = [
    ('PhaseA1', True), ('PhaseA2', True),
    ('PhaseB1', True), ('PhaseB2', True),
    ('PhaseC1', True), ('PhaseC2', True),
]
COIL_NAMES = [c[0] for c in COIL_CONFIG]


def id_iq_to_abc(Id, Iq, theta=0):
    Ialpha = Id * np.cos(theta) - Iq * np.sin(theta)
    Ibeta  = Id * np.sin(theta) + Iq * np.cos(theta)
    Ia = Ialpha
    Ib = -0.5 * Ialpha + np.sqrt(3)/2 * Ibeta
    Ic = -0.5 * Ialpha - np.sqrt(3)/2 * Ibeta
    return Ia, Ib, Ic


def park_abc_to_dq(psi_a, psi_b, psi_c, theta=0):
    c = np.cos(theta)
    s = np.sin(theta)
    psi_d = 2/3 * (psi_a*c + psi_b*np.cos(theta-2*np.pi/3) + psi_c*np.cos(theta+2*np.pi/3))
    psi_q = -2/3 * (psi_a*s + psi_b*np.sin(theta-2*np.pi/3) + psi_c*np.sin(theta+2*np.pi/3))
    return psi_d, psi_q


def set_coil_currents(m2d, boundaries, Ia, Ib, Ic):
    currents = [TURNS*Ia, TURNS*Ia, TURNS*Ib, TURNS*Ib, TURNS*Ic, TURNS*Ic]
    for bnd, (name, is_pos), val in zip(boundaries, COIL_CONFIG, currents):
        bnd.props['Current'] = f'{val}A'
        bnd.props['IsPositive'] = is_pos
        bnd.update()


def parse_matrix_flux_linkage(text):
    in_flux = False
    result = {}
    for line in text.splitlines():
        if line.strip().startswith('Flux Linkage'):
            in_flux = True
            continue
        if in_flux and line.strip() and not line.strip().startswith('Flux'):
            parts = line.strip().split()
            if len(parts) >= 2 and parts[0] in COIL_NAMES:
                result[parts[0]] = float(parts[1])
    return result


def extract_flux_linkages(m2d):
    m2d.analyze('Setup1')
    out_file = os.path.join(tempfile.gettempdir(), 'pmsm_matrix_export.txt')
    m2d.export_matrix('Matrix1', out_file)
    with open(out_file, 'r') as f:
        text = f.read()
    return parse_matrix_flux_linkage(text)


# FEA truth (reliable values)
FEA = {
    20: 296.1,   # accurate 3-period
    50: 371.1,   # accurate 3-period
}

P = 4
IMAX = 350.0
THETA_R_MS = -120.0  # degrees (measured from zero-current calibration)

# MS_angle = TR_Thet + THETA_R_MS
# But THETA_R_MS = -120°, so MS_angle = TR_Thet - 120°

# Test: TR Thet=50° and Thet=20° scanned from MS
test_angles = [
    # (MS_angle, TR_Thet_equiv, FEA_T)
    (-70, 50, FEA[50]),    # -70 = 50 - 120
    (-100, 20, FEA[20]),   # -100 = 20 - 120
]

print("=" * 80)
print("Verification: MS→TR frame rotation for Ld/Lq")
print(f"  Rotor offset: theta_r(MS) = {THETA_R_MS} deg, theta_r(TR@t=0) = 0 deg")
print(f"  MS_angle = TR_Thet + {THETA_R_MS} deg")
print("=" * 80)

project_path = 'pmsm_projects/2026-07-19_dq_fix'
template = str(Path(project_path) / 'Prius_2D_Practice.aedt')

m2d = Maxwell2d(
    project=template, design='4_Partial_motor_MS2',
    solution_type='MagnetostaticXY',
    non_graphical=False, new_desktop=False, close_on_exit=False
)

boundaries = []
for name, _ in COIL_CONFIG:
    b = [b for b in m2d.boundaries if b.name == name][0]
    boundaries.append(b)

# Zero-current calibration
print("\n[1] Zero-current calibration:")
set_coil_currents(m2d, boundaries, 0, 0, 0)
psi_coil_zero = extract_flux_linkages(m2d)
psi_a_z = TURNS * (psi_coil_zero['PhaseA1'] + psi_coil_zero['PhaseA2'])
psi_b_z = TURNS * (psi_coil_zero['PhaseB1'] + psi_coil_zero['PhaseB2'])
psi_c_z = TURNS * (psi_coil_zero['PhaseC1'] + psi_coil_zero['PhaseC2'])
psi_alpha = 2/3*(psi_a_z - 0.5*psi_b_z - 0.5*psi_c_z)
psi_beta = 2/3*(np.sqrt(3)/2*psi_b_z - np.sqrt(3)/2*psi_c_z)
theta_r_meas = np.arctan2(psi_beta, psi_alpha)
phi_pm = np.sqrt(psi_alpha**2 + psi_beta**2)
print(f"  theta_r = {np.degrees(theta_r_meas):.1f} deg, Phi = {phi_pm:.6f} Wb")

# Phi effective from torque
Phi_eff = 207.6 / (1.5 * P * IMAX)
print(f"  Phi_eff = {Phi_eff:.6f} Wb (from T=T@0deg=207.6)")

# Run test points
print(f"\n[2] Test points:")
theta_r = np.radians(theta_r_meas)  # use measured value

for ms_angle, tr_thet, T_fea in test_angles:
    print(f"\n  MS angle={ms_angle}deg → TR Thet={tr_thet}deg (FEA T={T_fea:.1f} Nm)")

    # Set currents in MS
    Id_ms = IMAX * np.sin(np.radians(ms_angle))
    Iq_ms = IMAX * np.cos(np.radians(ms_angle))
    print(f"  Id_ms={Id_ms:.1f}, Iq_ms={Iq_ms:.1f}")

    Ia, Ib, Ic = id_iq_to_abc(Id_ms, Iq_ms, theta=theta_r)
    set_coil_currents(m2d, boundaries, Ia, Ib, Ic)

    psi_coil = extract_flux_linkages(m2d)
    psi_a = TURNS * (psi_coil['PhaseA1'] + psi_coil['PhaseA2'])
    psi_b = TURNS * (psi_coil['PhaseB1'] + psi_coil['PhaseB2'])
    psi_c = TURNS * (psi_coil['PhaseC1'] + psi_coil['PhaseC2'])

    # Park transform in MS frame
    psi_d, psi_q = park_abc_to_dq(psi_a, psi_b, psi_c, theta=theta_r)
    print(f"  psi_d(MS)={psi_d:.6f}, psi_q(MS)={psi_q:.6f}")

    # ---- CORRECTED APPROACH ----
    # The Ld, Lq computed in MS frame ARE the correct Ld, Lq for the motor,
    # because they represent the actual flux response to current in the dq frame.
    # BUT we need to be careful about WHICH (Id, Iq) we use.
    #
    # The key: when we set Id_ms, Iq_ms in the MS frame and measure psi_d, psi_q,
    # these ARE at the physical (Id_ms, Iq_ms) point.
    # The flaw was using Id_tr, Iq_tr (TR frame values) with MS-frame Ld/Lq.
    #
    # FIX: Compute Ld/Lq using MS-frame (Id_ms, Iq_ms) and psi_d, psi_q.
    # Then PREDICT TORQUE also using MS-frame (Id_ms, Iq_ms).
    # The TORQUE IS FRAME-INVARIANT! T = 1.5*P*(ψ_d*Iq - ψ_q*Id) works in any frame.

    # Method A: Co-energy formula (frame-invariant)
    T_coenergy = 1.5 * P * (psi_d * Iq_ms - psi_q * Id_ms)

    # Method B: L-based formula in MS frame
    Ld_ms = (psi_d - phi_pm) / Id_ms if abs(Id_ms) > 1e-6 else 0.0
    Lq_ms = psi_q / Iq_ms if abs(Iq_ms) > 1e-6 else 0.0
    T_L_ms = 1.5 * P * (phi_pm * Iq_ms + (Ld_ms - Lq_ms) * Id_ms * Iq_ms)

    print(f"  Ld(MS frame)={Ld_ms:.6f}, Lq(MS frame)={Lq_ms:.6f}")
    print(f"  T_coenergy = {T_coenergy:.1f} Nm  (frame-invariant)")
    print(f"  T_L_ms     = {T_L_ms:.1f} Nm  (MS-frame Ld/Lq)")
    print(f"  T_FEA      = {T_fea:.1f} Nm")
    print(f"  Error (coenergy) = {abs(T_coenergy - T_fea):.1f} Nm")

m2d.close_project()
print("\nDone.")
