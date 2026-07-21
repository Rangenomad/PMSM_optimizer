"""PMSM dq convention fix - joint 2D grid search for correct Park transform parameters.

Combines:
- theta_r search: rotor d-axis angle used in Park transform
- MS-to-TR offset search: relationship between MS scan angle and TR Thet_deg
- Also allows Phi to be re-calibrated

Uses raw flux linkage data from 4 test points to find the set of parameters
that makes torque formula match FEA across ALL points.
"""

import numpy as np

# ── Raw data from Magnetostatic ──
# Zero-current PM flux
PSI_PM = {
    'a': 0.128260, 'b': 0.128220, 'c': 0.250542,
}

# Test point data: {angle_deg: (psi_a, psi_b, psi_c, Id_set, Iq_set)}
RAW = {
    0:  (0.5816, -0.3583, 0.1276,   0.0, 350.0),
    30: (0.5732, -0.3092, 0.1859, 175.0, 303.1),
    50: (0.5111, -0.2195, 0.2675, 268.1, 225.0),
    60: (0.3334, -0.1205, 0.3595, 303.1, 175.0),
}

# ── FEA torque truth data (1-period rough values) ──
FEA_TRUTH = {
    0:  207.6,
    10: 254.5,
    20: 296.0,
    30: 331.7,
    40: 358.2,
    50: 370.4,
    60: 357.1,
    70: 300.2,
}

# ── Standard Park transform ──
def park_abc_to_dq(psi_a, psi_b, psi_c, theta):
    c = np.cos(theta)
    s = np.sin(theta)
    psi_d = 2/3 * (psi_a*c + psi_b*np.cos(theta-2*np.pi/3) + psi_c*np.cos(theta+2*np.pi/3))
    psi_q = -2/3 * (psi_a*s + psi_b*np.sin(theta-2*np.pi/3) + psi_c*np.sin(theta+2*np.pi/3))
    return psi_d, psi_q


def id_iq_to_abc(Id, Iq, theta):
    Ialpha = Id * np.cos(theta) - Iq * np.sin(theta)
    Ibeta  = Id * np.sin(theta) + Iq * np.cos(theta)
    Ia = Ialpha
    Ib = -0.5 * Ialpha + np.sqrt(3)/2 * Ibeta
    Ic = -0.5 * Ialpha - np.sqrt(3)/2 * Ibeta
    return Ia, Ib, Ic


# ── Phi calibration from zero-current point ──
psi_a0, psi_b0, psi_c0 = PSI_PM['a'], PSI_PM['b'], PSI_PM['c']
psi_alpha = 2/3 * (psi_a0 - 0.5*psi_b0 - 0.5*psi_c0)
psi_beta  = 2/3 * (np.sqrt(3)/2*psi_b0 - np.sqrt(3)/2*psi_c0)
theta_r_zero = np.arctan2(psi_beta, psi_alpha)
phi_pm_raw = np.sqrt(psi_alpha**2 + psi_beta**2)

print(f"Zero-current calibration:")
print(f"  psi_alpha={psi_alpha:.6f}, psi_beta={psi_beta:.6f}")
print(f"  theta_r_zero = {np.degrees(theta_r_zero):.1f} deg")
print(f"  phi_pm_raw   = {phi_pm_raw:.6f} Wb")

# Phi implied by T@0deg:  T = 1.5*P*Φ*350  →  Φ = T/(6*350)
P = 4
IMAX = 350.0
Phi_eff = FEA_TRUTH[0] / (1.5 * P * IMAX)
print(f"  Phi_eff (from T@0deg) = {Phi_eff:.6f} Wb")
print()

# ── Joint 2D search grid ──
# Parameters: theta_r, offset (MS angle → TR Thet_deg), use_Phi_eff

print("=" * 90)
print("Joint 2D search: theta_r × MS-to-TR offset")
print("  theta_r: rotor d-axis angle for Park transform (varied -180 to 180)")
print("  offset:  ms_angle + offset = TR Thet_deg (varied -90 to 90)")
print("  Score = mean |T_pred - T_FEA| across 4 test points")
print("=" * 90)

best_score = 1e9
best_config = None
best_details = []

for tr_deg in np.linspace(-180, 175, 72):  # 5-deg steps
    tr = np.radians(tr_deg)
    for offset in range(-90, 95, 5):
        for phi_use, phi_label in [(phi_pm_raw, 'raw_Phi'), (Phi_eff, 'eff_Phi')]:
            score = 0
            n_pts = 0

            for ms_angle, (psi_a, psi_b, psi_c, Id_set, Iq_set) in RAW.items():
                # Park transform with candidate theta_r
                psi_d, psi_q = park_abc_to_dq(psi_a, psi_b, psi_c, theta=tr)

                # Map MS angle to TR Thet_deg
                Thet_deg = ms_angle + offset
                if Thet_deg not in FEA_TRUTH:
                    continue

                # Compute Id, Iq in TR frame (standard Thet convention)
                Id_tr = IMAX * np.sin(np.radians(Thet_deg))
                Iq_tr = IMAX * np.cos(np.radians(Thet_deg))

                # BUT: the flux pattern in MS was set with Id_set, Iq_set
                # The question is: does the flux ψ_d, ψ_q correspond to
                # the same Id_set, Iq_set in the Park-transformed frame?

                # Approach 1: Use Id_set, Iq_set directly (same physical point)
                # This is what Sub-flow A does
                Ld_set = (psi_d - phi_use) / Id_set if abs(Id_set) > 1e-6 else 0.0
                Lq_set = psi_q / Iq_set if abs(Iq_set) > 1e-6 else 0.0
                T_pred_set = 1.5 * P * (phi_use * Iq_set + (Ld_set - Lq_set) * Id_set * Iq_set)

                # Approach 2: Use Id_tr, Iq_tr (TR frame decomposition)
                Ld_tr = (psi_d - phi_use) / Id_tr if abs(Id_tr) > 1e-6 else 0.0
                Lq_tr = psi_q / Iq_tr if abs(Iq_tr) > 1e-6 else 0.0
                T_pred_tr = 1.5 * P * (phi_use * Iq_tr + (Ld_tr - Lq_tr) * Id_tr * Iq_tr)

                T_fea = FEA_TRUTH[Thet_deg]
                err_set = abs(T_pred_set - T_fea)
                err_tr = abs(T_pred_tr - T_fea)
                score += min(err_set, err_tr)
                n_pts += 1

            if n_pts == len(RAW) and score < best_score:
                best_score = score
                best_config = (tr_deg, offset, phi_label)
                # Compute detailed results
                best_details = []
                for ms_angle, (psi_a, psi_b, psi_c, Id_set, Iq_set) in RAW.items():
                    psi_d, psi_q = park_abc_to_dq(psi_a, psi_b, psi_c, theta=tr)
                    Thet_deg = ms_angle + offset
                    Id_tr = IMAX * np.sin(np.radians(Thet_deg))
                    Iq_tr = IMAX * np.cos(np.radians(Thet_deg))

                    Ld_set = (psi_d - phi_use) / Id_set if abs(Id_set) > 1e-6 else 0.0
                    Lq_set = psi_q / Iq_set if abs(Iq_set) > 1e-6 else 0.0
                    T_pred_set = 1.5 * P * (phi_use * Iq_set + (Ld_set - Lq_set) * Id_set * Iq_set)

                    Ld_tr = (psi_d - phi_use) / Id_tr if abs(Id_tr) > 1e-6 else 0.0
                    Lq_tr = psi_q / Iq_tr if abs(Iq_tr) > 1e-6 else 0.0
                    T_pred_tr = 1.5 * P * (phi_use * Iq_tr + (Ld_tr - Lq_tr) * Id_tr * Iq_tr)

                    T_fea = FEA_TRUTH.get(Thet_deg, None)
                    best_details.append({
                        'ms_angle': ms_angle, 'Thet_deg': Thet_deg,
                        'psi_d': psi_d, 'psi_q': psi_q,
                        'Id_set': Id_set, 'Iq_set': Iq_set,
                        'Id_tr': Id_tr, 'Iq_tr': Iq_tr,
                        'Ld_set': Ld_set, 'Lq_set': Lq_set, 'T_set': T_pred_set,
                        'Ld_tr': Ld_tr, 'Lq_tr': Lq_tr, 'T_tr': T_pred_tr,
                        'T_fea': T_fea,
                    })

print(f"\nBest configuration:")
print(f"  theta_r (Park)  = {best_config[0]:.1f} deg")
print(f"  MS→TR offset     = {best_config[1]:+d} deg")
print(f"  Phi used         = {best_config[2]}")
print(f"  Avg abs error    = {best_score/len(RAW):.1f} Nm")
print(f"  Total error      = {best_score:.1f} Nm")

print(f"\nDetailed results for best configuration:")
print(f"  {'MS':>3s} {'TR':>4s} {'psi_d':>8s} {'psi_q':>8s} {'Ld_set':>10s} {'Lq_set':>10s} {'T_set':>8s} {'T_fea':>8s} {'err':>7s}")
print(f"  {'-'*75}")
for d in best_details:
    err = abs(d['T_set'] - d['T_fea'])
    print(f"  {d['ms_angle']:3d} {d['Thet_deg']:4d} "
          f"{d['psi_d']:8.4f} {d['psi_q']:8.4f} "
          f"{d['Ld_set']:+10.6f} {d['Lq_set']:+10.6f} "
          f"{d['T_set']:8.1f} {d['T_fea']:8.1f} {err:7.1f}")

# Also show best configurations sorted by score
print(f"\n--- Top 10 configurations ---")
# Re-do the top-N collection
configs = []
for tr_deg in np.linspace(-180, 175, 72):
    tr = np.radians(tr_deg)
    for offset in range(-90, 95, 5):
        for phi_use, phi_label in [(phi_pm_raw, 'raw'), (Phi_eff, 'eff')]:
            score = 0
            n_pts = 0
            for ms_angle, (psi_a, psi_b, psi_c, Id_set, Iq_set) in RAW.items():
                psi_d, psi_q = park_abc_to_dq(psi_a, psi_b, psi_c, theta=tr)
                Thet_deg = ms_angle + offset
                if Thet_deg not in FEA_TRUTH:
                    continue
                Id_tr = IMAX * np.sin(np.radians(Thet_deg))
                Iq_tr = IMAX * np.cos(np.radians(Thet_deg))
                Ld_set = (psi_d - phi_use)/Id_set if abs(Id_set)>1e-6 else 0.0
                Lq_set = psi_q/Iq_set if abs(Iq_set)>1e-6 else 0.0
                T_pred_set = 1.5*P*(phi_use*Iq_set+(Ld_set-Lq_set)*Id_set*Iq_set)
                Ld_tr = (psi_d - phi_use)/Id_tr if abs(Id_tr)>1e-6 else 0.0
                Lq_tr = psi_q/Iq_tr if abs(Iq_tr)>1e-6 else 0.0
                T_pred_tr = 1.5*P*(phi_use*Iq_tr+(Ld_tr-Lq_tr)*Id_tr*Iq_tr)
                T_fea = FEA_TRUTH[Thet_deg]
                score += min(abs(T_pred_set-T_fea), abs(T_pred_tr-T_fea))
                n_pts += 1
            if n_pts == len(RAW):
                configs.append((score, tr_deg, offset, phi_label))

configs.sort()
for i, (s, tr, off, pl) in enumerate(configs[:10]):
    # Get detailed Ld/Lq for this config
    tr_rad = np.radians(tr)
    phi = phi_pm_raw if pl == 'raw' else Phi_eff
    ld_summary = []
    for ms_angle, (psi_a, psi_b, psi_c, Id_set, Iq_set) in RAW.items():
        psi_d, psi_q = park_abc_to_dq(psi_a, psi_b, psi_c, theta=tr_rad)
        Ld = (psi_d-phi)/Id_set if abs(Id_set)>1e-6 else 0.0
        Lq = psi_q/Iq_set if abs(Iq_set)>1e-6 else 0.0
        ld_summary.append(f"Ld={Ld:.5f}")
    ld_str = ', '.join(ld_summary)
    print(f"  #{i+1}: theta_r={tr:6.1f}deg  offset={off:+3d}deg  Phi={pl:>3s}  score={s:6.1f}  {ld_str}")

# ── Additional test: try negative Id convention ──
print(f"\n--- Testing NEGATIVE Id convention (Id = -Imax*sin(theta)) ---")
# Standard: Id = Imax*sin(Thet), Iq = Imax*cos(Thet)
# Some conventions: Id = -Imax*sin(Thet), Iq = Imax*cos(Thet)  (Thet>0 gives negative Id)
best_neg_score = 1e9
best_neg_config = None
for tr_deg in np.linspace(-180, 175, 72):
    tr = np.radians(tr_deg)
    for offset in range(-90, 95, 5):
        for phi_use, phi_label in [(phi_pm_raw, 'raw'), (Phi_eff, 'eff')]:
            score = 0
            n_pts = 0
            for ms_angle, (psi_a, psi_b, psi_c, Id_set, Iq_set) in RAW.items():
                psi_d, psi_q = park_abc_to_dq(psi_a, psi_b, psi_c, theta=tr)
                Thet_deg = ms_angle + offset
                if Thet_deg not in FEA_TRUTH:
                    continue
                # NEGATIVE Id convention
                Id_tr = -IMAX * np.sin(np.radians(Thet_deg))  # flipped sign
                Iq_tr = IMAX * np.cos(np.radians(Thet_deg))
                # Use Id_set sign-matched
                Id_ms = -Id_set  # also flip MS Id convention
                Ld = (psi_d - phi_use)/Id_ms if abs(Id_ms)>1e-6 else 0.0
                Lq = psi_q/Iq_set if abs(Iq_set)>1e-6 else 0.0
                T_pred = 1.5*P*(phi_use*Iq_set+(Ld-Lq)*Id_ms*Iq_set)
                T_fea = FEA_TRUTH[Thet_deg]
                score += abs(T_pred - T_fea)
                n_pts += 1
            if n_pts == len(RAW) and score < best_neg_score:
                best_neg_score = score
                best_neg_config = (tr_deg, offset, phi_label)

print(f"  Best: theta_r={best_neg_config[0]:.1f}deg  offset={best_neg_config[1]:+d}deg  "
      f"Phi={best_neg_config[2]:>3s}  score={best_neg_score:.1f}")
