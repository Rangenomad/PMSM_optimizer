"""Fit Ld, Lq, Phi from FEA torque data and compute MTPA analytically.

Approach:
1. Use FEA torque at Thet=0 to calibrate Phi
2. For each Thet, compute DeltaL = Ld - Lq from torque formula
3. Use Lq from MS (Id=0, most reliable) as anchor
4. Fit Ld(Id) linear + quadratic saturation models
5. Compare constant-Ld, linear saturation, quadratic saturation MTPA vs FEA

Results validated against brute-force FEA scan (350A, 2000rpm):
  Thet=0  -> 207.6 Nm, Thet=10 -> 254.5 Nm, Thet=20 -> 296.0 Nm
  Thet=30 -> 331.7 Nm, Thet=40 -> 358.2 Nm, Thet=50 -> 370.4 Nm
  Thet=60 -> 357.1 Nm, Thet=70 -> 300.2 Nm

Key findings:
- Motor has REVERSE SALIENCY (Ld > Lq), MTPA at POSITIVE current angle
- Phi_eff (0.0989 Wb) differs from MS Phi (0.0815 Wb) due to load saturation
- Quadratic saturation model: angle error 1.0deg, torque error ~3% (vs 24.5deg constant)
"""

import numpy as np

# --- Motor constants ---
P = 4           # pole pairs
IMAX = 350.0    # A

# --- FEA torque data (1-period rough values, validated against 3-period accurate) ---
FEA_DATA = [
    (0,   207.6),
    (10,  254.5),
    (20,  296.0),   # 3-period accurate: 296.1
    (30,  331.7),
    (40,  358.2),
    (50,  370.4),   # 3-period accurate: 371.1
    (60,  357.1),
    (70,  300.2),
]

# --- Step 1: Calibrate Phi from Thet=0 ---
T0 = FEA_DATA[0][1]
Phi_eff = T0 / (1.5 * P * IMAX)
print("=" * 70)
print("Parameter Fitting from FEA Torque Data")
print("=" * 70)
print(f"\n[1] Phi calibration from Thet=0 deg:")
print(f"    T0 = {T0:.1f} Nm")
print(f"    Phi_eff = T0 / (1.5*P*{IMAX}) = {T0:.1f} / {1.5*P*IMAX:.0f} = {Phi_eff:.6f} Wb")
print(f"    (cf. MS measured Phi = ~0.0815 Wb)")

# --- Step 2: Compute DeltaL = Ld - Lq from each point ---
print(f"\n[2] DeltaL = Ld - Lq from each FEA point:")
print(f"    {'Thet':>5s}  {'Id':>7s}  {'Iq':>7s}  {'T_FEA':>7s}  {'(Ld-Lq)':>10s}  {'Ld (if Lq=0.00155)':>15s}")
print(f"    {'-'*60}")

delta_L_data = []
Lq_anchor = 0.00155  # from MS at Id=0 (most reliable MS data point)

for thet_deg, T in FEA_DATA:
    if thet_deg == 0:
        delta_L_data.append((0, 0, 0))
        continue

    Id = IMAX * np.sin(np.radians(thet_deg))
    Iq = IMAX * np.cos(np.radians(thet_deg))
    dL = (T / (1.5 * P) - Phi_eff * Iq) / (Id * Iq)
    Ld = Lq_anchor + dL
    delta_L_data.append((thet_deg, Id, dL))
    print(f"    {thet_deg:5.1f}  {Id:7.1f}  {Iq:7.1f}  {T:7.1f}  {dL:+10.6f}  {Ld:15.6f}")

# --- Step 3: Analyze Ld saturation trend ---
print(f"\n[3] Ld vs Id trend (assuming Lq={Lq_anchor:.5f} H constant):")
ids = [d[1] for d in delta_L_data if d[0] != 0]
lds = [Lq_anchor + d[2] for d in delta_L_data if d[0] != 0]
for i in range(len(ids)):
    print(f"    Id={ids[i]:6.1f} A  ->  Ld={lds[i]:.6f} H")

# Linear fit
from numpy import polyfit
coeffs_lin = polyfit(ids, lds, 1)
Ld_slope_lin, Ld_intercept = coeffs_lin
print(f"\n    Linear fit:  Ld = {Ld_slope_lin:.8f}*Id + {Ld_intercept:.6f}")

# Quadratic fit
coeffs_quad = polyfit(ids, lds, 2)
c2, c1q, c0q_ld = coeffs_quad
print(f"    Quadratic fit: Ld = {c2:.10f}*Id^2 + {c1q:.8f}*Id + {c0q_ld:.6f}")

# --- Step 4: Constant-Ld MTPA ---
print(f"\n[4] Constant-Ld MTPA (Ld={Ld_intercept:.6f} H, intercept):")
dL_const = Ld_intercept - Lq_anchor
print(f"    dL = Ld - Lq = {dL_const:.6f} H ({'Ld > Lq (reverse saliency)' if dL_const > 0 else 'Ld < Lq (normal IPM)'})")

if abs(dL_const) > 1e-12:
    a = 2 * dL_const * IMAX
    b = Phi_eff
    c = -dL_const * IMAX
    disc = b*b - 4*a*c
    if disc >= 0:
        s1 = (-b + np.sqrt(disc)) / (2*a)
        s2 = (-b - np.sqrt(disc)) / (2*a)
        beta_candidates = []
        for s in [s1, s2]:
            if -1 <= s <= 1:
                beta = np.degrees(np.arcsin(s))
                Id_val = IMAX * np.sin(np.radians(beta))
                Iq_val = IMAX * np.cos(np.radians(beta))
                T_val = 1.5 * P * (Phi_eff*Iq_val + dL_const*Id_val*Iq_val)
                beta_candidates.append((beta, s, Id_val, Iq_val, T_val))
        for beta, s, Id_val, Iq_val, T_val in beta_candidates:
            print(f"    Solution: sinB={s:.4f}, B={beta:.1f}deg, T_pred={T_val:.1f} Nm")
        best = max(beta_candidates, key=lambda x: x[4])
        beta_const, _, _, _, T_const = best
        print(f"    *** MTPA: Thet={beta_const:.1f}deg, T={T_const:.1f} Nm ***")

# --- Step 5: Linear saturation MTPA (sweep-based) ---
print(f"\n[5] Linear saturation MTPA (Ld_slope={Ld_slope_lin:.8f}):")
best_T_lin, best_a_lin = 0, 0
for a in np.arange(0, 70, 0.25):
    Id = IMAX * np.sin(np.radians(a))
    Iq = IMAX * np.cos(np.radians(a))
    dL_eff = Ld_slope_lin * Id + (Ld_intercept - Lq_anchor)
    T = 1.5 * P * (Phi_eff*Iq + dL_eff*Id*Iq)
    if T > best_T_lin:
        best_T_lin, best_a_lin = T, a
print(f"    *** MTPA: Thet={best_a_lin:.1f}deg, T={best_T_lin:.1f} Nm ***")

# --- Step 6: Quadratic saturation MTPA (sweep-based) ---
print(f"\n[6] Quadratic saturation MTPA:")
print(f"    Ld(Id) = {c0q_ld:.6f} + ({c1q:.2e})*Id + ({c2:.2e})*Id^2")
best_T_q, best_a_q = 0, 0
for a in np.arange(0, 70, 0.1):
    Id = IMAX * np.sin(np.radians(a))
    Iq = IMAX * np.cos(np.radians(a))
    dL_eff = c2*Id*Id + c1q*Id + (c0q_ld - Lq_anchor)
    T = 1.5 * P * (Phi_eff*Iq + dL_eff*Id*Iq)
    if T > best_T_q:
        best_T_q, best_a_q = T, a
print(f"    *** MTPA: Thet={best_a_q:.1f}deg, T={best_T_q:.1f} Nm ***")

# --- Step 7: Compare all models vs FEA ---
print(f"\n[7] MTPA comparison (all models vs FEA brute-force):")
print(f"    Brute-force FEA:         Thet=50.0deg, T=370.4 Nm")
print(f"    Constant Ld (intercept): Thet={beta_const:.1f}deg, T={T_const:.1f} Nm  (err: {abs(beta_const-50):.1f}deg)")
print(f"    Linear saturation:       Thet={best_a_lin:.1f}deg, T={best_T_lin:.1f} Nm  (err: {abs(best_a_lin-50):.1f}deg)")
print(f"    Quadratic saturation:    Thet={best_a_q:.1f}deg, T={best_T_q:.1f} Nm  (err: {abs(best_a_q-50):.1f}deg)")

# --- Step 8: Torque verification for each model ---
print(f"\n[8] Torque verification vs FEA at all points:")
print(f"    {'Thet':>5s}  {'T_FEA':>7s}  {'T_const':>8s}  {'T_lin':>8s}  {'T_quad':>8s}  {'err_q':>6s}")
print(f"    {'-'*52}")
total_const = 0; total_lin = 0; total_quad = 0
for thet_deg, T_fea in FEA_DATA:
    Id = IMAX * np.sin(np.radians(thet_deg))
    Iq = IMAX * np.cos(np.radians(thet_deg))
    T_c = 1.5*P*(Phi_eff*Iq + dL_const*Id*Iq)
    T_l = 1.5*P*(Phi_eff*Iq + (Ld_slope_lin*Id + Ld_intercept - Lq_anchor)*Id*Iq)
    T_q = 1.5*P*(Phi_eff*Iq + (c2*Id*Id + c1q*Id + c0q_ld - Lq_anchor)*Id*Iq)
    e_c = abs(T_c-T_fea); e_l = abs(T_l-T_fea); e_q = abs(T_q-T_fea)
    total_const += e_c; total_lin += e_l; total_quad += e_q
    print(f"    {thet_deg:5.1f}  {T_fea:7.1f}  {T_c:8.1f}  {T_l:8.1f}  {T_q:8.1f}  {e_q:6.1f}")
print(f"    Total abs error: const={total_const:.1f}  linear={total_lin:.1f}  quadratic={total_quad:.1f} Nm")

# --- Step 9: Summary ---
print(f"\n[9] Recommended parameters for subflow_d_efficiency_map:")
print(f"    DEFAULT_Ld = {c0q_ld:.6f}      # quadratic intercept")
print(f"    DEFAULT_Lq = {Lq_anchor:.6f}      # MS anchor")
print(f"    DEFAULT_Phi = {Phi_eff:.6f}     # FEA torque-calibrated")
print(f"    DEFAULT_Ld_slope = {c1q:.8f}   # linear saturation coef")
print(f"    DEFAULT_Ld_slope2 = {c2:.10f}   # quadratic saturation coef")
print(f"    MTPA accuracy: ~{abs(best_a_q-50):.1f}deg angle, ~{abs(best_T_q-370.4):.1f} Nm torque")
print(f"    Improvement: {total_const/total_quad:.1f}x over constant-Ld model")

print("\nDone.")
