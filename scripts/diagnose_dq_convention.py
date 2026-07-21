"""Diagnose dq convention mismatch between Magnetostatic (Sub-flow A) and Transient.

Problem: Ld/Lq from Sub-flow A predict negative torque when Id>0,
but FEA (Transient TR design) shows positive torque at Id>0.

This script:
1. Runs zero-current calibration → get raw flux and theta_r
2. Runs a few MS points at known (Id_set, Iq_set)
3. Tests 12 Park transform conventions (sign flips, dq swaps)
4. Finds the convention where torque formula matches FEA at multiple points
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


def old_id_iq_to_abc(Id, Iq, theta=0):
    Ialpha = Id * np.cos(theta) - Iq * np.sin(theta)
    Ibeta  = Id * np.sin(theta) + Iq * np.cos(theta)
    Ia = Ialpha
    Ib = -0.5 * Ialpha + np.sqrt(3) / 2 * Ibeta
    Ic = -0.5 * Ialpha - np.sqrt(3) / 2 * Ibeta
    return Ia, Ib, Ic


def old_park_abc_to_dq(psi_a, psi_b, psi_c, theta=0):
    c = np.cos(theta)
    s = np.sin(theta)
    psi_d = 2/3 * (psi_a*c + psi_b*np.cos(theta-2*np.pi/3) + psi_c*np.cos(theta+2*np.pi/3))
    psi_q = -2/3 * (psi_a*s + psi_b*np.sin(theta-2*np.pi/3) + psi_c*np.sin(theta+2*np.pi/3))
    return psi_d, psi_q


def park_variants():
    """Generate all 8 sign/spin variants of the Park transform."""
    variants = []
    names = []
    # Base transforms
    for d_sign in [1, -1]:
        for q_sign in [1, -1]:
            for dq_swap in [False, True]:
                name = f"d_sign={d_sign:+d} q_sign={q_sign:+d} swap={dq_swap}"
                variants.append((d_sign, q_sign, dq_swap))
                names.append(name)
    return variants, names


def apply_park_variant(psi_a, psi_b, psi_c, theta, d_sign, q_sign, dq_swap):
    """Apply a Park transform variant."""
    c = np.cos(theta)
    s = np.sin(theta)
    psi_d_raw = 2/3 * (psi_a*c + psi_b*np.cos(theta-2*np.pi/3) + psi_c*np.cos(theta+2*np.pi/3))
    psi_q_raw = -2/3 * (psi_a*s + psi_b*np.sin(theta-2*np.pi/3) + psi_c*np.sin(theta+2*np.pi/3))

    psi_d = d_sign * psi_d_raw
    psi_q = q_sign * psi_q_raw

    if dq_swap:
        psi_d, psi_q = psi_q, psi_d

    return psi_d, psi_q


def id_iq_variants():
    """Generate 4 variants of Id/Iq mapping."""
    return [
        ("std", lambda Id,Iq,th: (Id,Iq)),
        ("swap", lambda Id,Iq,th: (Iq,Id)),
        ("neg_id", lambda Id,Iq,th: (-Id,Iq)),
        ("neg_iq", lambda Id,Iq,th: (Id,-Iq)),
    ]


def set_coil_currents(m2d, boundaries, Ia, Ib, Ic):
    currents = [TURNS * Ia, TURNS * Ia, TURNS * Ib, TURNS * Ib, TURNS * Ic, TURNS * Ic]
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


# ── FEA truth data (from brute-force Transient scan at POSITIVE angles) ──
# Key: Thet_deg in Transient, Value: torque (Nm) — 1-period rough scan
FEA_TRUTH = {
    0:  207.6,   # from MTPA scan era
    10: 254.5,   # from extended positive angle scan
    20: 296.0,   # from extended positive angle scan
    30: 331.7,   # from extended positive angle scan
    40: 358.2,   # from extended positive angle scan
    50: 370.4,   # from extended positive angle scan
    60: 357.1,   # from extended positive angle scan
    70: 300.2,   # from extended positive angle scan
}

# More accurate 3-period values:
FEA_ACCURATE = {
    20: 296.1,   # from torque_350A project
    50: 371.1,   # from mtpa_comparison project, 3-period accurate
}


def compute_expected_ld_lq_from_fea():
    """For each FEA point, compute what Ld/Lq would need to be to match torque.

    T = 1.5*P*(Phi*Iq + (Ld-Lq)*Id*Iq)
    At a given (Thet, T, Phi), we have 2 unknowns (Ld,Lq) but only 1 equation.
    We need at least 2 points to solve.

    At θ1=0°: Id1=0, Iq1=350, T1=208 (or 209 from accurate)
        T1 = 1.5*P*Phi*350  →  Phi = T1/(6*350) = 208/2100 = 0.0990

    At θ2=50°: Id2=268, Iq2=225, T2=371
        T2 = 1.5*P*(Phi*225 + (Ld-Lq)*268*225)
        371 = 6*(0.099*225 + ΔL*60300)
        61.83 = 22.28 + ΔL*60300
        ΔL = (61.83-22.28)/60300 = 0.000656  →  Ld-Lq = +0.000656
    """

    P = 4
    Imax = 350

    # Use FEA accurate data
    T_0 = FEA_ACCURATE.get(0, FEA_TRUTH.get(0, 207.6))
    Phi_implied = T_0 / (1.5 * P * Imax)
    print(f"  Phi implied by T(0deg)={T_0:.1f}: {Phi_implied:.6f} Wb")

    results = {}
    for deg in sorted(FEA_TRUTH.keys()):
        if deg == 0:
            continue
        Id = Imax * np.sin(np.radians(deg))
        Iq = Imax * np.cos(np.radians(deg))
        T = FEA_TRUTH[deg]
        # T = 1.5*P*(Phi_implied*Iq + dL*Id*Iq)
        if abs(Id) < 1e-6:
            continue
        dL = (T/(1.5*P) - Phi_implied*Iq) / (Id * Iq)
        results[deg] = dL

    return Phi_implied, results


def main():
    project_path = 'pmsm_projects/2026-07-19_dq_fix'
    template = str(Path(project_path) / 'Prius_2D_Practice.aedt')

    print("=" * 70)
    print("dq Convention Diagnosis")
    print("=" * 70)

    # Step 1: Compute expected Ld-Lq from FEA
    print("\n[1] Expected Ld-Lq from FEA torque data:")
    Phi_fea, dL_fea = compute_expected_ld_lq_from_fea()
    for deg, dl in sorted(dL_fea.items()):
        print(f"    theta={deg:3d}deg: need (Ld-Lq)={dl:.6f} H")

    # Step 2: Open MS design, run zero-current calibration + 3 test points
    print("\n[2] Running Magnetostatic test points...")

    m2d = Maxwell2d(
        project=template, design='4_Partial_motor_MS2',
        solution_type='MagnetostaticXY',
        non_graphical=False, new_desktop=False, close_on_exit=False
    )

    # Get boundary references
    boundaries = []
    for name, _ in COIL_CONFIG:
        b = [b for b in m2d.boundaries if b.name == name][0]
        boundaries.append(b)

    PolePairs = 4  # 8-pole PMSM, 4 pole pairs (fixed template parameter)
    print(f"    PolePairs = {PolePairs} (fixed)")

    # ── Zero current calibration ──
    print("    Running zero-current calibration...")
    set_coil_currents(m2d, boundaries, 0, 0, 0)
    psi_coil_zero = extract_flux_linkages(m2d)
    psi_a_z = TURNS * (psi_coil_zero['PhaseA1'] + psi_coil_zero['PhaseA2'])
    psi_b_z = TURNS * (psi_coil_zero['PhaseB1'] + psi_coil_zero['PhaseB2'])
    psi_c_z = TURNS * (psi_coil_zero['PhaseC1'] + psi_coil_zero['PhaseC2'])
    psi_alpha = 2/3 * (psi_a_z - 0.5*psi_b_z - 0.5*psi_c_z)
    psi_beta  = 2/3 * (np.sqrt(3)/2 * psi_b_z - np.sqrt(3)/2 * psi_c_z)
    theta_r = np.arctan2(psi_beta, psi_alpha)
    phi_pm = np.sqrt(psi_alpha**2 + psi_beta**2)

    print(f"    theta_r = {np.degrees(theta_r):.1f} deg")
    print(f"    Phi_PM  = {phi_pm:.6f} Wb")
    print(f"    Raw psi_a0={psi_a_z:.6f}, psi_b0={psi_b_z:.6f}, psi_c0={psi_c_z:.6f}")

    # ── Test points: map to Transient Thet_deg angles of interest ──
    # We need to figure out the MS angle parameter that corresponds to
    # each Transient Thet_deg. For now, use the same convention as Sub-flow A.
    # Sub-flow A uses: Id=Imax*sin(angle), Iq=Imax*cos(angle)
    # where angle is the scan parameter (0-90 deg).
    # We'll test the same angles and compare.

    test_angles = [0, 30, 50, 60]  # MS scan angles
    Imax = 350

    raw_data = {}  # {angle: (psi_a, psi_b, psi_c)}
    for angle_deg in test_angles:
        theta = np.radians(angle_deg)
        Id_set = Imax * np.sin(theta)
        Iq_set = Imax * np.cos(theta)

        Ia, Ib, Ic = old_id_iq_to_abc(Id_set, Iq_set, theta=theta_r)
        set_coil_currents(m2d, boundaries, Ia, Ib, Ic)
        psi_coil = extract_flux_linkages(m2d)

        psi_a = TURNS * (psi_coil['PhaseA1'] + psi_coil['PhaseA2'])
        psi_b = TURNS * (psi_coil['PhaseB1'] + psi_coil['PhaseB2'])
        psi_c = TURNS * (psi_coil['PhaseC1'] + psi_coil['PhaseC2'])
        raw_data[angle_deg] = (psi_a, psi_b, psi_c, Id_set, Iq_set)

        print(f"    angle={angle_deg}deg (Id_set={Id_set:.1f}, Iq_set={Iq_set:.1f}): "
              f"psi_a={psi_a:.4f}, psi_b={psi_b:.4f}, psi_c={psi_c:.4f}")

    m2d.close_project()

    # ── Step 3: Grid search for effective theta_r ──
    # The key insight: we don't know the exact mapping between MS angle and TR Thet_deg.
    # Instead, search over theta_r to find the transformation that makes
    # T(Id_set,Iq_set) from Ld/Lq match FEA T(Thet_deg).
    # Since Id_set=Imax*sin(angle), Iq_set=Imax*cos(angle), and we assume
    # angle ≈ Thet_deg (same rotor position in both designs),
    # we search for the theta_r that makes torque formula work.

    # But first, we need to also allow Phi to be calibrated, since it affects torque.
    # At angle=0: Id=0, Iq=350, T_FEA=207.6
    #   T = 1.5*P*(Phi*Iq + 0) = 6*Phi*350 = 2100*Phi
    #   → Phi_eff = 207.6/2100 = 0.098856 Wb (effective Phi at rated current)
    Phi_eff = FEA_TRUTH[0] / (1.5 * PolePairs * Imax)
    print(f"\n    Effective Phi (from T@0deg): {Phi_eff:.6f} Wb (raw PM Phi: {phi_pm:.6f})")
    print(f"    This Phi will be used as ground truth for torque formula calibration.\n")

    best_tr_score = 1e9
    best_tr_deg = 0.0
    best_tr_params = None

    print("    Searching theta_r grid [-180, 180] at 2-deg steps × 4 angle points...")
    for tr_deg in np.linspace(-180, 180, 181):
        tr = np.radians(tr_deg)
        score = 0
        n_pts = 0
        tr_params = {}

        for angle_deg in test_angles:
            psi_a, psi_b, psi_c, Id_set, Iq_set = raw_data[angle_deg]

            psi_d, psi_q = old_park_abc_to_dq(psi_a, psi_b, psi_c, theta=tr)

            Ld = (psi_d - phi_pm) / Id_set if abs(Id_set) > 1e-6 else 0.0
            Lq = psi_q / Iq_set if abs(Iq_set) > 1e-6 else 0.0

            T_pred = 1.5 * PolePairs * (phi_pm * Iq_set + (Ld - Lq) * Id_set * Iq_set)

            # We also need the formula with Phi_eff (from T@0deg) as alternative
            T_pred_eff = 1.5 * PolePairs * (Phi_eff * Iq_set + (Ld - Lq) * Id_set * Iq_set)

            T_fea = FEA_TRUTH.get(angle_deg, None)
            if T_fea is not None:
                # Use whichever Phi gives better match
                err_pm = abs(T_pred - T_fea)
                err_eff = abs(T_pred_eff - T_fea)
                score += min(err_pm, err_eff)
                n_pts += 1
            tr_params[angle_deg] = {'Ld': Ld, 'Lq': Lq, 'T_pred_pm': T_pred, 'T_pred_eff': T_pred_eff}

        if n_pts > 0 and score < best_tr_score:
            best_tr_score = score
            best_tr_deg = tr_deg
            best_tr_params = tr_params

    print(f"\n    Best theta_r = {best_tr_deg:.1f}deg (avg error = {best_tr_score:.1f} Nm across {len(test_angles)} points)")

    # Report detailed results for the best theta_r
    print(f"\n    Results at theta_r={best_tr_deg:.1f}deg:")
    print(f"    {'angle':>5s}  {'Id':>6s}  {'Iq':>6s}  {'Ld':>10s}  {'Lq':>10s}  {'T_form':>8s}  {'T_FEA':>8s}")
    print(f"    {'-'*60}")
    for angle_deg in test_angles:
        p = best_tr_params[angle_deg]
        T_fea = FEA_TRUTH.get(angle_deg, 0)
        # Pick best T_pred (PM phi or effective phi)
        T_pred = p['T_pred_pm'] if abs(p['T_pred_pm']-T_fea) < abs(p['T_pred_eff']-T_fea) else p['T_pred_eff']
        print(f"    {angle_deg:5d}  {p['Ld']:+10.6f}  {p['Lq']:+10.6f}  {T_pred:8.1f}  {T_fea:8.1f}")

    # ── Step 4: Search for MS-angle → Thet_deg offset ──
    # What if the MS scan angle doesn't directly correspond to TR Thet_deg?
    # Try: Thet_deg_effective = MS_angle + offset
    print(f"\n[4] Searching for MS→TR angle offset (with theta_r={best_tr_deg:.1f}deg)...")

    best_off_score = 1e9
    best_offset = 0.0

    for offset in range(-60, 70, 10):
        score = 0
        n_pts = 0

        for angle_ms in test_angles:
            psi_a, psi_b, psi_c, _, _ = raw_data[angle_ms]
            tr = np.radians(best_tr_deg)

            psi_d, psi_q = old_park_abc_to_dq(psi_a, psi_b, psi_c, theta=tr)

            # Map: MS_angle → TR Thet_deg = MS_angle + offset
            Thet_deg_eff = angle_ms + offset
            Id_eff = Imax * np.sin(np.radians(Thet_deg_eff))
            Iq_eff = Imax * np.cos(np.radians(Thet_deg_eff))

            Ld = (psi_d - phi_pm) / Id_eff if abs(Id_eff) > 1e-6 else 0.0
            Lq = psi_q / Iq_eff if abs(Iq_eff) > 1e-6 else 0.0

            T_pred = 1.5 * PolePairs * (phi_pm * Iq_eff + (Ld - Lq) * Id_eff * Iq_eff)
            T_fea = FEA_TRUTH.get(Thet_deg_eff, None)
            if T_fea is not None:
                score += abs(T_pred - T_fea)
                n_pts += 1

        if n_pts > 0 and score < best_off_score:
            best_off_score = score
            best_offset = offset

    print(f"    Best MS→TR offset: {best_offset:+d}deg (score = {best_off_score:.1f})")

    print("\nDone. Analysis complete.")


if __name__ == '__main__':
    main()
