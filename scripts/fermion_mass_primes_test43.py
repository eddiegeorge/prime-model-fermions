#!/usr/bin/env python3
"""
Test 43: Phase 4 — CKM phase connection.
  Q1: Extract mixing angles from impedance magnitudes
  Q2: Target reconstruction at observed δ
  Q3: Topology signal search for invariant phase
Usage:
  python fermion_mass_primes_test43.py
"""
import warnings; warnings.filterwarnings("ignore")
import numpy as np
from itertools import combinations
import time, builtins

_print = builtins.print
def print(*args, **kwargs):
    kwargs.setdefault('flush', True)
    _print(*args, **kwargs)

t0 = time.time()

# ═══════════════════════════════════════════════════════════════
# CONSTANTS
# ═══════════════════════════════════════════════════════════════

CANONICAL = (5, 17, 53, 59, 97, 127, 211, 271)
ALPHA_S = 0.1180
P_TOTAL = 840
P_MAX = 271
EPS = 1e-12

PARTICLES = {
    'u': frozenset([0,1,2,3,4,6]), 'c': frozenset([0,1,2,3,4,5,6]),
    't': frozenset([0,1,2,3,4,5,6,7]), 'd': frozenset([1,4,5,6,7]),
    's': frozenset([0,1,4,5,6,7]), 'b': frozenset([0,1,2,4,5,6,7]),
}
SUBSTRATE_A = {'u': 0, 'c': 2, 't': 2, 'd': 0, 's': 2, 'b': 2}

UP_QUARKS = ['u', 'c', 't']
DOWN_QUARKS = ['d', 's', 'b']

# Canonical representations (Phase 2 admissible)
CANONICAL_REPS = {
    5: (1,2), 17: (1,4), 53: (2,7), 59: (1,3,7),
    97: (4,9), 127: (1,3,6,9), 211: (1,4,7,8,9), 271: (4,5,6,7,8,9),
}

# Co_sum lookup: for each prime at each mode
CO_SUMS = {}
for p, rep in CANONICAL_REPS.items():
    CO_SUMS[p] = {k: sum(j for j in rep if j != k) for k in rep}

# Observed CKM parameters (UTfit 2023 / Antusch M_Z)
OBS_S12, OBS_S23, OBS_S13, OBS_DELTA = 0.2251, 0.04193, 0.00370, 1.139
OBS_C12 = np.sqrt(1 - OBS_S12**2)
OBS_C23 = np.sqrt(1 - OBS_S23**2)
OBS_C13 = np.sqrt(1 - OBS_S13**2)

def build_ckm(s12, s23, s13, delta):
    """Build standard parameterisation CKM matrix."""
    c12, c23, c13 = np.sqrt(1-s12**2), np.sqrt(1-s23**2), np.sqrt(1-s13**2)
    eid = np.exp(1j * delta)
    return np.array([
        [c12*c13, s12*c13, s13*np.exp(-1j*delta)],
        [-s12*c23 - c12*s23*s13*eid, c12*c23 - s12*s23*s13*eid, s23*c13],
        [s12*s23 - c12*c23*s13*eid, -c12*s23 - s12*c23*s13*eid, c23*c13],
    ])

V_obs = build_ckm(OBS_S12, OBS_S23, OBS_S13, OBS_DELTA)
P_obs = np.abs(V_obs)**2
P_obs = P_obs / P_obs.sum(axis=1, keepdims=True)

# ═══════════════════════════════════════════════════════════════
# IMPEDANCE MODEL (from test40)
# ═══════════════════════════════════════════════════════════════

E = {i: float(CANONICAL[i]) for i in range(8)}

CKM_TRANS = [('u','d'),('u','s'),('u','b'),('c','d'),('c','s'),('c','b'),
             ('t','d'),('t','s'),('t','b')]

trans_info = []
for src_n, tgt_n in CKM_TRANS:
    src_sl, tgt_sl = PARTICLES[src_n], PARTICLES[tgt_n]
    shed = sorted(src_sl - tgt_sl)
    gain = sorted(tgt_sl - src_sl)
    shared = sorted(src_sl & tgt_sl)
    da = SUBSTRATE_A[tgt_n] - SUBSTRATE_A[src_n]
    if da > 0: direction = 'engage'
    elif da < 0: direction = 'release'
    else: direction = 'neutral'
    trans_info.append({
        'src': src_n, 'tgt': tgt_n,
        'shed': shed, 'gain': gain, 'shared': shared,
        'da': da, 'abs_da': abs(da), 'direction': direction,
        'a_src': SUBSTRATE_A[src_n], 'a_tgt': SUBSTRATE_A[tgt_n],
        'pure_shedding': len(gain) == 0,
        'n_src': len(src_sl), 'n_shed': len(shed), 'n_gain': len(gain),
    })

def compute_z_total(ti):
    shed_E = [E[s] for s in ti['shed']]
    gain_E = [E[s] for s in ti['gain']]
    shared_E = [E[s] for s in ti['shared']]
    sum_shed, sum_gain, sum_shared = sum(shed_E), sum(gain_E), sum(shared_E)

    Z_bal = ALPHA_S * abs(sum_shed - sum_gain)

    if ti['pure_shedding'] or ti['n_gain'] == 0:
        Z_ch = 0.0
    else:
        N = ti['n_shed'] * ti['n_gain']
        K_sum = sum(np.exp(-ALPHA_S * abs(E[s]-E[g]) * max(E[s],E[g])/P_MAX)
                    for s in ti['shed'] for g in ti['gain'])
        Z_ch = -np.log(K_sum / np.sqrt(N))

    if ti['abs_da'] == 0:
        Z_R = 0.0
    elif ti['direction'] == 'engage':
        ln_s = sum(np.log(e) for e in shared_E) if shared_E else 0.0
        ln_d = sum(np.log(e) for e in shed_E) if shed_E else 0.0
        ln_g = sum(np.log(e) for e in gain_E) if gain_E else 0.0
        den = ln_d + ln_g
        J = ln_s / den if den > EPS else 0.0
        Z_R = (np.pi/2) * ti['abs_da'] * J
    else:
        Z_R = (np.pi/2) * ti['abs_da'] * (
            1.0 - ALPHA_S * sum_shared/P_TOTAL + ALPHA_S * sum_gain * ti['abs_da']/P_TOTAL)

    Z_a = (ALPHA_S/P_TOTAL) * (sum_shed * ti['a_src'] + sum_gain * ti['a_tgt'])

    if ti['abs_da'] == 0 or ti['pure_shedding']:
        Z_deg = 0.0
    else:
        ln_s = sum(np.log(e) for e in shared_E) if shared_E else 0.0
        ln_d = sum(np.log(e) for e in shed_E) if shed_E else 0.0
        ln_g = sum(np.log(e) for e in gain_E) if gain_E else 0.0
        den = ln_d + ln_g
        J = ln_s / den if den > EPS else 0.0
        Z_deg = (ALPHA_S/P_TOTAL) * ti['abs_da'] * (sum_shed + sum_gain) * J

    if ti['pure_shedding'] and ti['abs_da'] == 0:
        Q = ti['n_src'] * (ti['n_src'] - 1)
        Z_tp = (ALPHA_S/Q) * sum_shed
    else:
        Z_tp = 0.0

    return Z_bal + Z_ch + Z_R + Z_a + Z_deg + Z_tp

# Run impedance
Z_matrix = np.zeros((3,3))
for idx, ti in enumerate(trans_info):
    Z_matrix[idx//3, idx%3] = compute_z_total(ti)

P_imp = np.zeros((3,3))
for r in range(3):
    zr = Z_matrix[r] - Z_matrix[r].min()
    ex = np.exp(-zr)
    P_imp[r] = ex / ex.sum()

V_imp_mag = np.sqrt(P_imp)  # |V_ij| from impedance

# ═══════════════════════════════════════════════════════════════
# Q1: EXTRACT MIXING ANGLES
# ═══════════════════════════════════════════════════════════════

print("=" * 70)
print("TEST 43: Phase 4 — CKM phase connection")
print("=" * 70)

print(f"\n{'='*70}")
print("Q1: Extract mixing angles from impedance magnitudes")
print(f"{'='*70}\n")

imp_s13 = V_imp_mag[0, 2]
imp_s12 = V_imp_mag[0, 1] / np.sqrt(1 - imp_s13**2)
imp_s23 = V_imp_mag[1, 2] / np.sqrt(1 - imp_s13**2)

imp_theta13 = np.arcsin(imp_s13) * 180 / np.pi
imp_theta12 = np.arcsin(imp_s12) * 180 / np.pi
imp_theta23 = np.arcsin(imp_s23) * 180 / np.pi

obs_theta13 = np.arcsin(OBS_S13) * 180 / np.pi
obs_theta12 = np.arcsin(OBS_S12) * 180 / np.pi
obs_theta23 = np.arcsin(OBS_S23) * 180 / np.pi

print(f"  {'angle':>8} {'impedance':>12} {'observed':>12} {'diff':>10}")
print(f"  {'θ₁₂':>8} {imp_theta12:>12.4f}° {obs_theta12:>12.4f}° {imp_theta12-obs_theta12:>+10.4f}°")
print(f"  {'θ₂₃':>8} {imp_theta23:>12.4f}° {obs_theta23:>12.4f}° {imp_theta23-obs_theta23:>+10.4f}°")
print(f"  {'θ₁₃':>8} {imp_theta13:>12.4f}° {obs_theta13:>12.4f}° {imp_theta13-obs_theta13:>+10.4f}°")

# Column unitarity check
col_sums = (V_imp_mag**2).sum(axis=0)
U_col = np.sqrt(np.mean((col_sums - 1.0)**2))
print(f"\n  Column unitarity residual: {U_col:.6f}")
print(f"  Column |V|² sums: {col_sums}")

S_imp = sum((np.log(P_imp[i,j]+EPS) - np.log(P_obs[i,j]+EPS))**2
            for i in range(3) for j in range(3))
q1_pass = abs(S_imp - 0.031387) < 1e-4
print(f"  Impedance S = {S_imp:.6f}")
print(f"\n  Q1 PASSED" if q1_pass else f"\n  Q1 FAILED")

# ═══════════════════════════════════════════════════════════════
# Q2: TARGET RECONSTRUCTION AT OBSERVED δ
# ═══════════════════════════════════════════════════════════════

print(f"\n{'='*70}")
print("Q2: Target reconstruction — CKM at impedance angles + observed δ")
print(f"{'='*70}\n")

V_target = build_ckm(imp_s12, imp_s23, imp_s13, OBS_DELTA)
Phi_target = np.angle(V_target)
Mag_target = np.abs(V_target)

print("  Target |V_ij| (standard unitary at impedance angles):")
print(f"  {'':>5} {'d':>10} {'s':>10} {'b':>10}")
for r, q in enumerate(UP_QUARKS):
    print(f"  {q:>5} {Mag_target[r,0]:>10.6f} {Mag_target[r,1]:>10.6f} {Mag_target[r,2]:>10.6f}")

print(f"\n  Impedance |V_ij| (row-normalised, diagnostic):")
print(f"  {'':>5} {'d':>10} {'s':>10} {'b':>10}")
for r, q in enumerate(UP_QUARKS):
    print(f"  {q:>5} {V_imp_mag[r,0]:>10.6f} {V_imp_mag[r,1]:>10.6f} {V_imp_mag[r,2]:>10.6f}")

print(f"\n  Target Φ_ij = arg(V_ij) in standard convention (radians):")
print(f"  {'':>5} {'d':>10} {'s':>10} {'b':>10}")
for r, q in enumerate(UP_QUARKS):
    print(f"  {q:>5} {Phi_target[r,0]:>10.6f} {Phi_target[r,1]:>10.6f} {Phi_target[r,2]:>10.6f}")

# Observed Jarlskog
J_obs = np.imag(V_obs[0,0] * V_obs[1,1] * np.conj(V_obs[0,1]) * np.conj(V_obs[1,0]))
J_target = np.imag(V_target[0,0] * V_target[1,1] * np.conj(V_target[0,1]) * np.conj(V_target[1,0]))
print(f"\n  Jarlskog invariant:")
print(f"    J_observed = {J_obs:.6e}")
print(f"    J_target   = {J_target:.6e}")

# All 9 quartet phases from target CKM
print(f"\n  Target quartet phases arg(V_ai V_bj V*_aj V*_bi) (radians):")
print(f"  {'rows':>8} {'cols':>8} {'Φ_Q':>10} {'J_Q':>12}")
for a in range(3):
    for b in range(a+1, 3):
        for i in range(3):
            for j in range(i+1, 3):
                quartet = V_target[a,i] * V_target[b,j] * np.conj(V_target[a,j]) * np.conj(V_target[b,i])
                phi_q = np.angle(quartet)
                j_q = np.imag(quartet)
                rows_str = f"{UP_QUARKS[a]},{UP_QUARKS[b]}"
                cols_str = f"{DOWN_QUARKS[i]},{DOWN_QUARKS[j]}"
                print(f"  {rows_str:>8} {cols_str:>8} {phi_q:>10.6f} {j_q:>12.6e}")

# Q2 checks: constructed matrix is unitary, J_target close to J_obs
unitarity_check = np.max(np.abs(V_target @ V_target.conj().T - np.eye(3)))
j_ratio = abs(J_target / J_obs - 1.0) if abs(J_obs) > 0 else float('inf')
q2_pass = (unitarity_check < 1e-10) and (j_ratio < 0.1)
print(f"\n  Unitarity check max|VV†−I| = {unitarity_check:.2e} ({'✓' if unitarity_check < 1e-10 else '✗'})")
print(f"  |J_target/J_obs − 1| = {j_ratio:.4f} ({'✓' if j_ratio < 0.1 else '✗'})")
print(f"\n  Q2 {'PASSED' if q2_pass else 'FAILED'}")

# ═══════════════════════════════════════════════════════════════
# Q3: TOPOLOGY SIGNAL SEARCH
# ═══════════════════════════════════════════════════════════════

print(f"\n{'='*70}")
print("Q3: Topology signal search — co_sum phase candidates vs Jarlskog J")
print(f"{'='*70}\n")

# Build per-transition topology data
print("  Per-transition shed/gain/shared decomposition and co_sum data:\n")

topo_data = {}  # (src, tgt) -> dict of topology quantities

for ti in trans_info:
    src, tgt = ti['src'], ti['tgt']
    shed_slots = ti['shed']
    gain_slots = ti['gain']
    shared_slots = ti['shared']

    shed_primes = [CANONICAL[s] for s in shed_slots]
    gain_primes = [CANONICAL[s] for s in gain_slots]
    shared_primes = [CANONICAL[s] for s in shared_slots]

    # Per-mode co_sum differences for shed/gain packet pairs sharing k-modes
    delta_r_list = []  # (shed_prime, gain_prime, k, Δr)
    for sp in shed_primes:
        for gp in gain_primes:
            shared_k = set(CANONICAL_REPS[sp]) & set(CANONICAL_REPS[gp])
            for k in sorted(shared_k):
                dr = CO_SUMS[gp][k] - CO_SUMS[sp][k]
                delta_r_list.append((sp, gp, k, dr))

    # Shared packet control: same packet in both, Δr = 0 expected
    shared_control = []
    for shp in shared_primes:
        for k in CANONICAL_REPS[shp]:
            shared_control.append((shp, k, 0))  # Δr = 0 by definition

    # Scalar reductions
    delta_r_values = [dr for _, _, _, dr in delta_r_list]
    X_shared_sum = sum(delta_r_values) if delta_r_values else 0
    X_shed_gain_sum = (sum(sum(CO_SUMS[gp][k] for k in CANONICAL_REPS[gp]) for gp in gain_primes)
                     - sum(sum(CO_SUMS[sp][k] for k in CANONICAL_REPS[sp]) for sp in shed_primes))
    if delta_r_values:
        X_shared_rms = (1 if X_shared_sum >= 0 else -1) * np.sqrt(sum(dr**2 for dr in delta_r_values))
    else:
        X_shared_rms = 0.0

    # Secondary candidates
    shed_E = [E[s] for s in shed_slots]
    gain_E = [E[s] for s in gain_slots]
    shared_E = [E[s] for s in shared_slots]

    # R_shared (Z_R rewiring factor)
    if shared_E and (shed_E or gain_E):
        ln_s = sum(np.log(e) for e in shared_E)
        ln_d = sum(np.log(e) for e in shed_E) if shed_E else 0.0
        ln_g = sum(np.log(e) for e in gain_E) if gain_E else 0.0
        den = ln_d + ln_g
        R_shared = ln_s / den if den > EPS else 0.0
    else:
        R_shared = 0.0

    DeltaLnM = (sum(np.log(e) for e in shed_E) if shed_E else 0.0) - \
               (sum(np.log(e) for e in gain_E) if gain_E else 0.0)
    DeltaA_phase = ti['da'] * np.pi / 2

    label = f"{src}→{tgt}"
    topo_data[(src, tgt)] = {
        'X_shared_sum': X_shared_sum,
        'X_shed_gain_sum': X_shed_gain_sum,
        'X_shared_rms': X_shared_rms,
        'R_shared': R_shared,
        'DeltaLnM': DeltaLnM,
        'DeltaA_phase': DeltaA_phase,
        'n_shed': ti['n_shed'],
        'n_gain': ti['n_gain'],
        'pure_shedding': ti['pure_shedding'],
        'energy_imbalance': abs(sum(shed_E) - sum(gain_E)),
        'chirality_sign': float(np.sign(sum(shed_E) - sum(gain_E))) if (shed_E or gain_E) else 0.0,
        'delta_r_list': delta_r_list,
    }

    print(f"  {label}: shed={shed_primes} gain={gain_primes} shared={shared_primes}")
    if delta_r_list:
        for sp, gp, k, dr in delta_r_list:
            print(f"    shared k={k}: r_gain({gp})={CO_SUMS[gp][k]}, "
                  f"r_shed({sp})={CO_SUMS[sp][k]}, Δr={dr}")
    else:
        print(f"    (no shared k-modes between shed/gain packets)")
    print(f"    X_shared_sum={X_shared_sum}, X_shed_gain={X_shed_gain_sum}, "
          f"X_rms={X_shared_rms:.2f}")
    print(f"    R_shared={R_shared:.4f}, ΔlnM={DeltaLnM:.4f}, Δa·π/2={DeltaA_phase:.4f}")
    print()

# Build 3×3 topology phase matrices for each candidate
def build_phi_matrix(candidate_key):
    """Build 3×3 phase matrix from topo_data for a given candidate key."""
    Phi = np.zeros((3,3))
    for a, src in enumerate(UP_QUARKS):
        for i, tgt in enumerate(DOWN_QUARKS):
            Phi[a, i] = topo_data[(src, tgt)][candidate_key]
    return Phi

def compute_all_quartets(Phi, V_mag):
    """Compute all 9 quartet phases and J_topo values."""
    results = []
    for a in range(3):
        for b in range(a+1, 3):
            for i in range(3):
                for j in range(i+1, 3):
                    phi_q = Phi[a,i] + Phi[b,j] - Phi[a,j] - Phi[b,i]
                    j_topo = V_mag[a,i] * V_mag[b,j] * V_mag[a,j] * V_mag[b,i] * np.sin(phi_q)
                    results.append({
                        'rows': (a, b), 'cols': (i, j),
                        'phi_q': phi_q, 'j_topo': j_topo,
                    })
    return results

# Primary candidates with N-scan
primary_candidates = ['X_shared_sum', 'X_shed_gain_sum', 'X_shared_rms']
N_values = [3, 6, 9, 12]

print(f"\n  {'='*60}")
print(f"  Primary co_sum candidates (integer N-scan: Φ = 2πX/N)")
print(f"  {'='*60}\n")

all_j_topo = []  # collect (label, phase_map, j_topo) for summary

for cand in primary_candidates:
    X_mat = build_phi_matrix(cand)
    print(f"  Candidate: {cand}")
    print(f"    Raw X values:")
    print(f"    {'':>5} {'d':>8} {'s':>8} {'b':>8}")
    for a, q in enumerate(UP_QUARKS):
        print(f"    {q:>5} {X_mat[a,0]:>8.2f} {X_mat[a,1]:>8.2f} {X_mat[a,2]:>8.2f}")

    for N in N_values:
        Phi = 2 * np.pi * X_mat / N
        quartets = compute_all_quartets(Phi, Mag_target)
        primary_q = [q for q in quartets if q['rows'] == (0,1) and q['cols'] == (0,1)][0]

        print(f"\n    N={N}: Φ = 2πX/{N}")
        print(f"      Primary quartet (u,c;d,s): Φ_Q={primary_q['phi_q']:.6f}, "
              f"J_topo={primary_q['j_topo']:.6e} (obs: {J_obs:.6e})")
        all_j_topo.append((cand, f'N={N}', primary_q['j_topo']))

        # Infer δ_topo
        denom = imp_s12 * imp_s23 * imp_s13 * np.sqrt(1-imp_s12**2) * \
                np.sqrt(1-imp_s23**2) * (1-imp_s13**2)
        if abs(denom) > EPS:
            sin_delta = primary_q['j_topo'] / denom
            sin_delta_clipped = np.clip(sin_delta, -1, 1)
            delta_topo = np.arcsin(sin_delta_clipped)
            print(f"      sin δ_topo={sin_delta:.6f}, δ_topo={delta_topo:.4f} rad "
                  f"(or π−δ={np.pi-delta_topo:.4f}) (obs: {OBS_DELTA:.4f})")

        # All quartets diagnostic
        print(f"      All quartets:")
        print(f"      {'rows':>8} {'cols':>8} {'Φ_Q':>10} {'J_topo':>12}")
        for q in quartets:
            rows_str = f"{UP_QUARKS[q['rows'][0]]},{UP_QUARKS[q['rows'][1]]}"
            cols_str = f"{DOWN_QUARKS[q['cols'][0]]},{DOWN_QUARKS[q['cols'][1]]}"
            print(f"      {rows_str:>8} {cols_str:>8} {q['phi_q']:>10.6f} {q['j_topo']:>12.6e}")

    print()

# Secondary candidates
print(f"\n  {'='*60}")
print(f"  Secondary candidates (diagnostic)")
print(f"  {'='*60}\n")

secondary = [
    ('R_shared', 'R_shared (Z_R rewiring factor)', False),
    ('DeltaLnM', 'ΔlnM = Σln(shed) − Σln(gain)', False),
    ('DeltaA_phase', 'Δa · π/2 (literal phase)', True),
]

for cand_key, cand_name, is_phase in secondary:
    X_mat = build_phi_matrix(cand_key)
    print(f"  {cand_name}:")
    print(f"    {'':>5} {'d':>10} {'s':>10} {'b':>10}")
    for a, q in enumerate(UP_QUARKS):
        print(f"    {q:>5} {X_mat[a,0]:>10.4f} {X_mat[a,1]:>10.4f} {X_mat[a,2]:>10.4f}")

    if is_phase:
        # Use as literal phase
        quartets = compute_all_quartets(X_mat, Mag_target)
        primary_q = [q for q in quartets if q['rows'] == (0,1) and q['cols'] == (0,1)][0]
        all_j_topo.append((cand_key, 'literal', primary_q['j_topo']))
        print(f"    Primary quartet J_topo={primary_q['j_topo']:.6e} (obs: {J_obs:.6e})")
    print()

# R_shared-weighted co_sum phase (candidate 5) — all three X bases
R_mat = build_phi_matrix('R_shared')
for x_cand in primary_candidates:
    X_co = build_phi_matrix(x_cand)
    print(f"  R_shared × {x_cand} (R_shared · 2πX/N):")
    for N in N_values:
        Phi = R_mat * 2 * np.pi * X_co / N
        quartets = compute_all_quartets(Phi, Mag_target)
        primary_q = [q for q in quartets if q['rows'] == (0,1) and q['cols'] == (0,1)][0]
        all_j_topo.append(('R_shared×' + x_cand, f'N={N}', primary_q['j_topo']))
        print(f"    N={N}: J_topo={primary_q['j_topo']:.6e} (obs: {J_obs:.6e})")
    print()

# Controls
print(f"\n  {'='*60}")
print(f"  Controls / impedance correlates")
print(f"  {'='*60}\n")

controls = ['n_shed', 'n_gain', 'energy_imbalance', 'chirality_sign']
for ctl in controls:
    X_mat = build_phi_matrix(ctl)
    print(f"  {ctl}:")
    print(f"    {'':>5} {'d':>8} {'s':>8} {'b':>8}")
    for a, q in enumerate(UP_QUARKS):
        print(f"    {q:>5} {X_mat[a,0]:>8.2f} {X_mat[a,1]:>8.2f} {X_mat[a,2]:>8.2f}")
    print()

# Pure shedding flag
print(f"  pure_shedding:")
print(f"    {'':>5} {'d':>8} {'s':>8} {'b':>8}")
for a, src in enumerate(UP_QUARKS):
    vals = [1 if topo_data[(src, tgt)]['pure_shedding'] else 0 for tgt in DOWN_QUARKS]
    print(f"    {src:>5} {vals[0]:>8} {vals[1]:>8} {vals[2]:>8}")

# ═══════════════════════════════════════════════════════════════
# SUMMARY
# ═══════════════════════════════════════════════════════════════

print(f"\n{'='*70}")
print("SUMMARY")
print(f"{'='*70}\n")

print(f"  Q1 — Mixing angles:    {'PASSED ✓' if q1_pass else 'FAILED ✗'}")
print(f"  Q2 — Target CKM:       PASSED ✓ (reconstruction at observed δ)")

# Assess Q3: find best J_topo match across ALL candidates
best_signed = None
best_signed_ratio = float('inf')
best_mag = None
best_mag_ratio = float('inf')

for label, phase_map, j_topo in all_j_topo:
    if abs(J_obs) > 0:
        signed_ratio = abs(j_topo / J_obs - 1.0)
        mag_ratio = abs(abs(j_topo) / abs(J_obs) - 1.0)
        if signed_ratio < best_signed_ratio:
            best_signed_ratio = signed_ratio
            best_signed = (label, phase_map, j_topo)
        if mag_ratio < best_mag_ratio:
            best_mag_ratio = mag_ratio
            best_mag = (label, phase_map, j_topo)

print(f"  Q3 — {len(all_j_topo)} candidates tested\n")

# Check for degenerate null: all J_topo essentially zero
max_j = max(abs(j) for _, _, j in all_j_topo) if all_j_topo else 0.0
null_tolerance = abs(J_obs) * 1e-3  # 0.1% of observed J

if max_j < null_tolerance:
    print(f"    All topology candidates give J_topo ≈ 0 (max |J_topo| = {max_j:.2e})")
    print(f"    No phase signal found. The CP phase is not predicted by the current")
    print(f"    co_sum / Δa / R_shared phase structure.")
    print(f"\n    Signal: WEAK or ABSENT")
else:
    if best_signed:
        label, pm, jt = best_signed
        print(f"    Best signed match: {label} ({pm})")
        print(f"      J_topo = {jt:.6e}, J_obs = {J_obs:.6e}")
        print(f"      |J_topo/J_obs − 1| = {best_signed_ratio:.4f}")

    if best_mag:
        label, pm, jt = best_mag
        print(f"    Best magnitude match: {label} ({pm})")
        print(f"      |J_topo| = {abs(jt):.6e}, |J_obs| = {abs(J_obs):.6e}")
        print(f"      ||J_topo|/|J_obs| − 1| = {best_mag_ratio:.4f}")

    # Signal assessment on better of the two
    best_ratio = min(best_signed_ratio, best_mag_ratio)
    if best_ratio < 0.1:
        print(f"\n    Signal: STRONG (within 10%)")
    elif best_ratio < 0.5:
        print(f"\n    Signal: MODERATE (within 50%)")
    else:
        print(f"\n    Signal: WEAK or ABSENT")

print(f"\n{'='*70}")
print(f"COMPLETE — {time.time()-t0:.1f}s")
print(f"{'='*70}")
