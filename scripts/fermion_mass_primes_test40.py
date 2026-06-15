#!/usr/bin/env python3
"""
Test 40: Phase 1 — Cavity-grounded impedance model.
  Stage 1:  v2.0 arithmetic control (Gate A: S ≈ 0.031)
  Stage 2:  Ontological translation — E_p = Σk² confirms E_p = p (Gate B: identical)
Usage:
  python fermion_mass_primes_test40.py              # both stages
  python fermion_mass_primes_test40.py --control     # stage 1 only
  python fermion_mass_primes_test40.py --cavity      # stage 2 only
"""
import warnings; warnings.filterwarnings("ignore")
import numpy as np
from itertools import combinations
from sympy import isprime
import time, argparse

parser = argparse.ArgumentParser(description="Test 40: Cavity-grounded impedance")
parser.add_argument('--control', action='store_true', help='Stage 1 only')
parser.add_argument('--cavity', action='store_true', help='Stage 2 only')
args = parser.parse_args()

RUN_CONTROL = not args.cavity
RUN_CAVITY = not args.control

# ═══════════════════════════════════════════════════════════════
# CONSTANTS
# ═══════════════════════════════════════════════════════════════

CANONICAL = (5, 17, 53, 59, 97, 127, 211, 271)
ALPHA_S = 0.1180
P_TOTAL = 840  # Σp over codebook
P_MAX = 271    # max p over codebook

PARTICLES = {
    'u': frozenset([0,1,2,3,4,6]), 'c': frozenset([0,1,2,3,4,5,6]),
    't': frozenset([0,1,2,3,4,5,6,7]), 'd': frozenset([1,4,5,6,7]),
    's': frozenset([0,1,4,5,6,7]), 'b': frozenset([0,1,2,4,5,6,7]),
}
SUBSTRATE_A = {'u': 0, 'c': 2, 't': 2, 'd': 0, 's': 2, 'b': 2}

CKM_TRANS = [('u','d'),('u','s'),('u','b'),('c','d'),('c','s'),('c','b'),
             ('t','d'),('t','s'),('t','b')]
EXPECTED_SIGS = ['d>s>b', 's>d>b', 'b>s>d']
EPS = 1e-12

# ═══════════════════════════════════════════════════════════════
# OBSERVED CKM
# ═══════════════════════════════════════════════════════════════

s12, s23, s13, delta_ckm = 0.2251, 0.04193, 0.00370, 1.139
c12, c23, c13 = np.sqrt(1-s12**2), np.sqrt(1-s23**2), np.sqrt(1-s13**2)
V_obs = np.array([
    [c12*c13, s12*c13, s13*np.exp(-1j*delta_ckm)],
    [-s12*c23-c12*s23*s13*np.exp(1j*delta_ckm),
      c12*c23-s12*s23*s13*np.exp(1j*delta_ckm), s23*c13],
    [s12*s23-c12*c23*s13*np.exp(1j*delta_ckm),
     -c12*s23-s12*c23*s13*np.exp(1j*delta_ckm), c23*c13],
])
P_obs = np.abs(V_obs)**2
P_obs = P_obs / P_obs.sum(axis=1, keepdims=True)

# ═══════════════════════════════════════════════════════════════
# TRANSITION STRUCTURE
# ═══════════════════════════════════════════════════════════════

trans_info = []
for src_n, tgt_n in CKM_TRANS:
    src_sl, tgt_sl = PARTICLES[src_n], PARTICLES[tgt_n]
    shed = sorted(src_sl - tgt_sl)
    gain = sorted(tgt_sl - src_sl)
    shared = sorted(src_sl & tgt_sl)
    da = SUBSTRATE_A[tgt_n] - SUBSTRATE_A[src_n]
    a_src = SUBSTRATE_A[src_n]
    a_tgt = SUBSTRATE_A[tgt_n]
    if da > 0: direction = 'engage'
    elif da < 0: direction = 'release'
    else: direction = 'neutral'
    pure_shedding = len(gain) == 0
    n_src = len(src_sl)
    trans_info.append({
        'src': src_n, 'tgt': tgt_n,
        'shed': shed, 'gain': gain, 'shared': shared,
        'da': da, 'abs_da': abs(da), 'direction': direction,
        'a_src': a_src, 'a_tgt': a_tgt,
        'pure_shedding': pure_shedding,
        'n_src': n_src,
        'n_shed': len(shed), 'n_gain': len(gain),
    })

# ═══════════════════════════════════════════════════════════════
# Z-TERM COMPUTATION
# ═══════════════════════════════════════════════════════════════

def compute_z_terms(ti, E, P_norm, E_max):
    """
    Compute all six Z-terms for one transition.
    E: dict mapping slot index → packet energy
    P_norm: total energy normalisation
    E_max: max energy over codebook
    """
    shed_E = [E[s] for s in ti['shed']]
    gain_E = [E[s] for s in ti['gain']]
    shared_E = [E[s] for s in ti['shared']]

    sum_shed = sum(shed_E)
    sum_gain = sum(gain_E)
    sum_shared = sum(shared_E)

    # --- Z_bal: energy balance ---
    Z_bal = ALPHA_S * abs(sum_shed - sum_gain)

    # --- Z_ch: channel conductance ---
    if ti['pure_shedding'] or ti['n_gain'] == 0:
        Z_ch = 0.0
    else:
        N = ti['n_shed'] * ti['n_gain']
        K_sum = 0.0
        for s in ti['shed']:
            for g in ti['gain']:
                dist = abs(E[s] - E[g])
                pos = max(E[s], E[g]) / E_max
                K_sum += np.exp(-ALPHA_S * dist * pos)
        K = K_sum / np.sqrt(N)
        Z_ch = -np.log(K)

    # --- Z_R: substrate change ---
    if ti['abs_da'] == 0:
        Z_R = 0.0
    elif ti['direction'] == 'engage':
        ln_shared = sum(np.log(e) for e in shared_E) if shared_E else 0.0
        ln_shed = sum(np.log(e) for e in shed_E) if shed_E else 0.0
        ln_gain = sum(np.log(e) for e in gain_E) if gain_E else 0.0
        denom = ln_shed + ln_gain
        J = ln_shared / denom if denom > EPS else 0.0
        Z_R = (np.pi / 2.0) * ti['abs_da'] * J
    else:
        # releasing
        sum_rewired = sum_shared
        back_emf = sum_gain * ti['abs_da']
        Z_R = (np.pi / 2.0) * ti['abs_da'] * (
            1.0 - ALPHA_S * sum_rewired / P_norm
                + ALPHA_S * back_emf / P_norm
        )

    # --- Z_a: substrate throughput ---
    Z_a = (ALPHA_S / P_norm) * (sum_shed * ti['a_src'] + sum_gain * ti['a_tgt'])

    # --- Z_deg: channel degradation ---
    if ti['abs_da'] == 0 or ti['pure_shedding']:
        Z_deg = 0.0
    else:
        ln_shared = sum(np.log(e) for e in shared_E) if shared_E else 0.0
        ln_shed = sum(np.log(e) for e in shed_E) if shed_E else 0.0
        ln_gain = sum(np.log(e) for e in gain_E) if gain_E else 0.0
        denom = ln_shed + ln_gain
        J = ln_shared / denom if denom > EPS else 0.0
        Z_deg = (ALPHA_S / P_norm) * ti['abs_da'] * (sum_shed + sum_gain) * J

    # --- Z_tp: ring traversal (pure shedding at Δa=0 only) ---
    if ti['pure_shedding'] and ti['abs_da'] == 0:
        Q = ti['n_src'] * (ti['n_src'] - 1)
        Z_tp = (ALPHA_S / Q) * sum_shed
    else:
        Z_tp = 0.0

    Z_total = Z_bal + Z_ch + Z_R + Z_a + Z_deg + Z_tp
    return {
        'Z_bal': Z_bal, 'Z_ch': Z_ch, 'Z_R': Z_R,
        'Z_a': Z_a, 'Z_deg': Z_deg, 'Z_tp': Z_tp,
        'Z_total': Z_total,
    }

# ═══════════════════════════════════════════════════════════════
# SCORING
# ═══════════════════════════════════════════════════════════════

def compute_probabilities(Z_matrix):
    """Row-normalised softmax: P_ij = exp(-Z_ij) / Σ_j exp(-Z_ij).
    Numerically stable: shift each row by its minimum Z before exponentiating."""
    P = np.zeros((3, 3))
    for r in range(3):
        row_Z = Z_matrix[r]
        row_Z_shifted = row_Z - row_Z.min()
        row_exp = np.exp(-row_Z_shifted)
        row_sum = row_exp.sum()
        P[r] = row_exp / row_sum if row_sum > EPS else row_exp
    return P

def score_ckm(P):
    return sum((np.log(P[i,j]+EPS) - np.log(P_obs[i,j]+EPS))**2
               for i in range(3) for j in range(3))

def row_scores(P):
    u = sum((np.log(P[0,j]+EPS) - np.log(P_obs[0,j]+EPS))**2 for j in range(3))
    c = sum((np.log(P[1,j]+EPS) - np.log(P_obs[1,j]+EPS))**2 for j in range(3))
    t = sum((np.log(P[2,j]+EPS) - np.log(P_obs[2,j]+EPS))**2 for j in range(3))
    return u, c, t

def score_sig(P):
    tgts = ['d','s','b']
    sigs, exp_count = [], 0
    for r in range(3):
        order = sorted(range(3), key=lambda c: P[r,c], reverse=True)
        sig = '>'.join(tgts[c] for c in order)
        sigs.append(sig)
        if sig == EXPECTED_SIGS[r]: exp_count += 1
    return '|'.join(sigs), exp_count

# ═══════════════════════════════════════════════════════════════
# RUN ONE STAGE
# ═══════════════════════════════════════════════════════════════

def run_stage(label, E, P_norm, E_max):
    """Run the impedance model with given energy values and scales."""
    print(f"\n{'='*70}")
    print(f"{label}")
    print(f"{'='*70}")

    print(f"\n  Packet energies:")
    for i, p in enumerate(CANONICAL):
        print(f"    slot {i}: prime {p:>3}, E = {E[i]:.4f}")
    print(f"  P_norm = {P_norm:.4f}, E_max = {E_max:.4f}")

    Z_matrix = np.zeros((3, 3))
    all_z = []
    print(f"\n  Per-transition Z-term breakdown:")
    print(f"  {'trans':<5} {'Z_bal':>7} {'Z_ch':>7} {'Z_R':>7} {'Z_a':>7} "
          f"{'Z_deg':>7} {'Z_tp':>7} {'Z_tot':>7}")

    for ti_idx, ti in enumerate(trans_info):
        r, c = ti_idx // 3, ti_idx % 3
        z = compute_z_terms(ti, E, P_norm, E_max)
        Z_matrix[r, c] = z['Z_total']
        all_z.append(z)
        label_t = f"{ti['src']}→{ti['tgt']}"
        print(f"  {label_t:<5} {z['Z_bal']:>7.4f} {z['Z_ch']:>7.4f} {z['Z_R']:>7.4f} "
              f"{z['Z_a']:>7.4f} {z['Z_deg']:>7.4f} {z['Z_tp']:>7.4f} {z['Z_total']:>7.4f}")

    P = compute_probabilities(Z_matrix)
    S = score_ckm(P)
    u_sc, c_sc, t_sc = row_scores(P)
    sig, exp = score_sig(P)

    rows = ['u','c','t']
    print(f"\n  Predicted P_ij:")
    print(f"  {'':>5} {'d':>10} {'s':>10} {'b':>10}")
    for r in range(3):
        print(f"  {rows[r]:>5} {P[r,0]:>10.6f} {P[r,1]:>10.6f} {P[r,2]:>10.6f}")

    print(f"\n  Observed P_ij:")
    print(f"  {'':>5} {'d':>10} {'s':>10} {'b':>10}")
    for r in range(3):
        print(f"  {rows[r]:>5} {P_obs[r,0]:>10.6f} {P_obs[r,1]:>10.6f} {P_obs[r,2]:>10.6f}")

    print(f"\n  Per-element log-error² (P_pred vs P_obs):")
    print(f"  {'':>5} {'d':>10} {'s':>10} {'b':>10}")
    for r in range(3):
        errs = [(np.log(P[r,c]+EPS) - np.log(P_obs[r,c]+EPS))**2 for c in range(3)]
        print(f"  {rows[r]:>5} {errs[0]:>10.6f} {errs[1]:>10.6f} {errs[2]:>10.6f}")

    print(f"\n  S = {S:.6f}")
    print(f"  Per-row: u={u_sc:.4f}, c={c_sc:.4f}, t={t_sc:.4f}")
    print(f"  Signature: {sig} ({exp}/3)")

    return S, u_sc, c_sc, t_sc, sig, exp

# ═══════════════════════════════════════════════════════════════
# CAVITY E_p DERIVATION
# ═══════════════════════════════════════════════════════════════

def derive_cavity_energies():
    """Compute E_p = Σk² over occupied modes for each codebook prime.
    Confirms E_p = p identically."""
    prime_reps = {}
    for r in range(1, 10):
        for subset in combinations(range(1, 10), r):
            energy = sum(k*k for k in subset)
            if isprime(energy):
                prime_reps.setdefault(energy, []).append(subset)

    print(f"\n  Cavity E_p derivation (E_p = Σk² over occupied modes):")
    E_cavity = {}
    for i, p in enumerate(CANONICAL):
        reps = prime_reps.get(p, [])
        if not reps:
            raise ValueError(f"No cavity representation found for prime {p}")

        rep_energies = [sum(k*k for k in rep) for rep in reps]
        if any(E != p for E in rep_energies):
            raise ValueError(f"Bad representation for prime {p}: energies={rep_energies}")

        E_cavity[i] = float(rep_energies[0])

        modes_str = "; ".join(
            "{" + ",".join(str(k) for k in rep) + "}" +
            f" → Σk²={sum(k*k for k in rep)}"
            for rep in reps[:3]
        )

        print(f"    slot {i}: prime {p:>3}, reps: {modes_str}")
        print(f"             E_p = {E_cavity[i]:.0f} = p ✓")

    return E_cavity

# ═══════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════

t0 = time.time()

print("=" * 70)
print("TEST 40: Phase 1 — Cavity-grounded impedance model")
print("=" * 70)

results = []

# ─────────────────────────────────────────────────────────────
# STAGE 1: v2.0 arithmetic control
# ─────────────────────────────────────────────────────────────

if RUN_CONTROL:
    E_v20 = {i: float(CANONICAL[i]) for i in range(8)}
    S1, u1, c1, t1, sig1, exp1 = run_stage(
        "STAGE 1: v2.0 arithmetic control",
        E_v20, float(P_TOTAL), float(P_MAX))
    results.append(('v2.0 control', S1, u1, c1, t1, sig1, exp1))

# ─────────────────────────────────────────────────────────────
# STAGE 2: Ontological translation
# ─────────────────────────────────────────────────────────────

if RUN_CAVITY:
    print(f"\n{'='*70}")
    print(f"CAVITY ONTOLOGICAL GROUNDING")
    print(f"{'='*70}")

    E_cavity = derive_cavity_energies()

    E_cav_sum = sum(E_cavity[i] for i in range(8))
    E_cav_max = max(E_cavity[i] for i in range(8))
    print(f"\n  Σ E_p = {E_cav_sum:.0f} (should be {P_TOTAL})")
    print(f"  max E_p = {E_cav_max:.0f} (should be {P_MAX})")

    print(f"\n  Parameter grounding:")
    print(f"    P = {P_TOTAL} = Σ E_p over codebook (derived, not free)")
    print(f"    p_max = {P_MAX} = max E_p over codebook (derived, not free)")
    print(f"    α_s = {ALPHA_S} = strong coupling constant at M_Z (measured SM constant)")
    print(f"    π/2 = {np.pi/2:.4f} = geometric phase for orthogonal substrate (open: not yet derived from cavity)")

    S2, u2, c2, t2, sig2, exp2 = run_stage(
        "STAGE 2: Cavity ontological translation (E_p = Σk² = p)",
        E_cavity, E_cav_sum, E_cav_max)
    results.append(('Cavity E_p', S2, u2, c2, t2, sig2, exp2))

    # Confirm identity
    if RUN_CONTROL:
        S_diff = abs(results[0][1] - results[1][1])
        print(f"\n  Gate B check: |S_control - S_cavity| = {S_diff:.2e}")
        if S_diff < 1e-10:
            print(f"  GATE B PASSED: ontological translation reproduces v2.0 identically.")
        else:
            print(f"  GATE B FAILED: discrepancy detected.")

# ═══════════════════════════════════════════════════════════════
# COMPARISON TABLE
# ═══════════════════════════════════════════════════════════════

print(f"\n{'='*70}")
print("COMPARISON TABLE")
print(f"{'='*70}\n")

print(f"  {'Model':<25} {'S':>8} {'u_sc':>7} {'c_sc':>7} {'t_sc':>7} {'Sig':>5}")
for name, S, u, c, t, sig, exp in results:
    print(f"  {name:<25} {S:>8.4f} {u:>7.4f} {c:>7.4f} {t:>7.4f} {exp:>3}/3")

print(f"\n{'='*70}")
print(f"COMPLETE — {time.time()-t0:.1f}s")
print(f"{'='*70}")
