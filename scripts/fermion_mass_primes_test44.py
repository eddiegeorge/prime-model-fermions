#!/usr/bin/env python3
"""
Test 44: Phase 5 — Z-term ablation.
  Q1: Per-term contribution profile (9×6 matrix)
  Q2: Single-term ablation (6 removals)
  Q3: Pairwise ablation with interaction (15 pairs)
  Q4: Full subset ablation (63 subsets + 6 single-term-only)
Usage:
  python fermion_mass_primes_test44.py
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
CKM_TRANS = [('u','d'),('u','s'),('u','b'),('c','d'),('c','s'),('c','b'),
             ('t','d'),('t','s'),('t','b')]
EXPECTED_SIGS = ['d>s>b', 's>d>b', 'b>s>d']

TERM_NAMES = ['Z_bal', 'Z_ch', 'Z_R', 'Z_a', 'Z_deg', 'Z_tp']

# Observed CKM
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

E = {i: float(CANONICAL[i]) for i in range(8)}

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

# ═══════════════════════════════════════════════════════════════
# COMPUTE INDIVIDUAL Z-TERMS (returns 6-element array per transition)
# ═══════════════════════════════════════════════════════════════

def compute_z_terms(ti):
    """Return array of [Z_bal, Z_ch, Z_R, Z_a, Z_deg, Z_tp] for one transition."""
    shed_E = [E[s] for s in ti['shed']]
    gain_E = [E[s] for s in ti['gain']]
    shared_E = [E[s] for s in ti['shared']]
    sum_shed, sum_gain, sum_shared = sum(shed_E), sum(gain_E), sum(shared_E)

    # Z_bal
    Z_bal = ALPHA_S * abs(sum_shed - sum_gain)

    # Z_ch
    if ti['pure_shedding'] or ti['n_gain'] == 0:
        Z_ch = 0.0
    else:
        N = ti['n_shed'] * ti['n_gain']
        K_sum = sum(np.exp(-ALPHA_S * abs(E[s]-E[g]) * max(E[s],E[g])/P_MAX)
                    for s in ti['shed'] for g in ti['gain'])
        Z_ch = -np.log(K_sum / np.sqrt(N))

    # Z_R
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

    # Z_a
    Z_a = (ALPHA_S/P_TOTAL) * (sum_shed * ti['a_src'] + sum_gain * ti['a_tgt'])

    # Z_deg
    if ti['abs_da'] == 0 or ti['pure_shedding']:
        Z_deg = 0.0
    else:
        ln_s = sum(np.log(e) for e in shared_E) if shared_E else 0.0
        ln_d = sum(np.log(e) for e in shed_E) if shed_E else 0.0
        ln_g = sum(np.log(e) for e in gain_E) if gain_E else 0.0
        den = ln_d + ln_g
        J = ln_s / den if den > EPS else 0.0
        Z_deg = (ALPHA_S/P_TOTAL) * ti['abs_da'] * (sum_shed + sum_gain) * J

    # Z_tp
    if ti['pure_shedding'] and ti['abs_da'] == 0:
        Q = ti['n_src'] * (ti['n_src'] - 1)
        Z_tp = (ALPHA_S/Q) * sum_shed
    else:
        Z_tp = 0.0

    return np.array([Z_bal, Z_ch, Z_R, Z_a, Z_deg, Z_tp])

# Precompute all Z-terms for all 9 transitions
Z_all = np.zeros((9, 6))  # [transition_idx, term_idx]
for idx, ti in enumerate(trans_info):
    Z_all[idx] = compute_z_terms(ti)

# ═══════════════════════════════════════════════════════════════
# SCORING FUNCTIONS
# ═══════════════════════════════════════════════════════════════

def score_from_z_matrix(Z_mat):
    """Compute S, per-row S, P matrix, and row signatures from a 3×3 Z matrix."""
    P = np.zeros((3,3))
    for r in range(3):
        zr = Z_mat[r] - Z_mat[r].min()
        ex = np.exp(-zr)
        P[r] = ex / ex.sum()

    S_rows = np.zeros(3)
    for r in range(3):
        S_rows[r] = sum((np.log(P[r,c]+EPS) - np.log(P_obs[r,c]+EPS))**2 for c in range(3))

    S_total = S_rows.sum()

    sigs = []
    for r in range(3):
        order = np.argsort(-P[r])
        sig = '>'.join(DOWN_QUARKS[k] for k in order)
        sigs.append(sig)

    return S_total, S_rows, P, sigs

def evaluate_term_mask(mask):
    """Evaluate impedance model with only the terms indicated by mask (6-bit).
    mask[k] = True means term k is active."""
    Z_mat = np.zeros((3,3))
    for idx in range(9):
        r, c = idx // 3, idx % 3
        Z_mat[r, c] = sum(Z_all[idx, k] for k in range(6) if mask[k])
    return score_from_z_matrix(Z_mat)

# ═══════════════════════════════════════════════════════════════
# BASELINE
# ═══════════════════════════════════════════════════════════════

print("=" * 70)
print("TEST 44: Phase 5 — Z-term ablation")
print("=" * 70)

baseline_mask = [True] * 6
S_base, S_rows_base, P_base, sigs_base = evaluate_term_mask(baseline_mask)

print(f"\n  Baseline: all 6 terms active")
print(f"  S_baseline = {S_base:.6f}")
print(f"  Per-row: S_u={S_rows_base[0]:.6f}, S_c={S_rows_base[1]:.6f}, S_t={S_rows_base[2]:.6f}")
print(f"  Signatures: {sigs_base}")

assert abs(S_base - 0.031387) < 1e-4, f"Baseline drift: {S_base}"
assert sigs_base == EXPECTED_SIGS, f"Baseline signature drift: {sigs_base}"

# ═══════════════════════════════════════════════════════════════
# Q1: PER-TERM CONTRIBUTION PROFILE
# ═══════════════════════════════════════════════════════════════

print(f"\n{'='*70}")
print("Q1: Per-term contribution profile (9 transitions × 6 terms)")
print(f"{'='*70}\n")

# Header
print(f"  {'trans':<6}", end="")
for tn in TERM_NAMES:
    print(f" {tn:>8}", end="")
print(f" {'Z_total':>8}")
print(f"  {'─'*6}", end="")
for _ in TERM_NAMES:
    print(f" {'─'*8}", end="")
print(f" {'─'*8}")

for idx, ti in enumerate(trans_info):
    label = f"{ti['src']}→{ti['tgt']}"
    z_terms = Z_all[idx]
    z_total = z_terms.sum()
    print(f"  {label:<6}", end="")
    for k in range(6):
        print(f" {z_terms[k]:>8.4f}", end="")
    print(f" {z_total:>8.4f}")

# Signed shares
print(f"\n  Signed share Z_k / Z_total:")
print(f"  {'trans':<6}", end="")
for tn in TERM_NAMES:
    print(f" {tn:>8}", end="")
print()

for idx, ti in enumerate(trans_info):
    label = f"{ti['src']}→{ti['tgt']}"
    z_terms = Z_all[idx]
    z_total = z_terms.sum()
    print(f"  {label:<6}", end="")
    for k in range(6):
        share = z_terms[k] / z_total * 100 if abs(z_total) > EPS else 0
        print(f" {share:>7.1f}%", end="")
    print()

# Absolute shares
print(f"\n  Absolute share |Z_k| / Σ|Z_j|:")
print(f"  {'trans':<6}", end="")
for tn in TERM_NAMES:
    print(f" {tn:>8}", end="")
print()

for idx, ti in enumerate(trans_info):
    label = f"{ti['src']}→{ti['tgt']}"
    z_terms = Z_all[idx]
    abs_total = sum(abs(z) for z in z_terms)
    print(f"  {label:<6}", end="")
    for k in range(6):
        share = abs(z_terms[k]) / abs_total * 100 if abs_total > EPS else 0
        print(f" {share:>7.1f}%", end="")
    print()

# Structurally zero report
print(f"\n  Structurally zero terms:")
for idx, ti in enumerate(trans_info):
    label = f"{ti['src']}→{ti['tgt']}"
    zeros = [TERM_NAMES[k] for k in range(6) if abs(Z_all[idx, k]) < EPS]
    if zeros:
        print(f"    {label}: {', '.join(zeros)}")

print(f"\n  Q1 COMPLETE")

# ═══════════════════════════════════════════════════════════════
# Q2: SINGLE-TERM ABLATION
# ═══════════════════════════════════════════════════════════════

print(f"\n{'='*70}")
print("Q2: Single-term ablation (remove one at a time)")
print(f"{'='*70}\n")

single_results = {}

print(f"  {'removed':<8} {'S':>10} {'ΔS':>10} {'%degr':>8} {'sigs':>20} {'S_u':>8} {'S_c':>8} {'S_t':>8}")
print(f"  {'─'*8} {'─'*10} {'─'*10} {'─'*8} {'─'*20} {'─'*8} {'─'*8} {'─'*8}")

for k in range(6):
    mask = [True] * 6
    mask[k] = False
    S, S_rows, P, sigs = evaluate_term_mask(mask)
    dS = S - S_base
    pct = dS / S_base * 100 if S_base > EPS else 0
    sig_str = ' '.join(sigs)
    sig_ok = (sigs == EXPECTED_SIGS)

    single_results[k] = {'S': S, 'dS': dS, 'S_rows': S_rows, 'sigs': sigs, 'sig_ok': sig_ok}

    marker = '✓' if sig_ok else '✗'
    print(f"  {TERM_NAMES[k]:<8} {S:>10.6f} {dS:>+10.6f} {pct:>+7.1f}% "
          f"{sig_str:>20} {marker} {S_rows[0]:>8.4f} {S_rows[1]:>8.4f} {S_rows[2]:>8.4f}")

# Per-row deltas
print(f"\n  Per-row ΔS (vs baseline S_u={S_rows_base[0]:.4f}, S_c={S_rows_base[1]:.4f}, S_t={S_rows_base[2]:.4f}):")
print(f"  {'removed':<8} {'ΔS_u':>10} {'ΔS_c':>10} {'ΔS_t':>10}")
for k in range(6):
    sr = single_results[k]['S_rows']
    print(f"  {TERM_NAMES[k]:<8} {sr[0]-S_rows_base[0]:>+10.4f} {sr[1]-S_rows_base[1]:>+10.4f} {sr[2]-S_rows_base[2]:>+10.4f}")

# Ranking
ranked = sorted(range(6), key=lambda k: single_results[k]['dS'], reverse=True)
print(f"\n  Load-bearing ranking (most to least impact):")
for rank, k in enumerate(ranked):
    print(f"    {rank+1}. {TERM_NAMES[k]:<8} ΔS={single_results[k]['dS']:>+.6f}")

print(f"\n  Q2 COMPLETE")

# ═══════════════════════════════════════════════════════════════
# Q3: PAIRWISE ABLATION WITH INTERACTION
# ═══════════════════════════════════════════════════════════════

print(f"\n{'='*70}")
print("Q3: Pairwise ablation with interaction (15 pairs)")
print(f"{'='*70}\n")

print(f"  {'pair':<16} {'S':>10} {'ΔS':>10} {'interact':>10} {'sigs_ok':>8}")
print(f"  {'─'*16} {'─'*10} {'─'*10} {'─'*10} {'─'*8}")

pair_results = {}
for a in range(6):
    for b in range(a+1, 6):
        mask = [True] * 6
        mask[a] = False
        mask[b] = False
        S, S_rows, P, sigs = evaluate_term_mask(mask)
        dS = S - S_base
        # Interaction = actual ΔS - (ΔS_a + ΔS_b), all after softmax
        interaction = dS - single_results[a]['dS'] - single_results[b]['dS']
        sig_ok = (sigs == EXPECTED_SIGS)

        pair_label = f"{TERM_NAMES[a]}+{TERM_NAMES[b]}"
        pair_results[(a,b)] = {'S': S, 'dS': dS, 'interaction': interaction, 'sig_ok': sig_ok}

        marker = '✓' if sig_ok else '✗'
        print(f"  {pair_label:<16} {S:>10.6f} {dS:>+10.6f} {interaction:>+10.6f} {marker:>8}")

# Strongest interactions
sorted_pairs = sorted(pair_results.items(), key=lambda x: abs(x[1]['interaction']), reverse=True)
print(f"\n  Strongest interactions:")
for (a,b), res in sorted_pairs[:5]:
    label = f"{TERM_NAMES[a]}+{TERM_NAMES[b]}"
    sign = "coupled" if res['interaction'] > 0 else "compensating"
    print(f"    {label:<16} interaction={res['interaction']:>+.6f} ({sign})")

print(f"\n  Q3 COMPLETE")

# ═══════════════════════════════════════════════════════════════
# Q4: FULL SUBSET ABLATION (63 subsets + 6 single-term-only)
# ═══════════════════════════════════════════════════════════════

print(f"\n{'='*70}")
print("Q4: Full subset ablation (63 non-empty subsets)")
print(f"{'='*70}\n")

subset_results = []

for bits in range(1, 64):  # 1 to 63: all non-empty subsets
    mask = [(bits >> k) & 1 == 1 for k in range(6)]
    n_terms = sum(mask)
    active = [TERM_NAMES[k] for k in range(6) if mask[k]]
    S, S_rows, P, sigs = evaluate_term_mask(mask)
    dS = S - S_base
    sig_ok = (sigs == EXPECTED_SIGS)

    subset_results.append({
        'bits': bits, 'mask': mask, 'n_terms': n_terms,
        'active': active, 'S': S, 'dS': dS, 'sig_ok': sig_ok,
        'S_rows': S_rows, 'sigs': sigs,
    })

# Sort by S
subset_results.sort(key=lambda x: x['S'])

# Top 20 by S
print(f"  Top 20 subsets by S (lowest = best):\n")
print(f"  {'#':>3} {'n':>2} {'S':>10} {'ΔS':>10} {'sigs':>5} {'terms'}")
print(f"  {'─'*3} {'─'*2} {'─'*10} {'─'*10} {'─'*5} {'─'*40}")
for i, sr in enumerate(subset_results[:20]):
    marker = '✓' if sr['sig_ok'] else '✗'
    terms_str = '+'.join(sr['active'])
    print(f"  {i+1:>3} {sr['n_terms']:>2} {sr['S']:>10.6f} {sr['dS']:>+10.6f} {marker:>5} {terms_str}")

# Minimum term count at each tolerance
print(f"\n  Minimum effective term count (with 3/3 signature preservation):\n")
tolerances = [('S ≤ 0.05', 0.05), ('S ≤ 0.10', 0.10), ('S ≤ 2×S_base', 2 * S_base)]
for tol_name, tol_val in tolerances:
    candidates = [sr for sr in subset_results if sr['S'] <= tol_val and sr['sig_ok']]
    if candidates:
        best = min(candidates, key=lambda x: x['n_terms'])
        terms_str = '+'.join(best['active'])
        print(f"    {tol_name:<16}: {best['n_terms']} terms (S={best['S']:.6f}, {terms_str})")
    else:
        print(f"    {tol_name:<16}: no subset qualifies")

# Single-term-only models
print(f"\n  Single-term-only models (Z_total = Z_k only):\n")
print(f"  {'term':<8} {'S':>10} {'sigs':>20} {'sigs_ok':>8}")
print(f"  {'─'*8} {'─'*10} {'─'*20} {'─'*8}")
for k in range(6):
    mask = [False] * 6
    mask[k] = True
    S, S_rows, P, sigs = evaluate_term_mask(mask)
    sig_ok = (sigs == EXPECTED_SIGS)
    sig_str = ' '.join(sigs)
    marker = '✓' if sig_ok else '✗'
    print(f"  {TERM_NAMES[k]:<8} {S:>10.6f} {sig_str:>20} {marker:>8}")

print(f"\n  Q4 COMPLETE")

# ═══════════════════════════════════════════════════════════════
# SUMMARY
# ═══════════════════════════════════════════════════════════════

print(f"\n{'='*70}")
print("SUMMARY")
print(f"{'='*70}\n")

print(f"  Baseline S = {S_base:.6f} (all 6 terms)")

# Most load-bearing
print(f"\n  Load-bearing ranking:")
for rank, k in enumerate(ranked):
    sr = single_results[k]
    print(f"    {rank+1}. {TERM_NAMES[k]:<8} ΔS={sr['dS']:>+.6f} sigs={'✓' if sr['sig_ok'] else '✗'}")

# Minimum term count summary
print(f"\n  Minimum effective term count:")
for tol_name, tol_val in tolerances:
    candidates = [sr for sr in subset_results if sr['S'] <= tol_val and sr['sig_ok']]
    if candidates:
        best = min(candidates, key=lambda x: x['n_terms'])
        print(f"    {tol_name:<16}: {best['n_terms']} terms")
    else:
        print(f"    {tol_name:<16}: no subset qualifies")

# Terms with broken signatures when removed
sig_breakers = [TERM_NAMES[k] for k in range(6) if not single_results[k]['sig_ok']]
sig_safe = [TERM_NAMES[k] for k in range(6) if single_results[k]['sig_ok']]
print(f"\n  Signature-breaking on removal: {sig_breakers if sig_breakers else 'none'}")
print(f"  Signature-safe on removal:     {sig_safe if sig_safe else 'none'}")

print(f"\n  Q1 — Contribution profile:     COMPLETE")
print(f"  Q2 — Single ablation:          COMPLETE")
print(f"  Q3 — Pairwise interaction:     COMPLETE")
print(f"  Q4 — Full subset ablation:     COMPLETE")

print(f"\n{'='*70}")
print(f"COMPLETE — {time.time()-t0:.1f}s")
print(f"{'='*70}")
