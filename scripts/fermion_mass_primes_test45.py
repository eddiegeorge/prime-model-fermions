#!/usr/bin/env python3
"""
Test 45: Phase 6 — Consecutive-prime band overlap diagnostic.
  Q1: Band overlap profile (9 transitions)
  Q2: Band-overlap Z-term candidates (5 candidates × 6th term)
  Q3: Complement vs replacement of Z_ch (5 candidates × 2 modes)
Usage:
  python fermion_mass_primes_test45.py
"""
import warnings; warnings.filterwarnings("ignore")
import numpy as np
from itertools import combinations
from sympy import primerange, primepi
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
# CONSECUTIVE PRIME BANDS
# ═══════════════════════════════════════════════════════════════

BANDS = {
    5:   [2, 3],
    17:  [2, 3, 5, 7],
    53:  [5, 7, 11, 13, 17],
    59:  [17, 19, 23],
    97:  [29, 31, 37],
    127: [3, 5, 7, 11, 13, 17, 19, 23, 29],
    211: [67, 71, 73],
    271: [7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43],
}

# Verify bands sum to primes
for p, band in BANDS.items():
    assert sum(band) == p, f"Band for {p} sums to {sum(band)}"

# Build prime index lookup for gap calculation
all_primes = sorted(set(q for band in BANDS.values() for q in band))
max_prime = max(all_primes)
prime_list = list(primerange(2, max_prime + 10))
prime_to_idx = {p: i for i, p in enumerate(prime_list)}

# ═══════════════════════════════════════════════════════════════
# TRANSITION STRUCTURE AND Z-TERM COMPUTATION
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

def compute_z_terms(ti):
    """Return [Z_bal, Z_ch, Z_R, Z_deg, Z_tp] for one transition (5-term model)."""
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

    return np.array([Z_bal, Z_ch, Z_R, Z_deg, Z_tp])

# Precompute baseline Z-terms
Z_base = np.zeros((9, 5))
for idx, ti in enumerate(trans_info):
    Z_base[idx] = compute_z_terms(ti)

TERM_NAMES_5 = ['Z_bal', 'Z_ch', 'Z_R', 'Z_deg', 'Z_tp']

# ═══════════════════════════════════════════════════════════════
# SCORING
# ═══════════════════════════════════════════════════════════════

def score_from_z_matrix(Z_mat):
    P = np.zeros((3,3))
    for r in range(3):
        zr = Z_mat[r] - Z_mat[r].min()
        ex = np.exp(-zr)
        P[r] = ex / ex.sum()

    S_rows = np.zeros(3)
    for r in range(3):
        S_rows[r] = sum((np.log(P[r,c]+EPS) - np.log(P_obs[r,c]+EPS))**2 for c in range(3))

    sigs = []
    for r in range(3):
        order = np.argsort(-P[r])
        sigs.append('>'.join(DOWN_QUARKS[k] for k in order))

    return S_rows.sum(), S_rows, P, sigs

# Baseline
Z_mat_base = np.zeros((3,3))
for idx in range(9):
    r, c = idx // 3, idx % 3
    Z_mat_base[r, c] = Z_base[idx].sum()

S_baseline, S_rows_base, _, sigs_base = score_from_z_matrix(Z_mat_base)

# ═══════════════════════════════════════════════════════════════
# BAND OVERLAP COMPUTATION
# ═══════════════════════════════════════════════════════════════

def compute_band_overlap(ti):
    """Compute band overlap metrics for one transition."""
    shed_primes = [CANONICAL[s] for s in ti['shed']]
    gain_primes = [CANONICAL[s] for s in ti['gain']]

    # Set overlap
    B_shed = set()
    for sp in shed_primes:
        B_shed.update(BANDS[sp])
    B_gain = set()
    for gp in gain_primes:
        B_gain.update(BANDS[gp])

    # Multiset (with multiplicity)
    shed_multi = []
    for sp in shed_primes:
        shed_multi.extend(BANDS[sp])
    gain_multi = []
    for gp in gain_primes:
        gain_multi.extend(BANDS[gp])

    if not B_gain:  # pure shedding
        return {
            'B_shed': B_shed, 'B_gain': B_gain,
            'shared_set': set(), 'shared_count': 0,
            'jaccard': 0.0, 'shared_sum': 0,
            'gap': None, 'contiguous': False,
            'shed_multi': shed_multi, 'gain_multi': gain_multi,
            'multi_shared_count': 0,
            'pure_shedding': True,
        }

    shared = B_shed & B_gain
    union = B_shed | B_gain
    jaccard = len(shared) / len(union) if union else 0.0
    shared_sum = sum(shared)

    # Prime-index gap
    min_gap = float('inf')
    for qs in B_shed:
        for qg in B_gain:
            gap = abs(prime_to_idx[qs] - prime_to_idx[qg])
            if gap < min_gap:
                min_gap = gap
    if min_gap == 0:
        min_gap = 0  # direct overlap

    # Contiguity
    all_q = sorted(union)
    all_idx = [prime_to_idx[q] for q in all_q]
    contiguous = (max(all_idx) - min(all_idx) + 1) == len(all_idx)

    # Multiset shared count
    from collections import Counter
    shed_counts = Counter(shed_multi)
    gain_counts = Counter(gain_multi)
    multi_shared = sum(min(shed_counts[q], gain_counts[q]) for q in shed_counts if q in gain_counts)

    return {
        'B_shed': B_shed, 'B_gain': B_gain,
        'shared_set': shared, 'shared_count': len(shared),
        'jaccard': jaccard, 'shared_sum': shared_sum,
        'gap': min_gap, 'contiguous': contiguous,
        'shed_multi': shed_multi, 'gain_multi': gain_multi,
        'multi_shared_count': multi_shared,
        'pure_shedding': False,
    }

# Precompute all overlaps
overlaps = []
for ti in trans_info:
    overlaps.append(compute_band_overlap(ti))

# ═══════════════════════════════════════════════════════════════
# BAND Z-TERM CANDIDATES
# ═══════════════════════════════════════════════════════════════

BAND_NAMES = ['count', 'jaccard', 'sum/P', 'gap', 'ln_overlap']

def compute_z_band(idx, candidate):
    """Compute Z_band for one transition and one candidate."""
    ov = overlaps[idx]
    if ov['pure_shedding']:
        return 0.0

    if candidate == 0:  # shared count / max band size
        denom = max(len(ov['B_shed']), len(ov['B_gain']))
        return -ALPHA_S * ov['shared_count'] / denom if denom > 0 else 0.0

    elif candidate == 1:  # Jaccard
        return -ALPHA_S * ov['jaccard']

    elif candidate == 2:  # shared sum / P
        return -ALPHA_S * ov['shared_sum'] / P_TOTAL

    elif candidate == 3:  # gap (impedance, not conductance)
        if ov['gap'] is None or ov['gap'] == 0:
            return 0.0
        return ALPHA_S * ov['gap'] / P_MAX

    elif candidate == 4:  # ln(1 + shared count)
        return -ALPHA_S * np.log(1 + ov['shared_count'])

    return 0.0

# ═══════════════════════════════════════════════════════════════

print("=" * 70)
print("TEST 45: Phase 6 — Consecutive-prime band overlap diagnostic")
print("=" * 70)

print(f"\n  Baseline: 5-term model (Z_bal + Z_ch + Z_R + Z_deg + Z_tp)")
print(f"  S_baseline = {S_baseline:.6f}")
print(f"  Per-row: S_u={S_rows_base[0]:.6f}, S_c={S_rows_base[1]:.6f}, S_t={S_rows_base[2]:.6f}")
print(f"  Signatures: {sigs_base}")

assert abs(S_baseline - 0.026666) < 1e-4, f"Baseline drift: {S_baseline}"
assert sigs_base == EXPECTED_SIGS, f"Baseline signature drift: {sigs_base}"

# ═══════════════════════════════════════════════════════════════
# Q1: BAND OVERLAP PROFILE
# ═══════════════════════════════════════════════════════════════

print(f"\n{'='*70}")
print("Q1: Band overlap profile (9 transitions)")
print(f"{'='*70}\n")

print(f"  {'trans':<6} {'|B_sh|':>6} {'|B_gn|':>6} {'|shrd|':>6} "
      f"{'Jacc':>6} {'Σshrd':>6} {'gap':>4} {'cont':>5} {'multi':>6}")
print(f"  {'─'*6} {'─'*6} {'─'*6} {'─'*6} {'─'*6} {'─'*6} {'─'*4} {'─'*5} {'─'*6}")

for idx, ti in enumerate(trans_info):
    ov = overlaps[idx]
    label = f"{ti['src']}→{ti['tgt']}"
    if ov['pure_shedding']:
        print(f"  {label:<6} {len(ov['B_shed']):>6} {'—':>6} {'—':>6} "
              f"{'—':>6} {'—':>6} {'—':>4} {'—':>5} {'—':>6}  (pure shedding)")
    else:
        gap_str = str(ov['gap']) if ov['gap'] is not None else '—'
        cont_str = '✓' if ov['contiguous'] else '✗'
        print(f"  {label:<6} {len(ov['B_shed']):>6} {len(ov['B_gain']):>6} "
              f"{ov['shared_count']:>6} {ov['jaccard']:>6.3f} {ov['shared_sum']:>6} "
              f"{gap_str:>4} {cont_str:>5} {ov['multi_shared_count']:>6}")

# Detail for non-pure-shedding transitions
print(f"\n  Detailed shared constituents:")
for idx, ti in enumerate(trans_info):
    ov = overlaps[idx]
    if ov['pure_shedding']:
        continue
    label = f"{ti['src']}→{ti['tgt']}"
    shed_primes = [CANONICAL[s] for s in ti['shed']]
    gain_primes = [CANONICAL[s] for s in ti['gain']]
    shared_str = ','.join(str(q) for q in sorted(ov['shared_set']))
    print(f"    {label}: shed bands {shed_primes} ∩ gain bands {gain_primes} = {{{shared_str}}}")

print(f"\n  Q1 COMPLETE")

# ═══════════════════════════════════════════════════════════════
# Q2: BAND-OVERLAP Z-TERM CANDIDATES
# ═══════════════════════════════════════════════════════════════

print(f"\n{'='*70}")
print("Q2: Band-overlap Z-term candidates (added as 6th term)")
print(f"{'='*70}\n")

q2_results = {}

for cand in range(5):
    # Compute Z_band for all 9 transitions
    Z_band = np.zeros(9)
    for idx in range(9):
        Z_band[idx] = compute_z_band(idx, cand)

    # Add to baseline
    Z_mat = np.zeros((3,3))
    for idx in range(9):
        r, c = idx // 3, idx % 3
        Z_mat[r, c] = Z_base[idx].sum() + Z_band[idx]

    S, S_rows, P, sigs = score_from_z_matrix(Z_mat)
    dS = S - S_baseline
    sig_ok = (sigs == EXPECTED_SIGS)

    q2_results[cand] = {'S': S, 'dS': dS, 'S_rows': S_rows, 'sigs': sigs, 'sig_ok': sig_ok}

    improved = S < S_baseline
    marker = '✓' if sig_ok else '✗'
    arrow = '▼' if improved else '▲'

    print(f"  Candidate {cand+1} ({BAND_NAMES[cand]}):")
    print(f"    S = {S:.6f} (ΔS = {dS:+.6f}) {arrow} sigs={marker}")
    print(f"    Per-row: S_u={S_rows[0]:.6f} (Δ{S_rows[0]-S_rows_base[0]:+.6f}), "
          f"S_c={S_rows[1]:.6f} (Δ{S_rows[1]-S_rows_base[1]:+.6f}), "
          f"S_t={S_rows[2]:.6f} (Δ{S_rows[2]-S_rows_base[2]:+.6f})")

    # Diagnostic: is improvement physically interesting?
    if improved and sig_ok:
        rows_improved = sum(1 for r in range(3) if S_rows[r] < S_rows_base[r])
        rows_damaged = sum(1 for r in range(3) if S_rows[r] > S_rows_base[r] + abs(dS))
        print(f"    → {rows_improved}/3 rows improved, {rows_damaged}/3 rows damaged beyond net gain")
        if rows_damaged == 0:
            print(f"    → PHYSICALLY INTERESTING: genuine row-level improvement")
        else:
            print(f"    → CAUTION: improvement may be softmax rescaling artefact")
    print()

print(f"  Q2 COMPLETE")

# ═══════════════════════════════════════════════════════════════
# Q3: COMPLEMENT VS REPLACEMENT OF Z_ch
# ═══════════════════════════════════════════════════════════════

print(f"\n{'='*70}")
print("Q3: Complement vs replacement of Z_ch (all 5 candidates)")
print(f"{'='*70}\n")

# Also compute Z_ch-removed baseline for comparison
Z_mat_no_ch = np.zeros((3,3))
for idx in range(9):
    r, c = idx // 3, idx % 3
    Z_mat_no_ch[r, c] = sum(Z_base[idx, k] for k in range(5) if k != 1)  # skip Z_ch
S_no_ch, S_rows_no_ch, _, sigs_no_ch = score_from_z_matrix(Z_mat_no_ch)
print(f"  Reference: Z_ch removed only → S = {S_no_ch:.6f}, sigs = {sigs_no_ch}")
print()

print(f"  {'cand':<12} {'complement S':>13} {'replace S':>13} {'base':>13}")
print(f"  {'─'*12} {'─'*13} {'─'*13} {'─'*13}")

best_overall = {'S': S_baseline, 'label': 'baseline', 'mode': '5-term'}

for cand in range(5):
    Z_band = np.zeros(9)
    for idx in range(9):
        Z_band[idx] = compute_z_band(idx, cand)

    # Complement: 5-term + Z_band
    Z_mat_comp = np.zeros((3,3))
    for idx in range(9):
        r, c = idx // 3, idx % 3
        Z_mat_comp[r, c] = Z_base[idx].sum() + Z_band[idx]
    S_comp, S_rows_comp, _, sigs_comp = score_from_z_matrix(Z_mat_comp)

    # Replacement: remove Z_ch, add Z_band
    Z_mat_repl = np.zeros((3,3))
    for idx in range(9):
        r, c = idx // 3, idx % 3
        Z_mat_repl[r, c] = sum(Z_base[idx, k] for k in range(5) if k != 1) + Z_band[idx]
    S_repl, S_rows_repl, _, sigs_repl = score_from_z_matrix(Z_mat_repl)

    comp_marker = '✓' if sigs_comp == EXPECTED_SIGS else '✗'
    repl_marker = '✓' if sigs_repl == EXPECTED_SIGS else '✗'

    print(f"  {BAND_NAMES[cand]:<12} {S_comp:>10.6f} {comp_marker:>2} {S_repl:>10.6f} {repl_marker:>2} {S_baseline:>10.6f}")

    # Track best
    if S_comp < best_overall['S'] and sigs_comp == EXPECTED_SIGS:
        best_overall = {'S': S_comp, 'label': BAND_NAMES[cand], 'mode': 'complement'}
    if S_repl < best_overall['S'] and sigs_repl == EXPECTED_SIGS:
        best_overall = {'S': S_repl, 'label': BAND_NAMES[cand], 'mode': 'replacement'}

print(f"\n  Best overall: {best_overall['label']} ({best_overall['mode']}), "
      f"S = {best_overall['S']:.6f}")

# If best improves, run full subset ablation
if best_overall['S'] < S_baseline:
    print(f"\n  Improvement found. Running full subset ablation with best band term...")
    best_cand = BAND_NAMES.index(best_overall['label'])

    Z_band_best = np.zeros(9)
    for idx in range(9):
        Z_band_best[idx] = compute_z_band(idx, best_cand)

    # 6 terms: Z_bal, Z_ch, Z_R, Z_deg, Z_tp, Z_band
    TERM_NAMES_6 = TERM_NAMES_5 + ['Z_band']

    print(f"\n  Top 10 subsets of {{Z_bal, Z_ch, Z_R, Z_deg, Z_tp, Z_band}}:\n")
    print(f"  {'#':>3} {'n':>2} {'S':>10} {'ΔS':>10} {'sigs':>5} {'terms'}")
    print(f"  {'─'*3} {'─'*2} {'─'*10} {'─'*10} {'─'*5} {'─'*40}")

    subset_results = []
    for bits in range(1, 64):  # 2^6 - 1 = 63
        mask = [(bits >> k) & 1 == 1 for k in range(6)]
        n_terms = sum(mask)
        active = [TERM_NAMES_6[k] for k in range(6) if mask[k]]

        Z_mat = np.zeros((3,3))
        for idx in range(9):
            r, c = idx // 3, idx % 3
            z_val = 0.0
            for k in range(5):
                if mask[k]:
                    z_val += Z_base[idx, k]
            if mask[5]:
                z_val += Z_band_best[idx]
            Z_mat[r, c] = z_val

        S, S_rows, _, sigs = score_from_z_matrix(Z_mat)
        sig_ok = (sigs == EXPECTED_SIGS)
        subset_results.append({
            'bits': bits, 'n_terms': n_terms, 'active': active,
            'S': S, 'dS': S - S_baseline, 'sig_ok': sig_ok,
        })

    subset_results.sort(key=lambda x: x['S'])
    print(f"  Evaluated {len(subset_results)} non-empty subsets.\n")
    for i, sr in enumerate(subset_results[:10]):
        marker = '✓' if sr['sig_ok'] else '✗'
        terms_str = '+'.join(sr['active'])
        print(f"  {i+1:>3} {sr['n_terms']:>2} {sr['S']:>10.6f} {sr['dS']:>+10.6f} "
              f"{marker:>5} {terms_str}")

print(f"\n  Q3 COMPLETE")

# ═══════════════════════════════════════════════════════════════
# SUMMARY
# ═══════════════════════════════════════════════════════════════

print(f"\n{'='*70}")
print("SUMMARY")
print(f"{'='*70}\n")

print(f"  Baseline S = {S_baseline:.6f} (5-term model)")

# Q2 ranking
print(f"\n  Q2 — Band-overlap candidates as 6th term:")
ranked = sorted(q2_results.items(), key=lambda x: x[1]['S'])
for cand, res in ranked:
    improved = '▼' if res['S'] < S_baseline else '▲'
    marker = '✓' if res['sig_ok'] else '✗'
    print(f"    {BAND_NAMES[cand]:<12} S={res['S']:.6f} (ΔS={res['dS']:+.6f}) {improved} {marker}")

print(f"\n  Q3 — Best model: {best_overall['label']} ({best_overall['mode']}), "
      f"S = {best_overall['S']:.6f}")

if best_overall['S'] < S_baseline:
    pct = (S_baseline - best_overall['S']) / S_baseline * 100
    print(f"    Improvement: {pct:.1f}% over 5-term baseline")
    print(f"    Signal: STRONG — band structure is physically active")
elif best_overall['S'] < S_baseline * 1.01:
    print(f"    Signal: MODERATE — marginal effect")
else:
    print(f"    Signal: WEAK or ABSENT — band structure is stability condition only")

print(f"\n  Q1 — Band overlap profile:      COMPLETE")
print(f"  Q2 — Candidate 6th terms:       COMPLETE")
print(f"  Q3 — Complement vs replacement: COMPLETE")

print(f"\n{'='*70}")
print(f"COMPLETE — {time.time()-t0:.1f}s")
print(f"{'='*70}")
