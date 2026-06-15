#!/usr/bin/env python3
"""
Test 41: Phase 2 — Co_sum phase separation diagnostic.
  Q1: Co_sum distinctness at shared modes
  Q2: Per-fermion Pauli legality
  Q3: Energy identity E_p = Σk² = p
  Q4: Impedance identity (reproduces test40)
  Q5: Null pool selectivity (exhaustive enumeration)
Usage:
  python fermion_mass_primes_test41.py
"""
import warnings; warnings.filterwarnings("ignore")
import numpy as np
from itertools import combinations, product as iter_product
from sympy import isprime
from math import comb
from multiprocessing import Pool, cpu_count
import time, sys

# Force unbuffered output
print = lambda *args, **kwargs: __builtins__['print' if isinstance(__builtins__, dict) else 'print'](*args, **{**kwargs, 'flush': True})
import builtins
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
    'u': [0,1,2,3,4,6], 'c': [0,1,2,3,4,5,6],
    't': [0,1,2,3,4,5,6,7], 'd': [1,4,5,6,7],
    's': [0,1,4,5,6,7], 'b': [0,1,2,4,5,6,7],
}
SUBSTRATE_A = {'u': 0, 'c': 2, 't': 2, 'd': 0, 's': 2, 'b': 2}

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
# BUILD ALL A279857 REPRESENTATIONS
# ═══════════════════════════════════════════════════════════════

all_reps = {}  # prime -> [tuple of modes, ...]
for r in range(1, 10):
    for subset in combinations(range(1, 10), r):
        energy = sum(k*k for k in subset)
        if isprime(energy):
            all_reps.setdefault(energy, []).append(subset)

A279857_pool = sorted(all_reps.keys())

print("=" * 70)
print("TEST 41: Phase 2 — Co_sum phase separation diagnostic")
print("=" * 70)

print(f"\nA279857 pool: {len(A279857_pool)} primes")
print(f"Codebook: {CANONICAL}")
print(f"Multi-rep primes in pool:")
for p in A279857_pool:
    if len(all_reps[p]) > 1:
        rep_strs = ["{" + ",".join(str(k) for k in r) + "}" for r in all_reps[p]]
        print(f"  p={p:>3}: {len(all_reps[p])} reps: {', '.join(rep_strs)}")

# ═══════════════════════════════════════════════════════════════
# Q1: CO_SUM DISTINCTNESS AT SHARED MODES
# ═══════════════════════════════════════════════════════════════

print(f"\n{'='*70}")
print("Q1: Co_sum distinctness at shared modes (codebook, first reps)")
print(f"{'='*70}\n")

# Use first rep for codebook primes
codebook_first_rep = {p: all_reps[p][0] for p in CANONICAL}

q1_pass = True
for k in range(1, 10):
    primes_at_k = []
    for p in CANONICAL:
        rep = codebook_first_rep[p]
        if k in rep:
            cs = sum(j for j in rep if j != k)
            primes_at_k.append((p, cs))

    if len(primes_at_k) <= 1:
        print(f"  k={k}: {len(primes_at_k)} prime(s) — no overlap")
        continue

    cosums = [cs for _, cs in primes_at_k]
    distinct = len(set(cosums)) == len(cosums)

    print(f"  k={k}: {len(primes_at_k)} primes sharing this mode")
    for p, cs in primes_at_k:
        print(f"    p={p:>3}, co_sum={cs}")
    if distinct:
        print(f"    → ALL DISTINCT ✓")
    else:
        print(f"    → COLLISION ✗")
        q1_pass = False

print(f"\n  Q1 {'PASSED' if q1_pass else 'FAILED'}: "
      f"{'all' if q1_pass else 'not all'} raw co_sums distinct at shared modes")

# ═══════════════════════════════════════════════════════════════
# Q2: PER-FERMION PAULI LEGALITY
# ═══════════════════════════════════════════════════════════════

print(f"\n{'='*70}")
print("Q2: Per-fermion Pauli legality (all representation combinations)")
print(f"{'='*70}\n")

def check_collision(prime_rep_assignment):
    """Check if a set of (prime, rep) assignments has any (k, co_sum) collision."""
    occupied = {}
    for p, rep in prime_rep_assignment:
        for k in rep:
            cs = sum(j for j in rep if j != k)
            key = (k, cs)
            if key in occupied and occupied[key] != p:
                return True, key, occupied[key], p
            occupied[key] = p
    return False, None, None, None

q2_pass = True
for fname, slots in PARTICLES.items():
    primes = [CANONICAL[s] for s in slots]
    reps_per_prime = [all_reps[p] for p in primes]
    n_combos = 1
    for r in reps_per_prime:
        n_combos *= len(r)

    collision_free = 0
    first_collision = None

    for combo in iter_product(*reps_per_prime):
        assignment = list(zip(primes, combo))
        has_coll, key, p1, p2 = check_collision(assignment)
        if not has_coll:
            collision_free += 1
        elif first_collision is None:
            first_collision = (combo, key, p1, p2)

    status = "✓" if collision_free > 0 else "✗"
    print(f"  {fname}: {n_combos} rep combos, {collision_free} collision-free {status}")
    if collision_free == 0:
        q2_pass = False
        k, cs = first_collision[1]
        print(f"    Example collision: (k={k}, co_sum={cs}) between "
              f"primes {first_collision[2]} and {first_collision[3]}")
    elif collision_free < n_combos:
        k, cs = first_collision[1]
        print(f"    Some combos collide: e.g. (k={k}, co_sum={cs}) between "
              f"primes {first_collision[2]} and {first_collision[3]}")

print(f"\n  Q2 {'PASSED' if q2_pass else 'FAILED'}: "
      f"{'every' if q2_pass else 'not every'} fermion has at least one "
      f"collision-free representation assignment")

# ═══════════════════════════════════════════════════════════════
# Q3: ENERGY IDENTITY
# ═══════════════════════════════════════════════════════════════

print(f"\n{'='*70}")
print("Q3: Energy identity E_p = Σk² = p (all primes, all representations)")
print(f"{'='*70}\n")

q3_pass = True
for p in CANONICAL:
    for ri, rep in enumerate(all_reps[p]):
        E_p = sum(k*k for k in rep)
        match = E_p == p
        if not match:
            q3_pass = False
        modes_str = "{" + ",".join(str(k) for k in rep) + "}"
        print(f"  p={p:>3} rep {ri}: {modes_str} → Σk²={E_p} "
              f"{'✓' if match else '✗ MISMATCH'}")

print(f"\n  Q3 {'PASSED' if q3_pass else 'FAILED'}: "
      f"E_p = p for {'all' if q3_pass else 'not all'} representations")

# ═══════════════════════════════════════════════════════════════
# Q4: IMPEDANCE IDENTITY (REPRODUCE TEST40)
# ═══════════════════════════════════════════════════════════════

print(f"\n{'='*70}")
print("Q4: Impedance identity (reproduce test40, S = 0.031387)")
print(f"{'='*70}\n")

# Transition structure
trans_info = []
for src_n, tgt_n in CKM_TRANS:
    src_sl, tgt_sl = set(PARTICLES[src_n]), set(PARTICLES[tgt_n])
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

E = {i: float(CANONICAL[i]) for i in range(8)}

def compute_z_terms(ti):
    shed_E = [E[s] for s in ti['shed']]
    gain_E = [E[s] for s in ti['gain']]
    shared_E = [E[s] for s in ti['shared']]
    sum_shed = sum(shed_E)
    sum_gain = sum(gain_E)
    sum_shared = sum(shared_E)

    Z_bal = ALPHA_S * abs(sum_shed - sum_gain)

    if ti['pure_shedding'] or ti['n_gain'] == 0:
        Z_ch = 0.0
    else:
        N = ti['n_shed'] * ti['n_gain']
        K_sum = 0.0
        for s in ti['shed']:
            for g in ti['gain']:
                dist = abs(E[s] - E[g])
                pos = max(E[s], E[g]) / P_MAX
                K_sum += np.exp(-ALPHA_S * dist * pos)
        K = K_sum / np.sqrt(N)
        Z_ch = -np.log(K)

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
        sum_rewired = sum_shared
        back_emf = sum_gain * ti['abs_da']
        Z_R = (np.pi / 2.0) * ti['abs_da'] * (
            1.0 - ALPHA_S * sum_rewired / P_TOTAL
                + ALPHA_S * back_emf / P_TOTAL)

    Z_a = (ALPHA_S / P_TOTAL) * (sum_shed * ti['a_src'] + sum_gain * ti['a_tgt'])

    if ti['abs_da'] == 0 or ti['pure_shedding']:
        Z_deg = 0.0
    else:
        ln_shared = sum(np.log(e) for e in shared_E) if shared_E else 0.0
        ln_shed = sum(np.log(e) for e in shed_E) if shed_E else 0.0
        ln_gain = sum(np.log(e) for e in gain_E) if gain_E else 0.0
        denom = ln_shed + ln_gain
        J = ln_shared / denom if denom > EPS else 0.0
        Z_deg = (ALPHA_S / P_TOTAL) * ti['abs_da'] * (sum_shed + sum_gain) * J

    if ti['pure_shedding'] and ti['abs_da'] == 0:
        Q = ti['n_src'] * (ti['n_src'] - 1)
        Z_tp = (ALPHA_S / Q) * sum_shed
    else:
        Z_tp = 0.0

    return Z_bal + Z_ch + Z_R + Z_a + Z_deg + Z_tp

Z_matrix = np.zeros((3, 3))
for ti_idx, ti in enumerate(trans_info):
    r, c = ti_idx // 3, ti_idx % 3
    Z_matrix[r, c] = compute_z_terms(ti)

P = np.zeros((3, 3))
for r in range(3):
    row_Z = Z_matrix[r]
    row_Z_shifted = row_Z - row_Z.min()
    row_exp = np.exp(-row_Z_shifted)
    P[r] = row_exp / row_exp.sum()

S = sum((np.log(P[i,j]+EPS) - np.log(P_obs[i,j]+EPS))**2
        for i in range(3) for j in range(3))

q4_pass = abs(S - 0.031387) < 1e-4
print(f"  S = {S:.6f} (target: 0.031387)")
print(f"  Q4 {'PASSED' if q4_pass else 'FAILED'}: "
      f"{'matches' if q4_pass else 'does not match'} test40")

# ═══════════════════════════════════════════════════════════════
# Q5: NULL POOL SELECTIVITY (EXHAUSTIVE ENUMERATION)
# ═══════════════════════════════════════════════════════════════

print(f"\n{'='*70}")
print("Q5: Null pool selectivity (exhaustive enumeration)")
print(f"{'='*70}\n")

pool = A279857_pool
pool_size = len(pool)
total_subsets = comb(pool_size, 8)

print(f"  Pool: {pool_size} primes")
print(f"  Total 8-prime subsets: C({pool_size},8) = {total_subsets:,}")

# Multi-rep primes in pool
multi_rep_primes = {p for p in pool if len(all_reps[p]) > 1}
print(f"  Primes with multiple representations: {len(multi_rep_primes)}")
for p in sorted(multi_rep_primes):
    print(f"    p={p}: {len(all_reps[p])} reps")

# Precompute pairwise collision status for all pool pairs
# For each pair (p,q): check all rep combinations
#   "always": every combo collides
#   "never": no combo collides
#   "sometimes": mixed
print(f"\n  Precomputing pairwise collision status for C({pool_size},2) = {comb(pool_size,2)} pairs...")

def pair_shares_kr(rep_p, rep_q):
    """Check if two specific representations share any (k, co_sum)."""
    for k in rep_p:
        if k in rep_q:
            cs_p = sum(j for j in rep_p if j != k)
            cs_q = sum(j for j in rep_q if j != k)
            if cs_p == cs_q:
                return True
    return False

pair_status = {}  # (p,q) -> "always" | "never" | "sometimes"
always_edges = set()
sometimes_pairs = set()

for i, p in enumerate(pool):
    for j, q in enumerate(pool):
        if j <= i:
            continue
        collides_any = False
        free_any = False
        for rp in all_reps[p]:
            for rq in all_reps[q]:
                if pair_shares_kr(rp, rq):
                    collides_any = True
                else:
                    free_any = True
                if collides_any and free_any:
                    break
            if collides_any and free_any:
                break

        if collides_any and not free_any:
            pair_status[(p, q)] = "always"
            always_edges.add((p, q))
        elif free_any and not collides_any:
            pair_status[(p, q)] = "never"
        else:
            pair_status[(p, q)] = "sometimes"
            sometimes_pairs.add((p, q))

n_always = sum(1 for v in pair_status.values() if v == "always")
n_never = sum(1 for v in pair_status.values() if v == "never")
n_sometimes = sum(1 for v in pair_status.values() if v == "sometimes")
print(f"  Pairwise status: {n_always} always-colliding, {n_never} never-colliding, "
      f"{n_sometimes} sometimes-colliding")

# Enumerate all C(pool,8) subsets

# --- Move viability checker and supporting data to module level for multiprocessing ---
# (Workers inherit module globals via fork on Linux)

_always_edges = set()
_sometimes_pairs = set()
_all_reps_global = all_reps  # already computed above

def _check_collision_simple(prime_rep_list):
    """Check if a list of (prime, rep) has any (k, co_sum) collision."""
    occupied = {}
    for p, rep in prime_rep_list:
        for k in rep:
            cs = sum(j for j in rep if j != k)
            key = (k, cs)
            if key in occupied and occupied[key] != p:
                return True
            occupied[key] = p
    return False

def subset_is_viable(subset):
    """Check if an 8-prime subset has a collision-free representation assignment."""
    for i_idx in range(len(subset)):
        for j_idx in range(i_idx + 1, len(subset)):
            p, q = subset[i_idx], subset[j_idx]
            key = (min(p, q), max(p, q))
            if key in _always_edges:
                return False

    subset_set = set(subset)
    sometimes_in_subset = [(p, q) for p, q in _sometimes_pairs
                           if p in subset_set and q in subset_set]

    if not sometimes_in_subset:
        return True

    involved = set()
    for p, q in sometimes_in_subset:
        involved.add(p)
        involved.add(q)

    involved_list = sorted(involved)
    involved_reps = [_all_reps_global[p] for p in involved_list]
    fixed_assignment = [(p, _all_reps_global[p][0]) for p in subset if p not in involved]

    for combo in iter_product(*involved_reps):
        assignment = fixed_assignment + list(zip(involved_list, combo))
        if not _check_collision_simple(assignment):
            return True

    return False

def _worker_batch(first_idx):
    """Process all 8-subsets starting with _pool[first_idx]."""
    remaining = _pool[first_idx + 1:]
    first_prime = _pool[first_idx]
    viable = 0
    checked = 0
    total_for_task = comb(len(remaining), 7)
    for rest in combinations(remaining, 7):
        subset = (first_prime,) + rest
        checked += 1
        if subset_is_viable(subset):
            viable += 1
        if checked % 1_000_000 == 0:
            _print(f"    [p={first_prime:>3}] {checked:,}/{total_for_task:,} "
                   f"({viable:,} viable)", flush=True)
    return first_idx, checked, viable

_pool = None  # set before forking

print(f"\n  Enumerating {total_subsets:,} subsets (parallel)...")

# Set module globals before fork
_always_edges = always_edges
_sometimes_pairs = sometimes_pairs
_pool = pool

n_workers = cpu_count()
print(f"  Using {n_workers} CPUs, {len(pool)} tasks")

t_enum = time.time()

# Each task handles subsets starting with pool[i]
tasks = list(range(len(pool)))
viable_count = 0
checked_total = 0

with Pool(n_workers) as p:
    for first_idx, checked, viable in p.imap_unordered(_worker_batch, tasks):
        viable_count += viable
        checked_total += checked
        elapsed = time.time() - t_enum
        print(f"    task {first_idx+1:>2}/{len(pool)} (first={pool[first_idx]:>3}): "
              f"{checked:>10,} subsets, {viable:>10,} viable  "
              f"[cumul: {checked_total:,}/{total_subsets:,}, {elapsed:.0f}s]")

elapsed_total = time.time() - t_enum
print(f"\n  Enumeration complete: {elapsed_total:.1f}s")

# Verify canonical codebook using the same full check
canon_viable = subset_is_viable(CANONICAL)
rival_viable_count = viable_count - int(canon_viable)

print(f"\n  Results:")
print(f"    Total 8-prime subsets:    {total_subsets:,}")
print(f"    Collision-free (viable):  {viable_count:,}")
print(f"    Fraction viable:          {viable_count/total_subsets:.4f} "
      f"({100*viable_count/total_subsets:.1f}%)")
print(f"    Eliminated:               {100*(1-viable_count/total_subsets):.1f}%")
print(f"    Canonical codebook viable: {'✓' if canon_viable else '✗'}")
print(f"    Rival viable sets:         {rival_viable_count:,}")

# ═══════════════════════════════════════════════════════════════
# SUMMARY
# ═══════════════════════════════════════════════════════════════

print(f"\n{'='*70}")
print("SUMMARY")
print(f"{'='*70}\n")

print(f"  Q1 — Co_sum distinctness:    {'PASSED ✓' if q1_pass else 'FAILED ✗'}")
print(f"  Q2 — Pauli legality:         {'PASSED ✓' if q2_pass else 'FAILED ✗'}")
print(f"  Q3 — Energy identity:        {'PASSED ✓' if q3_pass else 'FAILED ✗'}")
print(f"  Q4 — Impedance identity:     {'PASSED ✓' if q4_pass else 'FAILED ✗'}")
print(f"  Q5 — Null selectivity:       {viable_count:,}/{total_subsets:,} "
      f"({100*viable_count/total_subsets:.1f}%) collision-free, "
      f"{rival_viable_count:,} rivals")

print(f"\n{'='*70}")
print(f"COMPLETE — {time.time()-t0:.1f}s")
print(f"{'='*70}")
