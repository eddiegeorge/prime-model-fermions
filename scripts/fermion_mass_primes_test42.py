#!/usr/bin/env python3
"""
Test 42: Phase 3 — Transfer determinant mass mechanism.
  Q1: Packet transfer gain g_p = E_p = p
  Q2: Fermion mass M_f = Πp_i (mode numbers, hierarchy)
  Q3: CKM shedding ratios from packet gain algebra
  Q4: Impedance identity (reproduces test40)
Usage:
  python fermion_mass_primes_test42.py
"""
import warnings; warnings.filterwarnings("ignore")
import numpy as np
from itertools import combinations
from sympy import isprime
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

# A279857 representations
prime_reps = {}
for r in range(1, 10):
    for subset in combinations(range(1, 10), r):
        energy = sum(k*k for k in subset)
        if isprime(energy):
            if energy not in prime_reps:
                prime_reps[energy] = subset

# ═══════════════════════════════════════════════════════════════

print("=" * 70)
print("TEST 42: Phase 3 — Transfer determinant mass mechanism")
print("=" * 70)

# ═══════════════════════════════════════════════════════════════
# Q1: PACKET TRANSFER GAIN
# ═══════════════════════════════════════════════════════════════

print(f"\n{'='*70}")
print("Q1: Packet transfer gain g_p = E_p = p")
print(f"{'='*70}\n")

print("  Define effective packet-space transfer operator G = diag(g_0, ..., g_7)")
print("  with g_p identified as the packet energy E_p = Σk² = p.\n")

q1_pass = True
for i, p in enumerate(CANONICAL):
    rep = prime_reps[p]
    E_p = sum(k*k for k in rep)
    g_p = E_p  # Phase 3 identification
    modes_str = "{" + ",".join(str(k) for k in rep) + "}"
    match = (E_p == p)
    if not match:
        q1_pass = False
    print(f"  slot {i}: modes {modes_str:>15}  E_p = Σk² = {E_p:>3}  "
          f"g_p = {g_p:>3}  {'✓' if match else '✗'}")

print(f"\n  G = diag({', '.join(str(p) for p in CANONICAL)})")
print(f"  ln G = diag({', '.join(f'{np.log(p):.4f}' for p in CANONICAL)})")

print(f"\n  Q1 {'PASSED' if q1_pass else 'FAILED'}: "
      f"g_p = E_p = p for {'all' if q1_pass else 'not all'} 8 packets")

# ═══════════════════════════════════════════════════════════════
# Q2: FERMION MASS FROM TRANSFER DETERMINANT
# ═══════════════════════════════════════════════════════════════

print(f"\n{'='*70}")
print("Q2: Fermion mass M_f = det(G|_f) = Π g_p")
print(f"{'='*70}\n")

fermion_masses = {}
for fname in ['u', 'c', 't', 'd', 's', 'b']:
    slots = PARTICLES[fname]
    primes = [CANONICAL[s] for s in sorted(slots)]
    M_f = 1
    for p in primes:
        M_f *= p
    ln_M = sum(np.log(p) for p in primes)
    fermion_masses[fname] = M_f

    primes_str = ' × '.join(str(p) for p in primes)
    print(f"  {fname}: slots {sorted(slots)}")
    print(f"     M_{fname} = {primes_str} = {M_f:,}")
    print(f"     ln M_{fname} = {ln_M:.6f}")
    print()

# Verify exact mode numbers
expected_masses = {
    'u': 5440026265,
    'c': 690883335655,
    't': 187229383962505,
    'd': 11975016563,
    's': 59875082815,
    'b': 3173379389195,
}

q2_mass_pass = all(fermion_masses[f] == expected_masses[f] for f in expected_masses)
if not q2_mass_pass:
    for f in expected_masses:
        if fermion_masses[f] != expected_masses[f]:
            print(f"  MISMATCH: M_{f} = {fermion_masses[f]:,} (expected {expected_masses[f]:,})")

# Verify hierarchy
hierarchy = sorted(fermion_masses.keys(), key=lambda f: fermion_masses[f], reverse=True)
hierarchy_str = ' > '.join(f"M_{f}({fermion_masses[f]:,})" for f in hierarchy)
print(f"  Hierarchy: {hierarchy_str}")

expected_hierarchy = ['t', 'b', 'c', 's', 'd', 'u']
q2_pass = q2_mass_pass and (hierarchy == expected_hierarchy)
print(f"  Expected:  M_t > M_b > M_c > M_s > M_d > M_u")
print(f"\n  Q2 {'PASSED' if q2_pass else 'FAILED'}: "
      f"mode numbers {'correct' if q2_mass_pass else 'incorrect'}, "
      f"hierarchy {'correct' if hierarchy == expected_hierarchy else 'incorrect'}")

# ═══════════════════════════════════════════════════════════════
# Q3: CKM SHEDDING RATIOS
# ═══════════════════════════════════════════════════════════════

print(f"\n{'='*70}")
print("Q3: CKM shedding ratios from packet gain algebra")
print(f"{'='*70}\n")

print(f"  {'trans':<5} {'shed primes':>25} {'gain primes':>25} {'M_src/M_tgt':>15} {'Π_shed/Π_gain':>15} {'match':>6}")

q3_pass = True
for src_n, tgt_n in CKM_TRANS:
    src_sl = PARTICLES[src_n]
    tgt_sl = PARTICLES[tgt_n]
    shed_slots = sorted(src_sl - tgt_sl)
    gain_slots = sorted(tgt_sl - src_sl)

    shed_primes = [CANONICAL[s] for s in shed_slots]
    gain_primes = [CANONICAL[s] for s in gain_slots]

    prod_shed = 1
    for p in shed_primes:
        prod_shed *= p
    prod_gain = 1
    for p in gain_primes:
        prod_gain *= p

    M_src = fermion_masses[src_n]
    M_tgt = fermion_masses[tgt_n]

    # The ratio M_src/M_tgt should equal prod_shed/prod_gain
    ratio_mass = M_src / M_tgt
    ratio_packets = prod_shed / prod_gain
    match = abs(ratio_mass - ratio_packets) < 1e-10

    if not match:
        q3_pass = False

    shed_str = '×'.join(str(p) for p in shed_primes) if shed_primes else '(none)'
    gain_str = '×'.join(str(p) for p in gain_primes) if gain_primes else '(none)'
    label = f"{src_n}→{tgt_n}"

    print(f"  {label:<5} {shed_str:>25} {gain_str:>25} {ratio_mass:>15.4f} {ratio_packets:>15.4f} {'✓' if match else '✗':>6}")

print(f"\n  Generation ratios (single packet gains):")
gen_ratios = [
    ('c', 'u', 'g_5=127'),
    ('t', 'c', 'g_7=271'),
    ('s', 'd', 'g_0=5'),
    ('b', 's', 'g_2=53'),
]
for heavier, lighter, label in gen_ratios:
    ratio = fermion_masses[heavier] / fermion_masses[lighter]
    print(f"    M_{heavier}/M_{lighter} = {ratio:.0f} = {label}")

print(f"\n  Q3 {'PASSED' if q3_pass else 'FAILED'}: "
      f"{'all' if q3_pass else 'not all'} 9 shedding ratios consistent")

# ═══════════════════════════════════════════════════════════════
# Q4: IMPEDANCE IDENTITY (REPRODUCE TEST40)
# ═══════════════════════════════════════════════════════════════

print(f"\n{'='*70}")
print("Q4: Impedance identity (reproduce test40, S = 0.031387)")
print(f"{'='*70}\n")

# Transition structure
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

E = {i: float(CANONICAL[i]) for i in range(8)}

def compute_z_total(ti):
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
        back_emf = sum_gain * ti['abs_da']
        Z_R = (np.pi / 2.0) * ti['abs_da'] * (
            1.0 - ALPHA_S * sum_shared / P_TOTAL
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
    Z_matrix[r, c] = compute_z_total(ti)

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
# SUMMARY
# ═══════════════════════════════════════════════════════════════

print(f"\n{'='*70}")
print("SUMMARY")
print(f"{'='*70}")

print(f"\n  Two layers of the topology:")
print(f"    Mass layer:       det(G|_f) = Πp_i        (multiplicative, transfer gain)")
print(f"    Transition layer: Z_ij = Σ Z_terms         (additive, impedance costs)")
print(f"    Both use:         E_p = Σk² = p            (packet energy as common input)")

print(f"\n  Q1 — Packet transfer gain:    {'PASSED ✓' if q1_pass else 'FAILED ✗'}")
print(f"  Q2 — Fermion mass hierarchy:  {'PASSED ✓' if q2_pass else 'FAILED ✗'}")
print(f"  Q3 — Shedding ratios:         {'PASSED ✓' if q3_pass else 'FAILED ✗'}")
print(f"  Q4 — Impedance identity:      {'PASSED ✓' if q4_pass else 'FAILED ✗'}")

print(f"\n{'='*70}")
print(f"COMPLETE — {time.time()-t0:.1f}s")
print(f"{'='*70}")
