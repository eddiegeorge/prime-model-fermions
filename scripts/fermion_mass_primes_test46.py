#!/usr/bin/env python3
"""
fermion_mass_primes_test46.py
Phase 10A: Exhaustive codebook null.

Exhaustive enumeration of all C(43,8) ~ 145M candidate 8-packet codebooks
from the full A279857 prime pool, with all structural filters, negative
controls, leave-one-out analysis, and canonical impedance ranking.
"""

import itertools
import math
import time
import hashlib
import json
from collections import defaultdict
from multiprocessing import Pool, cpu_count

# ─── Constants ────────────────────────────────────────────────────────
ALPHA_S = 0.118
PI_HALF = math.pi / 2
P_TOTAL = 840
P_MAX_EXPECTED = 271
CODEBOOK = (5, 17, 53, 59, 97, 127, 211, 271)
CODEBOOK_SET = frozenset(CODEBOOK)
EXPECTED_CODEBOOK_S = 0.014188  # canonical post-Test45 control value

# Observed CKM at M_Z — UTfit 2023 values via Antusch, Hinze & Saad (2025)
# arXiv:2510.01312v2, Eq. (2.5). Pure M_Z scale, not PDG mixed convention.
import cmath
_s12, _s23, _s13, _delta = 0.2251, 0.04193, 0.00370, 1.139
_c12 = math.sqrt(1 - _s12**2)
_c23 = math.sqrt(1 - _s23**2)
_c13 = math.sqrt(1 - _s13**2)
_ed  = cmath.exp(1j * _delta)
_edn = cmath.exp(-1j * _delta)
_V = [
    [_c12*_c13,                                    _s12*_c13,                          _s13*_edn],
    [-_s12*_c23 - _c12*_s23*_s13*_ed,   _c12*_c23 - _s12*_s23*_s13*_ed,    _s23*_c13],
    [_s12*_s23 - _c12*_c23*_s13*_ed,    -_c12*_s23 - _s12*_c23*_s13*_ed,   _c23*_c13],
]
CKM_OBS = [[0.0]*3 for _ in range(3)]
for _i in range(3):
    row_abs2 = [abs(_V[_i][_j])**2 for _j in range(3)]
    row_sum = sum(row_abs2)
    for _j in range(3):
        CKM_OBS[_i][_j] = row_abs2[_j] / row_sum
EPS = 1e-12

# Row signatures: descending order of P within each row
# u: d>s>b → (0,1,2), c: s>d>b → (1,0,2), t: b>s>d → (2,1,0)
ROW_SIGS = [(0,1,2), (1,0,2), (2,1,0)]

# Fermion substrate factors a
A_UP   = [0, 2, 2]  # u, c, t
A_DOWN = [0, 2, 2]  # d, s, b

# Slot masks: which of 8 slots (0-indexed) each fermion uses
# t = all 8. Chain strips: Up: 7,5; Down: 3,2,0; Lepton: 4,1,6
SLOT_MASKS = {
    'u': frozenset({0,1,2,3,4,6}),
    'c': frozenset({0,1,2,3,4,5,6}),
    't': frozenset({0,1,2,3,4,5,6,7}),
    'd': frozenset({1,4,5,6,7}),
    's': frozenset({0,1,4,5,6,7}),
    'b': frozenset({0,1,2,4,5,6,7}),
}

UP_QUARKS  = ['u','c','t']
DOWN_QUARKS = ['d','s','b']

# Pre-compute transition structure (slot-level, independent of primes)
TRANSITIONS = []
for ui, uq in enumerate(UP_QUARKS):
    for di, dq in enumerate(DOWN_QUARKS):
        TRANSITIONS.append({
            'name': f'{uq}->{dq}', 'row': ui, 'col': di,
            'shed_slots': sorted(SLOT_MASKS[uq] - SLOT_MASKS[dq]),
            'gain_slots': sorted(SLOT_MASKS[dq] - SLOT_MASKS[uq]),
            'shared_slots': sorted(SLOT_MASKS[uq] & SLOT_MASKS[dq]),
            'delta_a': A_DOWN[di] - A_UP[ui],
            'n_src': len(SLOT_MASKS[uq]),
        })

# Codebook band window profile (sorted)
CODEBOOK_WINDOW_PROFILE = tuple(sorted([2, 4, 5, 3, 3, 9, 3, 11]))

# ─── Module-level globals set by main() before Pool ──────────────────
BANDS_GLOBAL = {}

# ─── Utility functions ────────────────────────────────────────────────

def is_prime(n):
    if n < 2: return False
    if n < 4: return True
    if n % 2 == 0 or n % 3 == 0: return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i+2) == 0: return False
        i += 6
    return True

def primes_up_to(n):
    sieve = [True] * (n + 1)
    sieve[0] = sieve[1] = False
    for i in range(2, int(n**0.5) + 1):
        if sieve[i]:
            for j in range(i*i, n+1, i): sieve[j] = False
    return [i for i in range(2, n+1) if sieve[i]]

def is_7smooth(n):
    for pr in [2, 3, 5, 7]:
        while n % pr == 0: n //= pr
    return n == 1

def p(msg):
    print(msg, flush=True)

# ─── Pool and lookup generation ───────────────────────────────────────

def compute_a279857_pool():
    squares = [k*k for k in range(1, 10)]
    prime_reps = defaultdict(list)
    for r in range(1, 10):
        for mask in itertools.combinations(range(9), r):
            s = sum(squares[i] for i in mask)
            if is_prime(s):
                modes = frozenset(i+1 for i in mask)
                prime_reps[s].append(modes)
    pool = sorted(prime_reps.keys())
    return pool, dict(prime_reps)

def compute_consecutive_prime_bands(pool_set):
    small_primes = primes_up_to(300)
    bands = defaultdict(list)
    for start in range(len(small_primes)):
        s = 0
        constituents = []
        for end in range(start, len(small_primes)):
            s += small_primes[end]
            constituents.append(small_primes[end])
            if s > 285: break
            if len(constituents) >= 2 and s in pool_set:
                bands[s].append(tuple(constituents))
    return dict(bands)

def check_cosum_legality(combo, prime_reps):
    rep_options = [prime_reps[pr] for pr in combo]
    for rep_combo in itertools.product(*rep_options):
        mode_cosums = defaultdict(list)
        for pi, (prime, modes) in enumerate(zip(combo, rep_combo)):
            for k in modes:
                cs = sum(m for m in modes if m != k)
                mode_cosums[k].append(cs)
        collision = False
        for k, cosums in mode_cosums.items():
            if len(cosums) > 1 and len(cosums) != len(set(cosums)):
                collision = True; break
        if not collision: return True
    return False

def check_cosum_first_rep_only(combo, prime_reps):
    """N4 control: co_sum using only first representation per prime."""
    reps = [prime_reps[pr][0] for pr in combo]
    mode_cosums = defaultdict(list)
    for pi, (prime, modes) in enumerate(zip(combo, reps)):
        for k in modes:
            cs = sum(m for m in modes if m != k)
            mode_cosums[k].append(cs)
    for k, cosums in mode_cosums.items():
        if len(cosums) > 1 and len(cosums) != len(set(cosums)):
            return False
    return True

def count_legal_reps(combo, prime_reps):
    rep_options = [prime_reps[pr] for pr in combo]
    count = 0
    for rep_combo in itertools.product(*rep_options):
        mode_cosums = defaultdict(list)
        for pi, (prime, modes) in enumerate(zip(combo, rep_combo)):
            for k in modes:
                cs = sum(m for m in modes if m != k)
                mode_cosums[k].append(cs)
        collision = False
        for k, cosums in mode_cosums.items():
            if len(cosums) > 1 and len(cosums) != len(set(cosums)):
                collision = True; break
        if not collision: count += 1
    return count

# ─── Impedance scoring ───────────────────────────────────────────────

def compute_impedance_score(perm, bands_dict):
    """Compute CKM impedance score for a given 8-prime permutation.
    perm: tuple of 8 primes, perm[i] = prime at slot i.
    bands_dict: {prime: [band_tuple, ...]} — full bands with lists.
    Returns (S, row_sigs_ok, per_row_S).
    """
    P = sum(perm)
    p_max = max(perm)
    z_matrix = [[0.0]*3 for _ in range(3)]

    for tr in TRANSITIONS:
        shed_primes  = [perm[s] for s in tr['shed_slots']]
        gain_primes  = [perm[s] for s in tr['gain_slots']]
        shared_primes = [perm[s] for s in tr['shared_slots']]
        da    = tr['delta_a']
        n_src = tr['n_src']
        n_shed = len(shed_primes)
        n_gain = len(gain_primes)
        sum_shed = sum(shed_primes)
        sum_gain = sum(gain_primes)

        # Z_bal
        z_bal = ALPHA_S * abs(sum_shed - sum_gain)

        # Z_ch
        z_ch = 0.0
        if n_shed > 0 and n_gain > 0:
            N = n_shed * n_gain
            K = 0.0
            for ps in shed_primes:
                for pg in gain_primes:
                    K += math.exp(-ALPHA_S * abs(ps-pg) * max(ps,pg) / p_max)
            K /= math.sqrt(N)
            z_ch = -math.log(K) if K > 1e-30 else 30.0

        # Z_R
        z_R = 0.0
        if da != 0:
            abs_da = abs(da)
            if da > 0:  # engaging
                sum_ln_sh = sum(math.log(x) for x in shared_primes) if shared_primes else 0
                denom = (sum(math.log(x) for x in shed_primes)
                       + sum(math.log(x) for x in gain_primes))
                J = sum_ln_sh / denom if denom > 0 else 0
                z_R = PI_HALF * abs_da * J
            else:  # releasing
                sum_rewired = sum(shared_primes)
                z_R = PI_HALF * abs_da * (
                    1 - ALPHA_S * sum_rewired / P
                      + ALPHA_S * sum_gain * abs_da / P)

        # Z_deg
        z_deg = 0.0
        if da != 0 and n_shed > 0 and n_gain > 0:
            sum_ln_sh = sum(math.log(x) for x in shared_primes) if shared_primes else 0
            denom = (sum(math.log(x) for x in shed_primes)
                   + sum(math.log(x) for x in gain_primes))
            J = sum_ln_sh / denom if denom > 0 else 0
            z_deg = (ALPHA_S / P) * abs(da) * (sum_shed + sum_gain) * J

        # Z_tp
        z_tp = 0.0
        if n_gain == 0 and da == 0:
            Q = n_src * (n_src - 1)
            z_tp = (ALPHA_S / Q) * sum_shed if Q > 0 else 0

        # Z_band — uses full bands_dict with list of band tuples
        z_band = 0.0
        if n_shed > 0 and n_gain > 0:
            b_shed = set()
            for ps in shed_primes:
                if ps in bands_dict:
                    for band_tuple in bands_dict[ps]:
                        b_shed.update(band_tuple)
            b_gain = set()
            for pg in gain_primes:
                if pg in bands_dict:
                    for band_tuple in bands_dict[pg]:
                        b_gain.update(band_tuple)
            overlap = len(b_shed & b_gain)
            z_band = -ALPHA_S * math.log(1 + overlap)

        z_total = z_bal + z_ch + z_R + z_deg + z_tp + z_band
        z_matrix[tr['row']][tr['col']] = z_total

    # Softmax per row — subtract row min for numerical stability (matches test45)
    p_pred = [[0.0]*3 for _ in range(3)]
    for i in range(3):
        z_min = min(z_matrix[i])
        exp_vals = [math.exp(-(z_matrix[i][j] - z_min)) for j in range(3)]
        row_sum = sum(exp_vals)
        for j in range(3):
            p_pred[i][j] = exp_vals[j] / row_sum if row_sum > 0 else 1/3

    # Score S
    S = 0.0
    for i in range(3):
        for j in range(3):
            diff = math.log(p_pred[i][j] + EPS) - math.log(CKM_OBS[i][j] + EPS)
            S += diff * diff

    # Row signatures
    row_sigs_ok = True
    for i in range(3):
        order = tuple(sorted(range(3), key=lambda j: -p_pred[i][j]))
        if order != ROW_SIGS[i]:
            row_sigs_ok = False

    per_row_S = []
    for i in range(3):
        rs = sum((math.log(p_pred[i][j]+EPS) - math.log(CKM_OBS[i][j]+EPS))**2
                 for j in range(3))
        per_row_S.append(rs)

    return S, row_sigs_ok, per_row_S

# ─── Module-level worker for Stage 3 (picklable) ─────────────────────

def score_worker(combo):
    """Full 8! assignment search for one candidate set.
    Uses module-level BANDS_GLOBAL set by main()."""
    combo_sorted = tuple(sorted(combo))
    nat_S, nat_sigs, nat_pr = compute_impedance_score(combo_sorted, BANDS_GLOBAL)
    best_S, best_perm, best_sigs, best_pr = nat_S, combo_sorted, nat_sigs, nat_pr

    for perm in itertools.permutations(combo_sorted):
        S, sigs, pr = compute_impedance_score(perm, BANDS_GLOBAL)
        if S < best_S:
            best_S, best_perm, best_sigs, best_pr = S, perm, sigs, pr

    return {
        'combo': combo_sorted,
        'nat_S': nat_S, 'nat_sigs': nat_sigs, 'nat_per_row': nat_pr,
        'best_S': best_S, 'best_perm': best_perm,
        'best_sigs': best_sigs, 'best_per_row': best_pr,
        'lift': nat_S - best_S,
    }

# ─── Stage 1 worker (picklable, module-level) ────────────────────────

# Module-level lookups set before Pool
CONSEC_SET_G   = set()
UNIQUE_BAND_G  = set()
POOL_MODES_G   = {}
POOL_G         = []
PRIME_REPS_G   = {}
BANDS_WINDOWS_G = {}  # prime -> window length

# Filter bit indices
F1, F2, F3, F5, F6 = 0, 1, 2, 3, 4
N_CHEAP = 5  # number of cheap filters (F7 retired to diagnostic)

def stage1_worker(start_idx):
    """Process all 8-combos starting with POOL_G[start_idx]."""
    pool = POOL_G
    n = len(pool)
    remaining = pool[start_idx+1:]

    # Counters: bitmask-based
    alone = [0] * N_CHEAP
    # Leave-one-out: count of combos passing all cheap except filter i
    loo = [0] * N_CHEAP
    # Pairwise with F1
    pw_f1 = [0] * N_CHEAP
    total = 0
    all_pass = 0

    survivors = []  # combos passing all cheap filters
    f1_survivors_info = []  # (combo, bitmask) for all F1 survivors

    for rest in itertools.combinations(remaining, 7):
        combo = (pool[start_idx],) + rest
        total += 1

        s = sum(combo)
        mx = max(combo)

        # Compute cheap filter results
        pf1 = (s == 840)
        pf2 = all(x in CONSEC_SET_G for x in combo)
        pf3 = all(x in UNIQUE_BAND_G for x in combo)
        # F5: 9-mode coverage (union of all reps)
        modes = set()
        for x in combo:
            modes.update(POOL_MODES_G[x])
        pf5 = (len(modes) == 9)
        pf6 = (mx == P_MAX_EXPECTED)

        bits = [pf1, pf2, pf3, pf5, pf6]
        mask = 0
        for i, b in enumerate(bits):
            if b:
                mask |= (1 << i)
                alone[i] += 1

        # Pairwise with F1
        if pf1:
            for i in range(N_CHEAP):
                if bits[i]:
                    pw_f1[i] += 1

        # Leave-one-out: passes all except filter i
        all_mask = (1 << N_CHEAP) - 1  # all bits set
        for i in range(N_CHEAP):
            without_i = all_mask ^ (1 << i)
            if (mask & without_i) == without_i:
                loo[i] += 1

        if mask == all_mask:
            survivors.append(combo)
            all_pass += 1

        # Track F1 survivors for F4 analysis
        if pf1:
            f1_survivors_info.append((combo, mask))

    return {
        'total': total, 'alone': alone, 'loo': loo,
        'pw_f1': pw_f1, 'all_pass': all_pass,
        'survivors': survivors,
        'f1_survivors_info': f1_survivors_info,
    }

# ─── Main ─────────────────────────────────────────────────────────────

def main():
    global BANDS_GLOBAL, CONSEC_SET_G, UNIQUE_BAND_G, POOL_MODES_G
    global POOL_G, PRIME_REPS_G, BANDS_WINDOWS_G

    t0 = time.time()
    ncpu = cpu_count()
    p(f"=== TEST 46: EXHAUSTIVE CODEBOOK NULL (Phase 10A) ===")
    p(f"CPUs available: {ncpu}, using all")
    p("")

    # ── Pre-computation ──
    p("--- Pre-computation ---")
    pool, prime_reps = compute_a279857_pool()
    p(f"A279857 pool: {len(pool)} primes")
    p(f"  {pool}")

    pool_modes = {}
    for pr in pool:
        all_modes = set()
        for rep in prime_reps[pr]:
            all_modes.update(rep)
        pool_modes[pr] = all_modes

    pool_set_primes = set(pool)
    bands = compute_consecutive_prime_bands(pool_set_primes)
    consec_set = set(bands.keys()) & pool_set_primes
    p(f"Consecutive-prime-sum primes in pool: {len(consec_set)}")
    p(f"  {sorted(consec_set)}")

    unique_band = {pr for pr in pool if pr in bands and len(bands[pr]) == 1}
    p(f"Primes with unique consecutive-prime representation: {len(unique_band)}")

    # Band windows
    bands_windows = {}
    for pr in pool:
        if pr in bands and len(bands[pr]) >= 1:
            bands_windows[pr] = len(bands[pr][0])

    p(f"\nRepresentations per prime:")
    for pr in pool:
        reps = prime_reps[pr]
        bnd = bands.get(pr, [])
        wl = bands_windows.get(pr, 0)
        p(f"  {pr:>3}: {len(reps)} k²-reps, {len(bnd)} band{'s' if len(bnd)!=1 else ''}"
          f", window={wl}{'  UNIQUE' if pr in unique_band else ''}")

    # Config hash
    n_candidates = math.comb(len(pool), 8)
    config = {'pool': pool, 'n': n_candidates, 'codebook': list(CODEBOOK),
              'filters': ['F1-F6','N1-N5','R1-R2','F7-diag']}
    config_hash = hashlib.md5(json.dumps(config, sort_keys=True).encode()).hexdigest()[:12]
    p(f"\nConfig hash: {config_hash}")
    p(f"Candidates: C({len(pool)},8) = {n_candidates:,}")
    p(f"Codebook in pool: {CODEBOOK_SET <= pool_set_primes}")

    # Set module-level globals for workers
    BANDS_GLOBAL    = bands
    CONSEC_SET_G    = consec_set
    UNIQUE_BAND_G   = unique_band
    POOL_MODES_G    = pool_modes
    POOL_G          = pool
    PRIME_REPS_G    = prime_reps
    BANDS_WINDOWS_G = bands_windows

    # ── Calibration gate ──
    p(f"\n--- Calibration gate ---")
    cb_S, cb_sigs, cb_pr = compute_impedance_score(CODEBOOK, bands)
    p(f"  Codebook natural S = {cb_S:.6f}")
    p(f"  Expected S         = {EXPECTED_CODEBOOK_S:.6f}")
    p(f"  Row signatures OK  = {cb_sigs}")
    p(f"  Per-row S          = {[f'{x:.6f}' for x in cb_pr]}")
    if abs(cb_S - EXPECTED_CODEBOOK_S) > 1e-6 or not cb_sigs:
        p(f"  *** CALIBRATION FAILED ***")
        p(f"  Canonical Test45 control not reproduced. Do not proceed.")
        raise RuntimeError(
            f"Calibration failed: S={cb_S:.6f} (expected {EXPECTED_CODEBOOK_S}), "
            f"sigs={cb_sigs}"
        )
    p(f"  Calibration PASSED")

    # ── Smoke test ──
    intersection_pool = sorted(consec_set)
    n_smoke = math.comb(len(intersection_pool), 8)
    p(f"\n--- Smoke test: C({len(intersection_pool)},8) = {n_smoke:,} ---")
    smoke_surv = 0
    smoke_cb = False
    for combo in itertools.combinations(intersection_pool, 8):
        if sum(combo) != 840: continue
        if max(combo) != P_MAX_EXPECTED: continue
        if not all(x in unique_band for x in combo): continue
        smoke_surv += 1
        if frozenset(combo) == CODEBOOK_SET:
            smoke_cb = True
    p(f"  Survivors (F1+F3+F6): {smoke_surv}")
    p(f"  Codebook found: {smoke_cb}")
    if not smoke_cb:
        p("  *** WARNING: codebook not in smoke test survivors ***")
    p(f"  Smoke time: {time.time()-t0:.1f}s")

    # ── Stage 1: Full enumeration ──
    p(f"\n--- Stage 1: Full C({len(pool)},8) enumeration ---")
    t1 = time.time()

    with Pool(ncpu) as executor:
        results_iter = executor.imap_unordered(stage1_worker, range(len(pool)))
        agg_alone = [0]*N_CHEAP
        agg_loo   = [0]*N_CHEAP
        agg_pw_f1 = [0]*N_CHEAP
        agg_total = 0
        all_survivors = []
        all_f1_info = []
        done = 0
        for res in results_iter:
            done += 1
            agg_total += res['total']
            for i in range(N_CHEAP):
                agg_alone[i] += res['alone'][i]
                agg_loo[i]   += res['loo'][i]
                agg_pw_f1[i] += res['pw_f1'][i]
            all_survivors.extend(res['survivors'])
            all_f1_info.extend(res['f1_survivors_info'])
            if done % 5 == 0 or done == len(pool):
                el = time.time() - t1
                p(f"  Workers: {done}/{len(pool)}, combos: {agg_total:,}, "
                  f"F1 surv: {agg_alone[F1]:,}, cascade: {len(all_survivors)}, "
                  f"elapsed: {el:.1f}s")

    FILTER_NAMES = ['F1:sum840', 'F2:consec', 'F3:uniq_band',
                    'F5:9mode', 'F6:pmax271']

    p(f"\n  Stage 1 complete: {time.time()-t1:.1f}s")
    p(f"  Total combinations: {agg_total:,}")

    p(f"\n--- Stage 1: Filter report ---")
    p(f"  {'Filter':<20} {'Alone':>12} {'Leave-1-out':>12} {'F1-pair':>12}")
    p(f"  {'-'*58}")
    for i in range(N_CHEAP):
        p(f"  {FILTER_NAMES[i]:<20} {agg_alone[i]:>12,} {agg_loo[i]:>12,} {agg_pw_f1[i]:>12,}")
    p(f"\n  Cascade survivors (F1∩F2∩F3∩F5∩F6): {len(all_survivors)}")

    cb_in_survivors = any(frozenset(s) == CODEBOOK_SET for s in all_survivors)
    p(f"  Codebook in survivors: {cb_in_survivors}")

    # ── Stage 1b: F4 co_sum legality ──
    p(f"\n--- Stage 1b: Co_sum legality for {len(all_survivors)} cascade survivors ---")
    t1b = time.time()
    cosum_survivors = []
    for i, combo in enumerate(all_survivors):
        if check_cosum_legality(combo, prime_reps):
            cosum_survivors.append(combo)
        if (i+1) % 50 == 0 or i+1 == len(all_survivors):
            p(f"  Checked {i+1}/{len(all_survivors)}, legal: {len(cosum_survivors)}")
    p(f"  F4 survivors: {len(cosum_survivors)} of {len(all_survivors)}")
    p(f"  Time: {time.time()-t1b:.1f}s")

    # F4 alone on F1-survivors (practical scope)
    p(f"\n  F4 alone on F1-survivors ({len(all_f1_info):,}):")
    f4_alone_on_f1 = 0
    for ci, (combo, mask) in enumerate(all_f1_info):
        if check_cosum_legality(combo, prime_reps):
            f4_alone_on_f1 += 1
        if (ci+1) % 10000 == 0:
            p(f"    Checked {ci+1:,}/{len(all_f1_info):,}, legal: {f4_alone_on_f1:,}")
    p(f"  F4 alone (F1-scoped): {f4_alone_on_f1:,} of {len(all_f1_info):,} "
      f"({100*f4_alone_on_f1/max(len(all_f1_info),1):.1f}%)")
    p(f"  Note: F4 full-alone over all C(43,8) not computed (too expensive).")
    p(f"  F4 leave-one-out computed cascade-scoped: sets passing F1-F3+F5-F6 but failing F4 = "
      f"{len(all_survivors) - len(cosum_survivors)}")

    # N4: first-rep-only co_sum
    p(f"\n  N4: first-rep-only co_sum for cascade survivors:")
    n4_pass = sum(1 for c in all_survivors if check_cosum_first_rep_only(c, prime_reps))
    p(f"  N4 pass: {n4_pass} of {len(all_survivors)}")

    # Final Stage 1 survivors = F1-F6 + F4 (R1/R2/F7 are diagnostics only)
    final_stage1 = cosum_survivors  # R1/R2 NOT used as filters

    # ── Role-aware diagnostics R1, R2 ──
    p(f"\n--- Role-aware diagnostics (not filters) ---")
    r1_pass = 0
    for combo in final_stage1:
        low  = sum(1 for x in combo if x % 6 == 5)
        high = sum(1 for x in combo if x % 6 == 1)
        if low == 4 and high == 4:
            r1_pass += 1
    p(f"  R1 (mod 6 split 4+4): {r1_pass} of {len(final_stage1)}")

    # R2: sorted-adjacent mod-8 pair profile under natural assignment
    # Pairs slots [0,1],[2,3],[4,5],[6,7] in ascending prime order and checks
    # whether the mod-8 pair pattern matches the codebook's. This is a set-level
    # residue diagnostic, not a role-derived sector test.
    r2_match = 0
    cb_mod8_pairs = []
    for i in range(0, 8, 2):
        cb_mod8_pairs.append((CODEBOOK[i] % 8, CODEBOOK[i+1] % 8))
    cb_mod8_pairs_sorted = tuple(sorted(cb_mod8_pairs))
    for combo in final_stage1:
        cs = tuple(sorted(combo))
        pairs = []
        for i in range(0, 8, 2):
            pairs.append((cs[i] % 8, cs[i+1] % 8))
        if tuple(sorted(pairs)) == cb_mod8_pairs_sorted:
            r2_match += 1
    p(f"  R2 (sorted-adjacent mod-8 pair profile): {r2_match} of {len(final_stage1)}")
    p(f"  Codebook mod8 pairs: {cb_mod8_pairs}")

    # F7 diagnostic (retired from cascade — empirically motivated, not derived from cavity physics)
    f7_pass = 0
    for combo in final_stage1:
        mod8 = [0]*4
        for x in combo:
            mod8[(x % 8) // 2] += 1
        if all(c == 2 for c in mod8):
            f7_pass += 1
    p(f"  F7 (mod 8 balanced, diagnostic only): {f7_pass} of {len(final_stage1)}")

    # Negative controls N1-N3, N5
    p(f"\n--- Negative controls on cascade survivors ---")
    n1_pass = sum(1 for c in all_survivors
                  if all(is_7smooth(x-1) and is_7smooth(x+1) for x in c))
    p(f"  N1 (smooth-neighbour): {n1_pass} of {len(all_survivors)}")

    # N2: band window profile matches codebook
    n2_pass = 0
    for combo in all_survivors:
        profile = tuple(sorted(bands_windows.get(x, 0) for x in combo))
        if profile == CODEBOOK_WINDOW_PROFILE:
            n2_pass += 1
    p(f"  N2 (window profile match): {n2_pass} of {len(all_survivors)}")

    # N3: A279857 mask contiguity (at least one rep has contiguous k-modes per prime)
    n3_pass = 0
    for combo in all_survivors:
        all_contig = True
        for x in combo:
            prime_contig = any(
                max(rep) - min(rep) + 1 == len(rep)
                for rep in prime_reps[x]
            )
            if not prime_contig:
                all_contig = False; break
        if all_contig:
            n3_pass += 1
    p(f"  N3 (mask contiguity): {n3_pass} of {len(all_survivors)}")
    p(f"  N5 (window as running predictor): Phase 9 negative, not used as selector")

    p(f"\n=== STAGE 1 FINAL: {len(final_stage1)} survivors (F1-F6 + F4) ===")
    p(f"  R1/R2 reported as diagnostics, not used as filters")

    # ── Stage 2: Band-transition structure ──
    p(f"\n--- Stage 2: Band-transition structure ---")
    for ci, combo in enumerate(final_stage1):
        cs = tuple(sorted(combo))
        is_cb = '***' if frozenset(cs) == CODEBOOK_SET else '   '

        # Compute band structure under natural assignment
        overlap_counts = []
        contiguous_count = 0
        for tr in TRANSITIONS:
            shed = [cs[s] for s in tr['shed_slots']]
            gain = [cs[s] for s in tr['gain_slots']]
            if not gain:
                overlap_counts.append(0)
                continue
            b_shed = set()
            for ps in shed:
                if ps in bands:
                    for bt in bands[ps]: b_shed.update(bt)
            b_gain = set()
            for pg in gain:
                if pg in bands:
                    for bt in bands[pg]: b_gain.update(bt)
            ov = len(b_shed & b_gain)
            overlap_counts.append(ov)
            union = b_shed | b_gain
            if union:
                sp = primes_up_to(max(union))
                sp_in = [x for x in sp if x >= min(union) and x <= max(union)]
                if union == set(sp_in):
                    contiguous_count += 1

        n_legal = count_legal_reps(cs, prime_reps)
        p(f"  {is_cb} {cs}  overlaps={overlap_counts}  "
          f"contiguous={contiguous_count}  legal_reps={n_legal}")

    # ── Stage 3: Full 8! impedance ranking ──
    p(f"\n--- Stage 3: Full 8! ranking for {len(final_stage1)} survivors ---")
    t3 = time.time()

    results_list = []
    with Pool(ncpu) as executor:
        for i, result in enumerate(executor.imap_unordered(score_worker, final_stage1)):
            results_list.append(result)
            if (i+1) % 5 == 0 or i+1 == len(final_stage1):
                p(f"  Scored {i+1}/{len(final_stage1)}, "
                  f"elapsed: {time.time()-t3:.1f}s")

    p(f"  Stage 3 complete: {time.time()-t3:.1f}s")
    results_list.sort(key=lambda r: r['best_S'])

    # ── Stage 2b: Band diagnostics under best permutation (top candidates) ──
    p(f"\n--- Stage 2b: Band structure under best permutation (top 10 + codebook) ---")
    for r in results_list[:10]:
        perm = r['best_perm']
        is_cb = '***' if frozenset(r['combo']) == CODEBOOK_SET else '   '
        overlaps = []
        for tr in TRANSITIONS:
            shed = [perm[s] for s in tr['shed_slots']]
            gain = [perm[s] for s in tr['gain_slots']]
            if not gain:
                overlaps.append(0); continue
            b_s = set()
            for ps in shed:
                if ps in bands:
                    for bt in bands[ps]: b_s.update(bt)
            b_g = set()
            for pg in gain:
                if pg in bands:
                    for bt in bands[pg]: b_g.update(bt)
            overlaps.append(len(b_s & b_g))
        p(f"  {is_cb} {r['combo']}  S_best={r['best_S']:.4f}  overlaps={overlaps}")

    # ── Stage 4: Combined result ──
    p(f"\n{'='*70}")
    p(f"  STAGE 4: COMBINED RESULT")
    p(f"{'='*70}")

    p(f"\n--- Cascade table ---")
    p(f"  C({len(pool)},8) universe:           {n_candidates:>12,}")
    p(f"  + F1 (sum=840):                {agg_alone[F1]:>12,}")
    p(f"  + F2 (consec-prime):           {agg_pw_f1[F2]:>12,}  (F1∩F2)")
    p(f"  + F3+F5+F6:                    {len(all_survivors):>12,}  (all cheap)")
    p(f"  + F4 (co_sum legal):           {len(cosum_survivors):>12,}")
    p(f"  R1 diagnostic:                 {r1_pass:>12,}  (not a filter)")
    p(f"  R2 diagnostic:                 {r2_match:>12,}  (not a filter)")
    p(f"  F7 diagnostic:                 {f7_pass:>12,}  (not a filter)")

    p(f"\n--- Impedance ranking (top 20 by best-permutation) ---")
    p(f"{'Rk':>3} {'Set':>42} {'S_nat':>8} {'S_best':>8} "
      f"{'Lift':>7} {'Sig':>3} {'CB':>3}")
    p(f"  {'-'*80}")
    for rank, r in enumerate(results_list[:20], 1):
        is_cb = '***' if frozenset(r['combo']) == CODEBOOK_SET else ''
        p(f"  {rank:>2} {str(r['combo']):>42} {r['nat_S']:>8.4f} "
          f"{r['best_S']:>8.4f} {r['lift']:>7.4f} "
          f"{'Y' if r['best_sigs'] else 'N':>3} {is_cb}")

    # Codebook rank
    cb_result = None
    cb_rank_best = None
    for rank, r in enumerate(results_list, 1):
        if frozenset(r['combo']) == CODEBOOK_SET:
            cb_rank_best = rank
            cb_result = r
            break

    cb_rank_nat = None
    results_nat = sorted(results_list, key=lambda r: r['nat_S'])
    for rank, r in enumerate(results_nat, 1):
        if frozenset(r['combo']) == CODEBOOK_SET:
            cb_rank_nat = rank
            break

    p(f"\n--- Codebook result ---")
    if cb_result:
        p(f"  Natural score:  S = {cb_result['nat_S']:.6f}, sigs = {cb_result['nat_sigs']}")
        p(f"  Best score:     S = {cb_result['best_S']:.6f}, sigs = {cb_result['best_sigs']}")
        p(f"  Best perm:      {cb_result['best_perm']}")
        p(f"  Lift:           {cb_result['lift']:.6f}")
        p(f"  Rank (natural): {cb_rank_nat} of {len(results_list)}")
        p(f"  Rank (best):    {cb_rank_best} of {len(results_list)}")
    else:
        p(f"  *** CODEBOOK NOT FOUND IN RESULTS ***")

    # Nearest rivals
    p(f"\n--- Nearest rivals ---")
    for rank, r in enumerate(results_list[:10], 1):
        if frozenset(r['combo']) != CODEBOOK_SET:
            diff = set(r['combo']) - CODEBOOK_SET
            p(f"  Rank {rank}: {r['combo']}")
            p(f"    S_best={r['best_S']:.4f}, sigs={r['best_sigs']}, "
              f"non-codebook: {sorted(diff)}")

    # Non-codebook prime diagnostic
    non_cb_primes = sorted(pr for pr in consec_set if pr not in CODEBOOK_SET)
    p(f"\n--- Non-codebook prime diagnostic ({len(non_cb_primes)} primes) ---")
    for pr in non_cb_primes:
        n_in = sum(1 for s in final_stage1 if pr in s)
        best_s = None
        for r in results_list:
            if pr in r['combo']:
                best_s = r['best_S']
                break
        p(f"  {pr:>3}: in {n_in} survivors"
          f"{f', best S={best_s:.4f}' if best_s is not None else ''}")

    # 271 neighbourhood
    p(f"\n--- 271 neighbourhood ---")
    for pr in [263, 269, 271, 277, 281, 283]:
        in_pool = pr in pool_set_primes
        in_consec = pr in consec_set
        in_unique = pr in unique_band
        p(f"  {pr}: A279857={'Y' if in_pool else 'N'}, "
          f"consec={'Y' if in_consec else 'N'}, "
          f"unique={'Y' if in_unique else 'N'}")

    # ── Gate checks ──
    gate_a = len(results_list) > 0 and agg_total == n_candidates
    gate_b = cb_in_survivors
    gate_c = cb_rank_best is not None and cb_rank_nat is not None
    gate_d = len(results_list) > 0
    gate_e = True  # N1-N5 all reported above
    gate_f = True  # rivals listed
    gate_g = len(results_list) == len(final_stage1)

    p(f"\n--- Gate checks ---")
    p(f"  A  Cascade complete, combos={agg_total:,} (expected {n_candidates:,}): "
      f"{'PASS' if gate_a else 'FAIL'}")
    p(f"  B  Codebook in survivors: {'PASS' if gate_b else 'FAIL'}")
    p(f"  C  Codebook rank reported: {'PASS' if gate_c else 'FAIL'}")
    p(f"  D  Rivals identified: {'PASS' if gate_d else 'FAIL'}")
    p(f"  E  Negative controls reported: {'PASS' if gate_e else 'FAIL'}")
    p(f"  F  Non-unique rivals listed: {'PASS' if gate_f else 'FAIL'}")
    n_scored = len(results_list)
    n_expected = len(final_stage1)
    p(f"  G  Full 8! for all {n_expected} survivors ({n_scored} scored): "
      f"{'PASS' if gate_g else 'FAIL'}")

    p(f"\n=== TOTAL TIME: {time.time()-t0:.1f}s ===")

if __name__ == '__main__':
    main()
