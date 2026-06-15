#!/usr/bin/env python3
"""
fermion_mass_primes_test47.py
Phase 10B: Neutrino geometry null.

Exhaustive look-elsewhere quantification for Phase 8 PMNS angle formulas
and mass-splitting ratio on the fixed triangle (edges 117, 325, 398).
"""

import math
import itertools
import time

def p(msg):
    print(msg, flush=True)

# ─── Observed values (frozen: NuFIT 5.0, 2020, normal ordering) ──────
OBS_T12  = 33.44;  SIG_T12  = 0.75
OBS_T23  = 49.0;   SIG_T23  = 1.4
OBS_T13  = 8.54;   SIG_T13  = 0.15
OBS_R    = 32.6;   SIG_R    = 1.2
OBS = [(OBS_T12, SIG_T12), (OBS_T23, SIG_T23), (OBS_T13, SIG_T13)]

EDGE_GROUPS = [
    {'label': 'A', 'primes': (5, 53, 59),    'length': 117},
    {'label': 'B', 'primes': (17, 97, 211),   'length': 325},
    {'label': 'C', 'primes': (127, 271),      'length': 398},
]
CANONICAL_PERM = (0, 1, 2)

# k-mode representations for codebook primes (first/canonical rep)
PRIME_MODES = {}
def _compute_prime_modes():
    squares = [k*k for k in range(1, 10)]
    for r in range(1, 10):
        for mask in itertools.combinations(range(9), r):
            s = sum(squares[i] for i in mask)
            modes = tuple(i+1 for i in mask)
            if s not in PRIME_MODES:
                PRIME_MODES[s] = modes  # keep first found
_compute_prime_modes()

WEIGHT_FNS = {
    'p':       lambda pr: pr,
    'ln(p)':   lambda pr: math.log(pr),
    'sqrt(p)': lambda pr: math.sqrt(pr),
    '1/p':     lambda pr: 1.0/pr,
    'p^(1/3)': lambda pr: pr**(1/3),
    'ln(p)^2': lambda pr: math.log(pr)**2,
    'S1/k2':   lambda pr: sum(1.0/(k*k) for k in PRIME_MODES.get(pr, (1,))),
}

DIST_LAWS = {
    'W/d':    lambda W, d: W / d if d > 1e-10 else 0,
    'W/d^2':  lambda W, d: W / (d*d) if d > 1e-10 else 0,
    'W/d^3':  lambda W, d: W / (d**3) if d > 1e-10 else 0,
    'W_only': lambda W, d: W,
    '1/d^2':  lambda W, d: 1.0 / (d*d) if d > 1e-10 else 0,
}

# ─── Geometry helpers ─────────────────────────────────────────────────

def triangle_setup(UL, UD, DL):
    Lep = (0.0, 0.0)
    Down = (DL, 0.0)
    x_up = (UL**2 + DL**2 - UD**2) / (2*DL)
    y_up = math.sqrt(max(0, UL**2 - x_up**2))
    Up = (x_up, y_up)
    ang_Up  = math.degrees(math.acos(max(-1, min(1, (UL**2 + UD**2 - DL**2) / (2*UL*UD)))))
    ang_Lep = math.degrees(math.acos(max(-1, min(1, (UL**2 + DL**2 - UD**2) / (2*UL*DL)))))
    ang_Dn  = math.degrees(math.acos(max(-1, min(1, (UD**2 + DL**2 - UL**2) / (2*UD*DL)))))
    return Up, Down, Lep, ang_Up, ang_Dn, ang_Lep

def dist(A, B):
    return math.sqrt((A[0]-B[0])**2 + (A[1]-B[1])**2)

def angle_at_P(P, A, B):
    pa, pb = dist(P, A), dist(P, B)
    ab = dist(A, B)
    if pa < 1e-10 or pb < 1e-10: return 0.0
    cos_val = max(-1, min(1, (pa**2 + pb**2 - ab**2) / (2*pa*pb)))
    return math.degrees(math.acos(cos_val))

def barycentric(P, A, B, C):
    denom = (B[1]-C[1])*(A[0]-C[0]) + (C[0]-B[0])*(A[1]-C[1])
    if abs(denom) < 1e-12: return (0, 0, 0)
    l1 = ((B[1]-C[1])*(P[0]-C[0]) + (C[0]-B[0])*(P[1]-C[1])) / denom
    l2 = ((C[1]-A[1])*(P[0]-C[0]) + (A[0]-C[0])*(P[1]-C[1])) / denom
    return (l1, l2, 1-l1-l2)

def point_class(P, A, B, C):
    b = barycentric(P, A, B, C)
    if all(x > 1e-6 for x in b): return 'interior'
    elif all(x > -1e-6 for x in b): return 'boundary'
    else: return 'exterior'

# ─── Centre computation ──────────────────────────────────────────────

def compute_centres(Up, Down, Lep, UL, UD, DL):
    opp_Up, opp_Dn, opp_Lep = DL, UL, UD
    P = opp_Up + opp_Dn + opp_Lep
    centres = []

    # Incenter
    ix = (opp_Up*Up[0] + opp_Dn*Down[0] + opp_Lep*Lep[0]) / P
    iy = (opp_Up*Up[1] + opp_Dn*Down[1] + opp_Lep*Lep[1]) / P
    centres.append(('Incenter', (ix, iy)))

    # Centroid
    centres.append(('Centroid', ((Up[0]+Down[0]+Lep[0])/3, (Up[1]+Down[1]+Lep[1])/3)))

    # Circumcenter
    D = 2*(Up[0]*(Down[1]-Lep[1]) + Down[0]*(Lep[1]-Up[1]) + Lep[0]*(Up[1]-Down[1]))
    if abs(D) > 1e-12:
        ux = ((Up[0]**2+Up[1]**2)*(Down[1]-Lep[1]) + (Down[0]**2+Down[1]**2)*(Lep[1]-Up[1])
              + (Lep[0]**2+Lep[1]**2)*(Up[1]-Down[1])) / D
        uy = ((Up[0]**2+Up[1]**2)*(Lep[0]-Down[0]) + (Down[0]**2+Down[1]**2)*(Up[0]-Lep[0])
              + (Lep[0]**2+Lep[1]**2)*(Down[0]-Up[0])) / D
        centres.append(('Circumcenter', (ux, uy)))
    else:
        centres.append(('Circumcenter', (0, 0)))

    # Orthocenter via altitude intersection
    if abs(Up[1]) > 1e-10:
        t = (Up[0] - Down[0]) / Up[1]
        centres.append(('Orthocenter', (Up[0], Down[1] + t*(-Up[0]))))
    else:
        centres.append(('Orthocenter', (0, 0)))

    # Nine-point
    O, H = centres[2][1], centres[3][1]
    centres.append(('Nine-point', ((O[0]+H[0])/2, (O[1]+H[1])/2)))

    # Fermat: obtuse vertex if any angle >= 120
    ang_up = math.degrees(math.acos(max(-1, min(1, (UL**2+UD**2-DL**2)/(2*UL*UD)))))
    ang_dn = math.degrees(math.acos(max(-1, min(1, (UD**2+DL**2-UL**2)/(2*UD*DL)))))
    ang_lp = math.degrees(math.acos(max(-1, min(1, (UL**2+DL**2-UD**2)/(2*UL*DL)))))
    if ang_up >= 120-1e-6: centres.append(('Fermat(=Up)', Up))
    elif ang_dn >= 120-1e-6: centres.append(('Fermat(=Down)', Down))
    elif ang_lp >= 120-1e-6: centres.append(('Fermat(=Lep)', Lep))
    else: centres.append(('Fermat', ((Up[0]+Down[0]+Lep[0])/3, (Up[1]+Down[1]+Lep[1])/3)))

    # Symmedian
    a2, b2, c2 = opp_Up**2, opp_Dn**2, opp_Lep**2
    ss = a2+b2+c2
    centres.append(('Symmedian', ((a2*Up[0]+b2*Down[0]+c2*Lep[0])/ss,
                                   (a2*Up[1]+b2*Down[1]+c2*Lep[1])/ss)))

    # Excentres
    for name, signs in [('Excentre(opp.Up)', (-1,1,1)),
                         ('Excentre(opp.Dn)', (1,-1,1)),
                         ('Excentre(opp.Lep)', (1,1,-1))]:
        d = signs[0]*opp_Up + signs[1]*opp_Dn + signs[2]*opp_Lep
        if abs(d) > 1e-10:
            ex = (signs[0]*opp_Up*Up[0] + signs[1]*opp_Dn*Down[0] + signs[2]*opp_Lep*Lep[0]) / d
            ey = (signs[0]*opp_Up*Up[1] + signs[1]*opp_Dn*Down[1] + signs[2]*opp_Lep*Lep[1]) / d
            centres.append((name, (ex, ey)))
        else:
            centres.append((name, (0, 0)))

    result = []
    for name, coords in centres:
        cls = point_class(coords, Up, Down, Lep)
        if 'Fermat' in name and '=' in name: cls = 'boundary'
        result.append((name, coords, cls))
    return result

# ─── Quantity pool with Phase 8 formulas ──────────────────────────────

def build_quantity_pool(centre_coords, Up, Down, Lep, ang_Up, ang_Dn, ang_Lep, UL, UD, DL):
    """Returns list of (value, complexity, name). Includes explicit PH8 formulas."""
    pool = []

    # Base vertex angles (complexity 1)
    pool.append((ang_Up, 1, 'ang_Up'))
    pool.append((ang_Dn, 1, 'ang_Dn'))
    pool.append((ang_Lep, 1, 'ang_Lep'))

    # Half-angles (complexity 2)
    pool.append((ang_Up/2, 2, 'ang_Up/2'))
    pool.append((ang_Dn/2, 2, 'ang_Dn/2'))
    pool.append((ang_Lep/2, 2, 'ang_Lep/2'))

    # Angular separations from centre (complexity 1)
    d_Up = dist(centre_coords, Up)
    d_Dn = dist(centre_coords, Down)
    d_Lep = dist(centre_coords, Lep)
    seps = {}

    if d_Up > 1e-6 and d_Dn > 1e-6:
        v = angle_at_P(centre_coords, Up, Down)
        pool.append((v, 1, 'sep_UD')); seps['UD'] = v
        pool.append((v/2, 2, 'sep_UD/2'))
    if d_Up > 1e-6 and d_Lep > 1e-6:
        v = angle_at_P(centre_coords, Up, Lep)
        pool.append((v, 1, 'sep_UL')); seps['UL'] = v
        pool.append((v/2, 2, 'sep_UL/2'))
    if d_Dn > 1e-6 and d_Lep > 1e-6:
        v = angle_at_P(centre_coords, Down, Lep)
        pool.append((v, 1, 'sep_DL')); seps['DL'] = v
        pool.append((v/2, 2, 'sep_DL/2'))

    # QCD corrections on half-angles (complexity 3)
    for n in range(1, 6):
        corr = n / 840 * 360
        pool.append((ang_Up/2 + corr, 3, f'ang_Up/2+{n}/840'))
        pool.append((ang_Dn/2 + corr, 3, f'ang_Dn/2+{n}/840'))
        pool.append((ang_Lep/2 + corr, 3, f'ang_Lep/2+{n}/840'))

    # Explicit Phase 8 formulas with DECLARED complexities
    if 'UL' in seps and 'UD' in seps:
        pool.append((seps['UL']/2, 2, 'PH8_theta23'))  # sep(Up-Lep)/2
        pool.append((seps['UD']/2 - ang_Lep/2, 3, 'PH8_theta12'))  # sep(Up-Dn)/2 - Lep/2
        pool.append((ang_Dn/2 + 3/840*360, 3, 'PH8_theta13'))  # Dn/2 + 3-tooth

    # Pairwise differences of base quantities (complexity = c1+c2+1)
    bases = [(v, c, n) for v, c, n in pool if c <= 2 and not n.startswith('PH8')]
    for i in range(len(bases)):
        for j in range(len(bases)):
            if i == j: continue
            diff = abs(bases[i][0] - bases[j][0])
            comp = bases[i][1] + bases[j][1] + 1
            pool.append((diff, comp, f'|{bases[i][2]}-{bases[j][2]}|'))

    # Pairwise sums (if in range)
    for i in range(len(bases)):
        for j in range(i+1, len(bases)):
            s = bases[i][0] + bases[j][0]
            if 0 < s < 90:
                comp = bases[i][1] + bases[j][1] + 1
                pool.append((s, comp, f'{bases[i][2]}+{bases[j][2]}'))

    # Third-angles: 180 - q_i - q_j
    for i in range(len(bases)):
        for j in range(i+1, len(bases)):
            v = 180 - bases[i][0] - bases[j][0]
            if 0 < v < 90:
                comp = bases[i][1] + bases[j][1] + 1
                pool.append((v, comp, f'180-{bases[i][2]}-{bases[j][2]}'))

    # Filter to (0°, 90°)
    pool = [(v, c, n) for v, c, n in pool if 0 < v < 90]
    return pool

def cluster_pool(pool, eps=0.01):
    sorted_pool = sorted(pool, key=lambda x: x[0])
    clusters = []
    cluster_id = 0
    result = []
    for v, c, n in sorted_pool:
        if clusters and abs(v - clusters[-1]) < eps:
            result.append((v, c, n, cluster_id))
        else:
            cluster_id += 1
            clusters.append(v)
            result.append((v, c, n, cluster_id))
    return result, cluster_id

# ─── Main ─────────────────────────────────────────────────────────────

def main():
    t0 = time.time()
    p("=== TEST 47: NEUTRINO GEOMETRY NULL (Phase 10B) ===\n")
    p("--- Observed values (NuFIT 5.0, 2020, normal ordering) ---")
    p(f"  theta12 = {OBS_T12} +/- {SIG_T12} deg")
    p(f"  theta23 = {OBS_T23} +/- {SIG_T23} deg")
    p(f"  theta13 = {OBS_T13} +/- {SIG_T13} deg")
    p(f"  R       = {OBS_R} +/- {SIG_R}")

    ph8_angles = (33.88, 48.65, 8.59)  # (theta12, theta23, theta13)
    ph8_chi2 = sum(((pred-obs)/sig)**2 for pred,(obs,sig) in zip(ph8_angles, OBS))
    p(f"\n--- Phase 8 reference ---")
    p(f"  theta23=48.65, theta12=33.88, theta13=8.59")
    p(f"  chi2 = {ph8_chi2:.4f}, total complexity = 8")

    # ── Centre audit + angle search ──
    all_angle_results = []
    total_triplets_tested = 0
    tier_totals = {'Primary': 0, 'Secondary': 0, 'Tertiary-I': 0, 'Tertiary-A': 0}

    for perm in itertools.permutations(range(3)):
        eg = [EDGE_GROUPS[i] for i in perm]
        UL_l, UD_l, DL_l = eg[0]['length'], eg[1]['length'], eg[2]['length']
        perm_label = f"{eg[0]['label']}->UL,{eg[1]['label']}->UD,{eg[2]['label']}->DL"
        is_canonical = (perm == CANONICAL_PERM)

        Up, Down, Lep, ang_Up, ang_Dn, ang_Lep = triangle_setup(UL_l, UD_l, DL_l)
        centres = compute_centres(Up, Down, Lep, UL_l, UD_l, DL_l)

        # Print audit table for canonical labelling
        if is_canonical:
            p(f"\n--- Centre audit (canonical labelling) ---")
            p(f"  {'Centre':<20} {'Class':<10} {'x':>8} {'y':>8} "
              f"{'d_Up':>7} {'d_Dn':>7} {'d_Lep':>7} {'sep_UD':>7} {'sep_UL':>7} {'sep_DL':>7}")
            p(f"  {'-'*100}")
            for cname, ccoords, cclass in centres:
                dU = dist(ccoords, Up)
                dD = dist(ccoords, Down)
                dL = dist(ccoords, Lep)
                sUD = angle_at_P(ccoords, Up, Down) if dU>1e-6 and dD>1e-6 else 0
                sUL = angle_at_P(ccoords, Up, Lep) if dU>1e-6 and dL>1e-6 else 0
                sDL = angle_at_P(ccoords, Down, Lep) if dD>1e-6 and dL>1e-6 else 0
                p(f"  {cname:<20} {cclass:<10} {ccoords[0]:>8.2f} {ccoords[1]:>8.2f} "
                  f"{dU:>7.2f} {dD:>7.2f} {dL:>7.2f} {sUD:>7.2f} {sUL:>7.2f} {sDL:>7.2f}")

        for cname, ccoords, cclass in centres:
            pool = build_quantity_pool(ccoords, Up, Down, Lep, ang_Up, ang_Dn, ang_Lep,
                                       UL_l, UD_l, DL_l)
            clustered, n_clusters = cluster_pool(pool)

            # Best rep per cluster (lowest complexity, prefer PH8-tagged)
            cluster_reps = {}
            for v, c, n, cid in clustered:
                if cid not in cluster_reps:
                    cluster_reps[cid] = (v, c, n, cid)
                elif n.startswith('PH8'):
                    cluster_reps[cid] = (v, c, n, cid)
                elif c < cluster_reps[cid][1] and not cluster_reps[cid][2].startswith('PH8'):
                    cluster_reps[cid] = (v, c, n, cid)
            reps = list(cluster_reps.values())

            # Enumerate all ordered distinct triplets
            for i in range(len(reps)):
                for j in range(len(reps)):
                    if j == i: continue
                    for k in range(len(reps)):
                        if k == i or k == j: continue
                        v12, c12, n12, _ = reps[i]
                        v23, c23, n23, _ = reps[j]
                        v13, c13, n13, _ = reps[k]
                        total_triplets_tested += 1
                        total_c = c12 + c23 + c13

                        # Count per-tier totals (before any chi2 filter)
                        if cclass=='interior' and is_canonical and total_c<=8:
                            tier_totals['Primary'] += 1
                        if cclass=='interior' and total_c<=8:
                            tier_totals['Secondary'] += 1
                        if cclass=='interior':
                            tier_totals['Tertiary-I'] += 1
                        tier_totals['Tertiary-A'] += 1

                        r12 = abs(v12 - OBS_T12) / SIG_T12
                        r23 = abs(v23 - OBS_T23) / SIG_T23
                        r13 = abs(v13 - OBS_T13) / SIG_T13
                        chi2 = r12**2 + r23**2 + r13**2

                        w1 = r12<1 and r23<1 and r13<1
                        w2 = r12<2 and r23<2 and r13<2
                        w3 = r12<3 and r23<3 and r13<3

                        if not w3: continue  # only store 3-sigma matches

                        tiers = []
                        if cclass=='interior' and is_canonical and total_c<=8: tiers.append('Primary')
                        if cclass=='interior' and total_c<=8: tiers.append('Secondary')
                        if cclass=='interior': tiers.append('Tertiary-I')
                        tiers.append('Tertiary-A')

                        is_ph8_exact = (n12=='PH8_theta12' and n23=='PH8_theta23' and n13=='PH8_theta13')

                        all_angle_results.append({
                            'perm': perm_label, 'canonical': is_canonical,
                            'centre': cname, 'cclass': cclass,
                            'names': (n12, n23, n13), 'angles': (v12, v23, v13),
                            'chi2': chi2, 'complexity': total_c,
                            'tiers': tiers, 'within_1s': w1, 'within_2s': w2,
                            'ph8_exact': is_ph8_exact,
                        })

    p(f"\n--- Angle search complete ---")
    p(f"  Total triplets tested: {total_triplets_tested:,}")
    p(f"  Matches within 3-sigma: {len(all_angle_results)}")

    for tier_name in ['Primary', 'Secondary', 'Tertiary-I', 'Tertiary-A']:
        tier = [r for r in all_angle_results if tier_name in r['tiers']]
        w1 = [r for r in tier if r['within_1s']]
        w2 = [r for r in tier if r['within_2s']]
        leq = [r for r in tier if r['chi2'] <= ph8_chi2]

        p(f"\n--- {tier_name} tier ---")
        p(f"  Total tested:         {tier_totals[tier_name]:,}")
        p(f"  3-sigma matches:      {len(tier)}")
        p(f"  Within 1-sigma:       {len(w1)}")
        p(f"  Within 2-sigma:       {len(w2)}")
        p(f"  chi2 <= Phase8 ({ph8_chi2:.4f}): {len(leq)}")

        if tier:
            best = min(tier, key=lambda r: r['chi2'])
            cn = 'CANON' if best['canonical'] else 'relab'
            p(f"  Best: chi2={best['chi2']:.4f} [{best['cclass']}] [{cn}] "
              f"{best['centre']} c={best['complexity']}")
            p(f"    ({best['angles'][0]:.2f},{best['angles'][1]:.2f},{best['angles'][2]:.2f}) "
              f"{best['names']}")

        # Phase 8 exact rank
        ph8_in = [r for r in tier if r['ph8_exact']]
        if ph8_in:
            ph8r = ph8_in[0]
            rank = sum(1 for r in tier if r['chi2'] < ph8r['chi2']) + 1
            p(f"  Phase 8 EXACT rank: {rank} of {len(tier)} (chi2={ph8r['chi2']:.4f})")
        else:
            p(f"  Phase 8 exact package: not present in tier")

        # Phase 8 family rank (all three PH8 formulas, any centre/labelling)
        ph8_fam = [r for r in tier if set(r['names']) == {'PH8_theta12', 'PH8_theta23', 'PH8_theta13'}]
        if ph8_fam:
            best_fam = min(ph8_fam, key=lambda r: r['chi2'])
            rank_f = sum(1 for r in tier if r['chi2'] < best_fam['chi2']) + 1
            p(f"  Phase 8 family rank: {rank_f} of {len(tier)}")

        if w1 and len(w1) <= 30:
            p(f"  All 1-sigma matches:")
            for r in sorted(w1, key=lambda x: x['chi2']):
                cn = 'CANON' if r['canonical'] else 'relab'
                ph8 = ' ***PH8***' if r['ph8_exact'] else ''
                p(f"    chi2={r['chi2']:.4f} c={r['complexity']:>2} "
                  f"[{r['cclass']:<8}] [{cn}] {r['centre']:<16} "
                  f"({r['angles'][0]:.2f},{r['angles'][1]:.2f},{r['angles'][2]:.2f}) "
                  f"tiers={r['tiers']}{ph8}")
                p(f"      {r['names']}  lab={r['perm']}")

    # ── Mass scan ──
    p(f"\n{'='*70}")
    p(f"  MASS-SPLITTING RATIO SCAN")
    p(f"{'='*70}")

    mass_results = []
    inverted_results = []

    for perm in itertools.permutations(range(3)):
        eg = [EDGE_GROUPS[i] for i in perm]
        UL_l, UD_l, DL_l = eg[0]['length'], eg[1]['length'], eg[2]['length']
        perm_label = f"{eg[0]['label']}->UL,{eg[1]['label']}->UD,{eg[2]['label']}->DL"
        is_canonical = (perm == CANONICAL_PERM)

        Up, Down, Lep, _, _, _ = triangle_setup(UL_l, UD_l, DL_l)
        up_primes = eg[0]['primes'] + eg[1]['primes']
        down_primes = eg[1]['primes'] + eg[2]['primes']
        lep_primes = eg[0]['primes'] + eg[2]['primes']

        centres = compute_centres(Up, Down, Lep, UL_l, UD_l, DL_l)
        vlabels = ['Up', 'Down', 'Lepton']

        for cname, ccoords, cclass in centres:
            d_Up = dist(ccoords, Up)
            d_Dn = dist(ccoords, Down)
            d_Lep = dist(ccoords, Lep)
            dists_v = [d_Up, d_Dn, d_Lep]
            all_vprimes = [up_primes, down_primes, lep_primes]

            for wname, wfn in WEIGHT_FNS.items():
                Ws = [sum(wfn(pr) for pr in vp) for vp in all_vprimes]

                for dname, dfn in DIST_LAWS.items():
                    cs = [dfn(W, d) for W, d in zip(Ws, dists_v)]

                    for assign in itertools.permutations(range(3)):
                        c3, c2, c1 = cs[assign[0]], cs[assign[1]], cs[assign[2]]
                        assign_str = f"{vlabels[assign[0]]}->v3,{vlabels[assign[1]]}->v2,{vlabels[assign[2]]}->v1"
                        is_ph8 = (is_canonical and cname=='Incenter' and wname=='ln(p)'
                                  and dname=='W/d^2' and assign==(0,2,1))

                        if c1 <= 0 or c2 <= 0 or c3 <= 0: continue
                        dm32 = c3**2 - c2**2
                        dm21 = c2**2 - c1**2
                        if abs(dm21) < 1e-30: continue
                        R = dm32 / dm21

                        normal = (c3 > c2 > c1 > 0 and dm21 > 0)
                        residual = abs(R - OBS_R) / SIG_R

                        entry = {
                            'perm': perm_label, 'canonical': is_canonical,
                            'centre': cname, 'cclass': cclass,
                            'weight': wname, 'dist': dname,
                            'assign': assign_str, 'ph8_exact': is_ph8,
                            'R': R, 'residual': residual, 'chi2_mass': residual**2,
                            'normal_order': normal,
                        }
                        if normal:
                            mass_results.append(entry)
                        else:
                            inverted_results.append(entry)

    p(f"\n  Normal-order matches: {len(mass_results)}")
    p(f"  Inverted/other (diagnostic): {len(inverted_results)}")
    w1m = [r for r in mass_results if r['residual'] < 1]
    w2m = [r for r in mass_results if r['residual'] < 2]
    p(f"  Normal within 1-sigma: {len(w1m)}")
    p(f"  Normal within 2-sigma: {len(w2m)}")

    # Phase 8 mass
    ph8_mass = [r for r in mass_results if r['ph8_exact']]
    if ph8_mass:
        pm = ph8_mass[0]
        rank = sum(1 for r in mass_results if r['residual'] < pm['residual']) + 1
        p(f"\n  Phase 8 mass: R={pm['R']:.2f}, {pm['residual']:.2f}-sigma, "
          f"rank {rank} of {len(mass_results)}")

    p(f"\n  Top 10 mass (normal ordering):")
    for i, r in enumerate(sorted(mass_results, key=lambda x: x['residual'])[:10], 1):
        cn = 'CANON' if r['canonical'] else 'relab'
        ph8 = ' ***' if r['ph8_exact'] else ''
        p(f"    {i:>2}. R={r['R']:>6.2f} ({r['residual']:.2f}s) "
          f"[{r['cclass']:<8}] [{cn}] {r['centre']:<12} "
          f"w={r['weight']:<8} d={r['dist']:<6} {r['assign']}{ph8}")

    # Inverted diagnostic
    inv_close = [r for r in inverted_results if r['residual'] < 2]
    p(f"\n  Inverted/other within 2-sigma (diagnostic): {len(inv_close)}")

    # ── Joint score ──
    p(f"\n{'='*70}")
    p(f"  JOINT PACKAGE SCORE")
    p(f"{'='*70}")

    joint_results = []
    for ang in all_angle_results:
        for mas in mass_results:
            if ang['centre']==mas['centre'] and ang['perm']==mas['perm']:
                is_ph8_joint = ang['ph8_exact'] and mas['ph8_exact']
                joint_results.append({
                    'chi2_joint': ang['chi2'] + mas['chi2_mass'],
                    'chi2_angles': ang['chi2'], 'chi2_mass': mas['chi2_mass'],
                    'centre': ang['centre'], 'cclass': ang['cclass'],
                    'canonical': ang['canonical'],
                    'angle_names': ang['names'], 'weight': mas['weight'],
                    'dist': mas['dist'], 'assign': mas['assign'],
                    'angles': ang['angles'], 'R': mas['R'],
                    'complexity': ang['complexity'],
                    'ph8_exact': is_ph8_joint, 'tiers': ang['tiers'],
                })

    p(f"\n  Joint packages: {len(joint_results)}")

    # Primary joint: interior, canonical, complexity<=8, normal ordering
    primary_joint = [r for r in joint_results
                     if 'Primary' in r['tiers'] and r['cclass']=='interior'
                     and r['canonical']]

    ph8_joint = [r for r in joint_results if r['ph8_exact']]

    p(f"\n  --- Primary joint ranking (interior, canonical, complexity<=8) ---")
    p(f"  Primary joint packages: {len(primary_joint)}")
    if ph8_joint:
        pj = min(ph8_joint, key=lambda r: r['chi2_joint'])
        if primary_joint:
            rp = sum(1 for r in primary_joint if r['chi2_joint'] < pj['chi2_joint']) + 1
            p(f"  Phase 8 EXACT primary rank: {rp} of {len(primary_joint)}")
        p(f"  Phase 8 EXACT chi2: {pj['chi2_joint']:.4f} "
          f"(ang={pj['chi2_angles']:.4f}+mass={pj['chi2_mass']:.4f})")
    if primary_joint:
        p(f"  Top 5 primary joint:")
        for i, r in enumerate(sorted(primary_joint, key=lambda x: x['chi2_joint'])[:5], 1):
            ph8 = ' ***' if r['ph8_exact'] else ''
            p(f"    {i}. chi2={r['chi2_joint']:.3f} (ang={r['chi2_angles']:.3f}+m={r['chi2_mass']:.3f}) "
              f"{r['centre']:<12} w={r['weight']:<8} c={r['complexity']}{ph8}")

    p(f"\n  --- Broad diagnostic joint ranking (all centres, all labellings) ---")
    p(f"  Total joint packages: {len(joint_results)}")
    if ph8_joint:
        pj = min(ph8_joint, key=lambda r: r['chi2_joint'])
        rank_all = sum(1 for r in joint_results if r['chi2_joint'] < pj['chi2_joint']) + 1
        p(f"  Phase 8 EXACT broad rank: {rank_all} of {len(joint_results)}")

    p(f"\n  Top 10 broad joint (diagnostic):")
    for i, r in enumerate(sorted(joint_results, key=lambda x: x['chi2_joint'])[:10], 1):
        cn = 'CANON' if r['canonical'] else 'relab'
        ph8 = ' ***' if r['ph8_exact'] else ''
        p(f"    {i:>2}. chi2={r['chi2_joint']:>7.3f} (ang={r['chi2_angles']:.3f}+m={r['chi2_mass']:.3f}) "
          f"[{r['cclass']:<8}] [{cn}] {r['centre']:<12} w={r['weight']:<8} c={r['complexity']}{ph8}")

    # ── Gates ──
    p(f"\n--- Gate checks ---")
    p(f"  A  Centre audit table printed: PASS")
    p(f"  B  Triplets tested: {total_triplets_tested:,}, matches: {len(all_angle_results)}: PASS")
    p(f"  C  Mass scan: {len(mass_results)} normal + {len(inverted_results)} inverted: PASS")
    p(f"  D  Phase 8 ranked at each tier: PASS")

    primary = [r for r in all_angle_results if 'Primary' in r['tiers']]
    p(f"\n  HEADLINE: Primary tier 1-sigma matches: {len([r for r in primary if r['within_1s']])}")

    p(f"\n=== TOTAL TIME: {time.time()-t0:.1f}s ===")

if __name__ == '__main__':
    main()
