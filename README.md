# Prime Model of Fermions

Scripts and data supporting the paper:

**"A Prime Model of Fermions"**

Eddie George (2026)

**Version:** 1.0
**DOI:** [10.5281/zenodo.20698434](https://doi.org/10.5281/zenodo.20698434)
**License:** Creative Commons Attribution 4.0 International (CC-BY 4.0). Copyright © 2026 Eddie George.

**Series:** This is the supporting repository for Paper 2 of a triptych.

1. **Fermion Mass Hierarchies from Prime Stability Conditions.** Repository: [fermion-mass-primes](https://github.com/eddiegeorge/fermion-mass-primes)
2. **A Prime Model of Fermions** (this repository).
3. **A Prime Origin of Fermion Structure.**

## Contents

- `scripts/` — Python implementations of the eight tests cited in the paper
- `logs/` — output logs of the runs cited in the paper
- `figures/` — figures used in the paper
- `README.md` — this file

## Running the scripts

Requirements: Python 3.10+ with `numpy`, `scipy`.

To reproduce the results cited in the paper, from the repository root:

```
cd scripts
python fermion_mass_primes_test40.py
python fermion_mass_primes_test41.py
python fermion_mass_primes_test42.py
python fermion_mass_primes_test43.py
python fermion_mass_primes_test44.py
python fermion_mass_primes_test45.py
python fermion_mass_primes_test46.py
python fermion_mass_primes_test47.py
```

Each script prints results to standard output. To capture the output in the `logs/` format used in this repo:

```
python fermion_mass_primes_test46.py | tee ../logs/test46_$(date +%Y%m%d_%H%M%S).log
```

Runtimes (approximate, on a 4-core VM):

- test40–test45: sub-second to a few seconds
- test46: ~6 minutes (exhaustive C(43,8) ≈ 145M enumeration)
- test47: ~2–5 minutes (exhaustive 35.7M projection triplet enumeration)

## Scripts

### fermion_mass_primes_test40.py
**Phase 1: Cavity-grounded impedance model — v2.0 control and ontological translation.**

Iteration 3 returns the probability model to v2.0's additive Z-space form. Test 40 implements Phase 1 in two stages: an exact v2.0 control, then an ontological translation confirming that every v2.0 input is a cavity-derived quantity.

##### Stage 1: v2.0 arithmetic control

Reproduce the v2.0 impedance equation exactly. This validates the Z-term implementation.

```
P(i→j) = exp(−Z_ij) / Σ_j exp(−Z_ij)
Z_total = Z_bal + Z_ch + Z_R + Z_a + Z_deg + Z_tp
```

Exact v2.0 formulas:

```
Z_bal = α_s × |Σ p_shed − Σ p_gain|

Z_ch  = −ln(K)
  K = (1/√N) Σ_p∈S Σ_q∈G exp(−α_s × |p−q| × max(p,q)/p_max)
  N = n_shed × n_gain
  Pure shedding: Z_ch = 0

Z_R (engaging, Δa > 0):
  Z_R = (π/2) × |Δa| × J
  J = Σ ln(p_shared) / (Σ ln(p_shed) + Σ ln(p_gain))

Z_R (releasing, Δa < 0):
  Z_R = (π/2) × |Δa| × (1 − α_s × Σp_rewired/P + α_s × Σp_gain × |Δa| / P)
  P = 840

Z_R (Δa = 0): Z_R = 0

Z_a = (α_s/P) × (Σ p_shed × a_src + Σ p_gain × a_tgt)

Z_deg = (α_s/P) × |Δa| × (Σ p_shed + Σ p_gain) × J

Z_tp = (α_s / Q) × Σ p_shed
  Q = n_src × (n_src − 1)
  Applies to pure shedding at Δa = 0 only
```

Gate A: S ≈ 0.031, matching v2.0 published result.

##### Stage 2: Ontological translation

The same computation as Stage 1, but with every input defined through the cavity ontology. The numerical result must be identical — the cavity does not change the numbers, it explains them.

Ontological grounding for each v2.0 input:

```
v2.0 input          Cavity interpretation
──────────────      ───────────────────────────────────────────────
p                   E_p = Σk² over occupied modes of packet ψ_p
Σ p                 Σ E_p = total packet energy over shed/gain set
|p − q|             |E_p − E_q| = energy gap between cavity packets
max(p,q)/p_max      max(E_p, E_q)/E_max = relative packet energy position
ln(p)               ln(E_p) = log packet energy (rewiring load)
P = 840             Σ E_p over codebook = total codebook energy
p_max = 271         E_max = max codebook packet energy
a_src, a_tgt        substrate sector coupling levels
Δa                  substrate sector transition
n_shed, n_gain      packet configuration change counts
Q = n(n−1)          directed mode-pair count (ring quality factor)
α_s                 strong coupling constant
π/2                 geometric phase for orthogonal substrate dimensions
```

Stage 2 computes E_p = Σk² from the cavity mode representations for each codebook prime and confirms E_p = p identically. Then runs the full impedance model using E_p, confirming S matches Stage 1 exactly.

Gate B: Stage 2 reproduces Stage 1 identically.

##### Report format

For both stages, print per-transition Z-term breakdown:

```
trans  Z_bal  Z_ch   Z_R   Z_a  Z_deg  Z_tp  Z_total
u→d    ...    ...    ...   ...  ...    ...   ...
...
```

Then P_ij matrix, observed P_obs matrix, per-element errors, row scores, total S, and row signature.

Stage 2 additionally prints the E_p = Σk² derivation for each codebook prime, confirming E_p = p.

Comparison table at end:

```
Model              S      u_sc   c_sc   t_sc   Sig
v2.0 control       ...    ...    ...    ...    3/3
Cavity E_p         ...    ...    ...    ...    3/3
```

##### Usage

```
python fermion_mass_primes_test40.py              # both stages
python fermion_mass_primes_test40.py --control     # stage 1 only
python fermion_mass_primes_test40.py --cavity      # stage 2 only
```

##### Parameter status

```
Parameter   Value    Status
─────────   ─────    ──────────────────────────────────────────────
P           840      Derived: Σ E_p over codebook. Fixed by cavity codebook selection.
p_max       271      Derived: max E_p over codebook. Fixed by cavity codebook selection.
α_s         0.118    Identified: strong coupling constant at M_Z. Measured SM constant,
                     identified post-hoc in v2.0, not fitted.
π/2         1.5708   Open: v2.0 calls it "geometric phase for orthogonal substrate
                     dimensions." Not yet derived from cavity geometry. Candidate
                     derivation from C₃ vertex structure or ring quarter-wave condition
                     is a Phase 2 question.
```

No free parameters in the impedance equation. P and p_max are codebook consequences. α_s is a physical constant. π/2 is justified physically but not yet derived from the cavity.

##### Success criteria

| Level | Criteria |
|-------|----------|
| Gate A | v2.0 control reproduces S ≈ 0.031 |
| Gate B | Cavity E_p translation reproduces Stage 1 identically |

##### What this establishes

Phase 1 confirms that the v2.0 impedance model, reinterpreted through the cavity ontology, reproduces the same result. Every input is grounded: p is not an arbitrary prime label but the packet energy Σk² on the triangular loop cavity. The impedance model describes energy transport through the cavity's Kirchhoff network. The bridge between ansatz and ontology is built.

### fermion_mass_primes_test41.py
**Phase 2: Co_sum phase separation — collision diagnostic.**

The translated impedance equation (test40) multi-counts eigenvalues: multiple primes excite the same k-modes. This is physically illegal under Pauli exclusion. Phase 2 proposes that the existing phase structure in the packet states distinguishes shared modes via a sector label r = raw co_sum (sum of the other mode indices in the mask), promoting the basis from |k⟩ to |k, r⟩.

Test 41 is a diagnostic: it computes the phase separation structure and checks whether it resolves the multi-counting problem without disturbing the impedance equation.

##### Five questions

1. For each shared k, are raw co_sum labels distinct across codebook packets?
2. Does this make every fermion state Pauli-legal in |k, r⟩ space?
3. Does E_p = Σk² = p remain unchanged?
4. Does the test40 impedance calculation reproduce identically?
5. How selective is collision-free co_sum structure across the A279857 null pool?

##### Method

**Q1: Co_sum distinctness at shared modes.**

For each codebook prime, compute all A279857 representations (not just the first). For each representation, at each mode k in the mask, compute:

```
r_p(k) = Σ(modes in mask except k)
```

For each mode k shared by two or more codebook primes, report the (prime, representation, co_sum) tuples. Flag any collision: same k and same co_sum from different primes.

Report per-mode: number of primes sharing this mode, all co_sum values, and whether all distinct.

**Q2: Per-fermion Pauli legality.**

For each fermion (u, c, t, d, s, b), enumerate all representation combinations for its assigned primes (most primes have one representation; 53 has 2, 127 has 4, 211 has 4). For each combination, build the occupied set {(k, r)} and check for duplicate occupation.

Report per fermion: total representation combinations, number that are collision-free, and any collisions found with the specific (k, r) pair and the colliding primes.

**Q3: Energy identity.**

For each codebook prime and each of its representations, compute E_p = Σk² over the mask. Confirm E_p = p for all representations.

**Q4: Impedance identity.**

Using the collision-free representation assignment (first/fewest-mode rep for each prime), run the test40 impedance equation. Confirm:

```
S = 0.031387 (identical to test40)
All Z-terms identical to test40
All P_ij identical to test40
```

The phase labels do not enter the impedance equation — it uses E_p only. This gate confirms that phase separation fixes the ontology without disturbing the effective theory.

**Q5: Null pool selectivity.**

Enumerate all A279857 primes expressible on the 9-mode k² basis. For each prime, compute ALL representations.

For each of the C(pool, 8) 8-prime subsets, determine viability: a subset is viable if there exists a representation assignment (one rep per prime) such that no two primes share any (k, r) pair. Enumerate all representation combinations per subset to check this exhaustively.

Report:
```
Pool size (number of A279857 primes on 9 modes)
Number of primes with multiple representations
Total 8-prime subsets: C(pool, 8)
Collision-free (viable) subsets: exact count
Fraction viable
```

##### Usage

```
python fermion_mass_primes_test41.py
```

No flags. Single run computes all five questions.

##### Success criteria

| Gate | Criteria |
|------|----------|
| Q1 | All raw co_sums distinct at every shared mode (first-rep assignment) |
| Q2 | At least one collision-free rep assignment exists for every fermion |
| Q3 | E_p = p for all primes, all representations |
| Q4 | S = 0.031387, identical to test40 |
| Q5 | Reports exact collision-free count and fraction |

##### What this establishes

If Q1–Q4 pass: the phase-extended basis |k, r⟩ with r = raw co_sum resolves the eigenvalue multi-counting problem. The topology is Pauli-legal. The impedance equation is unchanged. Phase 2 succeeds.

Q5 quantifies how much structural selectivity the collision-free constraint provides, independent of the impedance score. This feeds directly into Phase 6 (exhaustive null) as an additional filter.

### fermion_mass_primes_test42.py
**Phase 3: Transfer determinant mass mechanism.**

Phase 2 established that the eight codebook packets are independent in the phase-extended |k, r⟩ basis. Phase 3 treats fermion mass as the determinant of an open-system packet transfer operator. Each phase-separated packet ψ_p carries a transfer gain with determinant p. Because independent packet channels combine multiplicatively, the fermion wedge state has mass factor Πp_i. The logarithmic form ln M = Σ ln p is the additive action representation of the same mechanism.

The key conceptual move: mass is gain, not stored energy. The cavity packet energy E_p = Σk² = p defines the packet label (additive eigenvalues). The transfer determinant det(T_p) = p defines the mass contribution (multiplicative gain). Same number, different physical role.

##### Consistency with established frameworks

The mass-as-transfer-gain interpretation connects to several established frameworks in physics:

- **Lattice QCD transfer matrix.** On the lattice, hadron masses are extracted from the exponential decay of the transfer matrix: T^n ~ exp(−m·n), giving m = −ln(eigenvalue of T). Mass is the logarithm of a transfer matrix eigenvalue. Standard lattice QCD.
- **Froggatt-Nielsen mechanism.** A BSM model where fermion mass hierarchies arise from products of scalar VEVs in a chain of couplings. Mass = product of gain factors along the chain. Explicitly a cascade-gain mechanism in published BSM literature.
- **Effective mass in condensed matter.** The electron's effective mass in a crystal is a transport property determined by band curvature, not intrinsic rest mass. Mass as a transport-derived quantity is standard.
- **Dynamical mass in QCD.** Most of the proton's mass (~99%) comes from gluon field dynamics, not Higgs couplings. Mass as an emergent transport phenomenon is the mainstream explanation for hadronic mass.

Phase 3 does not claim to derive det(T_p) = p from a Lagrangian. It establishes that the identification is formally consistent and that the fermion mass algebra works. The open-system derivation (specifying the gain medium and showing det(T_p) = E_p) is a downstream question.

##### Four questions

1. Does the packet transfer gain g_p = p hold for all 8 codebook packets?
2. Does the fermion mass M_f = Π g_p reproduce the correct mode numbers for all 6 quarks?
3. Do CKM shedding ratios M_parent / M_child equal the product of shed packet gains?
4. Does the test40 impedance calculation reproduce identically?

##### Method

**Q1: Packet transfer gain.**

Define the minimal effective packet-space transfer operator G with eigenvalues g_p = E_p = p. Test42 verifies that this identification makes the fermion mass algebra consistent. For each codebook prime, report:
```
E_p = Σk² over occupied modes = p       (packet energy, from Phase 1)
g_p = E_p = p                            (transfer gain, Phase 3 identification)
ln g_p = ln p                            (additive action contribution)
```

This is not a derivation of det(T_p) = p from open-system dynamics — that remains an open question for later work. It is a formal bridge: the packet energy is identified as the transfer gain.

**Q2: Fermion mass from transfer determinant.**

Construct the 8×8 diagonal transfer gain operator in packet space:
```
G = diag(g_0, g_1, ..., g_7) = diag(5, 17, 53, 59, 97, 127, 211, 271)
```

For each fermion f with active packet slots S_f, compute:
```
M_f = det(G|_{S_f}) = Π_{i ∈ S_f} g_i
ln M_f = Σ_{i ∈ S_f} ln g_i
```

Report the mass mode number M_f and its logarithm for all 6 quarks. Verify the mass hierarchy: M_t > M_b > M_c > M_s > M_d > M_u.

**Q3: Shedding ratios.**

For each CKM transition (src → tgt), the shed packets are S_shed = S_src \ S_tgt and the gained packets are S_gain = S_tgt \ S_src. The mass ratio is:
```
M_src / M_tgt = Π_{i ∈ S_shed} g_i / Π_{i ∈ S_gain} g_i
```

For pure-shedding transitions (no gain), this simplifies to:
```
M_src / M_tgt = Π_{i ∈ S_shed} g_i
```

Report the shedding/gain ratio for all 9 CKM transitions. Verify that generation ratios involve single packet gains:
```
m_c/m_u = g_5 = 127
m_t/m_c = g_7 = 271
m_s/m_d = g_0 = 5
m_b/m_s = g_2 = 53
```

**Q4: Impedance identity.**

Run the test40 impedance equation using E_p = p (unchanged). Confirm S = 0.031387. The transfer gain mechanism operates at the mass level, not the impedance level — the two are independent computations on the same packet energies.

##### Usage

```
python fermion_mass_primes_test42.py
```

##### Success criteria

| Gate | Criteria |
|------|----------|
| Q1 | g_p = E_p = p for all 8 codebook packets |
| Q2 | M_f = Πp_i gives correct mode numbers for all 6 quarks, hierarchy correct |
| Q3 | All 9 CKM shedding ratios consistent with packet gain algebra |
| Q4 | S = 0.031387, identical to test40 |

##### What this establishes

If Q1–Q4 pass: the multiplicative mass mechanism is formally consistent with the phase-extended topology. Packet energy (additive eigenvalues) and transfer gain (multiplicative mass) are identified as the same number p, in different physical roles. The fermion mass is the determinant of the transfer operator restricted to the fermion's active packets. The impedance equation (which handles CKM transitions) operates on the same packet energies but through a separate additive Z-space mechanism.

The two layers of the topology are:
```
Mass layer:       det(T_f) = Πp_i        (multiplicative, transfer gain)
Transition layer:  Z_ij = Σ Z_terms       (additive, impedance costs)
Both use:          E_p = Σk² = p          (packet energy as common input)
```

### fermion_mass_primes_test43.py
**Phase 4: CKM phase connection — impedance magnitudes meet topology phases.**

The v2.0 impedance model gives CKM magnitudes |V_ij|² through real-valued transition costs Z_ij. The full CKM matrix is a complex unitary rotation with three mixing angles and one CP phase δ. The magnitudes determine the angles; the topology must supply the phase.

The natural framework is an open-system complex transition amplitude:

```
A_ij = exp(−Z_ij / 2) · exp(i Φ_ij)

|A_ij|² = exp(−Z_ij)   → impedance probability (real part, Phases 1–3)
Φ_ij                    → CKM phase structure (imaginary part, Phase 4)
```

Phase 4 does not replace the impedance model. It completes it: the impedance gives magnitudes, the topology gives phases, and together they give the full complex CKM matrix.

##### Three questions

1. What mixing angles does the impedance model imply?
2. At the observed δ = 1.139, what phases Φ_ij does the full unitary CKM matrix require?
3. Do those required phases correlate with topology quantities (co_sum, Δa, π/2, winding)?

##### Method

**Q1: Extract mixing angles from impedance magnitudes.**

From the predicted |V_ij| (test40 impedance model), extract the three standard mixing angles:

```
s₁₃ = |V_ub|
s₁₂ = |V_us| / √(1 − s₁₃²)
s₂₃ = |V_cb| / √(1 − s₁₃²)
```

Report predicted vs observed angles. Also report the column unitarity residual as a consistency check (the predicted |V_ij|² are not exactly unitary because the impedance model row-normalises independently).

**Q2: Target reconstruction — build the CKM phase target using observed δ.**

Q2 constructs the standard unitary CKM matrix implied by the three impedance-derived angles and observed δ = 1.139 (UTfit 2023 / Antusch M_Z values, consistent with tests 21/40/41/42). This gives the target phase structure in the standard parameterisation convention, not an exact phase completion of the original nine impedance magnitudes.

For each element V_ij, extract:

```
|V_ij|     (magnitude — already matched by impedance)
Φ_ij = arg(V_ij)   (phase — what the topology must supply)
```

Report the 3×3 phase matrix Φ_ij. These are the standard-convention phases of the unitary CKM matrix implied by the three impedance-derived angles and observed δ. They are a target phase structure, not an exact completion of all nine row-normalised impedance magnitudes. Also compute the Jarlskog invariant:

```
J = Im(V_ud · V_cs · V*_us · V*_cd)
```

for both the impedance-derived matrix and the observed CKM.

**Q3: Topology signal search — test whether topology quantities predict the invariant phase without using observed δ as input.**

All co_sum phase candidates use the Phase 2 admissible canonical representation assignment:

```
5={1,2}, 17={1,4}, 53={2,7}, 59={1,3,7},
97={4,9}, 127={1,3,6,9}, 211={1,4,7,8,9}, 271={4,5,6,7,8,9}
```

Alternative collision-free representation assignments may be reported diagnostically only.

Q3 tests a tiered set of topology phase candidates. The primary candidates are raw co_sum sector quantities, especially per-mode co_sum differences and rephasing-invariant quartet combinations. Secondary candidates include R_shared (the Z_R rewiring factor), R_shared-weighted co_sum phase, ΔlnM, and Δa·π/2. Packet counts, pure-shedding status, winding direction, and energy imbalance are reported as controls.

**Primary phase candidates (co_sum native):**

For each transition, compute per-mode co_sum differences over packet pairs in the shed/gain/shared decomposition:

```
- shared packet slots: same packet appears in source and target;
  expected Δr = 0, used as a control
- shed/gain packet pairs: compare packets that share the same raw k-mode
  and record Δr = r_gain(k) − r_shed(k)
- source/target aggregate: compare the multiset of occupied |k,r⟩ states
  in source vs target
```

Scalar reductions from per-mode data:

```
1. X_shared_sum     = Σ shared-mode Δr (over shed/gain packet pairs)
2. X_shed_gain_sum  = Σ gain r − Σ shed r (aggregate)
3. X_shared_signed_rms = sign(X_shared_sum) · √(Σ shared-mode Δr²)
```

Quartet co_sum phase (topology-side analogue of the CKM rephasing-invariant phase):

```
4. Φ_loop = Φ_ud + Φ_cs − Φ_us − Φ_cd
```

Each X_ij is tested separately under Φ_ij = 2π·X_ij/N. The N-scan (N ∈ {3, 6, 9, 12}) applies to integer co_sum-derived candidates only.

**Secondary phase-coupled candidates:**

```
4. R_shared from Z_R: Σln(shared) / (Σln(shed) + Σln(gain))
5. R_shared-weighted co_sum phase: R_shared · 2π·X/N
6. ΔlnM = Σln(shed) − Σln(gain)
7. Δa · π/2 (literal phase, no N-scan)
```

For Δa·π/2, use the literal phase Δa·π/2. For R_shared and ΔlnM, report as correlation/control diagnostics unless explicitly combined with a co_sum phase (e.g., candidate 5).

**Control / impedance correlates:**

```
8. n_shed, n_gain
9. pure-shedding flag
10. energy imbalance |ΣE_shed − ΣE_gain|
11. winding direction / chirality sign
```

The decisive comparison is not raw arg(V_ij), but rephasing-invariant quartet phases and Jarlskog J. Report raw Φ_ij from Q2 only in the chosen standard convention. The main topology comparison uses the quartet phase Φ_loop and J.

Report the primary quartet Φ_ud + Φ_cs − Φ_us − Φ_cd, then also report all nine 2×2 row/column quartet phases as diagnostics. A genuine topology phase signal should not depend on one arbitrary quartet choice — it should be consistent across all plaquettes. This guards against false positives from a single Cabibbo-sector quartet.

Test whether any primary or secondary candidate predicts J or δ without using the observed δ as an input. For each co_sum candidate, test explicit no-free-parameter phase maps:

```
Φ = 2π · X / N
```

with N ∈ {3, 6, 9, 12} reported separately, plus raw unwrapped X only as a diagnostic. This prevents the script from silently choosing a phase scale, which would turn the test into an implicit fit.

For each topology phase candidate and each N, construct a transition phase matrix Φ_topo, then compute quartet phases using the general form:

```
For rows a,b and columns i,j:

Φ_Q(a,b;i,j) = Φ_ai + Φ_bj − Φ_aj − Φ_bi

J_topo(a,b;i,j) = |V_ai|_Q2 · |V_bj|_Q2 · |V_aj|_Q2 · |V_bi|_Q2 · sin(Φ_Q)
```

Report:

```
Primary:    rows u,c and columns d,s
Diagnostic: all 9 row/column quartets
```

A genuine topology signal should give a coherent sign/magnitude pattern across quartet diagnostics, not necessarily identical raw Φ_Q values. The invariant is J, not the bare quartet angle — the angle depends on the magnitude product used in that plaquette.

Use the standard-unitary magnitudes from Q2 (derived from the three impedance angles) when computing J_topo. Report the original row-normalised impedance magnitudes only as a diagnostic.

Then infer:

```
sin δ_topo = J_topo / (s₁₂ · s₂₃ · s₁₃ · c₁₂ · c₂₃ · c₁₃²)
```

Because this determines only sin δ, compare J_topo primarily. Report δ_topo on the principal branch and also note the π−δ ambiguity. J is the decisive invariant.

Compare J_topo with the observed J for each phase candidate and each N. This makes Q3 non-circular and testable: co_sum candidates produce an actual CKM invariant prediction rather than correlations against convention-dependent element phases.

This is a signal search, not a fit. Observed δ defines the target (Q2); the topology must predict or correlate with the invariant separately (Q3).

##### Usage

```
python fermion_mass_primes_test43.py
```

##### Success criteria

| Gate | Criteria |
|------|----------|
| Q1 | Mixing angles extracted, consistent with observed values within impedance model accuracy |
| Q2 | Target CKM matrix constructed at observed δ, phases Φ_ij tabulated, Jarlskog J computed |
| Q3 | Reports correlation (or absence) between topology quantities and rephasing-invariant phase, without using observed δ as input |

Phase 4 success is graduated:
- **Weak:** no correlation found — the CP phase is a separate layer, document and move on
- **Medium:** correlation exists but with free coefficients — suggestive but not predictive
- **Strong:** the topology predicts J (or δ) from co_sum / Δa / π/2 with no free parameters

Any of these outcomes is a valid result. Phase 4 is a test question, not a premise.

##### What this establishes

If Q3 finds a signal: the topology supplies both magnitudes (through impedance) and phases (through co_sum / winding geometry), completing the full complex CKM matrix from a single physical structure.

If Q3 finds no signal: the impedance model handles magnitudes and the CP phase requires additional structure beyond what Phases 1–3 provide. This is still a clean result — it maps the boundary of what the current topology explains.

### fermion_mass_primes_test44.py
**Phase 5: Z-term ablation — structural X-ray of the impedance equation.**

The v2.0 impedance equation uses six additive Z-terms to compute CKM transition costs. Phase 5 removes each term in turn to measure its load-bearing contribution, identify redundancies, and determine the minimum effective term count. This is diagnostic, not repair — the equation already achieves S = 0.031 with zero fitted parameters.

The six Z-terms are:
```
Z_bal:  energy balance cost          α_s × |ΣE_shed − ΣE_gain|
Z_ch:   channel conductance          −ln(K), pairwise coupling kernel
Z_R:    substrate rewiring cost       π/2 × |Δa| × f(direction)
Z_a:    substrate alignment cost      (α_s/P) × (ΣE_shed·a_src + ΣE_gain·a_tgt)
Z_deg:  degeneracy cost              (α_s/P) × |Δa| × (ΣE_shed+ΣE_gain) × R_shared
Z_tp:   topology/pure-shedding cost  (α_s/Q) × ΣE_shed (pure-shedding neutral only)
```

##### Four questions

1. What is each term's individual contribution to the impedance score?
2. Which terms are load-bearing (removal degrades S significantly)?
3. Are any terms structurally coupled (pairwise interaction after softmax)?
4. What is the minimum effective term count?

##### Method

**Baseline:** all six Z-terms active, using the canonical test40/v2.0 cavity-grounded equation. S_baseline = 0.031387. All Q2–Q4 deltas are measured against this baseline.

**Q1: Per-term contribution profile.**

For each of the 9 CKM transitions, report the value of each Z-term and its contribution to Z_total. Because Z_ch can be negative (conductance relief), report both signed and absolute shares:

```
For each transition report:
- raw signed Z_k for each of 6 terms
- total Z_total
- signed share Z_k / Z_total (diagnostic only)
- absolute share |Z_k| / Σ|Z_j|
```

Also report which terms are structurally zero for which transitions:
```
Z_ch  = 0 for pure-shedding transitions (no gain packets)
Z_R   = 0 for Δa = 0 transitions
Z_deg = 0 for Δa = 0 or pure-shedding transitions
Z_tp  ≠ 0 only for pure-shedding neutral transitions (t-row)
```

**Q2: Single-term ablation.**

For each of the 6 terms, remove it (set to zero) and recompute the full impedance equation:
```
Z_total(−term_k) = Σ_{j≠k} Z_j
P_ij = softmax(−Z_total)     [row-wise renormalisation]
S(−term_k) = Σ (ln P_ij − ln P_obs_ij)²
```

Report:
```
Term removed    S(−term)    ΔS = S(−term) − S_baseline    % degradation
```

Rank terms by impact: the term whose removal causes the largest S degradation is the most load-bearing.

For each ablation also report per-row scores:
```
S_u = Σ_j (ln P_uj − ln P_obs,uj)²
S_c = Σ_j (ln P_cj − ln P_obs,cj)²
S_t = Σ_j (ln P_tj − ln P_obs,tj)²
```
and ΔS_u, ΔS_c, ΔS_t relative to baseline per-row scores. A term can preserve row ordering while badly damaging the quantitative depth.

Also report per-row CKM signature preservation: does removing a term change the row ordering (d>s>b, s>d>b, b>s>d)?

**Q3: Pairwise ablation and interaction.**

For each of the C(6,2) = 15 term pairs, remove both and recompute S after full row-wise softmax renormalisation. Report:
```
Terms removed    S(−pair)    ΔS(−pair)
```

Compute the interaction term:
```
interaction(a,b) = [S(−a,−b) − S_baseline] − [S(−a) − S_baseline] − [S(−b) − S_baseline]
```

All S values are computed after row-wise softmax renormalisation. The interaction measures whether the pair is structurally coupled through the softmax, not a simple additive property of raw Z-space. Positive interaction = coupled (removing both is worse than expected). Negative = compensating (removing both is better than expected).

**Q4: Full subset ablation.**

Evaluate all non-empty subsets of the six Z-terms: 2⁶ − 1 = 63 models. For each subset report:
```
- active terms (bitmask or list)
- S
- ΔS from baseline
- term count
- row signatures preserved (yes/no)
```

Identify the smallest subset whose S remains within tolerance bands of baseline AND preserves all 3 row signatures (d>s>b, s>d>b, b>s>d):
```
S ≤ 0.05       with 3/3 signatures preserved
S ≤ 0.10       with 3/3 signatures preserved
S ≤ 2 × S_baseline  with 3/3 signatures preserved
```

A subset must preserve the qualitative row hierarchy to count as an effective reduced model, not just have a numerically tolerable S.

Q4 also reports each single-term-only model (Z_total = Z_k only, for each k). Removal tells you whether a term is necessary in context; single-term-only tells you whether a term has standalone explanatory power. This is useful for physical derivation priority in Phase 6.

##### Usage

```
python fermion_mass_primes_test44.py
```

##### Success criteria

| Gate | Criteria |
|------|----------|
| Q1 | 9×6 contribution matrix reported with signed and absolute shares, structurally zero terms identified |
| Q2 | All 6 single ablations reported with per-row scores, terms ranked by load-bearing impact |
| Q3 | All 15 pairwise ablations reported, interaction terms computed after softmax renormalisation |
| Q4 | All 63 subsets evaluated, minimum effective term count identified at each tolerance band |

Phase 5 success is measured, not assumed: the minimum effective term count, load-bearing ranking, and interaction structure emerge from the data.

##### What this establishes

The ablation profile tells us:
- Which terms to prioritise for derivation (Phase 6) — derive the load-bearing terms first
- Whether the six-term structure is minimal or contains redundancy
- Whether any terms are structurally coupled (suggesting a shared physical origin)
- Whether the CKM row signatures are robust to single-term removal

### fermion_mass_primes_test45.py
**Phase 6: Consecutive-prime band overlap diagnostic.**

The codebook primes decompose uniquely as sums of contiguous primes (the "bands"). Phases 1–5 treat packet energy as generic (E_p = p, no internal structure). Test 45 asks whether the band structure carries active CKM transition information — does the overlap between shed and gain bands affect transition rates?

```
Bands:
  5  = [2, 3]
  17 = [2, 3, 5, 7]
  53 = [5, 7, 11, 13, 17]
  59 = [17, 19, 23]
  97 = [29, 31, 37]
 127 = [3, 5, 7, 11, 13, 17, 19, 23, 29]
 211 = [67, 71, 73]
 271 = [7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43]
```

##### Three questions

1. What is the band overlap structure across the 9 CKM transitions?
2. Does adding a band-overlap Z-term to the 5-term baseline improve S?
3. Does band overlap subsume or complement the existing Z_ch (k-mode channel conductance)?

##### Method

**Baseline:** the post-Test44 canonical 5-term model (Z_bal + Z_ch + Z_R + Z_deg + Z_tp), S_baseline = 0.026666.

**Q1: Band overlap profile.**

For each CKM transition (shed packets → gain packets):
- Compute the union of constituent prime bands for all shed packets: B_shed
- Compute the union of constituent prime bands for all gain packets: B_gain
- Report both:
  - set overlap: unique constituent primes in B_shed ∩ B_gain
  - multiset overlap: constituent primes counted with multiplicity across packets (diagnostic — shows channel weight/intensity)

Candidate Z_band terms use set overlap by default. Multiset results are diagnostic unless explicitly promoted.

Additional metrics:

```
|B_shed|, |B_gain|
|shared| (set overlap count)
Jaccard = |shared| / |B_shed ∪ B_gain|
Σ(shared constituent primes)
gap = minimum prime-index distance between any q in B_shed and any q in B_gain
      (0 means direct overlap; prime-index distance, not numerical distance)
contiguous = B_shed ∪ B_gain contains every prime in the prime-index interval
             from min(B_shed ∪ B_gain) to max(...), with no missing intermediate primes
```

For pure-shedding transitions (t-row), B_gain = ∅ and all overlap metrics are zero. Set all Z_band candidates to 0 for pure-shedding transitions, including gap-based terms. Test45 is about shed/gain band overlap; pure-shedding has no gain band, so any nonzero penalty would silently introduce a new t-row term unrelated to the stated hypothesis.

**Q2: Band-overlap Z-term candidates.**

Define candidate Z_band terms. Higher band overlap should facilitate transitions (lower impedance), so Z_band should decrease with increasing overlap:

```
Candidate 1: Z_band = -α_s × |shared| / max(|B_shed|, |B_gain|)
  (negative = conductance relief, normalised by larger band)

Candidate 2: Z_band = -α_s × Jaccard
  (Jaccard overlap as conductance)

Candidate 3: Z_band = -α_s × Σ(shared) / P
  (shared constituent energy, normalised by P=840)

Candidate 4: Z_band = α_s × gap / p_max
  (band gap as impedance, normalised by p_max=271)

Candidate 5: Z_band = -α_s × ln(1 + |shared|)
  (logarithmic overlap conductance)
```

For each candidate:
- Add Z_band as a 6th term to the 5-term baseline
- Recompute P_ij = softmax(-Z_total) with row-wise renormalisation
- Report S, ΔS, per-row S, signature preservation

For any improving candidate, report whether improvement comes from:
- row ordering/signature correction
- per-row S reduction (which rows improve, which degrade)
- only uniform row sharpening after softmax

A candidate counts as physically interesting only if it improves at least one row-specific residual without damaging another row more than the net gain. This protects against false positives where a negative term merely rescales one row under softmax without capturing genuine sub-packet structure.

**Q3: Complementarity with Z_ch.**

For all five candidates, test:
- Complement: 5-term baseline + Z_band (6 terms, Z_ch retained)
- Replacement: replace Z_ch with Z_band (5 terms, Z_ch removed)

Report whether band overlap captures different physics from k-mode coupling. A band term might fail as an extra 6th term but still partially substitute for Z_ch, since Q3 explicitly asks whether band overlap subsumes or complements Z_ch.

Also test the best Q2 candidate against the full 63-subset ablation framework from test44: evaluate all subsets of {Z_bal, Z_ch, Z_R, Z_deg, Z_tp, Z_band} to find the new minimum effective model.

##### Usage

```
python fermion_mass_primes_test45.py
```

##### Success criteria

| Gate | Criteria |
|------|----------|
| Q1 | Band overlap profile reported for all 9 transitions |
| Q2 | All candidates tested, S improvement (or absence) measured |
| Q3 | Complementarity with Z_ch assessed, best model identified |

Phase 6 success is graduated:
- **Strong:** a band-overlap term improves on the 5-term baseline, S < 0.026666, with 3/3 signatures preserved — the consecutive-prime structure is physically active in CKM transitions
- **Moderate:** improvement exists but is smaller than Z_ch — band structure modulates but doesn't dominate
- **Weak/absent:** no improvement — the consecutive-prime property is purely a stability condition, not an active transition term

Any outcome is a valid result. A weak result is not a failure — it cleanly separates the stability role (packet selection) from the transition role (CKM rates).

##### What this establishes

If band overlap improves S: the internal sub-packet structure (consecutive-prime bands) carries physically active information about CKM transition rates. The "greased channels" from overlapping bands are a real mechanism, and packets are not generic energy blobs.

If band overlap does not improve S: the consecutive-prime property is a stability/formation condition (which packets can exist and how they transform under strain) but does not directly set transition rates. The CKM equation operates at the packet level; the sub-packet band structure operates at the formation level. Different layers.

### fermion_mass_primes_test46.py
**Phase 10A: Exhaustive codebook null.**

The codebook {5, 17, 53, 59, 97, 127, 211, 271} was identified empirically and validated through Phases 1–9. Test 46 asks the definitive selectivity question: is this the only 8-prime subset of the A279857 pool that satisfies all structural constraints discovered so far, or are there viable alternatives?

The test is exhaustive — every candidate is evaluated, not sampled. The starting universe is the full C(43,8) ≈ 145 million 8-subsets of the 43 A279857 primes. No pre-filtering; intermediate survivor sets (the 23-prime intersection, the 46 candidates, the codebook itself) must be rediscovered by the filter cascade.

This null is exhaustive inside the declared 9-mode packet ontology. It does not test arbitrary prime sets outside A279857, because those do not correspond to physical packet states p = Σ distinct k² on this topology.

##### Five questions

1. How many 8-subsets survive each topology-only filter (alone, cumulative, leave-one-out, pairwise with F1)?
2. How many survivors have viable band-transition structure for CKM transitions?
3. What is the codebook's natural-assignment and best-permutation impedance rank among all survivors?
4. Which filters exclude each nearest rival or near-miss?
5. Do negative controls (smooth-neighbour, band window length, contiguity, first-rep-only) carry independent weight?

##### Definitions

**Natural assignment:** ascending prime order mapped to the eight canonical v2.0 packet slots A..H, where the paper's Table 8 defines how slots A..H appear in each fermion state. For the codebook this is A=5, B=17, C=53, D=59, E=97, F=127, G=211, H=271.

##### Method

**Universe definition.**

```
Pool:       43 primes in A279857 for k = 1..9
            Computed by enumeration: all p prime where p = Σ distinct k²
            for some subset of {1², 2², ..., 9²}
Candidates: C(43,8) = 145,008,513
```

All valid k² representations listed per prime. Primes with multiple representations (53, 127, 211, others) have all masks enumerated. Representation ambiguity is handled at the co_sum legality stage: a set passes if at least one complete representation assignment is collision-free.

**Stage 1: Topology-only filters.**

Applied to all 145M candidates. No filter uses observed CKM magnitudes or row signatures.

Each filter reported four ways:
- **Alone**: how many pass this filter with no others
- **Cumulative**: running survivor count as filters are added
- **Leave-one-out**: survivors when this filter is removed from the full set
- **Pairwise with F1**: marginal effect after sum=840, showing genuine independence

*Set-level filters (no role/slot assignment assumed):*

```
  F1  sum = 840
      Closure condition: total packet energy = LCM(1..8), the 8-channel
      transfer cycle period
  F2  all 8 are consecutive-prime sums
      Formation condition: each packet is a coherent contiguous band on
      the prime substrate
  F3  unique consecutive-prime representation for all 8
      Band uniqueness: Z_band is well-defined only if each packet has
      exactly one band decomposition. Report also the weaker F2-only
      survivor count, because Z_band can in principle be made well-defined
      by a deterministic representation rule even when multiple bands
      exist. F3 is the stronger clean-ontology version
  F4  raw co_sum legality (≥1 collision-free rep assignment)
      Pauli consistency: the k² eigenmode occupancy patterns must allow
      distinct phase labels at shared modes (Phase 2)
  F5  full 9-mode coverage (all k=1..9 occupied in at least one prime)
      Eigenvalue completeness: the codebook must span the full cavity
      mode spectrum
  F6  p_max = 271
      Spectrum-bound diagnostic: the largest packet energy equals the
      observed cavity endpoint. Report both independent survivor count
      and marginal effect after F1, as it may be partially redundant
      with the closure constraint
```

*Role/assignment-aware diagnostics (tested separately, dependency on slot structure labelled):*

```
  R1  mod 6 low/high split (4 primes ≡ 5 mod 6, 4 primes ≡ 1 mod 6)
      Weak-orientation balance: the mod 6 split encodes T₃ polarity
      across the codebook (Phase 7). Set-level count, but physical
      interpretation assumes role assignment
  R2  sorted-adjacent mod 8 pair profile
      Tests whether the mod-8 residue pairs of adjacent primes (in
      ascending order) match the codebook's pair pattern
  F7  mod 8 balanced (2 per odd residue class {1,3,5,7})
      Residue symmetry: reported as diagnostic only — empirically
      observed property of the codebook, not derived from cavity physics
```

*Negative controls (expected: weak or no effect):*

```
  N1  smooth-neighbour (7-smooth p±1 for all 8)
      Previously tested, ~80% of primes pass — weak discriminator
  N2  band window length profile matches codebook
      Phase 9 negative: band window does not predict scale walking
  N3  A279857 mask contiguity
      Only works for 5 and 271 — not a codebook-wide selector
  N4  first-representation-only co_sum legality
      Control diagnostic against the full representation-aware test (F4).
      First-rep-only was historically misleading
  N5  band window length as scale-walking predictor
      Expected negative from Phase 9; included to confirm the failed
      running hypothesis is not reused as a selector
```

Output: cascade table with survivor count at each step, plus leave-one-out and pairwise-with-F1 tables.

**Stage 2: Band-transition structure.**

For every set surviving Stage 1:

- Compute shed/gain prime sets for all 9 CKM transitions (using the natural slot assignment and the best-permutation assignment)
- Band overlap count |B_shed ∩ B_gain| for each exchange transition
- Contiguous union property for each transition
- Z_band contribution (non-zero for exchange transitions, zero for pure-shedding)
- Number of legal representation assignments (strength metric)

**Stage 3: Impedance ranking.**

For every set surviving Stage 1 set-level filters (F1–F6), run full 8! = 40,320 slot assignment search under the canonical equation:

```
Z_total = Z_bal + Z_ch + Z_R + Z_deg + Z_tp + Z_band
```

with Z_a retired (Phase 5).

For each candidate set report:
- Natural assignment score S_nat and row signatures
- Best-permutation score S_best and row signatures
- Permutation lift = S_nat − S_best
- Rank among all surviving sets (natural and best)

Controls (run for top-ranked sets):
- 5-term score without Z_band
- Historical v2.0 6-term score with Z_a
- Individual Z-term ablation scores

**Stage 4: Combined result.**

Cascade table:

```
Stage                              Survivors
C(43,8) universe                   145,008,513
+ F1: sum = 840                    ...
+ F2: consecutive-prime sums       ...
+ F3: unique band representation   ...
+ F5: 9-mode coverage              ...
+ F6: p_max = 271                  ...
+ F4: co_sum legality              ...
= topology-only admissible basin

Diagnostics (not filters):
  R1: mod 6 split 4+4
  R2: sorted-adjacent mod 8 pair profile
  F7: mod 8 balance (empirical, not derived)

Impedance ranking:
+ row signatures preserved         ...
+ impedance S ≤ codebook S         ...
```

Nearest-rival diagnostic:
- For rivals surviving all topology-only filters: report impedance score, row signatures, natural/best ranks, and permutation lift.
- For near-miss rivals excluded before impedance ranking: report the first excluding filter and all filters failed.

Nearby non-codebook diagnostic: for each of the 15 non-codebook primes in the 23-intersection pool {41, 71, 83, 101, 109, 131, 139, 181, 197, 199, 223, 233, 251, 269, 281}, report:
- number of F1–F6 survivor sets containing it
- best natural impedance score among those sets
- best-permutation impedance score among those sets
- first filter where it becomes impossible to include

271 neighbourhood: why 271 and not 263, 269, 277, 281, 283.

**CKM-input audit.**

- Set-level topology filters F1–F6 do not use observed CKM magnitudes, row signatures, or slot assignments
- Role-aware diagnostics R1–R2 do not use observed CKM magnitudes or row signatures, but their dependence on slot/role structure is reported separately
- Z_band functional form is fixed from Test45 before Test46
- Impedance scoring is performed only after topology-only filtering is complete
- Results reported both before and after impedance scoring

##### Usage

```
python fermion_mass_primes_test46.py
```

Progress reported at 1M-candidate intervals, flushed to stdout (`print(..., flush=True)`) so output is visible in real-time on console. The script prints the generated 43-prime pool, representation list, filter order, codebook index, and configuration hash before running the cascade. Smoke test on the 23-intersection pool (C(23,8) = 490,314) runs first as a gate check before the full 145M enumeration. Full run on local VM; Lambda GPU instance if runtime exceeds practical limits.

##### Success criteria

| Gate | Criteria |
|------|----------|
| A | Cascade table complete; cheap filters (F1–F3, F5, F6) reported alone/cumulative/leave-one-out/pairwise; F4 reported cascade-scoped and F1-scoped (full-alone over all C(43,8) not computed) |
| B | Codebook rediscovered as survivor of all topology-only filters |
| C | Codebook impedance rank reported (natural and best-permutation) |
| D | All nearest rivals and near-misses identified with exclusion reasons |
| E | Negative controls reported, including controls that do not help the codebook |
| F | If codebook is not unique, all surviving rivals and their distinguishing properties listed |
| G | Full 8! assignment search completed for every F1–F6 survivor; no sampling, shortlist, heuristic cutoff, adaptive cutoff, or score-based pruning used |

##### What this establishes

If the codebook is the unique survivor of all set-level structural filters F1–F6, with R1–R2 reported separately as role-aware diagnostics: the 8-prime set is a mathematical fixed point of the constraint basin — no alternative exists within the declared 9-mode packet ontology.

If a small number of rivals survive: the constraints narrow from 145M to a countable few, and the impedance score provides the final selection. The codebook's rank among survivors quantifies its selectivity.

If many rivals survive with comparable impedance scores: the current filter set is insufficient to uniquely select the codebook, and additional structural properties must exist. This would be the most informative outcome for guiding further work.

Any outcome is a valid result. The test establishes selectivity only within the declared 9-mode packet ontology. Alternative ontologies would require their own independent null tests.

### fermion_mass_primes_test47.py
**Phase 10B: Neutrino geometry null.**

The Phase 8 PMNS angle formulas are structural projection hypotheses: zero-continuous-parameter and numerically strong, but the specific projection rules were identified, not derived. Test 47 quantifies the look-elsewhere risk by exhaustively enumerating geometric projection rules on the fixed triangle (edges 117, 325, 398) and counting how many produce PMNS-quality angle triplets and mass-splitting ratios.

The triangle is fixed by the codebook, which Phase 10A showed is the unique rank-1 solution from 145 million candidates. Phase 10B does not re-test the codebook; it tests how special the geometric rules are on the unique triangle.

##### Four questions

1. How many triangle centres produce a mass-splitting ratio within 2σ of observed, and under which weight/distance-law combinations?
2. How many distinct angle triplets from any centre match all three PMNS angles, and what is the Phase 8 rule's χ² rank?
3. Does the incenter remain the unique best centre under the primary null (same centre, canonical labelling, low complexity)?
4. Is the w = ln(p), 1/d² combination uniquely selected, or do nearby alternatives also match?

##### Observed inputs (frozen before running)

```
Source: NuFIT 5.3 (2024), normal ordering

θ₁₂ = 33.44° ± 0.75°
θ₂₃ = 49.0°  ± 1.4°
θ₁₃ = 8.54°  ± 0.15°

Δm²₃₂/Δm²₂₁ = 32.6 ± 1.2 (1σ)
```

##### Definitions

**Triangle:** edges UL = 117, UD = 325, DL = 398. Perimeter = 840. Vertex angles from cosine rule: Up = 120.93°, Lepton = 44.46°, Down = 14.61°.

**Phase 8 rules (Tier 1):**

```
θ₂₃ = sep(Up-Lepton from incenter) / 2 = 48.65°
θ₁₂ = sep(Up-Down from incenter) / 2 − Lepton half-angle = 33.88°
θ₁₃ = Down half-angle + 3/840 × 360° = 8.59°
```

**Complexity:** number of primitive geometric quantities + number of arithmetic operations + QCD correction flag (0 or 1). Phase 8 complexities: θ₂₃ = 2 (one separation, one halving), θ₁₂ = 3 (one separation, one halving, one subtraction), θ₁₃ = 3 (one half-angle, one addition, one correction). Total complexity = 8. Derived pool quantities retain their construction complexity: |q_i − q_j| has complexity(q_i) + complexity(q_j) + 1, not complexity 1 merely because it appears in the pool.

##### Method

**Step 1: Enumerate triangle centres.**

Ten centres (interior/exterior status reported for this obtuse triangle):

```
  Incenter             interior
  Centroid             interior
  Circumcenter         interior or exterior (report which)
  Orthocenter          exterior (obtuse triangle)
  Nine-point centre    interior (always)
  Fermat point         boundary (for triangles with angle ≥120°, the Fermat
                       point coincides with the obtuse vertex; Up = 120.93°,
                       so Fermat point = Up vertex)
  Symmedian point      interior
  Excentre opp. Up     exterior
  Excentre opp. Down   exterior
  Excentre opp. Lepton exterior
```

For each centre, compute coordinates, spoke distances to all three vertices, and angular separations between all three spoke pairs.

Primary physical null uses interior centres only. Boundary centres (e.g. the obtuse-triangle Fermat point at the Up vertex) and exterior centres (orthocenter, excentres) are excluded from the primary null and included only in the diagnostic all-centres null.

**Step 2: Build the geometric quantity pool per centre.**

For each centre, enumerate all quantities in [0°, 90°]:

```
  3 vertex angles
  3 half-angles
  3 angular separations from centre
  3 half-separations
  6 pairwise differences |q_i − q_j|
  6 pairwise sums (if in range)
  3 third-angles: 180° − q_i − q_j
  Correction variants: half-angle + n/840 × 360° for n = 1..5
```

Cluster quantities within ε = 0.01° before triplet enumeration to prevent near-duplicate over-counting. Report both raw and clustered pool sizes. Triplets use three distinct clustered quantities; reuse of the same quantity for two PMNS angles is excluded from the primary null and reported separately as a reuse-control diagnostic.

**Step 3: Edge relabelling.**

Run two modes:
- Canonical labelling: UL = 117, UD = 325, DL = 398
- All 3! = 6 relabellings of {117, 325, 398} onto {UL, UD, DL}

**Step 4: Score all ordered triplets.**

For each (centre, labelling, triplet assignment to θ₁₂, θ₂₃, θ₁₃):

```
  residual_i = |predicted_i − observed_i| / σ_i
  χ² = Σ residual_i²
```

Report tiered results:

```
  Tier        Centre              Labelling      Complexity
  Primary:    same interior ctr   canonical      ≤ Phase 8
  Secondary:  same interior ctr   all 6          ≤ Phase 8
  Tertiary-I: all interior ctrs   all 6          any
  Tertiary-A: all ctrs incl.      all 6          any
              boundary/exterior
```

For each tier report:
- Total triplets tested
- Number matching at 1σ (all residuals < 1), 2σ, 3σ
- Phase 8 exact package rank: incenter + canonical labelling + the three Phase 8 formulas
- Best Phase 8-family rank: same formula family allowing centre/labelling substitutions
- Best χ² overall and its identity
- Number of triplets with χ² ≤ χ²(Phase 8)

**Step 5: Weight function × distance law scan.**

For each centre, test all combinations:

Weight functions:

```
  w = p, ln(p), sqrt(p), 1/p, p^(1/3), Σ(1/k²), ln(p)²
```

Distance laws:

```
  W/d, W/d², W/d³, W only, 1/d² only
```

For each (centre, weight, distance) triplet:
- Compute vertex weights W_V and coupling strengths c_V
- Test all 3! = 6 assignments of vertices to {ν₁, ν₂, ν₃}
- Compute mass-splitting ratio R for each assignment
- Report residual from observed R = 32.6
- Primary pass requires normal ordering (positive Δm²₂₁ and m₃ > m₂ > m₁). Inverted-ordering matches are reported separately as diagnostic failures for the Phase 8 normal-ordering prediction

Report:
- Phase 8 combination (incenter, ln(p), 1/d²) rank
- Number of combinations with R within 1σ, 2σ
- Best R and its identity
- Whether normal ordering is preserved

**Step 6: Joint package score.**

For packages that specify both a centre/projection triplet and a weight/distance/assignment rule, compute:

```
  χ²_joint = χ²_angles + χ²_mass_ratio
```

Joint packages must share the same centre, the same edge/vertex labelling, and the same vertex-to-neutrino assignment where applicable. The Phase 8 joint package is: incenter + canonical labelling + Phase 8 angle rules + ln(p)/d² + Up→ν₃, Lepton→ν₂, Down→ν₁.

Report Phase 8 joint rank under the primary tier. This tests the coherence of the full Phase 8 claim: one centre (incenter) with one projection rule set and one induction law (ln(p)/d²), not a mix of different centres for angles and masses.

**Headline statistic:** the clustered, same-centre, canonical-labelling, complexity-bounded count at the primary tier. Raw full-grammar triplet counts from the tertiary search are reported as look-elsewhere diagnostics only.

##### Usage

```
python fermion_mass_primes_test47.py
```

Console output with flush. No parallelisation needed.

##### Success criteria

| Gate | Criteria |
|------|----------|
| A | All centres computed with coordinates, spoke distances, angular separations, interior/exterior status |
| B | Full triplet search completed at all three tiers, χ² ranking reported |
| C | Weight × distance × assignment scan completed for all centres |
| D | Phase 8 rule's rank and selectivity quantified at each tier |

##### What this establishes

If the Phase 8 package is the only same-centre, canonical-labelling, equal-or-lower-complexity triplet matching all PMNS angles at 1σ, then the projection rules are strongly selected within the declared geometric grammar. They are still not derived; Test 47 quantifies selectivity, not derivation status.

If multiple formula sets match at the primary tier: the look-elsewhere factor is measured and reported. The PMNS formulas remain structural hypotheses with a quantified look-elsewhere risk.

If a non-incenter centre produces a tighter match: the Phase 8 centre choice requires revision.

Any outcome is a valid result.


