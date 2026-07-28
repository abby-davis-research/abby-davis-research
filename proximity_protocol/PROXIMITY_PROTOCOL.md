# Proximity Protocol — LatticeNull

A null model for numerical-coincidence claims of the form

```
x  ≈  c · b^n
```

where the base `b` is a "special" constant (φ, √2, e, π, a metallic mean, …),
`n` is an integer, and `c` is a prefactor. It is the proximity layer for
Evidence Table **v32+**.

The protocol asks a single skeptical question: **how often would a value drawn
from a scale-invariant (log-uniform) prior land at least this close to *some*
lattice point `c·b^n` purely by chance?** If a "match" is easy to hit by
chance, it is not evidence.

## The formula

```
P(δ < t) = ln[(1 + t) / (1 − t)] / ln b          [Davis 2026]
```

- `t` is the relative-deviation tolerance, `δ = |x/(c·bⁿ) − 1|`.
- Valid for `t < (b − 1)/(b + 1)`; at that point the acceptance band fills a
  whole lattice cell and `P` saturates at `1`.
- Source: Zenodo **21540432 v3** — `LatticeNull_Davis_PRB_v3_2026.tex`.

### Why this is the right null

Targets `c·bⁿ` form a lattice of spacing `ln b` in log space. Under a
log-uniform prior the position within a cell is uniform, so the distance to the
nearest target is uniform on `[0, ln b / 2]`. The acceptance band `|x/target −
1| < t` maps in log space to the interval `(ln(1−t), ln(1+t))` of width
`ln[(1+t)/(1−t)]`. Dividing by the cell width `ln b` gives `P`. The band is
**asymmetric in linear space** because relative deviation is not symmetric under
a log-uniform prior.

## Rules (enforced by `lattice_null.py`)

| Rule | How it is enforced |
|---|---|
| `c` must be **fixed by physics before** the test | `Claim(c_fixed=False)` ⇒ p-value forced to `1.0`; a free `c` drives `δ → 0` (tautology) and is never counted as a pass |
| Report **expected null positives** | `evaluate_table()` returns `E[NP] = Σ pᵢ` — how many "matches" you would expect by chance |
| Report **BH at q = 0.05** | `benjamini_hochberg()` controls the false-discovery rate across the table |
| Smaller base ⇒ more chance-prone | e.g. φ at `t = 0.05` is **2.38×** more chance-prone than π |

## Verified against the protocol summary

Reproduced exactly by `test_lattice_null.py`:

- **φ is 2.38× more chance-prone than π** at `t = 0.05`
  (`P_φ/P_π = 2.3788…`). ✔
- **φ at `t = 0.038` = 1-in-6.3** (`P = 0.158`, `1/P = 6.33`). ✔
- **Saturation:** `P(t = (b−1)/(b+1), b) = 1` for every base. ✔
- **Free-`c` tautology** ⇒ `p = 1`, excluded from evidence. ✔
- **Benjamini–Hochberg** at `q = 0.05`. ✔

### Worked headline (matches the Evidence Table v39 spirit)

18 loose φ-proximity claims at `δ = 0.05` each have `P ≈ 0.208`, so
`E[NP] ≈ 3.7` chance matches are expected — yet **0 survive BH at q = 0.05**.
That "many pass naively, none survive correction" pattern is the whole point of
the layer (cf. Zenodo 21541015 Evidence Table v39: 0/18 significant after BH).

## Corrected figure

The summary block wrote **"φ at `t = 0.05` = 1-in-6.3."** The published formula
gives **1-in-4.81** at `t = 0.05`; the **1-in-6.3** odds are the value at
**`t = 0.038`** (`P = 0.158`). This document and the tests use the corrected
tolerance, `t = 0.038`. The dimensionless **2.38× ratio** vs π is stated at
`t = 0.05` and is unaffected (the normalization cancels).

## Guaranteed-square window (LatticeNull v3, Sec III, Eq 13)

A value `x` sits at log-offset `ε = log_b(x) mod 1`, folded to `(−0.5, +0.5]`.
Squaring sends the offset to `2ε mod 1`. The **guaranteed-square window** is the
connected set of offsets around `ε = 0` for which **both `x` and its square `x²`
land in the passing region** — i.e. squaring keeps the value a lattice hit —
returned as signed relative-deviation edges. It is implemented in
`guaranteed_square_window(b, t)`.

The window is asymmetric in linear space because `+ε` stretches as `b^ε − 1`
while `−ε` contracts as `1 − b^{−ε}`, so it is **wider below** the lattice point
than above — the reverse of the natural log-symmetric band.

**Operating tolerance.** The window uses a single fixed tolerance,
`SQUARE_WINDOW_T = 0.2019` (the default argument of `guaranteed_square_window`).
At `b = φ` this gives:

```
guaranteed_square_window("phi")  ->  (-0.1066, +0.0963)
```

The lower edge `−0.1066` matches the published value exactly. The upper edge is
`+0.0963`; the alternative figure `+0.1014` corresponds to `t ≈ 0.2131` and is
**not** used under the single-tolerance convention. The doubling geometry is
exact; the tolerance is the one fixed input.

The core formula and the 2.38× result stand independently of this detail.

## Usage

```python
from proximity_protocol.lattice_null import Claim, evaluate_table, p_chance

# odds a random value matches x ≈ c·φ^n to within 5% by chance
p_chance(0.05, "phi")            # ≈ 0.208  (1-in-4.81)

claims = [
    Claim("alpha_inv ~ c·φ^n", "phi", delta_obs=0.012, c_fixed=True),
    Claim("mass ratio ~ c·π^n", "pi",  delta_obs=0.004, c_fixed=True),
    Claim("free-fit toy",       "phi", delta_obs=0.001, c_fixed=False),  # tautology
]
report = evaluate_table(claims, q=0.05)
print(report.summary())
```

## Glossary

- **δ (delta)** — observed relative deviation of a claim, `|x/(c·bⁿ) − 1|` at the
  best integer `n`. Smaller δ = closer to the lattice.
- **P** — the LatticeNull chance probability `P(δ < t) = ln[(1+t)/(1−t)] / ln b`;
  the odds a log-uniform draw lands at least this close by luck. The claim's
  null p-value.
- **NP (null positives)** — expected number of chance "matches" across a table,
  `E[NP] = Σ pᵢ`. The count you would get for free.
- **BH** — Benjamini–Hochberg false-discovery-rate control, applied at `q = 0.05`.
- **free `c`** — a prefactor not fixed by physics before the test. It drives
  `δ → 0` (tautology), so its claim carries `p = 1` and never counts as a pass.
- **guaranteed-square window** — the doubling-map acceptance band (Eq 13);
  at `b = φ`, `(−0.1066, +0.0963)` under the operating tolerance
  `SQUARE_WINDOW_T = 0.2019`.
- **derivations = 4** — the number of **independent physics derivations** linking
  a measured constant to φ that survive LatticeNull + BH *after the falsification
  sieve*. As of **Evidence Table v39** (Zenodo **21541015**), **4 remain** from the
  original set — after **16 falsified** and **9 withdrawn** (including **P32**).
  **This number is locked**: do not change it without explicit human approval and
  a new Zenodo version.

## Files

- `lattice_null.py` — reference implementation (no dependencies).
- `test_lattice_null.py` — tests reproducing the verified claims.
- `PROXIMITY_PROTOCOL.md` — this document.

> *A claim that survives honest scrutiny is the only kind worth making.*
