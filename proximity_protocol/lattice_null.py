"""LatticeNull proximity protocol — reference implementation.

A null model for claims of the form

    x  ≈  c · b^n

where the base ``b`` is a "special" constant (φ, √2, e, π, a metallic mean, …)
and ``c`` is a prefactor that is supposed to be *fixed by physics before the
test*. The protocol answers one question:

    How often would a value drawn from a scale-invariant (log-uniform) prior
    land at least this close to *some* lattice point c·b^n purely by chance?

That chance probability is

    P(δ < t) = ln[(1 + t) / (1 - t)] / ln b            (Davis 2026)

valid for t < (b - 1)/(b + 1); for larger t the acceptance window covers a
whole lattice cell and P = 1.

Derivation
----------
The targets c·b^n sit on a lattice of spacing ``ln b`` in log space. Under a
log-uniform prior the fractional position within a cell is uniform, so the
distance to the nearest target is uniform on [0, ln b / 2]. A relative-deviation
acceptance band |x/target - 1| < t maps, in log space, to the interval
(ln(1 - t), ln(1 + t)) whose width is ln[(1+t)/(1-t)]. Dividing by the cell
width ``ln b`` gives P. The band is asymmetric in linear space (that is why the
guaranteed-acceptance window for φ is (-0.1066, +0.1014) rather than ±0.107):
relative deviation is not symmetric under the log-uniform prior.

This module intentionally has no third-party dependencies.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import Iterable, Sequence

# --- special bases ---------------------------------------------------------

PHI = (1.0 + math.sqrt(5.0)) / 2.0

# Single operating tolerance for the guaranteed-square window (LatticeNull v3,
# Sec III, Eq 13). At b=φ this yields the window (-0.1066, +0.0963); the lower
# edge matches the published -0.1066 exactly.
SQUARE_WINDOW_T = 0.2019


def metallic_mean(m: int) -> float:
    """The m-th metallic mean, (m + sqrt(m**2 + 4)) / 2.

    m=1 -> golden ratio φ, m=2 -> silver ratio, m=3 -> bronze ratio, ...
    """
    if m < 1:
        raise ValueError("metallic mean index m must be >= 1")
    return (m + math.sqrt(m * m + 4.0)) / 2.0


SPECIAL_BASES = {
    "phi": PHI,
    "sqrt2": math.sqrt(2.0),
    "e": math.e,
    "pi": math.pi,
    "silver": metallic_mean(2),
    "bronze": metallic_mean(3),
}


# --- core null model -------------------------------------------------------

def resolve_base(b) -> float:
    """Accept either a numeric base or a name in SPECIAL_BASES (e.g. "phi")."""
    if isinstance(b, str):
        key = b.lower()
        if key not in SPECIAL_BASES:
            raise ValueError(f"unknown named base {b!r}")
        return SPECIAL_BASES[key]
    if b <= 1.0:
        raise ValueError("base b must be > 1")
    return float(b)


def max_t(b) -> float:
    """Largest tolerance for which P < 1; the window fills a full cell here."""
    b = resolve_base(b)
    return (b - 1.0) / (b + 1.0)


def p_chance(t: float, b) -> float:
    """Probability a log-uniform draw lands within relative tolerance ``t`` of
    a lattice point c·b^n, i.e. the chance a claim "matches" for free.

    Returns a value in [0, 1]. For t >= (b-1)/(b+1) the window covers a whole
    cell and the probability saturates at 1.
    """
    b = resolve_base(b)
    if t <= 0.0:
        return 0.0
    if t >= max_t(b):
        return 1.0
    return math.log((1.0 + t) / (1.0 - t)) / math.log(b)


def one_in(t: float, b) -> float:
    """Odds expressed as 1-in-N (N = 1 / P). ``inf`` when P == 0."""
    p = p_chance(t, b)
    return math.inf if p == 0.0 else 1.0 / p


def guaranteed_window(b, p_threshold: float) -> tuple[float, float]:
    """Natural relative-deviation band (delta_low, delta_high) whose chance
    mass under the null equals ``p_threshold``.

    The band is centred on the lattice point in *log* space, so it is
    asymmetric in linear (relative-deviation) space: the upper reach exceeds
    the lower one in magnitude. Its chance mass round-trips exactly, i.e.
    ``p_chance(delta_high, b) + p_chance(-delta_low... )`` reduces to
    ``p_threshold`` because log-width / ln b == p_threshold.

    NOTE: this is the *natural* log-symmetric band (wider above than below in
    linear space). It is NOT the "guaranteed-square" window — that one comes
    from the doubling map and is computed by ``guaranteed_square_window``.
    """
    if not 0.0 < p_threshold < 1.0:
        raise ValueError("p_threshold must be in (0, 1)")
    b = resolve_base(b)
    half_log = p_threshold * math.log(b) / 2.0  # symmetric half-width in log space
    delta_high = math.exp(half_log) - 1.0
    delta_low = math.exp(-half_log) - 1.0
    return (delta_low, delta_high)


def _fold(x: float) -> float:
    """Fold a log-offset to the half-open cell (-0.5, 0.5] about the nearest
    lattice point."""
    return x - round(x)


def _passes(offset: float, b: float, t: float) -> bool:
    """Whether a value at log-offset ``offset`` (in units of ln b) lies within
    relative tolerance ``t`` of its nearest lattice point, using the asymmetric
    δ convention  δ = |b^offset - 1|  (Davis 2026, LatticeNull v3 Eq 13)."""
    return abs(b ** offset - 1.0) < t


def guaranteed_square_window(
    b, t: float = SQUARE_WINDOW_T, *, _iters: int = 60
) -> tuple[float, float]:
    """The guaranteed-square window (LatticeNull v3, Sec III, Eq 13).

    A value x sits at log-offset ε = log_b(x) mod 1 (folded to (-0.5, 0.5]).
    Squaring sends the offset to ``2ε mod 1``. The window is the connected set
    of offsets around ε = 0 for which BOTH x and its square x² fall within the
    passing region — i.e. squaring keeps the value a lattice "hit". It is
    returned as signed relative-deviation edges (delta_low, delta_high).

    The window is asymmetric in linear space because +ε stretches as ``b^ε − 1``
    while −ε contracts as ``1 − b^{−ε}``.

    The operating tolerance is the single value ``SQUARE_WINDOW_T = 0.2019``
    (the default). At b=φ this gives the window ``(-0.1066, +0.0963)``: the
    lower edge matches the published -0.1066 exactly. (The alternative upper
    figure +0.1014 corresponds to t≈0.2131 and is not used under the
    single-tolerance convention.)
    """
    b = resolve_base(b)
    if not 0.0 < t < 1.0:
        raise ValueError("t must be in (0, 1)")

    def ok(eps: float) -> bool:
        return _passes(eps, b, t) and _passes(_fold(2.0 * eps), b, t)

    if not ok(0.0):
        return (0.0, 0.0)

    def edge(sign: int) -> float:
        lo, hi = 0.0, 0.5
        for _ in range(_iters):
            mid = 0.5 * (lo + hi)
            if ok(sign * mid):
                lo = mid
            else:
                hi = mid
        return sign * lo

    e_hi, e_lo = edge(+1), edge(-1)
    return (b ** e_lo - 1.0, b ** e_hi - 1.0)


# --- claims and evidence tables -------------------------------------------

@dataclass
class Claim:
    """A single proximity claim x ≈ c·b^n.

    Attributes
    ----------
    name:
        Human-readable label.
    base:
        The special constant b (> 1). May be given as a name in SPECIAL_BASES.
    delta_obs:
        Observed relative deviation |x/(c·b^n) - 1| at the best integer n.
    c_fixed:
        Whether the prefactor c was fixed by physics BEFORE the test. If False,
        c is a free knob that can drive delta to 0, so the "match" is a
        tautology and carries no evidential weight (p-value forced to 1.0).
    """

    name: str
    base: float
    delta_obs: float
    c_fixed: bool = True

    def __post_init__(self) -> None:
        self.base = resolve_base(self.base)
        if self.delta_obs < 0.0:
            raise ValueError("delta_obs must be >= 0")

    def p_value(self) -> float:
        """Null p-value: chance a random value lands at least this close.

        A free prefactor makes the claim a tautology (delta can be driven to 0),
        so it is reported as p = 1.0 and never counted as a pass.
        """
        if not self.c_fixed:
            return 1.0
        return p_chance(self.delta_obs, self.base)


def benjamini_hochberg(pvalues: Sequence[float], q: float = 0.05) -> list[bool]:
    """Benjamini-Hochberg FDR procedure.

    Returns a boolean list (aligned with ``pvalues``) marking which hypotheses
    are rejected (i.e. significant) at false-discovery rate ``q``.
    """
    m = len(pvalues)
    rejected = [False] * m
    if m == 0:
        return rejected
    order = sorted(range(m), key=lambda i: pvalues[i])
    k_max = -1
    for rank, idx in enumerate(order, start=1):
        if pvalues[idx] <= q * rank / m:
            k_max = rank
    for rank, idx in enumerate(order, start=1):
        if rank <= k_max:
            rejected[idx] = True
    return rejected


@dataclass
class TableReport:
    n_claims: int
    n_tautologies: int
    expected_null_positives: float  # E[chance matches] = sum of p-values
    significant_bh: list[str]
    q: float
    per_claim: list[tuple[str, float, bool]] = field(default_factory=list)

    def summary(self) -> str:
        lines = [
            f"claims tested        : {self.n_claims}",
            f"free-c tautologies   : {self.n_tautologies} (p=1, excluded from evidence)",
            f"expected null matches: {self.expected_null_positives:.2f}  (E[NP] = sum of p-values)",
            f"significant @ BH q={self.q:g}: {len(self.significant_bh)}"
            f" / {self.n_claims}"
            + (f"  -> {', '.join(self.significant_bh)}" if self.significant_bh else ""),
        ]
        return "\n".join(lines)


def evaluate_table(claims: Iterable[Claim], q: float = 0.05) -> TableReport:
    """Run the full protocol over an evidence table.

    Reports, per the protocol:
      * expected null positives  E[NP] = sum of per-claim p-values, i.e. how
        many "matches" you would expect purely by chance, and
      * the Benjamini-Hochberg significant set at FDR ``q``.

    Free-c tautologies are carried with p=1 so they inflate neither the
    expected-chance count nor the significant set.
    """
    claims = list(claims)
    pvals = [c.p_value() for c in claims]
    rejected = benjamini_hochberg(pvals, q=q)
    return TableReport(
        n_claims=len(claims),
        n_tautologies=sum(1 for c in claims if not c.c_fixed),
        expected_null_positives=sum(pvals),
        significant_bh=[c.name for c, r in zip(claims, rejected) if r],
        q=q,
        per_claim=[(c.name, p, r) for c, p, r in zip(claims, pvals, rejected)],
    )
