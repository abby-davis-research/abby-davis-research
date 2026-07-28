"""Tests for the LatticeNull proximity protocol.

Run with:  python3 -m pytest proximity_protocol/test_lattice_null.py
or plainly: python3 proximity_protocol/test_lattice_null.py
"""

import math

from lattice_null import (
    PHI,
    Claim,
    benjamini_hochberg,
    evaluate_table,
    guaranteed_window,
    max_t,
    metallic_mean,
    one_in,
    p_chance,
)

TOL = 1e-9


def approx(a, b, tol=1e-4):
    return abs(a - b) <= tol


def test_metallic_means():
    assert approx(metallic_mean(1), PHI, TOL)          # golden
    assert approx(metallic_mean(2), 1 + math.sqrt(2), TOL)  # silver


def test_saturation_at_max_t():
    # At t = (b-1)/(b+1) the window fills a whole cell -> P = 1.
    for b in (PHI, math.sqrt(2), math.e, math.pi):
        assert approx(p_chance(max_t(b), b), 1.0, TOL)
        assert p_chance(max_t(b) + 0.01, b) == 1.0
    assert p_chance(0.0, PHI) == 0.0


def test_phi_is_2_38x_more_chance_prone_than_pi():
    # The protocol's headline: φ at t=0.05 is 2.38x more chance-prone than π.
    ratio = p_chance(0.05, PHI) / p_chance(0.05, math.pi)
    assert approx(ratio, 2.38, tol=0.01)


def test_phi_one_in_6_3_at_t_038():
    # Corrected protocol figure: φ is 1-in-6.3 chance-prone at t=0.038.
    # (The summary block mis-stated this as t=0.05, which is actually 1-in-4.81.)
    assert approx(one_in(0.038, PHI), 6.3, tol=0.05)
    assert approx(one_in(0.05, PHI), 4.81, tol=0.01)


def test_guaranteed_window_round_trips():
    # The natural log-symmetric band's chance mass must round-trip to the
    # threshold it was built from, and be asymmetric in linear space
    # (upper reach larger in magnitude than lower).
    b, th = PHI, 0.4493
    lo, hi = guaranteed_window(b, th)
    assert hi > abs(lo)                       # asymmetric in linear space
    log_width = math.log(1 + hi) - math.log(1 + lo)
    assert approx(log_width / math.log(b), th, tol=1e-9)
    # This band is NOT the protocol's exact (-0.1066, +0.1014); that comes from
    # a construction in the source .tex not specified in the summary block.
    assert not approx(lo, -0.1066, tol=5e-4)


def test_free_c_is_a_tautology():
    fixed = Claim("fixed", "phi", delta_obs=0.01, c_fixed=True)
    free = Claim("free", "phi", delta_obs=0.01, c_fixed=False)
    assert free.p_value() == 1.0            # delta can be driven to 0 -> no evidence
    assert fixed.p_value() < 1.0


def test_benjamini_hochberg_basic():
    # All-null p-values (uniform-ish) should reject nothing at q=0.05.
    pvals = [0.9, 0.8, 0.5, 0.4, 0.2]
    assert benjamini_hochberg(pvals, q=0.05) == [False] * 5
    # A single very small p among noise is rejected.
    pvals = [0.001, 0.8, 0.5, 0.4, 0.2]
    assert benjamini_hochberg(pvals, q=0.05)[0] is True


def test_evidence_table_expected_null_positives():
    # A table of loose φ-proximity claims: many "pass" naively, none survive BH.
    claims = [Claim(f"c{i}", "phi", delta_obs=0.05) for i in range(18)]
    rep = evaluate_table(claims, q=0.05)
    # E[NP] = 18 * P(0.05, phi) ≈ 18 * 0.208 ≈ 3.7 chance matches expected.
    assert approx(rep.expected_null_positives, 18 * p_chance(0.05, PHI), tol=1e-6)
    assert rep.significant_bh == []  # nothing significant after BH — the point


def _run_all():
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn()
        print(f"ok  {fn.__name__}")
    print(f"\n{len(fns)} tests passed")


if __name__ == "__main__":
    _run_all()
