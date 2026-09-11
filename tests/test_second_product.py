"""The same chain on a second product, and the registered prediction for a third rate.

Lycopene is a real second product -- a lycopene strain is a beta-carotene strain with the
cyclase deleted -- but it is a WEAK generality test, because it shares the pathway, the
dataset and the fitted parameter. It shows the solver is not hard-wired to one product.
It does not show the architecture transfers to an unrelated pathway, and nothing here
should be read as showing that.
"""

import pandas as pd
import pytest

from ystwin import paths
from ystwin.kinetic.carotenoid import ELIZONDO2025, CyclaseCalibration, solve_branch_from_flux

CALIBRATION = paths.data_dir() / "carotenoid" / "elizondo2025_steady_states.tsv"


@pytest.fixture(scope="module")
def states():
    if not CALIBRATION.exists():
        pytest.skip(f"{CALIBRATION} not present")
    d = pd.read_csv(CALIBRATION, sep="\t")
    d["flux"] = d.q_lycopene + d.q_betacarotene
    return d


class TestLycopeneAsTheProduct:
    def test_it_reproduces_the_measured_pools(self, states):
        errs = []
        for r in states.itertuples():
            got = solve_branch_from_flux(r.flux, r.mu_per_h)
            want = r.q_lycopene / r.mu_per_h
            errs.append(abs(got.lycopene_content - want) / want)

        assert sorted(errs)[len(errs) // 2] < 0.10

    def test_deleting_the_cyclase_sends_every_carbon_to_lycopene(self):
        """The limit that must hold exactly: with no cyclase the pool is flux over
        dilution, since washout is the only sink. A solver that got this wrong would be
        wrong about the mass balance itself."""
        no_cyclase = CyclaseCalibration(
            **{**ELIZONDO2025.__dict__, "capacity_mmol_per_gdcw": 1e-12})
        flux, mu = 6.0e-4, 0.15

        got = solve_branch_from_flux(flux, mu, no_cyclase)

        assert got.lycopene_content == pytest.approx(flux / mu, rel=1e-6)
        assert got.beta_carotene_rate == pytest.approx(0.0, abs=1e-9)

    def test_a_stronger_cyclase_moves_carbon_off_lycopene(self, states):
        weak = CyclaseCalibration(**{**ELIZONDO2025.__dict__,
                                     "capacity_mmol_per_gdcw": ELIZONDO2025.capacity_ci95[0]})
        strong = CyclaseCalibration(**{**ELIZONDO2025.__dict__,
                                       "capacity_mmol_per_gdcw": ELIZONDO2025.capacity_ci95[1]})

        assert (solve_branch_from_flux(6e-4, 0.15, strong).lycopene_content
                < solve_branch_from_flux(6e-4, 0.15, weak).lycopene_content)

    def test_carbon_is_conserved_through_the_branch(self, states):
        """Whatever the desaturase makes leaves as one pool or the other."""
        for r in states.itertuples():
            got = solve_branch_from_flux(r.flux, r.mu_per_h)
            out = got.growth_rate * got.lycopene_content + got.beta_carotene_rate

            assert out == pytest.approx(r.flux, rel=1e-6)


class TestTheThirdRateIsWhatIsMissing:
    def test_the_calibration_has_exactly_two_dilution_rates(self, states):
        """The structural limit. The fitted law is vmax proportional to mu, and two rates
        are two points: no functional form on the growth-rate axis can be refuted by
        them. A third rate is the first observation that could distinguish this law from
        a saturating one, or a constant one."""
        assert states.mu_per_h.nunique() == 2

    def test_the_capacity_interval_is_wide_because_of_it(self):
        lo, hi = ELIZONDO2025.capacity_ci95

        assert hi / lo > 1.5

    def test_a_registered_prediction_exists_for_a_third_rate(self):
        """Written before any measurement at D = 0.18, so it can be scored rather than
        rationalised.

        Re-registered from the GENOTYPE. The first version took the pathway flux as the
        strain's own measured ``q_lycopene + q_betacarotene``, which made a prediction
        registered in advance rest on a measurement of the quantity being predicted --
        registering early does not make that not circular. The columns changed with it:
        ``ci95_*`` on a fitted capacity became ``flux_*``, the leave-one-strain-out band on
        the flux law, which is the spread a user of the number actually faces."""
        path = paths.outputs_dir() / "registered_prediction_D018.csv"
        if not path.exists():
            pytest.skip(f"{path} not present")
        registered = pd.read_csv(path)

        assert len(registered) == 3
        # `dilution_rate`, not `growth_rate_setpoint_per_h`: this reads the REGISTERED
        # table, whose schema is frozen at what was predicted under the chemostat plan.
        # The Environment field was renamed on 2026-09-02; the sealed record was not.
        assert (registered.dilution_rate == 0.18).all()
        assert (registered.flux_low < registered.predicted_pathway_flux).all()
        assert (registered.predicted_pathway_flux < registered.flux_high).all()

    def test_and_nothing_measured_about_the_product_went_into_it(self):
        """The property that makes it worth registering at all."""
        path = paths.outputs_dir() / "registered_prediction_D018.csv"
        if not path.exists():
            pytest.skip(f"{path} not present")
        registered = pd.read_csv(path)

        assert "entry_expression" in registered.columns
        assert not any("q_" in c or "measured" in c for c in registered.columns)
