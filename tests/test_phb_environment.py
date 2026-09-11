"""The environment layer scored against a real environment series, and failing.

Kocharin & Nielsen 2013 (PMID 23514405, PMC3610212, CC-BY) is the first dataset in this
repository that can test the environment layer at all: eleven aerobic chemostat steady
states, ONE genotype (SCKK006), four dilution rates crossed with three carbon feeds matched
at 0.666 Cmol/L. One genotype across eleven environments is what makes "flux is a strain
constant" falsifiable rather than fitted.

It falsifies it, and these tests pin the three facts that do the falsifying, so that a
later change claiming to fix the environment layer has something concrete to beat.

**What is NOT tested here, deliberately.** PHB carries no expression measurement, so the
flux law cannot be run against it. And the solver's prediction for a passthrough chain is
``content = flux / mu``, the identity Kocharin used to report their own numbers -- checking
it would test arithmetic. The one testable question is whether anything in the environment
predicts the flux with the genotype held fixed, and the answer is no.
"""

from __future__ import annotations

import importlib.util
import pathlib

import numpy as np
import pandas as pd
import pytest

from ystwin import paths

_REPO = pathlib.Path(__file__).resolve().parents[1]
_STATES = paths.data_dir() / "phb" / "kocharin2013_chemostat_states.tsv"


@pytest.fixture(scope="module")
def states():
    if not _STATES.exists():
        pytest.skip(f"{_STATES} not present")
    return pd.read_csv(_STATES, sep="\t")


@pytest.fixture(scope="module")
def script():
    path = _REPO / "scripts" / "score_phb_environment.py"
    spec = importlib.util.spec_from_file_location("score_phb_environment", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class TestTheDatasetIsTheOneClaimed:
    def test_eleven_states(self, states):
        assert len(states) == 11

    def test_four_dilution_rates_and_three_feeds(self, states):
        assert sorted(states.dilution_rate_per_h.unique()) == [0.05, 0.10, 0.15, 0.20]
        assert len(states.carbon_source.unique()) == 3

    def test_ethanol_washes_out_at_the_top_rate(self, states):
        """Which is why it is eleven states and not twelve, and is itself a physiology
        check: at D = 0.20 an ethanol-limited culture cannot hold."""
        ethanol = states[states.carbon_source == "ethanol"]

        assert 0.20 not in set(ethanol.dilution_rate_per_h)

    def test_the_steady_state_identity_holds_in_their_own_numbers(self, states):
        """mu = D exactly in a chemostat, so q = mu * content. If this failed, the table
        would be internally inconsistent and nothing below would mean anything."""
        residual = (states.q_phb_mmol_per_gdcw_h
                    - states.dilution_rate_per_h * states.phb_mmol_per_gdcw).abs()

        assert residual.max() < 1e-6


class TestTheEnvironmentMovesTheFluxAndTheModelCannotSaySo:
    def test_carbon_source_moves_it_fourfold_at_identical_growth_rate(self, states):
        """The single most damaging number for the current architecture. Environment reaches
        the product only through mu, and mu is identical across this row."""
        at_low = states[states.dilution_rate_per_h == 0.05]

        spread = at_low.q_phb_mmol_per_gdcw_h.max() / at_low.q_phb_mmol_per_gdcw_h.min()

        assert spread > 4.0

    def test_and_the_model_has_no_term_that_could_carry_it(self):
        """Stated as an import check because that is what it comes down to: nothing in the
        prediction chain takes a carbon source. If a term is ever added this fails, which is
        the intended signal."""
        from ystwin.pathway import flux

        source = (pathlib.Path(flux.__file__)).read_text()

        assert "carbon" not in source.lower()


class TestTheGrowthDependenceIsNotOneLaw:
    def test_the_exponent_differs_by_more_than_one_across_feeds(self, states, script):
        """Glucose gives roughly mu^1.2, ethanol and the mix roughly mu^0. A single fitted
        exponent -- which is what `Vmax = capacity * mu` asserts -- describes none of them."""
        exponents = script.growth_exponents(states)

        spread = exponents.exponent.max() - exponents.exponent.min()

        assert spread > 1.0

    def test_glucose_rises_with_growth_rate(self, states, script):
        exponents = script.growth_exponents(states).set_index("carbon_source")

        assert exponents.loc["glucose", "exponent"] > 0.8

    def test_while_the_other_two_are_flat(self, states, script):
        exponents = script.growth_exponents(states).set_index("carbon_source")

        for feed in ("ethanol", "glucose_ethanol_1to2"):
            assert abs(exponents.loc[feed, "exponent"]) < 0.3

    def test_which_contradicts_the_beta_carotene_reading(self, states, script):
        """Not fatal, but it is the reason the beta-carotene growth law cannot be assumed to
        transfer. There the exponent fitted to 1.04 with a 95% interval of [0.53, 1.78];
        here two of three feeds sit outside that interval entirely."""
        exponents = script.growth_exponents(states).set_index("carbon_source")
        outside = [f for f in exponents.index if not 0.53 <= exponents.loc[f, "exponent"] <= 1.78]

        assert len(outside) == 2


class TestNothingAvailablePredictsIt:
    def test_no_single_measured_variable_correlates_strongly(self, states):
        """The honest summary. Best is about +0.47, on eleven points, which is not a model."""
        log_flux = np.log(states.q_phb_mmol_per_gdcw_h)
        ethanol_fraction = (states.feed_ethanol_g_per_l
                            / (states.feed_ethanol_g_per_l + states.feed_glucose_g_per_l))
        candidates = [states.dilution_rate_per_h, np.log(states.dilution_rate_per_h),
                      states.ysx_cmol_per_cmol, ethanol_fraction]

        best = max(abs(np.corrcoef(c, log_flux)[0, 1]) for c in candidates)

        assert best < 0.6

    def test_a_per_feed_constant_barely_beats_a_global_one(self, states):
        """Grouping by carbon source is the most generous description available, and it
        removes under a fifth of the residual. There is no environment term hiding here that
        eleven states would identify."""
        log_flux = np.log(states.q_phb_mmol_per_gdcw_h)
        per_feed = states.groupby("carbon_source").q_phb_mmol_per_gdcw_h.transform(
            lambda s: np.log(s).mean())

        grouped = float(np.sqrt(np.mean((log_flux - per_feed) ** 2)))
        global_ = float(np.sqrt(np.mean((log_flux - log_flux.mean()) ** 2)))

        assert grouped > 0.7 * global_


class TestNothingBeatsAConstantOutOfSample:
    """The strongest form of the negative result, and a correction to a weaker one.

    An earlier version of `docs/EXTERNAL_PRODUCT_VALIDATION.md` reported that carbon source
    "removes 28.7% of the variance". That is an IN-SAMPLE number and on eleven points it
    means very little. Scored leave-one-out, carbon source is WORSE than predicting the
    mean -- and so is dilution rate, and so is both together.

    The permutation null is what makes that a result rather than an absence. Noise alone,
    with the same number of parameters on the same eleven points, reaches a 95th-percentile
    skill of about -0.05. Every candidate sits inside that. So the honest statement is not
    "we have not found the term yet" but "eleven states cannot identify one".
    """

    @staticmethod
    def _loo_skill(observed, design):
        n = len(observed)
        residuals = np.empty(n)
        for held_out in range(n):
            keep = np.arange(n) != held_out
            coefficients, *_ = np.linalg.lstsq(design[keep], observed[keep], rcond=None)
            residuals[held_out] = design[held_out] @ coefficients - observed[held_out]
        model = float(np.sqrt(np.mean(residuals ** 2)))
        baseline = float(np.sqrt(np.mean((observed - observed.mean()) ** 2)))
        return 1.0 - model / baseline

    def test_carbon_source_is_worse_than_a_constant_out_of_sample(self, states):
        observed = np.log(states.q_phb_mmol_per_gdcw_h.to_numpy())
        feeds = pd.get_dummies(states.carbon_source, drop_first=True).to_numpy(dtype=float)
        design = np.column_stack([np.ones(len(states)), feeds])

        assert self._loo_skill(observed, design) < 0.0

    def test_and_so_is_dilution_rate(self, states):
        observed = np.log(states.q_phb_mmol_per_gdcw_h.to_numpy())
        design = np.column_stack([np.ones(len(states)),
                                  np.log(states.dilution_rate_per_h.to_numpy())])

        assert self._loo_skill(observed, design) < 0.0

    def test_and_so_is_both_together(self, states):
        """Adding parameters makes it worse, which is the signature of fitting noise."""
        observed = np.log(states.q_phb_mmol_per_gdcw_h.to_numpy())
        feeds = pd.get_dummies(states.carbon_source, drop_first=True).to_numpy(dtype=float)
        design = np.column_stack([np.ones(len(states)), feeds,
                                  np.log(states.dilution_rate_per_h.to_numpy())])

        assert self._loo_skill(observed, design) < 0.0

    def test_the_in_sample_number_that_used_to_be_quoted_is_much_more_flattering(self,
                                                                                 states):
        """Both are true and only one is a result. Keeping the contrast in a test is how the
        distinction survives the next person who computes an R-squared on eleven points."""
        log_flux = np.log(states.q_phb_mmol_per_gdcw_h)
        in_sample = 1.0 - float(np.var(
            log_flux - states.groupby("carbon_source").q_phb_mmol_per_gdcw_h.transform(
                lambda s: np.log(s).mean()))) / float(np.var(log_flux))

        feeds = pd.get_dummies(states.carbon_source, drop_first=True).to_numpy(dtype=float)
        out_of_sample = self._loo_skill(
            log_flux.to_numpy(), np.column_stack([np.ones(len(states)), feeds]))

        assert in_sample > 0.25 and out_of_sample < 0.0
