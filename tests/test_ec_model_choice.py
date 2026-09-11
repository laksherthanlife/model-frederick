"""Why the enzyme-constrained model stays on yeast-GEM 8.3.4 while the main model is 9.0.2.

The mismatch is documented in `data/gem/MANIFEST.md` as load-bearing, and on 2026-09-03 it
was closed by building an ecYeastGEM on 9.0.2 with GECKO 3.2.5 -- the first on any 9.x, since
SysBioChalmers/ecModels still ships 8.3.4 and GECKO's own tutorial ships 8.6.2. That build is
vendored beside the old one and is NOT the default.

These tests exist because "the newer model is on the same base version and predicts growth
better" is a good argument that reaches the wrong conclusion, and nothing in the code would
otherwise stop someone acting on it. The 9.0.2 build is fully RESPIRATORY: it reaches the
same growth with no ethanol at all. Every culture this repository is calibrated against
ferments, so that is a qualitative failure a percentage point of growth accuracy does not
buy back.

Marked `integration`: each test solves a real GSMM.
"""

from __future__ import annotations

import pytest

from ystwin import paths

pytestmark = pytest.mark.integration

_NEW = "ecYeastGEM_yeast902.xml.gz"
_OLD = "ecYeastGEM_batch.xml.gz"


@pytest.fixture(scope="module")
def both():
    import cobra

    directory = paths.data_dir() / "gem"
    old, new = directory / _OLD, directory / _NEW
    if not (old.exists() and new.exists()):
        pytest.skip("both ec models must be vendored to compare them")
    return (cobra.io.read_sbml_model(str(old)), cobra.io.read_sbml_model(str(new)))


def _batch_fluxes(model):
    """Unlimited carbon and oxygen, maximise growth. The protein pool should be what binds."""
    from ystwin.fba.physiology import cap_uptake

    with model:
        cap_uptake(model, "r_1714", 1000.0)
        cap_uptake(model, "r_1992", 1000.0)
        for objective in ("r_2111", "r_4041"):
            if objective in {r.id for r in model.reactions}:
                model.objective = objective
                break
        solution = model.optimize()

        def flux(base):
            ids = {r.id for r in model.reactions}
            for candidate in (base, f"{base}_REV"):
                if candidate in ids and abs(solution.fluxes.get(candidate, 0.0)) > 1e-9:
                    return abs(solution.fluxes[candidate])
            return 0.0

        return {"growth": solution.objective_value, "glucose": flux("r_1714"),
                "oxygen": flux("r_1992"), "ethanol": flux("r_1761")}


class TestTheDefaultIsStillTheOlderModel:
    def test_ec_yeast_gem_resolves_to_the_8_3_4_build(self):
        resolved = paths.ec_yeast_gem()

        assert resolved is not None
        assert resolved.name == _OLD

    def test_the_9_0_2_build_is_vendored_and_reachable(self):
        """Kept, not deleted -- it is the better model on growth and the record of the
        build. One environment variable selects it."""
        assert (paths.data_dir() / "gem" / _NEW).exists()


class TestWhyTheNewerModelIsNotTheDefault:
    def test_the_9_0_2_build_produces_no_ethanol(self, both):
        """The disqualifying result. Elizondo's six calibration states secrete 2.9-16.1
        mmol/gDCW/h of ethanol at BOTH dilution rates."""
        _, new = both

        assert _batch_fluxes(new)["ethanol"] == pytest.approx(0.0, abs=1e-6)

    def test_the_8_3_4_model_does_ferment(self, both):
        old, _ = both

        assert _batch_fluxes(old)["ethanol"] > 20.0

    def test_the_newer_model_is_genuinely_better_on_growth(self, both):
        """Stated so the trade-off is not misremembered as the new model simply being worse.
        Against van Hoek 1998's 0.40 /h it is 4.5% out against the older model's 5.8%."""
        old, new = both

        assert abs(_batch_fluxes(new)["growth"] - 0.40) < abs(_batch_fluxes(old)["growth"] - 0.40)

    def test_it_reaches_that_growth_respiratively(self, both):
        """The mechanism behind the zero: it takes far less glucose and far more oxygen,
        so respiration is affordable and overflow never becomes optimal."""
        old, new = both
        old_flux, new_flux = _batch_fluxes(old), _batch_fluxes(new)

        assert new_flux["glucose"] < 0.5 * old_flux["glucose"]
        assert new_flux["oxygen"] > 2.0 * old_flux["oxygen"]


class TestBothConventionsRemainSupported:
    def test_cap_uptake_constrains_the_right_reaction_in_each(self, both):
        """GECKO 2 pins r_1714 at (0,0) and supplies through r_1714_REV; GECKO 3 leaves
        r_1714 unsplit. `cap_uptake` reads that from the model, which is why swapping the
        default needed no code change -- and why a future swap still will not."""
        from ystwin.fba.physiology import cap_uptake

        old, new = both
        with old:
            assert cap_uptake(old, "r_1714", 10.0) == "r_1714_REV"
        with new:
            assert cap_uptake(new, "r_1714", 10.0) == "r_1714"
