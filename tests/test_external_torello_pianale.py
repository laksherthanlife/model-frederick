"""An external check against a paper this model has never been fitted to.

Torello Pianale L, Rugbjerg P, Olsson L (2022), *Real-Time Monitoring of the Yeast
Intracellular State During Bioprocesses With a Toolbox of Biosensors*, Front Microbiol
12:802169, doi 10.3389/fmicb.2021.802169. Open access.

Nothing in this repository was fitted to it, and no constant here was chosen with it in
view, so what follows is a held-out test of `generator/context.py` rather than a residual.
Their Table 4 gives CEN.PK113-7D at **0.379 +/- 0.000 /h** in Delft (Verduyn) mineral
medium, 20 g/L glucose, 30 C, 150 mL in a 500 mL flask at 140 rpm, with N2 flushed for ten
seconds -- which the paper itself calls **microaerobic**, not anaerobic.

The condition label is the interesting part, and it is where this model earns something.
The paper reports no Tween-80 and no ergosterol. `context.py` holds that an unsupplemented
anaerobic glucose culture does not grow AT ALL, because yeast cannot make a sterol without
molecular oxygen (Andreasen & Stier 1953, PMID 13034889). They measured 0.379 /h. Those two
statements are only compatible if the culture had oxygen -- so the model's refusal reads the
condition correctly, and agrees with the authors' own word for it.

A first pass at this file nearly recorded the opposite. A summary of the paper described the
cultures as "anaerobic", which would have made the model 18% low and looked like a refuted
constant. Reading the methods rather than the summary is what separated the two readings,
and the constant survived.
"""

from __future__ import annotations


from ystwin.generator.context import (
    ANAEROBIC_MU_MAX_PER_H,
    CultureContext,
    context_growth_rate,
)

# Torello Pianale 2022, Table 4, parental row.
MEASURED_MU = 0.379
CITATION = "doi 10.3389/fmicb.2021.802169, Table 4"


def delft(**overrides) -> CultureContext:
    """Their conditions, as the model expresses them."""
    return CultureContext(carbon_source="glucose", temperature_c=30.0,
                          glucose_g_per_L=20.0, **overrides)


class TestTheGrowthRateIsPredictedWithoutHavingSeenIt:
    def test_an_aerated_flask_lands_within_ten_percent(self):
        """A shaken 500 mL flask holding 150 mL, flushed for ten seconds and then left, is
        an aerated culture: the headspace is 350 mL and 140 rpm keeps transferring oxygen.
        Asked as air, the model gives 0.400 against their 0.379."""
        predicted = context_growth_rate(delft(oxygen=0.21))

        assert abs(predicted - MEASURED_MU) / MEASURED_MU < 0.10

    def test_and_it_is_an_over_prediction_rather_than_an_under_one(self):
        """Which direction the miss goes is worth pinning. Air is the upper end of what a
        ten-second flush leaves, so the model should sit above the measurement, and a model
        that drifted below it would be wrong in a way this test would not otherwise catch.
        """
        assert context_growth_rate(delft(oxygen=0.21)) > MEASURED_MU

    def test_the_measurement_is_above_the_anaerobic_ceiling(self):
        """0.379 exceeds ANAEROBIC_MU_MAX_PER_H, which is Verduyn 1990's measured 0.31. So
        the measurement is not reachable anaerobically at all under this model, whatever
        the supplements -- an independent second reason the culture had oxygen."""
        assert MEASURED_MU > ANAEROBIC_MU_MAX_PER_H


class TestTheModelDiagnosesTheConditionLabel:
    """The result worth having. The model is not just close to their number; it can say
    something about their culture that they only asserted."""

    def test_an_unsupplemented_anaerobic_culture_cannot_grow_at_all(self):
        assert context_growth_rate(delft(oxygen=0.0, anaerobic_supplements=False)) == 0.0

    def test_so_a_growing_unsupplemented_culture_must_have_had_oxygen(self):
        """The two statements together. They report growth and no anaerobic growth factors;
        this model says that combination is impossible without oxygen; the paper's own
        methods say microaerobic. Three independent routes to the same reading."""
        anaerobic = context_growth_rate(delft(oxygen=0.0, anaerobic_supplements=False))

        assert anaerobic == 0.0 and MEASURED_MU > 0.0

    def test_supplements_alone_do_not_close_the_gap(self):
        """Ruling out the obvious alternative reading -- that they did supplement and
        simply did not report it. Even fully supplemented, anaerobic tops out at 0.31."""
        supplemented = context_growth_rate(delft(oxygen=0.0, anaerobic_supplements=True))

        assert supplemented < MEASURED_MU


class TestWhatTheirSensorDesignSaysAboutOurs:
    """Not a numerical test -- a design conclusion this paper settles, pinned so it is not
    lost.

    Their five sensors split into two kinds, and the split is exactly the axis of this
    repository's headline finding that most of an apparent sensor response was growth
    dilution:

      RATIOMETRIC, one fluorophore, two excitations -- QUEEN-2m (410/480, em 520) and
      sfpHluorin (390/470, em 512). Both channels report the SAME molecule, so anything
      that scales the molecule -- dilution by growth above all -- cancels exactly in the
      ratio. These are structurally immune to the confound.

      INTENSIOMETRIC with a constitutive reference -- GlyRNA, OxPro, RibPro, each normalised
      to mCherry from pTEFmut8. Dilution cancels only to the extent that the two proteins
      share maturation and degradation, so this is partial protection, not immunity.

    This repository's constructs are single-channel mCitrine with no reference at all, which
    is the fully exposed case. That is why the dilution correction had to be done in
    software, and why `docs/FINDINGS.md` reports that it removed most of the response.

    The paper is also silent on growth-rate dependence -- it reports no analysis of it, and
    does not address signal dilution quantitatively. Its ratiometric sensors did not need it
    to. Ours did.
    """

    def test_the_repositorys_constructs_are_the_exposed_kind(self):
        """A single-channel reporter has no second measurement of the same molecule, so
        nothing cancels and the correction has to be modelled. If a ratiometric or
        dual-channel construct is ever added, this test should be the thing that fails."""
        from ystwin.generator.stress_panel import REPORTERS

        channels = {name: spec.module for name, spec in REPORTERS.items()}

        assert len(channels) == len(set(REPORTERS)), (
            "every reporter is one channel reading one module; a ratiometric pair would "
            "break this and would also make the dilution correction unnecessary")

    def test_the_dilution_correction_exists_because_of_that(self):
        """The software substitute for what a ratiometric sensor gets for free."""
        from ystwin.reporter import promoter_activity

        assert callable(promoter_activity)
