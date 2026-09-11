"""Whether a sensor has actually been demonstrated in this organism.

The panel carried Peredox as an available channel. It has never been used in
S. cerevisiae -- not once, across every query formulation a literature sweep could put to
PubMed and Europe PMC. Neither has SoNar. Recommending a build around a sensor nobody has
got working in yeast is recommending a research project, not an experiment, and the model
had no way to say so.

What exists in yeast is narrower than the sensor catalogue suggests. roGFP2-Grx1 is
thoroughly characterised, with cytosolic E_GSH reported between -300 and -320 mV across
several labs. QUEEN has a yeast paper. pHluorin is standard. The NADH channel in yeast is a
promoter-fusion reporter, not a ratiometric sensor, so it cannot report a ratio at all.

A build the wet lab cannot make is worse than a narrower one it can, so availability is a
property the selector has to see.
"""


from ystwin.analysis.sensor_selection import RECOMMENDED_BUILD, select_sensors
from ystwin.analysis.design import NoiseModel
from ystwin.generator.stress_panel import REPORTERS, Kind, demonstrated_in_yeast

PLATE = NoiseModel(relative_cv=0.05)


class TestEveryReporterDeclaresWhetherItWorksInYeast:
    def test_all_of_them_say(self):
        assert all(isinstance(r.demonstrated_in_yeast, bool) for r in REPORTERS.values())

    def test_a_promoter_fusion_is_available_by_construction(self):
        """Any promoter can be fused to any fluorophore; nothing has to be ported."""
        for name, reporter in REPORTERS.items():
            if reporter.kind is Kind.TRANSCRIPTIONAL:
                assert reporter.demonstrated_in_yeast, name

    def test_the_glutathione_sensor_is_established(self):
        assert REPORTERS["roGFP2-Grx1"].demonstrated_in_yeast

    def test_the_nadh_sensor_is_not(self):
        """Peredox has no published application in S. cerevisiae."""
        assert not REPORTERS["Peredox"].demonstrated_in_yeast

    def test_the_unavailable_ones_say_why_in_their_source(self):
        for name, reporter in REPORTERS.items():
            if not reporter.demonstrated_in_yeast:
                assert "yeast" in reporter.source.lower(), name


class TestTheHelper:
    def test_it_lists_only_what_has_been_shown_to_work(self):
        available = demonstrated_in_yeast()

        assert "roGFP2-Grx1" in available
        assert "Peredox" not in available

    def test_every_name_it_returns_is_a_real_reporter(self):
        assert set(demonstrated_in_yeast()) <= set(REPORTERS)

    def test_it_is_a_strict_subset_of_the_catalogue(self):
        assert len(demonstrated_in_yeast()) < len(REPORTERS)


class TestSelectionCanBeHeldToWhatExists:
    def test_it_can_be_restricted_to_demonstrated_sensors(self):
        chosen = select_sensors(n_channels=5, noise=PLATE, min_effect=0.5, growth_cv=0.05,
                                library=demonstrated_in_yeast())

        assert all(REPORTERS[r].demonstrated_in_yeast for r in chosen)

    def test_the_recommended_build_only_uses_what_exists(self):
        """The build is a recommendation to a wet lab, so it has to be buildable today."""
        assert all(REPORTERS[r].demonstrated_in_yeast
                   for r in RECOMMENDED_BUILD.stress_reporters)

    def test_the_restriction_turns_out_to_cost_nothing(self):
        """Worth pinning, because it says no conclusion rested on the unavailable sensor.

        From four channels to nine the selector picks the same set whether or not it is held
        to what works in yeast. Peredox is never chosen: it reads a pool no other channel
        reaches, but the modules it would add are ones nothing in the stressor panel drives
        strongly enough to be worth a channel. The availability constraint is therefore a
        guard against a future mistake rather than a correction to a past one.
        """
        for channels in range(4, 10):
            free = select_sensors(n_channels=channels, noise=PLATE, min_effect=0.5,
                                  growth_cv=0.053)
            held = select_sensors(n_channels=channels, noise=PLATE, min_effect=0.5,
                                  growth_cv=0.053, library=demonstrated_in_yeast())
            assert set(free) == set(held), channels
