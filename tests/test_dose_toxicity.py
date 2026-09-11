"""Dose response is biphasic, and the generator's monotone Hill could not produce it.

The uploaded plates settle this. Corrected promoter activity rises to a peak around 1 mM
and then falls -- UPRE1 goes 1.00, 1.15, 1.27, 1.38 and back to 0.91 at 5 mM DTT -- because
a dose high enough to induce strongly is also high enough to stop the cell transcribing.
A saturating Hill can only ever approach a ceiling, so it fits the rising limb and denies
the falling one.

The two agents fail differently and the difference is not cosmetic. At 2 mM H2O2 the wells
are dead: growth is near zero and blank-corrected fluorescence goes negative, so there is
no activity to model, only a measurement that does not exist. Under DTT the cells are alive
but slowed, and there the naive readout inverts the answer -- raw fluorescence per OD keeps
climbing to 2.17-fold while true activity is falling, because reporter accumulates in cells
that have stopped diluting it.
"""

import numpy as np
import pytest

from ystwin.generator.stress_panel import STRESSORS, module_response, viability


class TestViability:
    def test_an_undosed_culture_is_fully_viable(self):
        assert viability("DTT", 0.0) == pytest.approx(1.0)

    def test_a_mild_dose_barely_costs_viability(self):
        assert viability("DTT", 0.1) > 0.9

    def test_a_lethal_dose_costs_most_of_it(self):
        assert viability("H2O2", 4.0) < 0.3

    def test_it_falls_monotonically_with_dose(self):
        values = [viability("DTT", d) for d in (0.0, 0.5, 1.0, 2.0, 5.0, 20.0)]

        assert all(a >= b for a, b in zip(values, values[1:]))

    def test_it_stays_a_fraction(self):
        for name in STRESSORS:
            for dose in (0.0, 1.0, 1e3, 1e6):
                assert 0.0 <= viability(name, dose) <= 1.0

    def test_every_stressor_declares_where_it_becomes_lethal(self):
        for name, spec in STRESSORS.items():
            assert spec.lethal_dose > spec.ec50, name

    def test_an_unknown_stressor_is_refused(self):
        with pytest.raises(KeyError, match="unobtainium"):
            viability("unobtainium", 1.0)


class TestTheResponseIsBiphasic:
    def test_activity_peaks_and_then_falls(self):
        doses = np.array([0.0, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0])
        response = np.array([module_response("DTT", d)["UPR"] for d in doses])

        assert response.argmax() not in (0, len(doses) - 1)

    def test_the_peak_sits_between_the_inducing_and_lethal_doses(self):
        """Where exactly it sits is a measurement, and is checked against the plates in
        test_panel_calibration; here the claim is only that it sits between the two."""
        spec = STRESSORS["DTT"]
        doses = np.linspace(0.01, 4 * spec.lethal_dose, 400)
        response = np.array([module_response("DTT", d)["UPR"] for d in doses])
        peak = doses[response.argmax()]

        assert spec.ec50 * 0.5 < peak < spec.lethal_dose * 2.0

    def test_a_high_dose_gives_less_than_the_peak(self):
        peak = module_response("DTT", 1.0)["UPR"]

        assert module_response("DTT", 5.0)["UPR"] < peak

    def test_the_rising_limb_still_rises(self):
        assert module_response("H2O2", 0.5)["oxidative"] > module_response("H2O2", 0.1)["oxidative"]

    def test_it_still_returns_nothing_at_zero_dose(self):
        assert all(v == 0.0 for v in module_response("H2O2", 0.0).values())

    def test_viability_scales_every_module_together(self):
        """Toxicity stops the cell transcribing, so it cannot spare one regulon."""
        mild, harsh = module_response("DTT", 1.0), module_response("DTT", 20.0)
        driven = [m for m in mild if abs(mild[m]) > 1e-9]

        ratios = [harsh[m] / mild[m] for m in driven]
        assert max(ratios) - min(ratios) < 0.05
