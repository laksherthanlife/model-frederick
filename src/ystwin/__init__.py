"""Recover promoter activity from reporter and optical-density time series.

The measurement problem this solves is not specific to any product or any host. A
fluorescent reporter accumulates in a cell and is split between daughters at every
division, so its per-cell concentration relaxes to ``k_synth / (mu + k_deg)``: a culture
that slows down grows brighter with no change in transcription at all. Since almost
every perturbation worth studying also slows growth, raw ``RFU/OD`` confounds induction
with dilution, and the two have to be separated before anything downstream means
anything.

``reporter.py`` inverts that relation, ``growth.py`` supplies the rate it needs,
``observation.py`` maps state to instrument reading, and the gates refuse the wells where
the optics do not support the question. None of that is organism-specific.

What *is* specific sits behind it and is named as such: the seven-regulon stress panel in
``generator/`` is S. cerevisiae, the GSMM validation in ``fba/`` is Yeast9, the reader in
``plate/synergy.py`` is one manufacturer's export format, and the carotenoid pathway is a
worked application parked until a strain exists to test it (see ``docs/PARKED.md``).
Those are instances, not the subject.
"""

__version__ = "0.1.0"
