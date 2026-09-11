"""Where the data lives, resolved rather than hardcoded.

Every real-data path in this package used to be written inline, and two of them
pointed at one person's home directory. That made every real-data result in the
README unreproducible by anyone else, including this machine after a reorganise.

Each asset resolves to the environment variable named below if that variable is set, and
otherwise to the first location that exists out of

1. a path inside this repository, for assets small enough to vendor;
2. the sibling-directory convention under the IGEM project root, which is where
   the wet lab's exports and the genome-scale models actually live.

The variable is checked first and, when set, it is checked *instead* rather than *before*:
a variable pointing at a directory that does not exist resolves to nothing, not to the
default. Falling through would mean a machine explicitly configured to read one dataset
could silently read another and report the difference as a result, which is the failure
this module exists to prevent -- and it also makes the override useless for testing what
happens when the data is absent.

Nothing here downloads, creates, or guesses. A resolver returns ``None`` when the
asset is absent so callers can skip -- an analysis that silently substitutes a
different file is worse than one that does not run.
"""

from __future__ import annotations

import os
import pathlib

__all__ = [
    "PACKAGE_ROOT",
    "REPO_ROOT",
    "IGEM_ROOT",
    "atp_sensor_plates",
    "biosensor_plates",
    "data_dir",
    "display_path",
    "ec_yeast_gem",
    "crosstalk_workbook",
    "gen5_xpt_dir",
    "igem_results",
    "outputs_dir",
    "qpcr_dir",
    "qpcr_raw_dir",
    "require",
    "resolve_or_exit",
    "thermo_dir",
    "yeast_gem",
]

PACKAGE_ROOT = pathlib.Path(__file__).resolve().parent
"""``src/ystwin``."""

REPO_ROOT = PACKAGE_ROOT.parents[1]
"""The ``model-v2`` checkout."""

IGEM_ROOT = REPO_ROOT.parent
"""The directory holding ``model-v2`` and its sibling data directories."""


def _resolve(env_var: str, *defaults: pathlib.Path) -> pathlib.Path | None:
    """The override if it is set at all, otherwise the first default that exists.

    An override that is set wins outright, even when it points at nothing -- see the
    module docstring for why it must not fall through to a default.

    Args:
        env_var: Name of the variable that overrides this asset's location.
        defaults: Candidates in preference order, used only when the variable is unset.

    Returns:
        The resolved path, or ``None`` if the asset is nowhere the caller said it is.
    """
    override = os.environ.get(env_var)
    if override:
        path = pathlib.Path(override).expanduser()
        return path.resolve() if path.exists() else None
    for candidate in defaults:
        if candidate.exists():
            return candidate
    return None


def data_dir() -> pathlib.Path:
    """Vendored data that ships with the repository."""
    return REPO_ROOT / "data"


def display_path(path: pathlib.Path) -> str:
    """A path as a reader wants to see it: repo-relative when it is inside, absolute when not.

    Exists because ``path.relative_to(REPO_ROOT)`` RAISES on a path outside the repository,
    and :func:`outputs_dir` documents ``YSTWIN_OUTPUTS`` as the supported way to redirect a
    table into a scratch directory. Seven scripts printed their destination with
    ``relative_to`` and so crashed under the very override this module offers -- found on
    2026-09-04 by `scripts/audit_output_tables.py`, whose entire method is to run each
    generator with ``YSTWIN_OUTPUTS`` pointed elsewhere and diff the result.
    """
    try:
        return str(path.relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def outputs_dir() -> pathlib.Path:
    """Where scripts write their tables. ``YSTWIN_OUTPUTS`` overrides.

    The one path here that is created rather than resolved, because it is a destination
    and not a measurement -- so it never returns ``None`` and never goes through
    :func:`_resolve`. Redirecting it is how you regenerate a table without overwriting
    the tracked one.
    """
    out = pathlib.Path(os.environ.get("YSTWIN_OUTPUTS") or REPO_ROOT / "outputs")
    out.mkdir(parents=True, exist_ok=True)
    return out


def crosstalk_workbook() -> pathlib.Path | None:
    """`ER&OX-summary.xlsx`, needed only to re-derive the crosstalk tables.

    ``YSTWIN_CROSSTALK_WORKBOOK`` is the only way to find it, and the absence of a default
    is the point. The workbook carries ``dc:creator`` and ``cp:lastModifiedBy`` naming a
    private individual, so what is tracked is the derived table under ``data/crosstalk/``
    and never the binary -- and every consumer reads the derived table, so nothing is lost
    by this returning ``None``, which is the ordinary case.

    It was written with a ``~/Downloads`` default for about ten minutes, and both gates
    caught it: ``audit_claims.py`` for a home path outside this module, then
    ``audit_reproducibility.py`` for resolving outside the checkout by convention. The
    second is the one that matters -- a path found by convention works on one machine -- and
    it is the same finding that took ``igem_results`` its fallback away an hour earlier.
    """
    return _resolve("YSTWIN_CROSSTALK_WORKBOOK")


def gen5_xpt_dir() -> pathlib.Path | None:
    """The Gen5 ``.xpt`` instrument files. ``YSTWIN_GEN5_XPT`` is the only way to find them.

    No default, for the reason ``igem_results`` and ``crosstalk_workbook`` have none: a path
    found by convention resolves on one machine, and `audit_reproducibility.py` is right to
    fail on it.

    These are the RAW experiment files, upstream of every ``.xlsx`` this project has seen.
    They matter for one specific thing: plate 20260804's exports are all blank-subtracted, so
    the panel sits at n=3, and the ``.xpt`` still carries the background -- H1-H3 read 0.087
    OD and 347 RFU where the export has them at +-0.004 and +-15. ``plate/gen5.py`` reads
    them, and reproduces the committed text for every plate that has both.
    """
    return _resolve("YSTWIN_GEN5_XPT")


def igem_results() -> pathlib.Path | None:
    """The two 2026-07 Synergy **workbooks**, needed only to re-export or re-verify.

    ``YSTWIN_IGEM_RESULTS`` is the only way to find them, and the absence of a default is
    the point rather than an omission.

    This resolver used to fall back to ``../igem-results``, a sibling of the checkout, and
    `audit_reproducibility.py` failed on it for months with the right reason: a path found
    by convention resolves on one machine, so every result computed through it was
    reproducible in exactly one place with nothing in the repository to say so. The
    convention existed because the workbooks cannot be committed -- ``20260709`` carries a
    named private individual in ``cp:lastModifiedBy`` -- and so there was nowhere else for
    the numbers to come from.

    There is now. ``data/plates/`` holds both exports as text, value-for-value, verified by
    ``plate_readings.py --export`` and carrying each workbook's sha256 so a holder of the
    original can prove provenance without the original being here -- the same treatment the
    NewProtocol replicates already had, and the reason this fallback is no longer load
    bearing. Every consumer reads the committed text through
    :func:`ystwin.plate.replay.install` when this returns ``None``, which is now the
    ordinary case rather than the degraded one.
    """
    return _resolve("YSTWIN_IGEM_RESULTS")


def biosensor_plates() -> pathlib.Path | None:
    """The NewProtocol biosensor replicates. ``YSTWIN_PLATES`` overrides.

    These are the usable set: four biological replicates of seven doses by three
    technical replicates by four constructs.
    """
    return _resolve(
        "YSTWIN_PLATES",
        IGEM_ROOT / "plates" / "Biosensor Testing",
        pathlib.Path.home() / "Desktop/Result & Analysis/Plate Reader Result/Biosensor Testing",
        pathlib.Path.home() / "Desktop/iGEM 2026/Results (excel sheets)",
    )


def atp_sensor_plates() -> pathlib.Path | None:
    """The two ATP-sensor Synergy exports. ``YSTWIN_ATP_SENSOR`` overrides.

    Vendored under ``data/atp_sensor`` but not tracked, because they are wet-lab exports
    renamed by hand from the reader's own filenames. Anyone else has to be told where
    theirs are rather than left to discover an empty directory.
    """
    return _resolve(
        "YSTWIN_ATP_SENSOR",
        data_dir() / "atp_sensor",
        IGEM_ROOT / "plates" / "atp_sensor",
    )


def thermo_dir() -> pathlib.Path | None:
    """ModelSEED alias tables and the pytfa thermo database. ``YSTWIN_THERMO`` overrides.

    A third-party download rather than anyone's measurement, and large enough that it is
    not tracked. Resolved rather than read blind so that its absence produces the name of
    the variable to set instead of a bare missing-file error three frames down.
    """
    return _resolve("YSTWIN_THERMO", data_dir() / "thermo")


def qpcr_dir() -> pathlib.Path | None:
    """Repaired qPCR exports for the G4 anchor. ``YSTWIN_QPCR`` overrides."""
    return _resolve("YSTWIN_QPCR", data_dir() / "qpcr_repaired")


def qpcr_raw_dir() -> pathlib.Path | None:
    """Bio-Rad Cq exports as the machine wrote them. ``YSTWIN_QPCR_RAW`` overrides.

    Distinct from :func:`qpcr_dir`, which holds only the one replicate whose archive
    entries were lowercased and had to be rebuilt by
    ``scripts/repair_qpcr_export.py``. The other replicates open as exported and are read
    from wherever the lab keeps them.
    """
    return _resolve(
        "YSTWIN_QPCR_RAW",
        IGEM_ROOT / "qpcr" / "ER and Oxidative Stress Biosensors",
        pathlib.Path.home()
        / "Desktop/Result & Analysis/rt-qPCR results/ER and Oxidative Stress Biosensors",
    )


def yeast_gem() -> pathlib.Path | None:
    """Yeast9 v9.0.2 SBML. ``YSTWIN_YEAST_GEM`` overrides.

    The gzipped copy in ``data/gem`` comes first and is tracked: 11 MB of SBML is 330 kB
    compressed, cobrapy reads ``.xml.gz`` directly, and depending on a sibling checkout
    made every FBA result in this repository unreproducible by anyone else.
    """
    return _resolve(
        "YSTWIN_YEAST_GEM",
        data_dir() / "gem" / "yeast-GEM.xml.gz",
        data_dir() / "gem" / "yeast-GEM.xml",
        IGEM_ROOT / "dente" / "data" / "yeast-GEM.xml",
        IGEM_ROOT / "orchid-engine" / "data" / "yeast-GEM.xml",
    )


def ec_yeast_gem() -> pathlib.Path | None:
    """GECKO enzyme-constrained model, on the SAME yeast-GEM as :func:`yeast_gem`.

    STAYS ON 8.3.4 DELIBERATELY, and the reason is measured rather than inherited.
    ``ecYeastGEM_yeast902.xml.gz`` is also vendored -- a GECKO 3.2.5 build on this repo's
    own yeast-GEM 9.0.2, made by ``scripts/gecko/build_ecyeastgem.m`` because no ecYeastGEM
    on any 9.x is published anywhere. It is BETTER on growth and WORSE where it matters
    more here. Against van Hoek 1998's aerobic batch phenotype:

        model    growth        glucose          ethanol
        8.3.4    0.377 (+5.8%) 17.93 (+61.5%)   29.57 (+113%)
        9.0.2    0.382 (+4.5%)  4.39 (-60.5%)    0.00 (-100%)

    The 9.0.2 build is FULLY RESPIRATORY: it reaches the same growth on 4.4 glucose and
    10.0 oxygen with no ethanol at all, where the 8.3.4 model ferments. Every culture this
    repository is calibrated against ferments -- Elizondo's six states secrete 2.9-16.1
    mmol/gDCW/h of ethanol at BOTH dilution rates -- so a model that cannot produce
    overflow is qualitatively wrong for them, and a percentage point of growth accuracy
    does not buy that back. Overflow in an ecModel emerges when the protein pool is tight
    enough to trade protein-expensive respiration for protein-cheap fermentation; the
    DLKcat-derived kcats in the 9.0.2 build are evidently generous enough that respiration
    stays affordable. Fixing that is a sigma/pool calibration exercise, not a rebuild.

    THE TWO MODELS USE DIFFERENT CONVENTIONS and both remain resolvable, so code that
    touches uptake must not hard-code either. GECKO 2 pins ``r_1714`` at ``(0, 0)`` and
    supplies through ``r_1714_REV``; GECKO 3 leaves ``r_1714`` unsplit. The protein pool
    is ``(0, 0.1037)`` grams in one and ``(-125, 0)`` milligrams in the other.
    :func:`ystwin.fba.physiology.cap_uptake` already reads the split from the model
    rather than trusting a reaction id, which is why this swap needed no code change.

    ``YSTWIN_EC_YEAST_GEM`` overrides, so trying the 9.0.2 build is one environment
    variable: ``YSTWIN_EC_YEAST_GEM=data/gem/ecYeastGEM_yeast902.xml.gz``.
    """
    return _resolve(
        "YSTWIN_EC_YEAST_GEM",
        data_dir() / "gem" / "ecYeastGEM_batch.xml.gz",
        data_dir() / "gem" / "ecYeastGEM_batch.xml",
        IGEM_ROOT / "dente" / "data" / "ecYeastGEM_batch.xml",
    )


def require(resolved: pathlib.Path | None, what: str, env_var: str) -> pathlib.Path:
    """Return ``resolved`` or explain precisely how to supply it.

    Args:
        resolved: The output of one of the resolvers above.
        what: Human name of the asset, for the message.
        env_var: The environment variable that overrides its location.

    Raises:
        FileNotFoundError: naming the variable to set. Callers that should skip
            rather than fail must check for ``None`` instead of calling this.
    """
    if resolved is None:
        raise FileNotFoundError(
            f"{what} not found. Set {env_var} to its location, or place it at the "
            f"documented default under {IGEM_ROOT}."
        )
    return resolved


def resolve_or_exit(resolved: pathlib.Path | None, what: str, env_var: str) -> pathlib.Path:
    """:func:`require`, but for a command-line entry point rather than a library caller.

    A traceback is a reasonable thing to hand a caller that can catch it and a useless
    thing to hand someone who has just cloned the repository and wants a table. Scripts
    get the same message on stderr and a non-zero exit status, with no stack.

    Args:
        resolved: The output of one of the resolvers above.
        what: Human name of the asset, for the message.
        env_var: The environment variable that overrides its location.

    Raises:
        SystemExit: naming the variable to set. Exit status 1.
    """
    try:
        return require(resolved, what, env_var)
    except FileNotFoundError as exc:
        raise SystemExit(str(exc)) from exc
