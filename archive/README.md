# Archive

Things that are no longer part of the argument, kept because deleting them would destroy
the record of why the current design is what it is.

**The rule for what belongs here.** A file is archived when it is *superseded* — something
else now does its job — or when it is *orphaned*, meaning nothing outside its own tests
would ever call it and no plan exists to. A file is **not** archived merely for being
untested, unused today, or parked pending a measurement: those are live code with a stated
reason, and they live in the main tree with that reason written down.

**Nothing here is imported by anything in `src/`, `scripts/` or `tests/`.** That is the
property that makes archiving different from moving clutter around, and
`tests/test_archive_is_inert.py` enforces it. If something here becomes needed again, move
it back rather than importing across the boundary.

**What is deliberately NOT here:**

- `docs/superseded/` — a different thing, and better where it is. Those files record results
  that were *published and then corrected*, and each sits next to the current result it
  replaced. Their value is being findable from the claim they overturned.
- `scripts/parked/` and everything in `docs/PARKED.md` — parked is not stale. Those are
  complete, tested, and waiting on a measurement rather than on a decision.
- Anything only reachable from tests but part of a coherent library layer — the
  thermodynamics stack (`bridge/thermodynamic.py`, `bridge/tmfa.py`,
  `bridge/equilibrator.py`) has no caller today and is a real layer with real tests. Being
  called by nothing yet is not the same as being superseded by something.

**Docs still cite these files, and should.** A finding does not stop being a finding
because the code that produced it was retired -- `docs/SOURCE_AUDIT.md` quotes
`archive/tests/test_deconvolve.py` for the conclusion that mOrange2 is the wrong partner
for YFP, which is a wet-lab decision that outlives the module. References were repointed
here rather than deleted. One set was not: `docs/CITATION_AUDIT.md` names `docs/CLEANUP.md`
at four line numbers and that path no longer exists. The rows are still right about the
file; it is now `archive/docs/CLEANUP.md`.

**An archived test is a frozen document, not a live assertion.** `tests/test_module_admission.py`
here builds its fixtures from `REFERENCE_AEROBIC_BATCH`'s old glucose 21.3 and ethanol 27.4 --
numbers `docs/superseded/reference-phenotype.md` withdrew as unsourced, and which the live
reference replaced with 11.1 and 13.9. Nothing collects the file, so it asserts nothing. It is
left as written because rewriting it against today's reference would document neither the
module it tested nor the day it was retired.

## Contents

| Path | Superseded or orphaned by | Why |
| --- | --- | --- |
| `code/deconvolve.py` | orphaned | Spectral unmixing for a four-channel plate. Nothing outside its tests ever called it, and the panel it was written for is not the panel that was run. |
| `code/module_admission.py` | orphaned | A gate for admitting a latent module, whose only caller was a one-line pass-through in `bridge/latent_bridge.py` that was itself deleted as dead. |
| `outputs/atp_sensor_*.png` | orphaned | Three committed figures no script produced. Kept because the panels are cited; archived because a committed artefact nobody can regenerate is not a result. **Half of that reason has since expired**: `analyse_atp_sensor.py --figures` was added in `8809344` and is the writer they never had, drawing the plate-derived panels and reproducing their archived values to the digit (two more are plate data and still not redrawn -- one folded into another panel, one printed as a number). They stay here because the rest of each render -- a literature bar on a broken citation, a model comparison whose constant has moved, an analytic curve in a free parameter -- still cannot be regenerated. `docs/research/XPT_INVENTORY.md` has the split panel by panel. |
| `docs/CLEANUP.md` | superseded | An 865-line running log of a 2026-08 cleanup, restated in `AUDIT_2026_08.md`. Historical, and **not uniformly correct**: `docs/CITATION_AUDIT.md` rows 2 and 3 record two of its own resolution entries as wrong, one of them marking `PARAMETER_SOURCES.md:240` RESOLVED when the assertion is still there verbatim today. Archived rather than deleted for exactly that reason -- a resolution log that overstated itself is worth being able to read. |
