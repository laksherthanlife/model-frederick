"""Stop a run from replacing a tracked table with one built from fewer inputs.

Two scripts were found doing the same thing on the same afternoon, and neither was a
mistake anyone would spot by reading it:

* ``make_splits.py`` without ``YSTWIN_PLATES`` wrote a simulated-only manifest over the
  tracked one, losing 714 of 1650 rows, printing "wrote" and exiting 0. The failure
  surfaced three scripts away as ``'DataFrame' object has no attribute 'plate'``.
* ``run_g4.py`` without ``YSTWIN_QPCR_RAW`` read one qPCR replicate instead of three, and
  the anchor for both oxidative constructs went from a measured effect to ``NaN`` with a
  replicate requirement of 0 -- which reads exactly like "no anchor exists" rather than
  "this machine could not see two of them".

The shape is the same in both. A script that resolves its inputs individually, skips the
ones it cannot find, and writes whatever it managed to build, is a script whose output
silently encodes the configuration of the machine that ran it. The tables look fine. They
are smaller, and nothing says so.

The rule here is deliberately narrow. Refuse only where writing would **destroy**
something -- a tracked artefact that exists and was built from more. A fresh clone with no
artefact at all can build whatever it can reach, because there is nothing to lose and
refusing would strand the first run.
"""

from __future__ import annotations

import pathlib

__all__ = ["refuse_partial_rebuild"]


def refuse_partial_rebuild(destination: str | pathlib.Path, *, found: int, expected: int,
                           what: str, remedy: str) -> None:
    """Raise if ``destination`` exists and this run reached fewer inputs than it should.

    Args:
        destination: The tracked artefact the caller is about to overwrite.
        found: How many input sources this run actually reached.
        expected: How many it would reach with everything configured.
        what: What the inputs are, for the message ("qPCR replicates").
        remedy: The environment variable or step that would supply the rest.

    Raises:
        SystemExit: when ``found < expected`` and ``destination`` already exists. Not
            ``ValueError``: the caller is a script and the right outcome is a non-zero
            exit with a message, not a traceback a user has to read past.
    """
    if found >= expected:
        return
    path = pathlib.Path(destination)
    missing = expected - found
    if not path.exists():
        # Nothing to destroy. Say what is thin so the artefact is not later mistaken for
        # a complete one, and let the run proceed.
        print(f"  NOTE: building {path.name} from {found} of {expected} {what}; "
              f"{remedy}")
        return
    raise SystemExit(
        f"{path} exists and this run reached only {found} of {expected} {what}, so it "
        f"would replace a complete table with one missing {missing}. The result would "
        f"look like a finding about the data rather than about this machine. {remedy}"
    )
