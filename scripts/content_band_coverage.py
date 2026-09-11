"""Measure conditional channel and joint coverage on outer held-out strains.

    python3 scripts/content_band_coverage.py --draws 500 --output-dir /tmp/content-coverage

The previous table compared simulated widths and called the Monte Carlo width truth.
It did not measure coverage of an observed outcome. This table instead retains all
six condition predictions from three independent strain groups and scores both
beta-carotene and lycopene, including failed predictions or unavailable intervals.

Gene and entry/branch model selection occur inside each outer training fold. Two-stage
paired-strain resampling then refits entry and kinetic parameters jointly, using the
released expression and product observation intervals. Future assay precision comes
only from training intervals, never from held-out products. Marginal and simultaneous
strain coverage have separate nominal targets and actual numerators/denominators.
These are conditional bands, not proof that three strains calibrate 95% coverage.
"""

from __future__ import annotations

import argparse
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

import pandas as pd

from ystwin import paths
from ystwin.pathway.flux import (
    INTERVAL_ASSUMPTIONS,
    carotenoid_measurements,
    score_product_validation,
    summarize_product_validation,
)
from ystwin.pathway.spec import load_pathway


def coverage(draws: int = 500, *, seed: int = 0, nominal_coverage: float = 0.95,
             mode: str = "nested") -> pd.DataFrame:
    """Actual outer-held-out coverage, conditional on the declared resampling assumptions."""
    # Seeded resampling makes the conditional intervals reproducible, not calibrated.
    return score_product_validation(
        load_pathway("beta_carotene"), carotenoid_measurements(), mode=mode,
        draws=draws, seed=seed, nominal_coverage=nominal_coverage)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--draws", type=int, default=500)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--nominal-coverage", type=float, default=0.95)
    parser.add_argument("--mode", choices=("nested", "fixed_crte_comparison"), default="nested")
    parser.add_argument("--output-dir", type=pathlib.Path)
    args = parser.parse_args()

    frame = coverage(args.draws, seed=args.seed, nominal_coverage=args.nominal_coverage,
                     mode=args.mode)
    # Report nominal and actual coverage on the same conditions, not simulated widths.
    # A missing band or failed prediction stays in each channel's denominator.
    # The joint rectangle covers both channels and all conditions of one held-out strain.
    # Six predictions still arise from only three independent strain groups.
    summary = summarize_product_validation(frame, comparisons=False)
    columns = ["state", "strain", "selected_model", "status", "interval_status",
               "product_inside_band", "lycopene_inside_band", "joint_inside_band",
               "nominal_channel_coverage", "nominal_strain_joint_coverage",
               "interval_successful_draws", "interval_failed_draws"]
    print(frame[columns].to_string(index=False))
    print("\n" + summary.to_string(index=False))
    print("\n" + INTERVAL_ASSUMPTIONS)
    # Kinetic parameters and entry scale are refitted on the SAME bootstrap draws.
    # Source measurement errors are conditionally independent only by assumption;
    # neither their missing covariance nor new-strain model discrepancy is identified.
    directory = args.output_dir if args.output_dir is not None else paths.outputs_dir()
    directory.mkdir(parents=True, exist_ok=True)
    frame.to_csv(directory / "content_band_coverage.csv", index=False)
    summary.to_csv(directory / "content_band_coverage_summary.csv", index=False)
    # Same selection scores as predict_product's; see outputs/product_prediction_selection.csv.
    print(f"\n-> {directory}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
