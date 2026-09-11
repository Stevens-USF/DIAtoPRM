"""
RT recalibration using a HeLa bridge run.

If a historic MS run was acquired on a different LC column/gradient than the
one that will actually be used for PRM, its recorded RT values need to be
translated onto the new column's time scale before they're usable as a
scheduling window. HeLa QC runs on both setups give thousands of matched
peptides spread across the whole gradient, which is enough to fit that
mapping and check whether it's linear or curved (gradient-length changes in
particular often compress retention non-uniformly near the start/end of the
run, not just uniformly stretch it).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

MIN_CALIBRATION_POINTS = 10


@dataclass
class RtCalibration:
    degree: int
    coeffs: np.ndarray  # np.polyfit output, highest power first
    r_squared: float
    n_points: int
    old_rt: np.ndarray  # matched HeLa anchor points, for diagnostics/plotting
    new_rt: np.ndarray

    def predict(self, rt):
        return np.polyval(self.coeffs, rt)

    def describe(self) -> str:
        if self.degree == 1:
            slope, intercept = self.coeffs
            return f"new_RT = {slope:.4f} * old_RT + {intercept:.4f}"
        terms = " + ".join(f"{c:.5g}*old_RT^{self.degree - i}" for i, c in enumerate(self.coeffs[:-1]))
        return f"new_RT = {terms} + {self.coeffs[-1]:.4f}"


def _median_rt_by_sequence(df_hela: pd.DataFrame) -> pd.Series:
    return df_hela.groupby("Modified.Sequence")["RT"].median()


def fit_rt_calibration(df_hela_old: pd.DataFrame, df_hela_new: pd.DataFrame, degree: int = 1) -> RtCalibration:
    """Fit an old-column-RT -> new-column-RT mapping from two HeLa reports, matched
    on Modified.Sequence. degree=1 is a linear (shift+stretch) fit; degree=2 also
    captures non-uniform compression near the ends of the gradient, which a
    gradient-length change often introduces."""
    old_by_seq = _median_rt_by_sequence(df_hela_old)
    new_by_seq = _median_rt_by_sequence(df_hela_new)
    matched = pd.concat([old_by_seq.rename("old"), new_by_seq.rename("new")], axis=1).dropna()

    if len(matched) < MIN_CALIBRATION_POINTS:
        raise ValueError(
            f"Only {len(matched)} peptides matched between the two HeLa runs (need at least "
            f"{MIN_CALIBRATION_POINTS}) -- too few to fit a reliable RT calibration. Check that "
            "both files are HeLa DIA-NN reports with overlapping identifications."
        )

    old_rt = matched["old"].to_numpy()
    new_rt = matched["new"].to_numpy()

    coeffs = np.polyfit(old_rt, new_rt, degree)
    predicted = np.polyval(coeffs, old_rt)
    residuals = new_rt - predicted
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((new_rt - new_rt.mean()) ** 2)
    r_squared = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")

    return RtCalibration(
        degree=degree, coeffs=coeffs, r_squared=r_squared, n_points=len(matched), old_rt=old_rt, new_rt=new_rt
    )


def apply_rt_calibration(df_ms: pd.DataFrame, calibration: RtCalibration) -> pd.DataFrame:
    df_ms = df_ms.copy()
    df_ms["RT"] = calibration.predict(df_ms["RT"].to_numpy())
    return df_ms
