import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd
from rt_calibration import fit_rt_calibration, apply_rt_calibration

rng = np.random.default_rng(0)
n = 500
seqs = [f"PEPT{i}IDER" for i in range(n)]
old_rt = np.sort(rng.uniform(2, 58, n))

print("=== Case 1: clean linear shift+stretch (uniform gradient time-rescale) ===")
true_slope, true_intercept = 1.5, 2.0
new_rt_linear = true_slope * old_rt + true_intercept + rng.normal(0, 0.05, n)

df_old = pd.DataFrame({"Modified.Sequence": seqs, "RT": old_rt})
df_new = pd.DataFrame({"Modified.Sequence": seqs, "RT": new_rt_linear})

cal_lin = fit_rt_calibration(df_old, df_new, degree=1)
print("fit:", cal_lin.describe(), "R2=", round(cal_lin.r_squared, 5), "n=", cal_lin.n_points)
assert abs(cal_lin.coeffs[0] - true_slope) < 0.01
assert abs(cal_lin.coeffs[1] - true_intercept) < 0.2
assert cal_lin.r_squared > 0.99

print("\n=== Case 2: reshaped gradient -> curved relationship (compression near ends) ===")
# simulate compression at both ends of the gradient (common with a reshaped/different
# length gradient, not just a uniform stretch)
new_rt_curved = 1.2 * old_rt + 0.02 * (old_rt - 30) ** 2 / 10 + rng.normal(0, 0.05, n)
df_new2 = pd.DataFrame({"Modified.Sequence": seqs, "RT": new_rt_curved})

cal_lin2 = fit_rt_calibration(df_old, df_new2, degree=1)
cal_quad2 = fit_rt_calibration(df_old, df_new2, degree=2)
print("linear fit R2=", round(cal_lin2.r_squared, 5))
print("quadratic fit R2=", round(cal_quad2.r_squared, 5))
assert cal_quad2.r_squared > cal_lin2.r_squared, "quadratic should fit the curved case better than linear"

print("\n=== apply_rt_calibration on a mock MS dataframe ===")
df_ms = pd.DataFrame({"RT": [10.0, 20.0, 30.0]})
df_ms_cal = apply_rt_calibration(df_ms, cal_lin)
print(df_ms_cal)
expected = true_slope * df_ms["RT"].to_numpy() + true_intercept
assert np.allclose(df_ms_cal["RT"].to_numpy(), expected, atol=0.05)

print("\n=== too few matched points raises a clear error ===")
try:
    fit_rt_calibration(df_old.head(3), df_new.head(3))
    raise SystemExit("expected ValueError for too few points")
except ValueError as e:
    print("Correctly raised:", e)

print("\nAll RT calibration tests passed.")
