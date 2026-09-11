# DIAtoPRM

Shiny for Python app that generates a PRM (Parallel Reaction Monitoring)
transition list from an IPA gene list + a DIA-NN / timsDIA-NN search output.
Groups peptides per gene, computes m/z and isolation width, drops weak
charge-state duplicates and undigested fragments, and writes a
QE/timsTOF-style transition CSV.

Ported from `original_script/PRM_script_b0_6.py` (Sameer Varma & Stanley M.
Stevens Jr., USF) — same algorithm, minus the terminal `input()` prompts.

- `prm_logic.py` — core algorithm (pure functions, no UI dependency)
- `rt_calibration.py` — optional RT recalibration from a HeLa bridge run (see below)
- `app.py` — Shiny UI: file uploads, parameter inputs, results table, log, CSV download

## Run locally

```
uv run shiny run app.py
```

Then open the printed `http://127.0.0.1:8000` URL. Upload the IPA `.txt`
file (tab-delimited, 2 header rows to skip, needs a `Symbol` column) and the
MS search output as `.tsv` or `.parquet` (DIA-NN/timsDIA-NN report — both are
read into the same schema; parquet is preferred when available since it's
much smaller and faster to upload), set the MS file type (DIANN/TIMSDIANN)
and tolerances in the sidebar, and click **Run**.

## RT recalibration (optional)

If the MS file you're scheduling from was acquired on a different LC
column/gradient than the one you'll actually run PRM on, its recorded RT
values need to be translated onto the new time scale before they're a
trustworthy scheduling window. Open **RT recalibration** in the sidebar,
enable it, and upload two HeLa DIA-NN reports (`.tsv`/`.parquet`, same as
the MS file): one from around when the historic run was acquired (old
column/gradient) and one from the current setup (new column/gradient).
Peptides are matched between them by `Modified.Sequence` (median RT per
peptide), and a mapping is fit from old RT → new RT, then applied to the
main MS file's RT column before the rest of the pipeline runs.

Two fit types are available:
- **Linear** — good when the column swap is basically a uniform time
  shift/stretch (e.g. same gradient shape, just scaled).
- **Quadratic** — also captures non-uniform compression near the start/end
  of the gradient, which a *reshaped* gradient (not just a uniform
  time-rescale) commonly introduces.

Check the **RT Calibration** tab after running: it shows R², the number of
matched HeLa peptides, and a scatter plot of old vs. new RT with the fit
line. A low R² or a scatter that visibly bows away from the line means the
selected fit type isn't capturing the relationship well — try the other one.
