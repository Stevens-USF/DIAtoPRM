# DIAtoPRM

Shiny for Python app that generates a PRM (Parallel Reaction Monitoring)
transition list from an IPA gene list + a DIA-NN / timsDIA-NN search output.
Groups peptides per gene, computes m/z and isolation width, drops weak
charge-state duplicates and undigested fragments, and writes a
QE/timsTOF-style transition CSV.

Ported from `original_script/PRM_script_b0_6.py` (Sameer Varma & Stanley M.
Stevens Jr., USF) — same algorithm, minus the terminal `input()` prompts.

- `prm_logic.py` — core algorithm (pure functions, no UI dependency)
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
