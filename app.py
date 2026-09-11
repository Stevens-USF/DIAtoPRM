"""
DIAtoPRM — Shiny for Python front end for the PRM transition-list generator.

Wraps prm_logic.run_prm() (ported from PRM_script_b0_6.py) with file uploads,
parameter inputs, a results table, a log panel, and a CSV download — replacing
the original script's input() prompts.
"""

from pathlib import Path

import pandas as pd
from shiny import App, reactive, render, ui

from prm_logic import OUTPUT_HEADER, PrmParams, run_prm


def read_ms_file(path: str) -> pd.DataFrame:
    """DIA-NN/timsDIA-NN reports ship as either .tsv or .parquet with the same schema."""
    if Path(path).suffix.lower() == ".parquet":
        return pd.read_parquet(path)
    return pd.read_csv(path, delimiter="\t")


app_ui = ui.page_sidebar(
    ui.sidebar(
        ui.input_file("ipa_file", "IPA gene list (.txt, tab-delimited)", accept=[".txt"]),
        ui.input_file("ms_file", "MS search output (.tsv or .parquet)", accept=[".tsv", ".txt", ".parquet"]),
        ui.input_select("filetype_ms", "MS file type", choices=["DIANN", "TIMSDIANN"]),
        ui.input_numeric("rt_range", "RT range (seconds)", value=240),
        ui.input_numeric("im_tolerance", "IM tolerance (%)", value=4),
        ui.input_numeric("charge_tol", "Charge-state relative intensity tolerance (fraction)", value=0.2, step=0.05),
        ui.input_numeric("digestion_tol", "Undigested-fragment relative intensity tolerance (fraction)", value=0.1, step=0.05),
        ui.input_checkbox("top_n_all", "Output all sequences (no top-N cutoff)", value=True),
        ui.input_numeric("top_n", "Top N highest-intensity sequences", value=5),
        ui.input_checkbox("drop_m", "Drop sequences containing Methionine", value=False),
        ui.input_checkbox("split_genes", "Split multi-gene 'Genes' entries", value=False),
        ui.input_action_button("run", "Run", class_="btn-primary"),
        width=380,
    ),
    ui.navset_tab(
        ui.nav_panel(
            "Results",
            ui.download_button("download_csv", "Download CSV"),
            ui.output_data_frame("results_table"),
        ),
        ui.nav_panel("Log", ui.output_text_verbatim("log_output")),
    ),
    title="DIAtoPRM — PRM Transition List Generator",
)


def server(input, output, session):
    result = reactive.Value(None)

    @reactive.effect
    @reactive.event(input.run)
    def _run():
        ipa_info = input.ipa_file()
        ms_info = input.ms_file()
        if not ipa_info or not ms_info:
            ui.notification_show("Please upload both the IPA file and the MS file.", type="error")
            return
        try:
            df_ipa = pd.read_csv(ipa_info[0]["datapath"], delimiter="\t", skiprows=2)
            df_ms = read_ms_file(ms_info[0]["datapath"])
            params = PrmParams(
                filetype_ms=input.filetype_ms(),
                split_genes=input.split_genes(),
                intensity_charge_state_tolerance=input.charge_tol(),
                intensity_digestion_tolerance=input.digestion_tol(),
                rt_range=input.rt_range(),
                im_tolerance=input.im_tolerance(),
                top_n=-1 if input.top_n_all() else int(input.top_n()),
                drop_sequence_aa_m=input.drop_m(),
            )
            res = run_prm(df_ipa, df_ms, params)
            result.set(res)
            ui.notification_show(f"Done — {len(res.output_df)} transitions.", type="message")
        except Exception as e:
            ui.notification_show(f"Error: {e}", type="error", duration=None)
            raise

    @render.data_frame
    def results_table():
        res = result.get()
        df = res.output_df.copy() if res is not None else pd.DataFrame(columns=OUTPUT_HEADER)
        if res is not None:
            df.columns = OUTPUT_HEADER
        return df

    @render.text
    def log_output():
        res = result.get()
        return "\n".join(res.log) if res is not None else "Run the analysis to see the log."

    @render.download(filename="output.csv")
    def download_csv():
        res = result.get()
        if res is None:
            return
        df = res.output_df.copy()
        df.columns = OUTPUT_HEADER
        yield df.to_csv(index=False)


app = App(app_ui, server)
