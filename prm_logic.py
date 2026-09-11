"""
Core PRM transition-list generation logic.

Ported from PRM_script_b0_6.py (Sameer Varma & Stanley M. Stevens Jr., USF),
with input() prompts and print-to-stdout removed in favor of a log list and
a returned DataFrame, so it can be driven from a UI instead of a terminal.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

SUPPORTED_MS_FILETYPES = ["DIANN", "TIMSDIANN"]

# monoisotopic residue masses, https://education.expasy.org/student_projects/isotopident/htdocs/aa-list.html
AA_MASSES = {
    "A": 71.03711, "L": 113.08406, "V": 99.06841, "I": 113.08406, "Y": 163.06333,
    "P": 97.05276, "F": 147.06841, "W": 186.07931, "D": 115.02694, "E": 129.04259,
    "R": 156.10111, "K": 128.09496, "N": 114.04293, "Q": 128.05858, "C": 103.00919,
    "S": 87.03203, "T": 101.04768, "G": 57.02146, "M": 131.04049, "H": 137.05891,
}
MASS_HYDROGEN = 1.00727647
MASS_WATER = 18.01057

# https://unimod.org
UNIMOD_MASSES = {1: 42.010565, 4: 57.021464, 7: 0.984016, 35: 15.994915, 21: 79.966331}

OUTPUT_COLUMNS = [
    "Mass_over_Charge", "Precursor.Charge", "Isolation_width", "RT_median",
    "RT_range", "Start_IM", "End_IM", "CE", "External_ID", "Description",
]
OUTPUT_HEADER = [
    "Mass [m/z]", "Charge", "Isolation Width [m/z]", "RT [s]", "RT Range [s]",
    "Start IM [1/K0]", "End IM [1/K0]", "CE [eV]", "External ID", "Description",
]


def assign_mass(sequence: str, charge: int) -> float:
    mass = sum(AA_MASSES[c] for c in sequence)
    mass += charge * MASS_HYDROGEN
    mass += MASS_WATER
    return mass


def adjust_unimod(modified_sequence: str) -> float:
    mods = [int(x) for x in re.findall(r"\d+", modified_sequence)]
    return sum(UNIMOD_MASSES[m] for m in mods)


def number_kr(sequence: str) -> int:
    return sequence.count("R") + sequence.count("K")


@dataclass
class PrmParams:
    filetype_ms: str = "DIANN"
    split_genes: bool = False
    isolation_width_cutoff: float = 700
    isolation_width_low: float = 2
    isolation_width_high: float = 3
    intensity_charge_state_tolerance: float = 0.2
    intensity_digestion_tolerance: float = 0.1
    min_fragment_length: int = 0
    rt_range: float = 240
    im_tolerance: float = 4
    top_n: int = -1
    drop_sequence_aa_m: bool = False
    decimals: int = 5


@dataclass
class PrmResult:
    output_df: pd.DataFrame
    log: list[str] = field(default_factory=list)


def run_prm(df_ipa: pd.DataFrame, df_ms: pd.DataFrame, params: PrmParams) -> PrmResult:
    log: list[str] = []

    def say(*parts) -> None:
        log.append(" ".join(str(p) for p in parts))

    df_ipa = df_ipa.copy()
    df_ms = df_ms.copy()

    say("Genes in IPA file:", len(df_ipa))
    say("Rows in MS file:", len(df_ms))

    if params.split_genes:
        df_ms["Genes"] = df_ms["Genes"].str.split(";")
        df_ms = df_ms.explode("Genes").reset_index(drop=True)
        cols = list(df_ms.columns)
        cols.append(cols.pop(cols.index("File.Name")))
        df_ms = df_ms[cols]
        say("Adding rows in MS file to accommodate multiple 'Genes' entries:", len(df_ms))

    df_ipa["Symbol"] = df_ipa["Symbol"].str.lower()
    df_ms["Genes"] = df_ms["Genes"].str.lower()

    im_col = "IM" if params.filetype_ms.upper() == "DIANN" else "Exp.1/K0"

    output_rows: list[pd.DataFrame] = []

    say("\nGrouping...")
    for g_ipa in df_ipa["Symbol"]:
        df_ms_select = df_ms.loc[df_ms["Genes"] == g_ipa]
        say("-" * 50)
        say("Searching for Gene", g_ipa, ": found", len(df_ms_select), "entries")
        say("-" * 50)

        if len(df_ms_select) == 0:
            say("---------------WARNING: ZERO ENTRIES FOR GENE", g_ipa, "---------------")
            continue

        grouped = df_ms_select.groupby(["Modified.Sequence", "Precursor.Charge"]).agg(
            {"RT": ["median"], im_col: ["median"], "Precursor.Normalised": ["mean"], "Stripped.Sequence": ["first"]}
        )
        grouped.columns = ["RT_median", "IM_median", "Intensity_mean", "Stripped.Sequence"]
        grouped = grouped.reset_index()
        grouped["Gene"] = g_ipa
        grouped["RT_median"] = np.round(60 * grouped["RT_median"], decimals=params.decimals)
        grouped["RT_range"] = params.rt_range
        grouped["Start_IM"] = np.round(grouped["IM_median"] * (1 - params.im_tolerance / 100), decimals=params.decimals)
        grouped["End_IM"] = np.round(grouped["IM_median"] * (1 + params.im_tolerance / 100), decimals=params.decimals)
        grouped["Intensity_mean"] = np.round(grouped["Intensity_mean"], decimals=params.decimals)

        for i in grouped.index:
            grouped.at[i, "Description"] = grouped.at[i, "Gene"] + "_" + grouped.at[i, "Modified.Sequence"]
            mass = assign_mass(grouped.at[i, "Stripped.Sequence"], grouped.at[i, "Precursor.Charge"])
            mass += adjust_unimod(grouped.at[i, "Modified.Sequence"])
            grouped.at[i, "Mass_over_Charge"] = mass / grouped.at[i, "Precursor.Charge"]
            grouped.at[i, "Isolation_width"] = (
                params.isolation_width_low
                if grouped.at[i, "Mass_over_Charge"] < params.isolation_width_cutoff
                else params.isolation_width_high
            )
        grouped["Mass_over_Charge"] = np.round(grouped["Mass_over_Charge"], decimals=params.decimals)

        grouped["CE"] = ""
        grouped["External_ID"] = ""

        # drop weak charge-state duplicates
        say("Examining sequences with different charges")
        dropped: set[int] = set()
        for i in grouped.index:
            for j in grouped.index:
                if grouped.at[j, "Stripped.Sequence"] == grouped.at[i, "Stripped.Sequence"] and j != i:
                    intensity_sum = grouped.at[i, "Intensity_mean"] + grouped.at[j, "Intensity_mean"]
                    intensity_j_rel = grouped.at[j, "Intensity_mean"] / intensity_sum
                    if intensity_j_rel < params.intensity_charge_state_tolerance:
                        say(
                            g_ipa, ": Removing sequence", grouped.at[j, "Stripped.Sequence"],
                            "with charge", grouped.at[j, "Precursor.Charge"], "as its Intensity=",
                            grouped.at[j, "Intensity_mean"], "is less than",
                            params.intensity_charge_state_tolerance, "of total Intensity of", intensity_sum,
                        )
                        dropped.add(j)
        grouped = grouped.drop(labels=list(dropped))
        say("Total sequences dropped =", len(dropped))

        # drop incompletely-digested fragments
        say("Examining incomplete digestion")
        dropped = set()
        for i in grouped.index:
            if number_kr(grouped.at[i, "Stripped.Sequence"]) > 1:
                grouped.at[i, "Digestion"] = "Incomplete"
                for j in grouped.index:
                    if i == j:
                        continue
                    seq_i, seq_j = grouped.at[i, "Stripped.Sequence"], grouped.at[j, "Stripped.Sequence"]
                    if seq_j in seq_i and seq_i != seq_j and len(seq_i) - len(seq_j) > params.min_fragment_length:
                        if grouped.at[i, "Intensity_mean"] < params.intensity_digestion_tolerance * grouped.at[j, "Intensity_mean"]:
                            dropped.add(i)
                            say(
                                "Dropping undigested fragment", seq_i, "as its Intensity=",
                                grouped.at[i, "Intensity_mean"], "is", params.intensity_digestion_tolerance,
                                "of the Intensity of", grouped.at[j, "Intensity_mean"], "of its digested fragment", seq_j,
                            )
                        else:
                            dropped.add(i)
                            dropped.add(j)
                            say(
                                "Dropping undigested fragment", seq_i, "and its digested fragment", seq_j,
                                "as their Intensities of", grouped.at[i, "Intensity_mean"], "and",
                                grouped.at[j, "Intensity_mean"], "are within", params.intensity_digestion_tolerance, "of each other.",
                            )
            else:
                grouped.at[i, "Digestion"] = "Complete"
        grouped = grouped.drop(labels=list(dropped))
        say("Total sequences dropped =", len(dropped))

        if params.drop_sequence_aa_m:
            before = len(grouped)
            grouped = grouped[~grouped["Stripped.Sequence"].str.contains("M")]
            say("Sequences containing methionine dropped =", before - len(grouped))

        grouped = grouped.sort_values(by=["Intensity_mean"], ascending=False)
        if params.top_n != -1:
            grouped = grouped.head(params.top_n)

        output_rows.append(grouped[OUTPUT_COLUMNS])

    if output_rows:
        output_df = pd.concat(output_rows, ignore_index=True)
    else:
        output_df = pd.DataFrame(columns=OUTPUT_COLUMNS)

    say("\nNORMAL TERMINATION: processed", len(df_ipa), "genes,", len(output_df), "transitions written")

    return PrmResult(output_df=output_df, log=log)
