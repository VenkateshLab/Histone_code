#!/usr/bin/env python3
"""Standalone Figure 1 generator for the GitHub package."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.colors import ListedColormap
from matplotlib.patches import Rectangle


SCRIPT_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_INPUT = SCRIPT_ROOT / "input"
DEFAULT_OUTPUT = SCRIPT_ROOT / "outputs" / "figure1"

BG = "#FFFDF8"
TEXT = "#000000"

STAGE_ORDER_A = ["E 10.5", "E 11.5", "E 12.5", "E 13.5", "E 14.5", "E 15.5", "E 16.5", "P0", "ADULT"]
ASSAY_ORDER = [
    "H3K4me1",
    "H3K4me2",
    "H3K4me3",
    "H3K27ac",
    "H3K9ac",
    "H3K27me3",
    "H3K9me3",
    "H3K36me3",
    "Bulk ATAC Seq",
    "Poly A RNA Seq",
]
DISPLAY_ASSAY_LABELS = {
    "H3K4me1": "H3K4me1",
    "H3K4me2": "H3K4me2",
    "H3K4me3": "H3K4me3",
    "H3K27ac": "H3K27ac",
    "H3K9ac": "H3K9ac",
    "H3K27me3": "H3K27me3",
    "H3K9me3": "H3K9me3",
    "H3K36me3": "H3K36me3",
    "Bulk ATAC Seq": "ATAC",
    "Poly A RNA Seq": "RNA",
}
COLORS_A = {
    "histone": "#355C7D",
    "chromatin": "#C06C84",
    "rna": "#6C9A8B",
    "missing": "#ECE7E1",
    "grid": "#C8C1B8",
    "legend_edge": "#6B6259",
    "text_light": "#F9F7F4",
    "text_dark": "#342E29",
    "bg": BG,
}

COLORS_B = {
    1: "#C93A3A",
    2: "#E8C547",
    3: "#D97A7A",
    4: "#F2B4B4",
    5: "#E67E22",
    6: "#F3A65A",
    7: "#F7C996",
    8: "#4A90E2",
    9: "#2E7D32",
    10: "#8BCF7A",
    11: "#B39DDB",
    12: "#9DD9F3",
    13: "#7B4FA3",
}
GRID_B = "#D7D0C7"
TEXT_B = "#1F1A17"

STAGE_ORDER = ["E10", "E11", "E12", "E13", "E14", "E15", "E16", "P0"]
STAGE_COLORS = {
    "E10": "#B33C3C",
    "E11": "#D66B3D",
    "E12": "#E6A141",
    "E13": "#C9B458",
    "E14": "#7AA974",
    "E15": "#4F9A8D",
    "E16": "#4D7FB8",
    "P0": "#6A4C93",
}
STATE_COLORS = sns.color_palette("tab20", 18)
STATE_CMAP = ListedColormap(STATE_COLORS)

matplotlib.rcParams["font.family"] = "DejaVu Sans"
matplotlib.rcParams["font.size"] = 6
matplotlib.rcParams["axes.linewidth"] = 0.6


def build_paths(output_dir: Path) -> dict[str, Path]:
    return {
        "a": output_dir / "Figure_1_A_sequencing_availability.png",
        "a_legend": output_dir / "Figure_1_A_legend.png",
        "b": output_dir / "Figure_1_B_states.png",
        "c": output_dir / "Figure_1_C.png",
        "d": output_dir / "Figure_1_D.png",
        "cd_legend": output_dir / "Figure_1_CD_legend.png",
        "e_e10": output_dir / "Figure_1_E_E10_5.png",
        "e_p0": output_dir / "Figure_1_E_P0.png",
        "e_legend": output_dir / "Figure_1_E_legend.png",
        "f_emissions": output_dir / "Figure_1_F_emissions.png",
        "f_jaccard": output_dir / "Figure_1_F_jaccard.png",
    }


def clean_value(value: object) -> str:
    if pd.isna(value):
        return ""
    text = str(value).replace("\xa0", " ").strip()
    if text in {"---", "-", ""}:
        return ""
    if text == "P":
        return "Present"
    return text


def stage_label(value: str) -> str:
    return value.replace(" ", "")


def load_table_a(csv_path: Path) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    first_col = df.columns[0]
    df = df.rename(columns={first_col: "Stage"})
    df["Stage"] = df["Stage"].map(clean_value)
    for column in ASSAY_ORDER:
        df[column] = df[column].map(clean_value)
    df["Stage"] = pd.Categorical(df["Stage"], categories=STAGE_ORDER_A, ordered=True)
    return df.sort_values("Stage").reset_index(drop=True)


def assay_kind(assay: str) -> str:
    if assay == "Bulk ATAC Seq":
        return "chromatin"
    if assay == "Poly A RNA Seq":
        return "rna"
    return "histone"


def draw_plot_a(df: pd.DataFrame, output_path: Path, width: float, height: float, dpi: int) -> None:
    n_rows = len(df)
    n_cols = len(ASSAY_ORDER)
    fig, ax = plt.subplots(figsize=(width, height), facecolor=BG)
    ax.set_facecolor(BG)

    for row_idx, (_, row) in enumerate(df.iterrows()):
        for col_idx, assay in enumerate(ASSAY_ORDER):
            value = row[assay]
            kind = assay_kind(assay)
            facecolor = COLORS_A[kind] if value else COLORS_A["missing"]
            rect = Rectangle(
                (col_idx, n_rows - row_idx - 1),
                1,
                1,
                facecolor=facecolor,
                edgecolor=COLORS_A["grid"],
                linewidth=0.45,
            )
            ax.add_patch(rect)

    ax.set_xlim(0, n_cols)
    ax.set_ylim(0, n_rows)
    ax.set_xticks([idx + 0.5 for idx in range(n_cols)])
    ax.set_xticklabels(
        [DISPLAY_ASSAY_LABELS[assay] for assay in ASSAY_ORDER],
        rotation=45,
        ha="right",
        rotation_mode="anchor",
        fontsize=5.6,
        weight="bold",
    )
    ax.set_yticks([idx + 0.5 for idx in range(n_rows)])
    ax.set_yticklabels([stage_label(v) for v in reversed(df["Stage"].astype(str))], fontsize=5.8, weight="bold")
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)

    fig.subplots_adjust(left=0.14, right=0.995, top=0.995, bottom=0.37)
    fig.savefig(output_path, dpi=dpi, facecolor=BG, pad_inches=0.02)
    plt.close(fig)


def draw_legend_a(output_path: Path, dpi: int) -> None:
    legend_items = [
        ("Histone", COLORS_A["histone"]),
        ("ATAC", COLORS_A["chromatin"]),
        ("RNA", COLORS_A["rna"]),
        ("N/A", COLORS_A["missing"]),
    ]
    fig, ax = plt.subplots(figsize=(3.0, 0.42), facecolor=BG)
    ax.set_facecolor(BG)
    ax.axis("off")
    x_positions = [0.02, 0.28, 0.50, 0.72]
    for (label, color), x in zip(legend_items, x_positions):
        ax.add_patch(
            Rectangle(
                (x, 0.3),
                0.03,
                0.32,
                transform=ax.transAxes,
                facecolor=color,
                edgecolor=COLORS_A["legend_edge"],
                linewidth=0.7,
            )
        )
        ax.text(x + 0.04, 0.46, label, transform=ax.transAxes, va="center", ha="left", fontsize=10, color=COLORS_A["text_dark"])
    fig.savefig(output_path, dpi=dpi, facecolor=BG, bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)


def load_states_b(csv_path: Path) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    df["STATES"] = df["STATES"].astype(int)
    return df.sort_values("STATES").reset_index(drop=True)


def draw_figure_b(df: pd.DataFrame, output_path: Path, width: float, height: float, dpi: int) -> None:
    n_rows = len(df)
    fig, ax = plt.subplots(figsize=(width, height), facecolor=BG)
    ax.set_facecolor(BG)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, n_rows + 0.8)
    ax.axis("off")

    x0, x1, x2, x3, x4 = 0.015, 0.09, 0.26, 0.61, 0.985
    c_state = (x0 + x1) / 2
    c_mn = (x1 + x2) / 2
    c_desc = (x2 + x3) / 2
    c_marks = (x3 + x4) / 2

    header_y = n_rows + 0.35
    ax.text(c_state, header_y, "State", ha="center", va="center", fontsize=5.2, weight="bold", color=TEXT_B)
    ax.text(c_mn, header_y, "Mnemonic", ha="center", va="center", fontsize=5.2, weight="bold", color=TEXT_B)
    ax.text(c_desc, header_y, "Description", ha="center", va="center", fontsize=5.2, weight="bold", color=TEXT_B)
    ax.text(c_marks, header_y, "Histone Mark", ha="center", va="center", fontsize=5.2, weight="bold", color=TEXT_B)

    for idx, row in df.iterrows():
        y = n_rows - idx - 0.02
        ax.add_patch(
            Rectangle((x0, y - 0.9), x4 - x0, 0.9, facecolor=COLORS_B[int(row["STATES"])], edgecolor=GRID_B, linewidth=0.5)
        )
        ax.text(c_state, y - 0.45, str(int(row["STATES"])), ha="center", va="center", fontsize=4.8, color=TEXT_B)
        ax.text(c_mn, y - 0.45, str(row["MNEMONICS"]), ha="center", va="center", fontsize=4.7, color=TEXT_B)
        ax.text(c_desc, y - 0.45, str(row["DESCRIPTION"]), ha="center", va="center", fontsize=4.4, color=TEXT_B)
        ax.text(c_marks, y - 0.45, str(row["Histone Mark Associated"]), ha="center", va="center", fontsize=3.9, color=TEXT_B)

    for x in [x0, x1, x2, x3, x4]:
        ax.plot([x, x], [0.08, n_rows + 0.02], color=GRID_B, linewidth=0.6, zorder=2)

    fig.subplots_adjust(left=0.02, right=0.995, top=0.98, bottom=0.04)
    fig.savefig(output_path, dpi=dpi, facecolor=BG)
    plt.close(fig)


def load_segment_tables(input_dir: Path) -> pd.DataFrame:
    frames = []
    for stage in STAGE_ORDER:
        path = input_dir / f"{stage}_18_segments.bed"
        df = pd.read_csv(path, sep="\t", header=None, names=["chr", "start", "end", "state"])
        df["stage"] = stage
        df["length"] = df["end"] - df["start"]
        df["state_num"] = df["state"].str.replace("E", "", regex=False).astype(int)
        frames.append(df)
    return pd.concat(frames, ignore_index=True)


def plot_fig_1c(all_df: pd.DataFrame, output_path: Path) -> None:
    count_df = all_df.groupby(["state_num", "stage"]).size().reset_index(name="segments")
    pivot = count_df.pivot(index="state_num", columns="stage", values="segments").reindex(index=range(1, 19), columns=STAGE_ORDER).fillna(0)

    fig, ax = plt.subplots(figsize=(3.0, 2.0), facecolor=BG)
    ax.set_facecolor(BG)
    x = list(range(1, 19))
    width = 0.1
    offsets = [-0.35, -0.25, -0.15, -0.05, 0.05, 0.15, 0.25, 0.35]
    for stage, offset in zip(STAGE_ORDER, offsets):
        ax.bar([i + offset for i in x], pivot[stage].values / 1000.0, width=width, color=STAGE_COLORS[stage], edgecolor="none", linewidth=0.0)

    ax.set_xlim(0.45, 18.55)
    ax.set_xticks(x)
    ax.set_xlabel("ChromHMM state", fontsize=7.2, color=TEXT)
    ax.set_ylabel("Number of segment count (x1000)", fontsize=7.2, color=TEXT)
    ax.tick_params(axis="x", labelsize=5.9, colors=TEXT, width=0.8, length=4)
    ax.tick_params(axis="y", labelsize=6.3, colors=TEXT, width=0.8, length=4)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.spines["left"].set_color(TEXT)
    ax.spines["bottom"].set_color(TEXT)
    ax.spines["left"].set_linewidth(0.8)
    ax.spines["bottom"].set_linewidth(0.8)
    fig.subplots_adjust(left=0.18, right=0.985, top=0.96, bottom=0.18)
    fig.savefig(output_path, dpi=600, facecolor=BG)
    plt.close(fig)


def plot_fig_1d(all_df: pd.DataFrame, output_path: Path) -> None:
    state_pct = all_df.groupby(["stage", "state_num"])["length"].sum().reset_index()
    totals = state_pct.groupby("stage")["length"].transform("sum")
    state_pct["percent"] = state_pct["length"] / totals * 100
    pivot = state_pct.pivot(index="state_num", columns="stage", values="percent").reindex(index=range(1, 19), columns=STAGE_ORDER).fillna(0)

    fig, ax = plt.subplots(figsize=(3.0, 2.0), facecolor=BG)
    ax.set_facecolor(BG)
    x = list(range(1, 19))
    width = 0.1
    offsets = [-0.35, -0.25, -0.15, -0.05, 0.05, 0.15, 0.25, 0.35]
    for stage, offset in zip(STAGE_ORDER, offsets):
        ax.bar([i + offset for i in x], pivot[stage].values, width=width, color=STAGE_COLORS[stage], edgecolor="none", linewidth=0.0)

    ax.set_xlim(0.45, 18.55)
    ax.set_xlabel("ChromHMM state", fontsize=7.2, color=TEXT)
    ax.set_ylabel("Genome coverage (%)", fontsize=7.2, color=TEXT)
    ax.set_xticks(range(1, 19))
    ax.tick_params(axis="x", labelsize=5.9, colors=TEXT, width=0.8, length=4)
    ax.tick_params(axis="y", labelsize=6.3, colors=TEXT, width=0.8, length=4)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.spines["left"].set_color(TEXT)
    ax.spines["bottom"].set_color(TEXT)
    ax.spines["left"].set_linewidth(0.8)
    ax.spines["bottom"].set_linewidth(0.8)
    fig.subplots_adjust(left=0.16, right=0.985, top=0.96, bottom=0.18)
    fig.savefig(output_path, dpi=600, facecolor=BG)
    plt.close(fig)


def plot_cd_legend(output_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(3, 0.2), facecolor=BG)
    ax.set_facecolor(BG)
    ax.axis("off")
    x_positions = [0.002 + i * 0.124 for i in range(len(STAGE_ORDER))]
    for x, label in zip(x_positions, STAGE_ORDER):
        ax.add_patch(plt.Rectangle((x, 0.18), 0.018, 0.64, transform=ax.transAxes, facecolor=STAGE_COLORS[label], edgecolor="none", linewidth=0))
        ax.text(x + 0.022, 0.5, label, transform=ax.transAxes, va="center", ha="left", fontsize=5.3, color=TEXT)
    fig.subplots_adjust(left=0.002, right=0.998, top=0.98, bottom=0.02)
    fig.savefig(output_path, dpi=600, facecolor=BG)
    plt.close(fig)


def is_canonical_chr(chrom: str) -> bool:
    return bool(re.match(r"^chr([1-9]|1[0-9]|X|Y)$", str(chrom)))


def load_segment_summary(input_dir: Path) -> pd.DataFrame:
    frames = []
    for stage in STAGE_ORDER:
        path = input_dir / f"{stage}_18_segments.bed"
        df = pd.read_csv(path, sep="\t", header=None, names=["chr", "start", "end", "state"])
        df = df[df["chr"].map(is_canonical_chr)].copy()
        df["state_num"] = df["state"].str.replace("E", "", regex=False).astype(int)
        chr_state = df.groupby(["chr", "state_num"]).size().reset_index(name="count")
        totals = chr_state.groupby("chr")["count"].transform("sum")
        chr_state["percent"] = chr_state["count"] / totals * 100
        chr_state["stage"] = stage
        frames.append(chr_state)
    return pd.concat(frames, ignore_index=True)


def chromosome_sort_key(chrom: str) -> tuple[int, int]:
    suffix = chrom.replace("chr", "")
    if suffix.isdigit():
        return (0, int(suffix))
    return (1, 1000 if suffix == "X" else 1001)


def plot_1e_single(summary: pd.DataFrame, stage: str, title: str, output_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(3, 2), facecolor=BG)
    ax.set_facecolor(BG)
    sub = summary[summary["stage"] == stage]
    pivot = sub.pivot(index="chr", columns="state_num", values="percent").reindex(columns=range(1, 19)).fillna(0)
    chrom_order = sorted(pivot.index, key=chromosome_sort_key)
    pivot = pivot.reindex(chrom_order)
    labels = [chrom.replace("chr", "") for chrom in chrom_order]

    bottom = pd.Series(0.0, index=pivot.index)
    handles = []
    for state_num in range(1, 19):
        bars = ax.bar(
            range(len(pivot.index)),
            pivot[state_num],
            bottom=bottom,
            color=STATE_COLORS[state_num - 1],
            width=0.9,
            edgecolor="none",
            linewidth=0,
            label=str(state_num),
        )
        handles.append(bars[0])
        bottom += pivot[state_num]

    ax.set_title(title, fontsize=6.5, color=TEXT, pad=2)
    ax.set_ylim(0, 100)
    ax.set_ylabel("Percentage of segments", fontsize=6.2, color=TEXT)
    ax.set_xlabel("Chromosome", fontsize=6.2, color=TEXT, labelpad=2)
    ax.set_yticks([0, 50, 100])
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=90, fontsize=4.3, color=TEXT)
    ax.tick_params(axis="y", labelsize=4.8, colors=TEXT, width=0.6, length=2)
    ax.tick_params(axis="x", colors=TEXT, width=0.6, length=2, pad=1)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.spines["left"].set_color(TEXT)
    ax.spines["bottom"].set_color(TEXT)
    ax.spines["left"].set_linewidth(0.6)
    ax.spines["bottom"].set_linewidth(0.6)

    legend = ax.legend(
        handles=handles,
        labels=[str(i) for i in range(1, 19)],
        ncol=6,
        frameon=False,
        fontsize=4.0,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.36),
        columnspacing=0.7,
        handlelength=0.8,
        handletextpad=0.3,
        borderaxespad=0,
    )
    for text in legend.get_texts():
        text.set_color(TEXT)

    fig.subplots_adjust(left=0.18, right=0.995, top=0.92, bottom=0.40)
    fig.savefig(output_path, dpi=600, facecolor=BG)
    plt.close(fig)


def plot_1e_legend(output_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(6, 0.2), facecolor=BG)
    ax.set_facecolor(BG)
    ax.axis("off")
    x_positions = [0.002 + i * (0.996 / 18.0) for i in range(18)]
    for i, x in enumerate(x_positions, start=1):
        ax.add_patch(plt.Rectangle((x, 0.18), 0.018, 0.64, transform=ax.transAxes, facecolor=STATE_COLORS[i - 1], edgecolor="none", linewidth=0))
        ax.text(x + 0.022, 0.5, str(i), transform=ax.transAxes, va="center", ha="left", fontsize=6.2, color=TEXT)
    fig.subplots_adjust(left=0.002, right=0.998, top=0.98, bottom=0.02)
    fig.savefig(output_path, dpi=600, facecolor=BG)
    plt.close(fig)


def load_integrated_emissions(input_dir: Path) -> pd.DataFrame:
    df = pd.read_csv(input_dir / "Int_ChromHMM_emissions.txt", sep="\t")
    df.columns = [c.strip() for c in df.columns]
    df = df.set_index(df.columns[0])
    desired = ["H3K27ac", "H3K9ac", "H3K4me1", "H3K4me2", "H3K4me3", "H3K36me3", "H3K27me3", "H3K9me3", "ATAC"]
    existing = [c for c in desired if c in df.columns]
    return df[existing]


def plot_1f_emissions(input_dir: Path, output_path: Path) -> None:
    df = load_integrated_emissions(input_dir)
    fig, ax = plt.subplots(figsize=(3, 2), facecolor=BG)
    ax.set_facecolor(BG)
    hm = sns.heatmap(df, cmap=sns.color_palette("light:b", as_cmap=True), cbar_kws={"label": "Emission Probability", "shrink": 0.8}, ax=ax)
    cbar = hm.collections[0].colorbar
    cbar.ax.tick_params(labelsize=5, colors=TEXT, width=0.5, length=2)
    cbar.set_label("Emission Probability", fontsize=5.5, color=TEXT)
    ax.set_ylabel("ChromHMM States", fontsize=6, color=TEXT, labelpad=2)
    ax.set_xlabel("Histone Marks", fontsize=6, color=TEXT, labelpad=2)
    ax.tick_params(axis="x", rotation=90, labelsize=5, colors=TEXT, width=0.5, length=2)
    ax.set_yticks([i + 0.5 for i in range(len(df.index))])
    ax.set_yticklabels([str(int(v)) for v in df.index], rotation=0)
    ax.tick_params(axis="y", rotation=0, labelsize=4.6, colors=TEXT, width=0.5, length=0)
    plt.tight_layout(pad=0.4)
    fig.savefig(output_path, dpi=600, facecolor=BG)
    plt.close(fig)


def load_stage_emissions(input_dir: Path) -> dict[str, pd.DataFrame]:
    files = {"INT": input_dir / "Int_ChromHMM_emissions.txt"}
    for stage in STAGE_ORDER:
        files[stage] = input_dir / f"{stage}_emission.txt"

    dfs: dict[str, pd.DataFrame] = {}
    for stage, path in files.items():
        df = pd.read_csv(path, sep="\t")
        df.columns = df.columns.str.strip()
        if "H3K36me3.1" in df.columns:
            df = df.drop(columns=["H3K36me3.1"])
        dfs[stage] = df
    return dfs


def build_jaccard_matrix(input_dir: Path) -> pd.DataFrame:
    dfs = load_stage_emissions(input_dir)
    state_col = dfs["INT"].columns[0]
    all_marks = [c for c in dfs["INT"].columns if c != state_col]

    for stage, df in dfs.items():
        for mark in all_marks:
            if mark not in df.columns:
                df[mark] = 0.0
        dfs[stage] = df[[state_col] + all_marks]

    def get_enriched_sets(df: pd.DataFrame, threshold: float = 0.8) -> dict[str, set[str]]:
        enriched = {}
        for _, row in df.iterrows():
            key = str(int(row[state_col]))
            enriched[key] = {mark for mark in all_marks if float(row[mark]) >= threshold}
        return enriched

    integrated_sets = get_enriched_sets(dfs["INT"])
    stage_sets = {stage: get_enriched_sets(df) for stage, df in dfs.items() if stage != "INT"}
    jaccard_matrix = pd.DataFrame(index=[str(i) for i in range(1, 19)], columns=STAGE_ORDER, dtype=float)

    def jaccard(a: set[str], b: set[str]) -> float:
        if not a and not b:
            return 0.0
        return len(a & b) / len(a | b)

    for i in range(1, 19):
        key = str(i)
        for stage in STAGE_ORDER:
            jaccard_matrix.loc[key, stage] = round(jaccard(integrated_sets[key], stage_sets[stage][key]), 2)
    return jaccard_matrix


def plot_1f_jaccard(input_dir: Path, output_path: Path) -> None:
    jaccard_matrix = build_jaccard_matrix(input_dir)
    fig, ax = plt.subplots(figsize=(3, 2), facecolor=BG)
    ax.set_facecolor(BG)
    hm = sns.heatmap(
        jaccard_matrix.astype(float),
        annot=False,
        cmap=plt.cm.OrRd,
        vmin=0,
        vmax=1,
        linewidths=0.3,
        linecolor="white",
        cbar_kws={"label": "Jaccard Similarity", "shrink": 0.8},
        ax=ax,
    )
    cbar = hm.collections[0].colorbar
    cbar.ax.tick_params(labelsize=5, colors=TEXT, width=0.5, length=2)
    cbar.set_label("Jaccard Similarity", fontsize=5.5, color=TEXT)
    ax.set_xlabel("Developmental Stages", fontsize=6, color=TEXT, labelpad=2)
    ax.set_ylabel("ChromHMM States", fontsize=6, color=TEXT, labelpad=2)
    ax.tick_params(axis="x", rotation=45, labelsize=5, colors=TEXT, width=0.5, length=2)
    ax.set_yticks([i + 0.5 for i in range(len(jaccard_matrix.index))])
    ax.set_yticklabels(list(jaccard_matrix.index), rotation=0)
    ax.tick_params(axis="y", rotation=0, labelsize=4.6, colors=TEXT, width=0.5, length=0)
    plt.tight_layout(pad=0.4)
    fig.savefig(output_path, dpi=600, facecolor=BG)
    plt.close(fig)


def generate_all(
    figure1_a_csv: Path,
    figure1_b_csv: Path,
    seg_input_dir: Path,
    output_dir: Path,
    dpi: int,
    panel_width: float,
    panel_height: float,
) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    out = build_paths(output_dir)

    draw_plot_a(load_table_a(figure1_a_csv), out["a"], panel_width, panel_height, dpi)
    draw_legend_a(out["a_legend"], dpi)
    draw_figure_b(load_states_b(figure1_b_csv), out["b"], panel_width, panel_height, dpi)

    seg_df = load_segment_tables(seg_input_dir)
    plot_fig_1c(seg_df, out["c"])
    plot_fig_1d(seg_df, out["d"])
    plot_cd_legend(out["cd_legend"])

    seg_summary = load_segment_summary(seg_input_dir)
    plot_1e_single(seg_summary, "E10", "E10.5", out["e_e10"])
    plot_1e_single(seg_summary, "P0", "P0", out["e_p0"])
    plot_1e_legend(out["e_legend"])

    plot_1f_emissions(seg_input_dir, out["f_emissions"])
    plot_1f_jaccard(seg_input_dir, out["f_jaccard"])

    return out


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--figure1-a-csv", default=str(DEFAULT_INPUT / "Figure_1_A.csv"))
    parser.add_argument("--figure1-b-csv", default=str(DEFAULT_INPUT / "Figure_1_B.csv"))
    parser.add_argument("--seg-input-dir", default=str(DEFAULT_INPUT / "Fig_1"))
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT))
    parser.add_argument("--dpi", type=int, default=600)
    parser.add_argument("--panel-width", type=float, default=3.0)
    parser.add_argument("--panel-height", type=float, default=2.0)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outputs = generate_all(
        figure1_a_csv=Path(args.figure1_a_csv),
        figure1_b_csv=Path(args.figure1_b_csv),
        seg_input_dir=Path(args.seg_input_dir),
        output_dir=Path(args.output_dir),
        dpi=args.dpi,
        panel_width=args.panel_width,
        panel_height=args.panel_height,
    )
    for path in outputs.values():
        print(path)


if __name__ == "__main__":
    main()
