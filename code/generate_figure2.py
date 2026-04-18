#!/usr/bin/env python3
"""Standalone Figure 2 generator for the GitHub package."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.colors import BoundaryNorm, ListedColormap


ROOT = Path(__file__).resolve().parent.parent
FIG1_INPUT = ROOT / "input" / "Fig_1"
FIG2_INPUT = ROOT / "input" / "Fig_2"
DEFAULT_OUTPUT = ROOT / "outputs" / "figure2"
BG = "#FFFDF8"
TEXT = "#000000"

EMISSION_WINDOWS = {
    "Figure_2_B_E10_5_to_E11_5.png": ("E10.5 to E11.5", ["E10", "E11"]),
    "Figure_2_B_E12_5_to_E14_5.png": ("E12.5 to E14.5", ["E12", "E13", "E14"]),
    "Figure_2_B_E15_5_to_E16_5.png": ("E15.5 to E16.5", ["E15", "E16"]),
    "Figure_2_B_P0.png": ("P0", ["P0"]),
}

PANEL_GROUPS = {
    "C": [
        ("ActPro_FE.txt", "Active promoter"),
        ("FlnkPro_FE.txt", "Flanking promoter"),
        ("BiPro.txt", "Bivalent promoter"),
        ("WkPro.txt", "Weak promoter"),
    ],
    "D": [
        ("ActEnh.txt", "Active enhancer"),
        ("BiEnh.txt", "Bivalent enhancer"),
        ("WkEnh.txt", "Weak enhancer"),
    ],
    "E": [
        ("Act_Trns.txt", "Active transcription"),
        ("TrGB.txt", "Transcribed gene body"),
        ("WK_TrGB.txt", "Weak transcribed"),
    ],
    "F": [
        ("RpPC.txt", "Repressive PolyComb"),
        ("Wk_RpPC.txt", "Weak Repressive PolyComb"),
        ("Het.txt", "Heterochromatin"),
    ],
}

CATEGORY_LABELS = {
    "CpGIsland": "CpG island",
    "Promoter": "Promoter",
    "Enhancer": "Enhancer",
    "RefSeqExon": "RefSeq exon",
    "RefSeqGene": "RefSeq gene",
    "RefSeqTES": "RefSeq TES",
    "RefSeqTSS": "RefSeq TSS",
    "RefSeqTSS2kb": "RefSeq TSS2kb",
}
CATEGORY_LABELS_SHORT = {
    "CpG island": "CpG",
    "Promoter": "Prom",
    "Enhancer": "Enh",
    "RefSeq exon": "Exon",
    "RefSeq gene": "Gene",
    "RefSeq TES": "TES",
    "RefSeq TSS": "TSS",
    "RefSeq TSS2kb": "TSS2kb",
}

WINDOW_LABELS = {
    "E10.5 to E11.5": "E10.5-E11.5",
    "E12.5 to E14.5": "E12.5-E14.5",
    "E15.5 to E16.5": "E15.5-E16.5",
    "P0": "P0",
}


def emission_cmap() -> tuple[ListedColormap, BoundaryNorm]:
    colors = ["#FFFFFF", "#FFFFFF", "#90EE90", "#1B7837"]
    bounds = [0.0, 0.2, 0.5, 0.8, 1.0]
    return ListedColormap(colors), BoundaryNorm(bounds, ncolors=len(colors))


def fe_cmap() -> tuple[ListedColormap, BoundaryNorm]:
    colors = ["#FFFFFF", "#FFFFFF", "#ADD8E6", "#00008B"]
    bounds = [0.0, 0.2, 0.5, 0.8, 1.0]
    return ListedColormap(colors), BoundaryNorm(bounds, ncolors=len(colors))


def load_emission_file(stage: str) -> pd.DataFrame:
    path = FIG1_INPUT / f"{stage}_emission.txt"
    df = pd.read_csv(path, sep="\t", index_col=0)
    df = df.loc[:, ~df.columns.astype(str).str.contains(r"^Unnamed")]
    if "H3K36me3.1" in df.columns:
        df = df.drop(columns=["H3K36me3.1"])
    return df


def average_emissions(stages: list[str]) -> pd.DataFrame:
    frames = [load_emission_file(stage) for stage in stages]
    combined = pd.concat(frames).groupby(level=0).mean()
    desired_order = [
        "H3K27ac",
        "H3K9ac",
        "H3K4me1",
        "H3K4me2",
        "H3K4me3",
        "H3K36me3",
        "H3K27me3",
        "H3K9me3",
    ]
    cols = [c for c in desired_order if c in combined.columns]
    return combined[cols]


def clean_fe_table(filename: str) -> pd.DataFrame:
    raw = pd.read_csv(FIG2_INPUT / filename, sep="\t")
    raw.columns = [str(c).strip() for c in raw.columns]
    stage_col = raw.columns[1]
    value_cols = raw.columns[2:]
    df = raw[[stage_col] + list(value_cols)].copy()
    df = df.rename(columns={stage_col: "window"})
    df["window"] = df["window"].astype(str).str.strip()
    df = df[df["window"] != "nan"].reset_index(drop=True)
    df = df.set_index("window")
    df.index = [WINDOW_LABELS.get(v, v) for v in df.index]
    df.columns = [CATEGORY_LABELS.get(c, c) for c in df.columns]
    return df.astype(float)


def save_emission_heatmap(title: str, stages: list[str], output_path: Path, width: float, height: float, dpi: int) -> None:
    cmap, norm = emission_cmap()
    df = average_emissions(stages)

    fig, ax = plt.subplots(figsize=(width, height), dpi=dpi, facecolor=BG)
    ax.set_facecolor(BG)
    hm = sns.heatmap(
        df,
        cmap=cmap,
        norm=norm,
        ax=ax,
        linewidths=0.25,
        linecolor="#D0D0D0",
        cbar=True,
        cbar_kws={"label": "Emission probability", "shrink": 0.86},
    )
    cbar = hm.collections[0].colorbar
    cbar.ax.tick_params(labelsize=5.2, colors=TEXT, width=0.5, length=2)
    cbar.set_label("Emission probability", fontsize=5.6, color=TEXT)

    ax.set_title(title, fontsize=7.0, color=TEXT, pad=4)
    ax.set_xlabel("Histone marks", fontsize=6.0, color=TEXT, labelpad=2)
    ax.set_ylabel("ChromHMM states", fontsize=6.0, color=TEXT, labelpad=2)
    ax.set_yticks([i + 0.5 for i in range(len(df.index))])
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right", fontsize=5.0, color=TEXT)
    ax.set_yticklabels([str(int(v)) for v in df.index], rotation=0, fontsize=5.0, color=TEXT)
    ax.tick_params(axis="x", colors=TEXT, width=0.5, length=2)
    ax.tick_params(axis="y", colors=TEXT, width=0.5, length=0)

    plt.tight_layout(pad=0.5)
    fig.savefig(output_path, dpi=dpi, facecolor=BG)
    plt.close(fig)


def save_fe_heatmap(
    title: str,
    filename: str,
    output_path: Path,
    width: float,
    height: float,
    dpi: int,
    show_cbar: bool = True,
) -> None:
    cmap, norm = fe_cmap()
    df = clean_fe_table(filename)
    df_norm = (df - df.min()) / (df.max() - df.min())
    df_norm = df_norm.fillna(0.0)

    fig, ax = plt.subplots(figsize=(width, height), dpi=dpi, facecolor=BG)
    ax.set_facecolor(BG)
    hm = sns.heatmap(
        df_norm,
        cmap=cmap,
        norm=norm,
        ax=ax,
        linewidths=0.25,
        linecolor="#FFFFFF",
        cbar=show_cbar,
        cbar_kws={"label": "Normalised enrichment", "shrink": 0.86},
    )
    if show_cbar:
        cbar = hm.collections[0].colorbar
        cbar.ax.tick_params(labelsize=5.2, colors=TEXT, width=0.5, length=2)
        cbar.set_label("Normalised enrichment", fontsize=5.6, color=TEXT)

    # Keep category separators from the original code.
    for pos in [1, 3]:
        ax.axvline(pos, color="white", lw=1.8)

    ax.set_title("" if width <= 1.8 else title, fontsize=3.8 if width <= 1.8 else 7.0, color=TEXT, pad=1)
    ax.set_xlabel("" if width <= 1.8 else "Genomic category", fontsize=4.8 if width <= 1.8 else 6.0, color=TEXT, labelpad=1)
    ax.set_ylabel("" if width <= 1.8 else "Developmental window", fontsize=4.8 if width <= 1.8 else 6.0, color=TEXT, labelpad=1)
    ax.set_yticks([i + 0.5 for i in range(len(df_norm.index))])
    if width <= 1.8:
        ax.set_xticklabels(
            [CATEGORY_LABELS_SHORT.get(c, c) for c in df_norm.columns],
            rotation=45,
            ha="right",
            fontsize=2.2,
            color=TEXT,
        )
        ax.set_yticklabels(
            [str(v).replace(".5-", "-").replace(".5", "") for v in df_norm.index],
            rotation=0,
            fontsize=2.9,
            color=TEXT,
        )
    else:
        ax.set_xticklabels(list(df_norm.columns), rotation=45, ha="right", fontsize=5.0, color=TEXT)
        ax.set_yticklabels(list(df_norm.index), rotation=0, fontsize=5.0, color=TEXT)
    ax.tick_params(axis="x", colors=TEXT, width=0.5, length=0)
    ax.tick_params(axis="y", colors=TEXT, width=0.5, length=0)

    if width <= 1.8:
        fig.subplots_adjust(left=0.18, right=0.99, top=0.98, bottom=0.38)
    else:
        plt.tight_layout(pad=0.5)
    fig.savefig(output_path, dpi=dpi, facecolor=BG)
    plt.close(fig)


def save_fe_legend(output_path: Path, width: float, height: float, dpi: int) -> None:
    cmap, norm = fe_cmap()
    fig, ax = plt.subplots(figsize=(width, height), dpi=dpi, facecolor=BG)
    fig.subplots_adjust(left=0.46, right=0.90, top=0.985, bottom=0.035)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=ax)
    cbar.set_label("", fontsize=5.5, color=TEXT, labelpad=4)
    cbar.ax.yaxis.set_ticks_position("left")
    cbar.ax.tick_params(labelsize=2.8, colors=TEXT, width=0.45, length=1.5, pad=0.5)
    cbar.set_ticks([0.0, 0.2, 0.5, 0.8, 1.0])
    cbar.set_ticklabels(["0", "0.2", "0.5", "0.8", "1"])
    fig.savefig(output_path, dpi=dpi, facecolor=BG)
    plt.close(fig)


def save_enrichment_heatmap(output_path: Path, width: float, height: float, dpi: int) -> None:
    data = pd.read_csv(FIG2_INPUT / "Pro_enrichHM.txt", sep="\t", index_col=0)
    data = data.replace("", pd.NA)
    data = data.apply(pd.to_numeric, errors="coerce")

    fig, ax = plt.subplots(figsize=(width, height), dpi=dpi, facecolor=BG)
    ax.set_facecolor(BG)
    hm = sns.heatmap(
        data,
        annot=True,
        fmt=".1f",
        cmap="YlGn",
        cbar_kws={"label": "Enrichment"},
        linewidths=0.2,
        linecolor="gray",
        mask=data.isna(),
        annot_kws={"size": 5.5, "color": "black"},
        ax=ax,
    )
    cbar = hm.collections[0].colorbar
    cbar.ax.tick_params(labelsize=5.8, colors=TEXT, width=0.5, length=2)
    cbar.set_label("Enrichment", fontsize=6.5, color=TEXT)

    ax.set_xlabel("Developmental stage", fontsize=7.5, color=TEXT)
    ax.set_ylabel("Chromatin state", fontsize=7.5, color=TEXT)
    ax.tick_params(axis="x", labelsize=6.2, rotation=45, colors=TEXT)
    ax.tick_params(axis="y", labelsize=5.8, rotation=0, colors=TEXT)

    fig.subplots_adjust(left=0.24, right=0.93, bottom=0.23, top=0.97)
    fig.savefig(output_path, dpi=dpi, facecolor=BG)
    plt.close(fig)


def generate_all(
    output_dir: Path,
    width: float,
    height: float,
    dpi: int,
    only_panel: str | None = None,
    b_height: float | None = None,
    include_enrichment: bool = True,
) -> list[Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    outputs: list[Path] = []

    if only_panel is None or only_panel.upper() == "B":
        current_b_height = b_height if b_height is not None else height
        for filename, (title, stages) in EMISSION_WINDOWS.items():
            out = output_dir / filename
            save_emission_heatmap(title, stages, out, width, current_b_height, dpi)
            outputs.append(out)

    for panel, specs in PANEL_GROUPS.items():
        if only_panel is not None and panel != only_panel.upper():
            continue
        panel_dir = output_dir / f"Figure_2_{panel}"
        panel_dir.mkdir(parents=True, exist_ok=True)
        for filename, title in specs:
            safe_name = title.lower().replace(" ", "_")
            out = panel_dir / f"Figure_2_{panel}_{safe_name}.png"
            save_fe_heatmap(title, filename, out, width, height, dpi)
            outputs.append(out)

    if include_enrichment and only_panel is None:
        out = output_dir / "Figure_2_enrichment_heatmap.png"
        save_enrichment_heatmap(out, 3.2, 2.75, dpi)
        outputs.append(out)

    return outputs


def generate_compact_panel_set(
    output_dir: Path,
    panel_key: str,
    panel_width: float,
    panel_height: float,
    legend_width: float,
    legend_height: float,
    dpi: int,
) -> list[Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    panel_key = panel_key.upper()
    panel_dir = output_dir / f"Figure_2_{panel_key}"
    panel_dir.mkdir(parents=True, exist_ok=True)
    outputs: list[Path] = []

    for filename, title in PANEL_GROUPS[panel_key]:
        safe_name = title.lower().replace(" ", "_")
        out = panel_dir / f"Figure_2_{panel_key}_{safe_name}.png"
        save_fe_heatmap(title, filename, out, panel_width, panel_height, dpi, show_cbar=False)
        outputs.append(out)

    legend_out = output_dir / f"Figure_2_{panel_key}_legend.png"
    save_fe_legend(legend_out, legend_width, legend_height, dpi)
    outputs.append(legend_out)
    return outputs


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT), help="Directory for Figure 2 outputs.")
    parser.add_argument("--width", type=float, default=3.2, help="Panel width in inches.")
    parser.add_argument("--height", type=float, default=2.2, help="Panel height in inches.")
    parser.add_argument("--dpi", type=int, default=600, help="Output DPI.")
    parser.add_argument("--only-panel", choices=["B", "C", "D", "E", "F", "b", "c", "d", "e", "f"], help="Generate only one Figure 2 panel set.")
    parser.add_argument("--b-height", type=float, help="Optional custom height in inches for Figure 2B heatmaps.")
    parser.add_argument("--compact-panel", choices=["C", "D", "E", "F", "c", "d", "e", "f"], help="Generate one compact panel set plus a shared legend.")
    parser.add_argument("--panel-width", type=float, default=1.75, help="Compact Figure 2C panel width in inches.")
    parser.add_argument("--panel-height", type=float, default=0.7, help="Compact Figure 2C panel height in inches.")
    parser.add_argument("--legend-width", type=float, default=0.2, help="Standalone Figure 2C legend width in inches.")
    parser.add_argument("--legend-height", type=float, default=2.0, help="Standalone Figure 2C legend height in inches.")
    parser.add_argument("--skip-enrichment", action="store_true", help="Skip the promoter enrichment heatmap export.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.compact_panel:
        outputs = generate_compact_panel_set(
            Path(args.output_dir),
            args.compact_panel,
            args.panel_width,
            args.panel_height,
            args.legend_width,
            args.legend_height,
            args.dpi,
        )
    else:
        outputs = generate_all(
            Path(args.output_dir),
            args.width,
            args.height,
            args.dpi,
            only_panel=args.only_panel,
            b_height=args.b_height,
            include_enrichment=not args.skip_enrichment,
        )
    for path in outputs:
        print(path)


if __name__ == "__main__":
    main()
