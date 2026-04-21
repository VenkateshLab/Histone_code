# -*- coding: utf-8 -*-
"""
Unified generator for the current Figure 3 panels:

1. Fig3_AD_Emission.png
2. Fig3_AD_FE.png
3. Fig3D_gcover.png

All required inputs are expected under Input/Fig_3.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm, to_rgba
from matplotlib.cm import ScalarMappable
from matplotlib.patches import Rectangle
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list
from scipy.spatial.distance import pdist
from sklearn.preprocessing import StandardScaler


BASE_DIR = Path(__file__).resolve().parent
INPUT_DIR = BASE_DIR / "Input" / "Fig_3"
OUTPUT_DIR = BASE_DIR / "outputs" / "figure3"

STATE_STYLE_MAP = {
    1: {"label": "Act Pro", "color": "#DD7E79"},
    2: {"label": "Wk Pro", "color": "#F3CDCC"},
    3: {"label": "Act Enh", "color": "#F4B183"},
    4: {"label": "Wk Enh", "color": "#F4E6AA"},
    5: {"label": "Acc", "color": "#A9D18E"},
    6: {"label": "Q/UnRp", "color": "#BFBFBF"},
}


def apply_common_matplotlib_style(font_size=7):
    plt.rcParams["font.family"] = "Arial"
    plt.rcParams["font.size"] = font_size
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42


def build_annotated_matrix_rgba(data_df, state_ids, cmap, norm):
    n_rows, n_data_cols = data_df.shape
    n_total_cols = n_data_cols + 1
    rgba = np.ones((n_rows, n_total_cols, 4), dtype=float)

    for row, state_id in enumerate(state_ids):
        rgba[row, 0, :] = to_rgba(STATE_STYLE_MAP.get(state_id, {"color": "#FFFFFF"})["color"])

    rgba[:, 1:, :] = cmap(norm(data_df.to_numpy()))
    return rgba


def draw_annotated_single_panel(
    ordered_df,
    state_ids,
    output_path,
    colors,
    xlabel,
    title,
    figsize,
    x_tick_rotation,
    x_tick_fontsize,
    left_text_fontsize,
    cbar_ticks,
    x_tick_ha="center",
    x_tick_rotation_mode=None,
    transparent=False,
):
    apply_common_matplotlib_style(font_size=7)

    bounds = [0.0, 0.2, 0.5, 0.8, 1.0]
    cmap = ListedColormap(colors)
    norm = BoundaryNorm(bounds, ncolors=len(colors))

    n_rows, n_data_cols = ordered_df.shape
    n_total_cols = n_data_cols + 1
    rgba = build_annotated_matrix_rgba(ordered_df, state_ids, cmap, norm)

    fig = plt.figure(figsize=figsize, dpi=600)
    gs = fig.add_gridspec(nrows=1, ncols=2, width_ratios=[5.2, 0.35], wspace=0.05)
    ax = fig.add_subplot(gs[0, 0])
    cbar_ax = fig.add_subplot(gs[0, 1])

    ax.imshow(rgba, aspect="auto", interpolation="none", origin="upper")

    for row, state_id in enumerate(state_ids):
        label = STATE_STYLE_MAP.get(state_id, {"label": f"State {state_id}"})["label"]
        ax.text(
            0,
            row,
            label,
            ha="center",
            va="center",
            fontsize=left_text_fontsize,
            color="black",
        )

    ax.set_xticks(np.arange(n_total_cols))
    ax.set_yticks(np.arange(n_rows))

    ax.set_xticklabels(
        [""] + list(ordered_df.columns),
        rotation=x_tick_rotation,
        ha=x_tick_ha,
        rotation_mode=x_tick_rotation_mode,
        fontsize=x_tick_fontsize,
    )
    ax.set_yticklabels([])
    ax.tick_params(axis="x", which="major", bottom=True, top=True, labelbottom=True, length=0, pad=2)
    ax.tick_params(axis="y", which="major", left=True, right=True, labelleft=False, length=0)

    ax.set_xlim(-0.5, n_total_cols - 0.5)
    ax.set_ylim(n_rows - 0.5, -0.5)

    ax.set_title(title, fontsize=8, weight="bold", pad=4)
    ax.set_xlabel(xlabel, fontsize=6.5, labelpad=2)
    ax.set_ylabel("State", fontsize=7, labelpad=10, rotation=90, va="center")
    ax.yaxis.set_label_coords(-0.09, 0.5)

    for row in range(n_rows):
        for col in range(1, n_total_cols):
            ax.add_patch(
                Rectangle(
                    (col - 0.5, row - 0.5),
                    1,
                    1,
                    fill=False,
                    edgecolor="gray",
                    linewidth=0.5,
                    linestyle=":",
                )
            )

    ax.spines["left"].set_position(("data", 0.5))
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.6)
        spine.set_color("black")

    ax.spines["bottom"].set_bounds(0.5, n_total_cols - 0.5)
    ax.spines["top"].set_bounds(0.5, n_total_cols - 0.5)

    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cbar_ax, ticks=cbar_ticks)
    cbar.ax.tick_params(labelsize=5, length=2, pad=2)

    if figsize[0] <= 2.8:
        fig.subplots_adjust(left=0.12, right=0.93, bottom=0.23, top=0.84)
    else:
        fig.subplots_adjust(left=0.14, right=0.94, bottom=0.33, top=0.84)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, format="png", dpi=600, transparent=transparent)
    plt.close(fig)


def generate_emission():
    input_path = INPUT_DIR / "Adult_EO.txt"
    output_path = OUTPUT_DIR / "Fig3_AD_Emission.png"

    df = pd.read_csv(input_path, sep="\t")
    df.set_index("State (Emission order)", inplace=True)

    desired_order = [
        "H3K27ac", "H3K9ac", "H3K4me1", "H3K4me2",
        "H3K4me3", "H3K36me3", "H3K27me3", "H3K9me3", "ATAC",
    ]
    ordered_cols = [col for col in desired_order if col in df.columns]
    ordered_df = df[ordered_cols]
    state_ids = [int(idx) for idx in ordered_df.index]

    draw_annotated_single_panel(
        ordered_df=ordered_df,
        state_ids=state_ids,
        output_path=output_path,
        colors=["#FFFFFF", "#FFFFFF", "#90EE90", "#1B7837"],
        xlabel="Histone Marks",
        title="ChromHMM - Adult",
        figsize=(2.8, 2.2),
        x_tick_rotation=0,
        x_tick_fontsize=5.5,
        left_text_fontsize=5.2,
        cbar_ticks=[0.0, 0.2, 0.5, 0.8, 1.0],
        transparent=False,
    )
    print(f"Saved figure to {output_path}")


def generate_fold_enrichment():
    input_path = INPUT_DIR / "Ad_overlap.txt"
    output_path = OUTPUT_DIR / "Fig3_AD_FE.png"

    df = pd.read_csv(input_path, sep="\t", index_col=0)
    df = df.loc[:, ~df.columns.str.contains("^Unnamed")]
    df = df.loc[[idx for idx in df.index if str(idx) in {"1", "2", "3", "4", "5", "6"}]]

    normalized_df = (df - df.min()) / (df.max() - df.min())
    normalized_df = normalized_df.fillna(0)

    desired_order = ["CpGIsland", "Promoter", "Enhancer"]
    ordered_cols = desired_order + [col for col in normalized_df.columns if col not in desired_order]
    ordered_df = normalized_df[ordered_cols]
    state_ids = [int(str(idx)) for idx in ordered_df.index]

    draw_annotated_single_panel(
        ordered_df=ordered_df,
        state_ids=state_ids,
        output_path=output_path,
        colors=["#FFFFFF", "#FFFFFF", "#ADD8E6", "#00008B"],
        xlabel="Genomic Features",
        title="ChromHMM - Adult",
        figsize=(3.0, 2.2),
        x_tick_rotation=45,
        x_tick_fontsize=5.3,
        left_text_fontsize=4.8,
        cbar_ticks=[0.0, 0.2, 0.5, 0.8, 1.0],
        x_tick_ha="right",
        x_tick_rotation_mode="anchor",
        transparent=False,
    )
    print(f"Saved figure to {output_path}")


def generate_genome_coverage():
    apply_common_matplotlib_style(font_size=6)
    matplotlib.rcParams["axes.linewidth"] = 1.0
    matplotlib.rcParams["xtick.major.width"] = 1.0
    matplotlib.rcParams["ytick.major.width"] = 1.0

    input_path = INPUT_DIR / "EtoA_genome.txt"
    output_path = OUTPUT_DIR / "Fig3D_gcover.png"

    df = pd.read_csv(input_path, sep="\t", index_col=0)
    df = df.apply(pd.to_numeric, errors="coerce")
    df = df.dropna(how="all", axis=0).dropna(how="all", axis=1)

    scaler = StandardScaler()
    df_scaled = pd.DataFrame(
        scaler.fit_transform(df),
        index=df.index,
        columns=df.columns,
    ).fillna(0)

    state_colors = {
        "Act Pro": (0.753, 0.220, 0.188),
        "Flnk Pro": (0.918, 0.925, 0.039),
        "Bi Pro": (0.878, 0.498, 0.502),
        "Wk Pro": (0.953, 0.804, 0.800),
        "Act Enh": (0.906, 0.529, 0.169),
        "Bi Enh": (0.949, 0.702, 0.259),
        "Wk Enh": (0.949, 0.702, 0.259),
        "Wk Bi Enh": (0.957, 0.902, 0.667),
        "Poised": (0.282, 0.624, 0.655),
        "Acc": (0.722, 0.886, 0.871),
        "Alt Splicing": (0.353, 0.667, 0.275),
        "Strong Transcribed": (0.580, 0.769, 0.490),
        "Transcribed": (0.784, 0.859, 0.784),
        "Rp PC": (0.729, 0.706, 0.847),
        "Wk Rp PC": (0.871, 0.859, 0.933),
        "Het": (0.510, 0.361, 0.682),
        "Quiescent": (0.949, 0.949, 0.973),
    }

    df_dist = pdist(df_scaled, metric="euclidean")
    df_linkage = linkage(df_dist, method="ward")
    order = leaves_list(df_linkage)
    df_display = df_scaled.iloc[order]
    n_rows, n_cols = df_display.shape

    color_array = np.zeros((n_rows, n_cols, 3))
    min_val = df_display.min().min()
    max_val = df_display.max().max()

    for i, state in enumerate(df_display.index):
        base_color = np.array(state_colors.get(state, (0.5, 0.5, 0.5)))
        z_scores = df_display.iloc[i].values
        for j, z_val in enumerate(z_scores):
            if z_val <= 0:
                intensity = (z_val - min_val) / (0 - min_val) if min_val < 0 else 0
                color_array[i, j] = base_color * intensity + (1 - intensity)
            else:
                intensity = z_val / max_val if max_val > 0 else 0
                color_array[i, j] = base_color * (0.3 + 0.7 * intensity)

    fig = plt.figure(figsize=(3.4, 2.0), dpi=600)
    gs = fig.add_gridspec(
        1,
        3,
        width_ratios=[0.26, 0.22, 0.85],
        wspace=0.03,
        left=0.12,
        right=0.96,
        top=0.84,
        bottom=0.14,
    )

    ax_dend = fig.add_subplot(gs[0])
    dendrogram(
        df_linkage,
        orientation="left",
        no_labels=True,
        color_threshold=0,
        above_threshold_color="black",
        ax=ax_dend,
    )
    ax_dend.invert_yaxis()
    ax_dend.set_xticks([])
    ax_dend.set_yticks([])
    for spine in ax_dend.spines.values():
        spine.set_visible(False)

    ax_labels = fig.add_subplot(gs[1])
    ax_labels.set_xlim(0, 1)
    ax_labels.set_ylim(0, n_rows)
    ax_labels.axis("off")
    for i, state in enumerate(df_display.index):
        ax_labels.text(0.0, i + 0.5, state, ha="left", va="center", fontsize=5.5, family="sans-serif")

    ax_heat = fig.add_subplot(gs[2])
    for i in range(n_rows):
        for j in range(n_cols):
            ax_heat.add_patch(
                Rectangle(
                    (j, i),
                    1,
                    1,
                    facecolor=color_array[i, j],
                    edgecolor="white",
                    linewidth=0.8,
                )
            )

    ax_heat.set_xticks(np.arange(n_cols) + 0.5)
    ax_heat.set_xticklabels(df_display.columns, rotation=90, ha="right", fontsize=5.5)
    ax_heat.set_yticks([])
    ax_heat.set_xlim(0, n_cols)
    ax_heat.set_ylim(0, n_rows)
    ax_heat.tick_params(left=False, bottom=False, pad=2)
    for spine in ["top", "right", "left", "bottom"]:
        ax_heat.spines[spine].set_visible(False)
    ax_heat.set_xlabel("Developmental Age", fontsize=6.5, fontweight="bold", labelpad=3)

    for i in range(n_rows):
        for j in range(n_cols):
            value = df_display.iloc[i, j]
            ax_heat.text(j + 0.5, i + 0.5, f"{value:.2f}", ha="center", va="center", fontsize=4.5, color="black")

    fig.suptitle("Hierarchical Clustering of Chromatin States", fontsize=7, fontweight="bold", y=0.95)
    fig.text(0.03, 0.5, "Chromatin States", va="center", rotation="vertical", fontsize=6.5, fontweight="bold")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=600, format="png", transparent=True)
    plt.close(fig)
    print(f"Saved figure to {output_path}")


def main():
    generate_emission()
    generate_fold_enrichment()
    generate_genome_coverage()


if __name__ == "__main__":
    main()
