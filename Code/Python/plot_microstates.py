"""
This script loads and combines:
      (1) microstate data
      (2) behavioral data
      (3) demographic data
      (4) microstate data from Thiele et al. (2023)

and plots the mean and standard deviation of microstate measures.

The original script was written by Jonas A. Thiele (2022) for the paper:
"Multimodal Brain Signal Complexity Predicts Human Intelligence" [1].
Orignial script is available at:
https://github.com/jonasAthiele/BrainComplexity_Intelligence/tree/main

The script was modified by Christoph Frühlinger (2025) using Python 3.13.3
in Visual Studio Code.

References:
[1] Thiele, J. A., Richter, A. & Hilger, K. (2023). Multimodal brain signal
complexity predicts human intelligence. eNeuro, 10(2), ENEURO.0345-22.2022.
https://doi.org/10.1523/eneuro.0345-22.2022
"""

# %% Imports

import os
import re
from pathlib import Path

import mne
import numpy as np
import pandas as pd
import scipy.stats
import seaborn as sns
from matplotlib import pyplot as plt

# %% Load data

cwd = Path.cwd()
dir_Root = cwd.parent.parent
dir_Data = dir_Root + "/Data"
dir_microstates = dir_Data + "/Microstates"
dir_preprocessed = dir_Data + "/Snipplet"  # Epoched data
dir_results = dir_Root + "/Results"

Participants = os.listdir(dir_Data + "/FluidData/task-IST")

gf_data = pd.DataFrame()
gf_data["ID"] = Participants

fluid_correct = pd.read_excel(
    dir_Root + "/Data/FluidData/IST_fluid_A.xlsx", header=None
)
fluid_correct = fluid_correct[1]

gf_scores = []

for idx, i_sub in enumerate(Participants):

    filepath = Path(dir_Root) / "Data" / "FluidData" / "task-IST" / i_sub / "beh"
    file = list(filepath.glob("*Fluid_beh.csv"))

    if not file:
        gf_scores.append(np.nan)
        continue

    gf = pd.read_csv(file[0], usecols=["ratings"])
    gf = gf.iloc[1:-1, 0].astype(str).reset_index(drop=True)

    if pd.isna(gf).all():
        gf_scores.append(np.nan)
        continue

    min_len = min(len(fluid_correct), len(gf))
    score = np.sum(gf[:min_len].to_numpy() == fluid_correct[:min_len].to_numpy())
    gf_scores.append(score)

gf_data["gf_score"] = gf_scores
gf_data["gf_score"] = scipy.stats.zscore(gf_data["gf_score"], nan_policy="omit", ddof=1)

demographics = pd.read_csv(
    dir_Root + "/Data/SocioDemographics.txt", encoding="unicode_escape"
)
demographics = demographics.loc[:, ["ID", "Gender", "Age"]]
demographics["Gender"] = (
    demographics["Gender"].replace({1: "female", 2: "male"}).astype("category")
)

beh_data = pd.merge(gf_data, demographics, on="ID", how="outer")
# beh_data = beh_data[beh_data['Age'].notna()]

intell_main = beh_data["gf_score"].to_numpy()
age_main = beh_data["Age"].to_numpy()  # age
subs_beh = beh_data["ID"].to_numpy()  # subject IDs corresponding to behavioral data

# Load group microstate data and check order
ms_path = Path(dir_microstates)
maps_group_files = list(ms_path.glob("microstates_group*.npy"))

fname = "/40_seconds_first_run_eyes_closed_sub-AA06WI11_task-Resting_run-1_eeg.set"

preproc = mne.io.read_raw_eeglab(
    dir_preprocessed + fname, preload=True
)  # raw for topomap plotting

for maps_group_file in maps_group_files:
    maps_group_data = np.load(maps_group_file)
    fig, axes = plt.subplots(nrows=1, ncols=5)
    cnt = 0
    print(maps_group_file)
    for axes_row in axes:

        mne.viz.plot_topomap(
            maps_group_data[cnt], preproc.info, axes=axes_row, show=False
        )
        axes_row.spines["right"].set_visible(False)
        axes_row.spines["top"].set_visible(False)

        cnt += 1
    fig.tight_layout()

ms_order_run1EC = [0, 4, 1, 2, 3]
ms_order_run1EO = [4, 3, 1, 2, 0]
ms_order_run2EC = [2, 4, 3, 0, 1]
ms_order_run2EO = [4, 0, 1, 3, 2]
ms_order_run3EC = [1, 2, 3, 0, 4]
ms_order_run3EO = [1, 4, 0, 3, 2]

ms_order_list = [
    ms_order_run1EC,
    ms_order_run1EO,
    ms_order_run2EC,
    ms_order_run2EO,
    ms_order_run3EC,
    ms_order_run3EO,
]

for cnt, maps_group_file in enumerate(maps_group_files):
    maps_group_data = np.load(maps_group_file)  # shape: (5, n_channels)
    maps_group_data = maps_group_data[ms_order_list[cnt], :]  # reorder maps
    fig, axes = plt.subplots(nrows=1, ncols=5)
    print(maps_group_file)

    for i, ax in enumerate(axes):
        mne.viz.plot_topomap(maps_group_data[i], preproc.info, axes=ax, show=False)
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)

    fig.tight_layout()

# Rename microstate labels from digits to letters
ms_files = list(ms_path.glob("df_microstate*.pkl"))

ms_dfs = []
condition_labels = ["run1EC", "run1EO", "run2EC", "run2EO", "run3EC", "run3EO"]

alphabet_labels = ["A", "B", "C", "D", "F"]
for ms_file, order, condition in zip(ms_files, ms_order_list, condition_labels):
    df = pd.read_pickle(ms_file)

    label_map = {
        str(orig_idx): alphabet_labels[new_idx]
        for new_idx, orig_idx in enumerate(order)
    }

    def replace_microstate_digits(col_name):

        col_name = re.sub(
            r"_(\d)_(\d)",
            lambda m: f"_{label_map.get(m.group(1),
                                        m.group(1))}_{label_map.get(m.group(2),
                                                                    m.group(2))}",
            col_name,
        )

        col_name = re.sub(
            r"_(\d)\b", lambda m: f"_{label_map.get(m.group(1), m.group(1))}", col_name
        )
        return col_name

    df = df.rename(columns={col: replace_microstate_digits(col) for col in df.columns})

    if "Condition" not in df.columns:
        df["Condition"] = condition

    ms_dfs.append(df)

df_microstate_all = pd.concat(ms_dfs, ignore_index=True)

all_columns = set()
for df in ms_dfs:
    all_columns.update(df.columns)
df_microstate_all = df_microstate_all.reindex(columns=sorted(all_columns))

# Sort data
df_microstate_all = df_microstate_all.sort_values(by=["ID", "Condition"]).reset_index(
    drop=True
)

# Remove rows with no microstate data
df_microstate_all = df_microstate_all.dropna(subset=["Length"])

# Merge data with behavioral data
df_microstate_all = pd.merge(beh_data, df_microstate_all, on="ID", how="outer")

if not os.path.exists(dir_microstates + "/Microstates_data_full.csv"):
    df_microstate_all.to_csv(
        dir_microstates + "/Microstates_data_full.csv",
        index=False,
        header=True,
    )
# Load microstate data from Thiele et al. (2023)
df_microstate_thiele = pd.read_pickle(
    "C:/Users/Christoph Frühlinger/Nextcloud/PhD/Forschung/Complexity/Data/"
    "Microstates/Thiele Data/df_complexity_org.pkl"
)

# Rename microstates from digits to letters
ms_thiele_order = [4, 3, 1, 2, 0]
label_map = {
    str(orig_idx): alphabet_labels[new_idx] for new_idx, orig_idx in enumerate(order)
}
df_microstate_thiele = df_microstate_thiele.rename(
    columns={
        col: replace_microstate_digits(col) for col in df_microstate_thiele.columns
    }
)

# Prepare df for merging
df_microstate_thiele = df_microstate_thiele.rename(columns={"IDs_subjects_neuro": "ID"})
df_microstate_thiele["Condition"] = ["Thiele et al. (2023)"] * len(df_microstate_thiele)

common_cols = [
    col for col in df_microstate_all.columns if col in df_microstate_thiele.columns
]
combined_ms = pd.concat(
    [df_microstate_all[common_cols], df_microstate_thiele[common_cols]],
    ignore_index=True,
)

# Plot mean and standard deviation of microstate measures
feature_cols = [
    col
    for col in combined_ms.columns
    if col.startswith(("coverage_", "frequence_", "lifespan_"))
]

# Transform to long format for plotting
df_long = combined_ms.melt(
    id_vars=["ID", "Condition"],
    value_vars=feature_cols,
    var_name="Feature",
    value_name="Value",
)

df_long[["Feature", "Microstate"]] = df_long["Feature"].str.rsplit(
    "_", n=1, expand=True
)

n_features = df_long["Feature"].nunique()
col_wrap = 5

labels_feature = {
    "coverage": "Coverage",
    "frequence": "Frequency",
    "lifespan": "Lifespan",
    "lifespan_peaks": "Lifespan at GFP Peaks",
}

df_long["Feature"] = df_long["Feature"].map(labels_feature)

labels_condition = {
    "first_run_eyes_open": "Run 1 EO",
    "first_run_eyes_closed": "Run 1 EC",
    "second_run_eyes_open": "Run 2 EO",
    "second_run_eyes_closed": "Run 2 EC",
    "third_run_eyes_open": "Run 3 EO",
    "third_run_eyes_closed": "Run 3 EC",
    "Thiele et al. (2023)": "Thiele et al. (2023)",
}

df_long["Condition"] = df_long["Condition"].map(labels_condition)

plot_order = ["Run 1 EO", "Run 1 EC", "Run 2 EO", "Run 3 EO", "Thiele et al. (2023)"]

# Plot
sns.set_style("whitegrid")
g = sns.FacetGrid(
    df_long,
    row="Feature",
    col="Microstate",
    hue="Condition",
    palette="viridis",
    margin_titles=True,
    height=3.5,
    sharey=False,
)

g.map_dataframe(
    sns.barplot, x="Condition", y="Value", errorbar="sd", alpha=0.9, order=plot_order
)

# Achsen & Layout anpassen
g.set_axis_labels("", "Mean", fontsize=22)
for ax in g.axes.flat:
    ax.set_xticklabels([])
    ax.set_xlabel("")
    ax.tick_params(axis="y", labelsize=18)
g.set_titles(
    row_template="{row_name}", col_template="{col_name}", size=22, weight="bold"
)
g.add_legend(title="Condition", label_order=plot_order, fontsize=18)
g._legend.get_title().set_fontsize(22)
plt.subplots_adjust(top=0.92)
plt.show()

g.savefig(
    dir_results + "/microstate_mean_sd.tiff",
    format="tiff",
    dpi=600,
    bbox_inches="tight",
)


# Plot Transition Probabilities
feature_cols = [col for col in combined_ms.columns if col.startswith(("transition_"))]

# Transform to long format for plotting
df_long = combined_ms.melt(
    id_vars=["ID", "Condition"],
    value_vars=feature_cols,
    var_name="Feature",
    value_name="Value",
)

df_long[["Feature", "From", "To"]] = df_long["Feature"].str.rsplit(
    "_", n=2, expand=True
)

labels_feature = {
    "transition_probability": "Transition Probability",
    "transition_probability_peaks": "Transition Probability at GFP Peaks",
}

df_long["Feature"] = df_long["Feature"].map(labels_feature)

labels_condition = {
    "first_run_eyes_open": "Run 1 EO",
    "first_run_eyes_closed": "Run 1 EC",
    "second_run_eyes_open": "Run 2 EO",
    "second_run_eyes_closed": "Run 2 EC",
    "third_run_eyes_open": "Run 3 EO",
    "third_run_eyes_closed": "Run 3 EC",
    "Thiele et al. (2023)": "Thiele et al. (2023)",
}

df_long["Condition"] = df_long["Condition"].map(labels_condition)

plot_order = ["Run 1 EO", "Run 1 EC", "Run 2 EO", "Run 3 EO", "Thiele et al. (2023)"]

features = df_long["Feature"].unique()
n_rows = len(features)
n_cols = len(plot_order)

fig, axes = plt.subplots(n_rows, n_cols, figsize=(7 * n_cols, 7.2 * n_rows))

for i, feature in enumerate(features):
    for j, condition in enumerate(plot_order):
        ax = axes[i][j]

        subset = df_long[
            (df_long["Feature"] == feature) & (df_long["Condition"] == condition)
        ]

        pivot = subset.pivot_table(
            index="From", columns="To", values="Value", aggfunc="mean"
        )

        # Plot als Heatmap
        sns.heatmap(
            pivot,
            ax=ax,
            cmap="RdBu_r",
            cbar=False,
            vmin=0,
            vmax=1,
            square=True,
            annot=False,
        )

        if i == 0:
            ax.set_title(condition, fontsize=22, weight="bold")
        else:
            ax.set_title("")

        if j == 0:
            ax.set_ylabel(feature, fontsize=22, weight="bold")
        else:
            ax.set_ylabel("")

        ax.set_xlabel("")
        ax.tick_params(rotation=0, labelsize=18)

sm = plt.cm.ScalarMappable(cmap="RdBu_r")
sm.set_array([])
cbar_ax = fig.add_axes([0.92, 0.15, 0.015, 0.7])
cbar = fig.colorbar(sm, cax=cbar_ax)
cbar.ax.tick_params(labelsize=18)
plt.show()

fig.savefig(
    dir_results + "/microstate_trans_prob_mean.tiff",
    format="tiff",
    dpi=600,
    bbox_inches="tight",
)
