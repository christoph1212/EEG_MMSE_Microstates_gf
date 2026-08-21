"""
This script plots the correlation results between entropy measures and gf.
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

# %%
import pickle
from pathlib import Path

import matplotlib.pyplot as plt
import mne
import numpy as np
import pandas as pd
from matplotlib import colors, ticker

# %% Load data
cwd = Path.cwd()
dir_root = cwd.parent.parent
dir_data = dir_root / "Data"
dir_preprocessed = dir_data / "Snipplet"
dir_entropy = dir_data / "Entropy_Thiele"
dir_microstates = dir_data / "Microstates"

df_main = pd.read_pickle(dir_entropy / "df_entropy_results_first_run_eyes_open.pkl")

column_names = np.array(df_main.columns)

with open(dir_microstates / "channel_names", "rb") as fp:
    channel_names = pickle.load(fp)

fname = "40_seconds_first_run_eyes_closed_sub-AA06WI11_task-Resting_run-1_eeg.set"

raw = mne.io.read_raw_eeglab(dir_preprocessed / fname, preload=True)

# %% Associations between single complexity measures and intelligence

# Shannon entropy and intelligence
idx_vars = []
for ch in channel_names:
    var_name = "shannon_" + ch
    idx_vars.append(np.where(column_names == var_name)[0][0])
idx_vars = np.array(idx_vars)

fig, a = plt.subplots(figsize=(12, 6))
plt.bar(
    range(idx_vars.size),
    np.array(df_main)[0, idx_vars],
    color="red",
    edgecolor="black",
    alpha=0.5,
)
plt.ylim([-0.15, 0.15])
plt.grid(color="gray", linestyle="dashed", linewidth=0.5)
plt.yticks(fontsize="9")
plt.xticks(np.arange(idx_vars.size), np.array(channel_names), fontsize="9", rotation=55)
# plt.xlabel('Channel')
a.spines["right"].set_visible(False)
a.spines["top"].set_visible(False)
plt.locator_params(axis="y", nbins=5)
plt.savefig(
    dir_entropy / "shannon_RAPM.jpg", format="jpg", dpi=1200, bbox_inches="tight"
)


# Fuzzy entropy and intelligence
idx_vars = []
for ch in channel_names:
    var_name = "fuzzy_" + ch
    idx_vars.append(np.where(column_names == var_name)[0][0])
idx_vars = np.array(idx_vars)

fig, a = plt.subplots(figsize=(12, 6))
plt.bar(
    range(idx_vars.size),
    np.array(df_main)[0, idx_vars],
    color="red",
    edgecolor="black",
    alpha=0.5,
)
plt.ylim([-0.15, 0.15])
plt.grid(color="gray", linestyle="dashed", linewidth=0.5)
plt.yticks(fontsize="9")
plt.xticks(np.arange(idx_vars.size), np.array(channel_names), fontsize="9", rotation=55)
a.spines["right"].set_visible(False)
a.spines["top"].set_visible(False)
plt.locator_params(axis="y", nbins=5)
plt.savefig(dir_entropy / "fuzzy_RAPM.jpg", format="jpg", dpi=1200, bbox_inches="tight")


# MSE of GFP and intelligence
idx_vars = []
for sc in range(20):
    var_name = "mse_gfp_scale_" + str(sc)
    idx_vars.append(np.where(column_names == var_name)[0][0])
idx_vars = np.array(idx_vars)

fig, a = plt.subplots()
plt.bar(
    range(idx_vars.size),
    np.array(df_main)[0, idx_vars],
    color="red",
    edgecolor="black",
    alpha=0.5,
)
plt.ylim([-0.15, 0.15])
plt.grid(color="gray", linestyle="dashed", linewidth=0.5)
plt.yticks(fontsize="9")
plt.xticks([0, 4, 9, 14, 19], [1, 5, 10, 15, 20], fontsize="9")
plt.xlabel("scale")
a.spines["right"].set_visible(False)
a.spines["top"].set_visible(False)
plt.locator_params(axis="y", nbins=5)
plt.savefig(
    dir_entropy / "mse_gfp_RAPM.jpg", format="jpg", dpi=1200, bbox_inches="tight"
)

# MSE and intelligence
idx_vars = []
for ch in range(59):
    for sc in range(20):
        var_name = "mse_ch_" + channel_names[ch] + "_sc_" + str(sc)
        idx_vars.append(np.where(column_names == var_name)[0][0])
idx_vars = np.array(idx_vars)

mse = np.array(df_main)[0, idx_vars]
mse = np.reshape(mse, (59, 20))

fig, a = plt.subplots(figsize=(18, 6))
divnorm = colors.TwoSlopeNorm(vmin=-0.31, vcenter=0, vmax=0.36)
plt.imshow(mse, cmap="bwr", norm=divnorm)
plt.colorbar()
plt.yticks(np.arange(mse.shape[0]), np.array(channel_names), fontsize="7")
plt.xticks([0, 4, 9, 14, 19], [1, 5, 10, 15, 20], fontsize="7")
# plt.grid(color='gray', linestyle='dashed', linewidth=0.5)
a.spines["right"].set_visible(False)
a.spines["top"].set_visible(False)
plt.savefig(
    dir_entropy / "mse_main_RAPM.jpg", format="jpg", dpi=1200, bbox_inches="tight"
)

# plot map of significant cluster MSE
p_adj = np.array(df_main)[2, idx_vars]
p_adj = np.reshape(p_adj, (59, 20))
p_sig = np.zeros((59, 20))
idx, idy = np.where(p_adj < 0.05)
p_sig[idx, idy] = 1
fig, a = plt.subplots()
plt.imshow(p_sig, cmap="bwr", norm=divnorm)

# %% Means and standard deviations of measures

# Shannon entropy
idx_vars = []
for ch in channel_names:
    var_name = "shannon_" + ch
    idx_vars.append(np.where(column_names == var_name)[0][0])
idx_vars = np.array(idx_vars)

fig, (ax1) = plt.subplots(ncols=1)
im, cm = mne.viz.plot_topomap(
    np.array(df_main)[3, idx_vars],
    raw.info,
    axes=ax1,
    show=False,
    vlim=(0.1, 0.9),
    contours=12,
)
cbar = fig.colorbar(im)
cbar.ax.tick_params(labelsize=20)
ax1.spines["right"].set_visible(False)
ax1.spines["top"].set_visible(False)
# plt.savefig("shannon_mean.jpg", format="jpg", dpi=1200, bbox_inches="tight")

fig, (ax1) = plt.subplots(ncols=1)
im, cm = mne.viz.plot_topomap(
    np.array(df_main)[4, idx_vars], raw.info, axes=ax1, show=False, vlim=(0.1, 0.25)
)
cbar = fig.colorbar(im)
cbar.ax.tick_params(labelsize=20)
ax1.spines["right"].set_visible(False)
ax1.spines["top"].set_visible(False)
# plt.savefig("shannon_SD.jpg", format="jpg", dpi=1200, bbox_inches="tight")


# Fuzzy entropy
idx_vars = []
for ch in channel_names:
    var_name = "fuzzy_" + ch
    idx_vars.append(np.where(column_names == var_name)[0][0])
idx_vars = np.array(idx_vars)

fig, (ax1) = plt.subplots(ncols=1)
im, cm = mne.viz.plot_topomap(
    np.array(df_main)[3, idx_vars],
    raw.info,
    axes=ax1,
    show=False,
    vlim=(0.1, 0.9),
    contours=12,
)
cbar = fig.colorbar(im)
cbar.ax.tick_params(labelsize=20)
ax1.spines["right"].set_visible(False)
ax1.spines["top"].set_visible(False)
# plt.savefig("fuzzy_mean.jpg", format="jpg", dpi=1200, bbox_inches="tight")

fig, (ax1) = plt.subplots(ncols=1)
im, cm = mne.viz.plot_topomap(
    np.array(df_main)[4, idx_vars],
    raw.info,
    axes=ax1,
    show=False,
    vlim=(0.1, 0.25),
    contours=12,
)
cbar = fig.colorbar(im)
cbar.ax.tick_params(labelsize=20)
ax1.spines["right"].set_visible(False)
ax1.spines["top"].set_visible(False)
# plt.savefig("fuzzy_SD.jpg", format="jpg", dpi=1200, bbox_inches="tight")


# MSE of GFP
idx_vars = []
for sc in range(20):
    var_name = "mse_gfp_scale_" + str(sc)
    idx_vars.append(np.where(column_names == var_name)[0][0])
idx_vars = np.array(idx_vars)

fig, ax = plt.subplots()
plt.bar(
    range(idx_vars.size),
    np.array(df_main)[3, idx_vars],
    color="red",
    edgecolor="black",
    alpha=0.5,
)
# plt.ylim([-0.25, 0.1])
plt.grid(color="gray", linestyle="dashed", linewidth=0.5)
plt.yticks(fontsize="20")
plt.xticks([0, 4, 9, 14, 19], [1, 5, 10, 15, 20], fontsize="20")
plt.xlabel("scale", fontsize="20")
plt.ylabel("mean MSE", fontsize="20")
plt.ylim((0, 0.9))
plt.locator_params(axis="y", nbins=5)
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
# plt.savefig("mse_gfp_mean.jpg", format="jpg", dpi=1200, bbox_inches="tight")

fig, ax = plt.subplots()
plt.bar(
    range(idx_vars.size),
    np.array(df_main)[4, idx_vars],
    color="red",
    edgecolor="black",
    alpha=0.5,
)
# plt.ylim([-0.25, 0.1])
plt.grid(color="gray", linestyle="dashed", linewidth=0.5)
plt.yticks(fontsize="20")
plt.xticks([0, 4, 9, 14, 19], [1, 5, 10, 15, 20], fontsize="20")
plt.xlabel("scale", fontsize="20")
plt.ylabel("SD MSE", fontsize="20")
plt.ylim((0, 0.25))
plt.locator_params(axis="y", nbins=4)
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
# plt.savefig("mse_gfp_SD.jpg", format="jpg", dpi=1200, bbox_inches="tight")


# MSE
idx_vars = []
for ch in range(59):
    for sc in range(20):
        var_name = "mse_ch_" + channel_names[ch] + "_sc_" + str(sc)
        idx_vars.append(np.where(column_names == var_name)[0][0])
idx_vars = np.array(idx_vars)

var = np.array(df_main)[3, idx_vars]
var = np.reshape(var, (59, 20))

scale = np.array([0, 4, 9, 14, 19])
fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(ncols=5)
vmin = 0.1
vmax = 0.9
im, cm = mne.viz.plot_topomap(
    var[:, scale[0]], raw.info, axes=ax1, show=False, vlim=(vmin, vmax)
)
im, cm = mne.viz.plot_topomap(
    var[:, scale[1]], raw.info, axes=ax2, show=False, vlim=(vmin, vmax)
)
im, cm = mne.viz.plot_topomap(
    var[:, scale[2]], raw.info, axes=ax3, show=False, vlim=(vmin, vmax)
)
im, cm = mne.viz.plot_topomap(
    var[:, scale[3]], raw.info, axes=ax4, show=False, vlim=(vmin, vmax)
)
im, cm = mne.viz.plot_topomap(
    var[:, scale[4]], raw.info, axes=ax5, show=False, vlim=(vmin, vmax)
)
cbar_ax = fig.add_axes([0.21, 0.35, 0.6, 0.05])
cb = fig.colorbar(im, cax=cbar_ax, orientation="horizontal")
cb.ax.tick_params(labelsize=12)
tick_locator = ticker.MaxNLocator(nbins=5)
cb.locator = tick_locator
cb.update_ticks()
# plt.savefig("MSE_mean.jpg", format="jpg", dpi=1200, bbox_inches="tight")


var = np.array(df_main)[4, idx_vars]
var = np.reshape(var, (59, 20))

scale = np.array([0, 4, 9, 14, 19])
fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(ncols=5)
vmin = 0.1
vmax = 0.25
im, cm = mne.viz.plot_topomap(
    var[:, scale[0]], raw.info, axes=ax1, show=False, vlim=(vmin, vmax)
)
im, cm = mne.viz.plot_topomap(
    var[:, scale[1]], raw.info, axes=ax2, show=False, vlim=(vmin, vmax)
)
im, cm = mne.viz.plot_topomap(
    var[:, scale[2]], raw.info, axes=ax3, show=False, vlim=(vmin, vmax)
)
im, cm = mne.viz.plot_topomap(
    var[:, scale[3]], raw.info, axes=ax4, show=False, vlim=(vmin, vmax)
)
im, cm = mne.viz.plot_topomap(
    var[:, scale[4]], raw.info, axes=ax5, show=False, vlim=(vmin, vmax)
)
cbar_ax = fig.add_axes([0.21, 0.35, 0.6, 0.05])
cb = fig.colorbar(im, cax=cbar_ax, orientation="horizontal")
cb.ax.tick_params(labelsize=12)
tick_locator = ticker.MaxNLocator(nbins=7)
cb.locator = tick_locator
cb.update_ticks()
# plt.savefig("MSE_SD.jpg", format="jpg", dpi=1200, bbox_inches="tight")
