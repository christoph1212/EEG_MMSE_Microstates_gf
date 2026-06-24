"""
This script calculates entropy measures from resting state EEG data.
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

import pickle
import traceback
from pathlib import Path

import mne
import neurokit2 as nk
import numpy as np
import pandas as pd

# https://github.com/ixjlyons/entro-py/tree/b371cbc5f72197ca9f59506266de8cb927a2f0d8
import entropy

# %% Define directories
cwd = Path.cwd()
dir_root = cwd.parent.parent
dir_data = dir_root / "Data"
dir_preprocessed = dir_data / "Snipplet"  # Epoched data
dir_entropy = dir_data / "Entropy_Thiele"
dir_log = dir_data / "Log/Entropy_Thiele"

if not dir_entropy.exists():
    dir_entropy.mkdir(parents=True, exist_ok=True)

if not dir_log.exists():
    dir_log.mkdir(parents=True, exist_ok=True)

# read names of files
files = list(dir_preprocessed.rglob("*.set"))

print(f"Loading {len(files)} files for microstate analysis\n")

# initialize lists
IDs_subjects = []
conditions = []
epoch_lengths = []
data_all = []

# loop over all subjects
for f in files:

    # get filename
    fname = f.name

    epoch_length = fname.split("_")[0]
    if epoch_length != "40":
        continue
    epoch_lengths.append(epoch_length)

    # store info in lists
    ID = fname.split("_")[6]
    IDs_subjects.append(ID)

    condition = [
        fname.split("_")[2]
        + "_"
        + fname.split("_")[3]
        + "_"
        + fname.split("_")[4]
        + "_"
        + fname.split("_")[5]
    ]
    conditions.append(condition)

    # read file
    preproc = mne.io.read_raw_eeglab(dir_preprocessed / fname, preload=True)

    # remap non-eeg channels
    data = preproc.get_data()
    data_all.append(data)

epoch_lengths = np.array(epoch_lengths).astype(int).tolist()
conditions = np.array(conditions).flatten()

# %% Compute entropies

for cond in list(set(conditions)):
    try:
        # select subset data
        indices = np.where(conditions == cond)[0]
        data_subset = [data_all[i] for i in indices]

        # Initial checks
        if len(data_subset) != len(np.array(IDs_subjects)[indices]):
            raise ValueError(
                "Mismatch between number of subjects and data subset length."
            )
        if any(np.array(epoch_lengths)[indices] != 40):
            raise ValueError("All epochs must have a length of 40 seconds.")

        print(f"\n***Computing microstates for condition: {cond}***\n")

        print(f"Number of subjects in this condition: {len(data_subset)}")

        # MSE on GFP signal
        gfp_all = []
        MSE_gfp_all = []
        for data in data_subset:

            # calculate whole brain MSE on GFP
            gfp = np.std(data, axis=0)
            gfp_all.append(gfp)
            MSE_gfp, info_MSE_gfp = nk.entropy_multiscale(gfp, show=False, scale=21)
            MSE_gfp_vals = info_MSE_gfp["Values"]
            MSE_gfp_all.append(MSE_gfp_vals)

        # MSE for each channel
        MSE_all = []
        for data in data_subset:
            MSE = []
            # loop over channels and calculate MSE
            for channel in range():

                MSE_auc, info_MSE = nk.entropy_multiscale(
                    data[channel, :], show=False, scale=21
                )
                MSE_vals = info_MSE["Values"]
                MSE.append(MSE_vals)

            MSE_all.append(np.array(MSE))

        # Fuzzy entropy of GFP
        divident = 16
        fuzzy_gfp_all = []
        for s in range(len(gfp_all)):

            fe_all_div = []
            gfp_signal = gfp_all[s]

            len_signal = int(np.ceil(gfp_signal.shape[0] / divident))
            for div in range(divident):
                signal_div = gfp_signal[
                    div * len_signal : div * len_signal + len_signal
                ]
                fe_div = entropy.fuzzyen(signal_div, dim=2, r=0.2, n=1)
                fe_all_div.append(fe_div)

            fuzzy_gfp_all.append(np.mean(np.array(fe_all_div)))

        # Fuzzy entropy all channels
        fuzzy_all = []
        for data in data_subset:

            len_signal = int(np.ceil(data.shape[1] / divident))

            fe_ch_all = []

            # loop over channels
            for channel in range(28):
                fe_ch = []

                for div in range(divident):
                    signal_div = data[
                        channel, div * len_signal : div * len_signal + len_signal
                    ]
                    fe = entropy.fuzzyen(signal_div, dim=2, r=0.2, n=1)
                    fe_ch.append(fe)
                fe_ch_avg = np.mean(np.array(fe_ch))
                fe_ch_all.append(fe_ch_avg)

            fuzzy_all.append(np.array(fe_ch_all))

        # Shannon entropy of GFP
        shannon_gfp_all = []
        for gfp_s in gfp_all:

            data_rounded = np.around(gfp_s, decimals=7)
            shannon_gfp = nk.entropy_shannon(data_rounded)
            shannon_gfp_all.append(shannon_gfp[0])

        # Shannon entropy of each channel
        shannon_all = []
        for data in data_subset:

            data_rounded = np.around(data, decimals=7)
            se_channel = []
            # loop over channels and calculate se
            for channel in range(28):
                se = nk.entropy_shannon(data_rounded[channel, :])
                se_channel.append(se[0])

            shannon_all.append(np.array(se_channel))

        # %% Save variables in tables

        channel_names = preproc.info.ch_names

        # Shannon entropy of GFP
        variables = np.array(shannon_gfp_all)
        table_shannon_gfp = pd.DataFrame(variables, columns=["shannon_gfp"])

        # Fuzzy entropy of GFP
        variables = np.array(fuzzy_gfp_all)
        table_fuzzy_gfp = pd.DataFrame(variables, columns=["fuzzy_gfp"])

        # Shannon entropy of all channels
        variables = np.array(shannon_all)
        names = []
        for i in range(variables.shape[1]):
            names.append("shannon_ch_" + str(channel_names[i]))
        names = np.array(names)
        table_shannon = pd.DataFrame(variables, columns=names)

        # Fuzzy entropy of all channels
        variables = np.array(fuzzy_all)
        names = []
        for i in range(variables.shape[1]):
            names.append("fuzzy_ch_" + str(channel_names[i]))
        names = np.array(names)
        table_fuzzy = pd.DataFrame(variables, columns=names)

        # MSE all channels all scales
        variables = np.array(MSE_all)
        column_names = []
        for ch in range(variables.shape[1]):
            names = []
            for sc in range(variables.shape[2]):
                names.append("MSE_" + channel_names[ch] + "_" + str(sc))
            column_names.append(np.array(names))
        column_names = np.array(column_names)

        data = np.reshape(
            variables, (variables.shape[0], variables.shape[1] * variables.shape[2])
        )
        columns = np.reshape(column_names, (variables.shape[1] * variables.shape[2]))
        table_MSE = pd.DataFrame(data, columns=columns)

        # MSE in spatial and temporal cluster for prediction models
        spatial_cluster = [
            ["Fp1", "Fp2"],
            ["FC1", "FC2", "FC5", "FC6"],
            ["F3", "F7", "Fz", "F4", "F8"],
            ["T7", "T8"],
            ["CP5", "CP1", "CP6", "CP2"],
            ["C3", "Cz", "C4"],
            ["P3", "P7", "P4", "P8", "Pz"],
            ["O1", "Oz", "O2"],
        ]
        scale_cluster = (
            np.arange(0, 5),
            np.arange(5, 10),
            np.arange(10, 15),
            np.arange(15, 20),
        )  # 4 clusters with 5 consecutive time steps each

        idx_spatial_cluster = []
        # get indexes for channels of spatial clusters
        for c in spatial_cluster:
            ind_c = []
            for e in c:

                ind_c.append(channel_names.index(e))

            idx_spatial_cluster.append(np.array(ind_c))

        MSE_cluster = []
        for s in range(len(MSE_all)):
            # averaging MSE within each spatial and temporal cluster
            MSE_s = MSE_all[s]
            MSE_cluster_s = []
            for idx in idx_spatial_cluster:
                MSE_cluster_ch = []
                for sc in scale_cluster:
                    Mx = MSE_s[idx, :]
                    Mx = Mx[:, sc]
                    MSE_cluster_ch.append(np.mean(Mx))

                MSE_cluster_s.append(np.array(MSE_cluster_ch))
            MSE_cluster.append(np.array(MSE_cluster_s))

        variables = np.array(MSE_cluster)
        column_names = []
        for ch in range(variables.shape[1]):
            names = []
            for sc in range(variables.shape[2]):
                names.append("MSE_cluster_channel_" + str(ch) + "_scales_" + str(sc))
            column_names.append(np.array(names))
        column_names = np.array(column_names)

        data = np.reshape(
            variables, (variables.shape[0], variables.shape[1] * variables.shape[2])
        )
        columns = np.reshape(column_names, (variables.shape[1] * variables.shape[2]))
        table_MSE_cluster = pd.DataFrame(data, columns=columns)

        # MSE of GFP
        variables = np.array(MSE_gfp_all)
        names = []
        for sc in range(variables.shape[1]):
            names.append("MSE_gfp_scale" + str(sc))
        names = np.array(names)
        table_MSE_gfp = pd.DataFrame(variables, columns=names)

        df_complexity = pd.concat(
            [
                pd.DataFrame(np.array(IDs_subjects)[indices], columns=["ID"]),
                pd.DataFrame(np.array(conditions)[indices], columns=["Condition"]),
                pd.DataFrame(np.array(epoch_lengths)[indices], columns=["Length"]),
                table_shannon_gfp,
                table_fuzzy_gfp,
                table_shannon,
                table_fuzzy,
                table_MSE,
                table_MSE_cluster,
                table_MSE_gfp,
            ],
            axis=1,
        )

        save_name = f"df_entropy_{cond}.pkl"
        df_complexity.to_pickle(dir_entropy / save_name)

    except Exception as e:
        print(f"An error occurred while processing condition {cond}: {e}")
        # Write log
        with open(dir_log / f"Error_entropy_{cond}.txt", "a") as log_file:
            log_file.write(f"Error processing condition {cond}: {e}\n")
            log_file.write(traceback.format_exc())

with open("channel_names", "wb") as fp:
    pickle.dump(channel_names, fp)

np.savetxt("channel_names.csv", np.array(channel_names), delimiter=",", fmt="%s")
