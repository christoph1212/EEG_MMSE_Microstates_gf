"""
This script calculates microstate measures from resting state EEG data.
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
from pathlib import Path

import mne
import numpy as np
import pandas as pd
from scipy.signal import find_peaks

import microstates_group
import microstates_subject

# %% Setup and load data

# Define directories
cwd = Path.cwd()
dir_root = cwd.parent.parent
dir_data = dir_root / "Data"
dir_preprocessed = dir_data / "Snipplet"  # Epoched data
dir_microstates = dir_data / "Microstates"

if not dir_microstates.exists():
    dir_microstates.mkdir(parents=True, exist_ok=True)

# read names of files
files = list(dir_preprocessed.rglob("*.set"))

print(f"Loading {len(files)} files for microstate analysis\n")

# initialize lists
IDs_subjects = []
runs = []
conditions = []
epoch_lengths = []
data_all = []

# loop over all subjects (filenames saved in filesEC)
for f in files:

    # get filename
    fname = f.name

    epoch_length = fname.split("_")[0]
    if epoch_length == "0":
        continue
    epoch_lengths.append(epoch_length)

    # store info in lists
    ID = fname.split("_")[6]
    IDs_subjects.append(ID)

    run = fname.split("_")[8]
    runs.append(run)

    condition = [fname.split("_")[4] + " " + fname.split("_")[5]]
    conditions.append(condition)

    # read file
    preproc = mne.io.read_raw_eeglab(dir_preprocessed / fname, preload=True)

    # remap non-eeg channels
    data = preproc.get_data()
    data_all.append(data)

epoch_lengths = np.array(epoch_lengths).astype(int).tolist()

# %% Compute microstates

# initialize lists
gfp_peaks_all = []
n_gfp_peaks_all = []
maps_subject_all = []
for data in data_all:

    # compute GFP and peaks of GFP
    gfp = np.std(data, axis=0)
    peaks, _ = find_peaks(gfp, distance=2)
    gfp_peaks_all.append(peaks)
    n_peaks = len(peaks)
    n_gfp_peaks_all.append(n_peaks)

    """
    get subject-specific microstates
        segment the data into 5 microstates that maximize GEV (GEV = global
        explained variance = percentage of total variance explained by a
        given microstate)

    only time points that are GFP peaks are used (within microstates_subject)
    """
    maps, segmentation = microstates_subject.segment(
        data, n_states=5, n_inits=1000, thresh=1e-10, max_iter=10000
    )
    maps_subject_all.append(maps)

maps_subjects_all_arr = np.array(maps_subject_all)

# reshape for clustering all subject-specific microstates
maps_subjects_all_arr = np.reshape(
    maps_subjects_all_arr,
    (
        (maps_subjects_all_arr.shape[0] * maps_subjects_all_arr.shape[1]),
        maps_subjects_all_arr.shape[2],
    ),
)
maps_subjects_all_arr = maps_subjects_all_arr.T

"""
get group specific microstates
    segment the data into 5 microstates that maximize GEV (GEV = global
    explained variance = percentage of total variance explained by a given
    microstate)
"""
maps_group, segmentation = microstates_group.segment(
    maps_subjects_all_arr, n_states=5, n_inits=1000, thresh=1e-10, max_iter=10000
)
np.save(dir_microstates / "microstates_group.npy", maps_group)  # plot group microstates

# plot maps_group
microstates_group.plot_maps(maps_group, preproc.info)

# backfitting to whole eeg signal
maps_group = maps_group.T
map_sequence_all = []
corr_backmap = []
for data in data_all:

    corr_subj = []
    for m in range(maps_group.shape[1]):
        map_corr = microstates_group._corr_vectors(
            data, maps_group[:, m].reshape(-1, 1)
        )
        corr_subj.append(map_corr)

    corr_subj = np.array(corr_subj)
    corr_subj_abs = abs(corr_subj)
    corr_backmap.append(corr_subj_abs)
    map_sequence = np.argmax(corr_subj_abs, axis=0)
    map_sequence_all.append(map_sequence)

# backfitting to GFP peaks of eeg signal only
map_sequence_peaks_all = []
for data, gfp_peaks_s in zip(data_all, gfp_peaks_all):
    corr_subj = []

    for m in range(maps_group.shape[1]):
        map_corr_ind = microstates_group._corr_vectors(
            data[:, gfp_peaks_s], maps_group[:, m].reshape(-1, 1)
        )
        corr_subj.append(map_corr_ind)

    corr_subj = np.array(corr_subj)
    corr_subj_abs = abs(corr_subj)
    map_sequence_peaks = np.argmax(corr_subj_abs, axis=0)
    map_sequence_peaks_all.append(map_sequence_peaks)

# %% Compute microstate measures

# compute coverage time per microstate all time points
n_total_occurences_states_all = []
coverage = []
for seq in map_sequence_all:
    n_total_occurences_states_all.append(np.bincount(seq))
    coverage.append(np.bincount(seq) / len(seq))

# compute coverage time per microstate at GFP peaks
n_total_occurences_states_peaks_all = []
coverage_peak = []
for seq in map_sequence_peaks_all:
    n_total_occurences_states_peaks_all.append(np.bincount(seq))
    coverage_peak.append(np.bincount(seq) / len(seq))

# compute occurences of microstates (times it is transitioned into a
# microstate (no matter its duration))

# occurences at all time points
states_single_all = []
for st in map_sequence_all:
    st = np.array(st)
    diff_st = np.diff(st)
    pos = np.where(diff_st != 0)
    pos = pos[0] + 1
    pos = np.insert(pos, 0, 0)
    states_single_all.append(st[pos])

n_single_occurences_states_all = []
for st in states_single_all:
    n_single_occurences_states_all.append(np.bincount(st))

# occurences at all GFP peaks
states_single_peaks_all = []
for st in map_sequence_peaks_all:
    st = np.array(st)
    diff_st = np.diff(st)
    pos = np.where(diff_st != 0)
    pos = pos[0] + 1
    pos = np.insert(pos, 0, 0)
    states_single_peaks_all.append(st[pos])

n_single_occurences_states_peaks_all = []
for st in states_single_peaks_all:
    n_single_occurences_states_peaks_all.append(np.bincount(st))

# frequency of microstates
frequency = []
for st, e in zip(n_single_occurences_states_all, epoch_lengths):
    freq = st / e
    frequency.append(freq)

# lifespan
lifespan = np.array(n_total_occurences_states_all) / np.array(
    n_single_occurences_states_all
)

# lifespan at GFP peaks
lifespan_peaks = np.array(n_total_occurences_states_peaks_all) / np.array(
    n_single_occurences_states_peaks_all
)

# Transition probabilities
# Function t_empirical used from:
# https://github.com/Frederic-vW/eeg_microstates/blob/78283c71fb82d80704ba954d29bf88b830f2e416/eeg_microstates3.py
# Copyright (c) 2017 Frederic von Wegner
# MIT License


# transition matrix
def t_empirical(data, n_clusters):
    T = np.zeros((n_clusters, n_clusters))
    n = len(data)
    for i in range(n - 1):
        T[data[i], data[i + 1]] += 1.0
    p_row = np.sum(T, axis=1)
    for i in range(n_clusters):
        if p_row[i] != 0.0:
            for j in range(n_clusters):
                T[i, j] /= p_row[i]  # normalize row sums to 1.0
    return T


# transition matrix of microstates all time points
trans_mat = []
for seq in map_sequence_all:

    trans_mat.append(t_empirical(seq, 5))

# transition matrix of microstates GFP peaks only
trans_mat_peak = []
for seq in map_sequence_peaks_all:

    trans_mat_peak.append(t_empirical(seq, 5))

# %% Post-hoc analyses microstate measures

# similarity of subject specific microstates
similarity_subject_microstates_all = []
for maps_subject in maps_subject_all:

    correlation_maps_subject = abs(np.corrcoef(maps_subject))
    mean_correlation_maps_subject = np.mean(correlation_maps_subject)
    similarity_subject_microstates_all.append(mean_correlation_maps_subject)

# explained variance of individual signals by subject-specific microstates
gev_subject_all = []
for data, maps_subject in zip(data_all, maps_subject_all):

    activation = maps_subject.dot(data)
    segmentation = np.argmax(np.abs(activation), axis=0)
    map_corr = microstates_group._corr_vectors(data, maps_subject[segmentation].T)
    gfp = np.std(data, axis=0)
    gfp_sum_sq = np.sum(gfp**2)
    gev = sum((gfp * map_corr) ** 2) / gfp_sum_sq
    gev_subject_all.append(gev)

# explained variance of individual signals by group states
gev_group_all = []
for data in data_all:

    activation = maps_group.T.dot(data)
    segmentation = np.argmax(np.abs(activation), axis=0)
    map_corr = microstates_group._corr_vectors(data, maps_group.T[segmentation].T)
    gfp = np.std(data, axis=0)
    gfp_sum_sq = np.sum(gfp**2)
    gev = sum((gfp * map_corr) ** 2) / gfp_sum_sq
    gev_group_all.append(gev)

# %% Save variables in tables

channel_names = preproc.info.ch_names

# number of GFP peaks
variables = np.array(n_gfp_peaks_all)
table_n_gfp_peaks = pd.DataFrame(variables, columns=["n_gfp_peaks"])

# coverage of micorsates
variables = np.array(coverage)
names = []
for i in range(variables.shape[1]):
    names.append("coverage_" + str(i))

names = np.array(names)
table_coverage = pd.DataFrame(variables, columns=names)

# lifespan of microstates
variables = np.array(lifespan)
names = []
for i in range(variables.shape[1]):
    names.append("lifespan_" + str(i))

names = np.array(names)
table_lifespan = pd.DataFrame(variables, columns=names)

# lifespan of microstates GFP peaks only
variables = np.array(lifespan_peaks)
names = []
for i in range(variables.shape[1]):
    names.append("lifespan_peaks_" + str(i))

names = np.array(names)
table_lifespan_peaks = pd.DataFrame(variables, columns=names)

# frequency of microstates
variables = np.array(frequency)
names = []
for i in range(variables.shape[1]):
    names.append("frequence_" + str(i))

names = np.array(names)
table_frequence = pd.DataFrame(variables, columns=names)

# transition probabilities of microstates
variables = np.array(trans_mat)
column_names = []
for ms1 in range(variables.shape[1]):
    names = []
    for ms2 in range(variables.shape[2]):
        names.append("transition_probability_" + str(ms1) + "_" + str(ms2))
    column_names.append(np.array(names))
column_names = np.array(column_names)

data = np.reshape(
    variables, (variables.shape[0], variables.shape[1] * variables.shape[2])
)
columns = np.reshape(column_names, (variables.shape[1] * variables.shape[2]))
table_transmat = pd.DataFrame(data, columns=columns)

# transition probabilities of microstates GFP peaks only
variables = np.array(trans_mat_peak)
column_names = []
for ms1 in range(variables.shape[1]):
    names = []
    for ms2 in range(variables.shape[2]):
        names.append("transition_probability_peaks_" + str(ms1) + "_" + str(ms2))
    column_names.append(np.array(names))
column_names = np.array(column_names)

data = np.reshape(
    variables, (variables.shape[0], variables.shape[1] * variables.shape[2])
)
columns = np.reshape(column_names, (variables.shape[1] * variables.shape[2]))
table_transmat_peak = pd.DataFrame(data, columns=columns)

# Post-hoc analyses
variable_1 = np.array(similarity_subject_microstates_all)
variable_2 = np.array(gev_subject_all)
variable_3 = np.array(gev_group_all)
variables = np.vstack((variable_1, variable_2, variable_3)).T
table_posthoc = pd.DataFrame(
    variables, columns=["similarity_subject_micosates", "gev_subject", "gev_group"]
)

df_microstate = pd.concat(
    [
        pd.DataFrame(IDs_subjects, columns=["ID"]),
        pd.DataFrame(runs, columns=["run"]),
        pd.DataFrame(conditions, columns=["condition"]),
        pd.DataFrame(epoch_lengths, columns=["epoch_length"]),
        table_n_gfp_peaks,
        table_coverage,
        table_lifespan,
        table_lifespan_peaks,
        table_frequence,
        table_transmat,
        table_transmat_peak,
        table_posthoc,
    ],
    axis=1,
)


df_microstate.to_pickle(dir_microstates / "df_microstate.pkl")

with open(dir_microstates / "channel_names", "wb") as fp:
    pickle.dump(channel_names, fp)

np.savetxt(
    dir_microstates / "channel_names.csv",
    np.array(channel_names),
    delimiter=",",
    fmt="%s",
)
