"""
This script calculates correlations between entropy measures and fluid intelligence.
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

import numpy as np
import pandas as pd
import pingouin as pg
from scipy.ndimage import measurements
from sklearn.preprocessing import minmax_scale
from statsmodels.stats.multitest import fdrcorrection

# %% Load data
cwd = Path.cwd()
dir_root = cwd.parent.parent
dir_data = dir_root / "Data"
dir_entropy = dir_data / "Entropy_Thiele"
dir_microstates = dir_data / "Microstates"

gf_table = pd.read_csv(
    dir_data / "IST_table.csv", delimiter=";", encoding="unicode_escape"
)
gf_table = gf_table.assign(
    gf_score=lambda gf_table: np.sum(gf_table.iloc[:, 1:], axis=1)
)

age_table = pd.read_csv(
    dir_data / "SocioDemographics.txt", delimiter=",", encoding="unicode_escape"
)
age_table = age_table[["ID", "Age"]]

dv_table = pd.merge(gf_table[["ID", "gf_score"]], age_table, on="ID", how="inner")
dv_table = dv_table[dv_table["gf_score"].notna()]
dv_table = dv_table[dv_table["Age"].notna()]

intell_main = dv_table["gf_score"].to_numpy()  # gf scores
age_main = dv_table["Age"].to_numpy()  # age
subs_beh = dv_table["ID"].to_numpy()  # subject IDs corresponding to behavioral data

df_complexity_run2_EO = pd.read_csv(
    dir_entropy / "df_entropy_first_run_eyes_open.csv"
)  # load complexity measures
subs_neuro_run2_EO = np.array(
    df_complexity_run2_EO.ID
)  # subject IDs corresponding to complexity data

# channel names
with open(dir_microstates / "channel_names", "rb") as fp:
    channel_names = pickle.load(fp)

# %% Sort behavioral data according to order of subjects in neuro data

# keep only subjects that exist in both the neuro and behavioral datasets
behavior_ids = set(subs_beh.tolist())
keep_mask = np.array([s in behavior_ids for s in subs_neuro_run2_EO], dtype=bool)
df_complexity = df_complexity_run2_EO.loc[keep_mask].copy()
subs_neuro_run2_EO = np.array(df_complexity.ID)

# get indexes where subjects from behavioral data are found in subjects of neuro data
idx_neuro_run2_EO = []
for s in list(subs_neuro_run2_EO):

    idx = np.where(subs_beh == s)[0]
    if idx.shape[0] == 0:
        raise ValueError(f"No behavioral data found for subject {s}")
    idx_neuro_run2_EO.append(int(idx[0]))

idx_neuro_run2_EO = np.array(idx_neuro_run2_EO)

# sort age and intelligence (RAPM) scores according to order of neuro data
age_sorted = age_main[idx_neuro_run2_EO]
intel_sorted = intell_main[idx_neuro_run2_EO]

# %% Compute associations between complexity values and RAPM scores


def get_association(X, y, C):

    C = np.array(C)
    X = minmax_scale(
        X, feature_range=(0, 1)
    )  # scaling between 0 and 1 for better comparision between mean and SD
    data = np.vstack((X, y, C)).T

    if C.shape[0] == 1:
        df_data = pd.DataFrame(data, columns=["X", "y", "conf1"])
        res = pg.partial_corr(
            data=df_data, x="X", y="y", covar=["conf1"]
        )  # Pearson partial correlation

    elif C.shape[0] == 2:
        df_data = pd.DataFrame(data, columns=["X", "y", "conf1", "conf2"])
        res = pg.partial_corr(
            data=df_data, x="X", y="y", covar=["conf1", "conf2"]
        )  # Pearson partial correlation

    elif C.shape[0] == 3:
        df_data = pd.DataFrame(data, columns=["X", "y", "conf1", "conf2", "conf3"])
        res = pg.partial_corr(
            data=df_data, x="X", y="y", covar=["conf1", "conf2", "conf3"]
        )
    else:
        raise ValueError(
            "get_associations: not implemented for your number of confounds"
        )

    return res, X


C = [age_sorted]  # confounds
y = intel_sorted  # RAPM scores

# lists to save results in
results_p_raw = []
results_p_adj = []
results_corr = []
results_mean = []
results_std = []
results_min = []
results_max = []
results_names = []

# association between Shannon entropy of GFP signals and RAPM scores
idx_cols_vars = np.where(df_complexity.columns == "shannon_gfp")[0]
X = np.array(df_complexity.iloc[:, idx_cols_vars]).ravel()
res, X = get_association(X, y, C)
p = res["p_val"].to_numpy()
corr = res["r"].to_numpy()
results_p_raw.append(p)
results_p_adj.append(p)
results_corr.append(corr)
results_mean.append(np.mean(X))
results_std.append(np.std(X))
results_min.append(np.min(X))
results_max.append(np.max(X))
results_names.append("shannon_gfp")

# association betweem Fuzzy entropy of GFP signals and RAPM scores
idx_cols_vars = np.where(df_complexity.columns == "fuzzy_gfp")[0]
X = np.array(df_complexity.iloc[:, idx_cols_vars]).ravel()
res, X = get_association(X, y, C)
p = res["p_val"].to_numpy()
corr = res["r"].to_numpy()
results_p_raw.append(p)
results_p_adj.append(p)
results_corr.append(corr)
results_mean.append(np.mean(X))
results_std.append(np.std(X))
results_min.append(np.min(X))
results_max.append(np.max(X))
results_names.append("fuzzy_gfp")


# assosiation between Shannon entropy of each channel and RAPM scores
corr_ch = np.zeros((59, 1))
p_ch = np.zeros((59, 1))
for ch in range(59):

    idx_cols_vars = np.where(
        df_complexity.columns == "shannon_ch_" + str(channel_names[ch])
    )[0]
    X = np.array(df_complexity.iloc[:, idx_cols_vars]).ravel()
    res, X = get_association(X, y, C)
    p_ch[ch] = res["p_val"].to_numpy()
    corr_ch[ch] = res["r"].to_numpy()
    results_mean.append(np.mean(X))
    results_std.append(np.std(X))
    results_min.append(np.min(X))
    results_max.append(np.max(X))
    results_names.append("shannon_" + channel_names[ch])

reject, p_val_adj = fdrcorrection(p_ch.flatten(), alpha=0.05)
results_corr.append(corr_ch.flatten())
results_p_raw.append(p_ch.flatten())
results_p_adj.append(p_val_adj)


# associations between Fuzzy entropy of each channel and RAPM scores
corr_ch = np.zeros((59, 1))
p_ch = np.zeros((59, 1))
for ch in range(59):

    idx_cols_vars = np.where(
        df_complexity.columns == "fuzzy_ch_" + str(channel_names[ch])
    )[0]
    X = np.array(df_complexity.iloc[:, idx_cols_vars]).ravel()
    res, X = get_association(X, y, C)
    p_ch[ch] = res["p_val"].to_numpy()
    corr_ch[ch] = res["r"].to_numpy()
    results_mean.append(np.mean(X))
    results_std.append(np.std(X))
    results_min.append(np.min(X))
    results_max.append(np.max(X))
    results_names.append("fuzzy_" + channel_names[ch])

reject, p_val_adj = fdrcorrection(p_ch.flatten(), alpha=0.05)
results_corr.append(corr_ch.flatten())
results_p_raw.append(p_ch.flatten())
results_p_adj.append(p_val_adj)


# associations between MSE and RAPM scores
corr_ch_sc = np.zeros((59, 20))
p_ch_sc = np.zeros((59, 20))
for ch in range(59):
    for sc in range(20):

        idx_cols_vars = np.where(
            df_complexity.columns == "MSE_" + channel_names[ch] + "_" + str(sc)
        )[0]
        X = np.array(df_complexity.iloc[:, idx_cols_vars]).ravel()
        res, X = get_association(X, y, C)

        p_ch_sc[ch, sc] = res["p_val"].to_numpy()
        corr_ch_sc[ch, sc] = res["r"].to_numpy()

        results_mean.append(np.mean(X))
        results_std.append(np.std(X))
        results_min.append(np.min(X))
        results_max.append(np.max(X))
        results_names.append("mse_ch_" + channel_names[ch] + "_sc_" + str(sc))


results_corr.append(corr_ch_sc.flatten())
results_p_raw.append(p_ch_sc.flatten())

# measure clustersize of significant correlations
signi_ch_sc = np.zeros((59, 20))
signi_ch_sc[np.where(p_ch_sc < 0.05)] = 1
lw_real, num_real = measurements.label(signi_ch_sc)
area_real = measurements.sum(
    signi_ch_sc, lw_real, index=np.arange(1, lw_real.max() + 1)
)  # clustersizes

# permutation test (null models) for computing p-value of clustersizes in real data
areas_all = []
n_permutations = 100
for p in range(n_permutations):

    print(p)
    corr_ch_sc = np.zeros((59, 20))
    p_ch_sc = np.zeros((59, 20))
    intel_sorted_perm = np.random.permutation(intel_sorted)  # Permuted RAPM scores

    for ch in range(59):
        for sc in range(20):

            idx_cols_vars = np.where(
                df_complexity.columns == "MSE_" + channel_names[ch] + "_" + str(sc)
            )[0]
            X = np.array(df_complexity.iloc[:, idx_cols_vars]).ravel()
            res, X = get_association(X, intel_sorted_perm, C)

            p_ch_sc[ch, sc] = res["p_val"].to_numpy()
            corr_ch_sc[ch, sc] = res["r"].to_numpy()

    signi_ch_sc = np.zeros((59, 20))
    signi_ch_sc[np.where(p_ch_sc < 0.05)] = 1
    lw, num = measurements.label(signi_ch_sc)
    area = measurements.sum(signi_ch_sc, lw, index=np.arange(1, lw.max() + 1))
    areas_all.append(area)

# cluster sizes of permutated data
size_cluster_perm = np.concatenate(areas_all).astype(int)
count_cluster_sizes = np.bincount(size_cluster_perm)

CI95 = np.percentile(
    size_cluster_perm, 95
)  # 95th percentile of cluster sizes of null models

# compute p-values for each cluster size of real data
p_cluster_size = []
for n in range(1, int(area_real.max()) + 1):
    permute = size_cluster_perm
    real = n
    sum_bad = np.sum(
        permute >= real
    )  # instances of null models have larger cluster size than real data
    p = sum_bad / permute.shape[0]
    p_cluster_size.append(p)

p_cluster_size = np.array(p_cluster_size)

# assign adjusted p-values to clusters of real data
p_adj = np.ones(p_ch_sc.shape)
for c in range(num_real):

    c_idx = c + 1
    area_c = area_real[c]
    p_c = p_cluster_size[int(area_c) - 1]
    p_adj[np.where(lw_real == c_idx)] = p_c

p_val_adj = np.array(p_adj).flatten()
results_p_adj.append(p_val_adj)


# associations between clustered MSE and RAPM scores
n_cluster = 8
n_scale = 4
corr_ch_sc = np.zeros((n_cluster, n_scale))
p_ch_sc = np.zeros((n_cluster, n_scale))
for ch in range(n_cluster):
    for sc in range(n_scale):

        idx_cols_vars = np.where(
            df_complexity.columns
            == "MSE_cluster_channel_" + str(ch) + "_scales_" + str(sc)
        )[0]
        X = np.array(df_complexity.iloc[:, idx_cols_vars]).ravel()
        res, X = get_association(X, y, C)

        p_ch_sc[ch, sc] = res["p_val"].to_numpy()
        corr_ch_sc[ch, sc] = res["r"].to_numpy()

        results_mean.append(np.mean(X))
        results_std.append(np.std(X))
        results_min.append(np.min(X))
        results_max.append(np.max(X))
        results_names.append("mse_clustered_cluster_" + str(ch) + "_" + str(sc))


results_corr.append(corr_ch_sc.flatten())
results_p_raw.append(p_ch_sc.flatten())
reject, p_val_adj = fdrcorrection(p_ch_sc.flatten(), alpha=0.05)
results_p_adj.append(p_val_adj)


# associations between MSE of GFP and RAPM scores
rho_sc = np.zeros((20))
p_sc = np.zeros((20))
for sc in range(20):

    idx_cols_vars = np.where(df_complexity.columns == "MSE_gfp_scale" + str(sc))[0]
    X = np.array(df_complexity.iloc[:, idx_cols_vars]).ravel()
    res, X = get_association(X, y, C)

    p_sc[sc] = res["p_val"].to_numpy()
    rho_sc[sc] = res["r"].to_numpy()
    results_mean.append(np.mean(X))
    results_std.append(np.std(X))
    results_min.append(np.min(X))
    results_max.append(np.max(X))
    results_names.append("mse_gfp_scale_" + str(sc))

reject, p_val_adj = fdrcorrection(p_sc.flatten(), alpha=0.05)
results_corr.append(rho_sc.flatten())
results_p_raw.append(p_sc.flatten())
results_p_adj.append(p_val_adj)

# combine all results in data frame and save
res1 = np.concatenate(results_corr)
res2 = np.concatenate(results_p_raw)
res3 = np.concatenate(results_p_adj)
res4 = np.array(results_mean)
res5 = np.array(results_std)
res6 = np.array(results_min)
res7 = np.array(results_max)

data = np.vstack((res1, res2, res3, res4, res5, res6, res7))
df_results_main = pd.DataFrame(
    data,
    columns=np.array(results_names),
    index=["corr", "p", "p_adj", "M", "SD", "min", "max"],
)

df_results_main.to_pickle(dir_entropy / "df_entropy_results_first_run_eyes_open.pkl")

# %%
