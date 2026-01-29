"""Decoy FDR calculation and validation."""

import os.path
import numpy as np
import pandas as pd
from novoboard.accuracy import WorkerTest


def read_denovo(denovo_csv, selected_features=None):
    """Read de novo sequencing results from CSV file."""
    denovo_psm = pd.read_csv(denovo_csv, keep_default_na=False)
    denovo_psm['feature_id'] = denovo_psm.apply(lambda row: row['Source File'].split('.mgf')[0] + '.mgf' + '||' + str(row['Scan']), axis=1)
    if selected_features:
        denovo_psm['selected_features'] = denovo_psm.apply(lambda row: row['feature_id'] in selected_features, axis=1)
        denovo_psm = denovo_psm[denovo_psm['selected_features']]
    return denovo_psm


def calculate_FDR(target_csv, decoy_csv, engine_score, fdr_list, selected_features=None):
    """Calculate FDR using target-decoy approach."""
    print("target_csv =", target_csv)
    print("decoy_csv =", decoy_csv)
    target_psm = read_denovo(target_csv, selected_features)
    decoy_psm = read_denovo(decoy_csv, selected_features)
    print("len(target_psm) =", len(target_psm)); print("len(decoy_psm) =", len(decoy_psm))
    dfs = pd.concat([target_psm, decoy_psm], keys=['target', 'decoy']).reset_index().rename(columns={'level_0': 'spectrum'})
    dfs['is_target'] = dfs.apply(lambda row: row['spectrum']=='target', axis=1)

    # target-decoy competition
    dfs.sort_values(by=[engine_score, 'is_target'], ascending=[False, False], inplace=True)
    # 1-1 competition on each scan id
#     dfs_fdr = dfs.drop_duplicates(subset=['feature_id'])
    # competition on whole dataset
    dfs_fdr = dfs
    dfs_fdr['feature_id'] = dfs_fdr.apply(lambda x: x['feature_id'] if x['is_target'] else x['feature_id']+'||decoy', axis=1)
    print("len(dfs) =", len(dfs))
    print("len(dfs_fdr) =", len(dfs_fdr))
    print("sum(dfs_fdr['is_target']) =", sum(dfs_fdr['is_target']))
    
    # fdr estimation
    cumsum = range(1, len(dfs_fdr) + 1)
    cumsum_target = np.cumsum(np.array(dfs_fdr['is_target'].astype(int)))
    cumsum_decoy = cumsum - cumsum_target
    estimated_fdr = cumsum_decoy / cumsum_target
    dfs_fdr['estimated_fdr'] = estimated_fdr

    score_list = []
    count_list = []
    for fdr in fdr_list:
        fdr_index = np.flatnonzero(estimated_fdr <= fdr)
        fdr_index = fdr_index[-1] if len(fdr_index) > 0 else 0
        score_list.append(dfs_fdr.iloc[fdr_index][engine_score])
        count_list.append(fdr_index + 1)

    return dfs, dfs_fdr, score_list, count_list


def validate_FDR(target_csv, decoy_csv, engine_score, db_csv, spectrum_file, p_decoy, T_pct, col_score, col_aa_score):
    """Validate FDR estimation against known database matches."""
    db_psm = pd.read_csv(db_csv, keep_default_na=False)
    db_psm['feature_id'] = db_psm.apply(lambda row: row['Source File'].split('.mgf')[0] + '.mgf' + '||' + str(row['Scan']), axis=1)
    selected_features = None # set(db_psm['feature_id'])

    dfs, dfs_fdr, score_list, count_list = calculate_FDR(target_csv, decoy_csv, engine_score, p_decoy, selected_features)

    target_decoy_csv = target_csv + '-' + decoy_csv.split('/')[-1]
    dfs_fdr.to_csv(target_decoy_csv, index=False)
    accuracy_file = target_decoy_csv + '.accuracy'
    if not os.path.isfile(accuracy_file):
        worker_test = WorkerTest(db_csv, target_decoy_csv, spectrum_file, col_score, col_aa_score)
        worker_test.test_accuracy()

    denovo_df = dfs_fdr
    denovo_df = denovo_df.set_index('feature_id')
    accuracy_df = pd.read_csv(accuracy_file, delimiter ='\t', index_col='feature_id')
    denovo_df['db_peptide'] = accuracy_df.apply(
        lambda x: x['target_sequence'].replace('C(Carbamidomethylation)', 'C(+57.02)').replace('M(Oxidation)', 'M(+15.99)').replace(',', ''), 
        axis=1)
    denovo_df['recall_AA'] = accuracy_df['recall_AA']
    denovo_df['predicted_len'] = accuracy_df['predicted_len']
    denovo_df['recall_peptide'] = accuracy_df.apply(lambda x: x['recall_AA'] == x['predicted_len'], axis=1)
    denovo_df['target_ion'] = accuracy_df['target_ion']
    denovo_df['matched_ion'] = accuracy_df['matched_ion']
    denovo_df['recall_peptide_I'] = accuracy_df.apply(lambda x: x['matched_ion'] == x['target_ion'], axis=1)
    denovo_df['recall_peptide_T'] = accuracy_df.apply(lambda x: x['matched_ion'] >= x['target_ion']*T_pct, axis=1)
    print("len(denovo_df) =", len(denovo_df))
    print("  with recall_AA =", len(denovo_df[~denovo_df['recall_AA'].isna()]))
    print("    is_target =", len(denovo_df[~denovo_df['recall_AA'].isna()][denovo_df['is_target']]))

    # calculate true FDR on annotated target spectra
    df = denovo_df
    df = df[~df['recall_AA'].isna()]
    df = df[df['is_target']]
    cumsum = range(1, len(df) + 1)
    cumsum_correct = np.cumsum(np.array(df['recall_peptide'].astype(int)))
    cumsum_false = cumsum - cumsum_correct
    true_fdr = cumsum_false / cumsum
    cumsum_correct = np.cumsum(np.array(df['recall_peptide_I'].astype(int)))
    cumsum_false = cumsum - cumsum_correct
    true_fdr_I = cumsum_false / cumsum
    cumsum_correct = np.cumsum(np.array(df['recall_peptide_T'].astype(int)))
    cumsum_false = cumsum - cumsum_correct
    true_fdr_T = cumsum_false / cumsum
    # only report entries with decreasing fdr from the bottom to avoid bumps
    min_est, min_true = 1, 1
    reported = []
    for x, y, z, v, w in list(zip(df['estimated_fdr'], cumsum, true_fdr, true_fdr_I, true_fdr_T))[::-1]:
        if x <= min_est and w <= min_true:
            min_est = x
            min_true = w
            reported.append((x, y, z, v, w))
    estimated_fdr, cumsum, true_fdr, true_fdr_I, true_fdr_T = zip(*reported)
    
    # calculate estimated FDR and #PSMs on all target spectra
    df = denovo_df
    df = df[df['is_target']]
    cumsum_full = range(1, len(df) + 1)
    # only report entries with decreasing fdr from the bottom to avoid bumps
    min_est = 1
    reported = []
    for x, y in list(zip(df['estimated_fdr'], cumsum_full))[::-1]:
        if x <= min_est:
            min_est = x
            reported.append((x, y))
    estimated_fdr_full, cumsum_full = zip(*reported)
    
    return {
        'denovo_df': denovo_df, 
        'df': df, 
        'estimated_fdr': estimated_fdr,
        'cumsum': cumsum,
        'true_fdr': true_fdr,
        'true_fdr_I': true_fdr_I,
        'true_fdr_T': true_fdr_T,
        'estimated_fdr_full': estimated_fdr_full,
        'cumsum_full': cumsum_full,
    }

