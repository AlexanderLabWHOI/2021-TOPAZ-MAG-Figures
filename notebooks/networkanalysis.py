# Python file to recreate correlation calculation.
import os
import sys
import pandas as pd
import numpy as np    
import networkx as nx
import scipy
from scipy import stats
from joblib import Parallel, delayed
import multiprocessing
num_cores = multiprocessing.cpu_count()
import matplotlib
import matplotlib.pyplot as plt
import plotnine
from plotnine import ggplot, aes, geom_line, geom_bar
from plotnine import *

def calculate_spearman(col1, col2):
    return stats.spearmanr(col1, col2).pvalue

def calculate_pvalues(df):
    df = df.dropna()._get_numeric_data()
    pvalues = Parallel(n_jobs=num_cores)(delayed(calculate_spearman)(df.iloc[:,i], df.iloc[:,j]) \
        for i in range(len(df.columns)) \
        for j in range(i,len(df.columns)))
    return pvalues

def calc_correlations(numericmatrix, strong_threshold = 0.5):
    correlation_spearman = numericmatrix.T.corr("spearman")
    p_values_spearman = calculate_pvalues(numericmatrix.T)
    pvalues = pd.DataFrame(np.nan, index=np.arange(len(numericmatrix.T.columns)), columns=numericmatrix.T.columns)
    i = 0
    while (i < len(numericmatrix.T.columns)):
    #for j in range(len(numericmatrix.T.columns)):
        pvalues.iloc[j,j:len(numericmatrix.T.columns)] = p_values_spearman[i:(i+len(numericmatrix.T.columns) - j)]
        i = i + (len(numericmatrix.T.columns) - j)
    pvalues["mag2"] = pvalues.columns
    correlation_spearman["mag2"] = correlation_spearman.columns,
    pvals_melted = pvalues.melt(id_vars = "mag2", var_name = "mag1", value_name = "pval")
    corr_melted = correlation_spearman.melt(id_vars = "mag2", var_name = "mag1", value_name = "corrcoef")
    pvals_melted = pvals_melted[(pvals_melted.pval == pvals_melted.pval)]
    combined_graph = pd.merge(corr_melted, pvals_melted, how = right)
    # Bonferroni: signifp / (0.5*(numgenes-1)*numgenes)
    # Sidak:
    filtered_graph = combined_graph[combined_graph.pval < (1-(1-0.05)**(1/len(combined_graph.index)))]
    filtered_strong_connections = combined_graph[abs(combined_graph.corrcoef) > strong_threshold]
    return combined_graph, filtered_strong_connections

datamatrix_tpm = pd.read_csv(os.path.join("..","input","MAG_tpm.csv"))
datamatrix_tpm = datamatrix_tpm[datamatrix_tpm.Genome != "unmapped"]
datamatrix_tpm.index = datamatrix_tpm.Genome
datamatrix_tpm = datamatrix_tpm.drop(columns=['Genome'])
qualify_cols = datamatrix_tpm.sum(axis=0)
numericmatrix_tpm = datamatrix_tpm.loc[:,qualify_cols > 0]
datamatrix_tpm["Method"] = "TPM"
combined_graph_tpm, filtered_strong_connections_tpm = calc_correlations(numericmatrix_tpm)
filtered_strong_connections_tpm.to_csv("filtered_strong_connections_tpm_fullpval.csv")
combined_graph_tpm.to_csv("combined_graph_tpm_fullpval.csv")
