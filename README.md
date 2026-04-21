# drought-propagation-code

MATLAB code for drought propagation analysis based on gridded drought-index time series over South Asia.

## Overview

This repository contains MATLAB scripts used for drought propagation analysis, non-stationary cluster-based analysis, and wavelet-based propagation-time analysis.

The input gridded drought-index time series are **not included in this repository** because of their large size. The datasets used in this study were preprocessed from the monthly Global Land Data Assimilation System (GLDAS) products **GLDAS_NOAH025_M_2.0** and **GLDAS_NOAH025_M_2.1**. Detailed preprocessing procedures are described in the associated manuscript.

## Repository structure

drought-propagation-code/
├─ scripts/
│  ├─ RUN1_DP_S.m
│  ├─ RUN2_PCA_CLUSTER.m
│  ├─ RUN3_DP_NS.m
│  └─ RUN4_WAVELET.m
├─ functions/
│  └─ select_best_copula_builtin.m
├─ data_sample/
└─ README.md
