# nichenet-learning-reproduction

This repository contains my learning-oriented reproduction of a NicheNet analysis workflow using **R**, **Seurat**, and **NicheNet**.

The goal of this repository is **not** to build a production-ready pipeline, but to:
- understand the main logic of the NicheNet workflow
- practice reading and organizing bioinformatics analysis code
- reproduce a representative ligand–receptor / ligand–target analysis example
- improve my ability to structure analysis scripts in a cleaner GitHub repository

---

## Project overview

In this project, I reorganized a longer exploratory script into several smaller scripts with clearer responsibilities:

- package installation and loading
- data downloading and preprocessing
- main NicheNet analysis
- result visualization

This repository is intended as a **learning and practice project**, showing my current progress in:
- understanding Seurat object handling
- running a basic NicheNet workflow
- checking intermediate results
- separating installation, preparation, analysis, and plotting steps

---

## Repository structure

```text
scripts/
├─ 00_install_packages.R
├─ 01_prepare_data.R
├─ 02_run_nichenet_analysis.R
└─ 03_plot_results.R
