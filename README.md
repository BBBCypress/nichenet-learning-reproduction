# nichenet-learning-reproduction

This repository is a learning-based reproduction of a NicheNet workflow using **R**, **Seurat**, and **nichenetr**.

The main purpose of this project is to:
- follow the logic of a typical NicheNet analysis
- practice reading and reorganizing bioinformatics code
- split a long exploratory script into clearer analysis steps
- make the workflow easier to read and rerun

This is still a **learning / practice repository**, not a production-level pipeline.

---

## What is included

The current workflow includes these parts:

1. package installation and loading  
2. data download and preprocessing  
3. receiver / sender definition  
4. differential expression analysis in receiver cells  
5. ligand activity ranking  
6. ligand-target and ligand-receptor prediction  
7. heatmap-based visualization of selected results  

---

## Repository structure

```text
scripts/
├─ 00_install_packages.R
├─ 01_prepare_data.R
├─ 02_run_nichenet_analysis.R
└─ 03_plot_results.R
