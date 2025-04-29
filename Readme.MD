# Nanoparticle growth analysis

### In directory `functions`
- contains all the codes to perform pre-processing, optimization, post-processing.
- run `addpath("functions")` to use the function in this directory.

### In each directory(Au, CdSe, Iron_oxide, Pt(acac)2_monomeric, Pt(COD)Cl2, Pt_ex_situ)
- `*.csv` file contains raw data of nanoparticle size in unit of pixel.
- `Preprocessing*.m` performs data processing of the raw data contained in csv file and generates *.mat file containing pre-processed experimental data.
- `Fitting*.m` performs parameter optimization and generates data file (Ga*.mat) containing optimized parameter as well as hyperparameters.

### In root directory
- `summarize_results.m` post-processes optimization and experimental data, generating:
  - `fd.mat`: optimization data for all nanoparticle systems.
  - `re.mat`: post-processed data for all nanoparticle systems.
  - `ed.mat`: experimental data for all nanoparticle systems.
  - System-specific files (Pt_acac_2.mat, Pt_COD_Cl_2.mat, AuHCl_4.mat, Pt_ex.mat, Fe_xO_y.mat, CdSe.mat) containing experimental, post-processed, and optimization data for each system.

### In directory `Figure_codes`
- contains all the codes to plot figures
- run `addpath("Figure_codes")` to use the function in this directory.
- run `Figures_main.m` to generate all the figures in the Main manuscript.
- run `Figures_supp.m` to generate some of the figures in the Supporting Information.


# How to run
### In each directory(Au, CdSe, Iron_oxide, Pt(acac)2_monomeric, Pt(COD)Cl2, Pt_ex_situ)**
1. If there is `Preprocessing*.m` in the directory, run it to generate `*.mat` file. 
Skip this if there is no `Preprocessing*.m` in the directory.
2. Run `Fitting*.m` to perform optimization process. 
This will generate `*.mat` file containing optimized parameter along with hyperparameters. 
This takes very long time (>6hours with 64 core processor CPU).
### In root directory
3. Run `Summarize_Result.m`.
This performs post-processing of optimized parameters. 
Use the right file name generated in step 2 in the script `Summarize_Result.m`.
This generates the data file (`*.mat`).
### In directory `Figure_codes`
4. Run `Figures_main.m` to generate all the figures in the Main manuscript.
5. Run `Figures_supp.m` to generate some of the figures in the Supporting Information.
