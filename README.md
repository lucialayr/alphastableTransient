# Disturbances and long-transients give the illusion of alternative stable states

This repository contains code and data needed to reproduce the data analysis and figures of Chapter 3.3 of my thesis. data will be released once the paper is published.

<img align="left" src="figures/potential_landscape.png" style="width: 50%;">

## Folder structure

├── **code** &#x1F4C1;

│&nbsp; &nbsp; &nbsp; &nbsp;└── cluster &#x1F4C1;  *scripts to perform numerical simulations on High Performance Computing infrastructure*

│&nbsp; &nbsp; &nbsp; &nbsp;└──  `lpj_densities.R`&#x1F4C4; *Estimates and plots* $\hat{p}(\chi^C_{\text{BNE}})$ *and* $\hat{p}(\chi^C_{\text{BNE}}, k_{T_G})$

│&nbsp; &nbsp; &nbsp; &nbsp;└──  `potential_estimation.R`&#x1F4C4; *Estimates and plots* $\hat{U}(\chi^C_{\text{BNE}})$ *and* $\hat{U}(X)$

│&nbsp; &nbsp; &nbsp; &nbsp;└──  `trajectories.R`&#x1F4C4; *Plots trajectories of* $\chi^C_{\text{i}}(t)$ *and* $X(t)$

│&nbsp; &nbsp; &nbsp; &nbsp;└──  `utils.R`&#x1F4C4; *Various helper functions for plotting*

├── **data** &#x1F4C1;  *model output on various steps of processing*

│&nbsp; &nbsp; &nbsp; &nbsp;└──  `final_states_a<alpha>_k<k>.csv`&#x1F4C4; *Data for potential estimations*

│&nbsp; &nbsp; &nbsp; &nbsp;└──  `trajectoriess_a<alpha>_k<k>.csv`&#x1F4C4; *Data for trajectories*

│&nbsp; &nbsp; &nbsp; &nbsp;└── processed &#x1F4C1;  *intermediate data*

│&nbsp; &nbsp; &nbsp; &nbsp;└── ext &#x1F4C1;  *External helper files*
  
├── **figures** &#x1F4C1; *Contains all the plots of the paper.

├── `*patches2.duckdb` &#x1F986; *DuckDB database which contains patch LPJ-GUESS data. See https://github.com/lucialayr/borealRecovery for details on methodlogy* 


*The figures is published under a CC-BY-SA license. To reuse please cite:*

Layritz, L. S. (2024). *Illustrations for 'Disturbances in the evergreen boreal forest and their impact on 21st century vegetation and climate dynamics - A stochastic modeling approach' (Doctoral thesis)*. Zenodo. https://doi.org/10.5281/zenodo.13731735
 