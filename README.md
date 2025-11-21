# Disturbances and long-transients give the illusion of alternative stable states

This repository contains code and data needed to reproduce the data analysis and figures of Chapter 3.3 of my thesis. Data will be released once the paper is published.

<img src="figures/potential_landscape.png" width="600" />

## Folder structure

├── **code** &#x1F4C1;

│&nbsp; &nbsp; &nbsp; &nbsp;└── numerics &#x1F4C1;  *scripts to perform numerical simulations*

│&nbsp; &nbsp; &nbsp; &nbsp;└── `Introduction_splitting_scheme.ipynb` &#x1F4C1; *An interactive notebook to give a friendly introduction to the splitting scheme*

│&nbsp; &nbsp; &nbsp; &nbsp;└──  `lpj_densities.R`&#x1F4C4; *Estimates and plots* $\hat{p}(\chi^C_{\text{BNE}})$ *and* $\hat{p}(\chi^C_{\text{BNE}}, k_{T_G})$

│&nbsp; &nbsp; &nbsp; &nbsp;└──  `potential_estimation.R`&#x1F4C4; *Estimates and plots* $\hat{U}(\chi^C_{\text{BNE}})$ *and* $\hat{U}(X)$

│&nbsp; &nbsp; &nbsp; &nbsp;└──  `trajectories.R`&#x1F4C4; *Plots trajectories of* $\chi^C_{\text{i}}(t)$ *and* $X(t)$

├── **data** &#x1F4C1;  *The .csv files will not be in the repository but need to be created with the scripts in `data/numerics`*

│&nbsp; &nbsp; &nbsp; &nbsp;└── processed &#x1F4C1;  *LPJ-GUESS simulation data and helper files*

│&nbsp; &nbsp; &nbsp; &nbsp;└──  `final_states_a<alpha>_k<k>.csv`&#x1F4C4; *Data for potential estimations*

│&nbsp; &nbsp; &nbsp; &nbsp;└──  `trajectoriess_a<alpha>_k<k>.csv`&#x1F4C4; *Data for trajectories*
  
├── **figures** &#x1F4C1; *Contains all the plots of the paper.*

$^*$*See https://github.com/lucialayr/borealRecovery for details on methodology* 

</p>


*The conceptual figure is published under a CC-BY-SA license. To reuse please cite:*

Layritz, L. S. (2024). *Illustrations for 'Disturbances in the evergreen boreal forest and their impact on 21st century vegetation and climate dynamics - A stochastic modeling approach' (Doctoral thesis)*. Zenodo. https://doi.org/10.5281/zenodo.13731735
 