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

## Reproducing the analysis

You have two options to reproduce the analysis:

### Option 1: Using Conda Environment (Recommended for Native Performance)

This approach provides better performance, especially on Apple Silicon Macs where Docker runs in emulation mode (30-50% slower).

#### Step 1: Create Python environment

```bash
# Navigate to the repository
cd /path/to/alphastableTransient

# Create conda environment from specification
conda env create -f environment.yml

# Activate the environment
conda activate alphastableTransient
```

#### Step 2: Install R packages

Open R and run:

```r
# Install remotes if not already installed
if (!require("remotes")) install.packages("remotes")

# Install required packages with specific versions
remotes::install_version('ggplot2', version = '3.5.1', repos = 'https://cloud.r-project.org')
remotes::install_version('dplyr', version = '1.1.4', repos = 'https://cloud.r-project.org')
remotes::install_version('tidyr', version = '1.3.1', repos = 'https://cloud.r-project.org')
remotes::install_version('readr', version = '2.1.5', repos = 'https://cloud.r-project.org')
remotes::install_version('purrr', version = '1.0.2', repos = 'https://cloud.r-project.org')
remotes::install_version('duckdb', version = '1.1.3', repos = 'https://cloud.r-project.org')
remotes::install_version('scico', version = '1.5.0', repos = 'https://cloud.r-project.org')
remotes::install_version('cowplot', version = '1.2.0', repos = 'https://cloud.r-project.org')
remotes::install_version('ggnewscale', version = '0.5.1', repos = 'https://cloud.r-project.org')
```

#### Step 3: Run the analysis

```bash
# Make sure conda environment is activated
conda activate alphastableTransient

# Run simulations (adjust number of simulations per batch as needed)
# For full analysis: 10000 simulations per batch × 10 batches = 100,000 total per parameter combination
bash code/numerics/start_all_runs_parallel_splitting.sh 10000

# For quick test: 100 simulations per batch
bash code/numerics/start_all_runs_parallel_splitting.sh 100

# Generate figures with R
Rscript code/Figure2_potential_estimations.R
Rscript code/Figure3_trajectories.R
Rscript code/Figure5_bifurcation_diagram.R
```

---

### Option 2: Using Docker Container

The repository is wrapped in a Docker container to ensure all libraries are installed in the correct version and provide reproducibility. 

**Note for Apple Silicon Macs**: Docker runs in emulation mode, which is 30-50% slower than native execution. Consider using Option 1 (conda environment) for better performance. If using Docker, configure Docker Desktop to use at least 10 CPU cores and 12GB RAM (Settings → Resources).

#### Building the Docker image

```bash
# For Intel/AMD systems (Linux, Windows):
docker build -t alphastable-transient:latest .

# For Apple Silicon Macs (M1/M2/M3):
docker build --platform linux/amd64 -t alphastable-transient:latest .
```

This step only needs to be done once and will take 10-15 minutes to install all dependencies.

#### Running analyses in the container

Once built, you can run scripts from the command line. The container will mount your local `data/` and `figures/` directories so results are saved to your host machine:

**Generate numerical simulations:**
```bash
# Full analysis (10000 simulations per batch)
docker run --rm -v $(pwd)/data:/workspace/data alphastable-transient:latest bash -c "bash code/numerics/start_all_runs_parallel_splitting.sh 10000"

# Quick test (100 simulations per batch)
docker run --rm -v $(pwd)/data:/workspace/data alphastable-transient:latest bash -c "bash code/numerics/start_all_runs_parallel_splitting.sh 100"
```

**Create figures:**
```bash
docker run --rm -v $(pwd)/data:/workspace/data -v $(pwd)/figures:/workspace/figures alphastable-transient:latest Rscript code/potential_estimations.R

docker run --rm -v $(pwd)/data:/workspace/data -v $(pwd)/figures:/workspace/figures alphastable-transient:latest Rscript code/lpj_densities.R

docker run --rm -v $(pwd)/data:/workspace/data -v $(pwd)/figures:/workspace/figures alphastable-transient:latest Rscript code/trajctories.R
```

**Interactive session:**
```bash
# R session
docker run --rm -it -v $(pwd)/data:/workspace/data -v $(pwd)/figures:/workspace/figures alphastable-transient:latest R

# Bash shell
docker run --rm -it -v $(pwd)/data:/workspace/data -v $(pwd)/figures:/workspace/figures alphastable-transient:latest /bin/bash
```

### Environment specifications

- **Base**: `rocker/r-ver:4.4.2` (lightweight R installation)
- **R**: 4.4.2
- **Python**: 3.13
- **R packages**: Specific versions installed via `remotes` (see `Dockerfile` for details)
  - ggplot2, dplyr, tidyr, readr, purrr (core tidyverse packages)
  - duckdb, scico, cowplot, ggnewscale
- **Python packages**: Specified in `environment.yml`

---

*The conceptual figure is published under a CC-BY-SA license. To reuse please cite:*

Layritz, L. S. (2024). *Illustrations for 'Disturbances in the evergreen boreal forest and their impact on 21st century vegetation and climate dynamics - A stochastic modeling approach' (Doctoral thesis)*. Zenodo. https://doi.org/10.5281/zenodo.13731735
 