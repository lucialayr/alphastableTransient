# Reproducible Docker container for alpha-stable transient analysis
# Based on rocker/r-ver which includes R 4.4.2 (lightweight, no RStudio)
# Note: For Apple Silicon (ARM64) Macs, build with: docker build --platform linux/amd64 -t alphastable-transient:latest .

FROM --platform=linux/amd64 rocker/r-ver:4.4.2

LABEL maintainer="Lucia Layritz"
LABEL description="Reproducible environment for 'Disturbances and long-transients give the illusion of alternative stable states'"
LABEL org.opencontainers.image.title="alphastable-transient"

# Install system dependencies for R packages and conda
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    ca-certificates \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    && rm -rf /var/lib/apt/lists/*

# Install Miniforge3 (includes mamba) for faster conda package resolution
ENV CONDA_DIR=/opt/conda
ENV PATH=${CONDA_DIR}/bin:${PATH}

# Auto-detect platform (aarch64 for ARM, x86_64 for Intel)
RUN ARCH=$(uname -m) && \
    wget --quiet https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-${ARCH}.sh -O /tmp/miniforge.sh && \
    /bin/bash /tmp/miniforge.sh -b -p ${CONDA_DIR} && \
    rm /tmp/miniforge.sh && \
    ${CONDA_DIR}/bin/conda clean -afy

# Copy environment specification
COPY environment.yml /tmp/environment.yml

# Create conda environment and clean up to save space
RUN mamba env create -f /tmp/environment.yml && \
    mamba clean -afy && \
    rm /tmp/environment.yml

# Activate conda environment for all shell sessions
ENV PATH=${CONDA_DIR}/envs/alphastableTransient/bin:${PATH}
ENV CONDA_DEFAULT_ENV=alphastableTransient

# Install R packages using remotes
# Only installing specific tidyverse packages that are actually used instead of the full tidyverse
RUN install2.r --error --skipinstalled remotes && \
    R -e "remotes::install_version('ggplot2', version = '3.5.1', repos = 'https://cloud.r-project.org')" && \
    R -e "remotes::install_version('dplyr', version = '1.1.4', repos = 'https://cloud.r-project.org')" && \
    R -e "remotes::install_version('tidyr', version = '1.3.1', repos = 'https://cloud.r-project.org')" && \
    R -e "remotes::install_version('readr', version = '2.1.5', repos = 'https://cloud.r-project.org')" && \
    R -e "remotes::install_version('purrr', version = '1.0.2', repos = 'https://cloud.r-project.org')" && \
    R -e "remotes::install_version('duckdb', version = '1.1.3', repos = 'https://cloud.r-project.org')" && \
    R -e "remotes::install_version('scico', version = '1.5.0', repos = 'https://cloud.r-project.org')" && \
    R -e "remotes::install_version('cowplot', version = '1.2.0', repos = 'https://cloud.r-project.org')" && \
    R -e "remotes::install_version('ggnewscale', version = '0.5.1', repos = 'https://cloud.r-project.org')" && \
    rm -rf /tmp/downloaded_packages /tmp/Rtmp*

# Set working directory
WORKDIR /workspace

# Copy project files
COPY code/ /workspace/code/
COPY data/processed/ /workspace/data/processed/
COPY README.md /workspace/

# Create data directory for generated outputs
RUN mkdir -p /workspace/data/generated /workspace/figures

# Set default command to R 
CMD ["R"]
