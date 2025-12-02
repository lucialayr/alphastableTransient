# Reproducible Docker container for alpha-stable transient analysis
# Based on rocker/r-ver which includes R 4.4.2 (lightweight, no RStudio)
# Note: For Apple Silicon (ARM64) Macs, build with: docker build --platform linux/amd64 -t alphastable-transient:latest .

FROM rocker/binder 

LABEL maintainer="Lucia Layritz"
LABEL description="Reproducible environment for 'Disturbances and long-transients give the illusion of alternative stable states'"
LABEL org.opencontainers.image.title="alphastable-transient"

# Copy environment specification
COPY environment.yml /tmp/environment.yml

# Create conda environment and clean up to save space
RUN mamba env update -f /tmp/environment.yml && \
    mamba clean -afy

# Install R packages using remotes
# Only installing specific tidyverse packages that are actually used instead of the full tidyverse
RUN install2.r --error --skipinstalled remotes tidyverse scico cowplot ggnewscale

