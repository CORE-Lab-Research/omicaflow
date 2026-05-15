# Use micromamba for significantly faster environment builds
FROM mambaorg/micromamba:1.5.1

# Metadata
LABEL maintainer="OmicaFlow Team"
LABEL version="1.0.0"
LABEL description="Optimized Docker image for OmicaFlow Multi-Omics Pipeline"

USER root
# Install basic system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    procps \
    && rm -rf /var/lib/apt/lists/*

USER $MAMBA_USER
WORKDIR /app

# Copy environment files to leverage Docker caching
COPY --chown=$MAMBA_USER:$MAMBA_USER envs/r_base.yml /tmp/r_base.yml
COPY --chown=$MAMBA_USER:$MAMBA_USER envs/snakemake.yml /tmp/snakemake.yml

# Create a single unified environment for the container
# This avoids the overhead of managing multiple environments inside Docker
RUN micromamba create -y -n omicaflow -f /tmp/r_base.yml && \
    micromamba install -y -n omicaflow -f /tmp/snakemake.yml && \
    micromamba clean --all --yes

# Set the environment as active
ARG MAMBA_DOCKERFILE_ACTIVATE=1
ENV PATH="/opt/conda/envs/omicaflow/bin:$PATH"

# Copy the rest of the project
COPY --chown=$MAMBA_USER:$MAMBA_USER . .

# Ensure logs and results directories exist and are writable
RUN mkdir -p logs results data/raw && \
    chmod -R 777 logs results data/raw

# Default command: run dry-run to verify setup
CMD ["snakemake", "--cores", "1", "-n"]
