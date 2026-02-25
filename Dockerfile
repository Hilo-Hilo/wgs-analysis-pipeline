# WGS Analysis Pipeline - 16GB RAM Optimized Container
# This container includes all dependencies for the WGS analysis pipeline

FROM ubuntu:22.04

# Metadata
LABEL maintainer="WGS Analysis Pipeline"
LABEL description="16GB RAM optimized whole genome sequencing analysis pipeline"
LABEL version="2.0-16GB"

# Prevent interactive prompts during installation
ARG DEBIAN_FRONTEND=noninteractive
ARG MINICONDA_INSTALLER=Miniconda3-py39_24.3.0-0-Linux-x86_64.sh

# Set working directory
WORKDIR /opt/wgs-pipeline

# Install system dependencies
RUN apt-get update && apt-get install -y \
    # Basic utilities
    wget \
    curl \
    git \
    unzip \
    python3 \
    python3-pip \
    bc \
    time \
    htop \
    # Build essentials
    build-essential \
    cmake \
    autoconf \
    automake \
    # Libraries
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    # Perl and Java for bioinformatics tools
    perl \
    default-jre \
    # Core HTS tools from apt (stable, modern samtools on arm64)
    samtools \
    bcftools \
    # Clean up
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install miniconda for bioinformatics tools
RUN wget "https://repo.anaconda.com/miniconda/${MINICONDA_INSTALLER}" -O miniconda.sh \
    && bash miniconda.sh -b -p /opt/miniconda \
    && rm miniconda.sh \
    && /opt/miniconda/bin/conda clean -ya

# Add conda to PATH
ENV PATH="/opt/miniconda/bin:$PATH"

# Create conda environment with bioinformatics tools
# Use mamba + retry/backoff for transient channel/DNS failures in CI.
RUN bash -lc 'set -euo pipefail; \
    retry_with_backoff() { \
      local max_attempts="$1"; \
      local delay_seconds="$2"; \
      shift 2; \
      local attempt=1; \
      local exit_code=0; \
      while true; do \
        if "$@"; then \
          return 0; \
        else \
          exit_code="$?"; \
        fi; \
        if [ "$attempt" -ge "$max_attempts" ]; then \
          echo "ERROR: command failed after ${max_attempts} attempts: $*" >&2; \
          return "$exit_code"; \
        fi; \
        echo "WARN: command failed (attempt ${attempt}/${max_attempts}, exit ${exit_code}): $*" >&2; \
        echo "Retrying in ${delay_seconds}s..." >&2; \
        sleep "$delay_seconds"; \
        attempt=$((attempt + 1)); \
        delay_seconds=$((delay_seconds * 2)); \
      done; \
    }; \
    conda config --set always_yes true; \
    conda config --set channel_priority flexible; \
    retry_with_backoff 4 5 conda install -n base --override-channels -c conda-forge mamba; \
    retry_with_backoff 4 5 mamba create -n wgs_analysis --override-channels -c conda-forge -c bioconda \
      python=3.11 \
      fastqc \
      fastp \
      bwa \
      "samtools>=1.18" \
      "bcftools>=1.18"; \
    conda clean -ya'

# Ensure tools from the analysis environment are on PATH
ENV PATH="/opt/miniconda/envs/wgs_analysis/bin:$PATH"

# Guardrail: fail image build if legacy samtools slipped in
RUN bash -lc 'set -e; v=$(samtools 2>&1 | awk "/^Version:/{print \$2; exit}"); \
  if [ -z "$v" ]; then v=$(samtools --version 2>/dev/null | awk "NR==1{print \$2}"); fi; \
  echo "Detected samtools: $v"; \
  case "$v" in 0.*) echo "ERROR: legacy samtools detected ($v)"; exit 1;; esac'

# Install Python packages for testing and utilities
RUN pip3 install --no-cache-dir \
    pytest \
    pandas \
    matplotlib \
    seaborn \
    pyyaml

# Copy pipeline files
COPY . /opt/wgs-pipeline/

# Set proper permissions
RUN chmod +x scripts/*.sh tests/*.sh tests/*.py

# Create necessary directories
RUN mkdir -p data/{raw,processed,reference} results logs temp

# Create a non-root user for security
RUN useradd -m -s /bin/bash wgsuser && \
    chown -R wgsuser:wgsuser /opt/wgs-pipeline
USER wgsuser

# Set environment variables
ENV WGS_PIPELINE_HOME="/opt/wgs-pipeline"
ENV PATH="/opt/wgs-pipeline/scripts:$PATH"

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD check_requirements.sh --min-ram 1 --min-disk 10 --skip-conda --skip-env --skip-tools || exit 1

# Run all commands inside the wgs_analysis conda environment
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "wgs_analysis"]

# Default command
CMD ["/bin/bash"]

# Expose no ports (this is a compute container, not a service)

# Note: VEP is intentionally not installed in the base image to avoid
# incompatible legacy samtools resolution on ARM platforms.

# Build information
ARG BUILD_DATE
ARG BUILD_VERSION
ARG VCS_REF
LABEL org.label-schema.build-date=$BUILD_DATE
LABEL org.label-schema.version=$BUILD_VERSION
LABEL org.label-schema.vcs-ref=$VCS_REF
LABEL org.label-schema.schema-version="1.0"