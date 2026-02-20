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
# Use mamba + flexible priority for more reliable solves in CI.
RUN conda config --set always_yes true && \
    conda config --set channel_priority flexible && \
    conda install -n base -c conda-forge mamba && \
    mamba create -n wgs_analysis -c conda-forge -c bioconda \
    python=3.11 \
    fastqc \
    fastp \
    bwa \
    samtools \
    bcftools \
    ensembl-vep && \
    conda clean -ya

# Ensure tools from the analysis environment are on PATH
ENV PATH="/opt/miniconda/envs/wgs_analysis/bin:$PATH"

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

# Build information
ARG BUILD_DATE
ARG BUILD_VERSION
ARG VCS_REF
LABEL org.label-schema.build-date=$BUILD_DATE
LABEL org.label-schema.version=$BUILD_VERSION
LABEL org.label-schema.vcs-ref=$VCS_REF
LABEL org.label-schema.schema-version="1.0"