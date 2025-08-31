# WGS Analysis Pipeline - 16GB RAM Optimized Container
# This container includes all dependencies for the WGS analysis pipeline

FROM ubuntu:22.04

# Metadata
LABEL maintainer="WGS Analysis Pipeline"
LABEL description="16GB RAM optimized whole genome sequencing analysis pipeline"
LABEL version="2.0-16GB"

# Prevent interactive prompts during installation
ARG DEBIAN_FRONTEND=noninteractive

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
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh \
    && bash miniconda.sh -b -p /opt/miniconda \
    && rm miniconda.sh \
    && /opt/miniconda/bin/conda clean -ya

# Add conda to PATH
ENV PATH="/opt/miniconda/bin:$PATH"

# Create conda environment with bioinformatics tools
RUN conda create -n wgs_analysis python=3.9 -y && \
    conda install -n wgs_analysis -c bioconda -c conda-forge \
    fastqc=0.12.1 \
    fastp=0.23.4 \
    bwa=0.7.17 \
    samtools=1.17 \
    bcftools=1.17 \
    ensembl-vep=110.1 \
    -y && \
    conda clean -ya

# Activate environment by default
RUN echo "source activate wgs_analysis" > ~/.bashrc
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
    CMD check_requirements.sh --min-ram 1 --min-disk 10 || exit 1

# Default command
CMD ["/bin/bash"]

# Expose no ports (this is a compute container, not a service)

# Documentation
RUN echo "WGS Analysis Pipeline Container Ready!" > /tmp/container_ready.txt && \
    echo "Available commands:" >> /tmp/container_ready.txt && \
    echo "  check_requirements.sh  - Validate system" >> /tmp/container_ready.txt && \
    echo "  quality_control.sh     - Run QC analysis" >> /tmp/container_ready.txt && \
    echo "  data_cleaning.sh       - Clean reads" >> /tmp/container_ready.txt && \
    echo "  alignment.sh           - Align to genome" >> /tmp/container_ready.txt && \
    echo "  variant_calling.sh     - Call variants" >> /tmp/container_ready.txt && \
    echo "  vep_annotation.sh      - Annotate variants" >> /tmp/container_ready.txt && \
    echo "" >> /tmp/container_ready.txt && \
    echo "Example usage:" >> /tmp/container_ready.txt && \
    echo "  docker run -v /path/to/data:/opt/wgs-pipeline/data wgs-pipeline" >> /tmp/container_ready.txt

# Build information
ARG BUILD_DATE
ARG BUILD_VERSION
ARG VCS_REF
LABEL org.label-schema.build-date=$BUILD_DATE
LABEL org.label-schema.version=$BUILD_VERSION
LABEL org.label-schema.vcs-ref=$VCS_REF
LABEL org.label-schema.schema-version="1.0"