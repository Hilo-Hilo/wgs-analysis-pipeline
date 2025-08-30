# Whole Genome Sequencing Analysis Pipeline

A comprehensive pipeline for analyzing whole genome sequencing (WGS) data from raw FASTQ files to annotated variants. This repository contains scripts, documentation, and lessons learned from processing 30x coverage human WGS data.

## 🔬 Overview

This pipeline processes paired-end Illumina WGS data through quality control, alignment, variant calling, and annotation stages. It supports both CHM13 T2T and GRCh38 reference genomes and includes optimizations for large-scale data processing.

### Key Features
- Complete end-to-end WGS analysis workflow
- Support for CHM13 T2T v2.0 and GRCh38 references
- Optimized for 30x coverage whole genome data
- Parallel processing capabilities
- Comprehensive variant annotation with VEP
- Personal genomics analysis tools
- Cloud and local deployment options

## 📊 Pipeline Stages

1. **Quality Control**
   - FastQC for raw read assessment
   - Fastp for adapter trimming and quality filtering
   - MultiQC for aggregated reports

2. **Alignment**
   - BWA-MEM for GRCh38 (standard clinical reference)
   - BWA for CHM13 T2T (complete telomere-to-telomere assembly)
   - SAMtools for BAM processing

3. **Variant Calling**
   - BCFtools for traditional variant calling
   - DeepVariant for deep learning-based calling
   - Joint calling and GVCF generation

4. **Annotation**
   - VEP (Variant Effect Predictor) v110
   - gnomAD v4.0 population frequencies
   - ClinVar clinical significance
   - COSMIC cancer variants

5. **Analysis**
   - High-impact variant identification
   - Pharmacogenomics analysis
   - Ancestry and trait prediction
   - Clinical interpretation templates

## 🚀 Quick Start

### Prerequisites

```bash
# Create conda environment
conda create -n wgs_analysis python=3.9
conda activate wgs_analysis

# Install core tools
conda install -c bioconda bwa samtools bcftools
conda install -c bioconda fastp fastqc multiqc
```

### Basic Usage

```bash
# 1. Quality control
./scripts/quality_control.sh

# 2. Alignment to reference
./scripts/chm13_mapping.sh  # For CHM13
# OR
./scripts/grch38_mapping.sh  # For GRCh38

# 3. Variant calling
./scripts/variant_calling.sh

# 4. Annotation
./scripts/vep_annotation.sh
```

## 💡 Important Discoveries

### BWA-MEM2 Memory Requirements
**Critical Finding**: BWA-MEM2 requires **128GB RAM** for GRCh38 indexing, not the 28GB stated in documentation.

```bash
# Will fail with less than 128GB RAM
bwa-mem2 index GRCh38_latest_genomic.fna

# Solution: Use high-memory cloud instance or original BWA
```

See [BWA_MEM2_ROOT_CAUSE_ANALYSIS.md](documentation/BWA_MEM2_ROOT_CAUSE_ANALYSIS.md) for details.

### Performance Metrics
- **Input**: 131GB paired-end FASTQ files
- **Processing Time**: ~24 hours on 8 cores
- **Output**: ~5 million variants
- **Mapping Rate**: ~88.6% (typical for WGS)
- **Storage Required**: ~400GB total

## 📁 Repository Structure

```
wgs-analysis-pipeline/
├── scripts/              # Pipeline scripts
│   ├── alignment/       # Mapping scripts
│   ├── variant_calling/ # Variant calling scripts
│   ├── annotation/      # VEP annotation scripts
│   └── quality_control/ # QC scripts
├── analysis/            # Analysis tools
│   ├── pharmacogenomics.py
│   ├── trait_analysis.py
│   └── variant_impact.py
├── documentation/       # Detailed guides
│   ├── environment_setup.md
│   ├── cloud_infrastructure_guide.md
│   └── complete_analysis_workflow.md
├── templates/          # Example configurations
└── examples/          # Sample data and outputs
```

## 🧬 Analysis Tools

### Pharmacogenomics Analysis
Analyze genetic variants affecting drug metabolism:

```python
python analysis/pharmacogenomics.py --vcf input.vcf.gz
```

### Trait Prediction
Predict genetic traits and ancestry:

```python
python analysis/trait_analysis.py --vcf input.vcf.gz
```

### High-Impact Variant Analysis
Identify clinically significant variants:

```python
python analysis/variant_impact.py --vcf input.vcf.gz --impact HIGH
```

## ☁️ Cloud Deployment

### Google Cloud Platform
Optimized configuration for WGS analysis:

```bash
# Create high-memory instance for alignment
gcloud compute instances create wgs-analysis \
    --machine-type=n2d-highmem-16 \
    --boot-disk-size=500GB \
    --zone=us-central1-a
```

Cost optimization strategies:
- Use preemptible instances for non-critical stages
- Downgrade to n2-standard-4 after alignment
- Store results in Cloud Storage

See [cloud_infrastructure_guide.md](documentation/cloud_infrastructure_guide.md) for detailed setup.

## 📚 Documentation

- [Environment Setup](documentation/environment_setup.md) - Tool installation and configuration
- [Complete Analysis Workflow](documentation/complete_analysis_workflow.md) - Step-by-step pipeline guide
- [Cloud Infrastructure Guide](documentation/cloud_infrastructure_guide.md) - GCP deployment
- [BWA-MEM2 Analysis](documentation/BWA_MEM2_ROOT_CAUSE_ANALYSIS.md) - Memory requirement investigation

## 🔧 Troubleshooting

### Common Issues

1. **Out of Memory During Alignment**
   - Solution: Use original BWA instead of BWA-MEM2
   - Alternative: Use cloud instance with 128GB+ RAM

2. **Low Mapping Rate (<80%)**
   - Check read quality with FastQC
   - Verify correct reference genome version
   - Consider additional quality filtering

3. **VEP Installation Issues**
   - Use Docker container: `ensemblorg/ensembl-vep`
   - Or use web API for small variant sets

## 📊 Expected Results

From 30x WGS data, expect:
- **Total variants**: 4-5 million
- **SNVs**: ~4 million
- **Indels**: ~0.8 million
- **Novel variants**: 5-10%
- **High-impact variants**: 300-500

## 🤝 Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

### Areas for Contribution
- Additional variant callers (GATK, Strelka2)
- Structural variant detection
- Copy number variation analysis
- Improved visualization tools
- Additional population databases

## 📖 Citations

If you use this pipeline, please cite:

```bibtex
@software{wgs_analysis_pipeline,
  title = {Whole Genome Sequencing Analysis Pipeline},
  year = {2025},
  url = {https://github.com/yourusername/wgs-analysis-pipeline}
}
```

### Key Tools Used
- BWA: [Li and Durbin, 2009](https://doi.org/10.1093/bioinformatics/btp324)
- SAMtools: [Danecek et al., 2021](https://doi.org/10.1093/gigascience/giab008)
- BCFtools: [Danecek et al., 2021](https://doi.org/10.1093/gigascience/giab008)
- DeepVariant: [Poplin et al., 2018](https://doi.org/10.1038/nbt.4235)
- VEP: [McLaren et al., 2016](https://doi.org/10.1186/s13059-016-0974-4)

## 📝 License

This project is licensed under the MIT License - see [LICENSE](LICENSE) file for details.

## ⚠️ Disclaimer

This pipeline is for research and educational purposes. For clinical applications, please use validated clinical pipelines and consult with qualified geneticists and healthcare providers.

## 🙏 Acknowledgments

- CHM13 T2T Consortium for the complete human reference
- gnomAD team for population frequency data
- Ensembl VEP team for annotation tools
- The open-source bioinformatics community

## 📧 Contact

For questions or issues, please open a GitHub issue or see the documentation.

---

**Note**: This repository contains only pipeline code and documentation. No personal genomic data or identifying information is included.