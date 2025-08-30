#!/usr/bin/env python3
"""
Pharmacogenomics Analysis Script
Extracts drug response variants from VEP-annotated VCF
"""

import subprocess
import re
from collections import defaultdict

# Key pharmacogenomic genes and their important variants
PHARMACO_GENES = {
    'CYP2D6': {
        'chromosome': 'chr22',
        'region': '22:42522500-42555000',
        'drugs': ['Codeine', 'Tramadol', 'Antidepressants', 'Antipsychotics'],
        'star_alleles': {
            '*4': ['rs3892097', 'rs1065852'],  # Most common poor metabolizer
            '*3': ['rs35742686'],
            '*6': ['rs5030655'],
            '*10': ['rs1065852', 'rs1058164']
        }
    },
    'CYP2C19': {
        'chromosome': 'chr10', 
        'region': '10:94760681-94855547',
        'drugs': ['Clopidogrel', 'Omeprazole', 'Antidepressants'],
        'star_alleles': {
            '*2': ['rs4244285'],  # Poor metabolizer
            '*3': ['rs4986893'],  # Poor metabolizer  
            '*17': ['rs12248560']  # Ultrarapid metabolizer
        }
    },
    'CYP2C9': {
        'chromosome': 'chr10',
        'region': '10:94938683-94988859', 
        'drugs': ['Warfarin', 'Phenytoin'],
        'star_alleles': {
            '*2': ['rs1799853'],  # Intermediate metabolizer
            '*3': ['rs1057910']   # Poor metabolizer
        }
    },
    'SLCO1B1': {
        'chromosome': 'chr12',
        'region': '12:21178615-21239931',
        'drugs': ['Statins'],
        'star_alleles': {
            '*5': ['rs4149056']  # Increased statin myopathy risk
        }
    }
}

def extract_pharmaco_variants(vcf_path, sample_name='SAMPLE001'):
    """Extract pharmacogenomic variants from VCF"""
    results = {}
    
    for gene, info in PHARMACO_GENES.items():
        print(f"Analyzing {gene}...")
        
        # Extract variants in gene region
        cmd = [
            'bcftools', 'view', '-r', info['region'], vcf_path
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"Error extracting {gene}: {result.stderr}")
            continue
            
        # Parse variants
        gene_variants = []
        for line in result.stdout.split('\n'):
            if line.startswith('#') or not line.strip():
                continue
                
            fields = line.split('\t')
            if len(fields) < 10:
                continue
                
            chrom, pos, rsid, ref, alt = fields[0:5]
            genotype = fields[9].split(':')[0]  # GT field
            csq = fields[7] if 'CSQ=' in fields[7] else ''
            
            # Parse CSQ field for gene symbol
            if f'SYMBOL={gene}' in csq or gene in csq:
                variant_info = {
                    'position': f"{chrom}:{pos}",
                    'ref': ref,
                    'alt': alt,
                    'genotype': genotype,
                    'rsid': rsid if rsid != '.' else None,
                    'csq': csq
                }
                gene_variants.append(variant_info)
        
        results[gene] = {
            'variants': gene_variants,
            'drugs_affected': info['drugs'],
            'star_alleles': info['star_alleles']
        }
    
    return results

def interpret_pharmaco_results(results):
    """Interpret pharmacogenomic findings"""
    report = []
    
    for gene, data in results.items():
        report.append(f"\n=== {gene} Analysis ===")
        report.append(f"Drugs affected: {', '.join(data['drugs_affected'])}")
        report.append(f"Variants found: {len(data['variants'])}")
        
        if data['variants']:
            report.append("Key variants:")
            for variant in data['variants'][:5]:  # Show first 5
                gt = variant['genotype']
                pos = variant['position']
                ref_alt = f"{variant['ref']}>{variant['alt']}"
                report.append(f"  {pos} {ref_alt} ({gt})")
        
        # Basic interpretation
        variant_count = len([v for v in data['variants'] if v['genotype'] in ['0/1', '1/0', '1/1']])
        
        if gene == 'CYP2D6':
            if variant_count >= 2:
                report.append("⚠️  Multiple variants - may affect codeine/tramadol metabolism")
            else:
                report.append("✅ Likely normal codeine/tramadol metabolism")
                
        elif gene == 'CYP2C19':
            if variant_count >= 1:
                report.append("⚠️  May have reduced clopidogrel efficacy")
            else:
                report.append("✅ Likely normal clopidogrel response")
                
        elif gene == 'CYP2C9':
            if variant_count >= 1:
                report.append("⚠️  May require warfarin dose adjustment")
            else:
                report.append("✅ Likely normal warfarin sensitivity")
                
        elif gene == 'SLCO1B1':
            if variant_count >= 1:
                report.append("⚠️  Increased risk of statin-induced muscle problems")
            else:
                report.append("✅ Normal statin tolerance expected")
    
    return '\n'.join(report)

if __name__ == "__main__":
    import sys
    import os
    
    # Set up conda environment
    os.system("source ~/anaconda3/etc/profile.d/conda.sh && conda activate wgs_analysis")
    
    vcf_path = "../results/variants_fully_annotated.vcf.gz"
    
    if not os.path.exists(vcf_path):
        print(f"VCF file not found: {vcf_path}")
        sys.exit(1)
    
    print("🧬 Personal Pharmacogenomics Analysis")
    print("=" * 40)
    
    # Extract variants
    results = extract_pharmaco_variants(vcf_path)
    
    # Generate report
    report = interpret_pharmaco_results(results)
    print(report)
    
    # Save results
    with open('pharmacogenomics_report.txt', 'w') as f:
        f.write("Personal Pharmacogenomics Report\n")
        f.write("=" * 35 + "\n")
        f.write(report)
    
    print(f"\n📄 Report saved to: pharmacogenomics_report.txt")