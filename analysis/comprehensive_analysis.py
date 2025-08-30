#!/usr/bin/env python3
"""
Comprehensive Personal Genomics Analysis
Processes entire VCF for all insights in one pass
"""

import subprocess
import re
import sys
from collections import defaultdict, Counter

# All analysis targets in one place
ANALYSIS_CONFIG = {
    'pharmacogenomics': {
        'CYP2D6': {'drugs': ['Codeine', 'Tramadol', 'Antidepressants'], 'key_variants': ['rs3892097', 'rs1065852']},
        'CYP2C19': {'drugs': ['Clopidogrel', 'Omeprazole'], 'key_variants': ['rs4244285', 'rs12248560']},
        'CYP2C9': {'drugs': ['Warfarin'], 'key_variants': ['rs1799853', 'rs1057910']},
        'SLCO1B1': {'drugs': ['Statins'], 'key_variants': ['rs4149056']},
        'DPYD': {'drugs': ['5-FU', 'Capecitabine'], 'key_variants': ['rs3918290', 'rs55886062']},
        'TPMT': {'drugs': ['Thiopurines'], 'key_variants': ['rs1800462', 'rs1800460']},
        'UGT1A1': {'drugs': ['Irinotecan'], 'key_variants': ['rs8175347']},
    },
    'carrier_screening': {
        'CFTR': {'disease': 'Cystic Fibrosis', 'frequency': '1/2500'},
        'HBB': {'disease': 'Sickle Cell Disease', 'frequency': '1/365 (African)'},
        'HEXA': {'disease': 'Tay-Sachs Disease', 'frequency': '1/3500 (Ashkenazi)'},
        'SMN1': {'disease': 'Spinal Muscular Atrophy', 'frequency': '1/6000-10000'},
        'GBA': {'disease': 'Gaucher Disease', 'frequency': '1/450 (Ashkenazi)'},
        'ASPA': {'disease': 'Canavan Disease', 'frequency': '1/6400 (Ashkenazi)'},
        'MUTYH': {'disease': 'MUTYH-Associated Polyposis', 'frequency': '1/10000'},
        'HFE': {'disease': 'Hereditary Hemochromatosis', 'frequency': '1/200-400'},
    },
    'acmg_actionable': {
        'BRCA1': {'condition': 'Breast/Ovarian Cancer', 'action': 'Enhanced screening'},
        'BRCA2': {'condition': 'Breast/Ovarian Cancer', 'action': 'Enhanced screening'},
        'MLH1': {'condition': 'Lynch Syndrome', 'action': 'Colonoscopy screening'},
        'MSH2': {'condition': 'Lynch Syndrome', 'action': 'Colonoscopy screening'},
        'MSH6': {'condition': 'Lynch Syndrome', 'action': 'Colonoscopy screening'},
        'PMS2': {'condition': 'Lynch Syndrome', 'action': 'Colonoscopy screening'},
        'APC': {'condition': 'Familial Adenomatous Polyposis', 'action': 'Colonoscopy screening'},
        'VHL': {'condition': 'Von Hippel-Lindau', 'action': 'Multi-organ surveillance'},
        'RET': {'condition': 'Multiple Endocrine Neoplasia', 'action': 'Thyroid surveillance'},
        'TP53': {'condition': 'Li-Fraumeni Syndrome', 'action': 'Enhanced cancer screening'},
    }
}

def extract_all_variants_single_pass(vcf_path):
    """Extract all relevant variants in one bcftools query"""
    print("🔬 Extracting all variants (this may take 2-3 minutes)...")
    
    # Get all genes we care about
    all_genes = set()
    for category in ANALYSIS_CONFIG.values():
        all_genes.update(category.keys())
    
    gene_pattern = '|'.join(all_genes)
    
    # Single bcftools query for everything
    cmd = [
        'bcftools', 'query', 
        '-f', '%CHROM\t%POS\t%ID\t%REF\t%ALT\t[%GT]\t%INFO/CSQ\n',
        vcf_path
    ]
    
    print("Running bcftools query...")
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"Error: {result.stderr}")
        return {}
    
    print("📊 Processing variants...")
    
    # Parse all variants
    variant_data = {
        'pharmacogenomics': defaultdict(list),
        'carrier_screening': defaultdict(list), 
        'acmg_actionable': defaultdict(list),
        'high_impact': [],
        'statistics': Counter()
    }
    
    total_variants = 0
    for line in result.stdout.split('\n'):
        if not line.strip():
            continue
            
        fields = line.split('\t')
        if len(fields) < 7:
            continue
            
        chrom, pos, rsid, ref, alt, genotype, csq = fields
        total_variants += 1
        
        if total_variants % 100000 == 0:
            print(f"  Processed {total_variants:,} variants...")
        
        # Skip reference genotypes
        if genotype in ['0/0', './.']:
            continue
            
        # Parse CSQ field
        csq_fields = csq.split('|') if csq else []
        
        # Extract key info from CSQ
        symbol = ''
        consequence = ''
        impact = ''
        sift = ''
        polyphen = ''
        clin_sig = ''
        
        if len(csq_fields) > 60:  # VEP has many fields
            try:
                symbol = csq_fields[3] if len(csq_fields) > 3 else ''
                consequence = csq_fields[1] if len(csq_fields) > 1 else ''
                impact = csq_fields[2] if len(csq_fields) > 2 else ''
                sift = csq_fields[36] if len(csq_fields) > 36 else ''
                polyphen = csq_fields[37] if len(csq_fields) > 37 else ''
                clin_sig = csq_fields[64] if len(csq_fields) > 64 else ''
            except:
                pass
        
        variant_info = {
            'position': f"{chrom}:{pos}",
            'rsid': rsid if rsid != '.' else None,
            'ref_alt': f"{ref}>{alt}",
            'genotype': genotype,
            'consequence': consequence,
            'impact': impact,
            'sift': sift,
            'polyphen': polyphen,
            'clinical_significance': clin_sig
        }
        
        # Categorize variants
        if symbol in ANALYSIS_CONFIG['pharmacogenomics']:
            variant_data['pharmacogenomics'][symbol].append(variant_info)
            
        if symbol in ANALYSIS_CONFIG['carrier_screening']:
            if 'pathogenic' in clin_sig.lower() or impact == 'HIGH':
                variant_data['carrier_screening'][symbol].append(variant_info)
                
        if symbol in ANALYSIS_CONFIG['acmg_actionable']:
            if 'pathogenic' in clin_sig.lower() or impact == 'HIGH':
                variant_data['acmg_actionable'][symbol].append(variant_info)
        
        # High impact variants
        if impact == 'HIGH':
            variant_data['high_impact'].append(variant_info)
            
        # Statistics
        variant_data['statistics'][impact] += 1
        if 'pathogenic' in clin_sig.lower():
            variant_data['statistics']['pathogenic'] += 1
    
    print(f"✅ Processed {total_variants:,} total variants")
    return variant_data

def generate_comprehensive_report(variant_data):
    """Generate complete personal genomics report"""
    report = []
    report.append("🧬 COMPREHENSIVE PERSONAL GENOMICS ANALYSIS")
    report.append("=" * 50)
    
    # STATISTICS
    report.append(f"\n📊 VARIANT STATISTICS")
    report.append(f"High Impact: {variant_data['statistics']['HIGH']:,}")
    report.append(f"Moderate Impact: {variant_data['statistics']['MODERATE']:,}")
    report.append(f"Low Impact: {variant_data['statistics']['LOW']:,}")
    report.append(f"Pathogenic: {variant_data['statistics']['pathogenic']:,}")
    
    # PHARMACOGENOMICS
    report.append(f"\n💊 PHARMACOGENOMICS RESULTS")
    report.append("-" * 30)
    
    for gene, info in ANALYSIS_CONFIG['pharmacogenomics'].items():
        variants = variant_data['pharmacogenomics'][gene]
        drugs = ', '.join(info['drugs'])
        
        if variants:
            report.append(f"\n{gene} - {drugs}")
            report.append(f"  ⚠️  {len(variants)} variant(s) found")
            for v in variants[:3]:  # Show first 3
                report.append(f"    {v['position']} {v['ref_alt']} ({v['genotype']})")
                if v['consequence']:
                    report.append(f"      Impact: {v['consequence']}")
        else:
            report.append(f"\n{gene} - {drugs}")
            report.append(f"  ✅ No variants affecting metabolism")
    
    # CARRIER SCREENING  
    report.append(f"\n🧪 CARRIER SCREENING RESULTS")
    report.append("-" * 30)
    
    carrier_found = False
    for gene, info in ANALYSIS_CONFIG['carrier_screening'].items():
        variants = variant_data['carrier_screening'][gene]
        disease = info['disease']
        
        if variants:
            carrier_found = True
            report.append(f"\n{gene} - {disease}")
            report.append(f"  ⚠️  CARRIER: {len(variants)} pathogenic variant(s)")
            for v in variants:
                report.append(f"    {v['position']} {v['ref_alt']} ({v['genotype']})")
                if v['clinical_significance']:
                    report.append(f"      Clinical: {v['clinical_significance']}")
    
    if not carrier_found:
        report.append("\n✅ No carrier status detected for common recessive diseases")
    
    # ACMG ACTIONABLE FINDINGS
    report.append(f"\n🏥 ACMG ACTIONABLE FINDINGS")
    report.append("-" * 30)
    
    actionable_found = False
    for gene, info in ANALYSIS_CONFIG['acmg_actionable'].items():
        variants = variant_data['acmg_actionable'][gene]
        condition = info['condition']
        action = info['action']
        
        if variants:
            actionable_found = True
            report.append(f"\n⚠️  {gene} - {condition}")
            report.append(f"  Action: {action}")
            report.append(f"  Variants: {len(variants)}")
            for v in variants:
                report.append(f"    {v['position']} {v['ref_alt']} ({v['genotype']})")
    
    if not actionable_found:
        report.append("\n✅ No pathogenic variants in ACMG actionable genes")
    
    # HIGH IMPACT VARIANTS
    report.append(f"\n⚡ HIGH IMPACT VARIANTS")
    report.append("-" * 30)
    
    if variant_data['high_impact']:
        report.append(f"\nFound {len(variant_data['high_impact'])} high-impact variants:")
        for v in variant_data['high_impact'][:10]:  # Show first 10
            report.append(f"  {v['position']} {v['ref_alt']} ({v['genotype']})")
            if v['consequence']:
                report.append(f"    {v['consequence']}")
    else:
        report.append("\n✅ No high-impact variants detected")
    
    return '\n'.join(report)

if __name__ == "__main__":
    vcf_path = "../results/variants_fully_annotated.vcf.gz"
    
    print("🚀 Starting comprehensive genomics analysis...")
    print("This will analyze your entire genome for:")
    print("  • Drug response variants (pharmacogenomics)")
    print("  • Carrier status for genetic diseases")  
    print("  • ACMG actionable findings")
    print("  • High-impact variants")
    print()
    
    # Extract all variants
    variant_data = extract_all_variants_single_pass(vcf_path)
    
    # Generate report
    report = generate_comprehensive_report(variant_data)
    print(report)
    
    # Save report
    with open('comprehensive_personal_report.txt', 'w') as f:
        f.write(report)
    
    print(f"\n📄 Complete report saved to: comprehensive_personal_report.txt")
    print("🎉 Analysis complete!")