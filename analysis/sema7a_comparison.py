#!/usr/bin/env python3
"""
SEMA7A Variant Comparison Script
Compare lab results with VCF findings
"""

import subprocess
import re

def extract_sema7a_details():
    """Extract SEMA7A missense variants with full details"""
    print("🔬 SEMA7A Variant Analysis")
    print("=" * 40)
    
    # Search for SEMA7A missense variants
    cmd = [
        'bcftools', 'view', '-r', 'chr15:74000000-76000000',
        '../results/variants_fully_annotated.vcf.gz'
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    sema7a_variants = []
    for line in result.stdout.split('\n'):
        if 'SEMA7A' in line and 'missense' in line:
            fields = line.split('\t')
            if len(fields) >= 10:
                chrom, pos, rsid, ref, alt, qual, filter_field, info, format_field, genotype = fields[:10]
                
                # Parse CSQ field for detailed annotation
                csq_match = re.search(r'CSQ=([^;]+)', info)
                if csq_match:
                    csq_data = csq_match.group(1)
                    csq_fields = csq_data.split('|')
                    
                    if len(csq_fields) > 50:
                        consequence = csq_fields[1]
                        symbol = csq_fields[3]
                        transcript = csq_fields[6]
                        hgvsc = csq_fields[10]
                        hgvsp = csq_fields[11]
                        protein_pos = csq_fields[14]
                        amino_acids = csq_fields[15]
                        sift = csq_fields[36] if len(csq_fields) > 36 else ''
                        polyphen = csq_fields[37] if len(csq_fields) > 37 else ''
                        
                        variant_info = {
                            'position': f"{chrom}:{pos}",
                            'rsid': rsid if rsid != '.' else 'Not in dbSNP',
                            'ref_alt': f"{ref}>{alt}",
                            'genotype': genotype.split(':')[0],
                            'transcript': transcript,
                            'hgvsc': hgvsc,
                            'hgvsp': hgvsp,
                            'protein_position': protein_pos,
                            'amino_acids': amino_acids,
                            'sift': sift,
                            'polyphen': polyphen,
                            'consequence': consequence
                        }
                        sema7a_variants.append(variant_info)
    
    return sema7a_variants

def compare_with_lab_result(variants):
    """Compare findings with lab result"""
    print("\n📋 LAB RESULT:")
    print("Gene: SEMA7A")
    print("Variant: ENST00000261918.9:c.442C>T; p.(Arg148Trp)")
    print("Zygosity: Heterozygous")
    print("Type: missense variant")
    print("dbSNP: rs200895370")
    
    print("\n🔬 VCF FINDINGS:")
    if not variants:
        print("❌ No SEMA7A missense variants found in VCF")
        print("\nPOSSIBLE REASONS:")
        print("1. Variant not called by variant caller")
        print("2. Low coverage in this region")
        print("3. Quality filters removed the variant")
        print("4. Different reference genome coordinates")
        return
    
    print(f"✅ Found {len(variants)} SEMA7A missense variant(s)")
    
    for i, variant in enumerate(variants, 1):
        print(f"\n--- Variant {i} ---")
        print(f"Position: {variant['position']}")
        print(f"Change: {variant['ref_alt']}")
        print(f"Genotype: {variant['genotype']}")
        print(f"Transcript: {variant['transcript']}")
        print(f"cDNA: {variant['hgvsc']}")
        print(f"Protein: {variant['hgvsp']}")
        print(f"Amino Acids: {variant['amino_acids']}")
        print(f"dbSNP: {variant['rsid']}")
        print(f"SIFT: {variant['sift']}")
        print(f"PolyPhen: {variant['polyphen']}")
        
        # Check if this matches the lab result
        is_match = False
        if ('442' in variant['hgvsc'] and 'C>T' in variant['hgvsc']) or \
           ('Arg148Trp' in variant['hgvsp']) or \
           (variant['rsid'] == 'rs200895370'):
            is_match = True
            print("🎯 MATCH: This appears to be the lab-reported variant!")
        
        if variant['genotype'] in ['0/1', '1/0']:
            print("✅ Zygosity: Heterozygous (matches lab)")
        elif variant['genotype'] == '1/1':
            print("⚠️  Zygosity: Homozygous (differs from lab)")
        
        if not is_match:
            print("❓ No obvious match to lab result")

def check_coverage_region():
    """Check if we have coverage in the SEMA7A region"""
    print("\n📊 COVERAGE CHECK:")
    
    # Count total variants in SEMA7A region
    cmd = [
        'bcftools', 'view', '-r', 'chr15:74000000-76000000',
        '../results/variants_fully_annotated.vcf.gz',
        '--no-header'
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    variant_count = len([line for line in result.stdout.split('\n') if line.strip()])
    
    print(f"Total variants in SEMA7A region: {variant_count}")
    
    if variant_count > 0:
        print("✅ Region has good variant coverage")
    else:
        print("❌ No variants in region - possible coverage issue")

if __name__ == "__main__":
    import os
    
    # Set up conda environment
    os.system("source ~/anaconda3/etc/profile.d/conda.sh && conda activate wgs_analysis")
    
    # Extract SEMA7A variants
    variants = extract_sema7a_details()
    
    # Compare with lab results
    compare_with_lab_result(variants)
    
    # Check coverage
    check_coverage_region()
    
    print("\n" + "="*50)
    print("CONCLUSION:")
    
    if variants:
        print("Your VCF contains SEMA7A missense variants.")
        print("Manual review needed to confirm exact match with lab result.")
    else:
        print("Lab variant not found in your VCF data.")
        print("This could indicate:")
        print("- Coverage gaps")
        print("- Quality filtering")
        print("- Reference genome differences")
        print("- Variant calling sensitivity")