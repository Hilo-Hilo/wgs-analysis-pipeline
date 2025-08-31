#!/usr/bin/env python3
"""
Generate sample FASTQ data for testing the WGS pipeline.
Creates realistic but minimal test datasets.
"""

import argparse
import gzip
import random
import os
from pathlib import Path

def generate_realistic_read(read_id, length=150, quality_profile='good'):
    """Generate a realistic FASTQ read with varying quality patterns"""
    
    # Base composition (slightly GC-biased like real data)
    bases = ['A'] * 25 + ['T'] * 25 + ['G'] * 25 + ['C'] * 25
    
    # Generate sequence
    sequence = ''.join(random.choice(bases) for _ in range(length))
    
    # Generate quality scores based on profile
    if quality_profile == 'good':
        # High quality at start, slight decline toward end
        quality_scores = []
        for i in range(length):
            if i < length * 0.1:  # First 10% - very high quality
                base_q = random.randint(35, 40)
            elif i < length * 0.8:  # Middle 70% - good quality  
                base_q = random.randint(30, 38)
            else:  # Last 20% - declining quality
                decline = (i - length * 0.8) / (length * 0.2)
                base_q = random.randint(20, 35 - int(decline * 10))
            
            quality_scores.append(chr(base_q + 33))
    
    elif quality_profile == 'poor':
        # Lower quality throughout
        quality_scores = [chr(random.randint(15, 30) + 33) for _ in range(length)]
    
    else:  # mixed
        quality_scores = [chr(random.randint(20, 35) + 33) for _ in range(length)]
    
    quality = ''.join(quality_scores)
    
    return f"@{read_id}\n{sequence}\n+\n{quality}\n"

def generate_fastq_pair(output_dir, sample_name, num_reads=10000, read_length=150, quality_mix=True):
    """Generate paired FASTQ files"""
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    r1_file = output_dir / f"{sample_name}_R1.fastq.gz"
    r2_file = output_dir / f"{sample_name}_R2.fastq.gz"
    
    print(f"Generating {num_reads} paired reads...")
    print(f"Output: {r1_file}")
    print(f"Output: {r2_file}")
    
    with gzip.open(r1_file, 'wt') as f1, gzip.open(r2_file, 'wt') as f2:
        for i in range(num_reads):
            # Vary quality profiles for more realistic data
            if quality_mix:
                if i < num_reads * 0.7:
                    profile = 'good'
                elif i < num_reads * 0.9:
                    profile = 'mixed'
                else:
                    profile = 'poor'
            else:
                profile = 'good'
            
            # Generate paired reads
            read_id_1 = f"read_{i:06d}/1"
            read_id_2 = f"read_{i:06d}/2"
            
            r1_read = generate_realistic_read(read_id_1, read_length, profile)
            r2_read = generate_realistic_read(read_id_2, read_length, profile)
            
            f1.write(r1_read)
            f2.write(r2_read)
    
    # Report file sizes
    r1_size = r1_file.stat().st_size / (1024*1024)  # MB
    r2_size = r2_file.stat().st_size / (1024*1024)  # MB
    print(f"Generated files:")
    print(f"  {r1_file.name}: {r1_size:.1f} MB")
    print(f"  {r2_file.name}: {r2_size:.1f} MB")
    
    return r1_file, r2_file

def create_test_datasets():
    """Create multiple test datasets for different scenarios"""
    
    base_dir = Path(__file__).parent / "sample_data"
    
    datasets = [
        # Small dataset for quick tests
        {
            'name': 'small_test',
            'reads': 1000,
            'length': 100,
            'description': 'Small dataset for quick testing'
        },
        # Medium dataset for integration tests
        {
            'name': 'medium_test', 
            'reads': 10000,
            'length': 150,
            'description': 'Medium dataset for integration testing'
        },
        # Large dataset for performance tests
        {
            'name': 'large_test',
            'reads': 100000,
            'length': 150, 
            'description': 'Large dataset for performance testing'
        }
    ]
    
    for dataset in datasets:
        print(f"\n=== Creating {dataset['description']} ===")
        output_dir = base_dir / dataset['name']
        
        generate_fastq_pair(
            output_dir=output_dir,
            sample_name=dataset['name'],
            num_reads=dataset['reads'],
            read_length=dataset['length']
        )
        
        # Create a README for this dataset
        readme_file = output_dir / "README.md"
        with open(readme_file, 'w') as f:
            f.write(f"# {dataset['description']}\n\n")
            f.write(f"- **Reads**: {dataset['reads']:,} paired-end\n")
            f.write(f"- **Length**: {dataset['length']} bp\n")
            f.write(f"- **Purpose**: {dataset['description']}\n")
            f.write(f"- **Generated**: Synthetic test data\n\n")
            f.write("## Usage\n\n")
            f.write("```bash\n")
            f.write(f"# Use this dataset for testing:\n")
            f.write(f"./scripts/quality_control.sh --input-dir tests/sample_data/{dataset['name']}\n")
            f.write("```\n")

def main():
    parser = argparse.ArgumentParser(description='Generate sample FASTQ data for testing')
    parser.add_argument('--output-dir', '-o', default='sample_data',
                       help='Output directory (default: sample_data)')
    parser.add_argument('--sample-name', '-s', default='test_sample',
                       help='Sample name prefix (default: test_sample)')
    parser.add_argument('--num-reads', '-n', type=int, default=10000,
                       help='Number of paired reads (default: 10000)')
    parser.add_argument('--read-length', '-l', type=int, default=150,
                       help='Read length in bp (default: 150)')
    parser.add_argument('--create-all', action='store_true',
                       help='Create all standard test datasets')
    
    args = parser.parse_args()
    
    if args.create_all:
        create_test_datasets()
    else:
        generate_fastq_pair(
            output_dir=args.output_dir,
            sample_name=args.sample_name,
            num_reads=args.num_reads,
            read_length=args.read_length
        )
    
    print("\n✅ Sample data generation complete!")

if __name__ == '__main__':
    main()