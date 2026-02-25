#!/usr/bin/env python3
"""
Generate synthetic paired FASTQ files for smoke testing.

Creates deterministic paired-end reads with various quality profiles:
- good: High quality reads with few errors
- poor: Low quality reads for edge case testing
- adapter: Reads with adapter contamination
- mixed: Mix of quality profiles

Reads are generated to align to the synthetic reference genome.
"""

import argparse
import gzip
import hashlib
import random
import sys
from pathlib import Path
from typing import Iterator, Tuple


# Quality profiles
QUALITY_PROFILES = {
    "good": {
        "description": "High quality reads (Q30+)",
        "mean_quality": 35,
        "quality_std": 3,
        "error_rate": 0.001,
        "adapter_rate": 0.02,
    },
    "poor": {
        "description": "Low quality reads for edge case testing",
        "mean_quality": 18,
        "quality_std": 8,
        "error_rate": 0.05,
        "adapter_rate": 0.15,
    },
    "adapter": {
        "description": "Heavy adapter contamination",
        "mean_quality": 30,
        "quality_std": 4,
        "error_rate": 0.005,
        "adapter_rate": 0.40,
    },
    "mixed": {
        "description": "Mix of quality profiles",
        "mean_quality": 28,
        "quality_std": 10,
        "error_rate": 0.02,
        "adapter_rate": 0.08,
    },
}

# Illumina adapter sequences (partial)
ILLUMINA_ADAPTER_R1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
ILLUMINA_ADAPTER_R2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

# Default read parameters
DEFAULT_READ_LENGTH = 150
DEFAULT_INSERT_SIZE = 350
DEFAULT_INSERT_STD = 50
DEFAULT_NUM_READS = 2000  # Small for fast smoke tests


def load_reference(ref_path: Path) -> dict[str, str]:
    """Load reference genome FASTA into memory."""
    contigs = {}
    current_name = None
    current_seq = []
    
    with open(ref_path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_name:
                    contigs[current_name] = "".join(current_seq)
                current_name = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        
        if current_name:
            contigs[current_name] = "".join(current_seq)
    
    return contigs


def reverse_complement(seq: str) -> str:
    """Return reverse complement of DNA sequence."""
    complement = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
    return "".join(complement.get(b, "N") for b in reversed(seq))


def quality_char(phred: int) -> str:
    """Convert Phred quality score to ASCII character."""
    return chr(max(0, min(93, phred)) + 33)


def introduce_errors(seq: str, error_rate: float, rng: random.Random) -> str:
    """Introduce sequencing errors into a read."""
    if error_rate <= 0:
        return seq
    
    bases = list(seq)
    substitutes = {"A": "CGT", "T": "ACG", "G": "ACT", "C": "AGT", "N": "ACGT"}
    
    for i in range(len(bases)):
        if rng.random() < error_rate:
            original = bases[i]
            if original in substitutes:
                bases[i] = rng.choice(substitutes[original])
    
    return "".join(bases)


def add_adapter(seq: str, adapter: str, adapter_rate: float, rng: random.Random) -> str:
    """Potentially add adapter sequence to read end."""
    if rng.random() >= adapter_rate:
        return seq
    
    # Add partial adapter (variable length)
    adapter_len = rng.randint(10, min(40, len(adapter)))
    # Trim read and add adapter
    trim_pos = len(seq) - adapter_len
    if trim_pos < 20:
        trim_pos = 20  # Keep at least 20bp of read
    
    return seq[:trim_pos] + adapter[:len(seq) - trim_pos]


def generate_quality_string(
    length: int,
    mean_quality: float,
    quality_std: float,
    rng: random.Random
) -> str:
    """Generate a quality string with position-dependent degradation."""
    quals = []
    for i in range(length):
        # Quality typically degrades toward 3' end
        position_penalty = (i / length) * 5
        q = int(rng.gauss(mean_quality - position_penalty, quality_std))
        q = max(2, min(40, q))  # Clamp to valid range
        quals.append(quality_char(q))
    
    return "".join(quals)


def generate_read_pair(
    read_id: int,
    contigs: dict[str, str],
    profile: dict,
    read_length: int,
    insert_size: int,
    insert_std: int,
    rng: random.Random
) -> Tuple[str, str, str, str, str, str]:
    """Generate a paired-end read from reference.
    
    Returns: (r1_seq, r1_qual, r2_seq, r2_qual, chrom, pos)
    """
    # Select random contig weighted by length
    contig_names = list(contigs.keys())
    contig_lengths = [len(contigs[n]) for n in contig_names]
    total_len = sum(contig_lengths)
    weights = [l / total_len for l in contig_lengths]
    
    chrom = rng.choices(contig_names, weights=weights)[0]
    contig_seq = contigs[chrom]
    
    # Generate insert size
    actual_insert = max(read_length * 2, int(rng.gauss(insert_size, insert_std)))
    
    # Select position ensuring we don't go past contig end
    max_pos = len(contig_seq) - actual_insert - 1
    if max_pos < 0:
        max_pos = len(contig_seq) - read_length - 1
        actual_insert = read_length
    
    pos = rng.randint(0, max(0, max_pos))
    
    # Extract forward and reverse reads
    r1_seq = contig_seq[pos:pos + read_length]
    r2_start = pos + actual_insert - read_length
    r2_seq = contig_seq[r2_start:r2_start + read_length]
    
    # Pad if near contig end
    r1_seq = r1_seq.ljust(read_length, "N")
    r2_seq = r2_seq.ljust(read_length, "N")
    
    # R2 is reverse complement
    r2_seq = reverse_complement(r2_seq)
    
    # Introduce errors
    r1_seq = introduce_errors(r1_seq, profile["error_rate"], rng)
    r2_seq = introduce_errors(r2_seq, profile["error_rate"], rng)
    
    # Add adapters
    r1_seq = add_adapter(r1_seq, ILLUMINA_ADAPTER_R1, profile["adapter_rate"], rng)
    r2_seq = add_adapter(r2_seq, ILLUMINA_ADAPTER_R2, profile["adapter_rate"], rng)
    
    # Generate quality strings
    r1_qual = generate_quality_string(
        read_length, profile["mean_quality"], profile["quality_std"], rng
    )
    r2_qual = generate_quality_string(
        read_length, profile["mean_quality"], profile["quality_std"], rng
    )
    
    return r1_seq, r1_qual, r2_seq, r2_qual, chrom, str(pos)


def write_fastq_record(
    f,
    read_id: str,
    sequence: str,
    quality: str,
    extra_info: str = ""
) -> None:
    """Write a single FASTQ record."""
    header = f"@{read_id}"
    if extra_info:
        header += f" {extra_info}"
    f.write(f"{header}\n{sequence}\n+\n{quality}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Generate synthetic paired FASTQ files for smoke testing"
    )
    parser.add_argument(
        "-r", "--reference",
        type=Path,
        required=True,
        help="Reference genome FASTA file"
    )
    parser.add_argument(
        "-o", "--output-dir",
        type=Path,
        default=Path("tests/smoke/data/reads"),
        help="Output directory for FASTQ files"
    )
    parser.add_argument(
        "-n", "--num-reads",
        type=int,
        default=DEFAULT_NUM_READS,
        help=f"Number of read pairs to generate (default: {DEFAULT_NUM_READS})"
    )
    parser.add_argument(
        "-l", "--read-length",
        type=int,
        default=DEFAULT_READ_LENGTH,
        help=f"Read length (default: {DEFAULT_READ_LENGTH})"
    )
    parser.add_argument(
        "-p", "--profile",
        choices=list(QUALITY_PROFILES.keys()),
        default="good",
        help="Quality profile"
    )
    parser.add_argument(
        "--sample-name",
        default="smoke_sample",
        help="Sample name prefix"
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for reproducibility"
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Verbose output"
    )
    parser.add_argument(
        "--all-profiles",
        action="store_true",
        help="Generate reads for all quality profiles"
    )
    
    args = parser.parse_args()
    
    # Load reference
    if not args.reference.exists():
        print(f"Error: Reference file not found: {args.reference}", file=sys.stderr)
        return 1
    
    if args.verbose:
        print(f"Loading reference: {args.reference}")
    
    contigs = load_reference(args.reference)
    total_ref_len = sum(len(s) for s in contigs.values())
    
    if args.verbose:
        print(f"  Loaded {len(contigs)} contigs, {total_ref_len:,} bp total")
    
    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)
    
    # Determine profiles to generate
    profiles_to_run = list(QUALITY_PROFILES.keys()) if args.all_profiles else [args.profile]
    
    for profile_name in profiles_to_run:
        profile = QUALITY_PROFILES[profile_name]
        
        if args.verbose:
            print(f"\nGenerating {profile_name} reads: {profile['description']}")
        
        rng = random.Random(args.seed)
        
        # Output file paths
        sample_id = f"{args.sample_name}_{profile_name}"
        r1_path = args.output_dir / f"{sample_id}_R1.fastq.gz"
        r2_path = args.output_dir / f"{sample_id}_R2.fastq.gz"
        
        # Generate reads
        with gzip.open(r1_path, "wt") as f1, gzip.open(r2_path, "wt") as f2:
            for i in range(args.num_reads):
                r1_seq, r1_qual, r2_seq, r2_qual, chrom, pos = generate_read_pair(
                    read_id=i,
                    contigs=contigs,
                    profile=profile,
                    read_length=args.read_length,
                    insert_size=DEFAULT_INSERT_SIZE,
                    insert_std=DEFAULT_INSERT_STD,
                    rng=rng
                )
                
                read_name = f"{sample_id}:{i}:1"
                write_fastq_record(f1, f"{read_name}/1", r1_seq, r1_qual)
                write_fastq_record(f2, f"{read_name}/2", r2_seq, r2_qual)
        
        # Compute stats
        r1_size = r1_path.stat().st_size
        r2_size = r2_path.stat().st_size
        
        print(f"✓ Generated {profile_name} FASTQ pair:")
        print(f"    R1: {r1_path} ({r1_size:,} bytes)")
        print(f"    R2: {r2_path} ({r2_size:,} bytes)")
        print(f"    Reads: {args.num_reads:,} pairs")
        print(f"    Profile: {profile['description']}")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
