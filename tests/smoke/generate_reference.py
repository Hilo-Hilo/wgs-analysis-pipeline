#!/usr/bin/env python3
"""
Generate a tiny synthetic reference genome for smoke testing.

Creates a deterministic ~100kb genome with:
- Multiple "chromosomes" (contigs)
- Known variant positions for validation
- Predictable sequence patterns for debugging

Output: FASTA file ready for BWA-MEM2 indexing
"""

import argparse
import hashlib
import os
import sys
from pathlib import Path


# Deterministic seed for reproducibility
SEED = 42

# Contig configurations: (name, length, GC_content)
CONTIGS = [
    ("chr1_mini", 40000, 0.45),
    ("chr2_mini", 30000, 0.50),
    ("chr3_mini", 20000, 0.42),
    ("chrM_mini", 10000, 0.40),
]

# Known variant positions for validation (contig_idx, position, ref, alt)
# These will be embedded in the reference for test variant calling
KNOWN_VARIANTS = [
    (0, 1000, "A", "G"),    # SNP on chr1_mini
    (0, 5000, "T", "C"),    # SNP on chr1_mini
    (0, 10000, "GG", "G"),  # Deletion on chr1_mini
    (1, 2000, "C", "T"),    # SNP on chr2_mini
    (1, 8000, "A", "AT"),   # Insertion on chr2_mini
    (2, 5000, "G", "A"),    # SNP on chr3_mini
]


def seeded_random(seed_state: int) -> tuple[int, float]:
    """Simple LCG random number generator for reproducibility."""
    # Linear congruential generator parameters
    a = 1664525
    c = 1013904223
    m = 2**32
    new_state = (a * seed_state + c) % m
    return new_state, new_state / m


def generate_sequence(length: int, gc_content: float, seed: int) -> str:
    """Generate a random DNA sequence with specified GC content."""
    bases = []
    state = seed
    
    for _ in range(length):
        state, r = seeded_random(state)
        
        if r < gc_content / 2:
            bases.append("G")
        elif r < gc_content:
            bases.append("C")
        elif r < gc_content + (1 - gc_content) / 2:
            bases.append("A")
        else:
            bases.append("T")
    
    return "".join(bases)


def write_fasta(output_path: Path, contigs: list[tuple[str, str]]) -> None:
    """Write contigs to FASTA file with proper line wrapping."""
    with open(output_path, "w") as f:
        for name, sequence in contigs:
            f.write(f">{name}\n")
            # Wrap at 80 characters
            for i in range(0, len(sequence), 80):
                f.write(sequence[i:i+80] + "\n")


def write_variant_truth(output_path: Path, contigs: list[tuple[str, str]]) -> None:
    """Write known variant positions for test validation."""
    with open(output_path, "w") as f:
        f.write("# Known variant positions for smoke test validation\n")
        f.write("# CHROM\tPOS\tREF\tALT\n")
        for contig_idx, pos, ref, alt in KNOWN_VARIANTS:
            chrom = contigs[contig_idx][0]
            f.write(f"{chrom}\t{pos}\t{ref}\t{alt}\n")


def compute_checksum(filepath: Path) -> str:
    """Compute MD5 checksum for reproducibility verification."""
    md5 = hashlib.md5()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            md5.update(chunk)
    return md5.hexdigest()


def main():
    parser = argparse.ArgumentParser(
        description="Generate synthetic reference genome for smoke testing"
    )
    parser.add_argument(
        "-o", "--output-dir",
        type=Path,
        default=Path("tests/smoke/data/reference"),
        help="Output directory for reference files"
    )
    parser.add_argument(
        "--name",
        default="smoke_ref",
        help="Reference genome name prefix"
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=SEED,
        help="Random seed for reproducibility"
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Verbose output"
    )
    
    args = parser.parse_args()
    
    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)
    
    if args.verbose:
        print(f"Generating synthetic reference genome...")
        print(f"  Output directory: {args.output_dir}")
        print(f"  Seed: {args.seed}")
    
    # Generate contig sequences
    contigs = []
    state = args.seed
    total_length = 0
    
    for name, length, gc in CONTIGS:
        state, _ = seeded_random(state)
        seq = generate_sequence(length, gc, state)
        contigs.append((name, seq))
        total_length += length
        
        if args.verbose:
            print(f"  Generated {name}: {length}bp, GC={gc:.0%}")
    
    # Write FASTA
    fasta_path = args.output_dir / f"{args.name}.fa"
    write_fasta(fasta_path, contigs)
    
    # Write known variants truth file
    truth_path = args.output_dir / f"{args.name}_variants.truth"
    write_variant_truth(truth_path, contigs)
    
    # Compute checksum for reproducibility
    checksum = compute_checksum(fasta_path)
    
    # Write manifest
    manifest_path = args.output_dir / f"{args.name}.manifest"
    with open(manifest_path, "w") as f:
        f.write(f"# Smoke test reference genome manifest\n")
        f.write(f"name={args.name}\n")
        f.write(f"seed={args.seed}\n")
        f.write(f"total_length={total_length}\n")
        f.write(f"num_contigs={len(contigs)}\n")
        f.write(f"md5={checksum}\n")
        f.write(f"fasta={fasta_path.name}\n")
        f.write(f"truth={truth_path.name}\n")
    
    print(f"✓ Reference genome generated: {fasta_path}")
    print(f"  Total length: {total_length:,} bp")
    print(f"  Contigs: {len(contigs)}")
    print(f"  MD5: {checksum}")
    print(f"  Variants truth: {truth_path}")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
