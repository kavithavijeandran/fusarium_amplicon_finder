#!/usr/bin/env python

import argparse
import os
import csv
import re
import hashlib
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq

# IUPAC degenerate base codes
IUPAC_CODES = {
    'A': ['A'],
    'C': ['C'],
    'G': ['G'],
    'T': ['T'],
    'R': ['A', 'G'],
    'Y': ['C', 'T'],
    'S': ['G', 'C'],
    'W': ['A', 'T'],
    'K': ['G', 'T'],
    'M': ['A', 'C'],
    'B': ['C', 'G', 'T'],
    'D': ['A', 'G', 'T'],
    'H': ['A', 'C', 'T'],
    'V': ['A', 'C', 'G'],
    'N': ['A', 'C', 'G', 'T']
}

def rev_comp(seq):
    """Return reverse complement of a sequence."""
    return str(Seq(seq).reverse_complement())

def matches_with_degeneracy(base, primer_base):
    """Check if a genome base matches a primer base (accounting for degeneracy)."""
    primer_base = primer_base.upper()
    base = base.upper()
    
    if primer_base not in IUPAC_CODES:
        # If invalid code, fall back to exact match
        return base == primer_base
    
    return base in IUPAC_CODES[primer_base]

def find_approx_matches(sequence, primer, max_mismatches=5):
    """
    Find approximate matches of a primer in a sequence, accounting for degenerate bases.
    
    Args:
        sequence: Target DNA sequence (string)
        primer: Primer sequence with possible degenerate bases (string)
        max_mismatches: Maximum number of mismatches allowed (int)
    
    Returns:
        List of tuples: (position, mismatch_count)
    """
    matches = []
    primer = primer.upper()
    seq_len = len(sequence)
    primer_len = len(primer)
    
    for i in range(seq_len - primer_len + 1):
        window = sequence[i:i + primer_len]
        mismatches = 0
        
        for genome_base, primer_base in zip(window, primer):
            if not matches_with_degeneracy(genome_base, primer_base):
                mismatches += 1
                if mismatches > max_mismatches:
                    break  # Early exit if too many mismatches
        
        if mismatches <= max_mismatches:
            matches.append((i, mismatches))
    
    return matches

def extract_species_name(genome_file):
    """Extract species name from genome FASTA header."""
    try:
        with open(genome_file, "r") as f:
            for line in f:
                if line.startswith(">"):
                    parts = line.strip().split(" ")
                    if len(parts) >= 4:
                        return " ".join(parts[1:4])
                    elif len(parts) >= 2:
                        return " ".join(parts[1:])
                    else:
                        return parts[0].lstrip(">")
    except Exception as e:
        print(f"Warning: Could not extract species name from {genome_file}: {e}")
    
    return "Unknown_Species"

def extract_accession_number(genome_file):
    """Extract accession number from genome filename."""
    filename = os.path.basename(genome_file)
    parts = filename.split("_")
    if len(parts) >= 2:
        return parts[0] + "_" + parts[1]
    return filename.replace(".fna", "").replace(".fasta", "")

def detect_panel(genome_file):
    """Detect whether genome belongs to core or extended panel."""
    genome_file_normalized = genome_file.replace("\\", "/")
    if "/core/" in genome_file_normalized:
        return "core"
    elif "/extended/" in genome_file_normalized:
        return "extended"
    else:
        return "unknown"

def hash_sequence(seq):
    """Create a hash of a sequence for memory-efficient deduplication."""
    return hashlib.md5(seq.encode()).hexdigest()

def screen_genome(genome_file, forward_primer, reverse_primer, min_size, max_size, 
                  max_mismatches, deduplicate, seen):
    """
    Screen a single genome for amplicons matching the primer pair.
    
    Args:
        genome_file: Path to genome FASTA file
        forward_primer: Forward primer sequence (can contain degenerate bases)
        reverse_primer: Reverse primer sequence (can contain degenerate bases)
        min_size: Minimum amplicon size
        max_size: Maximum amplicon size
        max_mismatches: Maximum mismatches allowed per primer
        deduplicate: Whether to deduplicate results
        seen: Set for tracking seen amplicons (if deduplicating)
    
    Returns:
        List of result tuples
    """
    forward_primer_rc = rev_comp(forward_primer)
    reverse_primer_rc = rev_comp(reverse_primer)
    results = []

    species_name = extract_species_name(genome_file)
    accession = extract_accession_number(genome_file)
    panel = detect_panel(genome_file)

    # More intuitive orientation labels
    orientation_map = {
        0: "Plus_Strand_F_to_R",      # Forward primer on + strand, reverse on + strand
        1: "Plus_Strand_R_to_F",      # Reverse primer on + strand, forward on + strand
        2: "Minus_Strand_F_to_R",     # Forward primer on - strand, reverse on - strand
        3: "Minus_Strand_R_to_F"      # Reverse primer on - strand, forward on - strand
    }

    try:
        for record in SeqIO.parse(genome_file, "fasta"):
            seq = str(record.seq).upper()
            rev_seq = rev_comp(seq)
            contig = record.id
            contig_length = len(seq)

            # Check all four possible primer orientations on both strands
            primer_combinations = [
                # Plus strand: F ? R
                (find_approx_matches(seq, forward_primer, max_mismatches),
                 find_approx_matches(seq, reverse_primer_rc, max_mismatches),
                 seq, forward_primer, reverse_primer, 0),
                
                # Plus strand: R ? F
                (find_approx_matches(seq, reverse_primer, max_mismatches),
                 find_approx_matches(seq, forward_primer_rc, max_mismatches),
                 seq, reverse_primer, forward_primer_rc, 1),
                
                # Minus strand: F ? R
                (find_approx_matches(rev_seq, forward_primer_rc, max_mismatches),
                 find_approx_matches(rev_seq, reverse_primer, max_mismatches),
                 rev_seq, forward_primer_rc, reverse_primer, 2),
                
                # Minus strand: R ? F
                (find_approx_matches(rev_seq, reverse_primer_rc, max_mismatches),
                 find_approx_matches(rev_seq, forward_primer, max_mismatches),
                 rev_seq, reverse_primer_rc, forward_primer, 3)
            ]

            for fwd_positions, rev_positions, genome_seq, fwd_primer, rev_primer, orient_idx in primer_combinations:
                orientation = orientation_map[orient_idx]
                
                for f_pos, f_mismatches in fwd_positions:
                    for r_pos, r_mismatches in rev_positions:
                        # Check if primers are in correct order and within mismatch threshold
                        if r_pos > f_pos and f_mismatches <= max_mismatches and r_mismatches <= max_mismatches:
                            amplicon_end = r_pos + len(rev_primer)
                            amplicon_size = amplicon_end - f_pos
                            
                            if min_size <= amplicon_size <= max_size:
                                amplicon_seq = genome_seq[f_pos:amplicon_end]
                                forward_site_seq = genome_seq[f_pos:f_pos + len(fwd_primer)]
                                reverse_site_seq = genome_seq[r_pos:r_pos + len(rev_primer)]

                                # Memory-efficient deduplication using sequence hash
                                if deduplicate:
                                    seq_hash = hash_sequence(amplicon_seq)
                                    dedup_key = (species_name, accession, contig, f_pos + 1, amplicon_end, seq_hash)
                                    if dedup_key in seen:
                                        continue
                                    seen.add(dedup_key)

                                mismatch_tier = max(f_mismatches, r_mismatches)

                                results.append((
                                    species_name, accession, contig, contig_length,
                                    f_pos + 1, amplicon_end, amplicon_seq,
                                    amplicon_size,
                                    forward_site_seq, reverse_site_seq, orientation,
                                    f_mismatches, r_mismatches, mismatch_tier, panel
                                ))
    
    except Exception as e:
        print(f"Error processing {genome_file}: {e}")
    
    return results

def main():
    parser = argparse.ArgumentParser(
        description="Screen genomes for amplicons using primer pairs. Supports degenerate IUPAC codes.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s -d /path/to/genomes -f ATCGATCG -r GCTAGCTA
  %(prog)s -d genomes/ -f GARTGCCTAAAWTGG -r TCCTCCGCTTATTGATATGC --forward_name ITS1 --reverse_name ITS4

Degenerate codes supported: R, Y, S, W, K, M, B, D, H, V, N
        """
    )
    
    parser.add_argument("-d", "--directory", required=True, 
                        help="Directory containing genome FASTA files (searches core/ and extended/ subdirs)")
    parser.add_argument("-f", "--forward", required=True, 
                        help="Forward primer sequence (supports IUPAC degenerate codes)")
    parser.add_argument("-r", "--reverse", required=True, 
                        help="Reverse primer sequence (supports IUPAC degenerate codes)")
    parser.add_argument("--forward_name", default="", 
                        help="Forward primer name (e.g., ITS1)")
    parser.add_argument("--reverse_name", default="", 
                        help="Reverse primer name (e.g., ITS4)")
    parser.add_argument("--min", type=int, default=100, 
                        help="Minimum amplicon size in bp (default: 100)")
    parser.add_argument("--max", type=int, default=5000, 
                        help="Maximum amplicon size in bp (default: 5000)")
    parser.add_argument("--mismatches", type=int, default=5, 
                        help="Maximum mismatches allowed per primer (default: 5)")
    parser.add_argument("--no_dedup", action="store_true", 
                        help="Disable deduplication of identical amplicons")
    parser.add_argument("--outdir", default=None,
                        help="Output directory for results (default: ./amplicon_results)")
    parser.add_argument("--progress", action="store_true",
                        help="Show progress for each genome processed")
    
    args = parser.parse_args()

    # Set output directory
    if args.outdir is None:
        result_dir = os.path.join(os.getcwd(), "amplicon_results")
    else:
        result_dir = args.outdir
    
    amplicon_result_dir = os.path.join(result_dir, "amplicon_finder")
    os.makedirs(amplicon_result_dir, exist_ok=True)

    # Sanitize sequences for filename
    def sanitize(seq):
        return re.sub(r"[^A-Za-z0-9]", "_", seq)[:50]  # Limit length

    fwd_clean = sanitize(args.forward)
    rev_clean = sanitize(args.reverse)
    fwd_name_clean = sanitize(args.forward_name)
    rev_name_clean = sanitize(args.reverse_name)

    if args.forward_name and args.reverse_name:
        filename = f"amplicon_results_{fwd_name_clean}-{rev_name_clean}.csv"
    else:
        filename = f"amplicon_results_{fwd_clean}_{rev_clean}.csv"

    amplicon_csv = os.path.join(amplicon_result_dir, filename)

    # Find genome files (only in core/ and extended/ subdirectories)
    genome_files = []
    for root, _, files in os.walk(args.directory):
        root_normalized = root.replace("\\", "/")
        if "core" in root_normalized or "extended" in root_normalized:
            for f in files:
                if f.endswith(".fna") or f.endswith(".fasta"):
                    genome_files.append(os.path.join(root, f))

    print(f"Found {len(genome_files)} genome files to process.")
    print(f"Forward primer: {args.forward}")
    print(f"Reverse primer: {args.reverse}")
    print(f"Amplicon size range: {args.min}-{args.max} bp")
    print(f"Max mismatches per primer: {args.mismatches}")
    print(f"Deduplication: {'OFF' if args.no_dedup else 'ON'}")
    print()

    grouped_results = defaultdict(list)
    seen = set()  # Global deduplication tracker
    total_amplicons = 0

    for idx, genome_file in enumerate(genome_files, 1):
        if args.progress:
            print(f"Processing {idx}/{len(genome_files)}: {os.path.basename(genome_file)}")
        
        results = screen_genome(
            genome_file, args.forward, args.reverse,
            args.min, args.max, args.mismatches, not args.no_dedup, seen
        )
        
        if results:
            grouped_results[genome_file].extend(results)
            total_amplicons += len(results)

    print(f"\nFound {total_amplicons} total amplicons across {len(grouped_results)} genomes.")

    # Write results to CSV
    with open(amplicon_csv, "w", newline="") as csvfile:
        csv_writer = csv.writer(csvfile)
        
        # Write metadata
        fwd_label = f"Forward Primer ({args.forward_name})" if args.forward_name else "Forward Primer"
        rev_label = f"Reverse Primer ({args.reverse_name})" if args.reverse_name else "Reverse Primer"
        csv_writer.writerow([fwd_label, args.forward])
        csv_writer.writerow([rev_label, args.reverse])
        csv_writer.writerow(["Min Size", args.min])
        csv_writer.writerow(["Max Size", args.max])
        csv_writer.writerow(["Max Mismatches", args.mismatches])
        csv_writer.writerow([])
        
        # Write header
        csv_writer.writerow([
            "Species Name", "Accession", "Contig", "Contig Length",
            "Start", "End", "Amplicon Sequence", "Amplicon Size",
            "Forward Primer Match Site", "Reverse Primer Match Site", "Orientation",
            "Forward Primer Mismatches", "Reverse Primer Mismatches",
            "Mismatch Tier", "Panel"
        ])

        # Write data
        for genome_file in sorted(grouped_results.keys()):
            for row in grouped_results[genome_file]:
                csv_writer.writerow(row)

    print(f"\nAmplicon results saved to: {amplicon_csv}")

if __name__ == '__main__':
    main()
