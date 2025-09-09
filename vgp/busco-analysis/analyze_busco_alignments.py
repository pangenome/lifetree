#!/usr/bin/env python3

import sys
import subprocess
import argparse
from collections import defaultdict
import tempfile
import os
from multiprocessing import Pool
from functools import partial

def parse_bed_file(bed_file):
    """Parse BUSCO BED file and return dict of gene_id -> (chr, start, end, strand)"""
    busco_genes = {}
    with open(bed_file, 'r') as f:
        for line in f:
            if line.strip():
                parts = line.strip().split('\t')
                busco_genes[parts[3]] = (parts[0], int(parts[1]), int(parts[2]), 
                                         parts[5] if len(parts) > 5 else '+')
    return busco_genes

def run_impg_query(paf_file, target_chrom, target_start, target_end):
    """Run impg query for a specific region"""
    try:
        result = subprocess.run(
            ['impg', 'query', '-p', paf_file, '-r', f'{target_chrom}:{target_start}-{target_end}', 
             '-v', '1', '-d', '0', '-o', 'paf'],
            capture_output=True, text=True, check=True
        )
        return result.stdout
    except subprocess.CalledProcessError:
        return ""

def parse_impg_output(impg_output):
    """Parse impg PAF output and extract alignment information"""
    alignments = []
    for line in impg_output.strip().split('\n'):
        if not line or line.startswith('['):
            continue
        parts = line.split('\t')
        if len(parts) >= 12:
            alignments.append({
                'query_seq': parts[0],
                'query_start': int(parts[2]),
                'query_end': int(parts[3]),
                'target_start': int(parts[7]),
                'target_end': int(parts[8]),
            })
    return alignments

def analyze_gene(gene_id, target_busco, query_busco, paf_file):
    """Analyze a single BUSCO gene alignment"""
    target_chrom, target_start, target_end, _ = target_busco[gene_id]
    
    impg_output = run_impg_query(paf_file, target_chrom, target_start, target_end)
    alignments = parse_impg_output(impg_output)
    
    result = {
        'gene_id': gene_id,
        'target_name': target_chrom,
        'target_start': target_start,
        'target_end': target_end,
        'target_length': target_end - target_start,
        'gene_present_in_query': gene_id in query_busco,
        'query_n_total_alignments': len(alignments),
        'query_n_busco_alignments': 0,
        'query_n_non_busco_alignments': 0,
        'query_busco_completeness': 0.0
    }
    
    if alignments:
        # Check alignments overlapping with BUSCO location and calculate query coverage
        if gene_id in query_busco:
            query_chrom, query_start, query_end, _ = query_busco[gene_id]
            busco_alignments = []
            
            for aln in alignments:
                if (aln['query_seq'] == query_chrom and 
                    not (aln['query_end'] <= query_start or aln['query_start'] >= query_end)):
                    result['query_n_busco_alignments'] += 1
                    busco_alignments.append(aln)
                else:
                    result['query_n_non_busco_alignments'] += 1
            
            # Calculate coverage of query BUSCO region
            result['query_busco_completeness'] = calculate_coverage(busco_alignments, query_start, query_end, 'query')
        else:
            result['query_n_non_busco_alignments'] = len(alignments)
    
    return result

def calculate_coverage(alignments, gene_start, gene_end, coord_field='target'):
    """Calculate how much of the gene is covered by alignments"""
    if not alignments:
        return 0.0
    
    # Extract and clip intervals to gene boundaries
    intervals = []
    for aln in alignments:
        start = max(aln[f'{coord_field}_start'], gene_start)
        end = min(aln[f'{coord_field}_end'], gene_end)
        if start < end:
            intervals.append((start, end))
    
    if not intervals:
        return 0.0
    
    # Sort intervals by start position
    intervals.sort()
    
    # Merge overlapping intervals
    merged = [intervals[0]]
    for start, end in intervals[1:]:
        if start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    
    # Calculate total coverage
    total_covered = sum(end - start for start, end in merged)
    return total_covered / (gene_end - gene_start)

def main():
    parser = argparse.ArgumentParser(description='Analyze BUSCO gene alignments between genomes')
    parser.add_argument('--target-bed', required=True, help='BUSCO BED file for target genome')
    parser.add_argument('--query-bed', required=True, help='BUSCO BED file for query genome')
    parser.add_argument('--paf', required=True, help='PAF alignment file')
    parser.add_argument('--output', help='Output file for per-gene results (default: stdout)')
    parser.add_argument('--summary-output', help='Output file for summary statistics')
    parser.add_argument('--gene', help='Analyze only this specific BUSCO gene')
    parser.add_argument('--threads', type=int, default=4, help='Number of threads (default: 4)')
    
    args = parser.parse_args()
    
    # Parse BUSCO annotations
    print("Loading BUSCO annotations...", file=sys.stderr)
    target_busco = parse_bed_file(args.target_bed)
    query_busco = parse_bed_file(args.query_bed)
    
    print(f"Found {len(target_busco)} BUSCO genes in target", file=sys.stderr)
    print(f"Found {len(query_busco)} BUSCO genes in query", file=sys.stderr)
    
    # Determine genes to analyze (input gene or all common genes)
    genes_to_analyze = [args.gene] if args.gene else sorted(set(target_busco.keys()) & set(query_busco.keys()))

    # Process genes in parallel
    print(f"Processing {len(genes_to_analyze)} genes using {args.threads} threads...", file=sys.stderr)
    with Pool(args.threads) as pool:
        analyze_func = partial(analyze_gene, target_busco=target_busco, 
                              query_busco=query_busco, paf_file=args.paf)
        results = pool.map(analyze_func, genes_to_analyze)
    
    # Output per-gene results
    output = open(args.output, 'w') if args.output else sys.stdout
    
    # Write header
    output.write('\t'.join([
        'gene_id', 'target_name', 'target_start', 'target_end', 'target_length',
        'gene_present_in_query', 'gene_present_in_alignments', 'query_n_total_alignments',
        'query_n_busco_alignments', 'query_n_non_busco_alignments', 'query_busco_completeness'
    ]) + '\n')
    
    # Write results and collect stats
    stats = {'missing': 0, 'incomplete': 0, 'wrong': 0, 'good': 0}
    
    for r in results:
        output.write('\t'.join([
            r['gene_id'], r['target_name'], str(r['target_start']), str(r['target_end']),
            str(r['target_length']),
            'Yes' if r['gene_present_in_query'] else 'No', 'Yes' if r['query_busco_completeness'] > 0.0 else 'No',
            str(r['query_n_total_alignments']), str(r['query_n_busco_alignments']), str(r['query_n_non_busco_alignments']),
            f"{r['query_busco_completeness']:.4f}"
        ]) + '\n')
        
        if r['query_busco_completeness'] == 0.0:
            stats['missing'] += 1
        elif r['query_busco_completeness'] < 0.9:
            stats['incomplete'] += 1
        elif r['query_n_non_busco_alignments'] > 0:
            stats['wrong'] += 1
        else:
            stats['good'] += 1
    
    if args.output:
        output.close()
    
    # Write summary statistics to file if requested
    if args.summary_output:
        with open(args.summary_output, 'w') as summary_file:
            # Write header
            summary_file.write('\t'.join(['total_analyzed', 'missing', 'incomplete', 'wrong_location_too', 'good']) + '\n')
            # Write data
            summary_file.write('\t'.join([
                str(len(results)),
                str(stats['missing']),
                str(stats['incomplete']),
                str(stats['wrong']),
                str(stats['good'])
            ]) + '\n')
        print(f"Summary statistics written to {args.summary_output}", file=sys.stderr)

if __name__ == '__main__':
    main()