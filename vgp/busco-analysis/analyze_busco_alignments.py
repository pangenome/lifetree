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
    except subprocess.CalledProcessError as e:
        # Report the error and exit
        print(f"\nERROR: impg query failed for {paf_file}", file=sys.stderr)
        print(f"Command: {' '.join(e.cmd)}", file=sys.stderr)
        print(f"Return code: {e.returncode}", file=sys.stderr)
        if e.stderr:
            print(f"Error message:\n{e.stderr}", file=sys.stderr)
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

def calculate_covered_bases(alignments, gene_start, gene_end, coord_field='target'):
    """Calculate how many bases of the gene are covered by alignments"""
    if not alignments:
        return 0
    
    # Extract and clip intervals to gene boundaries
    intervals = []
    for aln in alignments:
        start = max(aln[f'{coord_field}_start'], gene_start)
        end = min(aln[f'{coord_field}_end'], gene_end)
        if start < end:
            intervals.append((start, end))
    
    if not intervals:
        return 0
    
    # Sort intervals by start position
    intervals.sort()
    
    # Merge overlapping intervals
    merged = [intervals[0]]
    for start, end in intervals[1:]:
        if start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    
    # Calculate total covered bases
    total_covered = sum(end - start for start, end in merged)
    return total_covered

def count_covering_alignments(alignments, gene_start, gene_end, coord_field='target'):
    """Count number of non-overlapping alignments needed to cover the gene"""
    if not alignments:
        return 0
    
    # Extract intervals that overlap with gene
    intervals = []
    for aln in alignments:
        start = max(aln[f'{coord_field}_start'], gene_start)
        end = min(aln[f'{coord_field}_end'], gene_end)
        if start < end:
            intervals.append((start, end))
    
    if not intervals:
        return 0
    
    # Sort intervals by start position
    intervals.sort()
    
    # Count non-overlapping intervals (fragmentation)
    count = 1
    last_end = intervals[0][1]
    
    for start, end in intervals[1:]:
        if start > last_end:  # Gap between alignments
            count += 1
        last_end = max(last_end, end)
    
    return count

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
        'target_covered_bases': 0,
        'query_covered_bases': 0,
        'n_covering_alignments_target': 0,
        'n_covering_alignments_query': 0,
        'has_other_alignments': False
    }
    
    # Add query length if gene is present in query
    if gene_id in query_busco:
        query_chrom, query_start, query_end, _ = query_busco[gene_id]
        result['query_length'] = query_end - query_start
    else:
        result['query_length'] = 0
    
    # Calculate target coverage from all alignments
    if alignments:
        result['target_covered_bases'] = calculate_covered_bases(alignments, target_start, target_end, 'target')
        result['n_covering_alignments_target'] = count_covering_alignments(alignments, target_start, target_end, 'target')
        
        # Check alignments overlapping with BUSCO location and calculate query coverage
        if gene_id in query_busco:
            query_chrom, query_start, query_end, _ = query_busco[gene_id]
            busco_alignments = []
            non_busco_alignments = []
            
            for aln in alignments:
                if (aln['query_seq'] == query_chrom and 
                    not (aln['query_end'] <= query_start or aln['query_start'] >= query_end)):
                    result['query_n_busco_alignments'] += 1
                    busco_alignments.append(aln)
                else:
                    result['query_n_non_busco_alignments'] += 1
                    non_busco_alignments.append(aln)
            
            # Calculate coverage of query BUSCO region
            result['query_covered_bases'] = calculate_covered_bases(busco_alignments, query_start, query_end, 'query')
            result['n_covering_alignments_query'] = count_covering_alignments(busco_alignments, query_start, query_end, 'query')
            result['has_other_alignments'] = len(non_busco_alignments) > 0
        else:
            # Gene not in query BUSCO annotations, all alignments are to "other" regions
            result['query_n_non_busco_alignments'] = len(alignments)
            result['has_other_alignments'] = len(alignments) > 0
    
    return result

def classify_gene(result):
    """Classify gene based on coverage and alignment patterns"""
    # Calculate coverage fractions
    target_cov = result['target_covered_bases'] / result['target_length'] if result['target_length'] > 0 else 0.0
    query_cov = result['query_covered_bases'] / result['query_length'] if result['query_length'] > 0 else 0.0
    
    has_others = result['has_other_alignments']
    
    # Determine if fragmented (need multiple alignments for coverage)
    is_fragmented = (result['n_covering_alignments_target'] > 1 or 
                    result['n_covering_alignments_query'] > 1)
    
    # No alignments at all
    if result['query_n_total_alignments'] == 0:
        return 'No-alignments'
    
    # Check if either coverage is 0
    if target_cov == 0 or query_cov == 0:
        # Missing category
        if has_others or result['query_n_total_alignments'] > 0:
            return 'Missing-others'
        else:
            return 'No-alignments'
    
    # Both coverages are > 0, classify based on thresholds
    # According to the requirements, BOTH must meet the threshold for each category
    if target_cov >= 0.9 and query_cov >= 0.9:
        coverage_level = 'good'
    elif target_cov >= 0.4 and query_cov >= 0.4:
        coverage_level = 'medium'
    elif target_cov > 0 and query_cov > 0:
        coverage_level = 'bad'
    else:
        # This shouldn't happen as we checked for 0 above
        return 'Missing-others' if has_others else 'No-alignments'
    
    # Build classification string
    classification = coverage_level
    if is_fragmented:
        classification += '-fragmented'
    if has_others:
        classification += '-others'
    
    return classification

def main():
    parser = argparse.ArgumentParser(description='Analyze BUSCO gene alignments between genomes')
    parser.add_argument('--target-bed', required=True, help='BUSCO BED file for target genome')
    parser.add_argument('--query-bed', required=True, help='BUSCO BED file for query genome')
    parser.add_argument('--paf', required=True, help='PAF alignment file')
    parser.add_argument('--output', help='Output file for per-gene results (default: stdout)')
    parser.add_argument('--summary-output', help='Output file for summary statistics')
    parser.add_argument('--coverage-output', help='Output file for coverage statistics')
    parser.add_argument('--gene', help='Analyze only this specific BUSCO gene')
    parser.add_argument('--threads', type=int, default=4, help='Number of threads (default: 4)')
    
    args = parser.parse_args()
    
    # Parse BUSCO annotations
    print("Loading BUSCO annotations...", file=sys.stderr)
    target_busco = parse_bed_file(args.target_bed)
    query_busco = parse_bed_file(args.query_bed)
    
    print(f"Found {len(target_busco)} BUSCO genes in target", file=sys.stderr)
    print(f"Found {len(query_busco)} BUSCO genes in query", file=sys.stderr)
    
    # Determine genes to analyze (all target genes, not just common ones)
    genes_to_analyze = [args.gene] if args.gene else sorted(target_busco.keys())

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
        'gene_present_in_query', 'query_n_total_alignments',
        'query_n_busco_alignments', 'query_n_non_busco_alignments', 
        'target_busco_completeness', 'query_busco_completeness',
        'n_covering_alignments_target', 'n_covering_alignments_query',
        'has_other_alignments', 'classification'
    ]) + '\n')
    
    # Write results and collect stats
    stats = defaultdict(int)
    stats_query_only = defaultdict(int)  # Stats for genes present in query
    
    # Coverage statistics
    total_target_length = 0
    total_query_length = 0
    total_target_covered_bases = 0
    total_query_covered_bases = 0
    
    for r in results:
        classification = classify_gene(r)
        r['classification'] = classification
        
        # Collect stats for all genes
        stats[classification] += 1
        
        # Collect stats only for genes present in query
        if r['gene_present_in_query']:
            stats_query_only[classification] += 1
        
        # Collect coverage statistics
        total_target_length += r['target_length']
        total_target_covered_bases += r['target_covered_bases']
        
        if r['gene_present_in_query']:
            total_query_length += r['query_length']
            total_query_covered_bases += r['query_covered_bases']
        
        # Calculate coverage fractions for output
        target_completeness = r['target_covered_bases'] / r['target_length'] if r['target_length'] > 0 else 0.0
        query_completeness = r['query_covered_bases'] / r['query_length'] if r['query_length'] > 0 else 0.0
        
        output.write('\t'.join([
            r['gene_id'], r['target_name'], str(r['target_start']), str(r['target_end']),
            str(r['target_length']),
            'Yes' if r['gene_present_in_query'] else 'No',
            str(r['query_n_total_alignments']), str(r['query_n_busco_alignments']), 
            str(r['query_n_non_busco_alignments']),
            f"{target_completeness:.4f}",
            f"{query_completeness:.4f}",
            str(r['n_covering_alignments_target']),
            str(r['n_covering_alignments_query']),
            'Yes' if r['has_other_alignments'] else 'No',
            classification
        ]) + '\n')
    
    if args.output:
        output.close()
    
    # Write summary statistics to file if requested
    if args.summary_output:
        with open(args.summary_output, 'w') as summary_file:
            # Define all categories in order
            all_categories = ['good', 'good-fragmented', 'good-others', 'good-fragmented-others',
                            'medium', 'medium-fragmented', 'medium-others', 'medium-fragmented-others',
                            'bad', 'bad-fragmented', 'bad-others', 'bad-fragmented-others',
                            'Missing-others', 'No-alignments']
            
            # Write header with all categories plus total
            header_columns = all_categories + ['total']
            summary_file.write('\t'.join(header_columns) + '\n')
            
            # Write data row with counts for each category (only genes present in query)
            n_genes_in_query = sum(1 for r in results if r['gene_present_in_query'])
            data_row = [str(stats_query_only[category]) for category in all_categories]
            data_row.append(str(n_genes_in_query))  # Add total count of genes in query
            summary_file.write('\t'.join(data_row) + '\n')

        print(f"Summary statistics written to {args.summary_output}, (include only {n_genes_in_query} genes present in query out of {len(results)} total)", file=sys.stderr)

    # Write coverage statistics to file if requested
    if args.coverage_output:
        with open(args.coverage_output, 'w') as coverage_file:
            # Write header
            coverage_file.write('\t'.join([
                'total_target_busco_length',
                'total_query_busco_length', 
                'target_covered_bases',
                'query_covered_bases',
                'target_coverage_fraction',
                'query_coverage_fraction'
            ]) + '\n')
            
            # Calculate coverage fractions
            target_coverage_fraction = total_target_covered_bases / total_target_length if total_target_length > 0 else 0.0
            query_coverage_fraction = total_query_covered_bases / total_query_length if total_query_length > 0 else 0.0
            
            # Write data row
            coverage_file.write('\t'.join([
                str(total_target_length),
                str(total_query_length),
                str(total_target_covered_bases),
                str(total_query_covered_bases),
                f"{target_coverage_fraction:.6f}",
                f"{query_coverage_fraction:.6f}"
            ]) + '\n')
        
        print(f"Coverage statistics written to {args.coverage_output}", file=sys.stderr)
    
    # # Print summary to stderr
    # print("\nSummary:", file=sys.stderr)
    # for category in sorted(stats.keys()):
    #     print(f"  {category}: {stats[category]}", file=sys.stderr)
    # print(f"  Total: {len(results)}", file=sys.stderr)

if __name__ == '__main__':
    main()