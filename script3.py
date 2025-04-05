#!/usr/bin/env python3
import argparse
import sys
import numpy

def read_indiv(indiv_file):
    indivs = []
    with open(indiv_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if parts:
                indivs.append((parts[0], parts[1], parts[2]))
    return indivs

def read_snp(snp_file):
    snps = []
    with open(snp_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if parts:
                snp_name = parts[0]
                chrom = parts[1]
                genetic_pos = parts[2]
                physical_pos = parts[3]
                ref = parts[4] if len(parts) >= 5 else "N"
                alt = parts[5] if len(parts) >= 6 else "N"
                snps.append((snp_name, chrom, genetic_pos, physical_pos, ref, alt))
    return snps

def read_bed_file(bed_file):
    regions = []
    with open(bed_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 3:
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                regions.append((chrom, start, end))
    return regions

def snp_in_regions(snp_chrom, snp_pos, regions):
    for region_chrom, start, end in regions:
        if snp_chrom == region_chrom and start <= snp_pos <= end:
            return True
    return False

def read_info_file(info_file):
    info_dict = {}
    with open(info_file, 'r') as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 6:
                try:
                    pos = int(parts[0])
                except ValueError:
                    continue
                ancestral = parts[3]
                outgroup = parts[4]
                archaic = parts[5]
                info_dict[pos] = (ancestral, outgroup, archaic)
    return info_dict

def process_geno_line(geno_line, indices):
    geno = list(geno_line.strip())
    return [geno[i] for i in indices]

def convert_genotype_to_haplotypes(g):
    if g == '0':
        return ['0', '0']
    elif g == '2':
        return ['1', '1']
    elif g == '1':
        return ['0', '1']
    elif g == '9':
        return ['9', '9']
    else:
        return ['?', '?']

def compute_missing_rates(geno_file, snps, target_indices, target_chrom):
    missing_counts = {i: 0 for i in target_indices}
    total_snps = 0

    with open(geno_file, 'r') as gf:
        for snp_index, geno_line in enumerate(gf):
            if snp_index >= len(snps):
                break
            _, chrom, _, _, _, _ = snps[snp_index]
            if chrom != target_chrom:
                continue
            total_snps += 1
            geno_values = process_geno_line(geno_line, target_indices)
            for i, g in zip(target_indices, geno_values):
                if g == '9':
                    missing_counts[i] += 1
    missing_rates = {i: (missing_counts[i] / total_snps if total_snps > 0 else 1.0)
                     for i in target_indices}
    return missing_rates, total_snps

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--geno', required=True, help=".eigenstratgeno")
    parser.add_argument('--snp', required=True, help=".snp")
    parser.add_argument('--indiv', required=True, help=".ind")
    parser.add_argument('--bed', required=True, help=".bed")
    parser.add_argument('--info', required=True, help="for additional snp info")
    parser.add_argument('--pop', required=True, help="population")
    parser.add_argument('--chrom', required=True, help="chromosome")
    parser.add_argument('--out', required=True, help="output")
    args = parser.parse_args()

    regions = read_bed_file(args.bed)
    info_dict = read_info_file(args.info)

    indivs = read_indiv(args.indiv)
    target_indices_all = [i for i, rec in enumerate(indivs) if rec[2] == args.pop]
    target_indices = target_indices_all[:10]

    if not target_indices:
        print("No individuals found for population:", args.pop)
        sys.exit(1)

    snps = read_snp(args.snp)

    missing_rates, snps_considered = compute_missing_rates(args.geno, snps, target_indices, args.chrom)
    print(f"Total SNPs on chromosome {args.chrom} considered: {snps_considered}")

    filtered_indices = target_indices
    print(f"Filtered individuals count: {len(filtered_indices)} (from {len(target_indices)})")
    with open(args.geno, 'r') as geno_f, open(args.out, 'w') as out_f:
        out_f.write("#POSITIONS\t#REF\t#ALT\tANCESTRAL\t#OUTGROUP\t#ARCHAIC\t#OBSERVATIONS\n")
        snp_index = 0
        for geno_line in geno_f:
            if snp_index >= len(snps):
                break
            snp_info = snps[snp_index]
            snp_index += 1
            snp_name, chrom, _, physical_pos, ref, alt = snp_info

            if chrom != args.chrom:
                continue

            try:
                physical_pos_int = int(physical_pos)
            except ValueError:
                continue
            if not snp_in_regions(chrom, physical_pos_int, regions):
                continue

            geno_values = process_geno_line(geno_line, filtered_indices)
            total = len(geno_values)
            missing_count = geno_values.count('9')
            missing_fraction = missing_count / total if total > 0 else 0

            if missing_fraction > 0:
                if missing_fraction < 0.85:
                    counts = {}
                    for g in geno_values:
                        if g != '9':
                            counts[g] = counts.get(g, 0) + 1
                    if counts:
                        majority_call = max(counts, key=counts.get)
                    else:
                        majority_call = '0'
                    geno_values = [g if g != '9' else majority_call for g in geno_values]
                else:
                    geno_values = [g if g != '9' else '0' for g in geno_values]
            haplotypes = []
            for g in geno_values:
                haplotypes.extend(convert_genotype_to_haplotypes(g))

            if physical_pos_int in info_dict:
                ancestral, outgroup, archaic = info_dict[physical_pos_int]
            else:
                ancestral, outgroup, archaic = '.', '0', '.'

            observations = " ".join(haplotypes)
            out_line = f"{physical_pos}\t{ref}\t{alt}\t{ancestral}\t{outgroup}\t{archaic}\t{observations}\n"
            out_f.write(out_line)

if __name__ == '__main__':
    main()
