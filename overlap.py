import ast

def interval_overlap(seg_start, seg_end, gene_start, gene_end):
    overlap_start = max(seg_start, gene_start)
    overlap_end = min(seg_end, gene_end)
    return max(0, overlap_end - overlap_start)

def load_genes(bed_file):
    genes = []
    with open(bed_file) as f:
        for line in f:
            if line.strip() == "":
                continue
            parts = line.strip().split()
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            gene_name = parts[3]
            genes.append({
                "chrom": chrom,
                "start": start,
                "end": end,
                "name": gene_name,
                "length": end - start
            })
    return genes

def load_archaic_segments(archaic_file):
    archaic_data = []
    with open(archaic_file) as f:
        for line in f:
            if line.strip() == "":
                continue
            parts = line.strip().split("\t")
            individual = parts[0]
            haplotype = parts[1]
            segments = ast.literal_eval(parts[2])
            archaic_data.append({
                "id": individual,
                "haplotype": haplotype,
                "segments": segments
            })
    return archaic_data

def compute_average_overlap(genes, archaic_data):
    results = {}
    n = len(archaic_data)
    for gene in genes:
        total_overlap = 0
        gene_start = gene["start"]
        gene_end = gene["end"]
        gene_length = gene["length"]
        for record in archaic_data:
            record_overlap = 0
            for seg in record["segments"]:
                seg_start, seg_end = seg
                record_overlap += interval_overlap(seg_start, seg_end, gene_start, gene_end)
            total_overlap += record_overlap
        avg_percentage = (total_overlap / (gene_length * n)) * 100
        results[gene["name"]] = avg_percentage
    return results

def main():
    genes_file = "genes.bed"
    archaic_file = "test.archaic.txt"

    genes = load_genes(genes_file)
    archaic_data = load_archaic_segments(archaic_file)
    avg_overlaps = compute_average_overlap(genes, archaic_data)
    
    print("Gene\tAverage_Archaic_Percentage")
    for gene_name, avg in avg_overlaps.items():
        print(f"{gene_name}\t{avg:.2f}")

if __name__ == "__main__":
    main()
