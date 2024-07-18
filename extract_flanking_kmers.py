import gzip
import argparse


def parse_gene_coordinates(filename):
    gene_data = {}

    with open(filename, 'r') as file:
        for line in file:
            parts = line.strip().split()
            gene = parts[0]
            chromosome = parts[2].split(':')[0]
            coordinates = []

            for region in parts[2:]:
                coord_range = region.split(':')[1]
                start, end = map(int, coord_range.split('-'))
                coordinates.append((start, end))

            min_coord = min(start for start, _ in coordinates)
            max_coord = max(end for _, end in coordinates)

            gene_data[gene] = (chromosome, min_coord, max_coord)

    return gene_data


def get_kmers_from_file(file_path, min_coord, max_coord, n_before, n_after):
    kmers_before = []
    kmers_after = []
    open_func = gzip.open if file_path.endswith('.gz') else open

    with open_func(file_path, 'rt') as file:
        for line in file:
            parts = line.strip().split()
            pos = int(parts[0].split(':')[1])
            kmer_data = line.strip()

            if pos < min_coord and len(kmers_before) < n_before:
                kmers_before.append((pos, kmer_data))
            elif pos > max_coord and len(kmers_after) < n_after:
                kmers_after.append((pos, kmer_data))

    return kmers_before, kmers_after


def write_kmers_to_file(kmers_before, kmers_after, output_file):
    with open(output_file, 'w') as file:
        for pos, kmer in kmers_before:
            file.write(f"{kmer}\n")
        for pos, kmer in kmers_after:
            file.write(f"{kmer}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Extract flanking k-mers for genes from reference k-mer database (Ref_kmers_chr<chr>.db.gz files)")
    parser.add_argument("input_file", help="Input file with gene coordinates")
    parser.add_argument("-r", "--ref_kmers", default="Kmer_db/Ref_kmers",
                        help="Path to reference k-mer files (default: Kmer_db/Ref_kmers)")
    parser.add_argument("-n", type=int, default=2000, help="Total number of k-mers for flanking regions")
    args = parser.parse_args()

    input_file = args.input_file
    ref_path = args.ref_kmers.rstrip("/\\") + "/"
    total_kmers = args.n
    n_before = (total_kmers // 2) + (total_kmers % 2)
    n_after = total_kmers // 2

    gene_data = parse_gene_coordinates(input_file)

    output_path = "kmer_db_files.txt"

    print(f"Extracting {total_kmers} flanking k-mers for regions in {input_file}\n")

    with open(output_path, 'w') as file:
        for gene, (chromosome, min_coord, max_coord) in gene_data.items():
            print(f"{gene}, chr {chromosome}")
            print(f"Min-max coordinates: {min_coord} - {max_coord}")
            ref_kmer_file = f"{ref_path}Ref_kmers_chr{chromosome}.db.gz"
            kmers_before, kmers_after = get_kmers_from_file(ref_kmer_file, min_coord, max_coord, n_before, n_after)

            num_kmers_before = len(kmers_before)
            num_kmers_after = len(kmers_after)

            # Adjust k-mers if there are not enough on one side
            if num_kmers_before < n_before:
                print("Not enough k-mers before the gene. Trying to get more k-mers after the gene...")
                additional_needed = n_before - num_kmers_before
                kmers_before, kmers_after = get_kmers_from_file(ref_kmer_file, min_coord, max_coord,
                                                                num_kmers_before, n_after + additional_needed)
            elif num_kmers_after < n_after:
                print("Not enough k-mers after the gene. Trying to get more k-mers before the gene...")
                additional_needed = n_after - num_kmers_after
                kmers_before, kmers_after = get_kmers_from_file(ref_kmer_file, min_coord, max_coord,
                                                                n_before + additional_needed, num_kmers_after)

            num_kmers_before = len(kmers_before)
            num_kmers_after = len(kmers_after)

            if kmers_before:
                before_min = min(kmers_before, key=lambda x: x[0])[0]
                before_max = max(kmers_before, key=lambda x: x[0])[0]
                print(f"Flanking k-mers before the gene: {before_min} - {before_max} ({num_kmers_before} kmers)")
            if kmers_after:
                after_min = min(kmers_after, key=lambda x: x[0])[0]
                after_max = max(kmers_after, key=lambda x: x[0])[0]
                print(f"Flanking k-mers after the gene: {after_min} - {after_max} ({num_kmers_after} kmers)")

            if num_kmers_before + num_kmers_after < total_kmers:
                print("WARNING: The number of total k-mers for these flanking regions is smaller than required")

            output_file = f"Kmer_db/{gene}_flanking_kmers.db"
            write_kmers_to_file(kmers_before, kmers_after, output_file)
            print(f"Flanking k-mers saved to {output_file}\n")

            kmer_path = f"Kmer_db/{gene}_kmers.db"
            flanking_path = f"Kmer_db/{gene}_flanking_kmers.db"
            file.write(f"{gene}: {kmer_path}, {flanking_path}\n")


if __name__ == "__main__":
    main()
