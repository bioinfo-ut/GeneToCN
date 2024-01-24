import argparse
import time
from glob import glob
from subprocess import Popen
from statistics import median, quantiles
import matplotlib.pyplot as plt
from os import mkdir, path
from contextlib import suppress


def parse_arguments():
    """ Defines the possible command line arguments for the program """
    parser = argparse.ArgumentParser(description='Program for k-mer based copy number estimation')
    # Options for k-mer db files
    # Separate files for each gene and one flanking region - .db files
    parser.add_argument('-g', '--gene_files', help='File(s) containing gene k-mers', nargs='+')
    parser.add_argument('-f', '--flanking_file', help='File containing k-mers of the single-copy reference region')
    # One big file containing all genes multiple flanking regions - .db or list of k-mers
    parser.add_argument('-db', '--kmer_db', help='File containing k-mer database (multiple genes and flanking regions')
    parser.add_argument('-kl', '--kmer_list',
                        help='File containing a list of k-mers (multiple genes and flanking regions')

    # Options for sequencing data of idividual - .fastq or .list
    parser.add_argument('-s', '--seq_files', help='Fastq file(s) containing sequencing data', nargs='+')
    parser.add_argument('-l', '--list_file', help='GenomeTester4 list file from fastq data')

    # Needed if one big k-mer file is used (to know which k-mers belong to which genes
    parser.add_argument('-kp', '--kmer_db_paths',
                        help='File containing all k-mers (multiple genes and flanking regions')

    # Used for drawing plots for gene region, if not provided, the first and last location of the k-mers are used
    parser.add_argument('-r', '--region_file', help='File with gene coordinates (for plot)')

    parser.add_argument('-p', '--max_proc', help='Maximum number of gmer_counters allowed to run in parallel', type=int,
                        default=5, choices=range(1, 20))
    parser.add_argument('-t', '--test', help='Testing (run without calling gmer_counter)', action='store_true')
    parser.add_argument('-gm', '--gmercounter',
                        help='Location for gmercounter tool (or glistquery if a list is provided)')
    parser.add_argument('-o', '--output', help='Prefix for output files')
    parser.add_argument('-d', '--out_dir', help='Directory for output files', default='./')
    parser.add_argument('-i', '--info', help='Show information about the process', action='store_true')
    parser.add_argument('-kcn', '--kmer_cn', help='Only print out copy numbers separately for each k-mer',
                        action='store_true')

    # Options that would be nice to add, but not possible to use yet
    # parser.add_argument('-e', '--exon_file', help='File with exon coordinates (for plot)')
    # parser.add_argument('-b', '--binary', help='The input k-mer files are in binary format', action='store_true')
    # parser.add_argument('-if', '--infofile', help='File for information (default stdout)')
    return parser


def parse_filenames(arg_list):
    """ Parses the filenames for wildcards """
    file_list = []
    for name in arg_list:
        file_list += glob(name)
    return file_list[:]


def get_complement(kmer):
    complement = {"A": "T", "C": "G", "T": "A", "G": "C", "a": "t", "c": "g", "t": "a", "g": "c"}
    kmer2 = ""
    for nt in kmer:
        kmer2 = complement[nt] + kmer2
    return kmer2


def convert_time(t):
    """ converts time from seconds (float) to hh:mm:ss (string) """
    t = int(round(t))
    hours = t // 3600
    if hours < 10:
        hh = "0" + str(hours)
    else:
        hh = str(hours)
    min = (t % 3600) // 60
    if min < 10:
        mm = "0" + str(min)
    else:
        mm = str(min)
    sec = ((t % 3600) % 60)
    if sec < 10:
        ss = "0" + str(sec)
    else:
        ss = str(sec)
    ans = hh + ":" + mm + ":" + ss
    return ans


def draw_plot(gene, locations, values, gene_median, is_flanking):
    """ creates a plot (.png) with k-mer copy values for gene region """
    min_loc = min(locations)
    max_loc = max(locations)
    if regions and gene in regions.keys():
        min_loc, max_loc = regions[gene]
    plt.plot([min_loc, max_loc], [gene_median, gene_median], color="#639e93")
    plt.plot(locations, values, "o", color="#5377a6", ms=2)
    plt.grid(axis="y")
    plt.title(gene)
    if is_flanking:
        plt.xlabel("k-mer nr")
        plt.ylabel("k-mer frequency")
    else:
        plt.xlabel("position")
        plt.ylabel("k-mer copy number")
    plt.ticklabel_format(axis='x', style='plain')
    plt.savefig(out_dir + "/plots/" + gene + "_" + out_prefix + "plot.png")
    plt.clf()
    return


def get_kmer_counts(gene):
    """ get k-mer locations and counts from count file """
    with open(out_dir + "counts/" + gene + "_" + out_prefix + "counts.txt", "r") as gr:
        values = []
        locations = []
        for line in gr.readlines()[2:]:
            parts = line.strip().split("\t")
            values += [int(x) for x in parts[2:]]
            if len(parts) == 3:
                locations.append(int(parts[0].split(":")[1]))
    return locations, values


def get_median_and_CN(gene, locations, values, flanking_median=-1):
    """ counts the median and copy number (if not flanking region) based on k-mer counts """
    gene_median = median(values)
    q_start = quantiles(values, n=100)[4]
    q_end = quantiles(values, n=100)[94]
    q_str = ", 95% percentile: (" + str(q_start) + "," + str(q_end) + ")"
    if info:
        print("Median for " + gene + " " + str(gene_median) + q_str)

    # if the region is a gene region
    if flanking_median >= 0:
        try:
            kmer_cns = [x / flanking_median * 2 for x in values]
            cn = round(gene_median / flanking_median * 2, 2)
        except ZeroDivisionError:
            print("Median k-mer frequency in flanking region is 0 - cannot estimate copy number (ZeroDivisionError).")
            return -1
        if locations:
            if print_kmer_cn:
                kmer_cn_file = out_dir + "kmer_cn/" + gene + "_" + out_prefix + "kmer_cn.txt"
                with open(kmer_cn_file, "w") as kmer_cn_w:
                    for i in range(len(locations)):
                        kmer_cn_w.write(str(locations[i]) + "\t" + str(kmer_cns[i]) + "\n")
            else:
                draw_plot(gene, locations, kmer_cns, cn, False)
        # ci_str = ", 95% CI: (" + str(ci_start) + "," + str(ci_end) + ")"
        if info:
            print("Copy number for " + gene + ": " + str(cn))
        out_str = gene + "\t" + str(cn) + "\t" + str(gene_median) + "\t" + str(len(values)) + "\n"
        return cn, out_str

    # If the region is a flanking region
    out_str = str(gene_median) + "\n# Nr of reference k-mers: " + str(len(values)) \
              + "\nGene\tCN\tMedian_freq\tN_kmers\n"
    if locations and not print_kmer_cn:
        draw_plot(gene, locations, values, gene_median, True)
    return gene_median, out_str


def get_regions_from_file(region_file):
    """ reads the region file and saves the start and end position (of the first region) for every gene/flanking region"""
    regions = dict()
    if region_file:
        with open(region_file, "r") as file:
            for line in file:
                region = line.strip()
                loc_strings = region.split(" ")
                gene = loc_strings[0]
                gene_type = loc_strings[1]
                locs = []
                for loc_string in loc_strings[2:]:
                    parts = loc_string.split(":")
                    positions = [int(x) for x in parts[1].split("-")]
                    locs += positions
                    if gene_type == "G":
                        break
                regions[gene] = [min(positions), max(positions)]
    return regions


def run_separate_counters(flanking_prefix):
    """ runs gmer_counter for separate k-mer databases for each gene/flanking region"""
    count_files = []
    procs = []
    errors = open(out_dir + "counter_errors/" + out_prefix + "counter_errors.txt", "w")

    # Flanking region
    if info:
        print("Starting gmercounter for " + flanking_prefix)
    flanking_count_file = out_dir + "counts/" + flanking_prefix + "_" + out_prefix + "counts.txt"
    gmercounter_command = [gm_loc + "gmer_counter", "-db", flanking_file] + seq_files
    w = open(flanking_count_file, "w")
    count_files.append(w)
    procs.append(Popen(gmercounter_command, stdout=w, stderr=errors))

    # Gene regions
    genes = []
    for filename in gene_files:
        gene = filename.split("/")[-1].split("_")[0]
        genes.append(gene)
        if info:
            print("Starting gmercounter for " + gene)
        gmercounter_command = [gm_loc + "gmer_counter", "-db", filename] + seq_files
        w = open(out_dir + "counts/" + gene + "_" + out_prefix + "counts.txt", "w")
        count_files.append(w)
        procs.append(Popen(gmercounter_command, stdout=w, stderr=errors))

        # If there are too many processes running, wait for one of them to finish before starting a new one
        if len(procs) > max_proc and info:
            print("More than " + str(max_proc) + " processes running, waiting for one of them to finish")
        while len(procs) > max_proc:
            for i in range(len(procs)):
                if procs[i].poll() is not None:
                    count_files.pop(i).close()
                    procs.pop(i)
                    break
            # Wait before checking again
            time.sleep(10)

    for i in range(len(procs)):
        procs[i].wait()
        count_files[i].close()

    errors.close()

    time_2 = time.time()
    print("TIME: Counting k-mer frequencies - %s" % convert_time(time_2 - time_start))

    return genes


def run_counter(count_file):
    """ Runs gmer_counter for one big k-mer database with k-mers for all af the genes/flanking regions"""
    errors = open(out_dir + "counter_errors/" + out_prefix + "counter_errors.txt", "w")
    if info:
        print("Starting gmercounter")
    gmercounter_command = [gm_loc + "gmer_counter", "-db", kmer_db] + seq_files
    w = open(count_file, "w")
    p = Popen(gmercounter_command, stdout=w, stderr=errors)
    p.wait()

    w.close()
    errors.close()

    time_2 = time.time()
    print("TIME: Counting k-mer frequencies - %s" % convert_time(time_2 - time_start))
    return


def run_glistquery(query_file):
    if info:
        print("Starting glistquery")
    glistquery_command = [gm_loc + "glistquery", list_file, "-f", kmer_list]
    with open(query_file, "w") as w, open(out_dir + "counter_errors/" + out_prefix + "query_errors.txt", "w") as errors:
        p = Popen(glistquery_command, stdout=w, stderr=errors)
        p.wait()

    time_2 = time.time()
    print("TIME: Counting k-mer frequencies - %s" % convert_time(time_2 - time_start))
    return


def one_flanking_region():
    """ KmerToCN for one set of genes and one flanking region (with separate db files)"""
    flanking_prefix = path.basename(flanking_file).split(".")[0]
    genes = run_separate_counters(flanking_prefix)

    locations, values = get_kmer_counts(flanking_prefix)
    flanking_median, out_str = get_median_and_CN(flanking_prefix, locations, values)
    cn_out_str = "# Median frequency of reference k-mers: " + out_str

    cn = dict()
    for gene in genes:
        locations, values = get_kmer_counts(gene)
        gene_cn, out_str = get_median_and_CN(gene, locations, values, flanking_median)
        cn_out_str += out_str
        cn[gene] = gene_cn

    if not print_kmer_cn:
        with open(cn_result_file, "w") as w_cn:
            w_cn.write(cn_out_str)
    return


def many_regions():
    """ KmerToCN for multiple genes and flanking regions at the same time (one big db file)"""
    if kmer_list is None:
        file_prefix = path.basename(kmer_db).split(".")[0]
    else:
        file_prefix = path.basename(kmer_list).split(".")[0]
    count_file = out_dir + "counts/" + file_prefix + "_" + out_prefix + "counts.txt"
    if not test:
        if kmer_list is None:
            run_counter(count_file)
        else:
            run_glistquery(count_file)
    all_counts = dict()
    # Read all the counts to a dictionary
    with open(count_file, "r") as r:
        comment_lines = 0
        if kmer_list is None:
            comment_lines = 2
        for line in r.readlines()[comment_lines:]:
            parts = line.strip().split("\t")
            if kmer_list is None:
                counts = [int(x) for x in parts[2:]]
            else:
                counts = int(parts[1])
            all_counts[parts[0]] = counts

    cn_out_str = ""
    with open(kmer_paths, "r") as r:
        is_flanking = True
        for line in r:
            parts = line.strip().split(" ")
            if not line.strip():
                cn_out_str += "\n"
                is_flanking = True
                continue
            gene = parts[0]
            filename = parts[1]
            values = []
            locations = []
            with open(filename, "r") as r2:
                for line2 in r2:
                    parts2 = line2.strip().split("\t")
                    counts = []
                    if kmer_list is None:
                        try:
                            counts = all_counts[parts2[0]]
                        except KeyError:
                            print("\"" + parts2[0] + "\" not found in count file (KeyError)")
                            # print(len(all_counts.keys()))
                            # quit()
                    else:
                        try:
                            for kmer in parts2[2:]:
                                if kmer in all_counts:
                                    counts.append(all_counts[kmer])
                                else:
                                    counts.append(all_counts[get_complement(kmer)])
                        except KeyError:
                            print("\"" + kmer + "\" not found in count file (KeyError)")

                    values += counts
                    if len(counts) == 1:
                        locations.append(int(parts2[0].split(":")[1]))
            if is_flanking:
                flanking_median, out_str = get_median_and_CN(gene, locations, values)
                cn_out_str += "# Median frequency of reference k-mers: " + out_str
                is_flanking = False
            else:
                print(len(locations), len(values))
                print(locations)
                print(values)
                gene_cn, out_str = get_median_and_CN(gene, locations, values, flanking_median)
                cn_out_str += out_str
        if not print_kmer_cn:
            with open(cn_result_file, "w") as w_cn:
                w_cn.write(cn_out_str)
    return


def main():
    if kmer_db or kmer_list:
        many_regions()
    else:
        one_flanking_region()

    time_finish = time.time()
    print("TIME: Total - %s" % convert_time(time_finish - time_start))


time_start = time.time()
# time_string = now.strftime('%Y%m%d_%H%M%S')

# Parse arguments

parser = parse_arguments()
args = parser.parse_args()
kmer_db = args.kmer_db
flanking_file = args.flanking_file
kmer_paths = args.kmer_db_paths
list_file = args.list_file
kmer_list = args.kmer_list
seq_files_arg = args.seq_files
gene_files_arg = args.gene_files
test = args.test
print_kmer_cn = args.kmer_cn

option1 = seq_files_arg is not None and gene_files_arg is not None and flanking_file is not None
option2 = seq_files_arg is not None and kmer_db is not None and kmer_paths is not None
option3 = list_file is not None and kmer_list is not None and kmer_paths is not None

if sum([option1, option2, option3]) != 1 and (test == False):
    print("Possible combinations of options:\n1) -g -f -s\n2) -db -kp -s\n3) -kl -kp -l")
    parser.print_usage()
    print("Use -h for more information")
    quit()

if (option1 or option2) and (test == False):
    seq_files = parse_filenames(seq_files_arg)
    if option1:
        gene_files = parse_filenames(gene_files_arg)

gm_loc = args.gmercounter
out_dir = args.out_dir.rstrip("/\\") + "/"
with suppress(FileExistsError):
    mkdir(out_dir + "plots/")
with suppress(FileExistsError):
    mkdir(out_dir + "counts/")
with suppress(FileExistsError):
    mkdir(out_dir + "counter_errors/")
if print_kmer_cn:
    with suppress(FileExistsError):
        mkdir(out_dir + "kmer_cn/")
out_prefix = args.output.strip("_") + "_"
info = args.info
max_proc = args.max_proc
region_file = args.region_file

regions = get_regions_from_file(region_file)

cn_result_file = out_dir + out_prefix + "cn.txt"

print(out_prefix.strip("_"))

if __name__ == "__main__":
    main()
