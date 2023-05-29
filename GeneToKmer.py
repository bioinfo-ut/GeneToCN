__author__ = "Fanny-Dhelia Pajuste, University of Tartu"
__version__ = "1.0.1"
__email__ = "fanny@ut.ee"

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse
import time
from subprocess import Popen


def parse_arguments():
    """ Defines the possible command line arguments for the program """
    parser = argparse.ArgumentParser(description='Program for creating k-mer lists for copy number estimation')
    parser.add_argument('loc_file', help='File containing the locations for the genes and flanking regions')
    parser.add_argument('reference', help='Fasta file with reference genome sequence')
    parser.add_argument('ref_list', help='List file for reference genome')
    parser.add_argument('-k', '--k', help='k-mer length', type=int, default=25, choices=range(16, 32))
    parser.add_argument('-gc', '--gc_window', help='Window size for GC content', type=int, default=250,
                        choices=range(100, 500))
    parser.add_argument('-gc_min', '--gc_min', help='Minimum GC content allowed', type=int, default=20,
                        choices=range(0, 35))
    parser.add_argument('-gc_max', '--gc_max', help='Maximum GC content allowed', type=int, default=65,
                        choices=range(40, 101))
    parser.add_argument('-mm', '--mismatches', help='Number of mismatches for glistquery', type=int, default=1,
                        choices=range(0,2))
    # parser.add_argument('-if', '--infofile', help='File for information (default stdout)')
    parser.add_argument('-gt', '--genometester', help='Location for GenomeTester4 tools')
    parser.add_argument('-o', '--output', help='Prefix for output files')
    parser.add_argument('-i', '--info', help='Show information about the process', action='store_true')
    return parser.parse_args()


def read_fasta(reference_file):
    """ Reads the reference sequence from fasta file """
    if info:
        print("Reading the reference sequence")
    fasta_sequences = list(SeqIO.parse(reference_file, 'fasta'))
    fasta = fasta_sequences[0]
    return fasta.seq


def get_regions_from_file(filename):
    """ Opens the file and adds the strings with locations of the regions to a list """
    regions = []
    with open(filename, "r") as file:
        for line in file:
            region_line = line.strip()
            regions.append(region_line)
    return regions


def get_complement(kmer):
    complement = {"A": "T", "C": "G", "T": "A", "G": "C", "a": "t", "c": "g", "t": "a", "g": "c"}
    kmer2 = ""
    for nt in kmer:
        kmer2 = complement[nt] + kmer2
    return kmer2


def get_gc_content(sequence):
    """ Calculates the GC content value for a sequence """
    c_count = sequence.count("C") + sequence.count("c")
    g_count = sequence.count("G") + sequence.count("g")
    gc = (c_count + g_count) / len(sequence) * 100
    return gc


def write_sequence(gene_file, start, end, region_id, region_desc):
    """ writes the gene sequence to a fasta file """
    with open(gene_file, "w") as w:
        out_seq = seq[start:end]
        seq_record = SeqRecord(out_seq, id=region_id, description=region_desc)
        SeqIO.write(seq_record, w, "fasta")


def filter_kmers(gene, loc_strings, glistquery_command, max_count, out_file):
    """ Runs glistquery to get unique k-mers, gets GC contents, filters the k-mers and writes to a file """
    time_1 = time.time()

    query_file = out + gene + "_mm" + str(mm) + "_max" + max_count + ".txt"
    with open(query_file, "w") as w:
        p = Popen(glistquery_command, stdout=w)

        # Meanwhile, save the locations of the k-mers with suitable GC content values (Filtering step: GC content)
        df_kmers = pd.DataFrame(columns=["loc", "nr", "kmer"])
        ch, loc_string2 = loc_strings[1].split(":")
        start, end = map(int, loc_string2.split("-"))
        for loc in range(start, end + 1 - k):
            gc_sequence = seq[loc:(loc + gc_window)]
            gc = get_gc_content(gc_sequence)
            # Filter out k-mers that are our of bounds for gc percentage value
            if gc > gc_max or gc < gc_min:
                print("Filtered out with GC value {} (limits: {}-{})".format(gc, gc_min, gc_max))
                continue
            location = ch + ":" + str(loc)
            kmer = seq[loc:(loc + k)]
            # Keep results in the right format for gmer_counter (location, nr of kmers following, kmers)
            df_kmers.loc[len(df_kmers.index)] = [location, 1, kmer]

        if info:
            print("GC contents are calculated")
            print("Kmers with GC content: {}".format(len(df_kmers.index)))
        p.wait()
    if info:
        print("Glistquery is done")
    time_2 = time.time()
    print("TIME: Running glistquery + calculating GC contents - %s" % convert_time(time_2 - time_1))

    # When glistquery has finished, keep those that are in the output (Filtering step: uniqueness, 1 mismatch)
    #removed = []
    with open(query_file, "r") as r:
        query_kmers = set()
        for line in r:
            kmer = line.strip().split("\t")[0]
            if kmer not in df_kmers["kmer"].values:
                kmer = get_complement(kmer)
            query_kmers.add(kmer)
    df_kmers_filtered = df_kmers[df_kmers["kmer"].isin(query_kmers)]
    if info:
        print("Kmers from query: {}".format(len(query_kmers)))
        print("Filtered k-mers: {}".format(len(df_kmers_filtered)))

    df_kmers_filtered.to_csv(out_file, index=False, sep="\t", header=False)

    time_6 = time.time()
    print("TIME: Filtering and writing to file - %.2f" % (time_6 - time_5))


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


# noinspection SpellCheckingInspection
def main():
    if info:
        print("Reading the location file")

    regions = get_regions_from_file(region_file)

    for region in regions:
        time_region_start = time.time()
        loc_strings = region.split(" ")
        region_type = loc_strings.pop(1)
        gene = loc_strings[0]
        out_file = out + gene + "_kmers.txt"
        max_count = "1"
        if info:
            print("Finding unique k-mers for " + gene)

        if len(loc_strings) > 2 and region_type == "G":
            max_count = str(len(loc_strings) - 1)
            commands = []
            compare_files = []
            for i in range(1, len(loc_strings)):
                if info:
                    print("Region " + loc_strings[i])
                ch, loc_string2 = loc_strings[i].split(":")
                start, end = map(int, loc_string2.split("-"))
                gene_file = out + gene + "_sequence_" + str(i) + ".fa"

                # write the gene sequence to a file
                write_sequence(gene_file, start, end, loc_strings[i], gene + " copy " + str(i))

                compare_file = out + gene + "_initial_kmers_" + str(i)
                # Run glistmaker (needed for getting the intersect of k-mers with glistcompare)
                glistmaker_command = [gt_loc + "glistmaker", gene_file, "-w", str(k), "-o", compare_file]
                compare_files.append(compare_file + "_" + str(k) + ".list")
                commands.append(glistmaker_command)

            time_2 = time.time()
            # Run glistmakers in parallel
            procs = [Popen(i) for i in commands]
            for p in procs:
                p.wait()
            time_3 = time.time()
            print("TIME: Running glistmakers - %s" % convert_time(time_3 - time_2))
            i = 1
            while len(compare_files) > 1:
                compare_files_2 = []
                glistcompare_commands = []
                while len(compare_files) > 1:
                    intersect_file = out + gene + "_" + str(i)
                    i += 1
                    glistcompare_command = [gt_loc + "glistcompare", compare_files.pop(0), compare_files.pop(0), "-i",
                                            "-o",
                                            intersect_file]
                    glistcompare_commands.append(glistcompare_command)
                    compare_files_2.append(intersect_file + "_" + str(k) + "_intrsec.list")
                compare_files = compare_files_2[:] + compare_files[:]
                procs = [Popen(i) for i in glistcompare_commands]
                for p in procs:
                    p.wait()

            time_4 = time.time()
            print("TIME: Running glistcompare(s) - %s" % convert_time(time_4 - time_3))

            list_file = compare_files[0]

            # Filter k-mers by gc_count and uniqueness
            glistquery_command = [gt_loc + "glistquery", ref_list, "-l", list_file, "-mm", str(mm), "-max", max_count]
            filter_kmers(gene, loc_strings, glistquery_command, max_count, out_file)

        else:
            gene_file = out + gene + "_sequence.fa"
            if len(loc_strings) > 2 and region_type == "F":
                with open(gene_file, "w") as w:
                    for i in range(1, len(loc_strings)):
                        if info:
                            print("Region " + loc_strings[i])
                        ch, loc_string2 = loc_strings[i].split(":")
                        start, end = map(int, loc_string2.split("-"))
                        out_seq = seq[start:end]
                        seq_record = SeqRecord(out_seq, id=loc_strings[i], description=gene + " part " + str(i))
                        SeqIO.write(seq_record, w, "fasta")
            else:
                # If there is only one gene copy in reference, write the sequence to a file and
                # get the unique k-mers using glistquery
                if info:
                    print("Region " + loc_strings[1])

                ch, loc_string2 = loc_strings[1].split(":")
                start, end = map(int, loc_string2.split("-"))
                gene_file = out + gene + "_sequence.fa"

                # write the gene sequence to a file
                write_sequence(gene_file, start, end, loc_strings[1], gene)

            # Filter k-mers by gc_count and uniqueness
            glistquery_command = [gt_loc + "glistquery", ref_list, "-s", gene_file, "-mm", str(mm), "-max", max_count]
            filter_kmers(gene, loc_strings, glistquery_command, max_count, out_file)

        time_region_end = time.time()
        print("TIME: %s - %s" % (gene, convert_time(time_region_end - time_region_start)))

    time_finish = time.time()
    print("TIME: Total - %s" % convert_time(time_finish - time_start))


time_start = time.time()
# time_string = now.strftime('%Y%m%d_%H%M%S')

# Parse arguments
args = parse_arguments()
info = args.info
reference_file = args.reference
ref_list = args.ref_list
gc_window = args.gc_window
gc_min = args.gc_min
gc_max = args.gc_max
k = args.k
out = args.output
region_file = args.loc_file
gt_loc = args.genometester
mm = args.mismatches

seq = read_fasta(reference_file)

if __name__ == "__main__":
    main()
