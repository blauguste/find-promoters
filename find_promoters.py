from Bio import SeqIO, SeqFeature
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import nt_search
import csv

gb_path = 'NC_004350_2.gb'
ref_record = SeqIO.read(open(gb_path), 'genbank')
ref_sequence = ref_record.seq

# List of sRNA TSSs. Csv with [0] uid, [1] threshold, [2] genome, [3] feature name, [4] strand, [5] TSS
TSS_list = []
with open('s_mutans_sample_pos.csv', 'r') as TSS_input:
    TSS_reader = csv.reader(TSS_input, delimiter = ',')
    next(TSS_reader, None) # ignore the header row
    for row in TSS_reader:
        TSS_list.append((row[0], row[1], row[2], row[3], row[4], int(row[5])))
print(TSS_list)

# Promoter sequences to look for: A2 and T6 always conserved plus 2 out of 4 variable sites (total of 4/6).
promoters_fwd = ['TATNNT', 'TANANT', 'TANNAT', 'NATANT', 'NATNAT', 'NANAAT']
# promoters_rev = ['ANNATA', 'ANTNTA', 'ATNNTA', 'ANTATN', 'ATNATN', 'ATTNTN']

# Length to extend the search window (upstream and downstream of sRNA TSS)
upstream_margin = 50
downstream_margin = 10

# sample_seq = Seq('CCCGGGGGCATTATAGGGGGGGG', IUPAC.unambiguous_dna)

def search_for_substring(promoter_list, seq_to_search):
    search_array = []
    for promoter in promoter_list:
        search_results = nt_search(str(seq_to_search), promoter)
        search_array.append(search_results)
    return search_array

with open('s_mutans_predicted_promoters.csv', 'w') as outfile:
    writer = csv.writer(outfile, delimiter = '\t')
    writer.writerow(['uid', 'threshold', 'genome', 'feature name', 'strand', 'TSS', 'sequence searched', 'promoter pos relative to TSS, promoter seq'])
    for entry in TSS_list:
        uid = entry[0]
        threshold = entry[1]
        genome = entry[2]
        feature_name = entry[3]
        strand = entry[4]
        TSS = entry[5]
        if strand == 'F':
            window_start_pos = TSS - 1 - upstream_margin
            window_end_pos = TSS + downstream_margin
            window_to_search = ref_sequence[window_start_pos:window_end_pos]
            print('feature name: ', feature_name)
            print('length of search window: ', len(window_to_search))
            results = search_for_substring(promoters_fwd, window_to_search)
            print(results)
            results_array = []
            for pos_list in results:
                num_results = len(pos_list) - 1
                if num_results >= 1:
                    for i, entry in enumerate(pos_list, start=1):
                        if i < len(pos_list):
                            # the promoter position via pythonic counting
                            promoter_pos = window_start_pos + (pos_list[i])
                            distance_to_TSS = promoter_pos + 1 - TSS
                            promoter_seq = ref_sequence[promoter_pos:promoter_pos + 6]
                            results_array.append((distance_to_TSS, str(promoter_seq)))
            if len(results_array) > 0:
                unique_results = set(results_array)
                writer.writerow([uid, threshold, genome, feature_name, strand, TSS, str(window_to_search), unique_results])
        elif strand == 'R':
            window_start_pos = TSS - 1 - downstream_margin
            window_end_pos = TSS + upstream_margin
            window_to_search = ref_sequence[window_start_pos:window_end_pos].reverse_complement()
            window_to_search_fwd = ref_sequence[window_start_pos:window_end_pos]
            print(window_to_search_fwd)
            print(window_to_search)
            print('feature name: ', feature_name)
            print('length of search window: ', len(window_to_search))
            results = search_for_substring(promoters_fwd, window_to_search)
            print(results)
            results_array = []
            for pos_list in results:
                num_results = len(pos_list) - 1
                if num_results >= 1:
                    for i, entry in enumerate(pos_list, start=1):
                        if i < len(pos_list):
                            promoter_pos = window_end_pos - pos_list[i]
                            distance_to_TSS = TSS - promoter_pos
                            promoter_seq = ref_sequence[(promoter_pos - 6):promoter_pos].reverse_complement()
                            results_array.append((distance_to_TSS, str(promoter_seq)))
            if len(results_array) > 0:
                unique_results = set(results_array)
                writer.writerow([uid, threshold, genome, feature_name, strand, TSS, str(window_to_search), unique_results])
        else: print('strand error! please input forward strand as F, reverse strand as R')

