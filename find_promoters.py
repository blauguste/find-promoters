from Bio import SeqIO, SeqFeature
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import nt_search
import csv

#gb_path = 'inputgb.gb'
#ref_record = SeqIO.read(open(gb_path), 'genbank')

# List of sRNA TSSs. Csv with [1] sRNA name, [2] strand, [3] TSS start
TSS_list = []
with open('sample_tss_list.csv', 'r') as TSS_input:
    TSS_reader = csv.reader(TSS_input, delimiter = ',')
    next(TSS_reader, None) # ignore the header row
    for row in TSS_reader:
        TSS_list.append((row[0], row[1], int(row[2])))
print(TSS_list)

# Promoter sequences to look for: A2 and T6 always conserved plus 2 out of 4 variable sites (total of 4/6).
promoters_fwd = ['TATNNT', 'TANANT', 'TANNAT', 'NATANT', 'NATNAT', 'NANAAT']
promoters_rev = ['ANNATA', 'ANTNTA', 'ATNNTA', 'ANTATN', 'ATNATN', 'ATTNTN']

# Length to extend the search window (upstream and downstream of sRNA TSS)
upstream_margin = 2
downstream_margin = 10

sample_seq = Seq('ATCGTGTCGATAGCTAGACTTAGGCTAGCTGATCGACTAGCACTATCTAGCTATATATTGCA', IUPAC.unambiguous_dna)

def search_for_substring(promoter_list, seq_to_search):
    search_array = []
    for promoter in promoter_list:
        search_results = nt_search(str(seq_to_search), promoter)
        search_array.append(search_results)
    return search_array

for entry in TSS_list:
    sRNA_name = entry[0]
    strand = entry[1]
    TSS = entry[2]
    if strand == 'F':
        window_start_pos = TSS - 1 - upstream_margin
        window_end_pos = TSS + downstream_margin
        print('window start: ', window_start_pos)
        window_to_search = sample_seq[window_start_pos:window_end_pos]
        results = search_for_substring(promoters_fwd, window_to_search)
        print(results)
        for pos_list in results:
            num_results = len(pos_list) - 1
            if num_results >= 1:
                for entry in pos_list[1:]:
                    promoter_pos = window_start_pos + pos_list[1]
                    print(promoter_pos)
                    print(sample_seq[promoter_pos:promoter_pos + 6])
    elif strand == 'R':
        window_start_pos = TSS - 1 - downstream_margin
        window_end_pos = TSS + upstream_margin
        window_to_search = sample_seq[window_start_pos:window_end_pos]
        results = search_for_substring(promoters_rev, window_to_search)
        print(results)
        for pos_list in results:
            num_results = len(pos_list) - 1
            if num_results >= 1:
                for entry in pos_list[1:]:
                    promoter_pos = window_start_pos + pos_list[1]
                    print(promoter_pos)
                    print(sample_seq[promoter_pos:promoter_pos + 6])
                      
#results_array.append((sRNA_name, strand, TSS, ))

