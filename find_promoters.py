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


# Promoter sequences to look for: A2 and T6 always conserved plus 2 out of 4 variable sites (total of 4/6).
promoters = ['TATNNT', 'TANANT', 'TANNAT', 'NATANT', 'NATNAT', 'NANAAT']

# Length to extend the search window (upstream and downstream of sRNA TSS)
upstream_window = 50
downstream_windown = 10

sample_seq = Seq('ATCGTGTCGATAGCTAGACTTAGGCTAGCTGATCGACTAGCTAGCTAGCA', IUPAC.unambiguous_dna)
print(sample_seq[5:20])
print(sample_seq[5:20].reverse_complement())

