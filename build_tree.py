'''
    Description :
    this script is used to evaluate the results of the second
    generation genome sequencing assembly.
    Date :
    2018/09/01
    author:
    lmz
    '''

import re
import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

# Define a list to accept contig sequences
all_contig = []
# Define a list to store the length of each contig sequence
contig_length = []

def get_all_seq(file):
    for seq_record in SeqIO.parse(file,"fasta"):
        all_contig.append(str(seq_record.seq))

def count_GC():
    long_seq = ''.join(all_contig)
    seq = Seq(long_seq, IUPAC.unambiguous_dna)
    GC_content = 100 * ((seq.count("G") + seq.count("C")) / len(long_seq))
    return '%.2f%%' % GC_content

#Get the length of each contig sequence
def get_length_contig(all_contig):
    for i in all_contig:
        contig_length.append(len(i))

def get_average_length():
    average_length = sum(contig_length)//len(contig_length)
    return average_length

#count contig length more than 500 or 2000
def count_length():
    num500 = 0
    num2000 = 0
    for i in contig_length:
        if i > 500:
            num500 += 1
    for i in contig_length:
        if i > 2000:
            num2000 += 1
    return num500,num2000

#count Nxx
def Nxx(percentage):
    contig_length_sort = list(reversed(sorted(contig_length)))
    N = int(sum(contig_length) * percentage)
    sum_base = 0
    i = 0
    while True:
        sum_base += contig_length_sort[i]
        if sum_base < N:
            i += 1
        else:
            return contig_length_sort[i], i+1
            break


def main():
    parser = argparse.ArgumentParser(description="Evaluate the results of the second \
                                     generation genome sequencing assembly")
    parser.add_argument("-i", '--input', help="for example genome.utg.fasta", required=1)
    args = parser.parse_args()
    print(args)
    get_all_seq(args.input)
    get_length_contig(all_contig)
    average_length = get_average_length()
    num500, num2000 = count_length()
    N50 = Nxx(0.5)
    N60 = Nxx(0.6)
    N70 = Nxx(0.7)
    N80 = Nxx(0.8)
    N90 = Nxx(0.9)
    GC = count_GC()
    print('Percent GC: %s' % GC)                                     
    print('Number of sequence: %d' % len(all_contig))
    print('Min length: %d' % min(contig_length))
    print('Max length: %d' % max(contig_length))
    print('Total size: %d' % sum(contig_length))
    print('N90: %d' % N90[0])
    print('N80: %d' % N80[0])
    print('N70: %d' % N70[0])
    print('N60: %d' % N60[0])
    print('N50: %d' % N50[0])
    print('L50: %d' % N50[1])
    print('Average length: %d' % average_length)
    print('Total number (>500bp): %d' % num500)
    print('Total number (>2000bp): %d' % num2000)

if __name__ == "__main__":
    main()
