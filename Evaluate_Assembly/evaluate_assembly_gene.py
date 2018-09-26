import re
import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from collections import OrderedDict
from Bio.Alphabet import IUPAC

Gene_all = []
GeneSeq_length = []

def get_AllGene_length(file):
    GeneDict = OrderedDict()
    re_name = []
    result = []
    for seq_record in SeqIO.parse(file,"fasta"):
        GeneDict[seq_record.id] = str(seq_record.seq)
    trans_str = ' '.join(list(GeneDict.keys()))
    transName = list(GeneDict.keys())
    for name in transName:
        re_name.append(name.split('i')[0])

    def filter_id(ID):
        final_name = ''
        p = ID+'i\d+'
        ls = re.findall(p, trans_str)
        if len(ls) == 1:
            final_name = ls[0]
        else:
            max_len = max([len(GeneDict[i]) for i in ls])
            for i in ls:
                if len(GeneDict[i]) == max_len:
                    final_name = i
        return final_name

    for ID in set(re_name):
        result.append(filter_id(ID))

    #Get the length of each gene sequence and length
    for i in result:
        GeneSeq_length.append(len(GeneDict[i]))
        Gene_all.append(GeneDict[i])

def count_GC():
    long_seq = ''.join(Gene_all)
    seq = Seq(long_seq, IUPAC.unambiguous_dna)
    GC_content = 100 * ((seq.count("G") + seq.count("C")) / len(long_seq))
    return '%.2f%%' % GC_content

def get_average_length():
    average_length = sum(GeneSeq_length)//len(GeneSeq_length)
    return average_length

#count contig length more than 500 or 2000
def count_length():
    num500 = 0
    num2000 = 0
    for i in GeneSeq_length:
        if i > 500:
            num500 += 1
    for i in GeneSeq_length:
        if i > 2000:
            num2000 += 1
    return num500,num2000

#count Nxx
def Nxx(percentage):
    GeneSeq_length_sort = list(reversed(sorted(GeneSeq_length)))
    N = int(sum(GeneSeq_length) * percentage)
    sum_base = 0
    i = 0
    while True:
        sum_base += GeneSeq_length_sort[i]
        if sum_base < N:
            i += 1
        else:
            return GeneSeq_length_sort[i], i+1
            break


def main():
    parser = argparse.ArgumentParser(description="Evaluate the results of the second \
                                     generation genome sequencing assembly")
    parser.add_argument("-i", '--input', help="for example genome.utg.fasta", required=1)
    args = parser.parse_args()
    get_AllGene_length(args.input)
    average_length = get_average_length()
    num500, num2000 = count_length()
    N50 = Nxx(0.5)
    N60 = Nxx(0.6)
    N70 = Nxx(0.7)
    N80 = Nxx(0.8)
    N90 = Nxx(0.9)
    GC = count_GC()
    print('Percent GC: %s' % GC)
    print('Number of sequence: %d' % len(Gene_all))
    print('Min length: %d' % min(GeneSeq_length))
    print('Max length: %d' % max(GeneSeq_length))
    print('Total size: %d' % sum(GeneSeq_length))
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

