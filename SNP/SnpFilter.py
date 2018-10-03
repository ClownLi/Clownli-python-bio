'''
Deletion criteria：
a.Samtools.vcf and gatk.vcf have overlapping parts, take gatk.vcf inside
b.5 <= DP <= 100
c.SNP distance (default 5 bp)
'''
import argparse


def get_line(file):
    list_all = []
    with open(file, 'r') as f1:
            con = f1.readlines()
            for line in con:
                line = line.strip()
                if line[0] != '#':
                    list_all.append(line)
    return list_all

def AddHeader(file1, file2):
    with open(file1, 'r') as fr:
        with open(file2, 'w') as fw:
            con = fr.readlines()
            for line in con:
                line = line.strip()
                if line[0] == "#":
                    fw.write(line + '\n')

def FirstDeal(samtool, gatk):
    first_deal = []
    for x in samtool:
        for y in gatk:
            if x.split('\t')[0] == y.split('\t')[0] and x.split('\t')[1] == y.split('\t')[1]:
                if y.split('\t')[7].split(';')[3].split('=')[0] == 'DP' \
                    and 5 <= int(y.split('\t')[7].split(';')[3].split('=')[1]) <= 100:
                    first_deal.append(y)
    return first_deal

def SecondDeal(file, List):
    result = []
    for i in range(len(List)-1):
        if int(List[i+1].split('\t')[1]) - int(List[i].split('\t')[1]) <= 5:
            result.append(List[i+1])
    with open(file, 'a+') as f:
        for i in List:
            if i not in result:
            	f.write(i + '\n')

def main():
    parser = argparse.ArgumentParser(description="For example: python filter_snp.py -i1 sample.samtools.raw.vcf -i2 sample.gatk.raw.vcf -o gatk_filter.vcf")
    parser.add_argument("-i1", "--input1", metavar='', type=str, required=True, help="sample.samtools.raw.vcf")
    parser.add_argument("-i2", "--input2", metavar='', type=str, required=True, help="sample.gatk.raw.vcf")
    parser.add_argument("-o", "--output", metavar='', type=str, required=True, help="gatk_filter.vcf")
    args = parser.parse_args()
    
    samtool = get_line(args.input1)
    gatk = get_line(args.input2)
    AddHeader(args.input2, args.output)
    firstdeal = FirstDeal(samtool, gatk)
    SecondDeal(args.output, firstdeal)

if __name__ == "__main__":
    main()
