"""
    Description: "I can't remember makeblastdb and blast paramters all the time
                  Therefore develop this script for simple blast"
    Author: LMZ
    Date: 2018-09-22
"""


import argparse
import subprocess


def run_command(cmd):
    return_code = subprocess.call(cmd, shell=True)
    if return_code != 0:
        print("ERROR: Return code {0} when running the following command: {1}".format(return_code, cmd))

def makeBlastDB(programe, database):
    dbtype = ''
    if programe == 'blastn':
        dbtype = 'nucl'
    else:
        dbtype = 'prot'

    cmd = "makeblastdb -in " + database + " -dbtype " + dbtype + " -out dbname"
    run_command(cmd)

def blastProgram(program, input, output, evalue, outfmt, cpu, NumAli):
    cmd = program + " -query " + input + " -db dbname -out " + output + \
          " -evalue " + str(evalue) + " -outfmt " + str(outfmt) + " -num_threads " + str(cpu) + \
          " -num_alignments " + str(NumAli)
    print(cmd)
    run_command(cmd)

def main():
    parser = argparse.ArgumentParser(description="Usage: python simple_blast.py -i query.fasta \
                                    -d database -p blastn or blastp or blastx")
    parser.add_argument("-i", "--input", metavar='', type=str, required=True, help="query.fasta")

    parser.add_argument("-d", "--database", metavar='', type=str , required=True, help="database.fasta")

    parser.add_argument("-p", "--programe", metavar='', type=str, required=True, help="program, could be blastn, blastp or blastx")

    parser.add_argument("-c", "--cpu", type=int, metavar='', choices=[1, 2, 3, 4], default= 1,
                        help="cpu (default 4)")

    parser.add_argument("-e", "--evalue", metavar='', type=str, default= '1e-3',
                        help="evalue (default 1e-3)")

    parser.add_argument("-o", "--output", metavar='', type=str, default='blast.out', help="default: blast.out")

    parser.add_argument("-f", "--outfmt", metavar='', type=int, choices=[6,7], default= 6,
                        help="output format (default 6)")

    parser.add_argument("-n", "--NumAli", metavar='', type=int, default= 250,
                        help="num of alignment (default 250)")

    
    args = parser.parse_args()
    makeBlastDB(args.programe, args.database)
    blastProgram(args.programe, args.input, args.output, args.evalue, args.outfmt, args.cpu, args.NumAli)
    

if __name__ == "__main__":
    main()
