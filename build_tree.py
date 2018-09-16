"""
Author: "LMZ"
Time: "2018-09-14"
Discripton: "This script was designed to easily construct the species tree based on the single-copy gene obtained from the OrthoFinder results."
@@@PS：Thank you my friend ～Dong Wei～ Adapted according to his script！！！
"""

import re
import os
from collections import OrderedDict
from Bio import SeqIO
import argparse
import subprocess

class BuildTree():
    
    orthofinder = '/home/research001/program/OrthoFinder-2.0.0/orthofinder'
    MAFFT = 'mafft'
    RAxML = '/home/research001/program/standard-RAxML-8.1.17/raxmlHPC-PTHREADS'
    ASTRAL = '/home/research001/program/Astral/astral.5.6.2.jar'
    TRIMAL = '/home/research001/program/trimAl/source/trimal'
    
    SpeciesID = []
    SingleOrtho = []
    OrthoDict = OrderedDict()
    SeqDict = OrderedDict()
    
    def __init__(self, *args):
        self.args = args
    
    def _run_command(self, cmd):
        return_code = subprocess.call(cmd, shell=True)
        if return_code != 0:
            print("ERROR: Return code {0} when running the following command: {1}".format(return_code, cmd))

    def _merge_protein(self, dir):
        os.chdir(dir)
        cmd = 'cat * > ../all.fas'
        self._run_command(cmd)
        os.chdir("../")

    def _orthofinder(self, dir):
        cmd = BuildTree.orthofinder + ' -f ' + dir
        self._run_command(cmd)

    def _get_SpeciesID(self, file):
        with open(file) as f:
            for line in f:
                BuildTree.SpeciesID.append(line.strip())

    def _get_SingleID(self, file):
        with open(file, 'r') as f:
            for line in f:
                BuildTree.SingleOrtho.append(line.strip())

    def _get_SingleGeneSeq(self, file1, file2):
        with open(file1, 'r') as f:
            for line in f:
                line_split = line.strip().split(":")
                if line_split[0] in BuildTree.SingleOrtho:
                    BuildTree.OrthoDict[line_split[0]] = line_split[1].split()

        for seq_record in SeqIO.parse(file2, 'fasta'):
            BuildTree.SeqDict[seq_record.id] = seq_record.seq

        os.mkdir("SingleGene")
        os.chdir("SingleGene")
        for single_id in BuildTree.SingleOrtho:
            file_name = single_id + ".fas"
            with open(file_name,'w') as f:
                for protein_id in BuildTree.OrthoDict[single_id]:
                    if protein_id in BuildTree.SeqDict.keys():
                        for id in BuildTree.SpeciesID:
                            if protein_id.startswith(id):
                                f.write(">"+id+"\n"+str(BuildTree.SeqDict[protein_id])+"\n")
        os.chdir("../")

    def _prepare_MSA(self, dir_name):
        os.chdir(dir_name)
        for dir in os.listdir():
            if dir[:6] == "Result":
                os.chdir(dir)
        if "SingleCopyOrthogroups.txt" in os.listdir():
            self._get_SingleID("SingleCopyOrthogroups.txt")
            self._get_SingleGeneSeq("Orthogroups.txt","../../all.fas")

    def _MSA_SingleCopyGene(self):
        os.mkdir("SingleGene_MSA")
        os.chdir("SingleGene")
        for file in os.listdir():
            seq_aln_file = os.path.splitext(file)[0] + "_aln.fas"
            seq_aln_trimed_file = os.path.splitext(file)[0] + "_aln_trimed.fas"
            cmd1 = BuildTree.MAFFT + ' ' + file + ' > ../SingleGene_MSA/' + seq_aln_file
            self._run_command(cmd1)
            cmd2 = BuildTree.TRIMAL + ' -in ' + '../SingleGene_MSA/' + seq_aln_file + ' -out ' + '../SingleGene_MSA/' + seq_aln_trimed_file + ' -automated1'
            self._run_command(cmd2)

    def _merge_SingleCopyGene(self):
        os.chdir("../SingleGene_MSA")
        D = {}
        for i in BuildTree.SpeciesID:
            D[i] = ""
        for file in os.listdir():
            if file.split('_')[-1] == "trimed.fas":
                for seq_record in SeqIO.parse(file, 'fasta'):
                    D[seq_record.id] += str(seq_record.seq)

        with open("merged_allSingleGenes.fas", "w") as f:
            for key,value in D:
                fh.write(">" + key + "\n" + value + "\n")
        os.chdir("../")

    def _MLtree_Concatenation(self):
        os.mkdir("Concatenation")
        os.chdir("Concatenation")
        cmd = BuildTree.RAxML + ' -T 10 -f a -N 100 -m PROTGAMMAJTT -x 123456 -p 123456 -s ../SingleGene_MSA/merged_allSingleGenes.fas -n concatenation_out.nwk'
        self._run_command(cmd)

    def main(self, args):
        parameter1 = args.input1
        parameter2 = args.input2
        self._merge_protein(parameter1)
        self._orthofinder(parameter1)
        self._get_SpeciesID(parameter2)
        self._prepare_MSA(parameter1)
        self._MSA_SingleCopyGene()
        self._merge_SingleCopyGene()
        self._MLtree_Concatenation()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
                                     prog="EasySpeciesTree",
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     description='''
                                         Easily construct the ML species tree with all single-copy gene's protein sequences

                                         ''')
    parser.add_argument(
                        '-in1', '--input1',
                        required=True,
                        help="offer the directory of protein sequence")
    parser.add_argument(
                        '-in2', '--input2',
                        required=True,
                        help="offer the name.txt")


    args = parser.parse_args()
    tree = BuildTree(args)
    tree.main(args)

