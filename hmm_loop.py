#!/usr/bin/python

"""
Identifying genes of interest in a batch of genomes using provided HMM profiles.
Requires genome files that have previously been converted to CDS amino acid format.
Must have HMMsearch loaded in path.


"""


import sys, os.path
import glob
import argparse
from collections import defaultdict
import subprocess


parser = argparse.ArgumentParser(description="Loop over a folder of genomes and run HMMSEARCH from a folder of HMM profiles")
parser.add_argument('-gen_path', help="Path to folder with ONLY the genomes you want searched - In translated amino acid format")
parser.add_argument('-hmm_path', help="Path to folder with ONLY the HMM profiles you want used for searching")
parser.add_argument('--cut_tc',action='store_true',help="This option only works of your HMM profiles contain a 'Trusted Cutoff', this is not true for all HMM models depending on the available data for the protein group")
parser.add_argument('-E',help='Set E-Value to be used in hmmsearch. Default: 1E-5',default='1E-5')
args = parser.parse_args()
arg_dict = vars(args)

fpath = arg_dict['gen_path']
fpaths = glob.glob(os.path.join(fpath, "*"))

hpath = arg_dict['hmm_path']
hpaths = glob.glob(os.path.join(hpath, "*"))

Evalue = arg_dict['E']

for faa in fpaths:
        (dir, file) = os.path.split(faa)
        gname = file.rsplit('.',1)[0]
        for hmm in hpaths:
                (dir, file) = os.path.split(hmm)
                hname = file.rsplit('.',1)[0]
                hmmer_log = open("hmmsearhc-log.txt","w")
                if arg_dict["cut_tc"] == "True":
                        subprocess.call(["hmmsearch","--cut_tc","--tblout", str(str(gname)+"_"+str(hname)+".HMM.txt"), "--notextw", hmm, str(faa)],stdout=hmmer_log)
                else:
                        subprocess.call(["hmmsearch","-E",Evalue,"--tblout", str(str(gname)+"_"+str(hname)+".HMM.txt"),"--notextw", hmm, str(faa)],stdout=hmmer_log)
