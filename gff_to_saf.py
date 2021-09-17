'''
Subreads featureCounts doesn't like GFF3 formatting, which is used by Prokka.  This script reformats Prokka GFF3s into the SAF format used by featureCounts. 
'''

import sys, os.path
import pandas as pd
import glob
import argparse

parser = argparse.ArgumentParser(description="Create SAF files from GFF files for using FeatureCounts")
parser.add_argument('-gff_dir', help="Directory of GFF files produced by Prokka")
parser.add_argument('-out_dir', help="Name of directory for resutling SAF files")

args = parser.parse_args()
arg_dict = vars(args)

gff_dir = arg_dict['gff_dir']
out_dir = arg_dict['out_dir']

abs_path = os.path.abspath(os.getcwd()) 
wo_dir = os.path.join(abs_path, out_dir)
if not os.path.exists(wo_dir):
    os.makedirs(wo_dir)
    
    
gffs = glob.glob(os.path.join(gff_dir, "*"))

for gff in gffs:
    SAF_temp = pd.DataFrame(columns = ['GeneID', 'Chr', 'Start', 'End', 'Strand'])
    rfile = open(gff, "r")
    (dir, file) = os.path.split(gff)
    outname = os.path.join(wo_dir, file.rsplit(".")[0])
    for line in rfile.readlines():
        if line.startswith("##"):
            continue
        else:
            tabs = line.rsplit("\t")
            if len(tabs) == 9:
                Chr = tabs[0]
                Start = tabs[3]
                End = tabs[4]
                Strand = tabs[6]
                GeneID = tabs[8].rsplit(";")[0].rsplit("=")[1]
                SAF_temp = SAF_temp.append({'GeneID' : GeneID, 'Chr' : Chr, 'Start': Start, 'End' : End, 'Strand' : Strand}, ignore_index=True)
    SAF_temp.to_csv("%s.SAF" % (outname), index = False, header = True, sep = "\t")
