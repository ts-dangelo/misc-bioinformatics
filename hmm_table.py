import sys, os.path
import glob
import argparse
from collections import defaultdict
import pandas as pd

parser = argparse.ArgumentParser(description="Create table of Positive CRISPR Hits based on a directory of HMMSearch result files")
parser.add_argument('-input', help="Path to HMMsearch results.  One file for every CRISPR hmm for every genome in the genomes file. Genome name occrus in the file name before the first underscore.")
parser.add_argument('-output', help="Output table file name")
parser.add_argument('-ext', help="File extension after unique CRISPR name")
parser.add_argument('-genomes', help= "Text file with genome name (what is before first underscore in CRISPR HMM result file)")
args = parser.parse_args()
arg_dict = vars(args)

gen_file = arg_dict['genomes']
cpath = arg_dict['input']
cpaths = glob.glob(os.path.join(cpath,  "*%s" % arg_dict['ext']))

crispr_dict = {} # dictionary of KO:COUNT dictionaries for each genome

for line in open(str(gen_file), "r"):
        line = line.rstrip()
        genome = line
        print(genome)
        crispr_count = defaultdict(int)
        for txt in cpaths:
                (dir, file) = os.path.split(txt)
                genome_name = file.rsplit('_',3)[0] #define genome name, to be the key in the resulting dictionary
                if genome_name != genome:
                        pass
                else:
                        #print(genome_name)
                        gen, crisp = file.split('_', 1)
                        crispr_name = crisp.rsplit('.', 2)[0] 
                        for line in open(str(txt), "r"):
                                line = line.rstrip()
                                if line[0] == "#":
                                        pass
                                else:
                                        crispr_count[crispr_name] += 1   
                                        crispr_dict[genome_name] = crispr_count
#print(crispr_count)            

        
                
        
out = arg_dict['output']

with open("%s.txt" % (out), "w") as outfile:
        df = pd.DataFrame.from_dict(crispr_dict, orient='index')
        df = df.fillna(0)
        print(df)          
        df.to_csv(outfile, index = True, header=True, sep='\t')
