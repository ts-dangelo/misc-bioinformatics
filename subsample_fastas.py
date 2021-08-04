import sys, os.path
from random import sample
from Bio import SeqIO
import glob
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse


parser = argparse.ArgumentParser(description="Subsample fastas to a given number of sequences")
parser.add_argument('-path', help="Path to directory containing the fastas you want subsampled")
parser.add_argument('-nseq', help="Number of sequences to subsample")
parser.add_argument('-out_dir', help="Directory for output files.")

args = parser.parse_args()
arg_dict = vars(args)

fpath = arg_dict['path']
fpaths = glob.glob(os.path.join(fpath, "*"))

nseq = int(arg_dict['nseq'])

o_dir = arg_dict['out_dir'] 
abs_path = os.path.abspath(os.getcwd()) 
wo_dir = os.path.join(abs_path, o_dir) 
if not os.path.exists(wo_dir):
	os.makedirs(wo_dir)

for faa in fpaths:
	(dir, file) = os.path.split(faa)
	clean_out = os.path.join(wo_dir, file)
	hdr_list = []
	with open(clean_out, "w") as outfile:
		for seq_record in SeqIO.parse(faa, "fasta"):
      			hdr = seq_record.id
			hdr_list.append(hdr)
		rsub = sample(hdr_list,nseq)
		for seq_record in SeqIO.parse(faa, "fasta"):
			hdr = seq_record.id
			if hdr in rsub:
				outfile.write(">%s\n" % (hdr) + str(seq_record.seq)+"\n")
			else:
				continue
		outfile.close()
