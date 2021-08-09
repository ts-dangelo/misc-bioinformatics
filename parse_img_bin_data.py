'''
For some perplexing reason, there is no way to download individual IMG produced bins from the website. You can only download the underlying assembled data
and then a file that summarizes which contigs (or scaffolds?) within the assembled .fna file belongs to each bin, which is pretty useless for most applications.
This script takes the .fna file containing all bins and the scaffoldsets_DATE-OF-DOWNLOAD.tab file and splits the one fasta into fastas for the individual bins.
'''


import sys, os.path
from collections import defaultdict
import glob
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse

parser = argparse.ArgumentParser(description="Split up IMG 'scaffold-sets' from metagenomic data into the bins.")
parser.add_argument('-fasta', help="Fasta file with mixed contigs from multiple bins")
parser.add_argument('-list', help="IMG contig to bin summary file. From IMG is has the name scaffoldsets_DATE-OF-DOWNLOAD.tab")

args = parser.parse_args()
arg_dict = vars(args)
fasta = arg_dict['fasta']
lis = arg_dict['list']

with open(lis, "r")as file:
	status = 0
	for i in file.readlines():
		if i.startswith("scaffold"):
			status = 1
			header_list = []
			if status == 1:
				idline = i.rstrip()
				id1 = idline.split('\t')[1]
				img = id1.split('_', 1)[0]
				img2 = id1.split('_', 2)[1]
				img_id = img + "_" + img2
				print(img_id)
		elif "assembled" in i:
			datline = i.rstrip()
			header = datline.split()[2]
			header_list.append(header)
            
		else:
			with open("%s%s" % (img_id, ".fna"), "w") as outfile:
				for seq_record in SeqIO.parse(fasta, "fasta"):
					header_id = str(seq_record.id)
					if header_id in header_list:
						outfile.write(">%s\n" % (header_id) + str(seq_record.seq)+"\n") 
					else:
						continue
				outfile.close()
