"""
This program runs KOFAMSCAN (-mapper option) on a set of input protein fasta files. It then clusters the genomes by similarity based on their KEGG annotations. 
Written by Tim D'Angelo (tdangelo@bigelow.org)
Fasta file format requirements:
Please format fasta files with a simple name that is one string with no spaces or special characters besides dashes (-) or underscores (_). The shorter and simpler the name, the better.
The fasta file name should be formatted so that everything before .faa is what you want your genomes to be called in downstream figures.
Please ensure the fasta headers are formatted in a simple way that is one string of text.  It is best that the individual ORFs are formatted in a way that is easily identifiable for downstream parsing.
Example of an ideal fasta header: 
>GenomeName_Orf1
 
If you are using genomes downloaded from a database like NCBI or IMG it is probably best to use the identifiers used on those databases for future reference:
>ImgGenomeeID_ImgOrfID 
The resulting KOFAMSCAN annotation files are used to build a Genome x KEGG ID count table, with KEGG ID's missing from a given genome counted as 0
The KEGG content table is an output (PREFIX__KO_Table.txt).  Additionally, the genomes are hierarchically clustered based on KEGG content using the euclidean distance of KEGG presence/absence data and the Wards linkage method.
This program requires that you have KOFAMSCAN installed in your path. This can be set up using Anaconda: https://anaconda.org/bioconda/kofamscan 
Be sure to load the Conda environment that KOFAMSCAN is active prior to running this script.
The necessary HMM profiles and KO IDs to run KOFAMSCAN can be downloaded from this link: ftp://ftp.genome.jp/pub/db/kofam/
Input arguments allow you to point to where you have stored those files for running KOFAMSCAN

For the Aminicenantes phylum, This workflow was tested against an annotation-independent workflow that clusters genomes based on protein families (links below). The sequence similarity calculation followed by MCL clustering used is similar to the pangenome workflow in the Anvio program.

https://github.com/raphael-upmc/proteinClusteringPipeline and in this publication: https://www.nature.com/articles/s41467-019-12171-z

The Mantel test in the R library Vegan was used to compare the data resulting from this KEGG annotation based genome-clustering, and annotation-independent protein-family based genome-clustering mentioned above. The comparison shows that the resulting data from both workflows is highly correlated (R2 = 0.9332 for Euclidean distance and 0.8625 for the Jaccard distance), and thus the intepretations will be similar. In the analysis of the Amincenantes Phylum it was found that an average of 46% of ORFs in the genomes were able to be assigned a KO Id by KOFAMSCAN. These Mantel test results show that although the annotation based analysis is only using ~half of the ORF data, the resulting relationships between the genomes are basically the same. Although ommitting large amounts of genome information may be a con to this workflow, those proteins are likely hypothetical or have no known function. The pro to this annotation based analsysis is that the underlying data (KEGG IDs) represent putative functions that have a defined confidence (the adaptive score threshold described in the KofamKOALA publication).


"""

def	id_table(kpaths):
	
	import sys, os.path
	from collections import defaultdict
	import pandas as pd
	
	gen_KOs = {} 

	for txt in kpaths:
		KO_count = defaultdict(int)
		(dir, file) = os.path.split(txt)
		genome_name = file.rsplit('_',2)[0] 
		clean_dat = {}	
		for line in open(str(txt), "r"):	
			line = line.rstrip()
			fields = line.split()
			if len(fields) > 1:
				clean_dat[fields[0]] = fields[1]
			else:
				continue

		kofamscan_output = pd.DataFrame.from_dict(clean_dat, orient = 'index')		 
		kofamscan_output.rename( columns={0 :'KO'}, inplace=True ) 
		ko_list = kofamscan_output['KO'].to_list()
		unique_ko = list(set(ko_list))
	       
		for ko in ko_list:
			if ko in unique_ko:
				KO_count[ko] += 1                   
			else:
				continue
		gen_KOs[genome_name] = KO_count

		kegg_table = pd.DataFrame.from_dict(gen_KOs, orient='index')
		kegg_table = kegg_table.fillna(0)   
                     
	return kegg_table



def plot(kegg_table, outfig):

	import seaborn as sns
	import matplotlib.pyplot as plt
	
	kegg_table[kegg_table > 0] = 1 # set table to binary data
	
	sns.set(font_scale=0.35)
	sns.set_style({"savefig.dpi": 300})
	g = sns.clustermap(kegg_table, cmap="binary",vmax =1 , col_cluster=True, figsize=(12, 8), method = "ward", metric = "euclidean", dendrogram_ratio=(.1, .15))
	g.cax.set_visible(False)
	g.ax_col_dendrogram.set_visible(False)
	g.ax_heatmap.xaxis.set_visible(False)
	plt.savefig(outfig, dpi= 300, orientation='landscape', format = "png")


def main():

	import glob
	import argparse
	import pandas as pd
	import subprocess
	import os,shutil
	from collections import defaultdict
	
	parser = argparse.ArgumentParser(description="Create table of KEGG ID count per genome using translated protein fasta files as input")
	parser.add_argument('-i', help="Path to a directory containing the protein fasta files and no other files. Please format fasta files with a simple name that is one string with no spaces or special characters besides dashes or underscores..")
	parser.add_argument('-annot_dir', help="Directory for KOMFAMSCAN output files")
	parser.add_argument('-prefix', help="Prefix for output file names")
	parser.add_argument('-db_dir', help="directory containing KOFAMSCAN profiles")
	parser.add_argument('-ko_list', help ="path ko_list file from KOFAMSCAN, inlcude file name: /path/to/ko_list ")
	args = parser.parse_args()
	arg_dict = vars(args)
	
	outname = arg_dict['prefix']
		
	fpath = arg_dict['i']
	fpaths = glob.glob(os.path.join(fpath, "*"))
	
	koloc = arg_dict['ko_list']
	profloc = arg_dict['db_dir']
	
	k_dir = arg_dict['annot_dir']
	abs_path = os.path.abspath(os.getcwd()) 
	annot_dir = os.path.join(abs_path, k_dir) 
	if not os.path.exists(annot_dir):
		os.makedirs(annot_dir)
	
	for faa in fpaths:
		(dir, file) = os.path.split(faa)
		genome = file.rsplit('.',1)[0]
		print ("Annotating "  + genome)		
		subprocess.call(["exec_annotation", "-f", "mapper", "--cpu=8", "-p", profloc, "-k", koloc, "-o", "%s_annotation.txt" % (genome), faa])
		for f in os.listdir(abs_path):
			if f.endswith("_annotation.txt"):
				annot_file = os.path.join(abs_path, f)
				shutil.move(annot_file, annot_dir)
   		
   	kpaths = glob.glob(os.path.join(annot_dir, "*"))
	out_table = id_table(kpaths)
	with open("%s_KO_Table.txt" % (outname), "w") as outfile: 
		out_table.to_csv(outfile, index = True, header=True, sep='\t')
	outfig = os.path.join(abs_path, outname + ".png")
	plot(out_table, outfig)         
			
	
if __name__ == "__main__":
	main()
