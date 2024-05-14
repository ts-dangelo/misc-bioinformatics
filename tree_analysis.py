''' 

This program calculate verticality and 'sister diversity' metrics for a set of rooted
phylogenetic trees.  The trees can be made via an ortholog identifying pipeline and a phylogeny
constructing tool.  Rooting by the Minimal Ancestor Deviation (M.A.D) method can be used. 

Tree tips must be formatted like this: GENOMENAME_GENEINFORMATION. The genome name must be
the first string of the tip ID, with the gene info being the last string in the tip ID seperated by an
underscore. 

Example: GCA_011056015_403 - GCA_011056015 is the genome name and 403 is the gene info.  

The genome name GCA_011056015 needs to be in the metadata file under column "ID" and the Phylum 
of the genome under a column named "Phylum"

'''


import sys, os.path
import pandas as pd
from collections import defaultdict
import glob
from scipy.stats import iqr
import numpy as np
from ete3 import Tree
import argparse
from skbio.diversity import alpha_diversity

parser = argparse.ArgumentParser(description="Calculating tree metrics")
parser.add_argument('-trees', help="Path of directory containing the trees, and nothing else - the trees must be rooted already.")
parser.add_argument('-metadata', help="Tab seperated text metadata file that needs to have at least the two columns 'ID' and 'Phylum'. See description")
parser.add_argument('-outname', help="Prefix for output files")
args = parser.parse_args()
arg_dict = vars(args)

tree_path = arg_dict['trees']
tree_paths = glob.glob(os.path.join(tree_path, "*"))

mdat = arg_dict['metadata']

prefix = arg_dict['outname']
abs_path = os.path.abspath(os.getcwd()) 
out_pref = os.path.join(abs_path, prefix)

mdata = pd.read_csv(mdat, sep = "\t")


def tree_stats(mdata, tree_paths):

	tree_df = pd.DataFrame(columns = ['gene_cluster_id', 'mean_rt_dist', 'median_rt_dist', 'rt_dist_std.dev', 'num_seqs', 'num_phyla', 'num_clades', "simpson_index", 'num_mono_phyla', 'mono_phyla', 'frac_genomes_in_mono', 'monophyly_score'])

	phyla_l = list(set(mdata.Phylum.to_list()))
	phyla_dicts = defaultdict(list)
	for f in phyla_l:
		phy_list = []
		for row in mdata.itertuples():
			gen_id = row.ID
			phy = row.Phylum
			if phy is f:
				phy_list.append(gen_id)
			else:
				continue
		phyla_dicts[f] = phy_list

	for tree in tree_paths:
		(dir, file) = os.path.split(tree)
		gc_id = file.rsplit('.',4)[0] #double check this
		print(gc_id)
		bl_list = []
		num_leaves = 0
		t = Tree(tree)
		madR = t.get_tree_root()
		for node in t.traverse():
			if node.is_leaf():
				num_leaves += 1
				name = node.name
				bl = madR.get_distance(node)
				bl_list.append(bl)
		arr = np.array(bl_list)
		mean = np.mean(arr)
		median = np.median(arr)
		std = np.std(arr)
		
		genome_list=[]
		for leaf in t.iter_leaves():
			header = leaf.name
			splits = header.rsplit("_")[0:-1]
			nhead = '_'.join(splits)
			genome_list.append(nhead)
			for key, value in phyla_dicts.items():
				if nhead in value:
					leaf.add_feature("Phyla", key)
		print(genome_list)
			
		phyla_list = []
		for genome in genome_list:
			for key, value in phyla_dicts.items():
				if genome in value:
					phyla_list.append(key)
					
		uniq_phyla = list(set(phyla_list))
		clade_min = len(uniq_phyla)
		
		print(uniq_phyla)

		num_mono_clades = 0
		for phyla in uniq_phyla:
			print(phyla)
			for node in t.get_monophyletic(values=[phyla], target_attr="Phyla"):
				#print(node)
				num_mono_clades += 1
				print(num_mono_clades)
				
		num_mono_phyla = 0
		num_mono_leaves = 0
		mono_phyla = []
		
		for phyla in uniq_phyla:
			mono_test = t.check_monophyly(values=[phyla], target_attr="Phyla")
			if str(mono_test[0]) == "True":
				mono_phyla.append(phyla)
				num_mono_phyla += 1
				for leaf in t.iter_leaves():
					nphy = str(leaf.get_ascii(attributes=["Phyla"], show_internal=False))
					nphy = nphy.replace("--", "")
					nphy = nphy.strip()
					if nphy == phyla:
						num_mono_leaves += 1
					else:
						continue
						
		mono_phyla_output = ', '.join([str(x) for x in mono_phyla])
		
		totals_list = []
		for phyla in uniq_phyla:
			num_leafs_per_phyla = 0 
			for leaf in t.iter_leaves():
				nphy = str(leaf.get_ascii(attributes=["Phyla"], show_internal=False))
				nphy = nphy.replace("--", "")
				nphy = nphy.strip()
				if nphy == phyla:
					num_leafs_per_phyla += 1
			totals_list.append(num_leafs_per_phyla)
			
		totals_list.sort(reverse=True)
		diversity_list = totals_list #copy of list before iteration for simpsons index
		
		min_denom = 0
		
		if len(totals_list) == 2:
		
			min_denom += (2 * (min(totals_list)) + 1)
			
		else:
		
			while len(totals_list) > 2:
			
				diff = totals_list[0] - totals_list[1]
				min_denom += (2 * totals_list[1])
				
				totals_list.remove(max(totals_list))
				totals_list.remove(max(totals_list))
				totals_list.append(diff)
				totals_list.sort(reverse=True)
				
			if len(totals_list) == 2:
				min_denom += (2 * (min(totals_list)) + 1)
				
		SI_1D = alpha_diversity('simpson', diversity_list)
		simpson = SI_1D[0]
		unscaled_mono_score = clade_min / num_mono_clades
		if min_denom == 0:
			min_mono = 1
		else:
			min_mono = clade_min / min_denom
		mono_score = (unscaled_mono_score * (1 - min_mono)) + min_mono #linear scaling monophyly score using theoretical lowest score for every given tree, puts them on all on the same 0-1 scale
		frac_leaves_in_mono = num_mono_leaves / num_leaves
		
		tree_df = tree_df.append({'gene_cluster_id' : gc_id, 'mean_rt_dist' : mean, 'median_rt_dist' : median, 'rt_dist_std.dev' : std, 'num_seqs' : num_leaves, 'num_phyla' : clade_min, 'num_clades' : num_mono_clades, 'simpson_index' : simpson, 'num_mono_phyla' : num_mono_phyla, 'mono_phyla' : mono_phyla_output, 'frac_leaves_in_mono' : frac_leaves_in_mono, 'monophyly_score' : mono_score}, ignore_index = True)

    
	return(tree_df, phyla_l, phyla_dicts)


def sister_analysis(phya_l, tree_paths, phyla_dicts):

	phyla_sister_counts = {}
	for tree in tree_paths:
		(dir, file) = os.path.split(tree)
		t = Tree(tree)
		gc_id = file.rsplit('.',4)[0] #double check this
		
		genome_list=[]
		for leaf in t.iter_leaves():
			header = leaf.name
			splits = header.rsplit("_")[0:-1]
			nhead = '_'.join(splits)
			genome_list.append(nhead)
			for key, value in phyla_dicts.items():
				if nhead in value:
					leaf.add_feature("Phyla", key)
					
		phyla_list = []
		for genome in genome_list:
			for key, value in phyla_dicts.items():
				if genome in value:
					phyla_list.append(key)
					
		uniq_phyla = list(set(phyla_list))
		
		sister_count = defaultdict(int)
		for phyla in uniq_phyla:
			num_mono_clades = 0
			sister_list = []
			for node in t.get_monophyletic(values=[phyla], target_attr="Phyla"):
				num_mono_clades += 1
				sisters = node.get_sisters()
				for sister_node in sisters:
					for leaf in sister_node.iter_leaves():
						taxonomy = leaf.Phyla
						if taxonomy == phyla: #skipping counting the same phyla as a sister of itself
							pass
						else:
							sister_list.append(taxonomy)
							
			count = len(list(set(sister_list)))
			sister_count[phyla] = count
			phyla_sister_counts[gc_id] = sister_count
	sister_df = pd.DataFrame.from_dict(phyla_sister_counts, orient='index')
	#sister_df = sister_df.fillna(method='ffill')
	sister_df.loc['mean'] = sister_df.mean()
	
	return(sister_df)

def main():

	tree_df, phya_l, phyla_dicts = tree_stats(mdata, tree_paths)
	sister_df =  sister_analysis(phya_l, tree_paths, phyla_dicts)
	
	with open("%s_tree-metrics.txt" % (out_pref), "w") as outfile:
		tree_df.to_csv(outfile, index = False, header=True, sep='\t')
		print("Tree metrics written to " + "%s_tree-metrics.txt" % (out_pref))
	with open("%s_sister-analysis.txt" % (out_pref), "w") as outfile:
		sister_df.to_csv(outfile, index = True, header=True, sep='\t', na_rep='NA')
		print("Sister analysis results written to " + "%s_sister-analysis.txt" % (out_pref))
	
if __name__ == "__main__":
	main()	
	
