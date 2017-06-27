from Bio import SeqIO
import importlib
import objects
import csv

importlib.reload(objects)

# Takes an antiSMASH cluster output gbk file and parses
# for use in creation of an orthogroups object
def parse_cluster_gbk(cluster_fname):
	types = []
	for seq_record in SeqIO.parse(cluster_fname, "genbank"):
		species = seq_record.description
		acc = seq_record.id
		features = seq_record.features
		pfams = []
		genes = []
		gene_index = -1
		for feature_index, feature in enumerate(features):
			if feature.type == 'gene':
				gene_index += 1
				genes.append(objects.gene())
				genes[gene_index].add_gene(feature)
				pfams.append([])
			elif feature.type == 'CDS':
				genes[gene_index].add_CDS(feature)
			elif feature.type == 'PFAM_domain':
				if gene_index != -1:
					pfam_obj = objects.pfam(feature)
					# print(gene_index, feature_index)
					pfams[gene_index].append(pfam_obj)
	return (species, acc, genes, pfams)
