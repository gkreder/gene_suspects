import os
import sys
import importlib
import objects
import mapping_funcs
import parsing_funcs

# Set antiSMASH data directories
bashi_data_dir = "/Users/student/Dropbox/UCSF/Fischbach/Bashi/data/"
lac_smash_dir = bashi_data_dir + 'lactobacillus_plantarum_AL935263.2_antiSMASH/'
clos_smash_dir = bashi_data_dir + 'clostridium_aerotolerans_NZ_JHWJ01000001_antiSMASH/'


# grab filenames for antiSMASH cluster .gbk files
lac_cluster_fnames = []
for fname in os.listdir(lac_smash_dir):
    if '.cluster' in fname and '.gbk' in fname:
        lac_cluster_fnames.append(lac_smash_dir + fname)

clos_cluster_fnames = []
for fname in os.listdir(clos_smash_dir):
    if '.cluster' in fname and '.gbk' in fname:
        clos_cluster_fnames.append(clos_smash_dir + fname)


# Make clusters objects from antiSMASH predicted BGCs
lac_fname = lac_smash_dir + 'geneclusters.txt'
lac_clusters = objects.clusters(lac_cluster_fnames)
clos_fname = clos_smash_dir + 'geneclusters.txt'
clos_clusters = objects.clusters(clos_cluster_fnames)
clusters_list = [lac_clusters, clos_clusters]

# Set orthoFinder data directory and filename
of_dir = bashi_data_dir + 'orthoFinder_Results_Jun15/'
groups_csv = of_dir + 'Orthogroups.csv'
orthogroups = objects.orthogroups(groups_csv)

common_genes = mapping_funcs.find_common_gene_set(orthogroups, clusters_list)
