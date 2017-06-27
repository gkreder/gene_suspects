import objects


def find_common_gene_set(orthogroups, clusters_list):
    """
    Meant to take a list of clusters objects (from all the species that
    were passed through to antiSMASH) and an orthogroups object (from an
    orthoFinder run) and find the common gene set between the two.
    Assuming here that antiSMASH gene sets will have
    more complete info (except for species maybe?) since the gbk file
    lookup was probably done automatically through the antiSMASH
    web server

    Inputs:
        orthgroups (orthogroups object)
        cluster_list (list of clusters objects)
    """


    # Grab total set of genes in orthogroups object
    orthogroups_genes = []
    # groups = orthogroups.groups
    for group in orthogroups.groups:
        genes = group.genes
        for gene in genes:
            orthogroups_genes.append(gene)

    # Grab total set of genes in clusters_list
    clusters_genes = []
    for clusters_obj in clusters_list:
        for cluster in clusters_obj.clusters:
            for gene in cluster.genes:
                clusters_genes.append(gene)

    # find the overlap genes and save them
    common_genes = []
    for o_gene in orthogroups_genes:
        for c_gene in clusters_genes:
            if o_gene.equals(c_gene):
                merged_gene = o_gene.merge(c_gene)
                common_genes.append(merged_gene)

    return common_genes
