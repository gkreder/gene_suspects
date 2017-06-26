import importlib
import parsing_funcs
import csv

importlib.reload(parsing_funcs)

class bgc(object):
    """
    A single biosynthetic gene cluster object

    Attributes:
        species (string) - name of species from which clusters came from
        acc (string) - accession number of species
        genes (list) - list of genes in cluster
        pfams (list) - list of pfam terms in cluster
    """

    def __init__(self, cluster_fname):
        parse_results = parsing_funcs.parse_cluster_gbk(cluster_fname)
        (species, acc, genes, pfams) = parse_results
        self.species = species
        self.acc = acc
        self.genes = genes
        self.pfams = pfams

    def get_genes(self):
        genes = []
        for gene in self.genes:
            genes.append(gene.protein_id)
        return genes

class clusters(object):
    """
    Stores relevant information concerning a species's antiSMASH
    found biosynthetic gene clusters

    Attributes:
        species (string) - name of species from which clusters came from
        acc (string) - accession number of species
        clusters (list) - list of bgc objects
    """
    def __init__(self, cluster_fnames):
        clusters = []
        for fname in cluster_fnames:
            bgc_obj = bgc(fname)
            clusters.append(bgc_obj)
        self.clusters = clusters
        if len(clusters) > 0:
            self.species = clusters[0].species
            self.acc = clusters[0].acc



class gene(object):
    """
    Stores information for a found gene

    Attributes:
        species (string)
        locus_tag (string)
        product (string)
        protein_id (string)
        protein_seq (string)
        db_xref (list)
        location (string)
    """

    def add_gene(self, gene):
        qualifiers = gene.qualifiers
        self.locus_tag = qualifiers['locus_tag'][0]
        self.location = gene.location
        if 'gene' in qualifiers:
            self.name = qualifiers['gene'][0]

    def add_CDS(self, CDS):
        qualifiers = CDS.qualifiers
        self.product = qualifiers['product'][0]
        self.protein_id = qualifiers['protein_id'][0]
        self.protein_seq = qualifiers['translation'][0]
        if 'db_xref' in qualifiers:
            self.db_xref = qualifiers['db_xref']
        self.locus_tag = qualifiers['locus_tag'][0]


class pfam(object):
    """pfam object to store relevant information for a given pfam hit taken
    from output .pfd files

    Attributes:
        location (string)
        sequence (string)
        locus_tag (string)
        evalue (float)
        domain (string)
        description (string)
        db_xref (list)
        notes (string)
        score (float)
    """

    def __init__(self, pfam_feature):
        q = pfam_feature.qualifiers
        self.location = pfam_feature.location
        self.sequence = q['translation'][0]
        self.locus_tag = q['locus_tag'][0]
        self.evalue = float(q['evalue'][0])
        self.domain = q['domain'][0]
        self.description = q['description'][0]
        if 'db_xref' in q:
            self.db_xref = q['db_xref']
        self.notes = q['note'][0]
        self.score = float(q['score'][0])


class orthogroup(object):
    """
    Attributes:
        raw (list)
        genes (list)
        species (list)
    """
    def __init__(self, row):
        species_list = []
        gene_list = []
        self.raw = row

        row = row.split('\t')[1 : ]
        # Clean out empty entries
        if row[0] == '':
            row = row[1 : ]

        # Flatten out any entries that are multiple genes
        # separated by a comma
        flattened_row = []
        for x in row:
            temp_split = x.split(',')
            for element in temp_split:
                flattened_row.append(element)
        for x in flattened_row:
            [species, gene_id] = x.strip().split('|')
            species_list.append(species)
            gene_obj = gene()
            gene_obj.protein_id = gene_id
            gene_obj.species = species
            gene_list.append(gene_obj)

        self.species = species_list
        self.genes = gene_list
        self.raw = row



class orthogroups(object):
    """
    Holds information about found orthogroups from clusterFinder output

    Attributes:
        groups (list)
        species (list)
    """

    def __init__(self, txt_fname):
        groups = []
        with open(txt_fname) as f:
            lines = f.readlines()
        for row_index, row in enumerate(lines):
            row = row.strip()
            if row_index == 0:
                species = row.split('\t')
                continue

            group_obj = orthogroup(row)
            groups.append(group_obj)

        self.species = species
        self.groups = groups
