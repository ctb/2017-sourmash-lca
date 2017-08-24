import gzip
import csv


want_taxonomy = ['superkingdom', 'phylum', 'order', 'class', 'family', 'genus', 'species']

class NCBI_TaxonomyFoo(object):
    def __init__(self):
        self.child_to_parent = None
        self.node_to_info = None
        self.taxid_to_names = None
        self.accessions = None

    def load_nodes_dmp(self, filename):
        self.child_to_parent, self.node_to_info = parse_nodes(filename)

    def load_names_dmp(self, filename):
        self.taxid_to_names = parse_names(filename)

    def load_accessions_csv(self, filename):
        self.accessions = load_genbank_accessions_csv(filename)

    # get taxid
    def get_taxid(self, acc):
        # @CTB hack hack split off NZ from accession
        if acc.startswith('NZ_'):
            acc = acc[3:]

        info = self.accessions.get(acc)
        if not info:
            return None

        taxid = info['taxid']
        taxid = int(taxid)
        return taxid

    # code to find the last common ancestor from the lineage string
    def find_lca(self, acc_set):
        # empty? exit.
        if not acc_set:
            return 1

        # get the first full path
        taxid = acc_set.pop()
        path = []
        while taxid != 1:
            path.insert(0, taxid)
            taxid = self.child_to_parent.get(taxid, 1)   # @CTB reexamine

        # find the first shared taxid in each follow-on path
        while acc_set:
            taxid = acc_set.pop()

            path2 = []
            while taxid != 1:
                if taxid in path:
                    path2.insert(0, taxid)
                taxid = self.child_to_parent[taxid]

            path = path2

        if path:
            return path[-1]
        return 1

    def get_lineage(self, taxid, want_taxonomy=None):
        taxid = int(taxid)

        lineage = []
        while taxid != 1:
            if taxid not in self.node_to_info:
                print('cannot find taxid {}; quitting.'.format(taxid))
                break
            rank = self.node_to_info[taxid][0]
            name = self.taxid_to_names[taxid][0]
            if not want_taxonomy or rank in want_taxonomy:
                lineage.insert(0, name)
            taxid = self.child_to_parent[taxid]

        return lineage


### internal utility functions
        

def parse_nodes(filename):
    "Parse the NCBI nodes_dmp file."
    child_to_parent = dict()
    node_to_info = dict()

    with open(filename) as fp:
        for n, line in enumerate(fp):
            x = line.split('\t|\t')

            node_id, parent_node_id, rank, embl, div_id, div_flag, gencode, mgc_inherit, mgc_flag, mgc_id, hidden_flag, subtree_flag, comments = x

            node_id = int(node_id)
            parent_node_id = int(parent_node_id)

            child_to_parent[node_id] = parent_node_id
            node_to_info[node_id] = rank, embl, div_id, div_flag, comments

    return child_to_parent, node_to_info




def parse_names(filename):
    taxid_to_names = dict()
    with open(filename) as fp:
        for n, line in enumerate(fp):
            line = line.rstrip('\t|\n')
            x = line.split('\t|\t')

            taxid, name, uniqname, name_class = x
            taxid = int(taxid)

            if name_class == 'scientific name':
                taxid_to_names[taxid] = (name, uniqname, name_class)

    return taxid_to_names


def load_genbank_accessions_csv(filename):
    xopen = open
    if filename.endswith('.gz'):
        xopen = gzip.open

    # load the genbank CSV (accession -> lineage)
    print('loading genbank accession -> lineage info')
    with xopen(filename, 'rt') as fp:
        accessions = {}
        for row in csv.DictReader(fp, fieldnames=['acc', 'taxid', 'lineage']):
            acc = row['acc']
            accessions[acc] = row

    return accessions
