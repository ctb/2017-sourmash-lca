import json
import gzip
import os
from pickle import load

from ncbi_taxdump_utils import NCBI_TaxonomyFoo


class LCA_Database(object):
    def __init__(self, filename=None):
        if filename:
            self.load(filename)
        else:
            self.lca = {}
            self.lca['version'] = 1
            self.lca['dblist'] = []

    def load(self, filename):
        with xopen(filename, 'rt') as json_fp:
            info = json.load(json_fp)
            assert info['version'] == 1

        info['basepath'] = os.path.dirname(filename)

        self.lca = info

    def save(self, filename):
        if 'basepath' in self.lca:
            del self.lca['basepath']

        with open(filename, 'wt') as json_fp:
            json.dump(self.lca, json_fp)

    def add_db(self, ksize, scaled, lca_path, nodes_path, names_path):
        db_info = {}
        db_info['ksize'] = int(ksize)
        db_info['scaled'] = int(scaled)
        db_info['lca_db'] = lca_path
        db_info['nodes'] = nodes_path
        db_info['names'] = names_path

        self.lca['dblist'].append(db_info)

    def get_database(self, ksize, scaled):
        lca_info = self.lca

        assert lca_info['version'] == 1
        basepath = lca_info['basepath']

        matching_ksizes = []
        for db in lca_info['dblist']:
            if db['ksize'] == ksize:
                matching_ksizes.append(db)

        # ignore scaled matching for now, take first one
        assert len(matching_ksizes) == 1
        entry = matching_ksizes[0]

        taxfoo = NCBI_TaxonomyFoo()

        # load the nodes_dmp file to get the tax tree
        nodes_file = os.path.join(basepath, entry['nodes'])
        print('loading taxonomic nodes from:', nodes_file)
        taxfoo.load_nodes_dmp(nodes_file)

        names_file = os.path.join(basepath, entry['names'])
        print('loading taxonomic names from:', names_file)
        taxfoo.load_names_dmp(names_file)

        lca_file = os.path.join(basepath, entry['lca_db'])
        print('loading k-mer DB from:', lca_file)
        with xopen(lca_file, 'rb') as hashval_fp:
            hashval_to_lca = load(hashval_fp)

        return taxfoo, hashval_to_lca, entry['scaled']


### utility functions


def xopen(filename, mode):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    return open(filename, mode)
