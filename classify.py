#! /usr/bin/env python
"""
Emulate the kraken k-mer classification step, but on banded MinHash signatures.

Briefly,

* load in the k-mer-to-lineage database produced by 'extract.py'.
* for every hash in the signature, find the lineage for that k-mer
  (here 'lineage' would be the computed last-common-ancestor from the NCBI
  taxonomy, based on GenBank genomes)
* count & summarize

Usage:

   kraken/classify.py nodes.dmp names.dmp foobar ecoli_many_sigs/ecoli-1.sig

where 'foobar' is the --savename from extract.py, and nodes/names.dmp are from
NCBI's taxdump.tar.gz --

    curl -O -L ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
    tar xzf taxdump.tar.gz nodes.dmp names.dmp

TODO:

* add classification of FASTA/FASTQ?
"""

import sourmash_lib, sourmash_lib.signature
import argparse
from pickle import load
import collections
from ncbi_taxdump_utils import NCBI_TaxonomyFoo
import os
import json
import gzip

SCALED=10000                              # should match the LCA compute @CTB

kraken_rank_code = {
    'genus' : 'G',
    'species': 'S',
    'phylum': 'P',
    'class': 'C',
    'order': 'O',
    'family': 'F',
    'genus': 'G',
    'kingdom': 'K',
    'domain': 'D' }

def xopen(filename, mode):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    return open(filename, mode)


def load_lca_json(filename):
    with xopen(filename, 'rt') as json_fp:
        info = json.load(json_fp)
        assert info['version'] == 1

    info['basepath'] = os.path.dirname(filename)

    return info


def get_lca_info(lca_info, ksize, scaled):
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
    hashval_to_lca = load(xopen(lca_file, 'rb'))

    return taxfoo, hashval_to_lca, entry['scaled']


def main():
    p = argparse.ArgumentParser()
    p.add_argument('lca_json')
    p.add_argument('sigfiles', nargs='+')
    p.add_argument('-k', '--ksize', default=31, type=int)
    args = p.parse_args()

    # load lca.json
    info = load_lca_json(args.lca_json)
    taxfoo, hashval_to_lca, scaled = get_lca_info(info, args.ksize, SCALED)
    
    # load signatures
    siglist = []
    print('loading signatures from {} signature files'.format(len(args.sigfiles)))
    for sigfile in args.sigfiles:
        sigs = sourmash_lib.signature.load_signatures(sigfile,
                                                      select_ksize=args.ksize)
        sigs = list(sigs)
        siglist.extend(sigs)

    print('loaded {} signatures total at k={}'.format(len(siglist), args.ksize))

    # downsample
    print('downsampling to scaled value: {}'.format(scaled))
    for sig in siglist:
        if sig.minhash.scaled < scaled:
            sig.minhash = sig.minhash.downsample_scaled(scaled)

    # now, extract hash values!
    hashvals = collections.defaultdict(int)
    for sig in siglist:
        for hashval in sig.minhash.get_mins():
            hashvals[hashval] += 1

    found = 0
    total = 0
    by_taxid = collections.defaultdict(int)

    # for every hash, get LCA of labels
    for hashval, count in hashvals.items():
        lca = hashval_to_lca.get(hashval)
        total += count

        if lca is None:
            by_taxid[0] += count
            continue

        by_taxid[lca] += count
        found += count

    print('found LCA classifications for', found, 'of', total, 'hashes')
    not_found = total - found

    # now, propogate counts up the taxonomic tree.
    by_taxid_lca = collections.defaultdict(int)
    for taxid, count in by_taxid.items():
        by_taxid_lca[taxid] += count

        parent = taxfoo.child_to_parent.get(taxid)
        while parent != None and parent != 1:
            by_taxid_lca[parent] += count
            parent = taxfoo.child_to_parent.get(parent)

    total_count = sum(by_taxid.values())

    # sort by lineage length
    x = []
    for taxid, count in by_taxid_lca.items():
        x.append((len(taxfoo.get_lineage(taxid)), taxid, count))

    x.sort()

    # ...aaaaaand output.
    print('{}\t{}\t{}\t{}\t{}\t{}'.format('percent', 'below', 'at node',
                                          'code', 'taxid', 'name'))
    for _, taxid, count_below in x:
        percent = round(100 * count_below / total_count, 2)
        count_at = by_taxid[taxid]

        rank = taxfoo.node_to_info.get(taxid)
        if rank:
            rank = rank[0]
            classify_code = kraken_rank_code.get(rank, '-')
        else:
            classify_code = '-'

        name = taxfoo.taxid_to_names.get(taxid)
        if name:
            name = name[0]
        else:
            name = '-'

        print('{}\t{}\t{}\t{}\t{}\t{}'.format(percent, count_below, count_at,
                                              classify_code, taxid, name))

    if not_found:
        classify_code = 'U'
        percent = round(100 * not_found / total_count, 2)
        count_below = not_found
        count_at = not_found
        taxid = 0
        name = 'not classified'

        print('{}\t{}\t{}\t{}\t{}\t{}'.format(percent, count_below, count_at,
                                              classify_code, taxid, name))


if __name__ == '__main__':
    main()
