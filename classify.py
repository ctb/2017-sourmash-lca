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
* worry about --scaled factors here and in extract.
"""

import sourmash_lib, sourmash_lib.signature
import argparse
from pickle import load
import collections
from ncbi_taxdump_utils import NCBI_TaxonomyFoo

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

def main():
    p = argparse.ArgumentParser()
    p.add_argument('nodes_dmp')
    p.add_argument('names_dmp')
    p.add_argument('saved_db')
    p.add_argument('sigfile')
    p.add_argument('-k', '--ksize', default=31)
    args = p.parse_args()

    taxfoo = NCBI_TaxonomyFoo()
    
    # load the nodes_dmp file to get the tax tree
    print('loading taxonomic nodes from:', args.nodes_dmp)
    taxfoo.load_nodes_dmp(args.nodes_dmp)
    print('loading taxonomic names from:', args.names_dmp)
    taxfoo.load_names_dmp(args.names_dmp)
    
    print('loading k-mer DB from:', args.saved_db)
    hashval_to_lca = load(open(args.saved_db, 'rb'))
    
    # load signature
    sig = sourmash_lib.signature.load_one_signature(args.sigfile,
                                                    select_ksize=args.ksize)
    sig.minhash = sig.minhash.downsample_scaled(SCALED)

    found = 0
    by_taxid = collections.defaultdict(int)

    # for every hash, print out LCA of labels
    for n, hashval in enumerate(sig.minhash.get_mins()):
        lca = hashval_to_lca.get(hashval)
        if lca is None:
            by_taxid[0] += 1
            continue

        by_taxid[lca] += 1
        found += 1

    print('found LCA classifications for', found, 'of', n + 1, 'hashes')
    not_found = n + 1 - found

    # now, propogate counts up the taxonomic tree.
    by_taxid_lca = collections.defaultdict(int)
    for taxid, count in by_taxid.items():
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
    for _, taxid, count_below in x:
        percent = round(100 * count_below / total_count, 2)
        count_at = by_taxid[taxid]

        rank = taxfoo.node_to_info[taxid][0]
        classify_code = kraken_rank_code.get(rank, '-')

        name = taxfoo.taxid_to_names[taxid][0]

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
