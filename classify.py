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

import khmer
import sourmash_lib, sourmash_lib.signature
import argparse
from pickle import load
import collections
from ncbi_taxdump_utils import NCBI_TaxonomyFoo


want_taxonomy = ['superkingdom', 'phylum', 'order', 'class', 'family', 'genus', 'species']

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
    print('loading nodes_dmp / taxonomic tree')
    taxfoo.load_nodes_dmp(args.nodes_dmp)
    taxfoo.load_names_dmp(args.names_dmp)
    
    print('loading k-mer DB')
    hashval_to_lca = load(open(args.saved_db, 'rb'))
    
    # load signature
    sig = sourmash_lib.signature.load_one_signature(args.sigfile,
                                                    select_ksize=args.ksize)

    cnt = collections.Counter()
    found = 0

    # for every hash, print out LCA of labels
    for n, hashval in enumerate(sig.minhash.get_mins()):
        lca = hashval_to_lca.get(hashval)
        if lca is None or lca == 1:
            continue

        # extract the full lineage & format ~nicely
        lineage = taxfoo.get_lineage(lca)
        lineage = ";".join(lineage)

        found += 1
        cnt[lineage] += 1

    print('found database classifications for', found, 'of', n + 1, 'hashes')
    for item, count in cnt.most_common():
        print(count, item)


if __name__ == '__main__':
    main()
