#! /usr/bin/env python
"""
Emulate the kraken database-building step, but using banded MinHash signatures
instead.

Briefly,

* load in a list of accessions-to-taxonomic-lineage information (e.g. such
  as the .csv.gz file from https://github.com/dib-lab/2017-ncbi-taxdump);
* for every hash in each signature, connect that hash to the lineage for that
  signature (signature genbank accession ID -> lineage);
* after loading in all signatures & hashes, find the last-common-ancestor
  lineage for each hash and associate in a dictionary;
* save dictionary.

Usage::

   kraken/extract.py genbank/*.csv.gz nodes.dmp ecoli_many_sigs/ecoli-*.sig --savename foobar

The E. coli signatures used in the command above can be downloaded like so:

   curl -O -L https://github.com/dib-lab/sourmash/raw/master/data/eschericia-sigs.tar.gz

'nodes.dmp' comes from 'taxdump.tar.gz' --

    curl -O -L ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
    tar xzf taxdump.tar.gz nodes.dmp names.dmp

----

This could be pretty easily broken down into multiple steps:

* step 1: link hashes to lineages and save those links to a file;
* step 1b: update this file as many times as you want;
* step 2: find last-common-ancestor taxid for each hash.

This would make adding in new organisms to the LCA database much
faster, though no less memory intensive.

----

TODO:

* refactor this into two stages?
* add --traverse-dir argument a la sourmash
"""

import sourmash_lib, sourmash_lib.signature
import argparse
from pickle import dump
from collections import defaultdict
from ncbi_taxdump_utils import NCBI_TaxonomyFoo


def main():
    p = argparse.ArgumentParser()
    p.add_argument('genbank_csv')
    p.add_argument('nodes_dmp')
    p.add_argument('sigs', nargs='+')
    p.add_argument('-k', '--ksize', default=31, type=int)
    p.add_argument('-s', '--savename', default=None, type=str)
    args = p.parse_args()

    if not args.savename:
        print('no savename, quitting')
        sys.exit(0)

    taxfoo = NCBI_TaxonomyFoo()

    # load the accesions->taxid info
    taxfoo.load_accessions_csv(args.genbank_csv)

    # load the nodes_dmp file to get the tax tree
    print('loading nodes_dmp / taxonomic tree')
    taxfoo.load_nodes_dmp(args.nodes_dmp)

    # track hashval -> set of taxids
    hashval_to_taxids = defaultdict(set)

    # for every minhash in every signature, link it to its NCBI taxonomic ID.
    print('loading signatures & traversing hashes')
    for n, filename in enumerate(args.sigs):
        if n % 10 == 1:
            print('... on signature', n, 'of', filename)
        sig = sourmash_lib.signature.load_one_signature(filename,
                                                      select_ksize=args.ksize)
        acc = sig.name().split(' ')[0]   # first part of sequence name
        acc = acc.split('.')[0]          # get acc w/o version

        taxid = taxfoo.get_taxid(acc)
        mins = sig.minhash.get_mins()

        for m in mins:
            hashval_to_taxids[m].add(taxid)
    print('...done')

    ####

    print('traversing tags and finding last-common-ancestor for {} tags'.format(len(hashval_to_taxids)))

    hashval_to_lca = {}
    found_root = 0

    # find the LCA for each hashval and store.
    for n, (hashval, taxid_set) in enumerate(hashval_to_taxids.items()):
        if n % 1000 == 0:
            print('...', n)

        # find associated least-common-ancestors.
        lca = taxfoo.find_lca(taxid_set)

        if lca == 1:
            found_root += 1

        # save!!
        hashval_to_lca[hashval] = lca
    print('done')
    print('found root {} times'.format(found_root))

    print('saving to', args.savename)
    dump(hashval_to_lca, open(args.savename, 'wb'))


if __name__ == '__main__':
    main()
