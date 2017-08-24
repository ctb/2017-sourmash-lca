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
import sys, os
import sourmash_lib, sourmash_lib.signature
import argparse
from pickle import dump
from collections import defaultdict
from ncbi_taxdump_utils import NCBI_TaxonomyFoo


def traverse_find_sigs(dirnames):
    for dirname in dirnames:
        for root, dirs, files in os.walk(dirname):
            for name in files:
                if name.endswith('.sig') or name.endswith('.sbt'):
                    fullname = os.path.join(root, name)
                    yield fullname


def main():
    p = argparse.ArgumentParser()
    p.add_argument('genbank_csv')
    p.add_argument('nodes_dmp')
    p.add_argument('sigs', nargs='+')
    p.add_argument('-k', '--ksize', default=31, type=int)
    p.add_argument('-s', '--savename', default=None, type=str)
    p.add_argument('--traverse-directory', action='store_true')
    p.add_argument('--scaled', default=None, type=int)
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
    if args.traverse_directory:
        inp_files = list(traverse_find_sigs(args.sigs))
    else:
        inp_files = list(args.sigs)

    print('loading signatures & traversing hashes')
    bad_input = 0
    for n, filename in enumerate(inp_files):
        if n % 100 == 0:
            print('... loading file #', n, 'of', len(inp_files), end='\r')

        try:
            sig = sourmash_lib.signature.load_one_signature(filename,
                                                      select_ksize=args.ksize)
        except (FileNotFoundError, ValueError):
            if not args.traverse_directory:
                raise

            bad_input += 1
            continue

        acc = sig.name().split(' ')[0]   # first part of sequence name
        acc = acc.split('.')[0]          # get acc w/o version

        taxid = taxfoo.get_taxid(acc)
        if taxid == None:
            continue

        if args.scaled:
            sig.minhash = sig.minhash.downsample_scaled(args.scaled)

        mins = sig.minhash.get_mins()

        for m in mins:
            hashval_to_taxids[m].add(taxid)
    print('\n...done')
    if bad_input:
        print('failed to load {} of {} files found'.format(bad_input,
                                                           len(inp_files)))

    with open(args.savename + '.hashvals', 'wb') as hashval_fp:
        dump(hashval_to_taxids, hashval_fp)

    ####

    print('traversing tags and finding last-common-ancestor for {} tags'.format(len(hashval_to_taxids)))

    hashval_to_lca = {}
    found_root = 0
    empty_set = 0

    # find the LCA for each hashval and store.
    for n, (hashval, taxid_set) in enumerate(hashval_to_taxids.items()):
        if n % 10000 == 0:
            print('...', n, end='\r')

        # find associated least-common-ancestors.
        lca = taxfoo.find_lca(taxid_set)

        if lca == 1:
            if taxid_set:
                found_root += 1
            else:
                empty_set += 1
            continue

        # save!!
        hashval_to_lca[hashval] = lca
    print('\ndone')

    if found_root:
        print('found root {} times'.format(found_root))
    if empty_set:
        print('found empty set {} times'.format(empty_set))

    print('saving to', args.savename)
    dump(hashval_to_lca, open(args.savename, 'wb'))


if __name__ == '__main__':
    main()
