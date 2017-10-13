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

   kraken/extract.py foobar.lca genbank/*.csv.gz nodes.dmp ecoli_many_sigs/ecoli-*.sig --lca-json=db.lca.json

The E. coli signatures used in the command above can be downloaded like so:

   curl -O -L https://github.com/dib-lab/sourmash/raw/master/data/eschericia-sigs.tar.gz

'nodes.dmp' comes from 'taxdump.tar.gz' --

    curl -O -L ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
    tar xzf taxdump.tar.gz nodes.dmp names.dmp
"""

import sys, os
import sourmash_lib, sourmash_lib.signature
import argparse
from pickle import dump, load
from collections import defaultdict

import lca_json
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
    p.add_argument('lca_output')
    p.add_argument('genbank_csv')
    p.add_argument('nodes_dmp')
    p.add_argument('sigs', nargs='+')
    p.add_argument('-k', '--ksize', default=31, type=int)
    p.add_argument('--scaled', default=10000, type=int)

    p.add_argument('--traverse-directory', action='store_true')

    p.add_argument('-s', '--save-hashvals', action='store_true')
    p.add_argument('-l', '--load-hashvals', action='store_true')

    p.add_argument('--lca-json')
    p.add_argument('--names-dmp', default='')
    args = p.parse_args()

    taxfoo = NCBI_TaxonomyFoo()

    # load the accessions->taxid info
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

    if args.load_hashvals:
        with open(args.lca_output + '.hashvals', 'rb') as hashval_fp:
            print('loading hashvals dict per -l/--load-hashvals...')
            hashval_to_taxids = load(hashval_fp)
            print('loaded {} hashvals'.format(len(hashval_to_taxids)))
    else:

        print('loading signatures & traversing hashes')
        bad_input = 0
        for n, filename in enumerate(inp_files):
            if n % 100 == 0:
                print('... loading file #', n, 'of', len(inp_files), end='\r')

            try:
                sig = sourmash_lib.signature.load_one_signature(filename,
                                                          ksize=args.ksize)
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

            sig.minhash = sig.minhash.downsample_scaled(args.scaled)

            mins = sig.minhash.get_mins()

            for m in mins:
                hashval_to_taxids[m].add(taxid)
        print('\n...done')
        if bad_input:
            print('failed to load {} of {} files found'.format(bad_input,
                                                               len(inp_files)))

        if args.save_hashvals:
            with open(args.lca_output + '.hashvals', 'wb') as hashval_fp:
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

    print('saving to', args.lca_output)
    with lca_json.xopen(args.lca_output, 'wb') as lca_fp:
        dump(hashval_to_lca, lca_fp)

    # update LCA DB JSON file if provided
    if args.lca_json:
        lca_db = lca_json.LCA_Database()
        if os.path.exists(args.lca_json):
            print('loading LCA JSON file:', args.lca_json)
            lca_db.load(args.lca_json)

        prefix = os.path.dirname(args.lca_json) + '/'

        lca_output = args.lca_output
        if lca_output.startswith(prefix):
            lca_output = lca_output[len(prefix):]

        nodes_dmp = args.nodes_dmp
        if nodes_dmp.startswith(prefix):
            nodes_dmp = nodes_dmp[len(prefix):]

        names_dmp = args.names_dmp
        if names_dmp:
            if names_dmp.startswith(prefix):
                names_dmp = names_dmp[len(prefix):]
        else:
            names_dmp = nodes_dmp.replace('nodes', 'names')

        lca_db.add_db(args.ksize, args.scaled, lca_output, nodes_dmp, names_dmp)
        print('saving LCA JSON file:', args.lca_json)
        lca_db.save(args.lca_json)


if __name__ == '__main__':
    main()
