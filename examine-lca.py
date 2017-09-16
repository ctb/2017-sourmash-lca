#! /usr/bin/env python
"""
...
"""
import os
import argparse
import collections
import sys
import csv
from pickle import load, dump

import sys
sys.path.insert(0, '/Users/t/dev/2017-sourmash-revindex')
import revindex_utils
from revindex_utils import HashvalRevindex

import sourmash_lib.signature

import lca_json                      # from github.com/ctb/2017-sourmash-lca

LCA_DBs = []
SCALED=10000

want_taxonomy = ['superkingdom', 'phylum', 'order', 'class', 'family', 'genus', 'species']


def traverse_find_sigs(dirnames):
    """
    Find all the filenames ending with .sig under given directories.
    """
    for dirname in dirnames:
        for root, dirs, files in os.walk(dirname):
            for name in files:
                if name.endswith('.sig'):
                    fullname = os.path.join(root, name)
                    yield fullname


def load_all_signatures(dirname, ksize):
    """
    Load all signatures under given dirname with given ksize, return
    dictionary d[name] -> signature.

    Because this can be slow for many hundreds of signatures, cache dict
    using pickle.
    """
    sigd = {}
    filename_d = {}

    pickle_cache = dirname.rstrip('/') + '.pickle'
    if os.path.exists(pickle_cache):
        print('loading from cache:', pickle_cache)
        with open(pickle_cache, 'rb') as fp:
            filename_d, sigd = load(fp)

    loaded_new = False

    n = 0
    for filename in traverse_find_sigs([dirname]):
        if filename not in filename_d:
            loaded_new = True
            sig = sourmash_lib.signature.load_one_signature(filename,
                                                            select_ksize=ksize)
            filename_d[filename] = 1
            sigd[sig.name()] = sig

        n += 1

    if loaded_new:
        print('saving to cache:', pickle_cache)
        with open(pickle_cache, 'wb') as fp:
            dump((filename_d, sigd), fp)

    return sigd


def main():
    p = argparse.ArgumentParser()
    p.add_argument('-k', '--ksize', default=31, type=int)
    p.add_argument('lca')
    p.add_argument('revindex')
    p.add_argument('accessions_csv')
    args = p.parse_args()

    # load the LCA databases from the JSON file(s)
    print('loading LCA database from {}'.format(args.lca))
    lca_db = lca_json.LCA_Database(args.lca)
    taxfoo, hashval_to_lca, _ = lca_db.get_database(args.ksize, SCALED)

    print('loading revindex:', args.revindex)
    revidx = HashvalRevindex(args.revindex)
    print('...loaded.')

    # load classification CSV
    print('loading classifications:', args.accessions_csv)
    taxfoo.load_accessions_csv(args.accessions_csv)
    print('...loaded.')

    ###

    # track number of classifications at various rank levels
    classified_at = collections.defaultdict(int)
    classified_samples = collections.defaultdict(int)

    hashval_troubles = collections.defaultdict(set)
    for hashval, lca in hashval_to_lca.items():
        rank = taxfoo.get_taxid_rank(lca)
        classified_at[rank] += 1

        if rank == 'superkingdom' and 0:
            n_sigids = len(revidx.hashval_to_sigids[hashval])
            classified_samples[n_sigids] += 1
            if n_sigids >= 4:
                for sigid in revidx.hashval_to_sigids[hashval]:
                    siginfo = revidx.sigid_to_siginfo[sigid]
                    hashval_troubles[hashval].add(siginfo)

    for hashval, siginfo_set in hashval_troubles.items():
        break
        print('getting {} sigs for {}'.format(len(siginfo_set), hashval))
        siglist = []
        for (filename, md5) in siginfo_set:
            sig = revindex_utils.get_sourmash_signature(filename, md5)
            siglist.append(sig)

        for sig in siglist:
            acc = sig.name().split()[0]
            taxid = taxfoo.get_taxid(acc)
            if taxid:
                print('\t', ";".join(taxfoo.get_lineage(taxid)), taxfoo.get_taxid_rank(taxid))


    print('')
    for rank in want_taxonomy:
        if classified_at.get(rank):
            print('\t{}: {}'.format(rank, classified_at.get(rank, 0)))

    print(classified_samples)


if __name__ == '__main__':
    main()
