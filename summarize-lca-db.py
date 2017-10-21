#! /usr/bin/env python
"""
Collect LCA rank statistics across several LCA databases.
"""
import sys
import argparse
import collections
import lca_json
import csv

from ncbi_taxdump_utils import want_taxonomy


def summarize_lca_db(taxfoo, hashval_to_lca):
    rank_counts = collections.defaultdict(int)

    print('iterating over {} hash vals'.format(len(hashval_to_lca)))

    n = 0
    for hashval, lca in hashval_to_lca.items():
        n += 1
        if n and n % 100000 == 0:
            print('... {}'.format(n), end='\r')

        rank = taxfoo.get_taxid_rank(lca)

        # pull rank back to next interesting taxonomic rank:
        while rank not in want_taxonomy:
            if lca == 1 or lca is None:
                break
            
            lca = taxfoo.get_taxid_parent(lca)
            rank = taxfoo.get_taxid_rank(lca)

        if lca and lca != 1:
            rank_counts[rank] += 1

    print('... done! {}'.format(n))

    return rank_counts


def main():
    p = argparse.ArgumentParser()
    p.add_argument('lca_filename')
    p.add_argument('-k', '--ksize-list', default="31", type=str)
    p.add_argument('-o', '--output', type=argparse.FileType('wt'))
    args = p.parse_args()

    lca_db = lca_json.LCA_Database(args.lca_filename)

    ksizes = list(map(int, args.ksize_list.split(',')))

    ksize_to_rank_counts = dict()
    
    for ksize in ksizes:
        #assert ksize not in ksize_to_rank_counts
        taxfoo, hashval_to_lca, scaled = lca_db.get_database(ksize, None)

        rank_counts = summarize_lca_db(taxfoo, hashval_to_lca)
        ksize_to_rank_counts[ksize] = rank_counts

    # this should be enforced by summarize_lca_db(...)
    all_ranks = set()
    for rank_counts in ksize_to_rank_counts.values():
        all_ranks.update(rank_counts.keys())

    assert all_ranks - set(want_taxonomy) == set()

    if args.output:
        w = csv.writer(args.output)
    else:
        w = csv.writer(sys.stdout)

    w.writerow(['rank'] + ksizes)
    for rank in want_taxonomy:
        count_list = [rank]
        for ksize in ksizes:
            rank_counts = ksize_to_rank_counts[ksize]
            count = rank_counts.get(rank, 0)
            count_list.append(str(count))

        w.writerow(count_list)


if __name__ == '__main__':
    sys.exit(main())
