#! /usr/bin/env python
"""
Load a genbank-free lineage, anchor with genbank.
"""
import sys
import argparse
import csv
from collections import defaultdict, Counter, OrderedDict
import itertools
import pprint
import sourmash_lib
import json

sys.path.insert(0, '../2017-sourmash-revindex')
import revindex_utils

taxlist = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus',
           'species']
null_names = set(['[Blank]', 'na', 'null'])


_print_debug = False
def debug(*args):
    if _print_debug:
        print(*args)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('csv')
    p.add_argument('lca_db_out')
    p.add_argument('genome_sigs', nargs='+')
    p.add_argument('--scaled', default=10000, type=float)
    p.add_argument('-k', '--ksize', default=31, type=int)
    p.add_argument('-d', '--debug', action='store_true')
    args = p.parse_args()

    if args.debug:
        global _print_debug
        _print_debug = True

    scaled = int(args.scaled)
    ksize = int(args.ksize)

    ### parse spreadsheet
    r = csv.reader(open(args.csv, 'rt'))
    row_headers = ['identifier'] + taxlist

    print('examining spreadsheet headers...', file=sys.stderr)
    first_row = next(iter(r))

    n_disagree = 0
    for (column, value) in zip(row_headers, first_row):
        if column.lower() != value.lower():
            print('** assuming {} == {} in spreadsheet'.format(column, value),
                  file=sys.stderr)
            n_disagree += 1
            if n_disagree > 2:
                print('whoa, too many assumptions. are the headers right?',
                      file=sys.stderr)
                sys.exit(-1)

    assignments = {}
    for row in r:
        lineage = list(zip(row_headers, row))

        ident = lineage[0][1]
        lineage = lineage[1:]

        # clean lineage of null names
        lineage = [(a,b) for (a,b) in lineage if b not in null_names]

        # convert to dictionary

        assignments[ident] = tuple(lineage)

    ## clean up with some indirection: convert lineages to numbers
    next_lineage_index = 0
    lineage_dict = {}

    assignments_idx = {}
    for (ident, lineage_tuple) in assignments.items():
        idx = lineage_dict.get(lineage_tuple)
        if idx is None:
            idx = next_lineage_index
            next_lineage_index += 1

            lineage_dict[idx] = lineage_tuple

        assignments_idx[ident] = idx

    # load signatures, construct index of hashvals to lineages
    hashval_to_lineage = defaultdict(list)
    md5_to_lineage = {}

    n = 0
    total_n = len(args.genome_sigs)
    for filename in args.genome_sigs:
        for sig in sourmash_lib.load_signatures(filename, ksize=args.ksize):
            print(u'\r\033[K', end=u'', file=sys.stderr)
            print('... loading signature {} (file {} of {})'.format(sig.name(), n, total_n), end='\r',file=sys.stderr)

            # is this one for which we have a lineage assigned?
            lineage_idx = assignments_idx.get(sig.name())
            if lineage_idx:
                # downsample to specified scaled; this has the side effect of
                # making sure they're all at the same scaled value!
                sig.minhash = sig.minhash.downsample_scaled(args.scaled)

                # connect hashvals to lineage
                for hashval in sig.minhash.get_mins():
                    hashval_to_lineage[hashval].append(lineage_idx)

                # store md5 -> lineage too
                md5_to_lineage[sig.md5sum()] = lineage_idx

    print(u'\r\033[K', end=u'', file=sys.stderr)
    print('...found {} genomes with lineage assignments!!'.format(len(md5_to_lineage)), file=sys.stderr)

    # now, save!
    print('saving to LCA DB v2: {}'.format(args.lca_db_out))
    with open(args.lca_db_out, 'wt') as fp:
        save_d = OrderedDict()
        save_d['version'] = '1.0'
        save_d['type'] = 'sourmash_lca'
        save_d['ksize'] = ksize
        save_d['scaled'] = scaled
        # convert lineage internals from tuples to dictionaries
        save_d['lineages'] = OrderedDict([ (k, OrderedDict(v)) \
                                           for k, v in lineage_dict.items() ])
        save_d['hashval_assignments'] = hashval_to_lineage
        save_d['signatures_to_lineage'] = md5_to_lineage
        json.dump(save_d, fp)


if __name__ == '__main__':
    sys.exit(main())
