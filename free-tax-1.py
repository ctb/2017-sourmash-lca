#! /usr/bin/env python
"""
Load a genbank-free lineage, anchor with genbank.
"""
import sys
import argparse
import csv
from collections import defaultdict, Counter
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
    p.add_argument('revindex')
    p.add_argument('lca_db')
    p.add_argument('--scaled', default=10000, type=float)
    p.add_argument('-d', '--debug', action='store_true')
    args = p.parse_args()

    if args.debug:
        global _print_debug
        _print_debug = True

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

        assignments[ident] = tuple(lineage)

    ## load revindex
    print('loading reverse index:', args.revindex, file=sys.stderr)
    custom_bins_ri = revindex_utils.HashvalRevindex(args.revindex)

    # load the signatures associated with each revindex.
    print('loading signatures for custom genomes...', file=sys.stderr)
    sigids_to_sig = {}
    ksizes = set()

    n = 0
    total_n = len(custom_bins_ri.sigid_to_siginfo)
    for sigid, (filename, md5) in custom_bins_ri.sigid_to_siginfo.items():
        n += 1
        print(u'\r\033[K', end=u'', file=sys.stderr)
        print('... loading from {} ({} of {})'.format(filename, n, total_n), end='\r',file=sys.stderr)
        sig = revindex_utils.get_sourmash_signature(filename, md5)

        # downsample to specified scaled; this has the side effect of
        # making sure they're all at the same scaled value!
        sig.minhash = sig.minhash.downsample_scaled(args.scaled)

        # check if there's more than one ksize...
        ksizes.add(sig.minhash.ksize)
        if len(ksizes) > 1:
            raise Exception('ERROR: multiple ksizes!? {}'.format(str(ksizes)))
            
        if sig.name() in assignments:
            sigids_to_sig[sigid] = sig
        else:
            debug('no assignment:', sig.name())

    if not len(sigids_to_sig):
        raise Exception('no custom genomes with assignments!?')

    print('...found {} custom genomes in revindex with assignments!!'.format(len(sigids_to_sig)), file=sys.stderr)

    ksize = ksizes.pop()
    scaled = int(args.scaled)

    ## now, connect the dots: hashvals to custom classifications
    hashval_to_custom = defaultdict(list)
    for hashval, sigids in custom_bins_ri.hashval_to_sigids.items():
        for sigid in sigids:
            sig = sigids_to_sig.get(sigid, None)
            if sig:
                assignment = assignments[sig.name()]
                hashval_to_custom[hashval].append(assignment)

    ## clean up with some indirection
    next_lineage_index = 0
    lineage_dict = {}
    hashval_to_custom_idx = defaultdict(list)
    
    for hashval, assignments in hashval_to_custom.items():
        for lineage_tuple in assignments:
            idx = lineage_dict.get(lineage_tuple)
            if idx is None:
                idx = next_lineage_index
                next_lineage_index += 1

                lineage_dict[idx] = lineage_tuple

            hashval_to_custom_idx[hashval].append(idx)

    # whew! done!! we can now go from a hashval to a custom assignment!!
    print('saving to LCA DB v2: {}'.format(args.lca_db))
    with open(args.lca_db, 'wt') as fp:
        save_d = dict(ksize=ksize, scaled=scaled, lineages=lineage_dict,
                      hashval_assignments=hashval_to_custom_idx)
        json.dump([save_d], fp)


if __name__ == '__main__':
    sys.exit(main())
