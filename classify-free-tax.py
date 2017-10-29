#! /usr/bin/env python
"""
Load a genbank-free lineage, anchor with genbank.
"""
import sys
import argparse
import csv
import traceback
import ncbi_taxdump_utils
from collections import defaultdict
import pprint
import sourmash_lib
import lca_json                      # from github.com/ctb/2017-sourmash-lca

LCA_DBs = ['db/genbank.lca.json']
SCALED=10000
THRESHOLD=5                               # how many counts of a taxid at min

sys.path.insert(0, '../2017-sourmash-revindex')
import revindex_utils

class TaxidNotFound(Exception):
    pass


taxlist = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus',
           'species']
null_names = set(['[Blank]', 'na', 'null'])


def get_taxids_for_name(taxfoo, names_to_taxids, srank, sname):

    # what taxids does the query name have?
    tid_set = names_to_taxids.get(sname, set())
    tid_list = list(tid_set)

    # none? that's a problem, quit out.
    if not tid_list:
        raise TaxidNotFound(sname)

    # collect taxids at the right rank
    taxid_at_rank = []
    for taxid in tid_list:
        rank = taxfoo.get_taxid_rank(taxid)
        if rank == srank:
            assert taxfoo.get_taxid_name(taxid) == sname
            taxid_at_rank.append(taxid)

    # there should only be one 
    if len(taxid_at_rank) == 1:
        return taxid_at_rank[0]

    # @CTB need to do something more here.

    return -1


def get_lca_taxid_for_lineage(taxfoo, names_to_taxids, lineage):
    """
    Given taxfoo and a list of lineage pairs (rank, identifier), find the least
    common ancestor in the lineage, and return that taxid with the rest of
    the lineage pairs.
    """
    lineage = list(lineage)               # make a copy
    while lineage:
        (rank, name) = lineage.pop(0)
        try:
            taxid = get_taxids_for_name(taxfoo, names_to_taxids, rank, name)
            if taxid == -1:
                raise TaxidNotFound
            last_taxid = taxid

            assert taxfoo.get_taxid_rank(taxid) == rank
            assert taxfoo.get_taxid_name(taxid) == name
        except TaxidNotFound:
            lineage.insert(0, (rank, name))   # add back in!
            break

    return last_taxid, lineage


def get_lowest_taxid_for_lineage(taxfoo, names_to_taxids, lineage):
    """
    Given taxfoo and a list of lineage pairs (rank, identifier), find
    the lowest rank that has a match in NCBI lineage, and return that
    taxid with the rest of the lineage pairs.
    """
    lineage = list(lineage)               # make a copy
    remainder = []
    while 1:
        (rank, ident) = lineage.pop()     # pop from end
        try:
            taxid = get_taxids_for_name(taxfoo, names_to_taxids, rank, ident)
            if taxid == -1:
                raise TaxidNotFound
        except TaxidNotFound:
            remainder.append((rank, ident))
            continue                      # super borked logic

        break                             # super borked logic

    return taxid, list(reversed(remainder))


def main():
    p = argparse.ArgumentParser()
    p.add_argument('csv')
    p.add_argument('revindex')
    p.add_argument('siglist', nargs='+')
    p.add_argument('--lca', nargs='+', default=LCA_DBs)
    p.add_argument('-k', '--ksize', default=31, type=int)
    #p.add_argument('-v', '--verbose', action='store_true')
    args = p.parse_args()

    ## load LCA databases
    lca_db_list = []
    for lca_filename in args.lca:
        print('loading LCA database from {}'.format(lca_filename))
        lca_db = lca_json.LCA_Database(lca_filename)
        taxfoo, hashval_to_lca, _ = lca_db.get_database(args.ksize, SCALED)
        lca_db_list.append((taxfoo, hashval_to_lca))
    
    # reverse index names -> taxids
    names_to_taxids = defaultdict(set)
    for taxid, (name, _, _) in taxfoo.taxid_to_names.items():
        names_to_taxids[name].add(taxid)

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

    confusing_lineages = defaultdict(list)
    incompatible_lineages = defaultdict(list)
    assignments = {}
    for row in r:
        lineage = list(zip(row_headers, row))

        ident = lineage[0][1]
        lineage = lineage[1:]

        # clean lineage of null names
        lineage = [(a,b) for (a,b) in lineage if b not in null_names]

        # ok, find the least-common-ancestor taxid...
        taxid, rest = get_lca_taxid_for_lineage(taxfoo, names_to_taxids,
                                                lineage)

        # and find the *lowest* identifiable ancestor taxid, just to see
        # if there are confusing lineages.
        lowest_taxid, lowest_rest = \
          get_lowest_taxid_for_lineage(taxfoo, names_to_taxids, lineage)

        # do they match? if not, report.
        if lowest_taxid != taxid:
            lowest_lineage = taxfoo.get_lineage(lowest_taxid, taxlist)
            lowest_str = ', '.join(lowest_lineage)

            # find last matching, in case different classification levels.
            match_lineage = [ b for (a, b) in lineage ]
            end = match_lineage.index(lowest_lineage[-1])
            assert end >= 0
            match_lineage = match_lineage[:end + 1]
            match_str = ', '.join(match_lineage)

            confusing_lineages[(match_str, lowest_str)].append(ident)
            continue

        # check! NCBI lineage should be lineage of taxid + rest
        ncbi_lineage = taxfoo.get_lineage(taxid, taxlist)
        assert len(ncbi_lineage)
        reconstructed = ncbi_lineage + [ b for (a,b) in rest ]

        # ...make a comparable lineage from the CSV line...
        csv_lineage = [ b for (a, b) in lineage ]

        # are NCBI-rooted and CSV lineages the same?? if not, report.
        if csv_lineage != reconstructed:
            csv_str = ", ".join(csv_lineage[:len(ncbi_lineage)])
            ncbi_str = ", ".join(ncbi_lineage)
            incompatible_lineages[(csv_str, ncbi_str)].append(ident)
            continue

        # all is well if we've reached this point! We've got NCBI-rooted
        # taxonomies and now we need to record. next:
        #
        # build a set of triples: (rank, name, taxid), where taxid can
        # be None.

        lineage_taxids = taxfoo.get_lineage_as_taxids(taxid)
        triples_info = []
        for taxid in lineage_taxids:
            name = taxfoo.get_taxid_name(taxid)
            rank = taxfoo.get_taxid_rank(taxid)

            if rank in taxlist:
                triples_info.append((rank, name, taxid))

        for (rank, name) in rest:
            assert rank in taxlist
            triples_info.append((rank, name, None))

        assignments[ident] = triples_info

        #pprint.pprint(triples_info)

    print('{} weird lineages; ignoring for now.'.format(len(confusing_lineages) + len(incompatible_lineages)))

    ## next phase: collapse lineages etc.

    ## load revindex
    print('loading reverse index:', args.revindex)
    custom_bins_ri = revindex_utils.HashvalRevindex(args.revindex)

    # load the signatures associated with each revindex.
    print('loading signatures for custom genomes...')
    custom_sigs = {}
    sigids_to_sig = {}
    for sigid, (filename, md5) in custom_bins_ri.sigid_to_siginfo.items():
        sig = revindex_utils.get_sourmash_signature(filename, md5)
        if sig.name() in assignments:
            custom_sigs[md5] = sig
            sigids_to_sig[sigid] = sig

    # figure out what ksize we're talking about here! (this should
    # probably be stored on the revindex...)
    random_db_sig = next(iter(sigids_to_sig.values()))
    ksize = random_db_sig.minhash.ksize

    print('...found {} signatures that also have assignments!!'.format(len(custom_sigs)))

    ## now, connect the dots: hashvals to custom classifications
    hashval_to_custom = {}
    for hashval, sigids in custom_bins_ri.hashval_to_sigids.items():
        for sigid in sigids:
            sig = sigids_to_sig.get(sigid, None)
            if sig:
                assignment = assignments[sig.name()]
                hashval_to_custom[hashval] = assignment

    # whew! done!! we can now go from a hashval to a custom assignment!!
    for query_filename in args.siglist:
        for query_sig in sourmash_lib.load_signatures(query_filename,
                                                      ksize=ksize):
            print('Query!', query_sig.name())
            names = set()
            for hashval in query_sig.minhash.get_mins():
                assignment = hashval_to_custom.get(hashval)
                if assignment:
                    lineage_str = "; ".join([ b for (a,b,c) in assignment ])
                    names.add(lineage_str)
            print("\t{}".format("\n\t".join(names)))


if __name__ == '__main__':
    sys.exit(main())
