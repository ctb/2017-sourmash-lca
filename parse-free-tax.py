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
    p.add_argument('names_dmp')
    p.add_argument('nodes_dmp')
    p.add_argument('csv')
    #p.add_argument('-v', '--verbose', action='store_true')
    args = p.parse_args()

    # ok, read in all the tax info.
    print('loading NCBI tax names/nodes', file=sys.stderr)
    taxfoo = ncbi_taxdump_utils.NCBI_TaxonomyFoo()
    taxfoo.load_names_dmp(args.names_dmp, False)
    taxfoo.load_nodes_dmp(args.nodes_dmp, False)

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

        #pprint.pprint(triples_info)

    print('## {} confusing lineages --'.format(len(confusing_lineages)))

    n = 0
    for (csv_str, ncbi_str), ident_list in confusing_lineages.items():
        n += 1
        print('confusing lineage #{}'.format(n))
        print('\tCSV: ', csv_str)
        print('\tNCBI:', ncbi_str)
        print('({} rows in spreadsheet)'.format(len(ident_list)))
        print('---')

    print('\n## {} incompatible lineages --'.format(len(confusing_lineages)))

    n = 0
    for (csv_str, ncbi_str), ident_list in incompatible_lineages.items():
        n += 1
        print('incompatible lineage #{}'.format(n))
        print('\tCSV: ', csv_str)
        print('\tNCBI:', ncbi_str)
        print('({} rows in spreadsheet)'.format(len(ident_list)))
        print('---')

    ## next phase: collapse lineages etc.


if __name__ == '__main__':
    sys.exit(main())
