#! /usr/bin/env python
import argparse
import csv
import traceback
import ncbi_taxdump_utils
from collections import defaultdict


class TaxidNotFound(Exception):
    pass


want_taxonomy = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']


def get_taxids_for_name(taxfoo, names_to_taxids, srank, sname):

    tid_set = names_to_taxids.get(sname, set())
    tid_list = list(tid_set)

    if not tid_list:
        raise TaxidNotFound

    # get all of the lineages for each of the taxids found
    lineages = [taxfoo.get_lineage_as_dict(tid, want_taxonomy) for tid in tid_list]

    # is there a unique lineage? if so, *done*! otherwise:
    if len(lineages) > 1:
        # if non-unique, collect (hopefully unique) one
        # that matches the desired rank
        tid_list2 = []
        lin2 = []
        for (tid, lin) in zip(tid_list, lineages):
            if lin.get(srank) == sname:
                tid_list2.append(tid)
                lin2.append(lin)

        lineages = lin2
        tid_list = tid_list2

        # this should never happen!
        if len(tid_list) > 1:
            assert 0, (tid_list, lineages)
            return -1

    if not lineages:
        return -1

    taxid = tid_list[0]
    return taxid


def main():
    p = argparse.ArgumentParser()
    p.add_argument('names_dmp')
    p.add_argument('nodes_dmp')
    p.add_argument('csv')
    p.add_argument('outcsv')
    p.add_argument('-v', '--verbose', action='store_true')
    args = p.parse_args()

    # ok, read in all the tax info.
    taxfoo = ncbi_taxdump_utils.NCBI_TaxonomyFoo()
    taxfoo.load_names_dmp(args.names_dmp, False)
    taxfoo.load_nodes_dmp(args.nodes_dmp, False)

    # reverse index names -> tasxids
    names_to_taxids = defaultdict(set)
    for taxid, (name, _, _) in taxfoo.taxid_to_names.items():
        names_to_taxids[name].add(taxid)

    # now, grab the spreadsheet & ignore the headers.
    r = csv.reader(open(args.csv, 'rt', encoding='utf8'))
    next(r)

    # track results here:
    results = []

    n = 0
    found = 0

    # for every row, try to figure it out.
    for row in r:
        n += 1
        last_rank = None
        last_name = None

        genome_name, row = row[0], row[1:]

        # walk through the ranks top down:
        ranknames = []
        for rank, name in zip(want_taxonomy, row):
            ranknames.append((rank, name))
            if name == 'na':
                break

        specified_lineage = ";".join([name for (_, name) in ranknames])

        taxid = 1
        while ranknames:
            last_rank, last_name = ranknames.pop()   # get last/most specific
            try:
                taxid = get_taxids_for_name(taxfoo, names_to_taxids, last_rank, last_name)
                found += 1
                break
            except TaxidNotFound:
                if args.verbose:
                    print('\ncan\'t find {} {}, punting to next'.format(last_rank, last_name))
                continue
            except AssertionError:
                #print('\nconfused by {}, skipping.'.format(genome_name))
                break

        if taxid == 1:
            print("\ncouldn't find taxid for {}, {}".format(genome_name, specified_lineage))

        lineage = taxfoo.get_lineage(taxid, want_taxonomy)
        results.append((genome_name, str(taxid), ";".join(lineage)))

        end_of_lineage='...' + (";".join(lineage))[-30:]
        print("\r\033[Kfound! {} of {}: {} - {} - {}".format(found, n, genome_name, taxid, end_of_lineage), end="")

    print('')

    with open(args.outcsv, 'wt', encoding='utf8') as outfp:
        w = csv.writer(outfp)
        for name, taxid, lineage in results:
            w.writerow([name, taxid, lineage])

if __name__ == '__main__':
    main()
