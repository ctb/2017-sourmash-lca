#! /usr/bin/env python
"""
Load a genbank-free lineage, anchor with genbank.
"""
import sys
import argparse
import csv
import traceback
import ncbi_taxdump_utils
from collections import defaultdict, Counter
import itertools
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


_print_debug = False
def debug(*args):
    if _print_debug:
        print(*args)


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


def build_tree(assignments, initial=None):
    """
    Builds a tree of dictionaries from lists of (rank, name) tuples
    in 'assignments'.  This tree can then be used to find least common
    ancestor agreements/confusion.
    """
    if initial is None:
        tree = {}
    else:
        tree = initial

    for assignment in assignments:
        node = tree

        for rank, name in assignment:
            child = node.get((rank, name), {})
            node[(rank, name)] = child

            # shift -> down in tree
            node = child

    return tree


def test_build_tree():
    tree = build_tree([[('rank1', 'name1'), ('rank2', 'name2')]])
    assert tree == { ('rank1', 'name1'): { ('rank2', 'name2') : {}} }


def test_build_tree_2():
    tree = build_tree([[('rank1', 'name1'), ('rank2', 'name2a')],
                       [('rank1', 'name1'), ('rank2', 'name2b')],
                      ])

    assert tree == { ('rank1', 'name1'): { ('rank2', 'name2a') : {},
                                           ('rank2', 'name2b') : {}} }


def test_build_tree_3():
    tree = build_tree([[('rank1', 'name1'), ('rank2', 'name2a')],
                      ])

    tree = build_tree([[('rank1', 'name1'), ('rank2', 'name2b')],
                      ], tree)

    assert tree == { ('rank1', 'name1'): { ('rank2', 'name2a') : {},
                                           ('rank2', 'name2b') : {}} }


def find_lca(tree):
    """
    Given a tree produced by 'find_tree', find the first node with multiple
    children, OR the only leaf in the tree.  Return ((rank, name), reason),
    where 'reason' is the number of children of the returned node, i.e.e
    0 if it's a leaf and > 1 if it's an internal node.
    """

    node = tree
    cur = ('root', 'root')
    while 1:
        if len(node) == 1:                # descend to only child
            cur = next(iter(node.keys()))
            node = node[cur]
        elif len(node) == 0:              # at leaf; end
            return cur, 0
        else:                             # len(node) > 1 => confusion!!
            return cur, len(node)


def test_find_lca():
    tree = build_tree([[('rank1', 'name1'), ('rank2', 'name2')]])
    lca = find_lca(tree)

    assert lca == (('rank2', 'name2'), 0)


def test_find_lca_2():
    tree = build_tree([[('rank1', 'name1'), ('rank2', 'name2a')],
                       [('rank1', 'name1'), ('rank2', 'name2b')],
                      ])
    lca = find_lca(tree)

    assert lca == (('rank1', 'name1'), 2)


def build_reverse_tree(assignments, initial=None):
    """
    Builds a child -> parent dictionary (a reverse DAG) from lists of
    (rank, name) tuples in 'assignments'.
    """
    if initial is None:
        parents = {}
    else:
        parents = initial

    for assignment in assignments:
        last_node = ('root', 'root')
        for rank, name in assignment:
            parents[(rank, name)] = last_node
            last_node = (rank, name)

    return parents


def test_build_reverse_tree():
    parents = build_reverse_tree([[('rank1', 'name1'), ('rank2', 'name2')]])

    print(parents)
    assert parents == { ('rank2', 'name2'): ('rank1', 'name1'),
                        ('rank1', 'name1'): ('root', 'root') }


def test_build_reverse_tree_2():
    parents = build_reverse_tree([[('rank1', 'name1'), ('rank2', 'name2a')],
                                 [('rank1', 'name1'), ('rank2', 'name2b')],
                                 ])

    assert parents == { ('rank2', 'name2a'): ('rank1', 'name1'),
                        ('rank2', 'name2b'): ('rank1', 'name1'),
                        ('rank1', 'name1'): ('root', 'root') }


def test_build_reverse_tree_3():
    parents = build_reverse_tree([[('rank1', 'name1'), ('rank2', 'name2a')],
                                 ])
    parents = build_reverse_tree([[('rank1', 'name1'), ('rank2', 'name2b')],
                                 ], parents)

    assert parents == { ('rank2', 'name2a'): ('rank1', 'name1'),
                        ('rank2', 'name2b'): ('rank1', 'name1'),
                        ('rank1', 'name1'): ('root', 'root') }


def main():
    p = argparse.ArgumentParser()
    p.add_argument('csv')
    p.add_argument('revindex')
    p.add_argument('siglist', nargs='+')
    p.add_argument('--lca', nargs='+', default=LCA_DBs)
    p.add_argument('-k', '--ksize', default=31, type=int)
    p.add_argument('-o', '--output', type=argparse.FileType('wt'),
                   help='output CSV to this file instead of stdout')
    #p.add_argument('-v', '--verbose', action='store_true')
    p.add_argument('-d', '--debug', action='store_true')
    args = p.parse_args()

    if args.debug:
        global _print_debug
        _print_debug = True

    ## load LCA databases
    lca_db_list = []
    for lca_filename in args.lca:
        print('loading LCA database from {}'.format(lca_filename),
              file=sys.stderr)
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

        # all is well if we've reached this point! We've got NCBI-rooted
        # taxonomies and now we need to record. next:
        #
        # build a set of triples: (rank, name, taxid), where taxid can
        # be None.

        lineage_taxids = taxfoo.get_lineage_as_taxids(taxid)
        tuples_info = []
        for taxid in lineage_taxids:
            name = taxfoo.get_taxid_name(taxid)
            rank = taxfoo.get_taxid_rank(taxid)

            if rank in taxlist:
                tuples_info.append((rank, name))

        for (rank, name) in rest:
            assert rank in taxlist
            tuples_info.append((rank, name))

        assignments[ident] = tuples_info

    print("{} weird lineages that maybe don't match with NCBI.".format(len(confusing_lineages) + len(incompatible_lineages)), file=sys.stderr)

    ## next phase: collapse lineages etc.

    ## load revindex
    print('loading reverse index:', args.revindex, file=sys.stderr)
    custom_bins_ri = revindex_utils.HashvalRevindex(args.revindex)

    # load the signatures associated with each revindex.
    print('loading signatures for custom genomes...', file=sys.stderr)
    sigids_to_sig = {}
    for sigid, (filename, md5) in custom_bins_ri.sigid_to_siginfo.items():
        sig = revindex_utils.get_sourmash_signature(filename, md5)
        if sig.name() in assignments:
            sigids_to_sig[sigid] = sig
        else:
            debug('no assignment:', sig.name())

    # figure out what ksize we're talking about here! (this should
    # probably be stored on the revindex...)
    random_db_sig = next(iter(sigids_to_sig.values()))
    ksize = random_db_sig.minhash.ksize

    print('...found {} custom genomes that also have assignments!!'.format(len(sigids_to_sig)), file=sys.stderr)

    ## now, connect the dots: hashvals to custom classifications
    hashval_to_custom = defaultdict(list)
    for hashval, sigids in custom_bins_ri.hashval_to_sigids.items():
        for sigid in sigids:
            sig = sigids_to_sig.get(sigid, None)
            if sig:
                assignment = assignments[sig.name()]
                hashval_to_custom[hashval].append(assignment)

    # whew! done!! we can now go from a hashval to a custom assignment!!

    # for each query, gather all the matches in both custom and NCBI, then
    # classify.
    csvfp = csv.writer(sys.stdout)
    if args.output:
        print("outputting classifications to '{}'".format(args.output.name))
        csvfp = csv.writer(args.output)
    else:
        print("outputting classifications to stdout")
    csvfp.writerow(['ID'] + taxlist)

    total_count = 0
    for query_filename in args.siglist:
        for query_sig in sourmash_lib.load_signatures(query_filename,
                                                      ksize=ksize):
            print(u'\r\033[K', end=u'', file=sys.stderr)
            print('... classifying {}'.format(query_sig.name()), end='\r',
                  file=sys.stderr)
            debug('classifying', query_sig.name())
            total_count += 1

            these_assignments = defaultdict(list)
            n_custom = 0
            for hashval in query_sig.minhash.get_mins():
                # custom
                assignment = hashval_to_custom.get(hashval, [])
                if assignment:
                    these_assignments[hashval].extend(assignment)
                    n_custom += 1

                # NCBI
                for (this_taxfoo, hashval_to_lca) in lca_db_list:
                    hashval_lca = hashval_to_lca.get(hashval)
                    if hashval_lca is not None and hashval_lca != 1:
                        lineage = this_taxfoo.get_lineage_as_dict(hashval_lca,
                                                                  taxlist)

                        tuple_info = []
                        for rank in taxlist:
                            if rank not in lineage:
                                break
                            tuple_info.append((rank, lineage[rank]))
                        these_assignments[hashval_lca].append(tuple_info)

            check_counts = Counter()
            for tuple_info in these_assignments.values():
                last_tup = tuple(tuple_info[-1])
                check_counts[last_tup] += 1

            debug('n custom hashvals:', n_custom)
            debug(pprint.pformat(check_counts.most_common()))

            # now convert to trees -> do LCA & counts
            counts = Counter()
            parents = {}
            for hashval in these_assignments:

                # for each list of tuple_info [(rank, name), ...] build
                # a tree that lets us discover least-common-ancestor.
                tuple_info = these_assignments[hashval]
                tree = build_tree(tuple_info)

                # also update a tree that we can ascend from leaves -> parents
                # for all assignments for all hashvals
                parents = build_reverse_tree(tuple_info, parents)

                # now find either a leaf or the first node with multiple
                # children; that's our least-common-ancestor node.
                lca, reason = find_lca(tree)
                counts[lca] += 1

            # ok, we now have the LCAs for each hashval, and their number
            # of counts. Now sum across "significant" LCAs - those above
            # threshold.

            tree = {}
            tree_counts = defaultdict(int)

            debug(pprint.pformat(counts.most_common()))

            n = 0
            for lca, count in counts.most_common():
                if count < THRESHOLD:
                    break

                n += 1

                xx = []
                parent = lca
                while parent:
                    xx.insert(0, parent)
                    tree_counts[parent] += count
                    parent = parents.get(parent)
                debug(n, count, xx[1:])

                # update tree with this set of assignments
                build_tree([xx], tree)

            if n > 1:
                debug('XXX', n)

            # now find LCA? or whatever.
            lca, reason = find_lca(tree)
            if reason == 0:               # leaf node
                debug('END', lca)
            else:                         # internal node
                debug('MULTI', lca)

            # backtrack to full lineage via parents
            lineage = []
            parent = lca
            while parent != ('root', 'root'):
                lineage.insert(0, parent)
                parent = parents.get(parent)

            # output!
            row = [query_sig.name()]
            for taxrank, (rank, name) in itertools.zip_longest(taxlist, lineage, fillvalue=('', '')):
                if rank:
                    assert taxrank == rank
                row.append(name)

            csvfp.writerow(row)

    print(u'\r\033[K', end=u'', file=sys.stderr)
    print('classified {} signatures total'.format(total_count), file=sys.stderr)


if __name__ == '__main__':
    sys.exit(main())
