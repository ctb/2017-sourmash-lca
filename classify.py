#! /usr/bin/env python
"""
Emulate the kraken k-mer classification step, but on banded MinHash signatures.

Briefly,

* load in the k-mer-to-lineage database produced by 'extract.py'.
* for every hash in the signature, find the lineage for that k-mer
  (here 'lineage' would be the computed last-common-ancestor from the NCBI
  taxonomy, based on GenBank genomes)
* count & summarize

Usage:

   kraken/classify.py db.lca.json ecoli_many_sigs/ecoli-1.sig

where 'db.lca.json' is created by 'extract.py'.

TODO:

* add classification of FASTA/FASTQ?
"""

import argparse
import collections

import sourmash_lib, sourmash_lib.signature
import lca_json

SCALED=10000                              # should match the LCA compute @CTB

kraken_rank_code = {
    'genus' : 'G',
    'species': 'S',
    'phylum': 'P',
    'class': 'C',
    'order': 'O',
    'family': 'F',
    'genus': 'G',
    'kingdom': 'K',
    'domain': 'D' }


def main():
    p = argparse.ArgumentParser()
    p.add_argument('lca_filename')
    p.add_argument('sigfiles', nargs='+')
    p.add_argument('-k', '--ksize', default=31, type=int)
    args = p.parse_args()

    # load lca info
    lca_db = lca_json.LCA_Database(args.lca_filename)
    taxfoo, hashval_to_lca, scaled = lca_db.get_database(args.ksize, SCALED)
    
    # load signatures
    siglist = []
    print('loading signatures from {} signature files'.format(len(args.sigfiles)))
    for sigfile in args.sigfiles:
        sigs = sourmash_lib.signature.load_signatures(sigfile,
                                                      ksize=args.ksize)
        sigs = list(sigs)
        siglist.extend(sigs)

    print('loaded {} signatures total at k={}'.format(len(siglist), args.ksize))

    # downsample
    print('downsampling to scaled value: {}'.format(scaled))
    for sig in siglist:
        if sig.minhash.scaled < scaled:
            sig.minhash = sig.minhash.downsample_scaled(scaled)

    # now, extract hash values!
    hashvals = collections.defaultdict(int)
    for sig in siglist:
        for hashval in sig.minhash.get_mins():
            hashvals[hashval] += 1

    found = 0
    total = 0
    by_taxid = collections.defaultdict(int)

    # for every hash, get LCA of labels
    for hashval, count in hashvals.items():
        lca = hashval_to_lca.get(hashval)
        total += count

        if lca is None:
            by_taxid[0] += count
            continue

        by_taxid[lca] += count
        found += count

    print('found LCA classifications for', found, 'of', total, 'hashes')
    not_found = total - found

    # now, propogate counts up the taxonomic tree.
    by_taxid_lca = collections.defaultdict(int)
    for taxid, count in by_taxid.items():
        by_taxid_lca[taxid] += count

        parent = taxfoo.child_to_parent.get(taxid)
        while parent != None and parent != 1:
            by_taxid_lca[parent] += count
            parent = taxfoo.child_to_parent.get(parent)

    total_count = sum(by_taxid.values())

    # sort by lineage length
    x = []
    for taxid, count in by_taxid_lca.items():
        x.append((len(taxfoo.get_lineage(taxid)), taxid, count))

    x.sort()

    # ...aaaaaand output.
    print('{}\t{}\t{}\t{}\t{}\t{}'.format('percent', 'below', 'at node',
                                          'code', 'taxid', 'name'))
    for _, taxid, count_below in x:
        percent = round(100 * count_below / total_count, 2)
        count_at = by_taxid[taxid]

        rank = taxfoo.node_to_info.get(taxid)
        if rank:
            rank = rank[0]
            classify_code = kraken_rank_code.get(rank, '-')
        else:
            classify_code = '-'

        name = taxfoo.taxid_to_names.get(taxid)
        if name:
            name = name[0]
        else:
            name = '-'

        print('{}\t{}\t{}\t{}\t{}\t{}'.format(percent, count_below, count_at,
                                              classify_code, taxid, name))

    if not_found:
        classify_code = 'U'
        percent = round(100 * not_found / total_count, 2)
        count_below = not_found
        count_at = not_found
        taxid = 0
        name = 'not classified'

        print('{}\t{}\t{}\t{}\t{}\t{}'.format(percent, count_below, count_at,
                                              classify_code, taxid, name))


if __name__ == '__main__':
    main()
