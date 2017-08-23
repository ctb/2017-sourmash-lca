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

   kraken/classify.py foobar ecoli_many_sigs/ecoli-1.sig
"""

import khmer
import sourmash_lib, sourmash_lib.signature
import argparse
import gzip, csv
from pickle import load
import collections


def main():
    p = argparse.ArgumentParser()
    p.add_argument('saved_db')
    p.add_argument('sigfile')
    p.add_argument('-k', '--ksize', default=31)
    args = p.parse_args()

    (tag_to_lid, id_to_lineage, lineage_to_id) = load(open(args.saved_db, 'rb'))
    # load signature
    sig = sourmash_lib.signature.load_one_signature(args.sigfile,
                                                    select_ksize=args.ksize)

    # for every hash, print out LCA of labels
    cnt = collections.Counter()

    found = 0
    for n, tag in enumerate(sig.minhash.get_mins()):
        lid = tag_to_lid.get(tag)
        if lid is None:
            continue
        
        lineage = id_to_lineage.get(lid)
        if lineage is None:
            continue

        found += 1
        cnt[lineage] += 1

    print('found database classifications for', found, 'of', n + 1, 'hashes')
    for item, count in cnt.most_common():
        print(count, item)


if __name__ == '__main__':
    main()
