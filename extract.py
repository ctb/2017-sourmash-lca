#! /usr/bin/env python
"""
Emulate the kraken database-building step, but using banded MinHash signatures
instead.

Briefly,

* load in a list of accessions-to-taxonomic-lineage information (e.g. such
  as the .csv.gz file from https://github.com/dib-lab/2017-ncbi-taxdump);
* for every hash in each signature, connect that hash to the lineage for that
  signature (signature genbank accession ID -> lineage);
* after loading in all signatures & hashes, find the last-common-ancestor
  lineage for each hash and associate in a dictionary;
* save dictionary.

Usage::

   kraken/extract.py genbank/*.csv.gz ecoli_many_sigs/ecoli-*.sig --savename foobar

The E. coli signatures used in the command above can be downloaded like so:

   curl -O -L https://github.com/dib-lab/sourmash/raw/master/data/eschericia-sigs.tar.gz

----

This could be pretty easily broken down into multiple steps:

* step 1: link hashes to lineages and save those links to a file;
* step 1b: update this file as many times as you want;
* step 2: find last-common-ancestor taxid for each hash.

This would make adding in new organisms to the LCA database much
faster, though no less memory intensive.
"""

import khmer
import sourmash_lib, sourmash_lib.signature
import argparse
import gzip, csv
from pickle import dump


def main():
    p = argparse.ArgumentParser()
    p.add_argument('genbank_csv')
    p.add_argument('sigs', nargs='+')
    p.add_argument('-k', '--ksize', default=31)
    p.add_argument('-s', '--savename', default=None, type=str)
    args = p.parse_args()

    if not args.savename:
        print('no savename, quitting')
        sys.exit(0)


    xopen = open
    if args.genbank_csv.endswith('.gz'):
        xopen = gzip.open

    # load the genbank CSV (accession -> lineage)
    print('loading genbank accession -> lineage info')
    with xopen(args.genbank_csv, 'rt') as fp:
        accessions = {}
        for row in csv.DictReader(fp, fieldnames=['acc', 'taxid', 'lineage']):
            acc = row['acc']
            accessions[acc] = row

    # code to find the last common ancestor from the lineage string
    def find_lca(acc_list):
        info_list = []
        for acc in acc_list:
            info = accessions.get(acc)
            if info:
                info_list.append(info['lineage'])

        while 1:
            lins = set(info_list)
            if len(lins) <= 1:
                if len(lins):
                    return lins.pop()
                else:
                    return 'XXX'

            new_list = []
            for x in info_list:
                new_list.append(x.rsplit(';', 1)[0])
            info_list = new_list


    # initialize a GraphLabels/labelhash data structure to track
    # hashes -> taxonomic IDs.  We could use a dictionary/list here
    # but I think it would not scale as well? Not sure.
    lh = khmer.GraphLabels(args.ksize, 1, 1)
    id_to_acc = {}

    # for every minhash in every signature, link it to the taxonomic lineage
    # it's from.
    print('loading signatures & traversing hashes')
    for n, filename in enumerate(args.sigs):
        sig = sourmash_lib.signature.load_one_signature(filename,
                                                      select_ksize=args.ksize)
        acc = sig.name().split(' ')[0]   # first part of sequence name
        acc = acc.split('.')[0]          # get acc w/o version

        # @CTB hack hack split off NZ from accession
        if acc.startswith('NZ_'):
            acc = acc[3:]
        
        mins = sig.minhash.get_mins()

        id_to_acc[n] = acc

        for m in mins:
            lh.add_tag(m)
            lh.link_tag_and_label(m, n)
    print('...done')

    print('traversing tags and finding last-common-ancestor')
    cur_id = 1
    lineage_to_id = {}
    id_to_lineage = {}
    tag_to_lid = {}
    for n, tag in enumerate(lh.get_tagset()):
        tag = lh.hash(tag)

        # get all of the accessions associated with this ID, & find associated
        # lineages.  We could probably move this into the loop above
        # and go directly from hashes to lineages.
        labels = lh.get_tag_labels(tag)
        acc_list = [ id_to_acc[x] for x in labels ]
        lca = find_lca(acc_list)

        # get the lineage_id of the lca if we know it --
        lineage_id = lineage_to_id.get(lca)

        # new lineage? add, and do bookkeeping.
        if lineage_id == None:
            lineage_id = cur_id
            lineage_to_id[lca] = lineage_id
            id_to_lineage[lineage_id] = lca
            cur_id += 1

        # save!!
        tag_to_lid[tag] = lineage_id
    print('done')

    print('saving to', args.savename)
    dump((tag_to_lid, id_to_lineage, lineage_to_id),
         open(args.savename, 'wb'))


if __name__ == '__main__':
    main()
