#! /usr/bin/env python
import argparse
import csv
import traceback
from ncbi_taxdump_utils import NCBI_TaxonomyFoo


want_taxonomy = ['superkingdom', 'phylum', 'order', 'class', 'family', 'genus', 'species']


def main():
    p = argparse.ArgumentParser()
    p.add_argument('names_dmp')
    p.add_argument('nodes_dmp')
    p.add_argument('csv')
    p.add_argument('outcsv')
    args = p.parse_args()

    taxfoo = NCBI_TaxonomyFoo()
    taxfoo.load_nodes_dmp(args.nodes_dmp, False)
    taxfoo.load_names_dmp(args.names_dmp, False)

    r = csv.reader(open(args.csv, 'rt', encoding='utf8'))
    next(r)

    results = []

    n = 0
    for row in r:
        n += 1
        genome_name = row[0]
        taxid = 131567

        lineage = taxfoo.get_lineage(taxid, want_taxonomy)
        results.append((genome_name, str(taxid), ";".join(lineage)))
        print(genome_name, taxid, ";".join(lineage))

    with open(args.outcsv, 'wt', encoding='utf8') as outfp:
        w = csv.writer(outfp)
        for name, taxid, lineage in results:
            w.writerow([name, taxid, lineage])

if __name__ == '__main__':
    main()
