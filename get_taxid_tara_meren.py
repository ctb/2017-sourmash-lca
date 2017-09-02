#! /usr/bin/env python
import lca_json
import argparse
import csv
import traceback


want_taxonomy = ['superkingdom', 'phylum', 'order', 'class', 'family', 'genus', 'species']


def main():
    p = argparse.ArgumentParser()
    p.add_argument('lca_db')
    p.add_argument('csv')
    p.add_argument('outcsv')
    args = p.parse_args()

    lca_db = lca_json.LCA_Database(args.lca_db)
    taxfoo = lca_db.get_taxonomy()

    r = csv.reader(open(args.csv, 'rt', encoding='utf8'))
    next(r)

    results = []

    n = 0
    for row in r:
        n += 1
        genome_name = row[0]
        taxid = 1

        lineage = taxfoo.get_lineage(taxid, want_taxonomy)
        results.append((genome_name, str(taxid), ";".join(lineage)))
        print(genome_name, taxid, ";".join(lineage))

    with open(args.outcsv, 'wt', encoding='utf8') as outfp:
        w = csv.writer(outfp)
        for name, taxid, lineage in results:
            w.writerow([name, taxid, lineage])

if __name__ == '__main__':
    main()
