#! /usr/bin/env python
"""
Convert `sourmash gather`'s CSV output into a lineage file.
"""
import argparse
import csv
import ncbi_taxdump_utils


want_taxonomy = ['superkingdom', 'phylum', 'order', 'class', 'family', 'genus', 'species']


def main():
    p = argparse.ArgumentParser()
    p.add_argument('nodes_dmp')
    p.add_argument('names_dmp')
    p.add_argument('accession_csv')
    p.add_argument('gather_csv')
    p.add_argument('-o', '--output', type=argparse.FileType('wt'))
    a = p.parse_args()

    taxfoo = ncbi_taxdump_utils.NCBI_TaxonomyFoo()
    taxfoo.load_nodes_dmp(a.nodes_dmp)
    taxfoo.load_names_dmp(a.names_dmp)
    taxfoo.load_accessions_csv(a.accession_csv)

    csv_w = None
    if a.output:
        csv_w = csv.writer(a.output)
        csv_w.writerow(['name', 'lineage'])

    with open(a.gather_csv, 'rt') as gather_fp:
        r = csv.DictReader(gather_fp)
        for row in r:
            name = row['name']
            acc = name.split(' ')[0]
            taxid = taxfoo.get_taxid(acc)
            if taxid:
                lineage = taxfoo.get_lineage(taxid, want_taxonomy)
                lineage = ";".join(lineage)
                print('For {}, found lineage {}'.format(acc, lineage))

                if a.output:
                    csv_w.writerow([name,lineage])


if __name__ == '__main__':
    main()
