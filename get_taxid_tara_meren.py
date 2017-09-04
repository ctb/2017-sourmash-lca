#! /usr/bin/env python
import lca_json
import argparse
import csv
import traceback


want_taxonomy = ['superkingdom', 'phylum', 'order', 'class', 'family', 'genus', 'species']

def get_taxids_for_name(taxfoo, **ranknames):
    assert len(ranknames) == 1
    srank, sname = list(ranknames.items())[0]

    tid_set = []
    for taxid, (name, _, _) in taxfoo.taxid_to_names.items():
        if name == sname:
            tid_set.append(taxid)

    lineages = [taxfoo.get_lineage_as_dict(tid, want_taxonomy) for tid in tid_set]

    if len(lineages) > 1:
        # if non-unique, collect one that matches the desired rank
        tid_set2 = []
        lin2 = []
        for (tid, lin) in zip(tid_set, lineages):
            if lin.get(srank) == sname:
                tid_set2.append(tid)
                lin2.append(lin)

        lineages = lin2
        tid_set = tid_set2
                
        if len(tid_set) > 2:
            assert 0, (tid_set, lineages)
            return -1

    if not lineages:
        print('not found:', srank, sname)
        return -1
    
    taxid = tid_set[0]
    return taxid


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
        last_rank = None
        last_name = None

        genome_name, row = row[0], row[1:]

        for rank, name in zip(want_taxonomy, row):
            if name == 'na':
                break
            last_rank, last_name = rank, name
            
        taxid = get_taxids_for_name(taxfoo, **{ last_rank: last_name })
        if taxid == -1:
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
