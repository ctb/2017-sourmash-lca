import khmer
import sourmash_lib, sourmash_lib.signature
import argparse
import gzip, csv
from pickle import load
import collections


def main():
    p = argparse.ArgumentParser()
    p.add_argument('genbank_csv')
    p.add_argument('saved_db')
    p.add_argument('sigfile')
    p.add_argument('-k', '--ksize', default=31)
    args = p.parse_args()

    lh = khmer.GraphLabels(args.ksize, 1, 1)
    print('loading labels and tags: {}'.format(args.saved_db))
    lh.load_labels_and_tags(args.saved_db)

    dict_file = args.saved_db + '.id_to_acc'
    id_to_acc = load(open(dict_file, 'rb'))

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


    # load signature
    sig = sourmash_lib.signature.load_one_signature(args.sigfile,
                                                    select_ksize=args.ksize)

    # for every hash, print out LCA of labels
    cnt = collections.Counter()

    print('n tags', len(sig.minhash.get_mins()))
    for n, tag in enumerate(sig.minhash.get_mins()):
        if n % 100 == 0:
            print('...', n)
        id_list = lh.get_tag_labels(tag)
        if not id_list:
            continue

        acc_list = [id_to_acc[i] for i in id_list]
        lca = find_lca(acc_list)
        cnt[lca] += 1

    for item, count in cnt.most_common():
        print(item, count)


if __name__ == '__main__':
    main()
