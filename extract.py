import khmer
import sourmash_lib, sourmash_lib.signature
import argparse
import gzip, csv
from pickle import dump


def main():
    p = argparse.ArgumentParser()
    p.add_argument('sigs', nargs='+')
    p.add_argument('-k', '--ksize', default=31)
    p.add_argument('-s', '--savename', default=None, type=str)
    args = p.parse_args()

    if not args.savename:
        print('no savename, quitting')
        sys.exit(0)

    lh = khmer.GraphLabels(args.ksize, 1, 1)

    id_to_acc = {}

    # for every minhash in every signature, link it to the signature it's from.
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

    lh.save_labels_and_tags(args.savename)
    dump(id_to_acc, open(args.savename + '.id_to_acc', 'wb'))


if __name__ == '__main__':
    main()
