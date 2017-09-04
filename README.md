# banded kraken - example usage

Prepare the database:

```
mkdir ecoli_many_sigs
cd ecoli_many_sigs
curl -O -L https://github.com/dib-lab/sourmash/raw/master/data/eschericia-sigs.tar.gz
tar xzf eschericia-sigs.tar.gz
cd ../

mkdir -p ecoli/genbank
cd ecoli/genbank
curl -O -L ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar xzf taxdump.tar.gz nodes.dmp names.dmp
cd ../..

curl -O -L https://github.com/dib-lab/2017-ncbi-taxdump/raw/master/genbank-genomes-accession%2Blineage-20170529.csv.gz

# do this once
extract.py ecoli/ecoli.lca genbank*.csv.gz ecoli/genbank/nodes.dmp ecoli_many_sigs/ecoli-*.sig --lca-json=ecoli/ecoli.lca.json
```

Now, run classification on any signatures you have lying around:

```
# do this many times
classify.py ecoli/ecoli.lca.json ecoli_many_sigs/ecoli-1.sig
```


## Using genbank LCA database

## (constructed for k=21, 31, and 51)

Download 'sourmash-lca-mark1/genbank-lca-2017.08.26.tar.gz' from [osf.io/zfmbd/][https://osf.io/zfmbd/] and unpack it --

```
curl -L https://osf.io/zfmbd/download?version=1 -o genbank-lca-2017.08.26.tar.gz
tar xzf genbank-lca-2017.08.26.tar.gz
```

This will give you a `genbank.lca.json` file along with a `genbank/` subdirectory.

Now execute:

```
classify.py -k 31 genbank.lca.json sigfile.sig
```

## Taxonomy file sources

tara_meren_taxids.csv from
https://www.dropbox.com/s/xplbkhnlim5wj1j/Table_S3.xlsx?dl=0
