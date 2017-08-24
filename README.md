# banded kraken - example usage

Prepare the database:

```
mkdir ecoli_many_sigs
cd ecoli_many_sigs
curl -O -L https://github.com/dib-lab/sourmash/raw/master/data/eschericia-sigs.tar.gz
tar xzf eschericia-sigs.tar.gz
cd ../

curl -O -L ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar xzf taxdump.tar.gz nodes.dmp names.dmp

curl -O -L https://github.com/dib-lab/2017-ncbi-taxdump/raw/master/genbank-genomes-accession%2Blineage-20170529.csv.gz

# do this once
extract.py genbank*.csv.gz nodes.dmp ecoli_many_sigs/ecoli-*.sig --savename foobar
```

Now, run classification on any signatures you have lying around:

```
# do this many times
classify.py nodes.dmp names.dmp foobar ecoli_many_sigs/ecoli-1.sig
```

