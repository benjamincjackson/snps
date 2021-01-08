## snps

Find nucleotide changes between each sequence in an alignment, relative to a reference sequence.

### usage

To build the binary, you need [go](https://golang.org/).

```
git clone https://github.com/benjamincjackson/snps.git
cd snps/
go build snps.go

./snps reference.fasta alignment.fasta > snps.csv
