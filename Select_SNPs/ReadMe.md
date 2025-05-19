# This is an  overview how to extract the SNPs 

### Download 1000 genomes data
./download.sh

### Download sampling location list from
Here, European samples are from CEU, IBS, GBR, FIN, TSI.

Download the samples from https://www.internationalgenome.org/data-portal/sample and call it igsr_samples.tsv.

### Extract European samples
./extract_european_samples.sh

Now we follow
https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/admixture_files/README.admixture_20141217

### Filter bi-allelic SNPs which are 10k apart
./filtered.sh


### Connection to the paper

In our paper, we basically used the same code. However, we adapted it to only consider individuals from GBR/IBS and from IBS, CEU, TSI, respectively.
