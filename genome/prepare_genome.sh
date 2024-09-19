#!/bin/bash

# GRCh38
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/reference/1KG_ONT_VIENNA_hg38.fa.gz
zcat 1KG_ONT_VIENNA_hg38.fa.gz | bgzip > hg38.fa.gz
samtools faidx hg38.fa.gz
rm 1KG_ONT_VIENNA_hg38.fa.gz

# CHM13
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/reference/1KG_ONT_VIENNA_t2t.fa.gz
zcat 1KG_ONT_VIENNA_t2t.fa.gz | bgzip > t2t.fa.gz
samtools faidx t2t.fa.gz
rm 1KG_ONT_VIENNA_t2t.fa.gz
