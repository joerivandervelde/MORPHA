# MORPHA
MOlecular effect Range of germline Pathogenic variation in global protein HAplotypes

To reproduce the small genes analysis, run these scripts:
- `0-all-small-genes.sh`
- `1-extract-small-gene-variants-from-clinvar.sh`
- `2-signal-peptide-lengths.sh`
- `3-reconstruct-1kg-patho-haplotypes.sh`
- `4-fold.sh`
- `5-repair.sh`
- `6-relax.sh`
- `7-lowest.sh`
- `8-analysis.R`

These make use of the following reusable components:
- `./genename2uniqproteinhaplo.sh <gene> <RefSeq protein ID> <front signal peptide cutoff, optional, default = 0>`
    - Using GRCh38 and 1000G data, reconstruct and group all 1KG haplotypes.
    - If front signal peptide cutoff is provided, the first N residues will be removed from each output protein sequence.
- `./genename2pathohaplo.sh <gene> <RefSeq protein ID> <localized vcf.gz> <front signal peptide cutoff, optional, default = 0>`
    - Spike pathogenic variants into 1KG haplotypes and group them.
    - A custom variant vcf.gz can be supplied here:
        - Made with `bgzip -c abc.vcf > abc.vcf.gz` followed by `tabix -p vcf abc.vcf.gz`.
        - BEWARE: the VCF `POS` field must be calculated as: `(hg38 location) - (GFF transcript start location) + 1`.
        - Also, `ID` field may NOT include `:`.
- `./genename2foldedhaplo.sh <gene>`
    - Fold all unique haplotype proteins, 1KG and spiked.
    - Expects [esmfold-server](https://github.com/joerivandervelde/esmfold-server) to run on your network, default address http://192.168.2.9:8000, adjust accordingly.

General script dependencies:
- samtools 1.23.1 (also tested with 1.13, seems OK)
- bcftools 1.23.1
- CrossMap 0.7.0
- gffread 0.12.9
- tabix 1.17 (also tested with 1.13, seems OK)

## Notes

### Chromosome X needs manual download and rename:
```
cd data
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz.tbi
mv 1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz 1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
mv 1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz.tbi 1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.vcf.gz.tbi
```

## BCFtools version must be >1.21
- 1.13 and earlier: Fails with error `consensus: unrecognized option '--regions-overlap'`
- 1.14 – 1.20: Ignores all alt alleles from liftover pathogenic VCF because it required exact sample matching to variant VCF to be applied
- 1.21 or later: **Needed to run succesfully**. Smart fallback was restored, defaults to applying alt alleles without exact sample match.


# Ubuntu WSL install and run

## Install the following
```
sudo apt install samtools
sudo apt install gffread
sudo apt install python3-pip
sudo apt install tabix
sudo apt install dos2unix # use to convert line endings of script you want to run, e.g. dos2unix *.sh
pip install CrossMap
```

BCFtools: download and compile 1.23.1 (or later)
```
wget https://github.com/samtools/bcftools/releases/download/1.23.1/bcftools-1.23.1.tar.bz2
tar -xjf bcftools-1.23.1.tar.bz2
cd bcftools-1.23.1
./configure
make
sudo make install
hash -r
bcftools --version
```
