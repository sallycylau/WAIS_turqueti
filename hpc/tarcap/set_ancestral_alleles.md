Lau_etal_seaway
================
Sally Lau
30/03/2022

## Introduction

## Set ancestral alleles for P. turqueti

See <https://github.com/jessstapley/Set_ancestral_allele_vcf> for
complete tutorial <br> Keep outgroup samples (Pareledone cornuta: CT931,
Pareledone aequipapillae: 44064_1)

``` bash
cd ./turqueti/tarcap/SNPfiltering_outgroup/1_keep_outgroup

/sw/containers/vcftools-0.1.16.sif vcftools --vcf ./turqueti/tarcap/SNPfiltering1/1_rmindel/rmindels.recode.vcf --keep outgroup.tsv --out out --recode --recode-INFO-all
```

Keep biallelic SNPs only

``` bash
cd ./turqueti/tarcap/SNPfiltering_outgroup/2_biallelic_only

/sw/containers/vcftools-0.1.16.sif vcftools --vcf ./turqueti/tarcap/SNPfiltering_outgroup/1_keep_outgroup/out.recode.vcf --min-alleles 2 --max-alleles 2 --out out_biallelic --recode --recode-INFO-all
```

Keep SNPs present in both species

``` bash
cd ./turqueti/tarcap/SNPfiltering_outgroup/3_maxmissing1

/sw/containers/vcftools-0.1.16.sif vcftools --vcf ./turqueti/tarcap/SNPfiltering_outgroup/2_biallelic_only/out_biallelic.recode.vcf --max-missing 1 --out out_maxmissing1 --recode --recode-INFO-all
```

Create a table of SNP positions and alleles using bcftools

``` bash
/sw/containers/bcftools-1.13.sif bcftools annotate --set-id '%CHROM\_%POS' out_maxmissing1.recode.vcf > out_maxmissing1_setID.vcf

/sw/containers/bcftools-1.13.sif bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' out_maxmissing1_setID.vcf  > out_maxmissing1.tab
```

Create a table with ancestral allele information. If the site is not
identical across both outgroups, then record it as missing, i.e. ‘.’

``` bash
awk '{OFS="\t";if($5=="0/0" && $6=="0/0"){print $1,$2,$3,$4,$3} \
    if($5=="0/1" && $6=="0/1"){print $1,$2,$3,$4,$4} \
    if($5=="1/1" && $6=="1/1"){print $1,$2,$3,$4,$4} \
    if($5=="0/0" && $6=="0/1"){print $1,$2,$3,$4,"./."} \
    if($5=="0/0" && $6=="1/1"){print $1,$2,$3,$4,"./."} \
    if($5=="0/1" && $6=="0/0"){print $1,$2,$3,$4,"./."} \
    if($5=="0/1" && $6=="1/1"){print $1,$2,$3,$4,"./."} \
    if($5=="1/1" && $6=="0/0"){print $1,$2,$3,$4,"./."} \
    if($5=="1/1" && $6=="0/1"){print $1,$2,$3,$4,"./."}}' out_maxmissing1.tab > out_maxmissing1_aa.tab
```

Compress and index \*\_aa.tab file

``` bash
bgzip out_maxmissing1_aa.tab

tabix -s1 -b2 -e2 out_maxmissing1_aa.tab.gz
```
