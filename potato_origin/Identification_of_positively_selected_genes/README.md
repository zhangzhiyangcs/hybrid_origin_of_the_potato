





```shell
protein_mafft=${1}.mafft
cds_mafft=${1}.pal
sampleID=${1}

pal2nal.pl ./04.prot_mafft/$protein_mafft ./02.cds_mafft/$cds_mafft  -output codon >${sampleID}.codon
pal2nal.pl ./04.prot_mafft/$protein_mafft ./02.cds_mafft/$cds_mafft  -output fasta >${sampleID}.fasta
python v2_check_syn_non-syn.py  ${sampleID}.codon ${sampleID}.syn DM

snp-sites -c ${sampleID}.fasta -v | sed 's/^1/'${sampleID}'/'  >${sampleID}.snp-sites1.vcf
python filter_missing.py ${sampleID}.snp-sites1.vcf ${sampleID}.snp-sites.vcf
python add_syn_in_snp-sites_vcf.py ${sampleID}.syn ${sampleID}.snp-sites.vcf ./snp_sites/${sampleID}.snp-sites.vcf


vcftools  --vcf ./snp_sites/${sampleID}.snp-sites.vcf --keep  P_E.txt --maf 0.05 --recode --recode-INFO-all  --max-missing 0.8 --out ./P_E_poly_sites/${sampl
vcftools  --vcf ./snp_sites/${sampleID}.snp-sites.vcf --keep  P_T.txt --maf 0.05 --recode --recode-INFO-all  --max-missing 0.8 --out ./P_T_poly_sites/${sampl

vcftools --vcf  ./snp_sites/${sampleID}.snp-sites.vcf  --weir-fst-pop T.txt  --weir-fst-pop P_E.txt   --max-missing 0.8 --out ${sampleID}.P_E.fst
vcftools --vcf  ./snp_sites/${sampleID}.snp-sites.vcf  --weir-fst-pop E.txt  --weir-fst-pop P_T.txt   --max-missing 0.8 --out ${sampleID}.P_T.fst

python  add_PT_fst.py  ${sampleID}.P_T.fst.weir.fst  ./snp_sites/${sampleID}.snp-sites.vcf  ${sampleID}.PT.vcf
python  add_PE_fst.py  ${sampleID}.P_E.fst.weir.fst  ./snp_sites/${sampleID}.snp-sites.vcf  ${sampleID}.PE.vcf
```



## stats the PSG

```shell
cat *.PT.vcf   |grep -v "#" |awk '{print $1"\t"$3"\t"$6}' |awk '{if($3>=0.95) print $_}' |grep -v "nan" >0.95_PT.stats 

cat *.PE.vcf   |grep -v "#" |awk '{print $1"\t"$3"\t"$7}' |awk '{if($3>=0.95) print $_}'|grep -v "nan" >0.95_PE.stats 

less -S 0.95_PE.stats |awk '{print $1"\t"$2}' |sort |uniq -c | awk '{if($3==1) print $2"\t"$1}'>PE_syn.stats 
less -S 0.95_PE.stats |awk '{print $1"\t"$2}' |sort |uniq -c | awk '{if($3==2) print $2"\t"$1}'>PE_nonsyn.stats 

less -S 0.95_PT.stats |awk '{print $1"\t"$2}' |sort |uniq -c | awk '{if($3==1) print $2"\t"$1}'>PT_syn.stats 
less -S 0.95_PT.stats |awk '{print $1"\t"$2}' |sort |uniq -c | awk '{if($3==2) print $2"\t"$1}'>PT_nonsyn.stats 

python append.py  PT_syn.stats  ID.txt  tmp1
python append.py  PT_nonsyn.stats  ID.txt tmp2
paste tmp1 tmp2 >PT_fix.stats

python append.py  PE_syn.stats  ID.txt  tmp1
python append.py  PE_nonsyn.stats  ID.txt tmp2
paste tmp1 tmp2 >PE_fix.stats

cat ./P_E_poly_sites/*.recode.vcf |grep -v "#"  | awk '{print $1"\t"$3}'|sort |uniq -c | awk '{if($3==1) print $2"\t"$1}' >./P_E_poly_sites/PE_syn_poly.stats
cat ./P_E_poly_sites/*.recode.vcf |grep -v "#"  | awk '{print $1"\t"$3}'|sort |uniq -c | awk '{if($3==2) print $2"\t"$1}' >./P_E_poly_sites/PE_nonsyn_poly.st
python append.py  ./P_E_poly_sites/PE_syn_poly.stats  ID.txt tmp1
python append.py  ./P_E_poly_sites/PE_nonsyn_poly.stats  ID.txt tmp2
paste tmp1 tmp2 >PE_poly.stats

cat ./P_T_poly_sites/*.recode.vcf |grep -v "#"  | awk '{print $1"\t"$3}'|sort |uniq -c | awk '{if($3==1) print $2"\t"$1}' >./P_T_poly_sites/PT_syn_poly.stats
cat ./P_T_poly_sites/*.recode.vcf |grep -v "#"  | awk '{print $1"\t"$3}'|sort |uniq -c | awk '{if($3==2) print $2"\t"$1}' >./P_T_poly_sites/PT_nonsyn_poly.st
python append.py  ./P_T_poly_sites/PT_syn_poly.stats  ID.txt tmp1
python append.py  ./P_T_poly_sites/PT_nonsyn_poly.stats  ID.txt tmp2
paste tmp1 tmp2 >PT_poly.stats
```

