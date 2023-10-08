## CAll SNP

```shell
#########seperate chr -L
export PATH=/public/agis/softwares/samtools-1.9/bin:$PATH
export PATH=/public/agis/softwares/gatk-4.1.9.0:$PATH


ref=/vol3/agis/work/1_reference/DM_v6.1_all_chr.fa
path=/vol3/agis/work/6_pangenome/01_ccs_data

for i in $(cat list);
do
    for j in 01 02 03 04 05 06 07 08 09 10 11 12; do
    echo '
for i in '$i';
do
export PATH=/public/agis/softwares/samtools-1.9/bin:$PATH
export PATH=/public/agis/softwares/gatk-4.1.9.0:$PATH
export PATH=/public/agis/softwares/jre1.8.0_251/bin/:$PATH


tembam=/vol3/agis/work/6_pangenome/02_snpbam
ref=/vol3/agis/work/1_reference/DM_v6.1_all_chr.fa
path=/vol3/agis/work/6_pangenome/01_ccs_data
#samtools addreplacerg -r 'ID:'$i'' -r 'PL:Illumina' -r 'SM:'$i'' -r 'LB:lib1' -o $tembam/$i.AR.bam 01_bam/$i.bam
#rm $tembam/$i.markdup.bam
#samtools index $tembam/$i.AR.bam
gatk --java-options "-Xmx20G" HaplotypeCaller -R $ref --ERC GVCF -I 01_bam/$i.AR.bam -O 02_gvcf/$i.chr'$j'.g.vcf -L chr'$j'
gatk --java-options "-Xmx20G" GenotypeGVCFs -R $ref -V 02_gvcf/$i.chr01.g.vcf -O 03_vcf/$i.chr01.vcf
done
' >gatk-$i.chr$j.sh
chmod +x gatk-$i.chr$j.sh
    done
done


## submit some sh to one queue##
#for i in $(cat list1 list2 list3 list4);
for i in $(cat list1);
do
    for j in 01 02 03 04 05 06 07 08 09 10 11 12; do
    #nohup ./gatk-$i.sh >> gatk-$i.log 2>&1 &
    quick_qsub {-q queue8 -l mem=10G,nodes=1:ppn=1} ./gatk-$i.chr$j.sh
    #qsub -cwd -pe smp 1 ./gatk-$i.chr$j.sh
    done
done

for i in $(cat list);
do
export PATH=/public/agis/softwares/samtools-1.9/bin:$PATH
export PATH=/public/agis/softwares/gatk-4.1.9.0:$PATH
export PATH=/public/agis/softwares/jre1.8.0_251/bin/:$PATH

gatk GatherVcfs -I 03_vcf/$i.chr01.vcf -I 03_vcf/$i.chr02.vcf -I 03_vcf/$i.chr03.vcf -I 03_vcf/$i.chr04.vcf -I 03_vcf/$i.chr05.vcf -I 03_vcf/$i.chr06.vcf -I 03_vcf/$i.chr07.vcf -I 03_vcf/$i.chr08.vcf -I 03_vcf/$i.chr09.vcf -I 03_vcf/$i.chr10.vcf -I 03_vcf/$i.chr11.vcf -I 03_vcf/$i.chr12.vcf -O 04_one_vcf/$i.gatk.vcf

done
```



## Hyde

```shell
species_num=$1
position_num=$2
python ~/software/hybrid/HyDe/scripts/individual_hyde_mp.py -i hyde.fa -m sample.ID.txt -o   Outgroup  -n $species_num -t 3 -s $position_num -q -tr triples.txt.for-ind --prefix Individual
```



## Loter

 SNPs were phased by BEAGLE software
```shell
vcf_gz=$1
vcf_gz2=$2
vcf_gz3=$3
zcat $1 | perl -pe "s/\s\.:/\t.\/.:/g" | bgzip -c >$vcf_gz2
wait
mv $vcf_gz2 $vcf_gz
java -Xmx888888m -jar  ~/software/hybrid/beagle.08Feb22.fa4.jar gt=$vcf_gz   out=$vcf_gz3
```

Generate input files for Loter software
```shell 
ID=$1
gzvcf=$2
out=$3
/public/software/env01/bin/vcftools --keep $ID  --gzvcf $gzvcf  --recode --recode-INFO-all --out $out


loter_cli  -r E01_ETB_snps.vcf.recode.vcf   E01_tomato_snps.vcf.recode.vcf   -a E01_potato_snps.vcf.recode.vcf -f vcf -o LOTER.out.txt -n 15 -v

perl 001.Loter_AncestryTractLength.pl  -V potato_chr_01.vcf  -R LOTER_chr01.out.txt  -O LOTER1.total.step1

perl 002.SummaryGenome_or_Chrs.TractLength.pl LOTER1.total.step1
```



## calculate fd

```shell
zcat parse_431.vcf.gz |grep -E "#|chr01" >parse_chr01.vcf
zcat parse_431.vcf.gz |grep -E "#|chr12" >parse_chr12.vcf

python  ~/software/scripts/ABBABABAwindows.py -g parse_chr01.vcf  -f phased  -o chr01_output1  -P1 Tomato  -P2 Potato   -P3 ETB  -O Outgroup  --popsFile 431_sample_ID.txt  --minData 0.5 -w 100000 -m 50  -T 30
python  ~/software/scripts/ABBABABAwindows.py -g parse_chr12.vcf  -f phased  -o chr12_output1  -P1 Tomato  -P2 Potato   -P3 ETB  -O Outgroup  --popsFile 431_sample_ID.txt  --minData 0.5 -w 100000 -m 50  -T 30
```



## calculate recombination rate

```shell
# step1
#!/bin/bash
~/miniconda3/bin/R --no-save <<EOF
library(FastEPRR)
#dir.create("chr${1}")
FastEPRR_VCF_step1(vcfFilePath="./potato_chr_${1}.vcf",srcOutputFilePath="./chr${1}",winLength="100",winDXThreshold=30)
EOF

# for submit
for i in {01..12}; do echo "qsub -q queue1,big CHB_step1.sh ${i}"; done

# step2
#!/bin/bash
~/miniconda3/bin/R --no-save <<EOF
library(FastEPRR)
FastEPRR_VCF_step2(srcFolderPath="chr12_rr", jobNumber= 100,currJob=${1},DXOutputFolderPath="chr12_output")
EOF

# for submit
#!/bin/bash
for((i=1;i<101;i=i+1));
do
    qsub -q queue1,big CHB_step2.sh ${i}
done

# step3
#!/bin/bash
~/miniconda3/bin/R --no-save <<EOF
library(FastEPRR)
FastEPRR_VCF_step3(srcFolderPath="chr2_rr",DXFolderPath="chr2_output",finalOutputFolderPath="step3_outputData")
EOF

# for submit

```





```shell
# extract the vcf from the specific region
/public/software/env01/bin/bcftools filter  population.vcf.gz   --regions  chr:pos1-pos2 >~/specific_region.vcf
bcftools query -f '%CHROM %POS %REF %ALT [ %GT]\n'  specific_region.vcf
```

