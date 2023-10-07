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

