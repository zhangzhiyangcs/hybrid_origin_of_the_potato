## Acquire the sequences from gff

```shell
# reformat genome fasta to 60 bp size at per line
~/software/seqkit seq genome.fasta -w 60 >re_genome.fasta

# get gene sequences
for i in `less sample.txt`; do echo "agat_sp_extract_sequences.pl -gff ${i}.reformat.gff3  --fasta ${i}.fasta -t gene -o  ${i}.gene.fa"; done 

# get CDS
for i in `less sample.txt`; do echo "agat_sp_extract_sequences.pl -gff ${i}.reformat.gff3 --fasta ${i}.fasta -t cds -o ${i}.cds.fa"; done

# get intron sequences
for i in `less sample.txt`; do echo "agat_sp_extract_sequences.pl -gff ${i}.intron.gff3 --fasta ${i}.fasta -t intron --merge  -o ${i}.intron.fa"; done

# get protein sequences
seqkit translate cds.fa >protein.fa
```



## Construct the pseudo Etuberosum genome

```shell
# Generate VCF files based on the Solanum etuberosum genome FASTA and following the previously established SNP calling pipeline.
BWA-GATK-SAMtools pipeline -> solanum_etuberosum.vcf

# Filter VCF file according to the depth, quality and missing rate, and bialleic alleles
/public/software/env01/bin/vcftools --minDP 4 --maxDP 120 --vcf solanum_etuberosum.vcf  --recode --recode-INFO-all --out ./filter_snp_near_indel --minQ 20  --max-missing 0.8 --min-alleles 2 --max-alleles 2

# Filter the SNPs around 5 bp of the Indel
/public/software/env01/bin/bcftools filter --threads 10 -g 5 filter_snp_near_indel.recode.vcf   --output-type z -o filter_snp_near_5bp_indel.gz

# Filter missing rate according to 3 lineages
awk '{if($2=="E") print $1}' popman >1.keep
awk '{if($2=="P") print $1}' popman >2.keep
awk '{if($2=="T") print $1}' popman >3.keep
awk '{if($2=="egg") print $1}' popman >4.keep

# Caculate missing rate via VCFtools
/public/software/env01/bin/vcftools --vcf filter_bad.recode.vcf --keep 1.keep --missing-site --out 1 &
/public/software/env01/bin/vcftools --gzvcf filter_snp_near_5bp_indel.gz --keep 2.keep --missing-site --out 2 &
/public/software/env01/bin/vcftools --gzvcf filter_snp_near_5bp_indel.gz --keep 3.keep --missing-site --out 3 &
/public/software/env01/bin/vcftools --gzvcf filter_snp_near_5bp_indel.gz --keep 4.keep --missing-site --out 4 &

# Filter the sites which missing rate >20%
cat 1.lmiss 2.lmiss 3.lmiss 4.lmiss | mawk '!/CHR/' | mawk '$6 > 0.2' | cut -f1,2 > badloci
/public/software/env01/bin/vcftools --gzvcf filter_snp_near_5bp_indel.gz --exclude-positions badloci --maf 0.05 --recode --recode-INFO-all --out filter_snp_near_5bp_missing_lineage_indel

# Filter the indel
/public/software/env01/bin/vcftools --remove-indels --recode --recode-INFO-all --vcf filter_snp_near_5bp_missing_lineage_indel.recode.vcf --stdout >filter_snp_near_5bp_missing_lineage_indel.recode.snp.vcf

# Calculate the missing rate of each individual 
for  i  in  {01..12}; do vcftools --vcf ${i}.filter_depth_quality.recode.vcf  --missing-indv;done

# Replace all occurrences of the string "\s.:" in the input text with "\t./.:"
for  i  in  {01..12}; do  perl -pe "s/\s\.:/\t.\/.:/g" chr${i}.filter.recode.vcf >chr${i}.vcf;done

# Construct the pseudo Etuberosum genome (https://github.com/vcflib/vcflib/blob/master/doc/vcf2fasta.md)
vcf2fasta -f solanum_etuberosum.fasta -p pseudo -P 2 -n N solanum_etuberosum.vcf


```



## Acquire the unique synteny-constrained orthogroups by GENESPACE 

```shell
# The concrete instruction for GENESPACE can be download from (https://github.com/jtlovell/GENESPACE)
# activate enviroment
1) source activate orthofinder
# split diamond parts from genespace
2) orthofinder -f 01.peptide -op| grep 'diamond blastp' > diamond.commond
# find orthogroup
3) orthofinder -b Results_Sep24/ -t 140 -a 140
# parallel commands

# If you run the orthofinder seperately
# need 3 files SequenceIDs.txt,Orthogroups.tsv,SpeciesIDs.txt
gpar$params$nCores<-100
gpar$paths$blastDir <-"/Results_Oct03/WorkingDirectory"
gpar$paths$orthogroupsDir <-"/Results_Oct03_1/Orthogroups"
gpar$paths$orthologuesDir<-"/Results_Oct03_1/Orthologues"


# Prepartion raw data for analysis
# construct rawGenomes format files for Genespace
for i in `less ID`; do 
mkdir -p rawGenomes/${i}/${i}/annotation
cd rawGenomes/${i}/${i}/annotation
mv ../../../../${i}_pep.fa ./
mv ../../../../${i}_gene.gff ./
gzip ./*
cd ../../../../
done


library(GENESPACE)
runwd <- file.path("/home/project/05.HHS_petota/02.Genespace")
ID<-read.csv("ID",header=F)
ID<-as.vector(ID$V1)
# Initialize the GENESPACE run
gpar <- init_genespace(
  genomeIDs = ID, 
  speciesIDs = ID, 
  versionIDs = ID, 
  ploidy = rep(1,83), ##ploid and sample number
  wd = runwd, 
  gffString = "gff", 
  pepString = "pep",
  nCores = 100,
  path2orthofinder = "orthofinder", 
  path2mcscanx = "/home/softwares/MCScanX-master",
  rawGenomeDir = file.path(runwd, "rawGenomes"))
  
  parse_annotations(
  gsParam = gpar, 
  gffEntryType = "mRNA",   ##the third column
  gffIdColumn = "ID",    ###the ninth coplum
  gffStripText = "Parent=", 
  headerEntryIndex = 1,
  headerSep = " ", 
  headerStripText = "Parent=")
gpar <- synteny(gsParam = gpar)

# stats gffWithOgs.txt.gz (select TRUE from this file)
# 83 species +- 10 species
# work dir:/home/project/05.HHS_petota/02.Genespace/stats
sh readme
less -S gffWithOgs.txt  |awk '{print $1"\t"$2"\t"$NF}' >TRUE.txt
less -S TRUE.out  |awk '{print $1"\t"$3}' |sort |uniq -d >dup.txt
python ./01.script/rm_dup.py  >rm_dup.txt
less -S rm_dup.txt |awk '{print $3}' |sort |uniq -c  >rm_dup_count.txt
less rm_dup_count.txt  |awk '{if ($1>73 && $1<93) print $_}' >rm_dup_ID.txt
python ./01.script/TRUE.py rm_dup_ID.txt  TRUE.txt >TRUE.out
python ./01.script/rm_dup_73_93_ID.py rm_dup_ID.txt  TRUE.txt>rm_dup_ID.out

rm rm_dup.txt rm_dup_count.txt rm_dup_73_93_ID.txt rm_dup_73_93_ID.out
less -S rm_dup_73_93.out | grep -v "Ipomoeaaquatica" >rm_Ipomoeaaquatica.out

# genespace id convert
less -S rm_dup_ID.out2|grep "DM" |awk '{print $2"\t"$3}' >DM_og.txt
python change_og_genespace.py  previous_og_DM.txt DM_og.txt gffWithOgs.txt TRUE.txt
 awk  '!a[$1,$3]++' select_TRUE.txt >rm_dup_select_TRUE.txt
```

