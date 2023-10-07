## sample 500 combination 

```shell
python sample.py   select_ID.txt  ../data
ls |grep "Nicotianalongiflora"  >dir_ID.txt

# protein
for i in `less og_ID.txt`; do echo "~/miniconda3/bin/seqkit grep -f ID.txt ../../prot_data/${i}.prot.fa >${i}.prot.fa";done >01.seqkit.sh
for i in `less dir_ID.txt`;  do echo "cp 01.seqkit.sh seqkit_work.sh  ./${i}; cd  ./${i}; sbatch seqkit_work.sh 01.seqkit.sh; cd ../";done>01.sbatch_seqkit.sh


for i in `less og_ID.txt`; do echo "mafft ${i}.prot.fa >${i}_prot.mafft";done>02.prot_mafft.sh
for i in `less dir_ID.txt`;  do echo "cp 02.prot_mafft.sh  prot_mafft_work.sh  ./${i}; cd  ./${i}; sbatch prot_mafft_work.sh  02.prot_mafft.sh; cd ../";done>01.sbatch_prot_mafft.sh


for i  in `less og_ID.txt`; do echo "iqtree -s og_${i}_prot.mafft  -B 1000   -m MFP -nt 1 -o Nicotianalongiflora";done >02.iqtree.sh
for i in `less dir_ID.txt`;  do echo "cp 02.iqtree.sh iqtree_work.sh ./${i}; cd  ./${i}; sbatch iqtree_work.sh 02.iqtree.sh; cd ../";done>02.sbatch_iqtree.sh


## CDS
for i in `less og_ID.txt `; do echo " ~/miniconda3/bin/seqkit grep -f ID.txt ../../cds_data/${i}.fa >${i}.cds.fa";done >01.seqkit_cds.sh

for i in `less dir_ID.txt`;  do echo "cp 01.seqkit_cds.sh seqkit_work.sh  ./${i}; cd  ./${i}; sbatch seqkit_work.sh 01.seqkit_cds.sh; cd ../";done>01.sbatch_cds_seqkit.sh

for i in `less og_ID.txt `; do echo "mafft ${i}.cds.fa >${i}_cds.mafft";done >02.cds_mafft.sh
for i in `less dir_ID.txt`;  do echo "cp 02.cds_mafft.sh  cds_mafft_work.sh  ./${i}; cd  ./${i}; sbatch cds_mafft_work.sh  02.cds_mafft.sh; cd ../";done>01.sbatch_cds_mafft.sh


## stats 90 bootstrap
for i in `less dir_ID.txt`;  do echo " cp tmp1 tmp2  boots.R  ./${i};cd  ./${i};  Rscript boots.R >90_bootstrap.ID; awk '{print \$3}' 90_bootstrap.ID |sed 's/\"//' |xargs >tmp; paste tmp1 tmp tmp2 >cat.sh; sh cat.sh; cd ../"; done

##pal2nal
for i in `less og_ID.txt `; do echo " ~/miniconda3/bin/pal2nal.pl ${i}_prot.mafft ${i}_cds.mafft  -output fasta >${i}.pal";done >03.pal2nal.sh
for i in `less dir_ID.txt`;  do echo "cp 03.pal2nal.sh pal2nal_work.sh ./${i};cd  ./${i}; sbatch pal2nal_work.sh  03.pal2nal.sh ; cd ../";done>03.sbatch_pal2nal.sh

## snpsites
for i in `less og_ID.txt `; do echo "~/miniconda3/bin/snp-sites -c  -b -o ${i}.beast  ${i}.pal";done > 04.snpsites.sh 
for i in `less dir_ID.txt`;  do echo "cp 04.snpsites.sh snpsites_work.sh ./${i};cd  ./${i}; sbatch snpsites_work.sh  04.snpsites.sh ; cd ../";done>04.sbatch_snpsites.sh

## generate PT ET PE ID
for i in `less dir_ID.txt`; do echo "cd ./${i} ;    less -S 90_bootstrap.ID|grep \"PE\" |awk '{print \$3}' |awk -F '_prot' '{print \$1\".beast\"}' |xargs -n 1 >PE.ID; cd ../";done >>ID.sh
for i in `less dir_ID.txt`; do echo "cd ./${i} ;    less -S 90_bootstrap.ID|grep \"PT\" |awk '{print \$3}' |awk -F '_prot' '{print \$1\".beast\"}' |xargs -n 1 >PT.ID; cd ../";done >>ID.sh
for i in `less dir_ID.txt`; do echo "cd ./${i} ;    less -S 90_bootstrap.ID|grep \"ET\" |awk '{print \$3}' |awk -F '_prot' '{print \$1\".beast\"}' |xargs -n 1 >ET.ID; cd ../";done >>ID.sh

## prepare PE.fa PE.tree 
python PE_contentate.py  ./ PE.ID PE.fa
for i in `less dir_ID.txt`; do echo "cp PE_contentate.py  ET_contentate.py  PT_contentate.py ./${i};  cd ./${i}  ; python PE_contentate.py  ./ PE.ID PE.fa; python PT_contentate.py  ./ PT.ID PT.fa; python ET_contentate.py  ./ ET.ID ET.fa; cd ../";done >contentate.sh
## prepare appear 200 times sample
for i in `less dir_ID.txt`; do echo "cp PE_contentate.py  ET_contentate.py  PT_contentate.py 200_ET.ID  200_PE.ID  200_PT.ID ./${i};  cd ./${i}  ; python PE_contentate.py  ./ 200_PE.ID PE.fa; python PT_contentate.py  ./ 200_PT.ID PT.fa; python ET_contentate.py  ./ 200_ET.ID ET.fa; cd ../";done >contentate.sh


## mcmctree
rm  *mcmc*sh.completed
rm cp_mcmctree.sh.completed
for i in `less dir_ID.txt`;do echo " sed  's/seqfile.*/seqfile=PT.fa/;s/treefile.*/treefile=PT.tree/' mcmctree.ctl >${i}/mcmctree.ctl;"; done>cp_mcmctree.sh
~/miniconda3/bin/ParaFly  -c cp_mcmctree.sh -CPU 128
for i in `less dir_ID.txt`;do echo " sed  's/seqfile.*/seqfile=PE.fa/;s/treefile.*/treefile=PE.tree/' mcmctree.ctl >${i}/mcmctree.ctl;"; done>cp_mcmctree.sh
for i in `less dir_ID.txt`;do echo " sed  's/seqfile.*/seqfile=ET.fa/;s/treefile.*/treefile=ET.tree/' mcmctree.ctl >${i}/mcmctree.ctl;"; done>cp_mcmctree.sh

for i in `less dir_ID.txt`; do echo "cd ${i}; nohup mcmctree ; cd ../";done >mcmc_step1.sh
~/miniconda3/bin/ParaFly  -c mcmc_step1.sh -CPU 128

for i in `less dir_ID.txt`; do echo "cd ${i}; mv out.BV in.BV; sed -i 's/usedata = 3/usedata = 2/' mcmctree.ctl;  mcmctree ; cd ../";done >mcmc_step2.sh
~/miniconda3/bin/ParaFly  -c mcmc_step2.sh -CPU 128


## stats tree
for i in `less dir_ID.txt`; do echo "cp ${i}/FigTree.tre ./${i}/PT_time.tre; cp ${i}/PT_time.tre ./PT_stats/${i}.tre ";done>cp_PT.sh
for i in `less dir_ID.txt`; do echo "cp ${i}/FigTree.tre ./${i}/PE_time.tre; cp ${i}/PE_time.tre ./PE_stats/${i}.tre ";done>cp_PE.sh
for i in `less dir_ID.txt`; do echo "cp ${i}/FigTree.tre ./${i}/ET_time.tre; cp ${i}/ET_time.tre ./ET_stats/${i}.tre ";done>cp_ET.sh


cat *tre | awk -F ':' '{print $2}' |awk -F ',' '{print $1}'  |grep "." >tmp
cat *tre | awk -F ':'  '{print $5}' | awk -F ')' '{print $1}' |grep "." >tmp1
cat *tre | awk -F ':'  '{print $7}' | awk -F ')' '{print $1}' |grep "." >tmp2
paste tmp tmp1 tmp2 |less -S
```

