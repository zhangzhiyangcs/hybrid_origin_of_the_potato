## Something need before analysis
Software: bedtools,seqkit,iqtree,samtools,parallel
Files: chr01.fa (cactus alignments), DM_chr01.repeat.bed (repeat bed for DM), Solanum_tuberosumDM.fa (DM genome),sampleID (species name)

### remove previous result 
rm -r chr01 chr02
rm new_chr*bed

### generate repeat bed for cactus alignment
ls *sed.sh | parallel -j1 'sh ${i}'

### mask the repeat for cactus alignment according to bed file such as new_chr01.bed
for i in {01..02}; do echo "bedtools maskfasta -fi chr${i}.fa -fo chr${i}_repeat.fa -bed new_chr${i}.bed";done >bed_mask.sh
sh bed_mask.sh

### generate 10kb window for each chromosom (set -w for bedtools)
samtools faidx Solanum_tuberosumDM.fa
bedtools  makewindows -g <(cut -f1-2 Solanum_tuberosumDM.fa.fai) -w 10000 > DM_chr_1k_win_bed
for chr in {01..02};do grep ^chr${chr} DM_chr_1k_win_bed> chr${chr}_win.id;done
sed -i "1i#" *_win.id
for chr in {01..02};do sed -i 's/chr'${chr}'/Solanum_tuberosumDM/g' chr${chr}_win.id;done
sh generate_bed.sh

### generate commands for extract sequence in chr01,chr02 file
sh seqkit_seq.sh

### extract sequence according to window bed
for  i in {01..02}; do echo "cd chr${i}; cp ../work.sh ./; sbatch work.sh para_seqkit.sh; cd .."; done >get_sequence.sh
sh get_sequence.sh

### change N to gap "-"
for  i in {01..02}; do echo "cd chr${i}/chr${i}.fa ; ls *fa |awk '{print \"sed -i 's/N/-/g' \"\$_}' >sed.sh; sbatch ../work.sh sed.sh; cd ../../"; done >sub_N_gap.sh
sh sub_N_gap.sh

### phylogeny 
for chr in {01..02};do
cd 'chr'${chr}
mkdir chr${chr}.iqtree
for tree in `less chr${chr}.list`; do
echo "iqtree -quiet -s chr${chr}.fa/${tree}.fa -m GTR -nt 4 -bb 1000 -o Solanum_incanum  -pre chr${chr}.iqtree/${tree}" >>iqtree.sh
done
sbatch work.sh iqtree.sh
cd ..
done

### use the msa_view software to filter cactus results (alternatively steps)
for chr in {01..12};do
cd 'chr'${chr}
mkdir chr${chr}.align 
for bed in `less chr${chr}.list`; do
echo "msa_view -G ANY chr${chr}.fa/${bed}.fa |sed 's/ //'>chr${chr}.align/${bed}.fa" >>msa_view.sh
done
sbatch work.sh msa_view.sh
cd ..
done
