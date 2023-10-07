### Phylogeny analysis for protein, intron, Gene sequences


Remove duplicate sequences which may interfere with phylogeny
```shell
for i in `less orthogroup_gene.txt`; do echo "python remove_duplicate_seq.py   01.data/${i}.fa > 02.data_rm_Ni_sea/${i}.fa";done >01.rm_dup.sh
```
Aligned sequences using the MAFFT software
```shell
for i in `less orthogroup_gene.txt`; do echo "mafft  --thread 5  --auto 02.data_rm_Ni_sea/${i}.fa  >03.mafft_Ni_sea/${i}.mafft";done>02.mafft.sh
```

Trim the alignment results using the trimal
```shell
for i in `less orthogroup_gene.txt`; do echo "trimal   -in 03.mafft_Ni_sea/${i}.mafft -out  04.trimal_Ni_sea/${i}.mafft -gappyout";done >03.trimal.sh
```
The corresponding tree was inferred by IQTREE software
```shell
for i in `less orthogroup_gene.txt`; do echo "iqtree -quiet -s 04.trimal_Ni_sea/${i}.mafft -m GTR -nt 4 -bb 1000 -o Outgroup  -pre 05.iqtree/${i}";done >04.iqtree.sh
```



### Phylogeny analysis for CDS


Remove duplicate sequences which may interfere with phylogeny
```shell 
for i in `less orthogroup_gene.txt`; do echo "python remove_duplicate_seq.py   01.data/${i}.fa > 02.data_rm_Ni_sea/${i}.fa";done >01.rm_dup.sh
```

Aligned sequences using the MAFFT software
```shell
for i in `less orthogroup_gene.txt`; do echo "mafft  --thread 5  --auto 02.data_rm_Ni_sea/${i}.fa  >03.mafft_Ni_sea/${i}.mafft";done>02.mafft.sh
```

CDSs were aligned according to the protein alignments using PAL2NAL
```shell
for i in `less orthogroup_gene.txt`; do echo "~/miniconda3/bin/pal2nal.pl  02.protein_mafft/${i}.mafft 01.cds_mafft/${i}.mafft  -output fasta >03.pal2nal/${i}.pal";done>03.pal2nal.sh
```

The corresponding tree was inferred by IQTREE software
```shell
for i in `less orthogroup_gene.txt`; do echo "iqtree -quiet -s 04.trimal_Ni_sea/${i}.pal -m GTR -nt 4 -bb 1000 -o Outgroup  -pre 05.iqtree/${i}";done >04.iqtree.sh
```



## Construct the window tree
Something need before analysis, and a demo can be download from this directory ([01.windowtree_data.zip](https://github.com/zhangzhiyangcs/hybrid_origin_of_the_potato/blob/main/windowtree/01.windowtree_data.zip))
```shell
Software: bedtools,seqkit,iqtree,samtools,parallel Files: chr01.fa (cactus alignments), DM_chr01.repeat.bed (repeat bed for DM), Solanum_tuberosumDM.fa (DM genome),sampleID (species name)
remove previous result
rm -r chr01 chr02 rm new_chr*bed
```

Generate repeat bed for cactus alignment

```shell
ls *sed.sh | parallel -j1 'sh ${i}'
```

Mask the repeat for cactus alignment according to bed file such as new_chr01.bed

```shell
for i in {01..02}; do echo "bedtools maskfasta -fi chr${i}.fa -fo chr${i}_repeat.fa -bed new_chr${i}.bed";done >bed_mask.sh sh bed_mask.sh
```

Generate 10-kb window for each chromosom (set -w for bedtools)

```shell
samtools faidx Solanum_tuberosumDM.fa bedtools makewindows -g <(cut -f1-2 Solanum_tuberosumDM.fa.fai) -w 10000 > DM_chr_1k_win_bed for chr in {01..02};do grep ^chr${chr} DM_chr_1k_win_bed> chr${chr}_win.id;done sed -i "1i#" *_win.id for chr in {01..02};do sed -i 's/chr'${chr}'/Solanum_tuberosumDM/g' chr${chr}_win.id;done sh generate_bed.sh
generate commands for extract sequence in chr01,chr02 file
sh seqkit_seq.sh
```

Extract sequence according to window bed

```shell
for i in {01..02}; do echo "cd chr${i}; cp ../work.sh ./; sbatch work.sh para_seqkit.sh; cd .."; done >get_sequence.sh sh get_sequence.sh
```

Change N to gap "-"

```shell
for i in {01..02}; do echo "cd chr${i}/chr${i}.fa ; ls *fa |awk '{print "sed -i 's/N/-/g' "$_}' >sed.sh; sbatch ../work.sh sed.sh; cd ../../"; done >sub_N_gap.sh sh sub_N_gap.sh
```

phylogeny

```shell
for chr in {01..02};do cd 'chr'${chr} mkdir chr${chr}.iqtree for tree in less chr${chr}.list; do echo "iqtree -quiet -s chr${chr}.fa/${tree}.fa -m GTR -nt 4 -bb 1000 -o Solanum_incanum -pre chr${chr}.iqtree/${tree}" >>iqtree.sh done sbatch work.sh iqtree.sh cd .. done
```

Use the msa_view software to filter cactus results (alternatively steps)

```shell
for chr in {01..12};do cd 'chr'${chr} mkdir chr${chr}.align for bed in less chr${chr}.list; do echo "msa_view -G ANY chr${chr}.fa/${bed}.fa |sed 's/ //'>chr${chr}.align/${bed}.fa" >>msa_view.sh done sbatch work.sh msa_view.sh cd .. done
```



## Coalescent-based method for phylogeny

Astral method
```shell
java -jar  ~/software/Astral/astral.5.15.5.jar  -i total.tree  -o gene_astral -C -T 128
```

MP-EST method

```shell
~/software/mp-est-master/linux/mpest2.0   control
```



## Simulate solely ILS scenario

simulate tree according to the astral result via DendroPy package
```shell
python3 02.simulate_tree-2.py simulated.trees.2; mv simulated.trees.2 bootstrap/
perl 04.bootstrap.pl  >simulate.sh
```

Stats the topology of simulated tree
```shell
python count_topology.py sample.ID simulated.trees.2 simulated.trees.2.txt
```




## Stats the tree

```R
library(ape)
library(phangorn)

args <- commandArgs(trailingOnly = TRUE)
e <- read.tree(args[1])
sample_ID <- read.csv("sample_ID.txt", header = FALSE, sep = "\t")

ETB_clade <- sample_ID[sample_ID$V2 == "ETB", "V1"]
Tomato_clade <- sample_ID[sample_ID$V2 == "Tomato", "V1"]
Potato_clade <- sample_ID[sample_ID$V2 == "Potato", "V1"]
Outgroup_clade <- sample_ID[sample_ID$V2 == "Outgroup", "V1"]

CLADE1 <- ETB_clade
CLADE2 <- Tomato_clade
CLADE3 <- Potato_clade
CLADE12 <- c(ETB_clade, Tomato_clade)
CLADE13 <- c(ETB_clade, Potato_clade)
CLADE23 <- c(Tomato_clade, Potato_clade)
CLADE_out <- Outgroup_clade

cla1 <- cla2 <- cla3 <- cla12 <- cla13 <- cla23 <- all <- other <- 0

for (gtree in e) {
  gtree$edge.length <- NULL
  if (length(gtree$tip) != 40) {
    print("no_40")
    next
  }
  all <- all + 1
  gtree <- root(gtree, "Outgroup", node = NULL, resolve.root = FALSE, interactive = FALSE)
  if (is.monophyletic(gtree, CLADE13) && is.monophyletic(gtree, CLADE2)) {
    cla13 <- cla13 + 1
    print("PE ")
    next
  }
  if (is.monophyletic(gtree, CLADE12) && is.monophyletic(gtree, CLADE3)) {
    cla12 <- cla12 + 1
    print("ET ")
    next
  }
  if (is.monophyletic(gtree, CLADE23) && is.monophyletic(gtree, CLADE1)) {
    cla23 <- cla23 + 1
    print("PT ")
    next
  }
  else {
    print("Other")
    next
  }
}
print("ETB_Tomato_clade")
print(cla12)
print("ETB_Potato_clade")
print(cla13)
print("Tomato_Potato_clade")
print(cla23)
```



