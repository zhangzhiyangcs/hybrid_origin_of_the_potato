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