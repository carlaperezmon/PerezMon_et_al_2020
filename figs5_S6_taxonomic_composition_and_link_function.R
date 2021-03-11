# load packages needed
ipak <- function(pkg) {
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) {
    install.packages(new.pkg, dependencies = TRUE)
  }
  sapply(pkg, require, character.only = TRUE)
}


packages <- c(
  "stringr", "reshape2", "digest", "AnnotationDbi", "ggplot2", "scales", "GO.db", "spaa", "xlsx", "FactoMineR", "factoextra", "xlsx", "dplyr", "zoo", "ggplot2", "gplots", "reshape2", "RColorBrewer",
  "colorspace", "colorRamps", "vegan", "gdata", "reshape", "stats", "cluster", "Rmisc", "cowplot",
  "indicspecies", "userfriendlyscience", "onewaytests", "RVAideMemoire", "agricolae", "scales",
  "pairwiseAdonis", "EcolUtils", "stringr", "PMCMRplus", "PMCMR", "ape", "Rmisc", "yarr", "pals"
)


ipak(packages)

path_data <- "F:/metagenomics_analyses_fgcz/files_user_interface/statistical_analyses/results_own/composition/"
setwd(path_data)


##### LOAD DATA #######

# for otu
otu <- read.csv2("arb-silva.de_align_resultlist_edited_silva.csv", header = T, sep = ";", check.names = FALSE)


# load genes_length for the TPM normalization

genes_length <- read.csv2("F:/metagenomics_analyses_fgcz/files_user_interface/statistical_analyses/results_own/genes_length.txt", header = TRUE, sep = "\t")


# load identity of samples
levels_soils <- c("N10", "N160", "S10", "S160")
metadata <- read.xlsx("metadata.xlsx", sheetName = "Sheet1", header = TRUE)
rownames(metadata) <- metadata[, "id"]
metadata <- metadata[order(factor(metadata$names, levels = levels_soils)), ]
metadata$id2 <- paste(metadata$id, "_function", sep = "")


# for depth of contigs
depth <- read.csv2("persample.depth.txt", header = T, sep = "\t", check.names = FALSE) #### for the depth of coverage of contigs
colnames(depth)[1] <- "contig"


otu_depth <- merge(otu, depth, by = c("contig")) # 11 contigs are lost
otu_depth_sub <- otu_depth[otu_depth$gene_bps > 300, ] # select a minimum size of bp for the classification, in this case 300bp, which is the minimal lenght of SSU in SILVA


otu_depth_sub <- otu_depth_sub[, -which(names(otu_depth_sub) %in% c(
  "job_id", "sequence_number",
  "cutoff_head", "cutoff_tail",
  "quality", "startpos", "stoppos", "ecolipos",
  "bps", "gene_bps", "turn"
))]


########### CREATION MATRICES BACT, ARCH AND EUK AND REPRESENTATION #####


# separate the classification

classi <- as.data.frame(str_split(otu_depth_sub$lca_tax_slv, ";", simplify = TRUE))
names(classi) <- c("domain", "phylum", "class", "order", "family", "genus")

classi <- classi[, -14]

classi <- apply(classi, 2, function(x) as.character(x))
classi[nchar(classi) == 0] <- "unclassified"


# create the classification table for representation at the contig level

otu_contig <- cbind(classi, otu_depth_sub[, 6:ncol(otu_depth_sub)])

# transform to relative abundance: depth of coverage of taxa/depth of coverage of all contigs classified at least to kingdom

# now do relative abundance for everything

otu_contig[, c(14:ncol(otu_contig))] <- apply(otu_contig[, c(14:ncol(otu_contig))], 2, function(x) as.numeric(as.character(x)))
otu_contig[, c(14:ncol(otu_contig))] <- apply(otu_contig[, c(14:ncol(otu_contig))], 2, function(x) x / sum(x))

colSums(otu_contig[, c(14:ncol(otu_contig))])


####### create the tables of % abundance and richness for the most common bacterial phyla

selection <- c("abundance", "richess")

for (i in 1:length(selection)) {
  print(selection[i])

  otu_contig_bac_phy <- otu_contig[otu_contig$domain == "Bacteria", ]
  otu_contig_bac_phy <- otu_contig_bac_phy[, which(names(otu_contig_bac_phy) %in% c("phylum", as.character(metadata$id)))]

  if (i == 2) {
    otu_contig_bac_phy[, -which(names(otu_contig_bac_phy) %in% c("phylum"))] <- +(otu_contig_bac_phy[, -which(names(otu_contig_bac_phy) %in% c("phylum"))] > 0)
  }

  otu_contig_bac_phy <- aggregate(. ~ phylum, otu_contig_bac_phy, "sum")

  otu_contig_bac_phy_melt <- melt(otu_contig_bac_phy)
  otu_contig_bac_phy_melt$variable <- metadata$names[match(otu_contig_bac_phy_melt$variable, metadata$id)]

  otu_contig_bac_phy_summary <- summarySE(otu_contig_bac_phy_melt, measurevar = "value", groupvars = c("variable", "phylum"))
  write.xlsx(otu_contig_bac_phy_summary, paste(selection[i], "_phyla_summary_abundances.xls", sep = ""))


  otu_contig_bac_phy$sum <- rowSums(otu_contig_bac_phy[, -1])
  otu_contig_bac_phy$rel <- otu_contig_bac_phy$sum / sum(otu_contig_bac_phy$sum)

  write.xlsx(otu_contig_bac_phy, paste(selection[i], "_phyla_total_abundances.xls", sep = ""))
}



# calculate the relative abundance of domains for Table 1

otu_contig_per <- otu_contig
otu_contig_per <- otu_contig_per[, c(1, 14:ncol(otu_contig_per))]
otu_contig_per <- otu_contig_per[otu_contig_per$domain != "Unclassified", ]

otu_contig_per <- aggregate(. ~ domain, otu_contig_per, "sum") # since I want to access % of total reads I sum
otu_contig_per[, -1] <- otu_contig_per[, -1] * 100
colSums(otu_contig_per[, -1])

otu_contig_per <- melt(otu_contig_per, id.vars = "domain")
otu_contig_per$variable <- metadata$names[match(otu_contig_per$variable, as.character(metadata$id))]

otu_contig_domain <- summarySE(otu_contig_per, measurevar = "value", groupvars = c("variable", "domain"))

write.xlsx(otu_contig_domain, "domain_abundances.xls")


####### REPRESENT TAXONOMY (fig. S6) ######


soils_color <- c("N160" = "steelblue4", "N10" = "darkorchid1", "S160" = "springgreen4", "S10" = "coral4")

# represent
theme_rep <- theme(
  panel.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.text = element_text(size = 30),
  axis.text.y = element_text(size = 15, color = "black"),
  axis.text.x = element_text(size = 15, color = "black", angle = 0, hjust = 0.95, vjust = 0.2),
  axis.title = element_text(size = 20),
  legend.text = element_text(size = 10),
  legend.title = element_text(size = 10),
  axis.line = element_line(colour = "black", size = 1, linetype = 1),
  axis.ticks = element_line(colour = "black", size = 1),
  axis.ticks.length = unit(0.2, "cm"),
  panel.spacing = unit(3, "lines"),
  # legend.position = "none",
  # legend.direction = 'horizontal',
  legend.key = element_rect(size = 2),
  legend.key.size = unit(1, "lines"),
  axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
  axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0))
)



# Prokaryotes

otu_contig_prok <- otu_contig

otu_contig_prok$domain <- paste("d", otu_contig_prok$domain, sep = "__")
otu_contig_prok$phylum <- paste("p", otu_contig_prok$phylum, sep = "__")
otu_contig_prok$class <- paste("c", otu_contig_prok$class, sep = "__")
otu_contig_prok$order <- paste("o", otu_contig_prok$order, sep = "__")
otu_contig_prok$family <- paste("f", otu_contig_prok$family, sep = "__")
otu_contig_prok$genus <- paste("g", otu_contig_prok$genus, sep = "__")


otu_contig_prok[, c(1:13)] <- apply(otu_contig_prok[, c(1:13)], 2, function(x) as.character(x))


# classify to the deepest taxonomical level possible

for (i in 2:ncol(otu_contig_prok[, c(1:13)])) {
  replacement <- grep("[a-z]__unclassified", otu_contig_prok[, i], value = TRUE)
  for (j in 1:length(replacement)) {
    otu_contig_prok[, i][otu_contig_prok[, i] == paste(replacement[j])] <- paste(otu_contig_prok[, i - 1][otu_contig_prok[, i] == paste(replacement[j])])
  }
}

for (i in 1:ncol(otu_contig_prok[, c(13:ncol(otu_contig_prok))])) {
  replacement <- grep("[a-z]__uncultured", otu_contig_prok[, i], value = TRUE)
  for (j in 1:length(replacement)) {
    otu_contig_prok[, i][otu_contig_prok[, i] == paste(replacement[j])] <- paste(otu_contig_prok[, i - 1][otu_contig_prok[, i] == paste(replacement[j])])
  }
}

otu_contig_prok[, c(1:13)] <- lapply(otu_contig_prok[, c(1:13)], factor)

tree <- as.phylo.formula(~ domain / phylum / class / order / family / genus, data = otu_contig_prok)


# Separate bacteria and archaea

# level of bacteria

otu_contig_bac <- otu_contig_prok[otu_contig_prok$domain == "d__Bacteria", ]

# Shannon diversity of bacteria

otu_contig_bac_div <- otu_contig_bac[, which(names(otu_contig_bac) %in% c(as.character(metadata$id)))]
shannon_bac <- t(otu_contig_bac_div)
shannon_bac_div <- cbind.data.frame(variable = row.names(shannon_bac), value = diversity(shannon_bac, index = "shannon"))
shannon_bac_div$variable <- as.character(metadata$names)[match(shannon_bac_div$variable, as.character(metadata$id))]
bac_div_summary <- summarySE(shannon_bac_div, measurevar = "value", groupvars = c("variable"))


# abundance of bacterial taxa classified at deepest taxonomical level possible


otu_contig_bac <- otu_contig_bac[, which(names(otu_contig_bac) %in% c("genus", as.character(metadata$id)))]
otu_contig_bac <- otu_contig_bac[otu_contig_bac$genus != "d__Bacteria", ]

otu_contig_bac <- aggregate(. ~ genus, otu_contig_bac, "sum")

otu_contig_bac$sum <- rowSums(otu_contig_bac[, -1])
otu_contig_bac$rel <- rowSums(otu_contig_bac[, -1]) / sum(otu_contig_bac$sum)


otu_contig_bac <- otu_contig_bac[otu_contig_bac$rel > 0.02, ] # filter for 0.02 relative abundance for the total across all soil samples

otu_contig_bac <- otu_contig_bac[, which(names(otu_contig_bac) %in% c("genus", as.character(metadata$id)))]

otu_contig_bac <- melt(otu_contig_bac, id.vars = "genus")
otu_contig_bac$genus <- factor(otu_contig_bac$genus)
otu_contig_bac$variable <- metadata$names[match(otu_contig_bac$variable, as.character(metadata$id))]


# graphical representation

labl <- unique(tree$tip.label[tree$tip.label %in% otu_contig_bac$genus])

# generate colors
{ # color contributor
  c25 <- c(
    "dodgerblue2", "#E31A1C", # red
    "green4",
    "#6A3D9A", # purple
    "#FF7F00", # orange
    "black", "gold1",
    "skyblue2", "#FB9A99", # lt pink
    "palegreen2",
    "#CAB2D6", # lt purple
    "#FDBF6F", # lt orange
    "gray70", "khaki2",
    "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
    "darkturquoise", "green1", "yellow4", "yellow3",
    "darkorange4", "brown"
  )
}

piratepal("all") # package yarr

pal.bands(polychrome) # package pals

colors_own <- unname(c(
  piratepal("basel"), piratepal("pony"), piratepal("eternal"),
  piratepal("rat"), piratepal("appletv"), piratepal("cars"),
  c25, primary.colors(10)
))

colors_own_ran <- sample(colors_own, length(unique(otu_contig_bac$genus)))


colors_own_ran <- as.vector(glasbey(length(unique(otu_contig_bac$genus))))

otu_contig_bac_sum <- summarySE(otu_contig_bac, measurevar = "value", groupvars = c("variable", "genus"))


ggplot(otu_contig_bac_sum, aes(factor(variable, levels = levels_soils), value, fill = ordered(genus, levels = labl))) +
  geom_histogram(stat = "identity", position = "fill", colour = "black", size = 0.05) +
  scale_fill_manual(values = colors_own_ran) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab(" ") +
  ylab("Relative abundance") +
  theme_rep +
  guides(fill = guide_legend(title = "Taxa", ncol = 1))






# Archaea

otu_contig_arch <- otu_contig_prok[otu_contig_prok$domain == "d__Archaea", ]


# Shannon diversity of Archaea

otu_contig_arch_div <- otu_contig_arch[, which(names(otu_contig_arch) %in% c(as.character(metadata$id)))]
shannon_arch <- t(otu_contig_arch_div)
shannon_arch_div <- cbind.data.frame(variable = row.names(shannon_arch), value = diversity(shannon_arch, index = "shannon"))
shannon_arch_div$variable <- as.character(metadata$names)[match(shannon_arch_div$variable, as.character(metadata$id))]
arch_div_summary <- summarySE(shannon_arch_div, measurevar = "value", groupvars = c("variable"))


# abundance of archaeal taxa classified at deepest taxonomical level possible

otu_contig_arch <- otu_contig_arch[, which(names(otu_contig_arch) %in% c("genus", as.character(metadata$id)))]
otu_contig_arch <- otu_contig_arch[otu_contig_arch$genus != "d__Archaea", ]

otu_contig_arch <- aggregate(. ~ genus, otu_contig_arch, "sum")

otu_contig_arch <- melt(otu_contig_arch, id.vars = "genus")
otu_contig_arch$genus <- factor(otu_contig_arch$genus)
otu_contig_arch$variable <- metadata$names[match(otu_contig_arch$variable, as.character(metadata$id))]


otu_contig_arch_sum <- summarySE(otu_contig_arch, measurevar = "value", groupvars = c("variable", "genus"))


# graphical representation

labl <- unique(tree$tip.label[tree$tip.label %in% otu_contig_arch$genus])

colors_own_ran <- as.vector(watlington(length(unique(otu_contig_arch$genus)))) # for now I decided in a pallete of pals


ggplot(otu_contig_arch_sum, aes(factor(variable, levels = levels_soils), value, fill = ordered(genus, levels = labl))) +
  geom_histogram(stat = "identity", position = "fill", colour = "black", size = 0.05) +
  scale_fill_manual(values = colors_own_ran) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab(" ") +
  ylab("Relative abundance") +
  theme_rep +
  guides(fill = guide_legend(title = "Taxa", ncol = 1))



# Eukaryotes

# redo the first steps of the otu_contig

otu_contig_euk <- otu_contig[otu_contig$domain == "Eukaryota", ]


# Shannon diversity of eukaryotic taxa

otu_contig_euk_div <- otu_contig_euk[, which(names(otu_contig_euk) %in% c(as.character(metadata$id)))]
shannon_euk <- t(otu_contig_euk_div)
shannon_euk_div <- cbind.data.frame(variable = row.names(shannon_euk), value = diversity(shannon_euk, index = "shannon"))
shannon_euk_div$variable <- as.character(metadata$names)[match(shannon_euk_div$variable, as.character(metadata$id))]
euk_div_summary <- summarySE(shannon_euk_div, measurevar = "value", groupvars = c("variable"))


# rearrengement of tables for graphical representations

names(otu_contig_euk)[1:13] <- c("domain", "kingdom", "phylum", "subphylum", "class", "subclass", "superorder", "order", "suborder", "family", "subfamily", "genus", "species")

otu_contig_euk[, c(1:13)] <- apply(otu_contig_euk[, c(1:13)], 2, function(x) as.character(x))
otu_contig_euk <- na.omit(otu_contig_euk)

# for substitution of unclassified and uncultured by known taxa
for (i in 2:ncol(otu_contig_euk[, c(1:13)])) {
  replacement <- grep("unclassified", otu_contig_euk[, i], value = TRUE)
  for (j in 1:length(replacement)) {
    otu_contig_euk[, i][otu_contig_euk[, i] == paste(replacement[j])] <- paste(otu_contig_euk[, i - 1][otu_contig_euk[, i] == paste(replacement[j])])
  }
}

for (i in 1:ncol(otu_contig_euk[, c(13:ncol(otu_contig_euk))])) {
  replacement <- grep("uncultured", otu_contig_euk[, i], value = TRUE)
  for (j in 1:length(replacement)) {
    otu_contig_euk[, i][otu_contig_euk[, i] == paste(replacement[j])] <- paste(otu_contig_euk[, i - 1][otu_contig_euk[, i] == paste(replacement[j])])
  }
}

otu_contig_euk[, c(1:13)] <- lapply(otu_contig_euk[, c(1:13)], factor)

tree <- as.phylo.formula(~ domain / kingdom / phylum / subphylum / class / subclass / superorder / order / suborder / family / subfamily / genus / species, data = otu_contig_euk)


otu_contig_euk <- otu_contig_euk[, 13:25]
names(otu_contig_euk)[1] <- "genus"


otu_contig_euk <- aggregate(. ~ genus, otu_contig_euk, "sum")


otu_contig_euk <- melt(otu_contig_euk, id.vars = "genus")
otu_contig_euk$genus <- factor(otu_contig_euk$genus)
otu_contig_euk$variable <- metadata$names[match(otu_contig_euk$variable, as.character(metadata$id))]

otu_contig_euk_sum <- summarySE(otu_contig_euk, measurevar = "value", groupvars = c("variable", "genus"))


# graphical representation

labl <- unique(tree$tip.label[tree$tip.label %in% otu_contig_euk$genus])

colors_own_ran <- as.vector(polychrome(length(unique(otu_contig_euk$genus)))) # for now I decided in a pallete of pals

ggplot(otu_contig_euk_sum, aes(factor(variable, levels = levels_soils), value, fill = ordered(genus, levels = labl))) +
  geom_histogram(stat = "identity", position = "fill", colour = "black", size = 0.05) +
  scale_fill_manual(values = colors_own_ran) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab(" ") +
  ylab("Relative abundance") +
  theme_rep +
  guides(fill = guide_legend(title = "Taxa", ncol = 1))





####### LINK MICROBES TO FUNCTION ####


### MICROBES LINK TO CAZY AND NCYC DATABASE ####

# load table mapping of genes into classified contigs

genes_contig <- read.csv2("genes_annotation_byTranscript.txt", header = T, sep = "\t", check.names = FALSE)

genes_contig2 <- genes_contig[(genes_contig$seqid %in% otu_depth_sub$contig), ]
genes_contig2 <- genes_contig2[, which(names(genes_contig2) %in% c("gene_id", "seqid"))]

genes_contig3 <- aggregate(gene_id ~ seqid, genes_contig2, paste, collapse = ",")
names(genes_contig3)[1] <- "contig"


otu_depth_sub_contig <- merge(otu_depth_sub, genes_contig3, by = "contig")
otu_depth_sub_contig2 <- separate_rows(otu_depth_sub_contig, gene_id, sep = ",")
names(otu_depth_sub_contig2)[ncol(otu_depth_sub_contig2)] <- "Identifier"



# load cazy and ncyc db

dbs <- c("cazy", "ncyc")
path_tax <- c("F:/metagenomics_analyses_fgcz/files_user_interface/cazy/aggr/CAZyDB.aggr.anno.S1.txt", "F:/metagenomics_analyses_fgcz/files_user_interface/ncyc/aggr/NCyc.aggr.anno.S1.txt")


# taxonomy
for (i in 1:length(path_tax)) {
  tax <- read.csv2(path_tax[i], header = TRUE, sep = "\t")
  colnames(tax) <- gsub(".\\.", "", colnames(tax))

  tax <- cbind.data.frame(FunID = tax[, c("FunID")], str_split(tax$taxonomy, ";", simplify = TRUE))
  colnames(tax)[-1] <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")

  tax <- as.data.frame(apply(tax, 2, function(x) gsub("[a-z]__", "", x)))
  tax <- tax[, colSums(tax != " ") != 0]

  assign(paste(dbs[i], "tax", sep = "_"), tax)
  rm(tax)
}



# raw gene counts with taxonomy
for (i in 1:length(dbs)) {
  print(dbs[i])
  get(paste(dbs[i], "tax", sep = "_")) -> tax

  path_files <- paste("F:/metagenomics_analyses_fgcz/files_user_interface/", dbs[i], sep = "")

  files <- list.files(paste(path_files, "/anno/", sep = ""), pattern = ".TXT", ignore.case = T)
  names <- str_split(files, "\\.", simplify = TRUE)[, 3]

  counts <- read.csv2(paste(path_files, "/anno/", files[1], sep = ""), header = TRUE, sep = "\t")
  counts <- counts[, c("Identifier", "FunID", "matchCounts")] # reorder columns

  for (j in 2:length(files)) {
    tb <- read.csv2(paste(path_files, "/anno/", files[j], sep = ""), header = TRUE, sep = "\t")

    counts$m2 <- tb$matchCounts[match(counts$Identifier, tb$Identifier)]
    colnames(counts) <- c(names(counts)[1:2], names[1:j])
  }

  counts_tax <- merge(counts, tax, by = "FunID", all = FALSE)

  counts_tax$class <- gsub(" ", "", counts_tax$class)

  assign(paste(dbs[i], "counts_tax", sep = "_"), counts_tax)
  rm(counts, counts_tax)
}



# filtering out genes with counts < 10 (standard of Deseq analyses)

cazy_ident <- cazy_counts_tax$Identifier[rowSums(cazy_counts_tax[, which(names(cazy_counts_tax) %in% as.character(metadata$id))]) > 9]
length(cazy_ident)

ncyc_ident <- ncyc_counts_tax$Identifier[rowSums(ncyc_counts_tax[, which(names(ncyc_counts_tax) %in% as.character(metadata$id))]) > 9]
length(ncyc_ident)


cazy_counts_tax <- na.omit(cazy_counts_tax[cazy_counts_tax$Identifier %in% cazy_ident, ])
cazy_counts_tax$class2 <- paste(cazy_counts_tax$class, cazy_counts_tax$order)
cazy_counts_tax$class2 <- trimws(cazy_counts_tax$class2, "b")
cazy_counts_tax$phylum2 <- cazy_counts_tax$phylum
cazy_counts_tax$phylum2 <- trimws(cazy_counts_tax$phylum2, "b")

ncyc_counts_tax <- na.omit(ncyc_counts_tax[ncyc_counts_tax$Identifier %in% ncyc_ident, ])
ncyc_counts_tax$class2 <- ncyc_counts_tax$class
ncyc_counts_tax$class2 <- trimws(ncyc_counts_tax$class2, "b")
ncyc_counts_tax$phylum2 <- ncyc_counts_tax$phylum
ncyc_counts_tax$phylum2 <- trimws(ncyc_counts_tax$phylum2, "b")




# create matrices of normalized TPM counts

for (i in 1:length(dbs)) {
  print(dbs[i])
  get(paste(dbs[i], "counts_tax", sep = "_")) -> tpm_tax
  get(paste(dbs[i], "gene_length", sep = "_")) -> gene_length

  tpm_tax$length <- gene_length$Length[match(tpm_tax$Identifier, gene_length$Geneid)]

  tpm_tax[, which(names(tpm_tax) %in% as.character(metadata$id))] <- tpm_tax[, which(names(tpm_tax) %in% as.character(metadata$id))] / tpm_tax$length # divide by length of gene

  # tpm: I do it with apply function because somehow is the only way in which i get same tpm counts as in count qc
  tpm_tax[, which(names(tpm_tax) %in% as.character(metadata$id))] <- apply(tpm_tax[, which(names(tpm_tax) %in% as.character(metadata$id))], 2, function(x) x / sum(x) * 10^6)

  colSums(tpm_tax[, which(names(tpm_tax) %in% as.character(metadata$id))])

  assign(paste(dbs[i], "tpm_tax", sep = "_"), tpm_tax)
}



# graphical representation of fig 5 (a and b)

palletes <- c("Tofino", "Berlin")


for (i in 1:length(dbs)) {
  print(dbs[i])
  get(paste(dbs[i], "tpm_tax", sep = "_")) -> tpm_tax

  tpm_tax <- tpm_tax[, which(names(tpm_tax) %in% c("Identifier", paste(metadata$id), "phylum2", "class2", "FunID"))]
  tpm_tax <- tpm_tax[tpm_tax$phylum != "GlycosylTransferases", ]

  names(tpm_tax)[3:14] <- metadata$id2[match(names(tpm_tax)[3:14], rownames(metadata))]

  tpm_tax2 <- merge(tpm_tax, otu_depth_sub_contig2, by = "Identifier")

  check_tax <- tpm_tax2[, which(names(tpm_tax2) %in% c("Identifier", "FunID", "lca_tax_slv"))]

  # trial=otu_depth_sub_contig2
  # trial$ncyc=tpm_tax$class2[match(trial$Identifier, tpm_tax$Identifier)]

  classi2 <- as.data.frame(str_split(tpm_tax2$lca_tax_slv, ";", simplify = TRUE))
  names(classi2) <- c("domain", "phylum", "class", "order", "family", "genus")

  classi2 <- classi2[, -ncol(classi2)]

  classi2 <- apply(classi2, 2, function(x) as.character(x))
  classi2[nchar(classi2) == 0] <- "unclassified"

  tpm_tax3 <- tpm_tax2[, -which(names(tpm_tax2) %in% c("contig", "sequence_identifier", "bp_score", "identity", "lca_tax_slv"))]
  tpm_tax3 <- cbind(classi2, tpm_tax3)

  # have graphs united
  tpm_tax3 <- tpm_tax3[tpm_tax3$phylum != "unclassified", ]


  # spliting the graphs

  tax_function <- tpm_tax3[, which(names(tpm_tax3) %in% c(
    "domain", "kingdom", "phylum", "class", "order", "family", "genus", "species",
    paste(metadata$id2), "class2", "phylum2"
  ))]

  tax_function2 <- tpm_tax3[, which(names(tpm_tax3) %in% c("phylum", paste(metadata$id2), "class2"))]

  tax_function_agg <- aggregate(. ~ phylum + class2, tax_function2, "sum")
  tax_function_agg$total <- rowSums(tax_function_agg[, which(names(tax_function_agg) %in% paste(metadata$id2))])

  rep <- ggplot(tax_function_agg, aes(phylum, total, fill = class2)) +
    # geom_errorbar(aes(ymin=value-sd, ymax=value+sd),stat = "identity", position =position_dodge(0.9), colour="black", width=0.5) +
    geom_col(position = "stack", color = "black", width = 0.8) +
    # scale_fill_discrete_qualitative(palette="Dark 2")+
    scale_fill_discrete_diverging(palette = palletes[i]) +
    geom_text(aes(label = class2), color = "white", size = 3, position = position_stack(vjust = 0.5)) +
    # scale_fill_brewer(palette = 3)+
    xlab("") +
    ylab("abundance [tpm]") +
    # ggtitle(paste(dbs[i]))+
    # ylab("richness") +
    theme_rep +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, vjust = 1)
    )

  print(rep)


  tax_abund <- tpm_tax3[, which(names(tpm_tax3) %in% c(
    "domain", "kingdom", "phylum", "class", "order", "family", "genus", "species",
    paste(metadata$id)
  ))]


  tax_abund2 <- tpm_tax3[, which(names(tpm_tax3) %in% c(
    "phylum",
    paste(metadata$id)
  ))]

  tax_abund2 <- unique(tax_abund2)
  tax_abund2[, -1] <- apply(tax_abund2[, -1], 2, function(x) as.numeric(x))

  tax_abund_agg <- aggregate(. ~ phylum, tax_abund2, "sum")

  assign(paste(dbs[i], "tpm_tax_abund_function", sep = "_"), tpm_tax3)
  assign(paste(dbs[i], "tpm_tax_function", sep = "_"), tax_function_agg)
  assign(paste(dbs[i], "tpm_tax_abund", sep = "_"), tax_abund_agg)

  rm(tpm_tax3, tax_function, tax_abund, tax_abund_agg)
}



## graphical representation of fig 5 (c)

#### abundances of selected phyla accross the soils

otu_contig_bac_phy_rel <- otu_contig_bac_phy
otu_contig_bac_phy_rel[, which(names(otu_contig_bac_phy_rel) %in% as.character(metadata$id))] <- otu_contig_bac_phy_rel[, which(names(otu_contig_bac_phy_rel) %in% as.character(metadata$id))] / colSums(otu_contig_bac_phy_rel[, which(names(otu_contig_bac_phy_rel) %in% as.character(metadata$id))])

otu_contig_bac_phy_rel <- otu_contig_bac_phy_rel[, -which(names(otu_contig_bac_phy_rel) %in% c("sum"))]

otu_contig_bac_phy_rel_melt <- melt(otu_contig_bac_phy_rel)
otu_contig_bac_phy_rel_melt$variable <- metadata$names[match(otu_contig_bac_phy_rel_melt$variable, metadata$id)]

otu_contig_bac_phy_rel_summary <- summarySE(otu_contig_bac_phy_rel_melt, measurevar = "value", groupvars = c("variable", "phylum"))


all_taxa <- unique(c(as.character(cazy_tpm_tax_abund$phylum), as.character(ncyc_tpm_tax_abund$phylum)))
all_taxa2 <- otu_contig_bac_phy_rel_summary[otu_contig_bac_phy_rel_summary$phylum %in% all_taxa, ]


rep <- ggplot(all_taxa2, aes(factor(variable, levels = levels_soils), factor(phylum, levels = rev(unique(all_taxa2$phylum))))) +
  geom_tile(aes(fill = value, width = 1), color = "white") +
  # scale_fill_viridis_c(limits=c(min(all_taxa_melt$value),max(all_taxa_melt$value), oob = scales::squish))+
  # scale_fill_gradient2(low = "green",high = "red", na.value = "grey50", midpoint =median(all_taxa_melt$value)) +
  scale_fill_gradient2(low = "slateblue4", high = "orange3", mid = "slateblue4") +
  # scale_fill_fermenter()+
  # coord_flip()+
  scale_y_discrete(expand = c(0, 0)) +
  ylab("") +
  xlab("") +
  scale_x_discrete(position = "top") +
  # ggtitle(paste(tables[i])) +
  # labs(fill = paste(names[j])) +
  theme_rep +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, vjust = 1)
  )

print(rep)
