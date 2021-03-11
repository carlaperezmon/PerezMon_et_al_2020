# to have the messages and labels in english
Sys.setlocale(category = "LC_ALL", locale = "english")


#### Genomic data analyses ####

# load packages needed
ipak <- function(pkg) {
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}


packages <- c("stringr", "reshape2", "digest", "ggplot2", "scales", "spaa", "xlsx", "viridis", "DESeq2", "Rmisc", "onewaytests", "userfriendlyscience", "vegan")

ipak(packages)


###########################
#########################   ARRANGE TABLES AND DATA ###############################
##########################


##### load data and create independent tables

path_data <- "F:/metagenomics_analyses_fgcz/files_user_interface/statistical_analyses/results_own/deseq/cazy_ncyc/"
path_fig <- "F:/metagenomics_analyses_fgcz/files_user_interface/statistical_analyses/results_own/deseq/graphs/cazy_ncyc/"

# path_data="Z:/bowiss/RPR/Carla/experiments/field/villum/villum_2019/transplants/light_sensors_transplant_villum_2019/analysed_data/"

setwd(path_data)

# for identity of samples
metadata <- read.xlsx("metadata.xlsx", sheetName = "Sheet1", header = TRUE)
rownames(metadata) <- metadata[, "id"]
metadata <- metadata[order(factor(metadata$names, levels = levels_soils)), ]



# prepare the cazy and ncyc count tables with classification of "classes" and "families"
dbs <- c("cazy", "ncyc")
path_tax <- c("F:/metagenomics_analyses_fgcz/files_user_interface/cazy/aggr/CAZyDB.aggr.anno.S1.txt", "F:/metagenomics_analyses_fgcz/files_user_interface/ncyc/aggr/NCyc.aggr.anno.S1.txt")


# load the classes and family identities
for (i in 1:length(path_tax)) {
  print(path_tax[i])
  tax <- read.csv2(path_tax[i], header = TRUE, sep = "\t")
  colnames(tax) <- gsub(".\\.", "", colnames(tax))

  tax <- cbind.data.frame(FunID = tax[, c("FunID")], str_split(tax$taxonomy, ";", simplify = TRUE))
  colnames(tax)[-1] <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")

  tax <- as.data.frame(apply(tax, 2, function(x) gsub("[a-z]__", "", x)))
  tax <- tax[, colSums(tax != " ") != 0]

  assign(paste(dbs[i], "tax", sep = "_"), tax)
  rm(tax)
}



# join counts with classes and families

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

  assign(paste(dbs[i], "counts", sep = "_"), counts)
  assign(paste(dbs[i], "counts_tax", sep = "_"), counts_tax)
}


# filtering out genes with less than 10 read counts, which is the standard use by Deseq

cazy_ident <- cazy_counts_tax$Identifier[rowSums(cazy_counts_tax[, 3:14]) > 9]
length(cazy_ident)

ncyc_ident <- ncyc_counts_tax$Identifier[rowSums(ncyc_counts_tax[, 3:14]) > 9]
length(ncyc_ident)


cazy_counts_tax <- na.omit(cazy_counts_tax[cazy_counts_tax$Identifier %in% cazy_ident, ])
ncyc_counts_tax <- na.omit(ncyc_counts_tax[ncyc_counts_tax$Identifier %in% ncyc_ident, ])

cazy_counts <- na.omit(cazy_counts[cazy_counts$Identifier %in% cazy_ident, ]) # can decide whether having or not the filtering on the counts for the clustering analyses
ncyc_counts <- na.omit(ncyc_counts[ncyc_counts$Identifier %in% ncyc_ident, ])



# nornmalize counts to TPM

# load genes_length for the normalization

genes_length <- read.csv2("F:/metagenomics_analyses_fgcz/files_user_interface/statistical_analyses/results_own/genes_length.txt", header = TRUE, sep = "\t")

cazy_gene_length <- genes_length[, which(names(genes_length) %in% c("Geneid", "Length"))][genes_length$Geneid %in% cazy_counts_tax$Identifier, ]
ncyc_gene_length <- genes_length[, which(names(genes_length) %in% c("Geneid", "Length"))][genes_length$Geneid %in% ncyc_counts_tax$Identifier, ]

# create TPM matrices

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




##### DESEQ ANALYSES ########

names_deseq1 <- c("N10", "N160")
names_deseq2 <- c("S160", "N160")
names_deseq3 <- c("S10", "N10")
names_deseq4 <- c("S160", "S10")



for (i in 1:length(dbs)) {
  print(dbs[i])
  get(paste(dbs[i], "counts_tax", sep = "_")) -> counts

  for (k in 1:4) {
    get(paste("names_deseq", k, sep = "")) -> names_deseq
    print(names_deseq)

    sample_meta <- metadata[(metadata$names == names_deseq[1] | metadata$names == names_deseq[2]), ]
    sample_meta <- sample_meta[order(factor(sample_meta$names, levels = names_deseq)), ]
    sample_meta$group <- c(rep("a", 3), rep("b", 3))
    sample_meta$names2 <- as.factor(paste(sample_meta$group, sample_meta$names, sep = "_"))
    sample_meta <- droplevels(sample_meta)

    counts_subset <- counts[, c("Identifier", paste(sample_meta$id))]
    rownames(counts_subset) <- counts[, c("Identifier")]
    counts_subset <- counts_subset[, -1]

    counts_subset <- counts_subset[, colnames(counts_subset) %in% rownames(sample_meta)] # make sure again
    counts_subset <- counts_subset[, order(factor(names(counts_subset), levels = rownames(sample_meta)))]

    print(all(rownames(sample_meta) %in% colnames(counts_subset)))
    print(ncol(counts_subset) == nrow(sample_meta))
    print(is.data.frame(sample_meta))
    print(is.data.frame(counts_subset))


    dds <- DESeqDataSetFromMatrix(
      countData = counts_subset,
      colData = sample_meta,
      design = ~names2
    )

    dds_test <- DESeq(dds)
    dds_results <- results(dds_test)
    print(resultsNames(dds_test))

    print(paste(dbs[i], k))
    print(summary(dds_results))


    dds_results2 <- results(dds_test, name = resultsNames(dds_test)[2])
    dds_results2 <- as.data.frame(dds_results2)
    dds_results2 <- na.omit(dds_results2[dds_results2$padj < 0.05, ])
    dds_results2$Identifier <- rownames(dds_results2)

    dds_results2_tax <- merge(dds_results2, counts, by = "Identifier", all = FALSE)


    assign(paste(dbs[i], resultsNames(dds_test)[2], "deseq_results", sep = "_"), dds_results)
    rm(dds_results)

    assign(paste(dbs[i], resultsNames(dds_test)[2], "deseq_tax", sep = "_"), dds_results2_tax)
    rm(dds_results2_tax)
  }
}



######## DESEQ REPRESENT #########


# arrange tables
cazy_counts_tax$class2 <- as.factor(gsub(" ", "", cazy_counts_tax$class))
cazy_counts_tax$class2 <- as.factor(paste(cazy_counts_tax$class, cazy_counts_tax$order))

cazy_tpm_tax$class2 <- as.factor(gsub(" ", "", cazy_tpm_tax$class))
cazy_tpm_tax$class2 <- as.factor(paste(cazy_tpm_tax$class, cazy_tpm_tax$order))

ncyc_counts_tax$class2 <- ncyc_counts_tax$class
ncyc_tpm_tax$class2 <- ncyc_tpm_tax$class

cazy_tpm_tax$total <- rowSums(cazy_tpm_tax[, which(names(cazy_tpm_tax) %in% as.character(metadata$id))]) # absolute abundance
ncyc_tpm_tax$total <- rowSums(ncyc_tpm_tax[, which(names(ncyc_tpm_tax) %in% as.character(metadata$id))]) # absolute abundance

cazy_tpm_tax$rel <- rowSums(cazy_tpm_tax[, which(names(cazy_tpm_tax) %in% as.character(metadata$id))]) / sum(cazy_tpm_tax[, which(names(tpm_tax) %in% as.character(metadata$id))]) # relative abundance
ncyc_tpm_tax$rel <- rowSums(ncyc_tpm_tax[, which(names(ncyc_tpm_tax) %in% as.character(metadata$id))]) / sum(ncyc_tpm_tax[, which(names(tpm_tax) %in% as.character(metadata$id))]) # relative abundance


# load complementary info

cazy_info <- read.xlsx("cazy_enzyme_descriptions_simplified2.xlsx", sheetName = "arranged2", header = TRUE, encoding = "UTF-8", stringsAsFactors = F)
cazy_info <- cazy_info[cazy_info$involved_in == "degradation", ]
cazy_info <- cazy_info[cazy_info$compound_associated != "nondescribed" | cazy_info$compound_associated != "varied", ] # delete the non-describe activities or the activities that groups that englobe too many activities

cazy_info$class3 <- gsub(" ", "", cazy_info$class)
droplevels(cazy_info)


cazy_lvls <- c("lignin", "pectin", "cellulose", "hemicellulose", "chitin", "starch") # excluding the bacterial and algi materials


ncyc_info <- ncyc_tax[, which(names(ncyc_tax) %in% c("class", "phylum"))]
ncyc_info$class3 <- gsub(" ", "", ncyc_info$class)
ncyc_info <- ncyc_info[!duplicated(ncyc_info), ]

write.xlsx(ncyc_info, "ncyc_info.xls")

ncyc_info <- read.xlsx("ncyc_info_edited.xls", sheetName = "Sheet1", header = TRUE, encoding = "UTF-8", stringsAsFactors = T)
ncyc_info <- ncyc_info[!duplicated(ncyc_info), ]
ncyc_info$phylum <- as.factor(trimws(ncyc_info$phylum, "b"))


ncyc_info$phylum[ncyc_info$phylum == "Denitrification"] <- "Denitrification, Dissimilatory nitrate reduction"
ncyc_info$phylum[ncyc_info$phylum == "Dissimilatory nitrate reduction"] <- "Denitrification, Dissimilatory nitrate reduction"


ncyc_info$class4 <- as.character(ncyc_info$class4)
ncyc_info$class4[ncyc_info$phylum == "Denitrification, Dissimilatory nitrate reduction" & ncyc_info$class4 == "nir"] <- "nir1"
ncyc_info$class4[ncyc_info$phylum == "Assimilatory nitrate reduction" & ncyc_info$class4 == "nir"] <- "nir2"
ncyc_info$class4[ncyc_info$phylum == "Denitrification, Dissimilatory nitrate reduction" & ncyc_info$class4 == "nar"] <- "nar1"
ncyc_info$class4[ncyc_info$phylum == "Assimilatory nitrate reduction" & ncyc_info$class4 == "nar"] <- "nar2"

ncyc_info <- droplevels(ncyc_info)

ncyc_lvls <- c("Nitrogen fixation", "Nitrification", "Denitrification, Dissimilatory nitrate reduction", "Anammox", "Assimilatory nitrate reduction", "Organic degradation and synthesis", "Others")


# for representation fig 3 and 4

text_size_1 <- 5
text_size <- 10

theme_rep <- theme(
  panel.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background = element_rect(fill = "white"),
  axis.text.y = element_text(size = text_size, color = "black"),
  axis.text.x = element_text(size = text_size, color = "black"),
  axis.title = element_text(size = text_size),
  legend.text = element_text(size = text_size),
  legend.title = element_text(size = text_size),
  legend.position = "top",
  legend.justification = "top",
  legend.direction = "horizontal",
  legend.key.size = unit(1, "cm"),
  legend.key.width = unit(1, "cm"),
  legend.key.height = unit(0.5, "cm"),
  axis.line = element_line(colour = "black", size = 1, linetype = 1),
  axis.ticks = element_line(colour = "black", size = 1),
  axis.ticks.length = unit(0.15, "cm"),
  axis.text.y.right = element_text(color = "red"),
  axis.title.y.right = element_text(color = "red"),
  plot.margin = grid::unit(c(0, 2, 0, 0), "cm")
)


theme_rep_1 <- theme(
  panel.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background = element_rect(fill = "white"),
  axis.text.y = element_text(size = text_size_1, color = "black"),
  axis.text.x = element_text(size = text_size, color = "black"),
  axis.title = element_text(size = text_size),
  legend.text = element_text(size = text_size),
  legend.title = element_text(size = text_size),
  legend.position = "top",
  legend.justification = "top",
  legend.direction = "horizontal",
  legend.key.size = unit(1, "cm"),
  legend.key.width = unit(1, "cm"),
  legend.key.height = unit(0.5, "cm"),
  axis.line = element_line(colour = "black", size = 1, linetype = 1),
  axis.ticks = element_line(colour = "black", size = 1),
  axis.ticks.length = unit(0.15, "cm"),
  axis.text.y.right = element_text(color = "red"),
  axis.title.y.right = element_text(color = "red"),
  plot.margin = grid::unit(c(0, 2, 0, 0), "cm")
)


color <- c("up" = "orange3", "down" = "slateblue4")


names_tables_ds <- c("names2_b_N160_vs_a_N10", "names2_b_N160_vs_a_S160", "names2_b_N10_vs_a_S10", "names2_b_S10_vs_a_S160") # add the deseq of south for a control
pair <- c("N160 vs N10", "N160 vs S160", "N10 vs S10", "S10 vs S160")


# selection of genes to and representation of CAZYs and Ncyc groups fig 3 and 4


for (i in 1:length(dbs)) {
  for (k in 1:length(names_tables_ds)) {
    print(dbs[i])
    print(names_tables_ds[k])

    get(paste(dbs[i], names_tables_ds[k], "deseq_tax", sep = "_")) -> dds_tax
    get(paste(dbs[i], "tpm_tax", sep = "_")) -> tpm
    get(paste(dbs[i], "info", sep = "_")) -> table
    get(paste(dbs[i], "lvls", sep = "_")) -> lvls

    tpm <- tpm[, which(names(tpm) %in% c("FunID", "Identifier", "rel", "total", as.character(metadata$id)))]

    table$class <- trimws(table$class, "b")
    table$class3 <- trimws(table$class3, "b")

    dds_tax$phylum <- trimws(dds_tax$phylum, "b")
    dds_tax$class <- trimws(dds_tax$class, "b")
    dds_tax$class2 <- trimws(dds_tax$class, "b")

    dds_tax <- dds_tax[dds_tax$phylum != "GlycosylTransferases", ]
    dds_tax <- dds_tax[dds_tax$phylum != "Carbohydrate-Binding Modules", ]

    dds_tax <- dds_tax[, which(names(dds_tax) %in% c("Identifier", "log2FoldChange", "padj", "phylum", "class", "order", "class2"))]
    dds_tax2 <- merge(dds_tax, tpm, by = "Identifier", all = FALSE)

    dds_tax2 <- dds_tax2[dds_tax2$total > quantile(dds_tax2$total, 0.20), ] # select based on abundance for now of 80%
    dds_tax2$class3 <- gsub(" ", "", dds_tax2$class2)

    dds_tax2 <- na.omit(dds_tax2[, -which(names(dds_tax2) %in% c("Identifier", "FunID", "kingdom", "phylum", "order", "class", "length"))])

    counts2 <- dds_tax2[, which(names(dds_tax2) %in% c("class3", "rel"))]
    counts2 <- aggregate(. ~ class3, counts2, "sum")


    # only visualizing data with all values overrepresented or underrepresented

    categ <- unique(dds_tax2$class3)
    categ_select1 <- c()
    categ_select2 <- c()

    for (l in 1:length(categ)) {
      if (all(dds_tax2$log2FoldChange[dds_tax2$class3 == categ[l]] > 0)) {
        # categ_select[[m]]<-categ[l]
        categ_select1 <- c(categ_select1, as.character(categ[l]))
        # print(categ[l])
      }
    }

    for (l in 1:length(categ)) {
      if (all(dds_tax2$log2FoldChange[dds_tax2$class3 == categ[l]] < 0)) {
        # categ_select[[m]]<-categ[l]
        categ_select2 <- c(categ_select2, as.character(categ[l]))
        # print(categ[l])
      }
    }


    dds_taxa_subset_up <- dds_tax2[dds_tax2$class3 %in% categ_select1, ]
    dds_taxa_subset_up$trend <- rep("up", nrow(dds_taxa_subset_up))

    dds_taxa_subset_down <- dds_tax2[dds_tax2$class3 %in% categ_select2, ]
    dds_taxa_subset_down$trend <- rep("down", nrow(dds_taxa_subset_down))

    dds_taxa_all <- rbind(dds_taxa_subset_up, dds_taxa_subset_down)
    dds_taxa_all <- merge(dds_taxa_all, table, by = "class3", all = F)

    print(setdiff(dds_taxa_all$class3, table$class3))

    dds_taxa_all$class_desc <- as.factor(paste(dds_taxa_all$activity_main, "(", dds_taxa_all$class2, ")", sep = ""))


    if (i == 1) {
      dds_taxa_all <- dds_taxa_all[(dds_taxa_all$compound_associated %in% lvls), ]
      dds_taxa_all$compound_associated <- factor(dds_taxa_all$compound_associated, levels = lvls)
      # dds_taxa_all=dds_taxa_all[order(factor(dds_taxa_all$compound_associated, levels = lvls), -dds_taxa_all$log2FoldChange),]
      dds_taxa_all <- dds_taxa_all[order(-dds_taxa_all$log2FoldChange), ]
    } else if (i == 2) {
      # dds_taxa_all=dds_taxa_all[order(factor(dds_taxa_all$phylum, levels = lvls), -dds_taxa_all$log2FoldChange),]
      dds_taxa_all$compound_associated <- dds_taxa_all$phylum
      dds_taxa_all$class3 <- dds_taxa_all$class4
    }


    dds_taxa_summary <- summarySE(dds_taxa_all, measurevar = "log2FoldChange", groupvars = c("class3", "compound_associated"))
    dds_taxa_summary$trend <- ifelse(dds_taxa_summary$log2FoldChange > 0, "up", "down")
    dds_taxa_summary[is.na(dds_taxa_summary)] <- 0
    dds_taxa_summary$count2 <- counts2$rel[match(dds_taxa_summary$class3, counts2$class3)]
    dds_taxa_summary <- dds_taxa_summary[dds_taxa_summary$compound_associated != "Others", ]
    dds_taxa_all <- dds_taxa_all[dds_taxa_all$compound_associated != "Others", ]

    if (i == 1) {
      dds_taxa_summary <- dds_taxa_summary[order(factor(dds_taxa_summary$compound_associated, levels = lvls), -dds_taxa_summary$log2FoldChange), ]
      order <- dds_taxa_summary$class3
    } else if (i == 2) {
      order <- rev(dds_taxa_summary$class3)
    }

    # represent like bars

    rep <- ggplot(dds_taxa_summary, aes(factor(class3, levels = rev(order)), log2FoldChange)) +
      geom_errorbar(aes(ymin = log2FoldChange - sd, ymax = log2FoldChange + sd),
        colour = "black", stat = "identity", position = position_dodge(width = 0.9), width = 0.5
      ) +
      geom_col(color = "black", aes(fill = trend), show.legend = FALSE) +
      geom_point(data = dds_taxa_all, aes(factor(class3, levels = rev(order)), y = log2FoldChange), pch = 21, color = "azure4", fill = NA) +
      scale_fill_manual(values = color, guide = "none") +
      facet_grid(factor(compound_associated, levels = lvls) ~ ., scales = "free", space = "free") +
      theme(panel.spacing = unit(0.8, "lines")) +
      coord_flip(clip = "on") +
      xlab("") +
      ggtitle(paste(names_tables_ds[k], dbs[i]), "class3") +
      ylab("Log2 fold-change") +
      theme_rep +
      guides(size = guide_legend(reverse = T)) +
      theme(
        legend.position = "top",
        legend.direction = "horizontal",
        axis.text.x = element_text(size = 8, color = "black")
      )

    print(rep)



    # path=paste(path_data,"tables_own_analyses/", sep="")

    # write.xlsx(dds_taxa_summary, file = paste(path, names_tables_ds[k], dbs[i], "class_selected.xlsx", sep="_"),
    # sheetName = "summary", append = FALSE)

    # write.xlsx(dds_tax2, file = paste(path, names_tables_ds[k], dbs[i], "class_all_significant.xlsx", sep="_"),
    # sheetName = "summary", append = FALSE)


    assign(paste(names_tables_ds[k], dbs[i], "dds_taxa_summary", sep = "_"), dds_taxa_summary)
  }
}









#########################   CLUSTERING ANALYSES (fig S2) ###############################


for (i in 1:length(dbs)) {
  print(dbs[i])
  get(paste(dbs[i], "tpm_tax", sep = "_")) -> tpm

  tpm <- tpm[, which(names(tpm) %in% metadata$id)]

  tpm_counts <- t(tpm)
  colSums(tpm)

  dist_mat <- vegdist(tpm_counts, method = "bray")

  write.csv2(tpm_counts, paste(dbs[i], "_counts_prime2.csv", sep = ""))


  ## PCOA

  pcoa <- cmdscale(dist_mat, eig = TRUE) # pcoa_ with vegan ## I sqrt because is recommended/ usual in the literature to equilibrate variance/ equilize weights between abundant and rare taxa

  barplot(pcoa$eig)

  pcoa_table <- as.data.frame(pcoa$points[, 1:2])

  Eigenvalues <- eigenvals(pcoa)
  Variance <- Eigenvalues / sum(Eigenvalues)
  Variance1 <- 100 * signif(Variance[1], 2)
  Variance2 <- 100 * signif(Variance[2], 2)

  pcoa_table$name <- metadata$names[match(rownames(pcoa_table), metadata$id)]
  pcoa_table$sample <- rownames(pcoa_table)


  pcoa_rep <- ggplot(pcoa_table, aes(V1, V2)) +
    geom_point(
      size = 5, aes(fill = name),
      color = "black", pch = 21
    ) +
    scale_fill_manual(values = soils_color, guide = "legend") +
    # scale_color_manual(values =soils_color) +
    xlab(paste("PC1 (", Variance1, "% )")) +
    ylab(paste("PC2 (", Variance2, "% )")) +
    ggtitle(paste(dbs[i], "annotated")) +
    theme_rep #+
  # guides(fill=guide_legend(title="samples", ncol = 4))
  # guides(fill=FALSE)


  print(pcoa_rep)
}









#########################   DIVERSITY AND ABUNDANCE, fig S4 and S5 ###############################

# for representation
soils_color <- c("N160" = "steelblue4", "N10" = "darkorchid1", "S160" = "springgreen4", "S10" = "coral4") # to match the colors of the previous paper
levels_soils <- c("N10", "N160", "S10", "S160")


for (i in 1:length(dbs)) {
  print(dbs[i])
  get(paste(dbs[i], "counts_tax", sep = "_")) -> richness

  if (i == 2) {
    richness$phylum <- ncyc_info$phylum[match(richness$class, as.character(ncyc_info$class3))]
    richness <- separate_rows(richness, phylum, sep = ",")
    richness$phylum <- trimws(richness$phylum, "b")
    richness <- richness[richness$class2 != "hcp", ]
  }


  richness <- richness[, -which(names(richness) %in% c("FunID", "Identifier", "kingdom", "class", "order", "class2"))]


  richness[, -which(names(richness) %in% c("phylum"))] <- +(richness[, -which(names(richness) %in% c("phylum"))] > 0)

  richness_agg <- aggregate(. ~ phylum, richness, "sum") # first sum the numbers to acquire the total richness for each sample

  richness_melt <- melt(richness_agg, by = "phylum")
  richness_melt$variable <- as.character(metadata$names)[match(richness_melt$variable, as.character(metadata$id))]

  richness_summary <- summarySE(richness_melt, measurevar = "value", groupvars = c("phylum", "variable"))
  write.xlsx(richness_melt, paste(path_data, "richness_stats/", dbs[i], "_richness.xlsx", sep = ""))
  write.xlsx(richness_summary, paste(path_data, "richness_stats/", dbs[i], "_richness_summary.xlsx", sep = ""))

  assign(paste(dbs[i], "richness", sep = "_"), richness_melt)
}




# abundance and shannon in the tpm counts


for (i in 1:length(dbs)) {
  print(dbs[i])
  get(paste(dbs[i], "tpm_tax", sep = "_")) -> abund

  if (i == 2) {
    abund$phylum <- ncyc_info$phylum[match(abund$class, as.character(ncyc_info$class3))]
    abund <- separate_rows(abund, phylum, sep = ",")
    abund$phylum <- trimws(abund$phylum, "b")
    abund <- abund[abund$class2 != "hcp", ]
  }

  abund <- abund[, -which(names(abund) %in% c("FunID", "Identifier", "kingdom", "class", "order", "length", "class2", "total", "rel"))]

  abund_agg <- aggregate(. ~ phylum, abund, "sum") # first sum the numbers to acquire the total abund for each sample


  abund_melt <- melt(abund_agg, by = "phylum")
  abund_melt$variable <- as.character(metadata$names)[match(abund_melt$variable, as.character(metadata$id))]

  abund_summary <- summarySE(abund_melt, measurevar = "value", groupvars = c("phylum", "variable"))

  write.xlsx(abund_melt, paste(path_data, "abund_stats/", dbs[i], "_abund.xlsx", sep = ""))
  write.xlsx(abund_summary, paste(path_data, "abund_stats/", dbs[i], "_abund_summary.xlsx", sep = ""))


  assign(paste(dbs[i], "abund", sep = "_"), abund_melt)


  #### shannon diversity

  shannon_all <- data.frame()

  for (j in 1:length(variables)) {
    abund_subset2 <- abund[abund$phylum == paste(variables[j]), ]
    abund_subset2 <- abund_subset2[, -ncol(abund_subset2)]
    abund_subset2 <- t(abund_subset2)

    if (variables[j] == "Nitrogen fixation") {
      abund_subset2 <- abund_subset2[abund_subset2 > 0, ]
      abund_subset_div <- cbind.data.frame(variable = names(abund_subset2), value = diversity(abund_subset2, index = "shannon"))
    } else {
      abund_subset_div <- cbind.data.frame(variable = row.names(abund_subset2), value = diversity(abund_subset2, index = "shannon"))
    }

    abund_subset_div$variable <- as.character(metadata$names)[match(abund_subset_div$variable, as.character(metadata$id))]

    abund_div_summary <- summarySE(abund_subset_div, measurevar = "value", groupvars = c("variable"))
    abund_div_summary$phylum <- rep(paste(variables[j]), nrow(abund_div_summary))

    shannon_all <- rbind(shannon_all, abund_div_summary)
    row1 <- c("N160", 3, 0, NA, NA, NA, "Nitrogen fixation")
    row2 <- c("S160", 3, 0, NA, NA, NA, "Nitrogen fixation")

    shannon_all <- rbind(shannon_all, row1, row2)
    shannon_all$value <- as.numeric(shannon_all$value)
    shannon_all$sd <- as.numeric(shannon_all$sd)

    # represent shannon

    rep <- ggplot(shannon_all, aes(factor(phylum), value, group = factor(variable, levels = levels_soils), fill = variable)) +
      geom_errorbar(aes(ymin = value - sd, ymax = value + sd), stat = "identity", position = position_dodge(0.9), colour = "black", width = 0.5) +
      geom_col(position = position_dodge(0.9), color = "black") +
      scale_fill_manual(values = soils_color) +
      scale_y_continuous(labels = function(x) format(x, big.mark = " ", decimal.mark = ".", scientific = FALSE)) +
      facet_wrap(~phylum, scales = "free", ncol = length(variables), drop = FALSE) +
      theme(panel.spacing = unit(2, "lines")) +
      xlab("") +
      ylab("") +
      theme_rep +
      theme(
        legend.position = "none",
        legend.direction = "horizontal",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )

    print(rep)
  }

  assign(paste(dbs[i], "shannon", sep = "_"), shannon_all)
}
