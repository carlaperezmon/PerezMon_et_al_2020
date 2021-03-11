# to have the messages and labels in english
Sys.setlocale(category = "LC_ALL", locale = "english")

options(digits = 3, scipen = -2) # for scientific notation with 2 decimals

#### Genomic data analyses ####

# load packages needed
ipak <- function(pkg) {
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) {
    install.packages(new.pkg, dependencies = TRUE)
  }
  sapply(pkg, require, character.only = TRUE)
}


packages <- c(
  "stringr", "reshape2", "digest", "ggplot2", "scales", "spaa", "xlsx", "viridis", "DESeq2", "Rmisc",
  "tidyr", "vegan", "onewaytests", "userfriendlyscience", "mclust"
)

ipak(packages)

path_data <- "F:/metagenomics_analyses_fgcz/files_user_interface/statistical_analyses/results_own/eggnog/"
path_fig <- "F:/metagenomics_analyses_fgcz/files_user_interface/statistical_analyses/results_own/eggnog/graphs/"

setwd(path_data)


# for sample identificaton
levels_soils <- c("N160", "N10", "S160", "S10")
metadata <- read.xlsx("metadata.xlsx", sheetName = "Sheet1", header = TRUE)
rownames(metadata) <- metadata[, "id"]
metadata <- metadata[order(factor(metadata$names, levels = levels_soils)), ]

cog_cat <- read.xlsx("cog_categories.xlsx", sheetName = "Sheet1", header = TRUE)


# load tables of genes and associated COGs

tax <- read.csv2("F:/metagenomics_analyses_fgcz/files_user_interface/eggnog/anno/eggNOG.anno.S1.txt", header = TRUE, sep = "\t")

tax_eggnog_raw <- tax

files <- list.files("F:/metagenomics_analyses_fgcz/files_user_interface/eggnog/anno/", pattern = ".TXT", ignore.case = T)
names <- str_split(files, "\\.", simplify = TRUE)[, 3]

check <- read.table("F:/metagenomics_analyses_fgcz/files_user_interface/eggnog/aggr/eggNOG.aggr.S1.txt", header = TRUE, sep = "\t", quote = "", fill = TRUE)


for (i in 1:length(files)) {
  tb <- read.csv2(paste("F:/metagenomics_analyses_fgcz/files_user_interface/eggnog/anno/", files[i], sep = ""), header = TRUE, sep = "\t")

  tax_eggnog_raw$m1 <- tb$matchCounts[match(tax_eggnog_raw$Identifier, tb$Identifier)]
  colnames(tax_eggnog_raw)[ncol(tax_eggnog_raw)] <- names[i]
}


tax_eggnog_raw$sum <- rowSums(tax_eggnog_raw[, 6:ncol(tax_eggnog_raw)]) # to apply the threshold of the FGGZ of min 10 counts per gene
tax_eggnog_raw <- tax_eggnog_raw[tax_eggnog_raw$sum > 9, ]
tax_eggnog_raw <- tax_eggnog_raw[, -which(names(tax_eggnog_raw) %in% c("sum"))]

# create normalized TPM count matrices
genes_length <- read.csv2("F:/metagenomics_analyses_fgcz/files_user_interface/statistical_analyses/results_own/genes_length.txt", header = TRUE, sep = "\t")

tax_eggnog_tpm <- tax_eggnog_raw
tax_eggnog_tpm$length <- genes_length$Length[match(tax_eggnog_tpm$Identifier, genes_length$Geneid)]

tax_eggnog_tpm[, 6:17] <- tax_eggnog_tpm[, 6:17] / tax_eggnog_tpm$length # divide by length of gene

# tpm: I do it with apply function because somehow is the only way in which i get same tpm counts as in count qc
tax_eggnog_tpm[, 6:17] <- apply(tax_eggnog_tpm[, 6:17], 2, function(x) x / sum(x) * 10^6)
tax_eggnog_tpm$sum <- rowSums(tax_eggnog_raw[, 6:17]) # just the sum
tax_eggnog_tpm$rel <- rowSums(tax_eggnog_raw[, 6:17]) / sum(tax_eggnog_raw[, 6:17]) # the relative abundance compared to the whole dataset




################## ORDINATION ANALYSES ##################################

# clustering analyses


counts <- tax_eggnog_tpm[, which(names(tax_eggnog_tpm) %in% metadata$id)]
rownames(counts) <- tax_eggnog_tpm$Identifier
tpm_counts <- t(counts)

eggnog_dist_mat <- vegdist(tpm_counts, method = "bray")


## PCOA

pcoa_eggnog <- cmdscale(eggnog_dist_mat, eig = TRUE) # pcoa_eggnog with vegan ## I sqrt because is recommended/ usual in the literature to equilibrate variance/ equilize weights between abundant and rare taxa

barplot(pcoa_eggnog$eig)

pcoa_eggnog_table <- as.data.frame(pcoa_eggnog$points[, 1:2])

Eigenvalues <- eigenvals(pcoa_eggnog)
Variance <- Eigenvalues / sum(Eigenvalues)
Variance1 <- 100 * signif(Variance[1], 2)
Variance2 <- 100 * signif(Variance[2], 2)

pcoa_eggnog_table$name <- metadata$names[match(rownames(pcoa_eggnog_table), metadata$id)]
pcoa_eggnog_table$sample <- rownames(pcoa_eggnog_table)


# pcoa representation fig S2

pcoa_eggnog_rep <- ggplot(pcoa_eggnog_table, aes(V1, V2)) +
  geom_point(
    size = 5, aes(fill = name),
    color = "black", pch = 21
  ) +
  scale_fill_manual(values = soils_color, guide = "legend") +
  # scale_color_manual(values =soils_color) +
  xlab(paste("PC1 (", Variance1, "% )")) +
  ylab(paste("PC2 (", Variance2, "% )")) +
  ggtitle("Eggnog annotated gene") +
  theme_rep #+
# guides(fill=guide_legend(title="samples", ncol = 4))
# guides(fill=FALSE)


print(pcoa_eggnog_rep)



#####################   DESEQ ANALYSES #########################


tables <- c("tax_eggnog_raw")

names_deseq1 <- c("N10", "N160")
names_deseq2 <- c("S160", "N160")
names_deseq3 <- c("S10", "N10")
names_deseq4 <- c("S160", "S10")


for (k in 1:4) {
  print(k)
  get(paste("names_deseq", k, sep = "")) -> names_deseq
  print(names_deseq)

  table <- tax_eggnog_raw[, -which(names(tax_eggnog_raw) %in% c("FunID", "FunGrp", "FunDesc", "matchCounts"))]
  rownames(table) <- table[, c("Identifier")]
  table <- table[, -1]

  sample_meta <- metadata[(metadata$names == names_deseq[1] | metadata$names == names_deseq[2]), ]
  sample_meta <- sample_meta[order(factor(sample_meta$names, levels = names_deseq)), ]
  sample_meta$group <- c(rep("a", 3), rep("b", 3))
  sample_meta$names2 <- as.factor(paste(sample_meta$group, sample_meta$names, sep = "_"))
  sample_meta <- droplevels(sample_meta)

  table_subset <- table[, colnames(table) %in% rownames(sample_meta)]
  table_subset <- table_subset[, order(factor(names(table_subset), levels = rownames(sample_meta)))]

  print(all(rownames(sample_meta) %in% colnames(table_subset)))
  print(ncol(table_subset) == nrow(sample_meta))
  print(is.data.frame(sample_meta))
  print(is.data.frame(table_subset))


  dds <- DESeqDataSetFromMatrix(
    countData = table_subset,
    colData = sample_meta,
    design = ~names2
  )

  dds_test <- DESeq(dds)
  dds_results <- results(dds_test)
  print(resultsNames(dds_test))

  print(paste(table[i], k))
  print(summary(dds_results))


  dds_results2 <- results(dds_test, name = resultsNames(dds_test)[2])
  dds_results2 <- as.data.frame(dds_results2)
  dds_results2 <- na.omit(dds_results2[dds_results2$padj < 0.05, ])
  dds_results2$Identifier <- rownames(dds_results2)

  assign(paste(tables[i], resultsNames(dds_test)[2], "deseq_results", sep = "_"), dds_results)
  rm(dds_results)

  assign(paste(tables[i], resultsNames(dds_test)[2], "deseq_results2", sep = "_"), dds_results2)
  rm(dds_results2)
}







################ REPRESENTATION OF DESEQ AND CHOSEN CATEGORIES fig 2 #####################

rel_group <- tax_eggnog_tpm[, which(names(tax_eggnog_tpm) %in% c("FunGrp", "rel", "sum"))]
rel_group <- aggregate(. ~ FunGrp, rel_group, "sum")
rel_group <- rel_group[order(-rel_group$rel, -rel_group$sum), ]
rel_group$cat <- cog_cat$categories[match(rel_group$FunGrp, cog_cat$FunGrp)]

rel_group <- rel_group[-grep(",", rel_group$FunGrp), ] # remove gene matching to multiple levels
rel_group$FunGrp <- as.factor(gsub(" ", "", rel_group$FunGrp))

rel_group <- rel_group[rel_group$FunGrp != "", ] # remove the genes that did not classified to any category
rel_group <- rel_group[rel_group$FunGrp != "S", ]

rel_group$rel2 <- rel_group$sum / sum(rel_group$sum) # relative abundance compared with selected categories

rel_group_levels <- unique(rel_group$cat)
rel_group_levels2 <- unique(rel_group$FunGrp)

rel_group_levels <- as.character(subset_all_sum_grp$categories)
rel_group_levels2 <- as.character(subset_all_sum_grp$FunGrp) # based on the order of number of gene of higher percentage for N160 vs N10


tables_deseq <- c("tax_eggnog_raw_names2_b_N160_vs_a_N10_deseq_results2", "tax_eggnog_raw_names2_b_N160_vs_a_S160_deseq_results2", "tax_eggnog_raw_names2_b_N10_vs_a_S10_deseq_results2")

# for graphical representation

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
  legend.key.width = unit(2, "cm"),
  legend.key.height = unit(0.5, "cm"),
  axis.line = element_line(colour = "black", size = 1, linetype = 1),
  axis.ticks = element_line(colour = "black", size = 1),
  axis.ticks.length = unit(0.15, "cm"),
  axis.text.y.right = element_text(color = "red"),
  axis.title.y.right = element_text(color = "red"),
  plot.margin = grid::unit(c(0, 2, 0, 0), "cm")
)


fill_color <- function(x) {
  if (x > 0) {
    "orange3"
  } else {
    "slateblue4"
  }
}

color2 <- c("up" = "orange3", "down" = "slateblue4")
color3 <- c("up" = "white", "down" = "white")


for (i in 1:length(tables_deseq)) {
  print(tables_deseq[i])
  get(paste(tables_deseq[i])) -> table_rep

  table_rep_tax <- merge(table_rep, tax_eggnog_tpm[, -which(names(tax_eggnog_tpm) %in% c("matchCounts"))], by = "Identifier", all = FALSE)
  table_rep_tax$FunGrp <- gsub(" ", "", table_rep_tax$FunGrp)

  table_rep_tax <- table_rep_tax[-grep(",", table_rep_tax$FunGrp), ] # remove gene matching to multiple levels
  table_rep_tax$FunGrp <- as.factor(gsub(" ", "", table_rep_tax$FunGrp))

  table_rep_tax <- table_rep_tax[table_rep_tax$FunGrp != "", ] # remove the genes that did not classified to any category
  table_rep_tax <- table_rep_tax[table_rep_tax$FunGrp != "S", ] # remove genes classified to unknown
  table_rep_tax <- table_rep_tax[table_rep_tax$FunGrp != "W", ] # keep only genes that classify to a function
  table_rep_tax <- table_rep_tax[table_rep_tax$FunGrp != "Y", ]
  table_rep_tax <- table_rep_tax[table_rep_tax$FunGrp != "Z", ]


  table_rep_tax <- table_rep_tax[table_rep_tax$FunID %in% tax_filtered$FunID, ]


  # preselection because there are many genes that gives significance

  table_rep_tax$sum <- rowSums(table_rep_tax[, 11:22]) # sum tpm abundance of selected gene
  table_rep_tax$rel2 <- table_rep_tax$sum / sum(table_rep_tax$sum) # relative abundance according to the selected gene. the column rel represent the relative abundance compared to all the annotated genes.
  table_rep_tax$categories <- cog_cat$categories[match(table_rep_tax$FunGrp, cog_cat$FunGrp)]

  table_rep_tax <- table_rep_tax[table_rep_tax$padj < 0.01, ] # select for significance
  table_rep_tax <- table_rep_tax[table_rep_tax$sum > quantile(table_rep_tax$sum, 0.20), ] # select for abundance

  table_rep_tax <- droplevels(table_rep_tax)


  # select for specific genes that only change positively or negatively

  table_rep_tax_subset <- table_rep_tax[, which(names(table_rep_tax) %in% c("log2FoldChange", "FunDesc", "FunGrp", "Identifier", "FunID", "sum", "rel", "categories"))]
  table_rep_tax_subset$FunID2 <- paste(table_rep_tax_subset$FunID, table_rep_tax_subset$FunGrp, sep = "_")

  categ <- as.character(unique(table_rep_tax_subset$FunID2))
  categ_select1 <- c()
  categ_select2 <- c()

  for (l in 1:length(categ)) {
    if (all(table_rep_tax_subset$log2FoldChange[table_rep_tax_subset$FunID2 == categ[l]] > 0)) {
      # categ_select[[m]]<-categ[l]
      categ_select1 <- c(categ_select1, as.character(categ[l]))
      print(categ[l])
    }
  }

  for (l in 1:length(categ)) {
    if (all(table_rep_tax_subset$log2FoldChange[table_rep_tax_subset$FunID2 == categ[l]] < 0)) {
      # categ_select[[m]]<-categ[l]
      categ_select2 <- c(categ_select2, as.character(categ[l]))
      print(categ[l])
    }
  }



  subset_up <- table_rep_tax_subset[table_rep_tax_subset$FunID2 %in% categ_select1, ]
  subset_up$trend <- rep("up", nrow(subset_up))
  subset_up_sum <- summarySE(subset_up, measurevar = "log2FoldChange", groupvars = c("FunGrp", "categories"))

  subset_down <- table_rep_tax_subset[table_rep_tax_subset$FunID2 %in% categ_select2, ]
  subset_down$trend <- rep("down", nrow(subset_down))
  subset_down_sum <- summarySE(subset_down, measurevar = "log2FoldChange", groupvars = c("FunGrp", "categories"))

  subset_all <- rbind(subset_up, subset_down)

  subset_all_agg_rel <- aggregate(rel ~ categories, subset_all, "sum") # relative abundance for the selected groups

  subset_all_sum_grp <- summarySE(subset_all, measurevar = "log2FoldChange", groupvars = c("FunGrp", "categories"))
  subset_all_sum_fun <- summarySE(subset_all, measurevar = "log2FoldChange", groupvars = c("FunID2", "categories"))

  subset_all_sum_fun_grp <- summarySE(subset_all, measurevar = "log2FoldChange", groupvars = c("FunGrp", "categories", "FunID2"))



  # represent grp
  {
    subset_all_sum_grp <- subset_all_sum_grp[order(-subset_all_sum_grp$log2FoldChange), ]

    subset_all_sum_grp$rel <- subset_all_agg_rel$rel[match(subset_all_sum_grp$categories, subset_all_agg_rel$categories)]
    subset_all_sum_grp$up_number <- subset_up_sum$N[match(subset_all_sum_grp$categories, subset_up_sum$categories)]
    subset_all_sum_grp$down_number <- subset_down_sum$N[match(subset_all_sum_grp$categories, subset_down_sum$categories)]
    subset_all_sum_grp$percentage_up <- (subset_all_sum_grp$up_number / subset_all_sum_grp$N) * 100 # percentage of groups going up

    subset_all_sum_grp[is.na(subset_all_sum_grp)] <- 0

    write.xlsx(rel_group_levels, "order_categories.xls")


    # chart for figure

    subset_all_agg_per <- aggregate(rel ~ categories + trend, subset_all, "sum")
    subset_all_agg_per$letter <- cog_cat$FunGrp[match(subset_all_agg_per$categories, cog_cat$categories)]
    subset_all_agg_per$N <- subset_all_sum_grp$N[match(subset_all_agg_per$categories, subset_all_sum_grp$categories)]
    subset_all_agg_per$labels <- paste(subset_all_agg_per$letter, "\n(", subset_all_agg_per$N, ")")

    # order_labels=subset_all_agg_per[order(factor(subset_all_agg_per$letter,levels=rel_group_levels2)),]
    order_labels <- subset_all_agg_per[order(subset_all_agg_per$letter), ]
    order_labels <- unique(subset_all_agg_per$labels)

    y.breaks <- cumsum(subset_all_agg_per$rel) - subset_all_agg_per$rel / 2

    rep3 <- ggplot(subset_all_agg_per, aes(factor(labels, levels = order_labels), rel, fill = trend)) +
      geom_bar(stat = "identity", color = "black", size = 0.01, position = position_fill()) +
      coord_polar(theta = "x", start = 0, direction = -1, clip = "on") +
      scale_fill_manual(values = color2) +
      theme_minimal() +
      ggtitle(paste(tables_deseq[i])) +
      theme(legend.position = "none") +
      xlab("") +
      ylab("") +
      theme(legend.position = "none") +
      theme(strip.background = element_blank(), strip.text = element_blank()) +
      theme(
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_blank()
      )


    print(rep3)
  }

  assign(paste(tables_deseq[i], "subset_all", sep = "_"), subset_all)
  assign(paste(tables_deseq[i], "subset_sum_grp", sep = "_"), subset_all_sum_grp)
  assign(paste(tables_deseq[i], "subset_sum_fun_grp", sep = "_"), subset_all_sum_fun_grp)

  # to create tables
  write.xlsx(subset_all_sum_fun_grp, paste(tables_deseq[i], "_subset_eggnog_grp_fun.xls", sep = ""))
  write.xlsx(subset_all, paste(tables_deseq[i], "_subset_eggnog_all.xls", sep = ""))

  subset_all2 <- subset_all
  subset_all2$FunID <- gsub("(.{7})", "\\1;", as.character(subset_all2$FunID))
  subset_all2 <- separate_rows(subset_all2, FunID, sep = ";")
  subset_all2 <- subset_all2[subset_all2$FunID != "", ]

  subset_all2$c <- tax_md$c[match(subset_all2$FunID, tax_md$FunID)]
  subset_all2$o <- tax_md$o[match(subset_all2$FunID, tax_md$FunID)]
}



###### combine EGGNOG with md5nr db #######

tax_md <- read.csv2("F:/metagenomics_analyses_fgcz/files_user_interface/statistical_analyses/results_own/md5nr/Md5nr_schafberg/Md5nr.aggr.anno.S1.txt", header = TRUE, sep = "\t")
tax_md <- separate(tax_md, "taxonomy", c("k", "p", "c", "o", "f", "g", "s"), sep = ";", extra = "drop", fill = "right") # select only for pathways
names(tax_md)[1] <- "FunID2"

tax_md$FunID2 <- gsub("_.*", "", tax_md$FunID2)

tables_subset <- c("tax_eggnog_raw_names2_b_N160_vs_a_N10_deseq_results2_subset_sum_fun_grp", "tax_eggnog_raw_names2_b_N160_vs_a_S160_deseq_results2_subset_sum_fun_grp", "tax_eggnog_raw_names2_b_N10_vs_a_S10_deseq_results2_subset_sum_fun_grp")
save2 <- c("N160_vs_N10_md", "N160_vs_S160_md", "N10_vs_S10_md")
soils_pair <- c("N160_N10", "N160_S160", "N10_S10")

# selection of groups
groups_selected <- c(
  "Carbohydrate transport and metabolism", "Defense mechanisms",
  "Energy production and conversion", "Inorganic ion transport and metabolism",
  "Posttranslational modification, protein turnover, chaperones", "Replication, recombination and repair",
  "Translation, ribosomal structure and biogenesis", "Cell motility"
)


# for the deseq groups

fill_trend <- function(x) {
  if (x > 0) {
    "up"
  } else {
    "down"
  }
}

for (j in 1:length(tables_deseq)) {
  print(tables_subset[j])
  print(save2[j])
  print(soils_pair[j])

  get(paste(tables_subset[j])) -> table
  # get(paste(tables_subset2[j])) -> table_desc

  table$FunID2 <- gsub("_.*", "", table$FunID2)

  table$c <- tax_md$c[match(table$FunID2, tax_md$FunID2)]
  table$o <- tax_md$o[match(table$FunID2, tax_md$FunID2)]

  # table_desc$o=tax_md$o[match(table_desc$FunID,tax_md$FunID2)]

  md_eggnog <- table[!is.na(table$c), ]

  # md_eggnog$o=gsub(",.*","",md_eggnog$o)


  md_eggnog$c <- gsub("c__", "", md_eggnog$c)
  md_eggnog$c <- trimws(md_eggnog$c, "b")

  md_eggnog$o <- gsub("o__", "", md_eggnog$o)
  md_eggnog$o <- trimws(md_eggnog$o, "b")

  md_eggnog <- md_eggnog[(md_eggnog$c %in% groups_selected), ] # stick to defined groups


  md_eggnog <- md_eggnog[!grepl("^Predicted", md_eggnog$o), ]
  md_eggnog <- md_eggnog[!grepl("^Uncharacterized", md_eggnog$o), ]
  md_eggnog <- md_eggnog[!grepl("^Putative", md_eggnog$o), ]
  md_eggnog <- md_eggnog[!grepl("^Truncated", md_eggnog$o), ]


  md_eggnog[grep("[0-9][\\.,][0-9]", md_eggnog$o), ]$o <- sub(",", ";", md_eggnog[grep("[0-9][\\.,][0-9]", md_eggnog$o), ]$o) # sub is to replace only first match, gsub is to replace all matches

  md_eggnog$o <- gsub("Ni,Fe", "Ni;Fe", md_eggnog$o) # sub is to replace only first match, gsub is to replace all matches
  md_eggnog$o <- gsub(",.*", "", md_eggnog$o) # sub is to replace only first match, gsub is to replace all matches
  md_eggnog$o <- gsub("systems", "system", md_eggnog$o)
  md_eggnog$o <- gsub("\\s*\\([^\\)]+\\)", "", md_eggnog$o)
  md_eggnog$o <- gsub("I", "", md_eggnog$o)
  md_eggnog$o <- gsub("large subunit", "", md_eggnog$o)
  md_eggnog$o <- gsub("small subunit", "", md_eggnog$o)
  md_eggnog$o <- gsub("alpha subunit", "", md_eggnog$o)
  md_eggnog$o <- gsub("beta subunit", "", md_eggnog$o)
  md_eggnog$o <- trimws(md_eggnog$o, "b")


  md_eggnog_sum <- summarySE(md_eggnog, measurevar = "log2FoldChange", groupvars = c("c", "o"))

  md_eggnog_sum$trend <- unlist(lapply(md_eggnog_sum$log2FoldChange, fill_trend))

  md_eggnog_sum$soils <- rep(soils_pair[j], nrow(md_eggnog_sum))

  assign(paste(save2[j], "eggnog_sum", sep = "_"), md_eggnog_sum)
}



md_eggnog_subset_all <- rbind(p_vs_n_md_eggnog_sum, p_vs_ds_md_eggnog_sum, n_vs_s_md_eggnog_sum)


write.xlsx(md_eggnog_subset_all, file = "md_eggnog_all_soils_summary.xlsx", append = FALSE) # table that I accomodated and edited a bit by hand



#### load table

eggnog_md_sel <- read.xlsx("md_eggnog_all_soils_summary_edited_1506.xlsx", sheetName = "Sheet1", header = TRUE, encoding = "UTF-8", stringsAsFactors = TRUE)

eggnog_md_sel_sum2 <- summarySE(eggnog_md_sel, measurevar = "log2FoldChange", groupvars = c("groups_1", "groups_2", "o", "soils"))
eggnog_md_sel_sum2$trend <- unlist(lapply(eggnog_md_sel_sum2$log2FoldChange, fill_trend))


groups_1_levels <- c(
  "cazys", "fermentation", "carbon_fixation", "TCA cycle", "respiration (aerobic)", "respiration (anaerobic)", "Stress response",
  "Defense", "Competition", "DNA, translation regulation", "protein modification and turnover respiration (aerobic)"
)

groups_2_levels <- c(
  "cazys", "fermentation", "Photosynthesis", "calvin cycle (aerobic)", "Arnon-Buchanan cycle (anaerobic)", "methanogenesis",
  "acetogenesis", "TCA cycle", "electron transport chains", "auxiliary", "anaerobic respiration",
  "anaerobic respiration precursors", "metal transporters",
  "Stress response", "Defense", "Competition", "DNA, translation regulation", "protein modification and turnover respiration (aerobic)"
)


eggnog_md_sel_sum2 <- eggnog_md_sel_sum2[order(
  factor(eggnog_md_sel_sum2$groups_1, levels = groups_1_levels),
  factor(eggnog_md_sel_sum2$groups_2, levels = groups_2_levels),
  factor(eggnog_md_sel_sum2$o)
), ]



order_rep <- unique(eggnog_md_sel_sum2$o)


# representation, fig 2B

ggplot(eggnog_md_sel_sum2, aes(factor(o, levels = rev(order_rep)), log2FoldChange, group = factor(soils, levels = soils_pair), fill = trend)) +
  geom_errorbar(aes(ymin = log2FoldChange - sd, ymax = log2FoldChange + sd), stat = "identity", position = position_dodge(0.9), colour = "black", width = 0.5) +
  geom_point(data = eggnog_md_sel, aes(x = factor(o, levels = rev(order_rep)), y = log2FoldChange), pch = 21, color = "azure4", fill = NA) +
  geom_col(position = position_dodge2(width = 0.8, preserve = "single"), color = "black") +
  scale_fill_manual(values = color2) +
  xlab("") +
  coord_flip() +
  facet_grid(cols = vars(factor(soils, levels = soils_pair)), rows = vars(factor(groups_1, levels = groups_1_levels)), scales = "free", space = "free") +
  theme(panel.spacing.x = unit(0.8, "lines"), panel.spacing.y = unit(0.3, "lines")) +
  theme_rep +
  theme(
    legend.position = "none",
    legend.direction = "horizontal"
  )






################### DIVERSITY AND ABUNDANCE, fig S3 #########################



soils_color <- c("N160" = "steelblue4", "N10" = "darkorchid1", "S160" = "springgreen4", "S10" = "coral4") # to match the colors of the previous paper
levels_soils <- c("N10", "N160", "S10", "S160")


#### richnesss

richness <- tax_eggnog_raw[, c(4, 6:ncol(tax_eggnog_raw))]
richness[, -which(names(richness) %in% c("FunGrp"))] <- +(richness[, -which(names(richness) %in% c("FunGrp"))] > 0)


richness_agg <- aggregate(. ~ FunGrp, richness, "sum") # first sum the numbers to acquire the total richness for each sample

richness_agg <- richness_agg[-grep(",", richness_agg$FunGrp), ]

richness_agg$FunGrp <- as.factor(gsub(" ", "", richness_agg$FunGrp))

richness_agg <- aggregate(. ~ FunGrp, richness_agg, "sum")
richness_agg <- richness_agg[richness_agg$FunGrp != "", ] # remove the genes that did not classified to any category
richness_agg <- richness_agg[richness_agg$FunGrp != "S", ]
richness_agg <- richness_agg[richness_agg$FunGrp != "W", ]
richness_agg <- richness_agg[richness_agg$FunGrp != "Y", ]
richness_agg <- richness_agg[richness_agg$FunGrp != "Z", ]


richness_melt <- melt(richness_agg, by = "FunGrp")
richness_melt$variable <- as.character(metadata$names)[match(richness_melt$variable, as.character(metadata$id))]
richness_melt$categories <- as.character(cog_cat$categories)[match(richness_melt$FunGrp, as.character(cog_cat$FunGrp))]


richness_summary <- summarySE(richness_melt, measurevar = "value", groupvars = c("FunGrp", "variable"))
write.xlsx(richness_melt, paste(path_data, "eggnog_rich/eggnog_richness.xlsx", sep = ""))
write.xlsx(richness_summary, paste(path_data, "eggnog_rich/eggnog_richness_summary.xlsx", sep = ""))

# abundance, using tpm normalized read counts

abund <- tax_eggnog_tpm[, c(4, 6:ncol(tax_eggnog_raw))]
abund_agg <- aggregate(. ~ FunGrp, abund, "sum") # first sum the numbers to acquire the total abund for each sample

abund_agg <- abund_agg[-grep(",", abund_agg$FunGrp), ]

abund_agg$FunGrp <- as.factor(gsub(" ", "", abund_agg$FunGrp))

abund_agg <- aggregate(. ~ FunGrp, abund_agg, "sum")
abund_agg <- abund_agg[abund_agg$FunGrp != "", ] # remove the genes that did not classified to any category
abund_agg <- abund_agg[abund_agg$FunGrp != "S", ]
abund_agg <- abund_agg[abund_agg$FunGrp != "W", ]
abund_agg <- abund_agg[abund_agg$FunGrp != "Y", ]
abund_agg <- abund_agg[abund_agg$FunGrp != "Z", ]


abund_melt <- melt(abund_agg, by = "FunGrp")
abund_melt$variable <- as.character(metadata$names)[match(abund_melt$variable, as.character(metadata$id))]
abund_melt$categories <- as.character(cog_cat$categories)[match(abund_melt$FunGrp, as.character(cog_cat$FunGrp))]


abund_summary <- summarySE(abund_melt, measurevar = "value", groupvars = c("FunGrp", "variable"))
write.xlsx(abund_melt, paste(path_data, "eggnog_abund/eggnog_abund.xlsx", sep = ""))
write.xlsx(abund_summary, paste(path_data, "eggnog_abund/eggnog_abund_summary.xlsx", sep = ""))


#### shannon diversity

shannon_all <- data.frame()

for (j in 1:length(variables)) {
  print(variables[j])

  if (variables[j] == "Nitrogen fixation") {
    j <- j + 1
  }


  print(variables[j])

  abund_subset2 <- abund[abund$FunGrp == paste(variables[j]), ]
  abund_subset2 <- abund_subset2[, -which(names(abund_subset2) %in% c("FunGrp"))]
  abund_subset2 <- t(abund_subset2)

  abund_subset_div <- cbind.data.frame(variable = row.names(abund_subset2), value = diversity(abund_subset2, index = "shannon"))

  abund_subset_div$variable <- as.character(metadata$names)[match(abund_subset_div$variable, as.character(metadata$id))]

  abund_div_summary <- summarySE(abund_subset_div, measurevar = "value", groupvars = c("variable"))
  abund_div_summary$FunGrp <- rep(paste(variables[j]), nrow(abund_div_summary))


  shannon_all <- rbind(shannon_all, abund_div_summary)

  assign("shannon_eggnog", shannon_all)
}


# represent shannon, fig S2

rep <- ggplot(shannon_eggnog, aes(factor(FunGrp), value, group = factor(variable, levels = levels_soils), fill = variable)) +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd), stat = "identity", position = position_dodge(0.9), colour = "black", width = 0.5) +
  geom_col(position = position_dodge(0.9), color = "black") +
  scale_fill_manual(values = soils_color) +
  scale_y_continuous(labels = function(x) format(x, big.mark = " ", decimal.mark = ".", scientific = FALSE)) +
  facet_wrap(~FunGrp, scales = "free", nrow = 4) +
  theme(panel.spacing = unit(2, "lines")) +
  xlab("") +
  ylab("Shannon's diversity index") +
  theme_rep +
  theme(
    legend.position = "none",
    legend.direction = "horizontal",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    strip.text.x = element_text(face = "bold")
  )

print(rep)
