#### Genomic data analyses ####

#load packages needed
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}


packages <- c("stringr", "reshape2","digest","AnnotationDbi", "ggplot2", 'scales',"GO.db", "spaa", "xlsx", "viridis")

ipak(packages)

path_data="F:/metagenomics_analyses_fgcz/files_user_interface/statistical_analyses/matchcounts/"
setwd(path_data)

#load data of identity of samples
metadata=read.xlsx("metadata.xlsx", sheetName = "Sheet1", header = TRUE)
rownames(metadata)=metadata[,"id"]
metadata=metadata[order(factor(metadata$names, levels=levels_soils)),]


##### load table of read counts of protein-coding genes


file.names <- dir(path_data, pattern =".txt")

all_samples=read.csv2(paste(path_data,file.names[1], sep=""), header=TRUE, sep="\t") 

for(i in 2:length(file.names)){
  a= read.csv(file.names[i], header=TRUE, sep="\t", stringsAsFactors=FALSE)
  a= a[order(a$Identifier), ]
  names(a)=c("Identifier",gsub(".txt","",file.names[i]))
  
  
  all_samples=merge(all_samples,a,by="Identifier")
  
  rm(a)
}

names(all_samples)[2]="S1"

all_samples_raw=all_samples
rownames(all_samples_raw)=all_samples_raw$Identifier

all_samples_raw=all_samples_raw[rowSums(all_samples_raw) > 9,] #filter for genes with at least 10 read counts

#normalize to TPM
genes_length=read.csv2("F:/metagenomics_analyses_fgcz/files_user_interface/statistical_analyses/results_own/genes_length.txt", header=TRUE, sep="\t") 

all_samples_tpm=all_samples_raw
all_samples_tpm$length=genes_length$Length[match(all_samples_tpm$Identifier, genes_length$Geneid)]

all_samples_tpm[,-ncol(all_samples_tpm)]=all_samples_tpm[,-ncol(all_samples_tpm)]/all_samples_tpm$length #divide by length of gene

#tpm: I do it with apply function because somehow is the only way in which i get same tpm counts as in count qc
all_samples_tpm[,-ncol(all_samples_tpm)]=apply(all_samples_tpm[,-ncol(all_samples_tpm)], 2, function (x) x/sum(x)*10^6)

all_samples_tpm=all_samples_tpm[,-ncol(all_samples_tpm)]
all_samples_tpm=all_samples_tpm[,-1]

colSums(all_samples_tpm)

tpm_counts=t(all_samples_tpm)



#### DIVERSITY AND ABUNDANCE OF PROTEIN-CODING GENES AND GRAPHICAL REPRESENTATION ######

# for graphical reppresentations

soils_color=c("N160"="steelblue4","N10"="darkorchid1","S160"="springgreen4","S10"="coral4")
levels_soils=c("N10","N160","S10","S160")

text_size=12
theme_rep= theme(panel.background= element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 strip.background =element_rect(fill="white"),
                 axis.text.y = element_text(size = text_size, color="black"),
                 axis.text.x = element_text(size = text_size, color="black"),
                 axis.title = element_text(size = text_size),
                 legend.text = element_text(size=text_size),
                 legend.title = element_text(size=text_size),
                 legend.position = "right", 
                 legend.justification = "top",
                 legend.direction = "vertical",
                 legend.key.size = unit(1, "cm"),
                 legend.key.width = unit(2,"cm"), 
                 legend.key.height = unit(0.5,"cm"), 
                 axis.line = element_line(colour="black", size= 1, linetype = 1), 
                 axis.ticks = element_line(colour="black", size=1),
                 axis.ticks.length = unit(0.15, "cm"),
                 axis.text.y.right = element_text(color = "red"),
                 axis.title.y.right = element_text(color= "red"),
                 plot.margin=grid::unit(c(0,2,0,0), "cm"))



#### richnesss


richness<- +(all_samples_tpm > 0)

richness=t(richness)
richness=as.data.frame(richness)

richness_total=rowSums(richness)

richness_total2=as.data.frame(richness_total)
richness_total2$variable=as.character(metadata$names)[match(rownames(richness_total2),as.character(metadata$id))]


richness_summary=summarySE(richness_total2, measurevar="richness_total", groupvars=c("variable"))

#representation of richness fig 1 a 

rep=ggplot(richness_summary, aes(factor(variable, levels=levels_soils), richness_total, fill=variable)) +
  geom_errorbar(aes(ymin=richness_total-sd, ymax=richness_total+sd),stat = "identity", position =position_dodge(0.9), colour="black", width=0.5) +
  geom_col(position=position_dodge(0.9), color="black") +
  scale_fill_manual(values=soils_color)+
  xlab("")+
  scale_y_continuous(labels=function(x) format(x, big.mark = " ", decimal.mark = ".", scientific = FALSE))+
  #coord_flip()+
  ggtitle("richness all functional gene")+
  ylab("richness") +
  theme_rep +
  theme(legend.position="top",
        legend.direction='horizontal',
        axis.text.x= element_blank(),
        axis.ticks.x =  element_blank())

print(rep)  


#statistical tests

print(barplot(tapply(richness_total2$richness_total, INDEX=richness_total2$variable, FUN=var), main="var_soil"))

qqnorm(residuals(aov(richness_total2$richness_total ~ richness_total2$variable)), main="qqnorm")
qqline(residuals(aov(richness_total2$richness_total ~ richness_total2$variable)))

welch.test(richness_total ~ variable, data = richness_total2)

gh=posthocTGH(richness_total2$richness_total, richness_total2$variable, method=c("games-howell"),
              conf.level = 0.95, digits=2, p.adjust="none",
              formatPvalue = TRUE) #Games-Howell test that is a recommended post hoc test after welch

gh_table=cbind.data.frame(sample=rownames(gh$output$games.howell), p=gh$output$games.howell$p)

gh_table=gh_table[order(gh_table$p),] 


#adundance and diversity representations

abund=t(all_samples_raw)
abund=as.data.frame(abund)
abund_total=rowSums(abund)

abund_total2=as.data.frame(abund_total)
abund_total2$variable=as.character(metadata$names)[match(rownames(abund_total2),as.character(metadata$id))]

abund_summary=na.omit(abund_summary)
abund_summary=summarySE(abund_total2, measurevar="abund_total", groupvars=c("variable"))

rep=ggplot(abund_summary, aes(factor(variable, levels=levels_soils), abund_total, fill=variable)) +
  geom_errorbar(aes(ymin=abund_total-sd, ymax=abund_total+sd),stat = "identity", position =position_dodge(0.9), colour="black", width=0.5) +
  geom_col(position=position_dodge(0.9), color="black") +
  scale_fill_manual(values=soils_color)+
  xlab("")+
  scale_y_continuous(labels=function(x) format(x, big.mark = " ", decimal.mark = ".", scientific = FALSE))+
  #coord_flip()+
  ggtitle("abundance all functional gene")+
  ylab("abundance") +
  theme_rep +
  theme(legend.position="top",
        legend.direction='horizontal',
        axis.text.x= element_blank(),
        axis.ticks.x =  element_blank())

print(rep)  

shannon=cbind.data.frame(variable=row.names(tpm_counts), value=diversity(tpm_counts, index="shannon"))
shannon$variable=as.character(metadata$names)[match(shannon$variable,as.character(metadata$id))]

shannon_summary=summarySE(shannon, measurevar="value", groupvars=c("variable"))


#statistical test
print(barplot(tapply(shannon$value, INDEX=shannon$variable, FUN=var), main="var_soil abund"))

qqnorm(residuals(aov(shannon$value ~ shannon$variable)), main=paste("qqnorm div"))
qqline(residuals(aov(shannon$value ~ shannon$variable)))

print(summary(aov(shannon$value ~ shannon$variable)))
TukeyHSD(aov(shannon$value ~ shannon$variable), conf.level=0.95)

print(paste(variables[j]))
welch.test(value ~ variable, data = shannon)

gh=posthocTGH(shannon$value, shannon$variable, method=c("games-howell"),
              conf.level = 0.95, digits=2, p.adjust="none",
              formatPvalue = TRUE) #Games-Howell test that is a recommended post hoc test after welch

gh_table=cbind.data.frame(sample=rownames(gh$output$games.howell), p=gh$output$games.howell$p)

gh_table=gh_table[order(gh_table$p),] 
#gh_significant=gh_table[(gh_table$p <= 0.05),]


#represent shannon (fig 1 b)

rep=ggplot(shannon_summary, aes(factor(variable, levels=levels_soils), value, fill=variable)) +
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd),stat = "identity", position =position_dodge(0.9), colour="black", width=0.5) +
  geom_col(position=position_dodge(0.9), color="black") +
  scale_fill_manual(values=soils_color)+
  xlab("")+
  coord_cartesian(ylim=c(12,16))+
  #coord_flip()+
  ggtitle("richness all functional gene")+
  ylab("Shannon's diversity index") +
  theme_rep +
  theme(legend.position="top",
        legend.direction='horizontal',
        axis.text.x= element_blank(),
        axis.ticks.x =  element_blank())

print(rep)  



## PCOA analyses and graphical representation (fig 1, c)


all_dist_mat=vegdist(tpm_counts,method="bray")

pcoa_all=cmdscale(all_dist_mat,eig=TRUE) #pcoa_all with vegan ## I sqrt because is recommended/ usual in the literature to equilibrate variance/ equilize weights between abundant and rare taxa

barplot(pcoa_all$eig)

pcoa_all_table=as.data.frame(pcoa_all$points[, 1:2])

Eigenvalues <- eigenvals(pcoa_all) 
Variance <- Eigenvalues / sum(Eigenvalues) 
Variance1 <- 100 * signif(Variance[1], 2)
Variance2 <- 100 * signif(Variance[2], 2)

pcoa_all_table$name  = metadata$names[match(rownames(pcoa_all_table), metadata$id)]
pcoa_all_table$sample=rownames(pcoa_all_table)


pcoa_all_rep=ggplot(pcoa_all_table, aes(V1,V2))+ 
  geom_point(size=5, aes(fill=name),
             color="black", pch=21
  ) +
  scale_fill_manual(values=soils_color, guide = "legend")+
  #scale_color_manual(values =soils_color) +
  xlab(paste("PC1 (", Variance1, "% )")) + 
  ylab(paste("PC2 (", Variance2, "% )")) +
  ggtitle("all annotated gene") +
  theme_rep #+
#guides(fill=guide_legend(title="samples", ncol = 4))
#guides(fill=FALSE)


print(pcoa_all_rep)


#export for Permanova analyses in PRIMER7

all_dist_mat2 <- as.matrix(all_dist_mat)
all_dist_mat2[upper.tri(all_dist_mat2, diag = FALSE)] <- ""
write.csv2(all_dist_mat2, "all_dist_mat.csv")







  



