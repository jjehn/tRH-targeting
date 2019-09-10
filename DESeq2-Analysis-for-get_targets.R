# load required packages
library("ggplot2")
library("ggrepel")


# get the results of the DESeq calculations (downloaded from Galaxy)
#res <- read.delim("Glu_Galaxy88-[DESeq2_result_file_on_data_75,_data_72,_and_others].tabular", header = F)
res <- read.delim("Gly_Galaxy71-[DESeq2_result_file_on_data_75,_data_72,_and_others].tabular", header = F)
names(res) <- c("GeneID",	"Base_mean",	"log2FoldChange",	"StdErr",	"Wald_Stats",	"P_value",	"P_adj")
head(res)
summary(res)
# remove rows containing NA values
row.has.na <- apply(res, 1, function(x){any(is.na(x))})
sum(row.has.na)
res <- res[!row.has.na,]


# get information of how many genes are differentially expressed
dim(res) # 17979/19522 genes are in the result object
dim(res[res$P_adj<0.01,]) # 6367/7039 genes are significantly differentially expressed
dim(res[res$log2FoldChange>0 & res$P_adj<0.01,]) # 2816/3246 are differentially upregulated in condition 1
dim(res[res$log2FoldChange<0 & res$P_adj<0.01,]) # 3551/3793 are differentially downregulated in condition 1

# make volcano plot
#prepare for colour after significance
res$sign <- ifelse(res$P_adj < 0.01, "significant", "not significant")


g = ggplot(data=res) +
  geom_point(aes(x=res$log2FoldChange, y=-log10(res$P_adj),colour=sign)) +
  scale_colour_manual(values=c("grey","skyblue")) + 
  xlab("log2 fold change") + ylab("-log10 p-value")+
  theme_bw()+
  theme(legend.title=element_blank()) + 
  ggtitle("Volcano plot: 5ptRh Glu/Gly vs NTC")
g



# cutoff genes that are very low expressed -> build count table of downloaded Galaxy featureCounts data
files <- list.files("./featureCountsGalaxy", pattern="_R.fq].tabular$", full.names = T)
counts <- lapply(files, read.delim, row.names=1, head=T)
counts <-do.call(cbind, counts)
head(counts)
tail(counts,10)
gene_length <- read.delim("./featureCountsGalaxy/feature_lengths__N1.tabular", row.names=1, header = T)
head(gene_length)
transcript_bp <- gene_length/1000
head(transcript_bp)
scaling_factor <- apply(counts,2,sum)/1E6
fpkm <- t(t(counts)/scaling_factor)
fpkm <- apply(fpkm,1,mean) # calculate mean of all as only wanted that fpkm is everywhere above 1
fpkm <- fpkm/transcript_bp
colnames(fpkm) <- c("fpkm_mean")

# volcano plot for genes with fpkm > 1 + print res_expressed as input for get_target_for_DESeq2.pl
res_expressed <- merge(res, fpkm, by.x = "GeneID", by.y = "row.names")
res_expressed <- res_expressed[res_expressed$fpkm_mean>1,]
dim(res) # 17979/19522
dim(res_expressed) # 12373/12444
dim(res_expressed[res_expressed$log2FoldChange>0 & res_expressed$P_adj < 0.01,]) # 2640/3016 sign. up regulated
dim(res_expressed[res_expressed$log2FoldChange<0 & res_expressed$P_adj < 0.01,]) # 3286/3426 sign. down regulated
res_expressed$sign <- ifelse(res_expressed$P_adj < 0.01, "significant", "not significant")
g = ggplot(data=res_expressed) +
  geom_point(aes(x=res_expressed$log2FoldChange, y=-log10(res_expressed$P_adj),colour=sign)) +
  scale_colour_manual(values=c("grey","skyblue")) + 
  xlab("log2 fold change") + ylab("-log10 p-value")+
  theme_bw()+
  theme(legend.title=element_blank()) + 
  ggtitle("Volcano plot: 5ptRh Glu/Gly (fpkm>1) vs NTC")
g

#write.table(res_expressed, file = "res_expressed_Glu.txt", row.names = F, col.names = T, quote = F, sep = "\t")
write.table(res_expressed, file = "res_expressed_Gly.txt", row.names = F, col.names = T, quote = F, sep = "\t")



#load libraries
library(ggplot2)
library(extrafontdb)
library(extrafont)
#font_import() just apply again if new fonts are installed on computer
fonts()
fonttable()
loadfonts(device="win")
library(Cairo)

# export plot in a better resolution (Glu)
g = ggplot(data=res_expressed) +
  geom_point(aes(x=res_expressed$log2FoldChange, y=-log10(res_expressed$P_adj),colour=sign)) +
  scale_colour_manual(values=c("grey","skyblue")) + 
  xlab("log2 fold change") + ylab("-log10 p-value")+
  theme_bw()+
  theme(legend.title=element_blank()) + 
  ggtitle("Volcano plot: Glu (fpkm>1) vs NTC")
g
Cairo(file="Glu_DE-Analysis_fpkm1_300_dpi.png", 
      type="png",
      units="cm", 
      width=15, 
      height=10, 
      pointsize=10, 
      dpi=300)
g
dev.off()

# export plot in a better resolution (Gly)
g = ggplot(data=res_expressed) +
  geom_point(aes(x=res_expressed$log2FoldChange, y=-log10(res_expressed$P_adj),colour=sign)) +
  scale_colour_manual(values=c("grey","skyblue")) + 
  xlab("log2 fold change") + ylab("-log10 p-value")+
  theme_bw()+
  theme(legend.title=element_blank()) + 
  ggtitle("Volcano plot: Gly (fpkm>1) vs NTC")
g
Cairo(file="Gly_DE-Analysis_fpkm1_300_dpi.png", 
      type="png",
      units="cm", 
      width=15, 
      height=10, 
      pointsize=10, 
      dpi=300)
g
dev.off()
