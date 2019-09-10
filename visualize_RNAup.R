# load required packages
library("ggplot2")
library("gridExtra")

file = "RNAup_in.hsap-3pUTR.genes.fas.tRF-5p-Glu-CTC.merge.out"
#file = "RNAup_in.hsap-5pUTR.genes.fas.tRF-5p-Glu-CTC.merge.out"
#file = "RNAup_in.hsap-CDS.genes.fas.tRF-5p-Glu-CTC.merge.out"

# get the results of the DESeq calculations (downloaded from Galaxy)
region <- read.delim(file, header = T)
head(region)
summary(region)

# histogram with density plot of alignment length
g1 <- ggplot(region, aes(x=alignment_length)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666") +
  xlab("alignment lenght") +
  xlim(1,33) +
  theme(text = element_text(size=20)) 
  

# histogram with density plot of start position within tRF
g2 <- ggplot(region, aes(x=tRF_start)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666") +
  xlab("alignment start within 5' tRH-Glu-CTC") +
  xlim(1,33) +
  theme(text = element_text(size=20)) 

# generate matrix with alignments
alignment_matrix <- matrix(0, nrow = nrow(region), ncol = 33)

for (i in c(1:nrow(region))) {
  
  alignment_matrix[i,c(region[i,]$tRF_start:region[i,]$tRF_end)] <- rep(1, c(region[i,]$alignment_length))
}

# identify positional alignment enrichment
alignment_matrix <- cbind(c(1:33), colSums(alignment_matrix))
alignment_df <- as.data.frame(alignment_matrix)
colnames(alignment_df) <- c("position", "counts")

g3 <- ggplot(alignment_df, aes(x=position, y=counts)) + 
  geom_bar(stat="identity", colour="black", fill="white") +
  xlim(1,33) +
  theme(text = element_text(size=20)) 

# arrange the 3 plots in 1 plot
grid.arrange(g1, g2, g3)
