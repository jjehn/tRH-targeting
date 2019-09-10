# load required packages
library("ggplot2")
library("caTools")
library("gridExtra")


# set file names
#file_UTR3_total <- c("DESeq2result_5ptRh-Glu-CTC_3pUTR.strict.txt")
#file_UTR3_EE <- c("DESeq2result_5ptRh-Glu-CTC_3pUTR.EE.strict.txt")
#file_UTR3_EJ <- c("DESeq2result_5ptRh-Glu-CTC_3pUTR.EJ.strict.txt")
#file_UTR3_ES <- c("DESeq2result_5ptRh-Glu-CTC_3pUTR.ES.strict.txt")

#file_UTR5_total <- c("DESeq2result_5ptRh-Glu-CTC_5pUTR.strict.txt")
#file_UTR5_EE <- c("DESeq2result_5ptRh-Glu-CTC_5pUTR.EE.strict.txt")
#file_UTR5_EJ <- c("DESeq2result_5ptRh-Glu-CTC_5pUTR.EJ.strict.txt")
#file_UTR5_ES <- c("DESeq2result_5ptRh-Glu-CTC_5pUTR.ES.strict.txt")

#file_CDS_total <- c("DESeq2result_5ptRh-Glu-CTC_CDS.strict.txt")
#file_CDS_EE <- c("DESeq2result_5ptRh-Glu-CTC_CDS.EE.strict.txt")
#file_CDS_EJ <- c("DESeq2result_5ptRh-Glu-CTC_CDS.EJ.strict.txt")
#file_CDS_ES <- c("DESeq2result_5ptRh-Glu-CTC_CDS.ES.strict.txt")


file_UTR3_total <- c("DESeq2result_5ptRh-Gly-GCC_3pUTR.strict.txt")
file_UTR3_EE <- c("DESeq2result_5ptRh-Gly-GCC_3pUTR.EE.strict.txt")
file_UTR3_EJ <- c("DESeq2result_5ptRh-Gly-GCC_3pUTR.EJ.strict.txt")
file_UTR3_ES <- c("DESeq2result_5ptRh-Gly-GCC_3pUTR.ES.strict.txt")

file_UTR5_total <- c("DESeq2result_5ptRh-Gly-GCC_5pUTR.strict.txt")
file_UTR5_EE <- c("DESeq2result_5ptRh-Gly-GCC_5pUTR.EE.strict.txt")
file_UTR5_EJ <- c("DESeq2result_5ptRh-Gly-GCC_5pUTR.EJ.strict.txt")
file_UTR5_ES <- c("DESeq2result_5ptRh-Gly-GCC_5pUTR.ES.strict.txt")

file_CDS_total <- c("DESeq2result_5ptRh-Gly-GCC_CDS.strict.txt")
file_CDS_EE <- c("DESeq2result_5ptRh-Gly-GCC_CDS.EE.strict.txt")
file_CDS_EJ <- c("DESeq2result_5ptRh-Gly-GCC_CDS.EJ.strict.txt")
file_CDS_ES <- c("DESeq2result_5ptRh-Gly-GCC_CDS.ES.strict.txt")


# function to get info of 5-mer mapping

gather_plot_data <- function(file_region) {
  region <- read.delim(file_region, header = T)
  skip_miR_info <- nrow(region)-1
  region <- read.delim(file_region, header = T, nrows = skip_miR_info)
  
  # subset for 5-mers (in case all kinds of kmers were tested)
  region <- region[region$length==5,]
  
  # introduce level for correct plotting order
  region$level <- row.names(region)
  region$level <- factor(region$level, levels = region$level)
  return(region)
  
}


#######################################################################################################
# 3' UTR
#######################################################################################################

# 3' UTR total
UTR3_total <- gather_plot_data(file_UTR3_total)

# 3' UTR EJ/ES/ESS
UTR3_EJ <- gather_plot_data(file_UTR3_EJ)
UTR3_ES <- gather_plot_data(file_UTR3_ES)
UTR3_EE <- gather_plot_data(file_UTR3_EE)

# combine 3' UTR info for 5-mers (total/EJ/ES/ESS)
UTR3_total$group <- c("total")
UTR3_EJ$group <- c("EJ")
UTR3_ES$group <- c("ES")
UTR3_EE$group <- c("EE")

UTR3 <- rbind(UTR3_total, UTR3_EJ, UTR3_EE, UTR3_ES)

# make length a factor for grouping by group (total/EJ/ES/ESS)
UTR3$group <- factor(UTR3$group, levels = c("total", "EJ", "EE", "ES"))

# plot down/total 
plot_UTR3_down <- ggplot(UTR3) + 
  geom_point(aes(x = level, y = runmean(x=UTR3$down.total, k=5)), colour = "blue") +
  geom_point(aes(x = level, y = runmean(x=UTR3$down.total.1, k=5)), colour = "red") + 
  geom_point(aes(x = level, y = runmean(x=UTR3$fraction.targeted, k=5)), colour = "black") +
  facet_grid(~group, scales = "free") +
  theme(axis.title = element_blank(), axis.text=element_blank(), axis.ticks=element_blank(), text = element_text(size=20)) + 
  scale_y_continuous(limits=c(0,1))
#plot_UTR3_down



#######################################################################################################
# 5' UTR
#######################################################################################################

# 5' UTR total/EJ/ES/ESS
UTR5_total <- gather_plot_data(file_UTR5_total)
UTR5_EJ <- gather_plot_data(file_UTR5_EJ)
UTR5_ES <- gather_plot_data(file_UTR5_ES)
UTR5_EE <- gather_plot_data(file_UTR5_EE)

# combine 5' UTR info for 5-mers (total/EJ/ES/ESS)
UTR5_total$group <- c("total")
UTR5_EJ$group <- c("EJ")
UTR5_ES$group <- c("ES")
UTR5_EE$group <- c("EE")

UTR5 <- rbind(UTR5_total, UTR5_EJ, UTR5_EE, UTR5_ES)

# make length a factor for grouping by group (total/EJ/ES/ESS)
UTR5$group <- factor(UTR5$group, levels = c("total", "EJ", "EE", "ES"))

# plot down/total 
plot_UTR5_down <- ggplot(UTR5) + 
  geom_point(aes(x = level, y = runmean(x=UTR5$down.total, k=5)), colour = "blue") +
  geom_point(aes(x = level, y = runmean(x=UTR5$down.total.1, k=5)), colour = "red") + 
  geom_point(aes(x = level, y = runmean(x=UTR5$fraction.targeted, k=5)), colour = "black") +
  facet_grid(~group, scales = "free") +
  theme(axis.title = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), text = element_text(size=20)) + 
  scale_y_continuous(limits=c(0,1))
#plot_UTR5_down



#######################################################################################################
# CDS
#######################################################################################################

# CDS total/EJ/ES/ESS
CDS_total <- gather_plot_data(file_CDS_total)
CDS_EJ <- gather_plot_data(file_CDS_EJ)
CDS_ES <- gather_plot_data(file_CDS_ES)
CDS_EE <- gather_plot_data(file_CDS_EE)

# combine CDS info for 5-mers (total/EJ/ES/ESS)
CDS_total$group <- c("total")
CDS_EJ$group <- c("EJ")
CDS_ES$group <- c("ES")
CDS_EE$group <- c("EE")

CDS <- rbind(CDS_total, CDS_EJ, CDS_EE, CDS_ES)

# make length a factor for grouping by group (total/EJ/ES/ESS)
CDS$group <- factor(CDS$group, levels = c("total", "EJ", "EE", "ES"))

# plot down/total 
plot_CDS_down <- ggplot(CDS) + 
  geom_point(aes(x = level, y = runmean(x=CDS$down.total, k=5)), colour = "blue") +
  geom_point(aes(x = level, y = runmean(x=CDS$down.total.1, k=5)), colour = "red") + 
  geom_point(aes(x = level, y = runmean(x=CDS$fraction.targeted, k=5)), colour = "black") +
  facet_grid(~group, scales = "free") +
  theme(axis.title = element_blank(), axis.text=element_blank(), axis.ticks=element_blank(), text = element_text(size=20)) + 
  scale_y_continuous(limits=c(0,1))
#plot_CDS_down



#######################################################################################################
# arrange plots in 1 figure to export
#######################################################################################################

grid.arrange(plot_UTR5_down, plot_UTR3_down, plot_CDS_down, nrow = 1)