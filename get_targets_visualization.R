# load required packages
library("ggplot2")
library("caTools")
library("gridExtra")


# set file names
file_UTR3 <- c("DESeq2result_tRF-5p-Glu-CTC_3pUTR.strict.txt")
file_UTR5 <- c("DESeq2result_tRF-5p-Glu-CTC_5pUTR.strict.txt")
file_CDS <- c("DESeq2result_tRF-5p-Glu-CTC_CDS.strict.txt")

#file_UTR3 <- c("DESeq2result_tRF-5p-Gly-GCC_3pUTR.strict.txt")
#file_UTR5 <- c("DESeq2result_tRF-5p-Gly-GCC_5pUTR.strict.txt")
#file_CDS <- c("DESeq2result_tRF-5p-Gly-GCC_CDS.strict.txt")


#######################################################################################################
# function gathering plot data
#######################################################################################################

gather_plot_data <- function(file_region) {
  # load get_target output for respective transcript region (remove miR info in last row)
  region <- read.delim(file_region, header = T)
  skip_miR_info <- nrow(region)-1
  region <- read.delim(file_region, header = T, nrows = skip_miR_info)
  
  # subset for interesting length of window
  length_max_plot <- region[which.max(region$length),1]-3
  length_max_plot <- ifelse(length_max_plot>21,21,length_max_plot)
  region <- region[region$length<length_max_plot,]
  
  # introduce level for correct plotting order
  region$level <- row.names(region)
  region$level <- factor(region$level, levels = region$level)
  return(region)
  
}



#######################################################################################################
# generating plots
#######################################################################################################

UTR5 <- gather_plot_data(file_UTR5)
CDS <- gather_plot_data(file_CDS)
UTR3 <- gather_plot_data(file_UTR3)


#######################################################################################################
# 5' UTR
#######################################################################################################

# plot down/total
plot_UTR5_down <- ggplot(UTR5) + 
  geom_point(aes(x = level, y = runmean(x=UTR5$down.total, k=5), colour = "targeted")) +
  geom_point(aes(x = level, y = runmean(x=UTR5$down.total.1, k=5), colour = "not targeted")) +
  scale_colour_manual(values=c("red", "blue")) +
  geom_point(aes(x = level, y = runmean(x=UTR5$fraction.targeted, k=5), shape = "fraction targeted"), colour = "black") +
  scale_shape_manual(values=c(16)) +
  facet_grid(~length, scales = "free") +
  labs(y = "downregulated / total", colour = NULL, shape = NULL) +
  theme(axis.title=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), text = element_text(size=20)) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(limits=c(0,1)) +
  theme(legend.justification = c(0.05,0.95), legend.position = c(0.05,0.95))
#plot_UTR5_down


# plot up/total
plot_UTR5_up <- ggplot(UTR5) + 
  geom_point(aes(x = level, y = runmean(x=UTR5$up.total, k=5)), colour = "blue") +
  geom_point(aes(x = level, y = runmean(x=UTR5$up.total.1, k=5)), colour = "red") + 
  geom_point(aes(x = level, y = runmean(x=UTR5$fraction.targeted, k=5)), colour = "black") +
  facet_grid(~length, scales = "free") +
  ylab("upregulated / total") + 
  theme(axis.title=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), text = element_text(size=20)) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(limits=c(0,1))
#plot_UTR5_up



#######################################################################################################
# CDS
#######################################################################################################

# plot down/total
plot_CDS_down <- ggplot(CDS) + 
  geom_point(aes(x = level, y = runmean(x=CDS$down.total, k=5)), colour = "blue") +
  geom_point(aes(x = level, y = runmean(x=CDS$down.total.1, k=5)), colour = "red") + 
  geom_point(aes(x = level, y = runmean(x=CDS$fraction.targeted, k=5)), colour = "black") +
  facet_grid(~length, scales = "free") +
  ylab("downregulated / total") + 
  theme(axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(), text = element_text(size=20)) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(limits=c(0,1))
#plot_CDS_down


# plot up/total
plot_CDS_up <- ggplot(CDS) + 
  geom_point(aes(x = level, y = runmean(x=CDS$up.total, k=5)), colour = "blue") +
  geom_point(aes(x = level, y = runmean(x=CDS$up.total.1, k=5)), colour = "red") + 
  geom_point(aes(x = level, y = runmean(x=CDS$fraction.targeted, k=5)), colour = "black") +
  facet_grid(~length, scales = "free") +
  ylab("upregulated / total") + 
  theme(axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(), text = element_text(size=20)) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(limits=c(0,1))
#plot_CDS_up



#######################################################################################################
# 3' UTR
#######################################################################################################

# plot down/total
plot_UTR3_down <- ggplot(UTR3) + 
  geom_point(aes(x = level, y = runmean(x=UTR3$down.total, k=5)), colour = "blue") +
  geom_point(aes(x = level, y = runmean(x=UTR3$down.total.1, k=5)), colour = "red") + 
  geom_point(aes(x = level, y = runmean(x=UTR3$fraction.targeted, k=5)), colour = "black") +
  facet_grid(~length, scales = "free") +
  ylab("downregulated / total") + 
  theme(axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(), text = element_text(size=20)) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(limits=c(0,1))
#plot_UTR3_down


# plot up/total
plot_UTR3_up <- ggplot(UTR3) + 
  geom_point(aes(x = level, y = runmean(x=UTR3$up.total, k=5)), colour = "blue") +
  geom_point(aes(x = level, y = runmean(x=UTR3$up.total.1, k=5)), colour = "red") + 
  geom_point(aes(x = level, y = runmean(x=UTR3$fraction.targeted, k=5)), colour = "black") +
  facet_grid(~length, scales = "free") +
  ylab("upregulated / total") + 
  theme(axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(), text = element_text(size=20)) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(limits=c(0,1))
#plot_UTR3_up


#######################################################################################################
# arrange plots in 1 figure to export
#######################################################################################################

grid.arrange(plot_UTR5_down, plot_CDS_down, plot_UTR3_down, plot_UTR5_up, plot_CDS_up, plot_UTR3_up, nrow = 2)



#######################################################################################################
# arrange plots up/down-separated in 1 figure to export
#######################################################################################################


grid.arrange(plot_UTR5_down, plot_CDS_down, plot_UTR3_down, nrow = 1)

plot_UTR5_up <- ggplot(UTR5) + 
  geom_point(aes(x = level, y = runmean(x=UTR5$up.total, k=5), colour = "targeted")) +
  geom_point(aes(x = level, y = runmean(x=UTR5$up.total.1, k=5), colour = "not targeted")) +
  scale_colour_manual(values=c("red", "blue")) +
  geom_point(aes(x = level, y = runmean(x=UTR5$fraction.targeted, k=5), shape = "fraction targeted"), colour = "black") +
  scale_shape_manual(values=c(16)) +
  facet_grid(~length, scales = "free") +
  labs(y = "upregulated / total", colour = NULL, shape = NULL) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), text = element_text(size=20)) + 
  ggtitle("5' UTR") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(limits=c(0,1)) +
  theme(legend.justification = c(0.05,0.95), legend.position = c(0.05,0.95))
#plot_UTR5_up

grid.arrange(plot_UTR5_up, plot_CDS_up, plot_UTR3_up, nrow = 1)
