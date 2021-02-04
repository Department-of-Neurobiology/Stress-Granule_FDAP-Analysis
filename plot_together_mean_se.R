###############################################
##Nataliya Trushina, 2020##
##Stress granule FDAP analysis##

#######
#Libraries
#######
library(tidyr)
library(ggplot2)
library(Rmisc)
library(stringr)
#######
#Set working directory to the file location in RStudio
#######
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 
#######
#Plot together multiple conditions
#######
all_filenames_plot_together <- Sys.glob(file.path("for_plotting*"))
data_merge <- data.frame()
for (i in all_filenames_plot_together){  
  x <- read.table(i, sep = ";",header = TRUE)
  data_merge <- rbind(data_merge, x)
}
data_merge$construct <- as.factor(data_merge$construct)
plot_mult <- ggplot(data_merge, aes(x = time, y = measured_intensity)) + 
  geom_line(aes(y=measured_intensity, color=construct)) +
  geom_ribbon(data = data_merge, aes(ymin=CI_lower,ymax=CI_upper,fill=construct),alpha=0.4) +
  theme_classic() +
  scale_colour_brewer(palette = "Dark2") + #max 8 conditions
  scale_fill_brewer(palette = "Dark2") +
  ggtitle("Plot of mean measured intensity +-SE") +
  xlab("Time (s)") + 
  ylab("Normalized intensity") +
  scale_x_continuous(breaks=seq(0, 300, 10)) +
  scale_y_continuous(limits = c(0,1.1),breaks=seq(0, 1.1, 0.1))
svg("mean_se.svg", width=5, height=3.5)
plot_mult
dev.off()
p <- ggplotly(
  p = plot_mult
)
htmlwidgets::saveWidget(as_widget(p), paste("mean_se.html"))