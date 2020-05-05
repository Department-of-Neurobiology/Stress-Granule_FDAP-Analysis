library(tidyr)
library(ggplot2)
library(Rmisc)
library(stringr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 

all_filenames_plot_together <- Sys.glob(file.path("for_plotting*"))
data_merge <- data.frame()
for (i in all_filenames_plot_together){  
  x <- read.table(i, sep = ";",header = TRUE)
  data_merge <- rbind(data_merge, x)
}

svg("mean_se.svg", width=10, height=6) 
ggplot(data_merge, aes(x = time, y = measured_intensity)) + 
  geom_line(aes(y=measured_intensity, color=construct)) +
  geom_ribbon(data = data_merge, aes(ymin=CI_lower,ymax=CI_upper,fill=construct),alpha=0.4) +
  theme_classic() +
  ggtitle("Plot of mean measured intensity +-se by time") +
  xlab("Time (s)") + 
  ylab("Normalized intensity") +
  scale_x_continuous(breaks=seq(0, 110, 10)) +
  scale_y_continuous(limits = c(-0.1,1.1),breaks=seq(-0.1, 1.1, 0.1))
dev.off() 