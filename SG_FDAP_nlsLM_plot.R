###############################################
##Nataliya Trushina, 2020##
##Stress granule FDAP analysis##

#Requirements: txt files from nikon software intensity measurement after photoactivation.

##Version 2##
#Plot all data (without extrapolation and fitting)
#Plot all biphasic exponential fits
#Plot fits with data points
###############################################

#######
#Libraries
#######
library(ggplot2)
library(minpack.lm)
library(tidyr)
library(Rmisc)

#######
#Set working directory to the file location
#######
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 

#######
#Read the txt files in the working directory
#######
filenames <- Sys.glob(file.path("*.txt"))  #parameters
filenames

#######
#Plotting raw data
#######
rawdatalist = list()
for (i in filenames){
  #i <- "20200217_D01_C07_HV60.txt"
  x <- read.table(i, sep = "\t", header = TRUE, fill = TRUE)
  Int0_background <- x[,6][1]
  x[,6] <- x[,6] - Int0_background
  x <- x[-c(1),]
  x
  name <- gsub(".txt","",i)
  print(name)
  #print(x[,6])
  dat <- data.frame(matrix(ncol=0,nrow=302))
  dat <- x[6]
  rawdatalist[name] <- dat
}
raw_data = as.data.frame(do.call(rbind, rawdatalist))
write.table(raw_data, "raw_data.csv", sep = ";",dec = '.', row.names = TRUE, col.names = FALSE)
raw_data_transpose <- as.data.frame(t(as.matrix(raw_data)))
plot_Int <- cbind(raw_data_transpose, "type"="Int","time"=1:nrow(raw_data_transpose))
long_plot_Int <- gather(plot_Int, cell, measured_intensity, -c(time:type))

pdf("raw_data.pdf",  width=14, height=6)
ggplot(long_plot_Int, aes(x = time, y = measured_intensity)) + 
  geom_line(aes(y=measured_intensity, color=cell)) +
  #geom_ribbon(data = plot_together, aes(ymin=CI_lower,ymax=CI_upper,fill=construct),alpha=0.4) +
  theme_classic() +
  ggtitle("Plot of raw data") +
  xlab("Time (s)") + 
  ylab("Intensity") +
  scale_x_continuous(breaks=seq(0, 302, 50)) #+
  #scale_y_continuous(limits = c(0,1.1),breaks=seq(0, 1.1, 0.1))
dev.off()


#######
#Set initial conditions
#######
datalist = list()
dir.create("jpeg_fits")
cell_counter=1
error_counter_1=0
error_counter=0
normdatalist = list()
normdatapoints = list()

#######
#Loop through txt files, extract the background intensity, 
#extrapolation for the first datapoint,
#normalization,
#biphasic exponential decay fit
#if not fitted - tries monophasic
#######
for (i in filenames){
  #i <- "20200211_D03_C02_HV65.txt"
  x <- read.table(i, sep = "\t", header = TRUE, fill = TRUE)
  name <- gsub(".txt","",i)
  print(name)
  Int0_background <- x[,6][1]
  x[,6] <- x[,6] - Int0_background
  time <- 1:nrow(x) -1
  y <- x[,6]
  #colnames(y) <- "y"
  str_1 <- cbind(time,y)[-c(1),]
  str_1 <- as.data.frame(str_1)
  n <- nrow(str_1)
  data_to_fit_1 <- structure(list(x = str_1$time, y=str_1$y), class = "data.frame", row.names = c(NA, -n), .Names = c("x", "y"))
  data_to_fit_1
  
  #https://stackoverflow.com/questions/33265467/nls-troubles-missing-value-or-an-infinity-produced-when-evaluating-the-model
  #one possible way to make different start for each measurement
  result_1 <- try(fit_1 <- nlsLM(y ~ y0 + a1*exp(-x/b1) + a2*exp(-x/b2), data=data_to_fit_1, 
                                 control = list(maxiter = 100),
                                 start=list(y0=max(y),
                                            a1=max(y), 
                                            b1=1, 
                                            a2=max(y), 
                                            b2=100)))
  #result_1 <- try(fit_1 <- nlsLM(y ~ y0 + a1*exp(-x/b1) + a2*exp(-x/b2), data=data_to_fit_1, start = list(y0=1000000, a1=100000, b1=1, a2=10000, b2=100), control = list(maxiter = 100)))
  
  if (class(result_1) == "try-error") {
    error_counter_1<-error_counter_1+1
    print(c("Could not be extrapolated:", name))
    next
  }
  else {
    Int0_extrap <- coef(fit_1)[1] + coef(fit_1)[2] + coef(fit_1)[4]
    
    print(c("Test activation: ", data_to_fit_1$y[1]/predict(fit_1)[1]))
    
    str <- str_1
    str[nrow(str) + 1,] = c(0,Int0_extrap)
    str$y <- str$y/Int0_extrap
    data_to_fit <- structure(list(x = str$time, y=str$y), class = "data.frame", row.names = c(NA, -n-1), .Names = c("x", "y"))
    
    #IMPORTANT - might need to change the initial parameters
    result <- try(fit <- nlsLM(y ~ offset + AF*exp(-(x/tf)) + AS*exp(-(x/ts)), data=data_to_fit, 
                               start = list(offset=0.1, AF=0.4, tf=10000, AS=0.4, ts=300), control = list(maxiter = 100), lower=c(offset=0, AF=0, tf=0, AS=0, ts=0)))
    if (class(result) == "try-error") {
      error_counter<-error_counter+1
      print(c("Could not be fitted with biphasic exponential decay:",name))
      
      fit_one_phase <- nlsLM(y ~ offset + AS*exp(-(x/ts)), data=data_to_fit, 
                       start = list(offset=0.1, AS=0.4, ts=30), control = list(maxiter = 100), lower=c(offset=0, AS=0, ts=0))
      
      
      summary(fit)
      summary(fit_one_phase)
      anova(fit, fit_one_phase)
      confint.default(fit)
     
      jpeg(paste("jpeg_fits/",paste(paste(name,'one_phase',sep="_"),"jpg",sep="."),sep=""))
      plot(y~x, data = data_to_fit, main =paste(i,sep="_"), xlab = "Time (s)", ylab = "Intensity", ylim=c(0,1))
      curve(predict(fit_one_phase, data.frame(x)), col = "blue", add = TRUE)
      dev.off()
      
      create_empty_table <- function(num_rows, num_cols) {
        frame <- data.frame(matrix(NA, nrow = num_rows, ncol = num_cols))
        return(frame)
      }
      
      empty_df <- create_empty_table(6,1)
      empty_df[1,] <- coef(fit_one_phase)[1]
      empty_df[4,] <- coef(fit_one_phase)[2]
      empty_df[5,] <- coef(fit_one_phase)[3]
      empty_df[6,] <- coef(fit_one_phase)[1] + coef(fit_one_phase)[2]
      empty_df
      datalist[name] <-  as.list(empty_df)
      next
    }
    else {
      jpeg(paste("jpeg_fits/",paste(paste(name,sep="_"),"jpg",sep="."),sep=""))
      plot(y~x, data = data_to_fit, main =paste(i,sep="_"), xlab = "Time (s)", ylab = "Intensity", ylim=c(-0.1,1))
      curve(predict(fit, data.frame(x)), col = "red", add = TRUE)
      dev.off()
      
      print(coef(fit))
      coefs <- as.data.frame(coef(fit))
      coefs[nrow(coefs) + 1,] <- coef(fit)[1] + coef(fit)[2] + coef(fit)[4]
      print(coef(fit)[1] + coef(fit)[2] + coef(fit)[4])
      row.names(coefs)[nrow(coefs)] <- "A_total"
      coefs
      colnames(coefs) <- name
      coefs
      as.list((coefs))
      datalist[name] <-  as.list(coefs)
      
      normdat <- data.frame(matrix(ncol=0,nrow=302))
      normdat <- as.data.frame(predict(fit))
      normdatalist[name] <- normdat
      
      normdatp <- data.frame(matrix(ncol=0,nrow=302))
      normdatp <- data_to_fit[2]
      normdatapoints[name] <- normdatp
      
    }
  }
}

norm_data_points = as.data.frame(do.call(rbind, normdatapoints)) 
write.table(norm_data_points, "norm_data_points.csv", sep = ";",dec = '.', row.names = TRUE, col.names = FALSE)
norm_data_points_transpose <- as.data.frame(t(as.matrix(norm_data_points)))
plot_Int_points <- cbind(norm_data_points_transpose, "type"="Int","time"=1:nrow(norm_data_points_transpose))
plot_Int_points$time[301] <- 0
plot_Int_points$time <- as.numeric(plot_Int_points$time) + 1
long_plot_Int_points <- gather(plot_Int_points, cell, measured_intensity, -c(time:type))

plot_final_summary <- summarySE(long_plot_Int_points, measurevar="measured_intensity", groupvars=c("type","time"))
### change se to sd for further plotting
plot_final_summary$CI_lower <- plot_final_summary$measured_intensity + plot_final_summary$se
plot_final_summary$CI_upper <- plot_final_summary$measured_intensity - plot_final_summary$se



ggplot(long_plot_Int_points, aes(x = time, y = measured_intensity)) + 
  geom_point(aes(y=measured_intensity, color=cell)) +
  theme_classic() +
  ggtitle("Plot of normalized intensities") +
  xlab("Time (s)") + 
  ylab("Normalized intensity") +
  scale_x_continuous(breaks=seq(0, 302, 50))

norm_data = as.data.frame(do.call(rbind, normdatalist)) 
write.table(norm_data, "norm_data.csv", sep = ";",dec = '.', row.names = TRUE, col.names = FALSE)
norm_data_transpose <- as.data.frame(t(as.matrix(norm_data)))
plot_Int <- cbind(norm_data_transpose, "type"="Int","time"=1:nrow(norm_data_transpose))
plot_Int$time[301] <- 0
plot_Int$time <- as.numeric(plot_Int$time) + 1
long_plot_Int <- gather(plot_Int, cell, measured_intensity, -c(time:type))



pdf("norm_data.pdf",  width=14, height=6)
ggplot(long_plot_Int, aes(x = time, y = measured_intensity)) + 
  geom_line(aes(y=measured_intensity, color=cell)) +
  geom_ribbon(data = plot_final_summary, aes(ymin=CI_lower,ymax=CI_upper),alpha=0.4) +
  theme_classic() +
  ggtitle("Plot of fitted normalized intensities and mean +-se by time") +
  xlab("Time (s)") + 
  ylab("Fitted normalized intensity") +
  scale_x_continuous(breaks=seq(0, 302, 50)) #+
scale_y_continuous(limits = c(-0.1,1.1),breaks=seq(-0.1, 1.1, 0.1))
dev.off()

error_counter_1
error_counter
big_data = as.data.frame(do.call(rbind, datalist))
big_data <- setNames(big_data,c("offset","AF","tf","AS","ts","A_total"))
dfFinal <- big_data #dfFinal <- rbind(big_data, colMeans(na.omit(big_data)))
#row.names(dfFinal)[nrow(dfFinal)] <- "Mean"
write.table(dfFinal, "all_fit_coefs.csv", sep = ";",dec = '.', row.names = TRUE, col.names = NA)


pdf("norm_data_with_points.pdf",  width=20, height=14)
ggplot(long_plot_Int, aes(x = time, y = measured_intensity)) + 
  geom_line(aes(y=measured_intensity, color=cell)) +
  geom_point(data = long_plot_Int_points, aes(x = time, y = measured_intensity, color=cell)) +
  geom_ribbon(data = plot_final_summary, aes(ymin=CI_lower,ymax=CI_upper),alpha=0.4) +
  theme_classic() +
  ggtitle("Plot of fitted normalized intensities and mean +-se by time") +
  xlab("Time (s)") + 
  ylab("Normalized intensity") +
  scale_x_continuous(breaks=seq(0, 302, 50)) +
  scale_y_continuous(limits = c(-0.1,1.1),breaks=seq(-0.1, 1.1, 0.1))
dev.off()
