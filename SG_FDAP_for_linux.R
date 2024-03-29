#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
print(paste("Condition name:", args[1], ", number of rows to subset: ", args[2], ", analysis type: ", args[3])) #https://www.r-bloggers.com/2015/09/passing-arguments-to-an-r-script-from-command-lines/


#Prerequisites for using Ubuntu subsystem on Windows:
# sudo apt install r-base-core
# sudo apt install pandoc
# sudo -i R
# install.packages(c("ggplot2", "minpack.lm", "tidyr", "Rmisc", "plotly", "RColorBrewer")) #https://stackoverflow.com/questions/44013764/unable-to-install-plotly-in-r


#Examples of running the script:
# Rscript --vanilla SG_FDAP_for_linux.R G3BP1wt 300 FIJI
# Rscript --vanilla SG_FDAP_for_linux.R G3BP1wt 300 NIS

###############################################
##Nataliya Trushina, 2020##
##Stress granule FDAP analysis##

#Requirements: output folder with csv files with measured intensities in the assigned ROI and archives with assigned ROIs coordinates.
#OR txt files with measured intensities in the photoactivation ROI from NIS Elements Software (Nikon) software.

#Plot data without extrapolation and fitting
#Plot biphasic exponential fits (also with data points - measured intensities in the ROI).
#Output the fitting parameters and the parameters for the goodness of fit.
###############################################

#######
#Libraries
#######
library(ggplot2)
library(minpack.lm)
library(tidyr)
library(Rmisc)
library(plotly)
library(RColorBrewer)

#######
#Additional settings
#######
construct_name <- args[1] #input the condition name (e.g., analyzed construct)
subset_nrow <-  as.numeric(args[2]) #input different value for easy subsetting by frame number (the number should not be longer than the shortest measurement acquisition!)
input_type <- as.character(args[3]) #choose the analysis type, write "NIS" for NIS-elements output (different versions of the programme might require small changes in file reading)
col_gradient <- colorRampPalette(c("#7CC17B", "#074c00")) #choose multiple colors for gradient

#######
#Functions
#######
create_empty_table <- function(num_rows, num_cols) {
  frame <- data.frame(matrix(NA, nrow = num_rows, ncol = num_cols))
  return(frame)
}

#######
#Read files with intensity measurements for plotting raw data without normalization
#######
rawdatalist = list() #create a list for data assembling
if (input_type == "NIS") {
  #### FOR NIS OUTPUT
  #find all txt files in the working directory
  filenames <- Sys.glob(file.path("*.txt"))
  for (i in filenames){
    x <- read.table(i, sep = "\t", header = TRUE, fill = TRUE)
    x <- x[!apply(x == "", 1, all),] #removes the empty row from some problematic files #to also remove NA: data <- data[!apply(is.na(data) | data == "", 1, all),]
    Int0_background <- x[,6][1] #the background intensity is in the first row of the last column
    x[,6] <- x[,6] - Int0_background #substract the background intensity from all the values in the last column
    x <- x[-c(1),] #delete the first row with zero value
    x <-x[1:subset_nrow,] #easy subsetting
    nrow_input <- nrow(x)
    name <- gsub(".txt","",i)
    print(name)
    dat <- data.frame(matrix(ncol=0,nrow=nrow_input+1)) #depends on subsetting
    dat <- x[6]
    rawdatalist[name] <- dat
  }
} else {
  #### FOR FIJI OUTPUT
  #find all csv files in the working directory
  filenames <- Sys.glob(file.path("*RawIntDen1.csv")) 
  for (i in filenames){
    #i <- "171023_D2_cell07_GG80_RawIntDen1.csv"
    x <- read.table(i, sep = ",", header = TRUE, fill = TRUE, row.names = "X")
    Int0_background <- x$RawIntDen1[1] #the background intensity is in the value corresponding to the first
    x$RawIntDen1 <- x$RawIntDen1 - Int0_background #substract the background intensity from all the values in the last column
    x <- x[-c(1),] #delete the first row with zero value
    x <-x[1:subset_nrow,] #easy subsetting
    nrow_input <- nrow(x) 
    name <- gsub("_RawIntDen1.csv","",i)
    print(name)
    dat <- data.frame(matrix(ncol=0,nrow=nrow_input+1)) #depends on subsetting
    dat <- x[2]
    rawdatalist[name] <- dat
  }
}
#assemble data
raw_data = as.data.frame(do.call(rbind, rawdatalist))
write.table(raw_data, "raw_data.csv", sep = ";",dec = '.', row.names = TRUE, col.names = FALSE)
raw_data_transpose <- as.data.frame(t(as.matrix(raw_data)))
nrow_plot <- nrow(raw_data_transpose) 
plot_Int <- cbind(raw_data_transpose, "type"="Int","time"=1:nrow_plot)
long_plot_Int <- gather(plot_Int, cell, measured_intensity, -c(time:type))
long_plot_Int$cell <- as.factor(long_plot_Int$cell)
#plot raw intensity measurements
cbp <- col_gradient(nlevels(long_plot_Int$cell))
plot_Int_raw <- ggplot(long_plot_Int, aes(x = time, y = measured_intensity)) + 
  geom_line(aes(y=measured_intensity, color=cell)) +
  scale_color_manual(values = cbp) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.position = "none") +
  ggtitle("Plot of raw data with background subtraction") +
  xlab("Time (s)") + 
  ylab("Intensity") +
  labs(color = "File name") +
  scale_x_continuous(breaks=seq(0, nrow_plot, 20)) +
  scale_y_continuous(breaks = pretty(long_plot_Int$measured_intensity, n = 5))
svg("raw_data.svg",  width=5, height=3.5)
plot_Int_raw
dev.off()
p <- ggplotly(
  p = plot_Int_raw
)
htmlwidgets::saveWidget(as_widget(p), paste("raw_data.html"))

#######
#Settings for analysis
#######
datalist = list()
dir.create("jpeg_fits")
dir.create("jpeg_fits_extrapolation_check")
cell_counter=1
error_counter_1=0
error_counter=0
normdatalist = list()
normdatapoints = list()

#######
#Loop through intensity measurement files, 
# extract the background intensity, 
# extrapolate the first data point,
# normalize,
# try biphasic exponential decay fit,
# if cannot be fitted - try monophasic exponential decay fit.
#Plot successful fits separately.
#Collect transformed data for further plotting.
#######
if (input_type == "NIS") {
  for (i in filenames){
    #### FOR NIS OUTPUT
    x <- read.table(i, sep = "\t", header = TRUE, fill = TRUE)
    x <- x[!apply(x == "", 1, all),]
    name <- gsub(".txt","",i)
    print(name)
    Int0_background <- x[,6][1]
    x[,6] <- x[,6] - Int0_background
    x <-x[1:subset_nrow,] #easy subsetting from additional settings
    time <- 1:nrow(x) -1
    y <- x[,6]

    str_1 <- cbind(time,y)[-c(1),] #delete the first row with zero values
    str_1 <- as.data.frame(str_1)
    n <- nrow(str_1)
    data_to_fit_1 <- structure(list(x = str_1$time, y=str_1$y), class = "data.frame", row.names = c(NA, -n), .Names = c("x", "y"))
    #try extrapolating the first data point
    result_1 <- try(fit_1 <- nlsLM(y ~ y0 + a1*exp(-x/b1) + a2*exp(-x/b2), data=data_to_fit_1, 
                                   control = list(maxiter = 100),
                                   lower=c(offset=0, a1=0, b1=0, a2=0, b2=0),
                                   start=list(y0=max(y),
                                              a1=max(y), 
                                              b1=1, 
                                              a2=max(y), 
                                              b2=100)))
    if (class(result_1) == "try-error") {
      error_counter_1<-error_counter_1+1
      print(c("Could not be extrapolated:", name))
      next
    }
    else {
      print(coef(fit_1))
      #plot fits for extrapolation check 
      jpeg(paste("jpeg_fits_extrapolation_check/",paste(paste(name,'extrapolation_check',sep="_"),"jpg",sep="."),sep=""))
      plot(y~time, data = str_1, main ="Extrapolation parameters check", xlab = "Time (s)", ylab = "Intensity",ylim=c(0,max(x[,6])))
      curve(predict(fit_1, newdata = data.frame(x)), col = "pink", add = TRUE)
      dev.off()
      
      Int0_extrap <- coef(fit_1)[1] + coef(fit_1)[2] + coef(fit_1)[4]
      str <- str_1
      str[nrow(str) + 1,] = c(0,Int0_extrap)
      str$y <- str$y/Int0_extrap
      data_to_fit <- structure(list(x = str$time, y=str$y), class = "data.frame", row.names = c(NA, -n-1), .Names = c("x", "y"))
      #try biphasic exponential decay fit with different initial parameters for each measurement
      result <- try(fit <- nlsLM(y ~ offset + AF*exp(-(x/tf)) + AS*exp(-(x/ts)), data=data_to_fit, 
                                 start = list(offset=max(y), AF=max(y), tf=100, AS=max(y), ts=1),
                                 control = list(maxiter = 100), lower=c(offset=0, AF=0, tf=0, AS=0, ts=0)))
      #OR try biphasic exponential decay fit with set initial parameters for each measurement, for example:
      #start = list(offset=0.1, AF=0.4, tf=1000, AS=0.4, ts=300)
      if (class(result) == "try-error") {
        error_counter<-error_counter+1
        print(c("Could not be fitted with biphasic exponential decay:",name))
        #try monophasic exponential decay fit with different initial parameters for each measurement
        fit_one_phase <- nlsLM(y ~ offset + AS*exp(-(x/ts)), data=data_to_fit, 
                               start = list(offset=max(y), AS=max(y), ts=30), control = list(maxiter = 100), lower=c(offset=0, AS=0, ts=0))
        #plot
        jpeg(paste("jpeg_fits/",paste(paste(name,'one_phase',sep="_"),"jpg",sep="."),sep=""))
        nf <- layout(matrix(c(1,2,3,3), ncol = 2, byrow = TRUE))
        layout.show(nf)
        par(cex=0.8, mai=c(0.1,0.1,0.2,0.1))
        par(mar=c(2,2,2,0.5))
        plot(residuals(fit_one_phase),main="Residuals", ylim=c(-0.3,0.3))
        abline(0, 0, col = "blue", lty = 2) 
        par(mar=c(2,2,2,1))
        qqnorm(residuals(fit_one_phase))
        par(mar=c(2,2,2,1))
        plot(y~x, data = data_to_fit, main =paste(i,sep="_"), xlab = "Time (s)", ylab = "Intensity", ylim=c(0,1.1))
        curve(predict(fit_one_phase, data.frame(x)), col = "blue", add = TRUE)
        dev.off()      
        #for assembling output table
        empty_df <- create_empty_table(12,1)
        empty_df[1,] <- coef(fit_one_phase)[1]
        empty_df[4,] <- coef(fit_one_phase)[2]
        empty_df[5,] <- coef(fit_one_phase)[3]
        empty_df[6,] <- coef(fit_one_phase)[1] + coef(fit_one_phase)[2]
        empty_df[7,] <- sigma(fit_one_phase)
        empty_df[8,] <-cor(str$y,predict(fit_one_phase))
        empty_df
        datalist[name] <-  as.list(empty_df)
        next
      }
      else {
        #if the biphasic fit was possible the monophasic fit is made only for comparison, for example, anova(fit, fit_one_phase_compare)
        fit_one_phase_compare <- nlsLM(y ~ offset + AS*exp(-(x/ts)), data=data_to_fit, 
                                       start = list(offset=max(y), AS=max(y), ts=30), control = list(maxiter = 100), lower=c(offset=0, AS=0, ts=0))
        #plot
        jpeg(paste("jpeg_fits/",paste(paste(name,sep="_"),"jpg",sep="."),sep=""))
        nf <- layout(matrix(c(1,2,3,3), ncol = 2, byrow = TRUE))
        layout.show(nf)
        par(cex=0.8) #make labels and margins smaller - #, mai=c(0.1,0.1,0.2,0.1)
        par(mar=c(2,2,2,0.5))
        plot(residuals(fit),main="Residuals", ylim=c(-0.3,0.3))
        abline(0, 0, col = "red", lty = 2) 
        par(mar=c(2,2,2,1))
        #qqnorm(residuals(fit))
        plot(residuals(fit_one_phase_compare),main="Residuals one-phase", ylim=c(-0.3,0.3))
        abline(0, 0, col = "blue", lty = 2) 
        par(mar=c(2,2,2,1))
        plot(y~x, data = data_to_fit, main =paste(i,sep="_"), xlab = "Time (s)", ylab = "Intensity", ylim=c(0,1))
        curve(predict(fit, data.frame(x)), col = "red", add = TRUE)
        dev.off()
        #check goodness of fit 
        df_chi <- data.frame(data_to_fit$y,predict(fit))
        print(chisq.test(df_chi$data_to_fit.y,df_chi$predict.fit., simulate.p.value = TRUE))
        print(coef(fit))
        coefs <- as.data.frame(coef(fit))
        coefs[nrow(coefs) + 1,] <- coef(fit)[1] + coef(fit)[2] + coef(fit)[4]
        print(coef(fit)[1] + coef(fit)[2] + coef(fit)[4])
        row.names(coefs)[nrow(coefs)] <- "A_total"
        colnames(coefs) <- name
        #for assembling output table
        empty_df <- create_empty_table(12,1)
        empty_df[1:6,] <- as.list(coefs)
        empty_df[7,] <- sigma(fit)
        empty_df[8,] <- cor(str$y,predict(fit))
        empty_df[9,] <- sigma(fit_one_phase_compare)
        empty_df[10,] <- cor(str$y,predict(fit_one_phase_compare))
        empty_df[11,] <- sigma(fit) - sigma(fit_one_phase_compare)
        empty_df[12,] <- cor(str$y,predict(fit)) - cor(str$y,predict(fit_one_phase_compare))
        datalist[name] <-  as.list(empty_df)
        #for plotting normalized fitting curves
        normdat <- data.frame(matrix(ncol=0,nrow=subset_nrow))
        normdat <- as.data.frame(predict(fit))
        normdatalist[name] <- normdat
        #for plotting normalized points
        normdatp <- data.frame(matrix(ncol=0,nrow=subset_nrow))
        normdatp <- data_to_fit[2]
        normdatapoints[name] <- normdatp
      }
    }
  }
} else {
  for (i in filenames){
    #### FOR FIJI OUTPUT
    x <- read.table(i, sep = ",", header = TRUE, fill = TRUE, row.names = "X")
    name <- gsub("_RawIntDen1.csv","",i)
    print(i)
    Int0_background <- x$RawIntDen1[1] #the background intensity is in the value corresponding to the first
    x$RawIntDen1 <- x$RawIntDen1 - Int0_background #substract the background intensity from all the values in the last column
    
    x <-x[1:subset_nrow,] #easy subsetting
    time <- 1:nrow(x) - 1
    y <- x$RawIntDen1 
  
    str_1 <- cbind(time,y)[-c(1),] #delete the first row with zero values
    str_1 <- as.data.frame(str_1)
    n <- nrow(str_1)
    data_to_fit_1 <- structure(list(x = str_1$time, y=str_1$y), class = "data.frame", row.names = c(NA, -n), .Names = c("x", "y"))
    #try extrapolating the first data point
    result_1 <- try(fit_1 <- nlsLM(y ~ y0 + a1*exp(-x/b1) + a2*exp(-x/b2), data=data_to_fit_1, 
                                   control = list(maxiter = 100),
                                   lower=c(offset=0, a1=0, b1=0, a2=0, b2=0),
                                   start=list(y0=max(y),
                                              a1=max(y), 
                                              b1=1, 
                                              a2=max(y), 
                                              b2=100)))
    if (class(result_1) == "try-error") {
      error_counter_1<-error_counter_1+1
      print(c("Could not be extrapolated:", name))
      next
    }
    else {
      print(coef(fit_1))
      #plot fits for extrapolation check 
      jpeg(paste("jpeg_fits_extrapolation_check/",paste(paste(name,'extrapolation_check',sep="_"),"jpg",sep="."),sep=""))
      plot(y~time, data = str_1, main ="Extrapolation parameters check", xlab = "Time (s)", ylab = "Intensity",ylim=c(0,max(x[,2])))
      curve(predict(fit_1, newdata = data.frame(x)), col = "pink", add = TRUE)
      dev.off()
      
      Int0_extrap <- coef(fit_1)[1] + coef(fit_1)[2] + coef(fit_1)[4]
      str <- str_1
      str[nrow(str) + 1,] = c(0,Int0_extrap)
      str$y <- str$y/Int0_extrap
      data_to_fit <- structure(list(x = str$time, y=str$y), class = "data.frame", row.names = c(NA, -n-1), .Names = c("x", "y"))
      #try biphasic exponential decay fit with different initial parameters for each measurement
      result <- try(fit <- nlsLM(y ~ offset + AF*exp(-(x/tf)) + AS*exp(-(x/ts)), data=data_to_fit, 
                                 start = list(offset=max(y), AF=max(y), tf=100, AS=max(y), ts=1),
                                 control = list(maxiter = 100), lower=c(offset=0, AF=0, tf=0, AS=0, ts=0)))
      #OR try biphasic exponential decay fit with set initial parameters for each measurement, for example:
      #start = list(offset=0.1, AF=0.4, tf=1000, AS=0.4, ts=300)
      if (class(result) == "try-error") {
        error_counter<-error_counter+1
        print(c("Could not be fitted with biphasic exponential decay:",name))
        #try monophasic exponential decay fit with different initial parameters for each measurement
        fit_one_phase <- nlsLM(y ~ offset + AS*exp(-(x/ts)), data=data_to_fit, 
                               start = list(offset=max(y), AS=max(y), ts=30), control = list(maxiter = 100), lower=c(offset=0, AS=0, ts=0))
        #plot
        jpeg(paste("jpeg_fits/",paste(paste(name,'one_phase',sep="_"),"jpg",sep="."),sep=""))
        nf <- layout(matrix(c(1,2,3,3), ncol = 2, byrow = TRUE))
        layout.show(nf)
        par(cex=0.8, mai=c(0.1,0.1,0.2,0.1))
        par(mar=c(2,2,2,0.5))
        plot(residuals(fit_one_phase),main="Residuals", ylim=c(-0.3,0.3))
        abline(0, 0, col = "blue", lty = 2) 
        par(mar=c(2,2,2,1))
        qqnorm(residuals(fit_one_phase))
        par(mar=c(2,2,2,1))
        plot(y~x, data = data_to_fit, main =paste(i,sep="_"), xlab = "Time (s)", ylab = "Intensity", ylim=c(0,1.1))
        curve(predict(fit_one_phase, data.frame(x)), col = "blue", add = TRUE)
        dev.off()      
        #for assembling output table
        empty_df <- create_empty_table(12,1)
        empty_df[1,] <- coef(fit_one_phase)[1]
        empty_df[4,] <- coef(fit_one_phase)[2]
        empty_df[5,] <- coef(fit_one_phase)[3]
        empty_df[6,] <- coef(fit_one_phase)[1] + coef(fit_one_phase)[2]
        empty_df[7,] <- sigma(fit_one_phase)
        empty_df[8,] <-cor(str$y,predict(fit_one_phase))
        empty_df
        datalist[name] <-  as.list(empty_df)
        next
      }
      else {
        #if the biphasic fit was possible the monophasic fit is made only for comparison, for example, anova(fit, fit_one_phase_compare)
        fit_one_phase_compare <- nlsLM(y ~ offset + AS*exp(-(x/ts)), data=data_to_fit, 
                                       start = list(offset=max(y), AS=max(y), ts=30), control = list(maxiter = 100), lower=c(offset=0, AS=0, ts=0))
        #plot
        jpeg(paste("jpeg_fits/",paste(paste(name,sep="_"),"jpg",sep="."),sep=""))
        nf <- layout(matrix(c(1,2,3,3), ncol = 2, byrow = TRUE))
        layout.show(nf)
        par(cex=0.8) #make labels and margins smaller - #, mai=c(0.1,0.1,0.2,0.1)
        par(mar=c(2,2,2,0.5))
        plot(residuals(fit),main="Residuals", ylim=c(-0.3,0.3))
        abline(0, 0, col = "red", lty = 2) 
        par(mar=c(2,2,2,1))
        #qqnorm(residuals(fit))
        plot(residuals(fit_one_phase_compare),main="Residuals one-phase", ylim=c(-0.3,0.3))
        abline(0, 0, col = "blue", lty = 2) 
        par(mar=c(2,2,2,1))
        plot(y~x, data = data_to_fit, main =paste(i,sep="_"), xlab = "Time (s)", ylab = "Intensity", ylim=c(0,1))
        curve(predict(fit, data.frame(x)), col = "red", add = TRUE)
        dev.off()
        #check goodness of fit 
        df_chi <- data.frame(data_to_fit$y,predict(fit))
        print(chisq.test(df_chi$data_to_fit.y,df_chi$predict.fit., simulate.p.value = TRUE))
        print(coef(fit))
        coefs <- as.data.frame(coef(fit))
        coefs[nrow(coefs) + 1,] <- coef(fit)[1] + coef(fit)[2] + coef(fit)[4]
        print(coef(fit)[1] + coef(fit)[2] + coef(fit)[4])
        row.names(coefs)[nrow(coefs)] <- "A_total"
        colnames(coefs) <- name
        #for assembling output table
        empty_df <- create_empty_table(12,1)
        empty_df[1:6,] <- as.list(coefs)
        empty_df[7,] <- sigma(fit)
        empty_df[8,] <- cor(str$y,predict(fit))
        empty_df[9,] <- sigma(fit_one_phase_compare)
        empty_df[10,] <- cor(str$y,predict(fit_one_phase_compare))
        empty_df[11,] <- sigma(fit) - sigma(fit_one_phase_compare)
        empty_df[12,] <- cor(str$y,predict(fit)) - cor(str$y,predict(fit_one_phase_compare))
        datalist[name] <-  as.list(empty_df)
        #for plotting normalized fitting curves
        normdat <- data.frame(matrix(ncol=0,nrow=subset_nrow))
        normdat <- as.data.frame(predict(fit))
        normdatalist[name] <- normdat
        #for plotting normalized points
        normdatp <- data.frame(matrix(ncol=0,nrow=subset_nrow))
        normdatp <- data_to_fit[2]
        normdatapoints[name] <- normdatp
      }
    }
  }
}
  

#assemble normalized points
norm_data_points = as.data.frame(do.call(rbind, normdatapoints)) 
write.table(norm_data_points, "norm_data_points.csv", sep = ";",dec = '.', row.names = TRUE, col.names = FALSE)
norm_data_points_transpose <- as.data.frame(t(as.matrix(norm_data_points)))
plot_Int_points <- cbind(norm_data_points_transpose, "type"="Int","time"=1:subset_nrow)
plot_Int_points$time[subset_nrow] <- 0
plot_Int_points$time <- as.numeric(plot_Int_points$time) + 1
long_plot_Int_points <- gather(plot_Int_points, cell, measured_intensity, -c(time:type))
long_plot_Int_points$cell <- as.factor(long_plot_Int_points$cell)
plot_final_summary <- summarySE(long_plot_Int_points, measurevar="measured_intensity", groupvars=c("type","time"))
plot_final_summary$CI_lower <- plot_final_summary$measured_intensity + plot_final_summary$se
plot_final_summary$CI_upper <- plot_final_summary$measured_intensity - plot_final_summary$se
#assemble normalized fitting curves
norm_data = as.data.frame(do.call(rbind, normdatalist)) 
write.table(norm_data, "norm_data.csv", sep = ";",dec = '.', row.names = TRUE, col.names = FALSE)
norm_data_transpose <- as.data.frame(t(as.matrix(norm_data)))
plot_Int <- cbind(norm_data_transpose, "type"="Int","time"=1:subset_nrow)
plot_Int$time[subset_nrow] <- 0
plot_Int$time <- as.numeric(plot_Int$time) + 1
long_plot_Int <- gather(plot_Int, cell, measured_intensity, -c(time:type))
long_plot_Int$cell <- as.factor(long_plot_Int$cell)
#plot normalized fitting curves
cbp <- col_gradient(nlevels(long_plot_Int$cell))
plot_Int <- ggplot(long_plot_Int, aes(x = time, y = measured_intensity)) + 
  geom_line(aes(y=measured_intensity, color=cell)) +
  geom_ribbon(data = plot_final_summary, aes(ymin=CI_lower,ymax=CI_upper),alpha=0.4, fill="#7CC17B") +
  scale_color_manual(values = cbp) +
  scale_fill_manual(values = cbp) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.position = "none") +
  ggtitle("Plot of fitted normalized intensity measurements +-SEM") +
  xlab("Time (s)") + 
  ylab("Fitted normalized intensity") +
  labs(color = "File name") +
  scale_x_continuous(breaks=seq(0, subset_nrow, 20)) +
  scale_y_continuous(limits = c(-0.1,1.1),breaks=seq(-0.1, 1.1, 0.1))
svg("norm_data.svg",  width=5, height=3.5)
plot_Int
dev.off()
p <- ggplotly(
  p = plot_Int
)
htmlwidgets::saveWidget(as_widget(p), paste("norm_data.html"))
#assemble fitting parameters
big_data = as.data.frame(do.call(rbind, datalist))
big_data <- setNames(big_data,c("offset","AS","ts","AF","tf","A_total","Residual standard error (RSE)","Correlation coefficient (cor)","One-phase RSE","One-phase cor","RSE difference", "Cor difference"))
df_final <- big_data 
#calculate mean values for only successful biphasic fits
df_final <- rbind(big_data, colMeans(na.omit(big_data)))
row.names(df_final)[nrow(df_final)] <- "Mean_biphas"
write.table(df_final, "all_fit_coefs.csv", sep = ";",dec = '.', row.names = TRUE, col.names = NA)
#plot normalized fitting curves with points
plot_Int_points <- ggplot(long_plot_Int, aes(x = time, y = measured_intensity)) + 
  geom_line(aes(y=measured_intensity, color=cell)) +
  geom_point(size = 0.5,data = long_plot_Int_points, aes(x = time, y = measured_intensity, color=cell)) +
  geom_ribbon(data = plot_final_summary, aes(ymin=CI_lower,ymax=CI_upper),alpha=0.4, fill="#7CC17B") +
  scale_color_manual(values = cbp) +
  scale_fill_manual(values = cbp) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.position = "none") +
  ggtitle("Plot of fitted normalized intensity measurements +-SEM") +
  xlab("Time (s)") + 
  ylab("Fitted normalized intensity") +
  labs(color = "File name") +
  scale_x_continuous(breaks=seq(0, subset_nrow, 20)) +
  scale_y_continuous(limits = c(-0.1,1.1),breaks=seq(-0.1, 1.1, 0.1))
svg("norm_data_with_points.svg",  width=5, height=3.5)
plot_Int_points
dev.off()
p <- ggplotly(
  p = plot_Int_points
)
htmlwidgets::saveWidget(as_widget(p), paste("norm_data_with_points.html"))

#final check
print(c("Number of curves that could not be extrapolated:", error_counter_1))
print(c("Number of curves that could not be fitted with biphasic exponential decay:", error_counter))

#######
#Plot mean curves when multiple conditions are analyzed (check file paths for input files)
#######
long_plot_Int_points$construct <- construct_name
plot_final_summary <- summarySE(long_plot_Int_points, measurevar="measured_intensity", groupvars=c("construct","time"))
plot_final_summary$CI_lower <- plot_final_summary$measured_intensity + plot_final_summary$se
plot_final_summary$CI_upper <- plot_final_summary$measured_intensity - plot_final_summary$se
write.table(plot_final_summary, paste("../../for_plotting_mean_and_se_", construct_name, ".csv", sep=""), sep = ";", row.names = FALSE)
all_filenames_plot_together <- Sys.glob(file.path("../../for_plotting*"))
data_merge <- data.frame()
for (i in all_filenames_plot_together){  
  x <- read.table(i, sep = ";",header = TRUE)
  data_merge <- rbind(data_merge, x)
}
data_merge$construct <- as.factor(data_merge$construct)
#plot together multiple conditions
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
svg("../../mean_se.svg", width=5, height=3.5)
plot_mult
dev.off()
p <- ggplotly(
  p = plot_mult
)
htmlwidgets::saveWidget(as_widget(p), paste("../../mean_se.html"))
