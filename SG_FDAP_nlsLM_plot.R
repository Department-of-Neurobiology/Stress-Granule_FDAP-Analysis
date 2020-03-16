setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 
library(ggplot2)
library(minpack.lm)
filenames <- Sys.glob(file.path("*.txt"))  #parameters
filenames
datalist = list()
dir.create("jpeg_fits")
cell_counter=1
error_counter_1=0
error_counter=0

for (i in filenames){
  #i <- "20191023_D3_C03_GG100.txt"
  x <- read.table(i, sep = "\t", header = TRUE, fill = TRUE)
  name <- gsub(".txt","",i)
  print(name)
  Int0_bg <- x[,6][1]
  x[,6] <- x[,6] - Int0_bg
  time <- 1:nrow(x) -1
  y <- x[,6]
  #colnames(y) <- "y"
  str_1 <- cbind(time,y)[-c(1),]
  str_1 <- as.data.frame(str_1)
  n <- nrow(str_1)
  data_to_fit_1 <- structure(list(x = str_1$time, y=str_1$y), class = "data.frame", row.names = c(NA, -n), .Names = c("x", "y"))
  data_to_fit_1
  
  #https://stackoverflow.com/questions/33265467/nls-troubles-missing-value-or-an-infinity-produced-when-evaluating-the-model
  #another possible way to make different start for each measurement
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
    next
  }
  else {
    Int0_extrap <- coef(fit_1)[1] + coef(fit_1)[2] + coef(fit_1)[4]
    str <- str_1
    str[nrow(str) + 1,] = c(0,Int0_extrap)
    str$y <- str$y/Int0_extrap
    data_to_fit <- structure(list(x = str$time, y=str$y), class = "data.frame", row.names = c(NA, -n-1), .Names = c("x", "y"))
  
    result <- try(fit <- nlsLM(y ~ offset + AF*exp(-(x/tf)) + AS*exp(-(x/ts)), data=data_to_fit, 
                               start = list(offset=0.1, AF=0.4, tf=100, AS=0.4, ts=300), control = list(maxiter = 100), lower=c(offset=0, AF=0, tf=0, AS=0, ts=0)))
    if (class(result) == "try-error") {
      error_counter<-error_counter+1
      
      #fit_one_phase <- nlsLM(y ~ offset + AS*exp(-(x/ts)), data=data_to_fit, 
      #                 start = list(offset=0.1, AS=0.4, ts=30), control = list(maxiter = 100), lower=c(offset=0, AS=0, ts=0))
      
      #jpeg(paste("jpeg_fits/",paste(paste(name,'one_phase',sep="_"),"jpg",sep="."),sep=""))
      #plot(y~x, data = data_to_fit, main =paste(i,sep="_"), xlab = "Time (s)", ylab = "Intensity", ylim=c(0,1))
      #curve(predict(fit_one_phase, data.frame(x)), col = "blue", add = TRUE)
      #dev.off()
      
      create_empty_table <- function(num_rows, num_cols) {
        frame <- data.frame(matrix(NA, nrow = num_rows, ncol = num_cols))
        return(frame)
      }
      empty_df <- create_empty_table(6,1)
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
    }
  }
}

error_counter_1
error_counter
big_data = as.data.frame(do.call(rbind, datalist))
big_data <- setNames(big_data,c("offset","AF","tf","AS","ts","A_total"))
dfFinal <- rbind(big_data, colMeans(na.omit(big_data)))
row.names(dfFinal)[nrow(dfFinal)] <- "Mean"
write.table(dfFinal, "all_fit_coefs.csv", sep = ";",dec = '.', row.names = TRUE, col.names = NA)
