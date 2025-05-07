library(Matrix)
library(igraph)
library(ggplot2)
library(RSpectra)
setwd("/share/lab_crd/M010_ALFAK_2023/ALFA-K")
source("utils/ALFA-KQ.R")
source("tests/prediction/prediction_functions.R")

procDir <- "data/salehi/alfak_outputs_V1a_procv0/"
files <- list.files(procDir, full.names = TRUE)
# Read all files and keep only those with 13 columns
x_list <- lapply(files, readRDS)
#x_list <- x_list[sapply(x_list, ncol) == 13]
x <- do.call(rbind, x_list)

# Filter data and select observation with maximum xv per file and per clean_id
x <- x[!is.na(x$xv) & x$xv > 0 & x$test_treat != "NaN", ]
x <- do.call(rbind, lapply(split(x, x$fi), function(df) df[df$xv == max(df$xv), ]))
rownames(x) <- NULL
x$clean_id <- sub("_l_\\d+_d1_\\d+_d2_\\d+$", "", x$fi)
x <- do.call(rbind, lapply(split(x, x$clean_id), function(df) head(df[df$xv == max(df$xv), ],1)))


df <- lapply(1:nrow(x),function(i){
  inpath <- paste0("data/salehi/alfak_outputs_V1a/minobs_",
                   x$min_obs[i],"/",x$fi[i],"/landscape.Rds")
  
  
  if(!file.exists(inpath)) return(NULL)
  lscape <- readRDS(inpath)
  
#  coords <- gen_coords(lscape$k)
#  dims <- rep(nrow(lscape), 2)
  
  p <- c(0.0008,0.0009,0.001,0.002,0.003,0.004,0.005,0.006)
  res <- do.call(cbind,pbapply::pblapply(p,function(p0){
    find_steady_state(lscape,p0)
  }))
  colnames(res) <- p
  res <- res[order(rowSums(res),decreasing=T),]
  print(c(i,x$fi[i]))
  print(head(round(res,digits = 2)))
  
  dfi <- data.frame(fi=x$fi[i],nmax=length(unique(apply(res,2,which.max))))
  print(dfi)
  return(res)
})

names(df) <- x$fi

saveRDS(df,"figures/misc/data/misseg_rate_test.Rds")

