setwd("~/projects/008_birthrateLandscape/ALFA-K/")
source("utils/ALFA-K.R")
x <- readRDS("data/cellLines/02_optim_data/hTERTa.Rds")
x$x <- x$x[,1:(ncol(x$x)-1)]
