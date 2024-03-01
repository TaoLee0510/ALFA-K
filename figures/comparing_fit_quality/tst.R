setwd("~/projects/008_birthrateLandscape/ALFA-K/")
source("utils/ALFA-K.R")
x <- readRDS("data/cellLines/02_optim_data/hTERTb.Rds")
y <- alfak(x)
