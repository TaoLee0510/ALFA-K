setwd("~/projects/008_birthrateLandscape/ALFA-K/")
library(ggplot2)
source("utils/landscape_functions.R")

#set.seed(42)
nchrom <- 2
Nwaves <- 10
wavelength <- 1

nsims <- 500
npoints <- 10

y <- do.call(rbind,lapply(1:nsims,function(i){
  pk <- gen_rf_landscape(founder = sample(1:10,nchrom,replace = T),Nwaves = Nwaves,wavelength = wavelength)
  do.call(rbind,lapply(1:npoints,function(j){
    k1 <- runif(2,min=0,max=10)
    k2 <- runif(2,min=0,max=10)
    
    pr <- get_rf_fitness(k1,pk)*get_rf_fitness(k2,pk)
    d <- sqrt(sum((k1-k2)^2))
    data.frame(d,pr)
  }))
}))
y$d <- round(y$d,digits=2)
y <- aggregate(list(pr=y$pr),by=list(d=y$d),mean)

plot(y$d,y$pr)
