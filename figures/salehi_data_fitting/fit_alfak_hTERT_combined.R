setwd("~/projects/008_birthrateLandscape/ALFA-K/")


min_obs <- 10

dir <- "data/salehi/alfak_inputs/"
ff <- list.files(dir)

ids <- c("p53--a_X57","p53--b_X55" )

x <- lapply(ids, function(idi) readRDS(paste0(dir,idi,".Rds")))

n_obs <- lapply(x, function(xi) table(rowSums(xi$x)))
n_obs <- as.numeric(names(unlist(n_obs)))

n_obs <- unique(n_obs[n_obs>min_obs])
n_obs <- n_obs[order(n_obs,decreasing=T)]

library(parallel)
cl <- makeCluster(getOption("cl.cores", 3))
clusterExport(cl, "x")

fits <- parLapplyLB(cl=cl,X = n_obs,fun =  function(ni){
  tryCatch({
  source("utils/ALFA-K.R")
    ##sometimes getting the error:
    ## Error in qr.default(Tmatrix) : 
    ## NA/NaN/Inf in foreign function call (arg 1)
    ## inside this call to alfak2
  fit <- alfak2(x,min_obs = ni)
  xfq <- lapply(fit$xo, function(xi) xi[xi$id=="fq",])
  
  fq_ids <- rbind(cbind(1,1:nrow(xfq[[1]])),cbind(2,1:nrow(xfq[[2]])))
  
  xv_res <- do.call(rbind,lapply(1:nrow(fq_ids), function(ix){
    ## the try catch is necessary because sometimes we get:
    ## Error in if (tr <= 0) break : missing value where TRUE/FALSE needed
    ## inside the call to Krig... 
    tryCatch({optim_loo2(fq_ids[ix,1],fq_ids[ix,2],x,xfq)},
             error=function(e) return(NULL))
  }))
  
  fit$min_obs <- ni
  fit$xv_res <- xv_res
  fit},error=function(e) return(NULL))
})



saveRDS(fits,"data/salehi/alfak_fits_htert_combined.Rds")






