setwd("~/projects/008_birthrateLandscape/ALFA-K/")
source("utils/ALFA-K.R")

min_obs <- 5

ff <- list.files("data/salehi/alfak_inputs/")
ff <- ff[c(4,6,8,9)]

lapply(ff, function(fi){
  print(fi)
  tryCatch({
    x <- readRDS(paste0("data/salehi/alfak_inputs/",fi))
    
    n_obs <- table(rowSums(x$x))
    n <- as.numeric(names(n_obs))
    n <- n[n>=min_obs]
    n <- n[order(n,decreasing=T)]
    n <- n[-c(1:2)]
    print(n)
    fits <- lapply(n,function(ni){
      tryCatch({
      fit <- alfak(x,min_obs = ni)
      xfq <- fit$xo[fit$xo$id=="fq",]
      fv <- unlist(sapply(1:nrow(xfq),optim_loo,xx=x,xo=xfq))
      xfq$f_xv <- fv
      fit$min_obs <- ni
      fit$xv_res <- xfq
      fit$vx_cor <- tryCatch({cor(xfq$f_est,xfq$f_xv,use="complete")},
                             error=function(e) return(-Inf))
      print(fit$vx_cor)
      fit
      },error=function(e) return(NULL))
      
    })
    names(fits) <- paste0("minobs_",n)
    saveRDS(fits,paste0("data/salehi/alfak_fits_minobs_adaptive/",fi))
  },error=function(e) print(e))
  
})



