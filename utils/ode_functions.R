library(Matrix)
library(igraph)

tmbuild <- function(p0,coords,dims){
  xx <- coords[,"ff"]*coords[,"nc"]*p0/2
  sources <- coords[,"ii"]==coords[,"jj"] 
  xx[sources] <-  coords[sources,"ff"]*(1-coords[sources,"nc"]*p0)
  N <- max(coords[,"jj"])
  sparseMatrix(i = coords[,"ii"], 
               j = coords[,"jj"], 
               x = xx,
               dims=dims)
}
screenR <- function(fi,p=c(0.001,0.005)){
  y <- readRDS(paste0("figures/ode_analysis/coords/",fi))
  dims <- rep(length(y$kary),2)
  
  res <- do.call(cbind,lapply(p,function(p0){
    tm <- tmbuild(p0,y$coords,dims)
    
    func <- function(x, extra=NULL) { as.vector(tm %*% x) } 
    res <- arpack(func, options=list(n=dims[1], nev=1, ncv=3, which="LM",maxiter=5000), 
                  sym=FALSE, complex = FALSE)
    res <- abs(res$vectors)
    res <- res/sum(res)
    names(res) <- y$kary
    res
  }))
  
  res <- res[order(rowSums(res),decreasing=T),]
  z <- res[apply(res,1,max)>0.001,]
  n0 <- z[,1]
  ni <- z[,2]
  n1 <- z[,ncol(z)]
  id <- n0>n1
  
  fit <- readRDS(paste0("data/salehi/alfak_fits/minobs_5/",fi))$fit
  pts <- do.call(rbind,lapply(rownames(z),s2v))
  f <- predict(fit,pts)
  
  df <- data.frame(n0,ni,n1,id,f)
  
  return(df)
}


