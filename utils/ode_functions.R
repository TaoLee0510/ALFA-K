library(Matrix)
library(igraph)

tmbuild <- function(p0,coords,dims){
  xx <- coords[,"ff"]*coords[,"nc"]*p0#/2
  sources <- coords[,"ii"]==coords[,"jj"] 
  xx[sources] <-  coords[sources,"ff"]*(1-2*coords[sources,"nc"]*p0)
  N <- max(coords[,"jj"])
  sparseMatrix(i = coords[,"ii"], 
               j = coords[,"jj"], 
               x = xx,
               dims=dims)
}
gen_coords <- function(x){
  ##x is an alfa-k fitted object
  ## function generates a list containing an ordered list of karyotypes (k)
  ## and a list of non- zero entries (neighbours) in a Sparse Transition matrix
  ## of karyotypes (coords)
  
  u <- rownames(x$xo)
  
  v <- rbind(do.call(rbind,lapply(u,s2v)),
             gen_all_neighbours(u))
  
  
  f <- c(predict(x$fit,v))
  v <- apply(v,1,paste,collapse=".")
  indices <- unlist(lapply(1:22,rep,2))
  
  coords <- do.call(rbind,pbapply::pblapply(1:length(v),function(i){
    v0 <- s2v(v[i])
    
    nn <- apply(gen_all_neighbours(v[i]),1,paste,collapse=".")
    nc <- v0[indices]
    names(nc) <- nn
    ii <- which(v%in%nn)
    nc <- as.numeric(nc[v[ii]])
    
    ii <- c(i,ii)
    nc <- c(sum(v0),nc)
    
    jj <- rep(i,length(ii))
    ff <- rep(f[i],length(ii))
    cbind(ii,jj,ff,nc)
  }))
  
  return(list(coords=coords,kary=v))
}
screenR <- function(fi,p=c(0.001,0.005)){
  y <- readRDS(paste0("figures/misseg_landscape_exploration/coords/",fi))
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


