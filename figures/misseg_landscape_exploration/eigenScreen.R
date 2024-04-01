setwd("~/projects/008_birthrateLandscape/ALFA-K/")
source("utils/ALFA-K.R")
library(Matrix)
library(igraph)
library(ggplot2)
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

screenR <- function(fi){
  y <- readRDS(paste0("figures/misseg_landscape_exploration/coords/",fi))
  
  dims <- rep(length(y$kary),2)
  
  p <- c(0.0001,0.001,0.003,0.005)
  res <- do.call(cbind,lapply(p,function(p0){
    tm <- tmbuild(p0,y$coords,dims)
    func <- function(x, extra=NULL) { as.vector(tm %*% x) } 
    res <- tryCatch(expr = {
      res <- arpack(func, options=list(n=dims[1], nev=1, ncv=3, which="LM",maxiter=5000), 
                    sym=FALSE, complex = FALSE)
      res <- abs(res$vectors)
    },error=function(e) return(rep(0,nrow(tm))))
    
    names(res) <- y$kary
    res
  }))
  
  res <- res[order(rowSums(res),decreasing=T),]
  z <- res[apply(res,1,max)>0.02,]
  n0 <- z[,1]
  n1 <- z[,2]
  n2 <- z[,3]
  n3 <- z[,4]
  
  fit <- readRDS(paste0("data/salehi/alfak_fits/minobs_5/",fi))$fit
  pts <- do.call(rbind,lapply(rownames(z),s2v))
  f <- predict(fit,pts)
  ump <- umap::umap(pts)
  
  df <- data.frame(n0,n1,n2,n3,f)
  df <- cbind(ump$layout,df)
  colnames(df)[1:2] <- c("u1","u2")
  
  winlo <- rownames(df)[which.max(n0)]
  winhi <- rownames(df)[which.max(n3)]
  
  winners <- list(winhi=winhi,winlo=winlo)
  saveRDS(winners,paste0("figures/misseg_landscape_exploration/screen_winners/",fi))
  
  z1 <- reshape2::melt(df,measure.vars=c("n0","n1","n2","n3"))
  
  p1 <- ggplot(z1,aes(x=u1,y=u2,color=value))+
    facet_wrap(~variable,nrow=2)+
    geom_point()+
    scale_color_viridis_c(trans="log")
  
  
  p2 <- ggplot(df,aes(x=u1,y=u2,color=f))+
    geom_point()+
    scale_color_viridis_c()
  
  
  p <- cowplot::plot_grid(p1,p2,nrow=2)
  imname <- head(unlist(strsplit(fi,split=".Rds")),1)
  ggsave(paste0("figures/misseg_landscape_exploration/screen_ims/",imname,".png"),plot=p,width=6,height=9,units="in")
}

ff <- list.files("figures/misseg_landscape_exploration/coords/")
lapply(ff,screenR)
