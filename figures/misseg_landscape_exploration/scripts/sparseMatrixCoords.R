
longest_cons <- function(ff,longest=T){
  ids <- sapply(ff, function(fi){
    fsplit <- unlist(strsplit(fi,split="_"))
    lsplit <- length(fsplit)
    fsplit <- fsplit[1:(lsplit-6)]
    paste(fsplit,collapse="_")
  })
  ff <- split(ff,f=ids)
  
  ff <- sapply(ff, function(fi){
    fi <- fi[order(fi)]
    if(longest) return(tail(fi,1))
    return(head(fi,1))
  })
  as.character(ff)
}

gen_coords <- function(fi){
  x <- readRDS(paste0("data/salehi/alfak_fits/minobs_5/",fi))
  
  u <- rownames(x$xo)
  
  ## ALL KARYOTYPES
  v <- rbind(do.call(rbind,lapply(u,s2v)),
             gen_all_neighbours(u))
  
  
  f <- c(predict(x$fit,v)) #KARYOTYPE FITNESS
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

setwd("~/projects/008_birthrateLandscape/ALFA-K/")
source("utils/ALFA-K.R")

x <- readRDS("figures/salehi_data_fitting/data/fit_summaries.Rds")
xl <- x[!x$has_descendents&!x$has_parents&x$min_obs==5&x$r2>0.3,]
xl <- xl[xl$filenames%in%longest_cons(xl$filenames),]

lapply(xl$filenames,function(fi){
  coords <- gen_coords(fi)
  saveRDS(coords,paste0("figures/misseg_landscape_exploration/data/coords/",fi))
})








