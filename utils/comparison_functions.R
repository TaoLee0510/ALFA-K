
load_sim_dat <- function(fo,g){
 foi <- list.files(fo) 
 x <- foi[!foi%in%c("log.txt","proc_data")]
 
 x <- unlist(strsplit(x,split=".csv"))
 x <- as.numeric(x)
 x <- foi[which.min(abs(g-x))]
 x <- read.csv(paste0(fo,"/",x),header = F)
 nchrom <- ncol(x)-2
 list(x=x[,1:nchrom],n=x[,1+nchrom],fitness=x[,2+nchrom])
}

summarise_sim_dat <- function(x,repr=NaN,id=NaN,g=NaN,landscape=NaN){
  df <- data.frame(n=x$n,fitness=x$fitness,rep=repr,id=id,g=g,landscape=landscape)
  rownames(df) <- apply(x$x,1,paste,collapse=".")
  df
}

wrap_sim_summary <- function(folders,ids,landscapes,g){
  xi <- lapply(1:length(folders), function(i){
    foi <- paste0(folders[[i]],list.files(folders[[i]]))
    do.call(rbind,lapply(1:length(foi), function(j){
      xij <- load_sim_dat(foi[j],g=g)
      xij <- summarise_sim_dat(xij,repr=j,id=ids[i],g=g,landscape=landscapes[i])
    }))
  })
  do.call(rbind,xi)
}

make_tsne <- function(x){
  xtsne <- Rtsne(x)
  xtsne <- xtsne$Y
  xtsne <- data.frame(xtsne)
  rownames(xtsne) <- apply(x,1,paste,collapse=".")
  xtsne
}

wrap_tsne <- function(folders,ids,landscapes,g=1500){
  xi <- lapply(1:length(folders), function(i){
    foi <- paste0(folders[[i]],list.files(folders[[i]]))
    lapply(foi,load_sim_dat,g=g)
  })
  x <- do.call(rbind,lapply(xi, function(xij){
    do.call(rbind,lapply(xij,function(xijk) xijk$x))
  }))
  x <- unique(x)
  xtsne <- make_tsne(x) 
  xtsne
}

ssd <- function(f1,f2,g=g,landscape=NaN){
  x1 <- do.call(rbind,lapply(f1, function(f1i){
    xi <- load_sim_dat(f1i,g=g)
    data.frame(kary = apply(xi$x,1,paste,collapse="."),
               n=xi$n)
  }))
  x1 <- aggregate(list(n=x1$n),by=list(kary=x1$kary),sum)
  
  n1 <- x1$n
  x1 <- do.call(rbind,lapply(x1$kary, function(xi){
    as.numeric(unlist(strsplit(xi,split="[.]")))
  }))  
  
  
  x2 <- do.call(rbind,lapply(f2, function(f2i){
    xi <- load_sim_dat(f2i,g=g)
    data.frame(kary = apply(xi$x,1,paste,collapse="."),
               n=xi$n)
  }))
  x2 <- aggregate(list(n=x2$n),by=list(kary=x2$kary),sum)
  
  n2 <- x2$n
  x2 <- do.call(rbind,lapply(x2$kary, function(xi){
    as.numeric(unlist(strsplit(xi,split="[.]")))
  }))  
  
  d <- as.matrix(dist(rbind(x1,x2)))^2
  n <- c(n1,n2)
  d <- t(t(d)*n)*n ##
  
  id1 <- 1:length(n1)
  id2 <- length(n1)+1:length(n2)
  
  d11 <- sum(d[id1,id1])/(2*sum(n1))
  d22 <- sum(d[id2,id2])/(2*sum(n2))
  d12 <- sum(d[id1,id2])/(sum(n1)+sum(n2))
  
  data.frame(d11=d11,d22=d22,d12=d12,g=g,landscape=landscape)
}

get_dwass <- function(s1,s2){
  #wpp fails if there is only one coordinate. 
  ## to avoid this, add a tiny amount of mass close to the single coordinate
  if(nrow(s1$x)<2){
    s1$x <- rbind(s1$x,0.001+s1$x)
    s1$n <- c(s1$n,0.01)
  }
  if(nrow(s2$x)<2){
    s2$x <- rbind(s2$x,0.001+s2$x)
    s2$n <- c(s2$n,0.01)
  }
  pp1 <- transport::wpp(s1$x,s1$n/sum(s1$n))
  pp2 <- transport::wpp(s2$x,s2$n/sum(s2$n))
  transport::wasserstein(pp1,pp2)
}