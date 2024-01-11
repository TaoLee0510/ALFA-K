
R2 <- function(obs,pred){
  1-sum((pred-obs)^2)/sum((obs-mean(obs))^2)
}
RMSE <- function(obs,pred){
  round(sqrt(mean((obs-pred)^2)),digits=2)
}

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

## gets wasserstein distance between two populations formatted in a specific way.
## better use wasserstein distance function as a wrapper
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
##formats objects for get_dwass
make_wass_object <- function(x,t,is.multiple.objects=F){
  k <- NULL
  if(!is.multiple.objects){
    k <- kary_vec(x,t)
  }else{
    k <- unlist(lapply(x,kary_vec,t=t))
    tmp <- aggregate(list(n=k),by=list(k=names(k)),sum)
    k <- tmp$n
    names(k) <- tmp$k
  }
  
  kx <- do.call(rbind,lapply(names(k),function(ki){
    as.numeric(unlist(strsplit(ki,split="[.]")))
  }))
  list(x=kx,n=as.numeric(k))
}

## gets wasserstein distance between two populations formatted as output from 
## proc_sim. Can use lists of populations. 
wasserstein_distance <- function(test,ref,t,is.test.multiple.objects=F,is.ref.multiple.objects=F){
  s1 <- make_wass_object(test,1200,is.multiple.objects = is.test.multiple.objects)
  s2 <- make_wass_object(ref,1200,is.multiple.objects = is.ref.multiple.objects)
  get_dwass(s1,s2)
}

##function takes output of proc_sim and returns karyotype frequency vector at a certain timepoint
kary_vec <- function(x,t){
  col_id <- which.min(abs(as.numeric(colnames(x$x))-t))
  k <- x$x[,col_id]
  k[k>0]
}

kary_df <- function(vt,vr){
  rx <- unique(c(names(vt),names(vr)))
  df <- data.frame(test=rep(0,length(rx)),ref=rep(0,length(rx)),row.names = rx)
  tmp <- aggregate(list(n=vr),by=list(k=names(vr)),sum)
  vr <- tmp$n
  names(vr) <- tmp$k
  df$test[rownames(df)%in%names(vt)] <- vt[rownames(df)[rownames(df)%in%names(vt)]]
  df$ref[rownames(df)%in%names(vr)] <- vr[rownames(df)[rownames(df)%in%names(vr)]]
  #df <- df[rowSums(df)>0,]
  return(df)
}

cna <- function(k){
  if(!is.numeric(k)) k <- as.numeric(unlist(strsplit(k,split="[.]")))
  cnas <- rep(0,length(k))
  cnas[which(k>median(k))] <- 1
  cnas[which(k<median(k))] <- -1
  # or (?) paste0(which(k<median(k)),"-")
  cnas
}

ll_cna <- function(test,ref,t){
  vt <- kary_vec(test,t)
  vr <- unlist(lapply(ref, kary_vec,t=t))
  df <- kary_df(vt,vr)
  cnas <- do.call(rbind,lapply(rownames(df),cna))
  cngs <- apply(cnas,2,function(ci) sum(as.numeric(ci>0)*df$ref))
  cnns <- apply(cnas,2,function(ci) sum(as.numeric(ci==0)*df$ref))
  cnls <-  apply(cnas,2,function(ci) sum(as.numeric(ci<0)*df$ref))
  nobs <- sum(df$ref)
  cngs <- cngs/nobs
  cnls <- cnls/nobs
  cnns <- cnns/nobs
  ll <- sapply(1:ncol(cnas),function(i){
    sum(-log(cngs[i])*df$test[cnas[,i]>0])+sum(-log(cnls[i])*df$test[cnas[,i]<0])+sum(-log(cnns[i])*df$test[cnas[,i]==0])
  })
  sum(ll)
}

##functions for computing angle metric
mag <- function(v) sqrt(sum(v^2))
getangle <- function(a,b) 180*acos(sum(a*b)/(mag(a)*mag(b)))/pi
get_mean <- function(x,t,is.multiple.objects=F){
  k <- NULL
  if(!is.multiple.objects){
    k <- kary_vec(x,t)
  }else{
    k <- unlist(lapply(x,kary_vec,t=t))
    tmp <- aggregate(list(n=k),by=list(k=names(k)),sum)
    k <- tmp$n
    names(k) <- tmp$k
  }
  n <- as.numeric(k)
  N <- sum(n)
  xi <- do.call(rbind,lapply(names(k),function(ki){
    as.numeric(unlist(strsplit(ki,split="[.]")))
  }))
  as.numeric(apply(xi,2,function(k){
    sum(k*n)/N
  }))
}
angle_metric <- function(test,ref,t=1200,is.test.multiple.objects=F,is.ref.multiple.objects=F,x0=2){
  xt <- get_mean(test,t,is.multiple.objects = is.test.multiple.objects)-x0
  xr <- get_mean(ref,t,is.multiple.objects = is.ref.multiple.objects)-x0
  getangle(xt,xr)
}
