stop("this file is deprecated!")
library(lhs)
##prepare ABM data for input to pipeline
proc_sim <- function(dir,times){
  dt <- read.csv(paste0(dir,"/log.txt"),sep=",")$dt
  
  filenames <- list.files(dir)
  filenames <- filenames[!filenames%in%c("proc_data","log.txt")]
  tx <- as.numeric(substr(filenames,1,5))
  filenames <- paste0(dir,"/",sapply(times, function(ti) filenames[which.min(abs(ti-tx))]))
  tx <- sapply(times, function(ti) tx[which.min(abs(ti-tx))])
  
  x <- lapply(filenames,function(fi){
    xi <- read.csv(fi,header=F)
    nchrom <- ncol(xi)-2
    k <- apply(xi[,1:nchrom],1,paste,collapse=".")
    n <- xi[,nchrom+1]
    fitness <- xi[,nchrom+2]
    data.frame(karyotype=k,n=n,fitness=fitness)
  })
  
  
  
  pop.fitness <- sapply(x,function(xi) sum(xi$n*xi$fitness)/sum(xi$n))
  tmp <- do.call(rbind,x)  
  tmp <- tmp[!duplicated(tmp$karyotype),]
  unique.karyotypes <- tmp$karyotype
  clone.fitness <- tmp$fitness
  
  x <- do.call(rbind,lapply(as.character(unique.karyotypes), function(k){
    sapply(x, function(xi) sum(xi$n[xi$karyotype==k]))
  }))

  
  colnames(x) <- tx
  rownames(x) <- unique.karyotypes  
  list(x=x,pop.fitness=pop.fitness,clone.fitness = clone.fitness,dt=dt)
}

## get fitness from GRF
getf <- function(pk,fobj){
  d <- apply(fobj$peaks,1,function(pki) sqrt(sum((pk-pki)^2)))
  sum(sin(d/fobj$wavelength)*fobj$scale)
}

## read a GRF landscape for use to calculate fitness
gen_fitness_object <- function(cfig_path,lscape_path){
  lscape <- read.table(lscape_path,sep=",")
  cfig <- readLines(cfig_path)
  pk_scale <- cfig[grepl("scale",cfig)]
  pk_wl <- cfig[grepl("wavelength",cfig)]
  pk_scale <- as.numeric(tail(unlist(strsplit(pk_scale,split=",")),1))
  pk_wl <- as.numeric(tail(unlist(strsplit(pk_wl,split=",")),1))
  list(peaks=lscape,scale=pk_scale,wavelength=pk_wl)
}


adj_pop_fitness <- function(f,u,n,pop.fitness) (pop.fitness - u/n*f)/(1-u/n)


optimx_v0 <- function(f,u,fit_dat){
  
  ux <- u/colSums(fit_dat$x)
  up <- 1-ux
  #ux <- pmax(0.001,pmin(ux,0.999))
  #up <- pmax(0.001,pmin(up,0.999))
  ux0 <- ux
  up0 <- up
  fp <- (fit_dat$pop.fitness-ux*f)/up#pmin(3,pmax(0,(fit_dat$pop.fitness-ux*f)/up))
  
  tx <- as.numeric(colnames(fit_dat$x))
  
  ## we would like to capture any fitness implications of the first observable timepoint
  i <- min(which(u>0))
  if(FALSE & i>1 & length(which(u>0))==1){
    ## assume that a clone emerged in the middle of the observation period"
    t <- fit_dat$dt*diff(tx)[i-1]/2
    ux[i] <- ux0[i]/1000*exp(t*f)
    up[i] <- max(0.0001,up0[i-1])*exp(t*mean(c(fp[i],fp[i-1])))
    sx <- ux[i]+up[i]
    ux[i]<- ux[i]/sx
    up[i]<- up[i]/sx
  }
  
  for(i in (1+min(which(u>0))):length(ux)){
    t <- fit_dat$dt*diff(tx)[i-1]
    ux[i] <- ux0[i-1]*exp(t*f)
    up[i] <- max(0.0001,up0[i-1])*exp(t*mean(c(fp[i],fp[i-1])))
    sx <- ux[i]+up[i]
    ux[i]<- ux[i]/sx
    up[i]<- up[i]/sx
  }
  #print(ux)
  negll <- -log(dbinom(u,colSums(fit_dat$x),prob=ux))
  negll[!is.finite(negll)] <- 10^9
  sum(negll)
}

opt_tp1_v0 <- function(ftp1,n1,x1,f1,ff,t){
  ux0 <- x1/sum(c(x1,n1))
  up0 <- 1-ux0
  
  fp <- (ff-ux0*ftp1)/up0
  
  ux <- ux0*exp(t*ftp1)
  up <- pmax(0.0001,up0)*exp(t*fp)
  sx <- ux+up
  ux <- ux/sx
  
  negll <- -log(dbinom(0,sum(c(n1,x1)),prob=ux))
  f_all <- c(ftp1,f1)
  negll[!is.finite(negll)] <- 10^9
  negd <- -log(dnorm(ftp1,mean(f_all),max(0.001,sd(f_all))))
  sum(negll)+sum(negd)
}



optimsimple_v0 <- function(x,min_obs=5){
  ntp <- apply(x$x,1,function(i) sum(i>0))
  tp1 <- x$x[,1]>0 & ntp==1
  n <- rowSums(x$x)
  to_fit <- which(n>min_obs & !tp1)
  xtp1 <- x$x[n>min_obs & tp1,,drop=FALSE]
  x$x <- x$x[to_fit,,drop=FALSE]
  
  
  if(length(to_fit)==1){
    x_opt <- data.frame(f_est=mean(x$pop.fitness),err=0,n=sum(x$x),ntp=ncol(x$x))
    rownames(x_opt) <- rownames(x$x)
    return(x_opt)
  }
  ix <- 1:nrow(x$x)
  
  peakest <- apply(x$x,1,which.max)
  init_guess <- x$pop.fitness[peakest]
  opt <- do.call(rbind,lapply(ix, function(i){
    u <- x$x[i,]
    fi <- u/colSums(x$x)
    fmax <- max(fi)
    opt <- NULL
    if(fmax>0.95){
      ## if there is only one clone present the optimx function breaks,
      ## so it is better to just use the pop growth rate to determine
      ## fitness of that clone
      opt <- list(minimum=x$pop.fitness[which.max(fi)],objective=0)
    }else{
      igi <- init_guess[i]
      val <- Inf
      ntries <- 1
      
      while(val>1e4 & ntries < 10){
        ntries <- ntries + 1
        interval <- c(igi-1/ntries,igi+1/ntries)
        opt <- optimise(optimx_v0,interval=interval,u=u,fit_dat=x)
        val <- opt$objective
      }
      if(opt$objective>1e6) opt$minimum <- igi
    }
    
    data.frame(f_est=opt$minimum,u0=NaN,ll=opt$objective,
               n=sum(u),ntp=sum(u>0))
    
  }))
  rownames(opt) <- rownames(x$x)
  opt
}


optimx <- function(par,u,fit_dat){
  f <- par[1]
  pf <- adj_pop_fitness(f,u,fit_dat$N, fit_dat$pop.fitness)
  xp <- c(0,cumsum(diff(as.numeric(names(u))*fit_dat$dt)*zoo::rollmean(pf,k=2)))
  xi <- c(0,cumsum(diff(as.numeric(names(u))*fit_dat$dt)*f))+par[2]
  ux <- exp(xi-xp)
  negll <- -log(dbinom(u,fit_dat$N,prob=ux/(1+ux)))
  negll[!is.finite(negll)] <- 10^9
  sum(negll)
}

clone_traj <- function(par,u,fit_dat){
  f <- par[1]
  pf <- adj_pop_fitness(f,u,fit_dat$N, fit_dat$pop.fitness)
  xp <- c(0,cumsum(diff(as.numeric(colnames(fit_dat$x))*fit_dat$dt)*zoo::rollmean(pf,k=2)))
  xi <- c(0,cumsum(diff(as.numeric(colnames(fit_dat$x))*fit_dat$dt)*f))+par[2]
  ux <- exp(xi-xp)
  ux/(1+ux)
}



clone_traj <- function(par,u,fit_dat,ntp){
  f <- par[1]
  pf <- adj_pop_fitness(f,u,fit_dat$N, fit_dat$pop.fitness)
  xp <- c(0,cumsum(diff(as.numeric(colnames(fit_dat$x))*fit_dat$dt)*zoo::rollmean(pf,k=2)))
  xi <- c(0,cumsum(diff(as.numeric(colnames(fit_dat$x))*fit_dat$dt)*f))+par[2]
  ux <- exp(xi-xp)
  data.frame(time=as.numeric(colnames(fit_dat$x))*fit_dat$dt,frequency=ux/(1+ux))
}
clone_traj <- function(par,u,fit_dat,ntp){
  
  tt <- as.numeric(names(u))*x$dt
  t <- seq(min(tt),max(tt),length.out=ntp)
  
  f <- par[1]
  pf <- adj_pop_fitness(f,u,fit_dat$N, fit_dat$pop.fitness)
  pf <- zoo::rollmean(pf,k=2)
  tmid <- zoo::rollmean(t,k=2)
  pf <- as.numeric(sapply(tmid,function(ti) pf[sum(ti>tt)]))
  xp <- c(0,cumsum(diff(t)*pf))
  xi <- c(0,cumsum(diff(t)*f))+par[2]
  ux <- exp(xi-xp)
  ux <- ux/(1+ux)
  
  data.frame(time=t,frequency=ux)
  
}

optimsimple <- function(x,min_obs=5,mintp=0){
  ntp <- apply(x$x,1,function(i) sum(i>0))
  n <- rowSums(x$x)
  N <- colSums(x$x)
  to_fit <- (n>min_obs) & ntp>mintp
  x$x <- x$x[to_fit,,drop=FALSE]
  
  fit_dat <- x
  fit_dat$N <- N
  
  opt <- do.call(rbind,lapply(1:nrow(x$x), function(i){
    u <- x$x[i,]
    fi <- u/fit_dat$N
    fmax <- max(fi)
    opt <- NULL
    
    ##initial guess of fitness (pop. fitness when clone freq. peaks)
    init_guess <- x$pop.fitness[which.max(u)]
    
    ##initial guess of starting frequency
    ## if clone is present at tp 1, its frequency there is our guess
    u0 <- log(fi[1])
    if(!is.finite(u0)){
      pktm <- which.max(u)
      xp <- c(0,cumsum(diff(as.numeric(names(u))*x$dt)*zoo::rollmean(x$pop.fitness,k=2)))[pktm]
      xi <- c(0,cumsum(diff(as.numeric(names(u))*x$dt)*init_guess))[pktm]
      u0 <- xp-xi
    }
    
    par <- c(init_guess,u0)
    
    lower <- c(init_guess-0.5,-Inf)
    upper <- c(init_guess+0.5,Inf)
    
    optres <- optim(par,optimx,method = "L-BFGS-B",
                    upper=upper,lower=lower,
                    u=u,fit_dat=fit_dat)
    
    data.frame(f_est=optres$par[1],u0=exp(optres$par[2]), ll = optres$value, n=sum(u),ntp=sum(u>0))
    
    
  }))
  rownames(opt) <- rownames(x$x)
  
  N <- colSums(x$x)
  ## add in values for the known bad cases
  f_est <- sapply(1:nrow(x$x), function(i){
    u <- x$x[i,]
    default_f_est <- x$pop.fitness[which.max(u/N)]
    if(max(u/N)>0.95) return(default_f_est)
    if(sum(u>0)==1) return(default_f_est)
    return(opt$f_est[i])
  })
  opt$f_est <- f_est
  opt
}
fitm <- function(par,xj,dt){
  N <- colSums(xj)
  x0i <- exp(head(par,length(par)/2))
  fi <- tail(par,length(par)/2)
  
  logx <- do.call(rbind,lapply(1:length(x0i), function(i){
    fi[i]*as.numeric(colnames(xj))*dt + log(x0i[i])
  }))
  pred <- apply(logx,2,function(lxi) exp(lxi)/sum(exp(lxi)))
  pred <- matrix(pred,ncol=ncol(xj))
  
  ll <- do.call(rbind,lapply(1:nrow(xj), function(i){
    -dbinom(xj[i,],size = N,prob=pred[i,],log=T)
  }))
  ll <- ll[,!N==0]
  ll[!is.finite(ll)] <- 1e+9
  sum(ll)
}

inifitm <- function(inipar, opar,xj,dt){
  N <- colSums(xj)
  x0i <- c(exp(head(opar,length(opar)/2)),exp(inipar[1]))
  fi <- c(tail(opar,length(opar)/2),inipar[2])
  
  logx <- do.call(rbind,lapply(1:length(x0i), function(i){
    fi[i]*as.numeric(colnames(xj))*dt + log(x0i[i])
  }))
  pred <- apply(logx,2,function(lxi) exp(lxi)/sum(exp(lxi)))
  
  ll <- do.call(rbind,lapply(1:nrow(xj), function(i){
    -dbinom(xj[i,],size = N,prob=pred[i,],log=T)
  }))
  ll <- ll[,!N==0]
  ll[!is.finite(ll)] <- 1e+9
  sum(ll)
  
}
opt_g_free <- function(x,min_obs=5,mintp=0){
  
  dt <- x$dt
  x <- x$x
  x <- x[rowSums(x)>min_obs,]
  ntp <- apply(x,1,function(xi) sum(xi>0))
  x <- x[ntp>mintp,]
  x <- x[order(rowSums(x),decreasing=T),]
  ##some initial parameter guesses
  x0 <- log(1/10^9)
  f0 <- 0.5
  
  i <- 1
  xj <- x[1:i,,drop=F]
  par <- c(x0,f0)
  
  opt <- optim(par,fitm,xj=xj,dt=dt)
  par <- opt$par
  err <- opt$value
  x0 <- head(par,length(par)/2)
  f0 <- tail(par,length(par)/2)
  
  for(i in 2:nrow(x)){
    print(i)
    guesses <- randomLHS(n=100,k=2)
    xj <- x[1:i,,drop=F]
    rx0 <- max(x0+25)-min(x0-25)
    guesses[,1] <- (guesses[,1]-0.5)*rx0+median(x0)
    
    inipar <- apply(guesses,1,function(gi){
      iniopt <- optim(gi,fn = inifitm,opar=par,xj=xj,dt=dt)
      c(iniopt$value,iniopt$par)
    })
    
    inipar <- inipar[2:3,which.min(inipar[1,])]
    par <- c( c(x0,inipar[1]),c(f0,inipar[2]))
    opt <- optim(par,fitm,xj=xj,dt=dt)
    par <- opt$par
    err <- c(err,opt$value)
    x0 <- head(par,length(par)/2)
    f0 <- tail(par,length(par)/2)
  }
  min_f0 <- min(f0)
  delta_f0 <- 0.2-min_f0
  f0 <- f0+delta_f0
  xopt <- data.frame(f_est=f0,u0=x0,ll=err,n=rowSums(x),ntp=apply(x,1,function(xi) sum(xi>0)))
  
  rownames(xopt) <- rownames(xj)[1:nrow(xopt)]
  xopt
}



plotm <- function(xo,tt=seq(0,300,10),nN,t0){
  
  nNt <- approx(x=t0,y=nN,xout=tt)$y
  
  N <- do.call(rbind,lapply(1:nrow(xo), function(i){
    exp(xo$u0[i])*exp(xo$f_est[i]*tt)
  }))
  
  N <- t(apply(N,2,function(Ni) Ni/sum(Ni)))
  #N <- apply(N,1, function(Ni) Ni*nNt)
  N <- data.frame(N)
  colnames(N) <- rownames(xo)
  N$time <- tt
  N <- reshape2::melt(N,id.vars="time")
  colnames(N) <- c("time","clone","frequency")
  N
}

optimise_frequent_clones <- function(x,min_obs=5,use_gdata=F,mintp=0){
  if(use_gdata) return(optimsimple(x,min_obs,mintp))
  else return(opt_g_free(x,min_obs,mintp))
}

plot_frequent_clones <- function(xo,x,use_gdata=F){
  if(use_gdata){
    x$N <- colSums(x$x)
    x$x <- x$x[rownames(xo),]
    
    
    
    df <- do.call(rbind,lapply(1:nrow(xo),function(i){
      pari <- c(xo$f_est[i],log(xo$u0[i]))
      dfi <- clone_traj(par=pari,u=x$x[i,],fit_dat = x,ntp=100)
      dfi$clone <- rownames(xo)[i]
      dfi
    }))
    
    df <- df[,c("time","clone","frequency")]
    return(df)
  }
  else{
    N <- colSums(x$x)
    n <- colSums(x$x[rownames(xo),])
    nN <- n/N
    t0 <- as.numeric(names(nN))*x$dt
    df <- plotm(xo,tt=(round(min(t0)):round(max(t0))),nN=nN,t0=t0)
    tt <- as.numeric(as.character(df$time))
    df$time <- tt
    return(df)
  }
}

plot_frequent_clone_data <- function(x,xo){
  y <- data.frame(x$x,check.names=F)
  y <- y[rownames(xo),]
  for(i in 1:ncol(y)) y[,i] <- y[,i]/sum(y[,i])
  y$id <- rownames(y)
  y <- reshape2::melt(y,id.vars='id')
  colnames(y) <- c("clone","tt","frequency")
  y$time <- as.numeric(as.character(y$tt))*x$dt
  return(y)
}

gen_all_neighbours <- function(ids,as.strings=T){
  if(as.strings) ids <- lapply(ids, function(ii) as.numeric(unlist(strsplit(ii,split="[.]"))))
  nkern <- do.call(rbind,lapply(1:length(ids[[1]]), function(i){
    x0 <- rep(0,length(ids[[1]]))
    x1 <- x0
    
    x0[i] <- -1
    x1[i] <- 1
    rbind(x0,x1)
  }))
  n <- do.call(rbind,lapply(ids, function(ii) t(apply(nkern,1,function(i) i+ii))))
  n <- unique(n)
  nids <- length(ids)
  n <- rbind(do.call(rbind,ids),n)
  n <- unique(n)
  n <- n[-(1:nids),]  
  n <- n[apply(n,1,function(ni) sum(ni<1)==0),]
  n
}

## this function generates unique neighboring clones of an input population, 
## until there is a ratio of fc between new/old pop size.
expand_pop <- function(x,fc=10,as.strings=T){
  if(as.strings) x <- do.call(rbind,lapply(x, function(ii) as.numeric(unlist(strsplit(ii,split="[.]")))))
  n0 <- nrow(x)
  n1 <- nrow(x)
  while((n1-n0)/n0<fc){
    x <- rbind(x,apply(x,2,function(xi) xi+sample(c(-1,0,1),size=nrow(x),replace = T,prob = c(0.05,0.9,0.05))))
    x <- unique(x)
    n1 <- nrow(x)
  }
  x[-(1:n0),]
}

find_grad_dist <- function(x_opt){
  xmat <- do.call(rbind,lapply(rownames(x_opt), function(xi){
    as.numeric(unlist(strsplit(xi,split="[.]")))
  }))
  
  dx <- as.matrix(dist(xmat))
  dy <- as.matrix(dist(x_opt$f_est))
  dy <- dy[which(dx==1&upper.tri(dx))]
  dy <- c(dy,-dy)
  
  sd(dy)
}
pij<-function(i, j, beta){
  qij <- 0
  if(abs(i-j)>i){ ## not enough copies for i->j
    return(qij)
  }
  # code fails for j = 0, but we can get this result by noting i->0 always
  ## accompanies i->2i
  if(j==0) j <- 2*i
  s <- seq(abs(i-j), i, by=2 )
  for(z in s){
    qij <- qij + choose(i,z) * beta^z*(1-beta)^(i-z) * 0.5^z*choose(z, (z+i-j)/2)
  }
  ## if there is a mis-segregation: ploidy conservation implies that both daughter cells will emerge simultaneously
  #if(i!=j) qij <- 2*qij
  
  return(qij)
}
optim_neighbor_fitness <- function(f,tp,iflux,ftp,sdy,u,tt,pop_size,f_est,dt){
  xtp <- rep(0,length(iflux))
  xtp[1] <- iflux[1]
  for(j in 2:length(xtp)){
    ## in the following statement we want c1 to evaluate to -Inf when 
    ## xtp[j-1] = 1. Due to precision we can have cases where xtp[j-1] > 1
    ## which evaluates to NaN. Fixed by bounding input to log at zero
    c1 <- log(max(1/xtp[j-1]-1,0)) 
    t <- diff(tp)[1]*dt
    xi <- exp(f*t)/(exp(f*t)+exp(c1+ftp[j]*t))
    xtp[j] <- iflux[j]+xi
  }
  ux <- approx(tp,xtp,xout=tt)$y
  ux <- pmin(ux,0.9999)
  deltaf <- abs(f-f_est)
  negll <- c(-log(dbinom(u,pop_size,prob=ux)), -log(dnorm(deltaf,mean=0,sd=sdy)))
  negll[!is.finite(negll)] <- 10^9## to prevent warnings when probability is zero
  sum(negll)
}
get_neighbor_fitness <- function(ni,x_opt,x,pm0,ntp=100){
  ni_s <- paste(ni,collapse=".")
  
  sdy <- find_grad_dist(x_opt)
  
  ##for each 1ms neighbour, identify likely parents
  sources <- rownames(x_opt)[rownames(x_opt)%in%apply(gen_all_neighbours(ni_s),1,paste,collapse=".")]
  x_opt_i <- x_opt[sources,,drop=FALSE]
  fu <- x$x[sources,,drop=FALSE]/colSums(x$x)
  
  ## get the longitudinal frequency data for neighbor of interest
  u <- rep(0,ncol(x$x))
  names(u) <- colnames(x$x)
  if(ni_s%in%rownames(x$x)) u <- x$x[ni_s,]
  
  ##misseg probability from parent to neighbor
  pm <- sapply(sources,function(si){
    si <- as.numeric(unlist(strsplit(si,split="[.]")))
    prod(sapply(1:length(si), function(k) pij(si[k],ni[k],pm0)))
  })
  
  ##interpolate fitness, parent frequency and time
  ## essentially we are approximating an ODE from here on
  mintp <- min(as.numeric(colnames(x$x)))
  maxtp <- max(as.numeric(colnames(x$x)))
  tp <- seq(mintp,maxtp,(maxtp-mintp)/ntp)
  
  
  ftp <- approx(as.numeric(colnames(x$x)),x$pop.fitness,xout = tp)$y
  
  iflux <- colSums(do.call(rbind,lapply(1:length(pm),function(j) {
    aj <- approx(as.numeric(colnames(x$x)),as.numeric(fu[j,]),xout = tp)$y
    aj*diff(tp)[1]*x$dt*pm[j]*x_opt_i$f_est[j]
  })))
  
  oni <- optimise(optim_neighbor_fitness,interval=c(0,1),tp=tp,iflux=iflux,ftp=ftp,sdy=sdy,
                  u=u,tt=colnames(x$x),pop_size=colSums(x$x),f_est=x_opt_i$f_est,dt=x$dt)
  data.frame(f_est=oni$minimum,u0=NaN,ll=oni$objective,n=sum(u),ntp=sum(u>0))
}

infer_population_fitness <- function(x,x_opt){
  m <- x$x[rownames(x_opt),] 
  for(i in 1:ncol(m)){
    m[,i]<- m[,i]*x_opt$f_est/sum(m[,i])
  }
  colSums(m)
}

wrap_neighbor_fitness <- function(x,x_opt,pm0=0.00005,ntp=100,use_gdata=T){
  if(!use_gdata) x$pop.fitness <- infer_population_fitness(x,x_opt)
  n <- gen_all_neighbours(rownames(x_opt))
  x2 <- do.call(rbind,lapply(1:nrow(n), function(i){
    get_neighbor_fitness(n[i,],x_opt,x,pm0,ntp)
  }))
  rownames(x2) <- apply(n,1,paste,collapse=".")
  x2
}

R2 <- function(obs,pred){
  1-sum((pred-obs)^2)/sum((obs-mean(obs))^2)
}
RMSE <- function(obs,pred){
  round(sqrt(mean((obs-pred)^2)),digits=2)
}
