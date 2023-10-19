
process_dist <- function(dir,target.fitness){
  ksim <- lapply(dir, function(di){
    x <- proc_sim(di,times=seq(0,3000,100))
    #plot(x$pop.fitness)
    dx <- x$x[,which.min(abs(x$pop.fitness-target.fitness))]
    dx <- dx[order(dx,decreasing=T)]
    return(dx)
  })
  
  ksim <- unlist(ksim)
  ksim <- ksim[ksim>0]
  dsim <- aggregate(list(n=ksim),by=list(k=names(ksim)),sum)
  tmp <- dsim$n
  names(tmp) <- dsim$k
  return(tmp)
}

load_dat <- function(data_sources){
  kdat <- lapply(data_sources, function(di){
    dd <- readRDS(di)
    dd$x[,ncol(dd$x)] ## always takes the last, better have as a parameter
  })
  kd <- unlist(kdat)
  kd <- aggregate(list(n=kd),by=list(k=names(kd)),sum)
  ddat <- kd$n
  names(ddat) <- kd$k
  return(ddat)
}

ll_cna <- function(df){
  cnas <- do.call(rbind,lapply(rownames(df),cna))
  cngs <- apply(cnas,2,function(ci) sum(as.numeric(ci>0)*df$sim))
  cnns <- apply(cnas,2,function(ci) sum(as.numeric(ci==0)*df$sim))
  cnls <-  apply(cnas,2,function(ci) sum(as.numeric(ci<0)*df$sim))
  nobs <- sum(df$sim)
  cngs <- cngs/nobs
  cnls <- cnls/nobs
  cnns <- cnns/nobs
  ll <- sapply(1:ncol(cnas),function(i){
    sum(-log(cngs[i])*df$dat[cnas[,i]>0])+sum(-log(cnls[i])*df$dat[cnas[,i]<0])+sum(-log(cnns[i])*df$dat[cnas[,i]==0])
  })
  sum(ll)
}



dist_cna <- function(df){
  cnas <- do.call(rbind,lapply(rownames(df),cna))
  fneut <- sum(df$dat*apply(cnas,1,function(i) sum(i==0)))/sum(df$dat*nchrom)
  cngs <- apply(cnas,2,function(ci) sum(as.numeric(ci>0)*df$sim))
  cnns <- apply(cnas,2,function(ci) sum(as.numeric(ci==0)*df$sim))
  cnls <-  apply(cnas,2,function(ci) sum(as.numeric(ci<0)*df$sim))
  nobs <- sum(df$sim)
  cngs <- cngs/nobs
  cnls <- cnls/nobs
  cnns <- cnns/nobs
  
  cngd <- apply(cnas,2,function(ci) sum(as.numeric(ci>0)*df$dat))
  cnnd <- apply(cnas,2,function(ci) sum(as.numeric(ci==0)*df$dat))
  cnld <-  apply(cnas,2,function(ci) sum(as.numeric(ci<0)*df$dat))
  nobs <- sum(df$dat)
  cngd <- cngd/nobs
  cnld <- cnld/nobs
  cnnd <- cnnd/nobs
  
  rbind(data.frame(chrom=1:22,fgain=cngs,floss=cnls,fneut=cnns,id="sim"),
        data.frame(chrom=1:22,fgain=cngd,floss=cnld,fneut=cnnd,id="data"),
        data.frame(chrom=1:22,fgain=0.5*(1-fneut),floss=0.5*(1-fneut),fneut=fneut,id="null"))
}


getLL <- function(id,root.dir,sweep_dir,cpp_source,fit,config){
  
  this.dir <- paste0(sweep_dir,"/",id)
  
  # generate random loguniform distributed parameters
  pars <- exp(runif(2,log(c(config$pgd_range[1],config$p_range[1])),
                    log(c(config$pgd_range[2],config$p_range[2]))))
  pgd <- pars[1]
  p <- pars[2]
  
  # run ABM
  cpp_cmd <- abm_from_krig(fit,dir = this.dir,pars=c(Nsteps=config$Nsteps,dt=config$dt,pgd=pgd,p=p),
                           cpp_cmd = cpp_source)
  for(i in 1:config$nreps) tmp <- system(cpp_cmd,intern = T)
  
  ##compute ll on summary stat
  target.dir <- paste0(this.dir,"/train/")
  target.dir <- paste0(target.dir,list.files(target.dir))
  dsim <- process_dist(target.dir,config$target.fitness)
  ddat <- load_dat(config$data_sources) ##slightly wasteful to load each time
  rx <- unique(c(names(ddat),names(dsim)))
  df <- data.frame(sim=rep(0,length(rx)),dat=rep(0,length(rx)))
  df$sim[rx%in%names(nsim)] <- nsim[rx[rx%in%names(nsim)]]
  df$dat[rx%in%names(ndat)] <- ndat[rx[rx%in%names(ndat)]]
  rownames(df) <- rx
  logl <- ll_cna(df)
  dcna <- dist_cna(df)
  result <- c(p,pgd,logl)
  names(result) <- c("p","pgd","ll")
  
  ##save result
  saveRDS(result,paste0(this.dir,"/result.Rds"))
  saveRDS(dcna,paste0(this.dir,"/cna_dist.Rds"))
  #delete sim output (takes a lot of space)
  unlink(target.dir,recursive = T)
}

