## generate a random field landscape
gen_rf_landscape <- function(founder,Nwaves,scalef=NULL,wavelength=1){
  
  if(is.null(scalef)) scalef <- 1/(pi*sqrt(Nwaves))
  pk <- lapply(1:Nwaves, function(i){
    #pk <- sample(0:10,length(founder),replace=T)
    pk <- sample((-10):20,length(founder),replace=T)
  })
  
  d <- sapply(pk,function(ci) {
    sqrt(sum((founder-ci)^2))
  })

  return(pk)
}

##get fitness from random field landscape
get_rf_fitness <- function(k,pk,scalef=NULL,wavelength=1){
  if(is.null(scalef)) scalef <- 1/(pi*sqrt(length(pk)))
  d <- sapply(pk,function(ci) {
    sqrt(sum((k-ci)^2))
  })
  f0=sum(sin(d/wavelength)*scalef)
  return(f0)
}

get_ruggedness <- function(p,pk,scalef=NULL,wavelength=1){
  Nchrom <- length(p)
  n <- expand.grid(lapply(1:Nchrom, function(novar) c(-1,1)))
  n <- n[!apply(n,1,function(ni) sum(abs(ni))==0),]
  n <- t(apply(n,1,function(ni) ni+p))
  fn <- apply(n,1,get_rf_fitness,pk=pk,scalef=scalef,wavelength=wavelength)
  fp <- get_rf_fitness(p,pk,scalef=scalef,wavelength=wavelength)
  mean(abs(fp-fn))
}

get_complexity <- function(v,pk,scalef=NULL,wavelength=1){
  f0 <- get_rf_fitness(v,pk,scalef=scalef,wavelength=wavelength)
  mean(sapply(1:length(v), function(i){
    vi <- v
    vi[i] <- vi[i]+1
    ff <- get_rf_fitness(vi,pk,scalef=scalef,wavelength=wavelength)
    
    vi <- v
    vi[i] <- vi[i]-1
    fr <- get_rf_fitness(vi,pk,scalef=scalef,wavelength=wavelength)
    
    sd(c(ff-f0,f0-fr))
  }))
}

assess_complexity <- function(wavelength,Nwaves,Nchrom=4,Nreps= 10){
  print(c(wavelength,Nwaves,Nchrom))
  x <- sapply(1:Nreps, function(novar){
    pk <- gen_rf_landscape(founder=rep(2,Nchrom),Nwaves = Nwaves,
                           wavelength=wavelength)
    p <- rep(2,Nchrom)
    get_complexity(p,pk,wavelength=wavelength)
  })
  
  data.frame(wavelength=wavelength,Nwaves=Nwaves,Nchrom=Nchrom,
             mean.complexity=mean(x),sd.complexity=sd(x))
}

assess_ruggedness <- function(wavelength,Nwaves,Nchrom=4,Nreps= 10){
  print(c(wavelength,Nwaves,Nchrom))
  x <- sapply(1:Nreps, function(novar){
    pk <- gen_rf_landscape(founder=rep(2,Nchrom),Nwaves = Nwaves,
                           wavelength=wavelength)
    p <- rep(2,Nchrom)
    get_ruggedness(p,pk,wavelength=wavelength)
  })
  
  data.frame(wavelength=wavelength,Nwaves=Nwaves,Nchrom=Nchrom,
             mean.ruggedness=mean(x),sd.ruggedness=sd(x))
}