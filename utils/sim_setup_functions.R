
gen_randscape <- function(founder,Nwaves,scalef,wavelength=1){
  if(is.null(scalef)) scalef <- 1/(pi*sqrt(Nwaves))
  f0 <- 0
  ## we want to have simulations where the diploid founder is fit enough to survive.
  ## if scalef <- 1/(pi*sqrt(30))
  ## then this would make the diploid cell in the top 10% fittest
  ## clones:
  while(f0<=0.4){
    pk <- lapply(1:Nwaves, function(i){
      pk <- sample((-10):20,length(founder),replace=T)
    })
    
    d <- sapply(pk,function(ci) {
      sqrt(sum((founder-ci)^2))
    })
    f0=sum(sin(d/wavelength)*scalef)
    
  }
  sapply(pk, paste,collapse="," )
}


gen_config <- function(output_dir="output",init_kary=rep(2,10),fitness_landscape_type="gaussian",
                       fitness_landscape_file="landscapes/landscape.txt",dt=0.1,
                       p=0.00005,Nsteps=3000,init_size=100000,scalef=1,wavelength=1){
  
  c(paste(c("init_kary", init_kary),collapse=","),
    paste(c("fitness_landscape_type",fitness_landscape_type),collapse=","),
    paste(c("fitness_landscape_file", fitness_landscape_file),collapse=","),
    paste(c("dt",dt),collapse=","),
    paste(c("p", p),collapse=","),
    paste(c("Nsteps", Nsteps),collapse=","),
    paste(c("output_dir", output_dir),collapse=","),
    paste(c("init_size", init_size),collapse=","),
    paste(c("scale", scalef),collapse=","),
    paste(c("wavelength", wavelength),collapse=","))
}

gen_replicate <- function(Nchrom,wavelength,sweep_dir,cpp_source,p=0.00005,Nwaves=10){
  options(scipen=999)
  ##setup all the filepaths and create necessary directories:
  sweep_id <- paste("N",Nchrom,"w",paste0(floor(wavelength),"p",round(10*(wavelength-floor(wavelength)),digits=1)),"m",p,sep="_")
  rep_id <- sum(grepl(sweep_id,list.files(sweep_dir)))
  sweep_id <- paste(sweep_id,"rep",stringr::str_pad(rep_id,2,pad=0),sep="_")
  
  cpp_output_dir <- paste0(sweep_dir,"/",sweep_id,"/train")
  cpp_landscape_path <- paste0(sweep_dir,"/",sweep_id,"/landscape.txt")
  cpp_config_path <- paste0(sweep_dir,"/",sweep_id,"/config.txt")
  
  dir.create(paste0(sweep_dir,"/",sweep_id))
  dir.create(cpp_output_dir)
  dir.create(paste0(sweep_dir,"/",sweep_id,"/test"))
  
  cpp_run_cmd <- paste(cpp_source,cpp_config_path)
  
  ## generate landscape
  founder <- rep(2,Nchrom)
  scalef <- 1/(pi*sqrt(Nwaves))
  lscape <- gen_randscape(founder,Nwaves=10,scalef=scalef,wavelength=wavelength)
  writeLines(lscape,cpp_landscape_path)
  ##generate config
  cfig <- gen_config(fitness_landscape_type = "random",
                     init_kary = rep(2,Nchrom),
                     fitness_landscape_file = cpp_landscape_path,
                     output_dir = cpp_output_dir,scalef=scalef,p=p,
                     wavelength=wavelength)
  writeLines(cfig,cpp_config_path)

  return(list(cpp_run_cmd=cpp_run_cmd, sweep_id=sweep_id))
}

modify_config <- function(parname,parval,config){
  new_entry <- paste(parname,parval,sep=",")
  to_modify <- which(grepl(paste0(parname,","),config))
  config[to_modify] <- new_entry 
  return(config)
}


setup_abm <- function(sim_dir,pars=NULL,pop=NULL,fit=NULL){
  options(scipen=999)
  
  config <- c(init_kary="2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2",
              fitness_landscape_type="krig",
              fitness_landscape_file=paste0(sim_dir,"/landscape.txt"),
              dt=0.1,
              pop_write_freq=1000000000,
              max_size=2000000,
              p=0.00005,
              pgd=0.0001,
              Nsteps=3000,
              output_dir=paste0(sim_dir,"/out"),
              init_size=100000,
              population_file=paste0(sim_dir,"/pop.txt")
  )
  dir.create(config["output_dir"],recursive=T)
  if(is.null(pop)){
    config <- config[!names(config)=="population_file"]
  }else{
    pop <- apply(pop,1,paste,collapse=",")
    writeLines(pop,config["population_file"])
  }
  if(is.null(fit)){
    stop("I didn't add this yet")
  }else{
    knots <- fit$knots
    cc <- fit$c
    d <- fit$d
    fscape <- rbind(cbind(knots,cc),c(d))
    write.table(fscape, config["fitness_landscape_file"],row.names = FALSE,col.names=FALSE,sep=",")
  }
  
  config[names(pars)] <- pars
  config <- paste0(names(config),",",config)
  config_path <- paste0(sim_dir,"/config.txt")
  writeLines(config,config_path)
  return(config_path)
}

gen_abm_landscape <- function(fit,deltaf=0){
  knots <- fit$knots
  cc <- fit$c
  d <- fit$d
  d[1] <- d[1]+deltaf
  fscape <- rbind(cbind(knots,cc),c(d))
  return(fscape)
}

##given an existing sweep, run new replicates of the sims in the directory 
## with neutral evolution
run_sims_in_dir_neutral <- function(d0,tstart=2000,Nsteps=1000,Nreps=10){
  tryCatch({
    cpp_cmd <- "ABM/bin/ABM"
    source("utils/ALFA-K.R")
    x <- proc_sim(paste0(d0,"train/00000"),times=tstart)
    x <- x$x
    pop0 <- cbind(do.call(rbind,lapply(rownames(x),s2v)),x)
    rownames(pop0) <- NULL
    colnames(pop0) <- NULL
    pop0[,ncol(pop0)] <- round(100000*pop0[,ncol(pop0)]/sum(pop0[,ncol(pop0)]))
    write.table(pop0,paste0(d0,"pop2000.txt"),sep=",",row.names = F,col.names = F)
    
    cfig0 <- readLines(paste0(d0,"config.txt"))
    odir <- paste0(d0,"flats/")
    dir.create(odir,recursive = T)
    cfig <- modify_config("fitness_landscape_type",parval = "flat",config = cfig0)
    cfig <- cfig[!grepl("fitness_landscape_file",cfig)]
    cfig <- modify_config("output_dir",parval = odir,config = cfig)
    cfig <- c(cfig,paste0("population_file,",d0,"pop2000.txt"))
    cfig <- modify_config("Nsteps",Nsteps,cfig)
    cfig_path <- paste0(d0,"flat_config.txt")
    writeLines(cfig,cfig_path)
    
    
    
    cmd <- paste(cpp_cmd,cfig_path)
    for(i in 1:Nreps) system(cmd)
  },error=function(e) return(NULL))
}

## a wrapper function for run_sims_in_dir_neutral
generate_neutral_reference <- function(root.dir="~/projects/ALFA-K",
                                       target.dir="data/main/",
                                       tstart=2000,Nsteps=1000,Nreps=10,
                                       ncores=50){
  setwd(root.dir)
  conditions <- list.files(target.dir)
  ids <- as.character(sapply(conditions, function(i) unlist(strsplit(i,split="_"))[6]))
  dirs <- paste0(target.dir,conditions,"/")
  library(parallel)
  cl <- makeCluster(getOption("cl.cores", min(ncores,length(dirs))))
  clusterExport(cl = cl,c("modify_config"),envir=environment())
  
  x <- parLapplyLB(cl=cl,X=dirs, fun=run_sims_in_dir_neutral,
                   tstart=tstart,Nsteps=Nsteps,Nreps=Nreps)
}

## Take the sampled output of a population (evolved on a GRF landscape),
## at a particular timepoint, then rerun on the same landscape (test contigency
## and possible bottlenecking caused by limited info & noise)
fork_sims_in_dir <- function(d0,tstart=2000,Nsteps=1000,Nreps=10,noise=0){
  tryCatch({
    cpp_cmd <- "ABM/bin/ABM"
    source("utils/ALFA-K.R")
    x <- proc_sim(paste0(d0,"train/00000"),times=tstart)
    x <- x$x
    pop0 <- cbind(do.call(rbind,lapply(rownames(x),s2v)),x)
    rownames(pop0) <- NULL
    colnames(pop0) <- NULL
    pop0[,ncol(pop0)] <- round(100000*pop0[,ncol(pop0)]/sum(pop0[,ncol(pop0)]))
    write.table(pop0,paste0(d0,"pop2000.txt"),sep=",",row.names = F,col.names = F)
    
    cfig0 <- readLines(paste0(d0,"config.txt"))
    cfig0 <- cfig0[!grepl("fitness_noise",cfig0)]
    cfig0 <- c(cfig0,paste0("fitness_noise,",noise))
    odir <- paste0(d0,"forks/n_",gsub("[.]","p",noise),"/")
    dir.create(odir,recursive = T)
    cfig <- modify_config("output_dir",parval = odir,config = cfig0)
    cfig <- c(cfig,paste0("population_file,",d0,"pop2000.txt"))
    cfig <- modify_config("Nsteps",Nsteps,cfig)
    cfig_path <- paste0(d0,"forks_cfig/")
    dir.create(cfig_path,recursive = T,showWarnings = F)
    cfig_path <- paste0(cfig_path,"flat_config.txt")
    writeLines(cfig,cfig_path)
    
    
    
    cmd <- paste(cpp_cmd,cfig_path)
    for(i in 1:Nreps) system(cmd)
  },error=function(e) return(NULL))
}

## a wrapper function for run_sims_in_dir_neutral
generate_fork_reference <- function(root.dir="~/projects/ALFA-K",
                                       target.dir="data/main/",
                                       tstart=2000,Nsteps=1000,Nreps=10,
                                        noise = c(0,0.1,0.01),
                                       ncores=50){
  setwd(root.dir)
  conditions <- list.files(target.dir)
  ids <- as.character(sapply(conditions, function(i) unlist(strsplit(i,split="_"))[6]))
  dirs <- paste0(target.dir,conditions,"/")
  library(parallel)
  cl <- makeCluster(getOption("cl.cores", min(ncores,length(dirs))))
  clusterExport(cl = cl,c("modify_config"),envir=environment())
  
  for(ni in noise){
    x <- parLapplyLB(cl=cl,X=dirs, fun=fork_sims_in_dir,
                     tstart=tstart,Nsteps=Nsteps,Nreps=Nreps,noise=ni)
  }
  
  
}