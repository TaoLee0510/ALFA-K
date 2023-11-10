
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
  to_modify <- which(grepl(parname,config))
  config[to_modify] <- new_entry 
  return(config)
}


setup_abm <- function(pars,pop=NULL){
  options(scipen=999)
  config <- c(init_kary="2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2",
              fitness_landscape_type="krig",
              fitness_landscape_file="landscape.txt",
              dt=0.1,
              pop_write_freq=1000000000,
              max_size=2000000,
              p=0.00005,
              pgd=0.0001,
              Nsteps=3000,
              output_dir="out",
              init_size=100000,
              population_file="custom_pop.txt"
  )
  if(!"population_file"%in%names(pars)) config <- config[!names(config)=="population_file"]
  dir.create(pars["output_dir"],recursive=T)
  if(!is.null(pop)){
    pop <- apply(pop,1,paste,collapse=",")
    writeLines(pop,pars["population_file"])
  }
  if("fitness_landscape_file"%in%names(pars)){
    knots <- fit$knots
    cc <- fit$c
    d <- fit$d
    fscape <- rbind(cbind(knots,cc),c(d))
    write.table(fscape, pars["fitness_landscape_file"],row.names = FALSE,col.names=FALSE,sep=",")
  }
  
  
  config[names(pars)] <- pars
  config <- paste0(names(config),",",config)
  setwd(pars["output_dir"])
  setwd("../")
  config_path <- paste0(getwd(),"/config.txt")
  writeLines(config,config_path)
  return(config_path)
    
}



