
setup_and_run <- function(fi){
  fit <- readRDS(paste0(fitDir,fi))$fit
  x <- readRDS(paste0(dataDir,fi))
  
  id <- gsub(".Rds","",fi)
  
  tarDir <- paste0(simDir,id,"/")
  outputDir <- paste0(tarDir,"output/")
  dir.create(outputDir,recursive = T)
  configPath <- paste0(tarDir,"config.txt")
  popPath <- paste0(tarDir,"pop.txt")
  lscapePath <- paste0(tarDir,"lscape.txt")
  
  pop <- x$x[,ncol(x$x)]
  pop <- pop[pop>0]
  pop <- round(pop*100000/sum(pop))
  pop <- gsub("[.]",",",paste(names(pop),pop,sep="."))
  writeLines(pop,popPath)
  
  lscape <- gen_abm_landscape(fit)
  write.table(lscape, lscapePath,row.names = FALSE,col.names=FALSE,sep=",")
  
  config <- modify_config("population_file",popPath,base_config)
  config <- modify_config("fitness_landscape_file",lscapePath,config)
  config <- modify_config("output_dir",outputDir,config)
  writeLines(config,configPath)
  
  cmd <- paste(cpp_cmd,cfig_path)
  for(i in 1:nRuns) system(cmd)
}

forwardSims <- function(alfak_dir="~/projects/ALFA-K",
                        fitDir = "data/salehi/alfak_fits/",
                        dataDir = "data/salehi/alfak_inputs/",
                        simDir = "data/salehi/forward_sims/",
                        nCores=60,nRuns=3){
  setwd(alfak_dir)
  ff <- list.files(fitDir)
  base_config <- c("init_kary,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2",
                   "fitness_landscape_type,krig",
                   "dt,0.1",
                   "p,0.00005",
                   "Nsteps,1000",
                   "init_size,100000",
                   "pop_write_freq,10",
                   "population_file,tbd",
                   "output_dir,tbd",
                   "fitness_landscape_file,tbd")
  cpp_cmd <- "ABM/bin/ABM"
  library(parallel)
  cl <- makeCluster(getOption("cl.cores", nCores))
  clusterExport(cl,varlist=c("cpp_cmd","base_config","simDir","dataDir","fitDir","nRuns"))
  clusterCall(cl,fun=function() source("utils/sim_setup_functions.R"))
  x <- parLapplyLB(cl=cl,X=ff,fun=setup_and_run)
}





