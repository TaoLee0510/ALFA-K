## this file generates and runs a sweep of ABM simulations with user specifed
## parameters. 
library("optparse")

option_list = list(
  make_option(c("-n", "--name"), type="character", default=NULL, 
              help="sweep folder name. If sweep already exists, new runs are added.", metavar="character"),
  make_option(c("-c", "--cores"), type="numeric", default=1, 
              help="number of parallel cores [default= %default]"),
  make_option(c("-m", "--min_obs"), type="numeric", default=5, 
              help="detection threshold for frequent clones [default= %default]")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$name)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (name)", call.=FALSE)
}

## https://stackoverflow.com/questions/1815606/determine-path-of-the-executing-script
thisFile <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    # Rscript
    return(normalizePath(sub(needle, "", cmdArgs[match])))
  } else {
    # 'source'd via R console
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
}

root.dir <- gsub("[\\]","/",thisFile())
root.dir <- unlist(strsplit(root.dir,split="/"))
root.dir <- root.dir[1:(length(root.dir)-2)]
root.dir <- paste0(paste(root.dir,collapse="/"),"/")

script.dir <- paste0(root.dir,"utils/")

## reads data from all sims in a sweep and creates an output folder for it in processed output
setwd(root.dir)
library(fields)
library(parallel)
optim_sweep <- function(fi,sweep.dir,min_obs){
  target.dir <- paste0(sweep.dir,"/",fi,"/")
  p <- as.numeric(unlist(strsplit(fi,split="_"))[6])
  save.dir <- paste0(target.dir,"fits")
  dir.create(save.dir)
  cfig_path <- paste0(target.dir,"config.txt")
  lscape_path <- paste0(target.dir,"landscape.txt")
  fobj <- gen_fitness_object(cfig_path, lscape_path)
  
  folders <- list.files(paste0(target.dir,"train"))
  
  for(fo in folders){
    x <- proc_sim(paste0(target.dir,"train/",fo),times=seq(400,2800,400))
    fit <- alfak(x,min_obs,misseg_rate = p)
    xo <- fit$xo
    fit <- fit$fit
    pks <- lapply(rownames(xo), function(pki) as.numeric(unlist(strsplit(pki,split="[.]"))))
    xo$f_tru <- sapply(pks,getf,fobj=fobj)
    save.subdir <- paste0(save.dir,"/",fo,"/")
    dir.create(save.subdir)
    saveRDS(x,paste0(save.subdir,"proc_data.Rds"))
    saveRDS(xo,paste0(save.subdir,"opt_data.Rds"))
    saveRDS(fit,paste0(save.subdir,"krig.Rds"))
  }
}

ff <- list.files(opt$name)

cl <- makeCluster(getOption("cl.cores", opt$cores))
clusterCall(cl, function(script.dir){
  source(paste0(script.dir,"ALFA-K.R"))
  library(fields)
},script.dir=script.dir)



x <- parLapplyLB(cl=cl,ff,optim_sweep,sweep.dir=opt$name,min_obs=opt$min_obs)