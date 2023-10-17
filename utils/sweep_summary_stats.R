## this file compiles and saves summary statistics from ABM sweep.
## First run ->generate_abm_sweep->process_abm_sweep -> this file
library("optparse")

option_list = list(
  make_option(c("-n", "--name"), type="character", default=NULL, 
              help="sweep folder name. If sweep already exists, new runs are added.", metavar="character"),
  make_option(c("-c", "--cores"), type="numeric", default=1, 
              help="number of parallel cores [default= %default]"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="output filename", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$name)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (I/O names)", call.=FALSE)
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

library(parallel)


sweep.dir <- paste0(opt$name,"/")

proc_run <- function(fi,sweep.dir){
  R2 <- function(obs,pred){
    1-sum((pred-obs)^2)/sum((obs-mean(obs))^2)
  }
  optimR2 <- function(offset,f_tru,f_est) R2(f_tru,f_est+offset)
  dir <- paste0(sweep.dir,fi,"/fits/00000/")
  if(!file.exists(paste0(dir,"krig.Rds"))) return(NULL)
  info <- unlist(strsplit(fi,split="_"))
  df <- data.frame(t(info[(1:length(info))%%2==0]))
  colnames(df) <- info[(1+1:length(info))%%2==0]
  xo <- readRDS(paste0(dir,"opt_data.Rds"))
  fit <- readRDS(paste0(dir,"krig.Rds"))
  cfig_path <- paste0(sweep.dir,fi,"/config.txt")
  lscape_path <- paste0(sweep.dir,fi,"/landscape.txt")
  fobj <- gen_fitness_object(cfig_path, lscape_path)
  
  n <- gen_all_neighbours(rownames(xo))
  rn <- apply(n,1,paste,collapse=".")
  n <- n[!rn%in%rownames(xo),]
  rn <- rn[!rn%in%rownames(xo)]
  xo <- rbind(xo,data.frame(f_est = predict(fit,n),u0=NaN,ll=NaN,n=0,ntp=0,id="d2n",
                   f_tru=apply(n,1,getf,fobj=fobj)))
  xo <- split(xo,f=xo$id)
  
  df$r_fq <- optimise(optimR2,interval=c(-1,1),f_tru=xo[["fq"]]$f_tru,
                      f_est=xo[["fq"]]$f_est,maximum=T)$objective
  df$r_nn <- optimise(optimR2,interval=c(-1,1),f_tru=xo[["nn"]]$f_tru,
                      f_est=xo[["nn"]]$f_est,maximum=T)$objective
  df$r_d2n <- optimise(optimR2,interval=c(-1,1),f_tru=xo[["d2n"]]$f_tru,
                       f_est=xo[["d2n"]]$f_est,maximum=T)$objective
  df$ll_fq <- mean(xo[["fq"]]$ll)
  df$ll_nn <- mean(xo[["nn"]]$ll)
  df$nfq <- nrow(xo$fq)
  df$nnn <- nrow(xo$nn)
  df$nd2n <- nrow(xo$d2n)
  saveRDS(df,paste0(dir,"summary_stats.Rds"))
  return(df)
}

cl <- makeCluster(getOption("cl.cores", opt$cores))
clusterCall(cl, function(script.dir){
  source(paste0(script.dir,"ALFA-K.R"))
  library(fields)
},script.dir=script.dir)

fo <- list.files(sweep.dir)
df <- parLapplyLB(cl=cl,fo,proc_run,sweep.dir=sweep.dir)
df <- do.call(rbind,df)
saveRDS(df,opt$output)