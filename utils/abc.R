## this file generates and runs a sweep of ABM simulations with user specifed
## parameters. 
library("optparse")

option_list = list(
  make_option(c("-n", "--name"), type="character", default=NULL, 
              help="output directory.", metavar="character"),
  make_option(c("-c", "--cores"), type="numeric", default=1, 
              help="number of parallel cores [default= %default]"),
  make_option(c("-f", "--filename"), type="character", default=NULL, 
              help="path to config file")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$name)){
  print_help(opt_parser)
  stop("At least two argument must be supplied (name,filename)", call.=FALSE)
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

script_dir <- paste0(root.dir,"utils/")
cpp_source <- paste0(root.dir,"ABM/bin/ABM")
if(Sys.info()["sysname"]=="Windows") cpp_source <- paste0(cpp_source,".exe")
sweep_dir <- paste0(root.dir,opt$name)
dir.create(sweep_dir,recursive = T)
source(paste0(script_dir,"sim_setup_functions.R"))
source(paste0(script_dir,"abc_functions.R"))




library(parallel)

##
config <- readRDS(paste0(root.dir,opt$filename))
config$data_sources <- paste0(root.dir,config$data_sources)
fit <- readRDS(paste0(root.dir,config$landscape))$fit

cl <- makeCluster(getOption("cl.cores", min(opt$reps,opt$cores)))
clusterCall(cl, function(script_dir){
  source(paste0(script_dir,"abc_functions.R"))
  source(paste0(script_dir,"landscape_functions.R"))
  source(paste0(script_dir,"ALFA-K.R"))
},script_dir=script_dir)

i <- 1:config$npars
ids <- stringr::str_pad(i,5,pad=0)
x <- parLapplyLB(cl=cl,ids,getLL,root.dir=root.dir,sweep_dir=sweep_dir,cpp_source=cpp_source,fit=fit,config=config)
