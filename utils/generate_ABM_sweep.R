## this file generates and runs a sweep of ABM simulations with user specifed
## parameters. 
library("optparse")

option_list = list(
  make_option(c("-n", "--name"), type="character", default=NULL, 
              help="sweep folder name. If sweep already exists, new runs are added.", metavar="character"),
  make_option(c("-c", "--cores"), type="numeric", default=1, 
              help="number of parallel cores [default= %default]"),
  make_option(c("-l", "--nchrom"), type="numeric", default=22, 
              help="number of chromosomes per cell [default= %default]"),
  make_option(c("-a", "--augment"), type="numeric", default=0, 
              help="run additional replicates of a sweep using existing landscapes"),
  make_option(c("-r", "--reps"), type="numeric", default=1, 
              help="number of ABM repeats per parameter combination [default= %default]"),
  make_option(c("-s", "--landscapereps"), type="numeric", default=1, 
              help="number of ABM repeats per idential landscape [default= %default]"),
  make_option(c("-w", "--wavelengths"), type="character", default="0.1,0.2,0.4,0.8,1.6", 
              help="wavelengths for GRF landscape [default= %default]"),
  make_option(c("-m", "--misrates"), type="character", default="0.001,0.0001,0.00001", 
              help="missegregation rates [default= %default]")
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

script_dir <- paste0(root.dir,"utils/")
cpp_source <- paste0(root.dir,"ABM/bin/ABM")
if(Sys.info()["sysname"]=="Windows") cpp_source <- paste0(cpp_source,".exe")
sweep_dir <- paste0(root.dir,opt$name)

library(parallel)
cl <- makeCluster(getOption("cl.cores", min(opt$reps,opt$cores)))
dir_success <- dir.create(sweep_dir,recursive = T)
if(opt$augment>0){
  if(dir_success) stop("parameter a intended for use with existing sweep")
  fo <- list.files(sweep_dir)
  cpp_cmds <- paste(cpp_source,paste0(sweep_dir,"/",fo,"/config.txt"))
  for(i in 1:opt$augment){
    parSapplyLB(cl,cpp_cmds,function(xx) system(xx))
  }
}

if(opt$augment==0){
  source(paste0(script_dir,"sim_setup_functions.R"))
  
  wavelengths <- as.numeric(unlist(strsplit(opt$wavelengths,split=",")))
  misrates <- as.numeric(unlist(strsplit(opt$misrates,split=",")))
  
  
  for(wavelength in wavelengths){
    for(misrate in misrates){
      setup_info <- lapply(1:opt$reps, function(dummyvar){
        gen_replicate(Nchrom=opt$nchrom,wavelength = wavelength,p = misrate,sweep_dir = sweep_dir,cpp_source = cpp_source)
      })
      for(i in 1:opt$landscapereps){
        parSapplyLB(cl,setup_info,function(xx) system(xx$cpp_run_cmd))
      }
      #print(setup_info)
    }
  }
}

