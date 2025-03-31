run_fits_in_subdir <- function(ci,dir="data/main/",
                               target_dir="alfak_eq_space_fits",
                               output_dir="alfak_eq_space_preds"){
  setwd("~/projects/ALFA-K")
  cpp_cmd <- "ABM/bin/ABM"
  source("utils/sim_setup_functions.R")
  source("utils/ALFA-K.R")
  init_pop_size <- 100000
  path2fits <- paste(dir,ci,target_dir,sep="/")
  path2preds <- paste(dir,ci,output_dir,sep="/")
  
  cfig_dir <- paste(path2preds,"config",sep="/")
  lscape_dir <- paste(path2preds,"landscape",sep="/")
  pop_path <- paste(path2preds,"pop.txt",sep="/")
  
  
  dir.create(cfig_dir,recursive=T)
  dir.create(lscape_dir,recursive=T)
  
  fits <- list.files(path2fits)
  fits
  
  fi <- fits[1]
  x <- readRDS(paste(path2fits,fi,sep="/"))
  
  
  
  pop0 <- x$input$x
  pop0 <- cbind(do.call(rbind,lapply(rownames(pop0),s2v)),pop0[,ncol(pop0),drop=F])
  rownames(pop0) <- NULL
  colnames(pop0) <- NULL
  pop0[,ncol(pop0)] <- round(init_pop_size*pop0[,ncol(pop0)]/sum(pop0[,ncol(pop0)]))
  write.table(pop0,pop_path,sep=",",row.names = F,col.names = F)
  
  cfig0 <- readLines(paste(dir,ci,"config.txt",sep="/"))
  for(fi in fits){
    fitname <- unlist(strsplit(fi,split=".Rds"))[1]
    odir <- paste(path2preds,"abm_output",fitname,sep="/")
    dir.create(odir,recursive = T)
    x <- readRDS(paste(path2fits,fi,sep="/"))
    pop.fitness <- x$input$pop.fitness
    fpred <- predict(x$fit$fit,pop0[,-ncol(pop0)])
    deltaf <- tail(pop.fitness,1)-mean(fpred)
    lscape <- gen_abm_landscape(x$fit$fit,deltaf)
    lscape_path <- paste0(lscape_dir,"/",fitname,".csv")
    write.table(lscape, lscape_path,row.names = FALSE,col.names=FALSE,sep=",")
    cfig <- modify_config("fitness_landscape_file",parval = lscape_path,config = cfig0)
    cfig <- modify_config("output_dir",parval = odir,config = cfig)
    cfig <- modify_config("fitness_landscape_type","krig",cfig)
    cfig <- c(cfig,paste0("population_file,",paste(path2preds,"pop.txt",sep="/")))
    cfig <- modify_config("Nsteps",1000,cfig)
    cfig_path <- paste0(cfig_dir,"/",fitname,".txt")
    writeLines(cfig,cfig_path)
    
    cmd <- paste(cpp_cmd,cfig_path)
    for(i in 1:3) system(cmd)
  }
}

#library(parallel)
#setwd("~/projects/ALFA-K")
#cl <- makeCluster(70)
#subdirs <- list.files("data/main")
#x <- parLapplyLB(cl,subdirs,run_fits_in_subdir)
#print("done")


 

