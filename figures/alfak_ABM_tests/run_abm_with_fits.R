setwd("~/projects/ALFA-K")
#setwd("~/projects/008_birthrateLandscape/ALFA-K")

dir <- "data/main/"
conditions <- list.files(dir)
ids <- as.character(sapply(conditions, function(i) unlist(strsplit(i,split="_"))[6]))
conditions <- conditions[ids=="0.00005"]
dirs <- paste0(dir,conditions,"/")

run_sims_in_dir <- function(d0){
  gen_abm_landscape <- function(fit,deltaf=0){
    knots <- fit$knots
    cc <- fit$c
    d <- fit$d
    d[1] <- d[1]+deltaf
    fscape <- rbind(cbind(knots,cc),c(d))
    return(fscape)
  }
  tryCatch({
  cpp_cmd <- "ABM/bin/ABM.exe"
  cpp_cmd <- "ABM/bin/ABM"
  source("utils/sim_setup_functions.R")
  source("utils/ALFA-K.R")
  x <- proc_sim(paste0(d0,"train/00000"),times=2000)
  pop.fitness <- x$pop.fitness
  x <- x$x
  pop0 <- cbind(do.call(rbind,lapply(rownames(x),s2v)),x)
  rownames(pop0) <- NULL
  colnames(pop0) <- NULL
  pop0[,ncol(pop0)] <- round(100000*pop0[,ncol(pop0)]/sum(pop0[,ncol(pop0)]))
  write.table(pop0,paste0(d0,"pop2000.txt"),sep=",",row.names = F,col.names = F)
  cfig_dir <- paste0(d0,"test_cfig/")
  lscape_dir <- paste0(d0,"test_lscape/")
  dir.create(cfig_dir)
  dir.create(lscape_dir)
  
  cfig0 <- readLines(paste0(d0,"config.txt"))
  swp_od <- paste0(d0,"test/")
  fits <- list.files(paste0(d0,"sweep_fits"))
  
  for(fi in fits){
    fitname <- unlist(strsplit(fi,split=".Rds"))[1]
    odir <- paste0(d0,"test_v2/",fitname)
    dir.create(odir,recursive = T)
    x <- readRDS(paste0(d0,"sweep_fits/",fi))
    fpred <- predict(x$fit,pop0[,-ncol(pop0)])
    deltaf <- pop.fitness-mean(fpred)
    lscape <- gen_abm_landscape(x$fit,deltaf)
    lscape_path <- paste0(lscape_dir,fitname,".csv")
    write.table(lscape, lscape_path,row.names = FALSE,col.names=FALSE,sep=",")
    cfig <- modify_config("fitness_landscape_file",parval = lscape_path,config = cfig0)
    cfig <- modify_config("output_dir",parval = odir,config = cfig)
    cfig <- modify_config("fitness_landscape_type","krig",cfig)
    cfig <- c(cfig,paste0("population_file,",d0,"pop2000.txt"))
    cfig <- modify_config("Nsteps",1000,cfig)
    cfig_path <- paste0(cfig_dir,fitname,".txt")
    writeLines(cfig,cfig_path)
    
    
    
    cmd <- paste(cpp_cmd,cfig_path)
    for(i in 1:3) system(cmd)
  }
  },error=function(e) return(NULL))
}


#lapply(dirs,run_sims_in_dir)
#stop()

library(parallel)
cl <- makeCluster(getOption("cl.cores", 8))

x <- parLapplyLB(cl=cl,X=dirs, fun=run_sims_in_dir)