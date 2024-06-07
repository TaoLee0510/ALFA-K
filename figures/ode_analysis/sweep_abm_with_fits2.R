#setwd("~/projects/ALFA-K")
setwd("~/projects/008_birthrateLandscape/ALFA-K")

ff <- c("SA535_CISPLATIN_CombinedH_X7_l_3_d1_0_d2_0.Rds",
        "SA535_CISPLATIN_CombinedH_X9_l_5_d1_0_d2_0.Rds",
        "SA906a_X57_l_7_d1_0_d2_0.Rds")


#library(parallel)
#cl <- makeCluster(getOption("cl.cores", 3))
lapply(ff, function(fi){
  
  Nreps <- 3
  Nsteps <- 10000
  pop_write_freq <- 100
  pgd <- 0
  max_size <- 1000000
  gen_abm_landscape <- function(fit){
    knots <- fit$knots
    cc <- fit$c
    d <- fit$d
    fscape <- rbind(cbind(knots,cc),c(d))
    return(fscape)
  }
  source("utils/sim_setup_functions.R")
  source("utils/ALFA-K.R")
  options(scipen=999)
  dir <- "data/salehi/alfak_fits/minobs_5/"
  out_dir <- "data/salehi/misseg_landscape_exploration/minobs_5/"
  dir.create(out_dir,recursive = T)
  
  cfig_template <- readLines("data/salehi/config_template.txt")  
  
  conditions <- c(0.0001)
  
  lapply(conditions, function(pmis){
    xi <- readRDS(paste0(dir,fi))
    
    cpp_cmd <- "ABM/bin/ABM.exe"
    
    od <- gsub("[.]","p",pmis)
    
    id <- unlist(strsplit(fi,split=".Rds"))[1]
    dir.create(paste0(out_dir,id,"/",od),recursive = T)
    abm_out_dir <- paste0(out_dir,id,"/",od,"/output/")
    dir.create(abm_out_dir)
    
    x <- readRDS(paste0("data/salehi/alfak_inputs_v2/",fi))
    dt <- x$dt
    x0 <- x$x[,1]
    x0 <- x0[x0>0]
    pop0 <- cbind(do.call(rbind,lapply(names(x0),s2v)),x0)
    pop0 <- do.call(rbind,lapply(readRDS(paste0("figures/ode_analysis/screen_winners/",fi)),s2v))
    pop0 <- cbind(pop0,10)
    rownames(pop0) <- NULL
    colnames(pop0) <- NULL
    
    pop0[,ncol(pop0)] <- round(100000*pop0[,ncol(pop0)]/sum(pop0[,ncol(pop0)]))
    
    lscape <- gen_abm_landscape(xi$fit)
    
    cfig_path <- paste0(out_dir,id,"/",od,"/cfig.txt")
    lscape_path <- paste0(out_dir,id,"/",od,"/lscape.txt")
    pop_path <- paste0(out_dir,id,"/",od,"/pop.txt")
    
    cfig <- modify_config("fitness_landscape_file",lscape_path,cfig_template)
    cfig <- modify_config("output_dir",abm_out_dir,cfig)
    cfig <- modify_config("pop_write_freq",pop_write_freq,cfig)
    cfig <- modify_config("Nsteps",Nsteps,cfig)
    cfig <- modify_config("p",pmis,cfig)
    cfig <- c(cfig,paste0("population_file,",pop_path))
    cfig <- c(cfig,paste0("pgd,",pgd))
    cfig <- c(cfig,paste0("max_size,",max_size))
    
    writeLines(cfig,cfig_path)
    write.table(pop0,pop_path,sep=",",quote = F,col.names = F,row.names = F)
    write.table(lscape,lscape_path,sep=",",quote = F,col.names = F,row.names = F)
    
    cmd <- paste(cpp_cmd,cfig_path)
    
    for(i in 1:Nreps) (system(cmd))
  })
  
})





