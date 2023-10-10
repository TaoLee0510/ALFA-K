## reads data from all sims in a sweep and creates an output folder for it in processed output
setwd("~/projects/008_birthrateLandscape/karyotype_evolution/")
library(fields)
library(parallel)
optim_sweep <- function(fi,sweep.dir,save.dir){
  p <- as.numeric(unlist(strsplit(fi,split="_"))[6])
  save.dir <- paste0(save.dir,fi,"/")
  dir.create(save.dir)
  cfig_path <- paste0(sweep.dir,fi,"/config.txt")
  lscape_path <- paste0(sweep.dir,fi,"/landscape.txt")
  fobj <- gen_fitness_object(cfig_path, lscape_path)
  
  x <- proc_sim(paste0(sweep.dir,fi,"/train/00000/"),times=seq(400,2800,400))
  xo <- optimise_frequent_clones(x,min_obs = 20)
  pks <- lapply(rownames(xo), function(pki) as.numeric(unlist(strsplit(pki,split="[.]"))))
  xo$f_tru <- sapply(pks,getf,fobj=fobj)
  xn <- wrap_neighbor_fitness(x,xo,use_gdata = F,pm0 = p)
  pks <- lapply(rownames(xn), function(pki) as.numeric(unlist(strsplit(pki,split="[.]"))))
  xn$f_tru <- sapply(pks,getf,fobj=fobj)
  
  xo$id <- "fq"
  xn$id <- "nn"
  xo <- rbind(xo,xn)
  
  xmat <- do.call(rbind,lapply(rownames(xo), function(i){
    as.numeric(unlist(strsplit(i,split="[.]")))
  }))
  y <- xo$f_est
  fit <- Krig(xmat,y,m=1,give.warnings = F)
  n <- gen_all_neighbours(rownames(xo))
  rn <- apply(n,1,paste,collapse=".")
  n <- n[!rn%in%rownames(xo),]
  rn <- rn[!rn%in%rownames(xo)]
  df <- data.frame(f_est = predict(fit,n),u0=NaN,ll=NaN,n=0,ntp=0,
                   f_tru=apply(n,1,getf,fobj=fobj),id="d2n")
  xo <- rbind(xo,df)
  saveRDS(xo,paste0(save.dir,"xo.Rds"))
  saveRDS(fit,paste0(save.dir,"fit.Rds"))
  list(fit=fit,xo=xo)
}

#print("enter path of sweep to be processed:")
#sweep.dir = readLines(con = "stdin", n = 1)

#print("enter path to save output:")
#save.dir = as.numeric(readLines(con = "stdin", n = 1))
sweep.dir <- "ABM/output/misrate/"
save.dir <- "proc_data/misrate/"
dir.create(save.dir)
script.dir="rscripts/"
ff <- list.files(sweep.dir)

cl <- makeCluster(getOption("cl.cores", 3))
clusterCall(cl, function(script.dir){
  source(paste0(script.dir,"analysis_functions_v2.R"))
  library(fields)
},script.dir=script.dir)



x <- parLapplyLB(cl=cl,ff,optim_sweep,sweep.dir=sweep.dir,save.dir=save.dir)