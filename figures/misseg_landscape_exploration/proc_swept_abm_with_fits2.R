setwd("~/projects/ALFA-K")
source("utils/ALFA-K.R")
#setwd("~/projects/008_birthrateLandscape/ALFA-K")

#library(parallel)
#cl <- makeCluster(getOption("cl.cores", 8))

proc_abm_ss <- function(fi,df){
  tt <- seq(0,10000,200)
  pks <- df[[fi]]
  mr <- list.files(paste0(dir,fi))
  if(length(mr)==0) return(NULL)
  z <- do.call(rbind,pbapply::pblapply(mr,function(mri){
    reps <- list.files(paste0(dir,fi,"/",mri,"/output/"))
    # print(reps)
    do.call(rbind,lapply(reps,function(ri){
      print(ri)
      x <- proc_sim(paste0(dir,fi,"/",mri,"/output/",ri,"/"),
                    times=tt)
      xhi <- c(apply(x$x,2,function(xii){
        sum(xii[rownames(x$x)%in%pks$winhi])
      }))
      xlo <- c(apply(x$x,2,function(xii){
        sum(xii[rownames(x$x)%in%pks$winlo])
      }))
      xtot <- xhi+xlo
      #print(xhi)
      #print(xlo)
      z <- data.frame(fhi=xhi/xtot,flo=xlo/xtot,ntot=xtot,time=tt)
      z$rep <- ri
      z$misrate <- mri
      z
    }))
  }))
  z$cellLine <- fi
  return(z)
}

dir <- "data/salehi/misseg_landscape_exploration/trial_minobs_5/"
ff <- c("SA535_CISPLATIN_CombinedH_X7_l_3_d1_0_d2_0",
        "SA906a_X57_l_7_d1_0_d2_0")

df <- readRDS("figures/misseg_landscape_exploration/peak_karyotypes.Rds")
z <- do.call(rbind,lapply(ff,proc_abm_ss,df=df))
saveRDS(z,"figures/misseg_landscape_exploration/processed_swept_abm_with_fits_v2.Rds")



