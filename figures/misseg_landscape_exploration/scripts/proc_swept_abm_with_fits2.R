setwd("~/projects/ALFA-K")

#setwd("~/projects/008_birthrateLandscape/ALFA-K")

library(parallel)
cl <- makeCluster(getOption("cl.cores", 8))

proc_abm_ss <- function(tmp,df,dir){
  fi <- tmp$ff
  mri <- tmp$mr 
  source("utils/ALFA-K.R")  
  pks <- df[[fi]]

    reps <- list.files(paste0(dir,fi,"/",mri,"/output/"))
    # print(reps)
    do.call(rbind,lapply(reps,function(ri){
      #print(ri)
      print(paste(fi,mri,ri))
      #print(tt)
      target_dir <- paste0(dir,fi,"/",mri,"/output/",ri,"/")
      tt <- list.files(target_dir)
      tt<- as.numeric(gsub(".csv","",tt))
      tt <- tt[!is.na(tt)]      
      x <- proc_sim(target_dir,
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
      z <- z[!z$ntot==0,]
      z$rep <- ri
      
      z
    }))
  z$misrate <- mri
  z$cellLine <- fi
  return(z)
}

dir <- "data/salehi/misseg_landscape_exploration/trial_minobs_5/"
ff <- c("SA535_CISPLATIN_CombinedH_X7_l_3_d1_0_d2_0",
        "SA906a_X57_l_7_d1_0_d2_0")

tmp <- do.call(rbind,lapply(ff,function(ffi) {
  tmp <- list.files(paste0(dir,ffi))
  cbind(ffi,tmp)
}))

conds <- lapply(1:nrow(tmp),function(i){
  list(ff=as.character(tmp[i,1]),
       mr=as.character(tmp[i,2]))
})

df <- readRDS("figures/misseg_landscape_exploration/data/peak_karyotypes.Rds")
z <- do.call(rbind,parLapplyLB(cl=cl,X=conds,proc_abm_ss,df=df,dir=dir))
saveRDS(z,"figures/misseg_landscape_exploration/data/processed_swept_abm_with_fits_v2.Rds")



