setwd("~/projects/ALFA-K/data/proc/sweep_results/")

conditions <- list.files()

x <- pbapply::pblapply(conditions,function(cj){
  fi <- list.files(cj)
  lapply(fi,function(fij){
    if(!file.exists(paste0(cj,"/",fij,"/result_summary.Rds"))) return(NULL)
    return(readRDS(paste0(cj,"/",fij,"/result_summary.Rds")))
  })
  
})

dfl <- do.call(rbind,lapply(x,function(xi){
  do.call(rbind,lapply(xi,function(xij){
    cbind(xij$summary,xij$info)
  }))
}))

reformat_df <- function(df_,metric){
  colnames(df_) <- c("prediction","baseline")
  df_$passage <- 1:nrow(df_)
  df_$metric <- metric
  df_
}

dfp <- do.call(rbind,lapply(x,function(xi){
  do.call(rbind,lapply(xi,function(xij){
    tryCatch({
      tmp <- rbind(cbind(reformat_df(xij$df_so,"overlap"),xij$df_a),
                   cbind(reformat_df(xij$df_sc,"cosine"),xij$df_a),
                   cbind(reformat_df(xij$df_de,"euclidean"),xij$df_a),
                   cbind(reformat_df(xij$df_dw,"wasserstein"),xij$df_a))
      
      cbind(tmp,xij$info,xij$summary)
    },error=function(e) return(NULL))
  }))
}))

dfp$win <- FALSE
dfp$win[dfp$metric%in%c("overlap","cosine")] <- (dfp$prediction>dfp$baseline)[dfp$metric%in%c("overlap","cosine")] 
dfp$win[!dfp$metric%in%c("overlap","cosine")] <- (dfp$prediction<dfp$baseline)[!dfp$metric%in%c("overlap","cosine")]
dfp$pos_xv <- dfp$Rxv>0

res <- list(dfl=dfl,dfp=dfp)
saveRDS(res,"../sweep_results_summaries.Rds")
