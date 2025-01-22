## various functions for aggregating ABM sweep output.

## this function is used when we have data structure
##  inDir
##    subdir
##      target
## where target is a data frame.
aggregate_fit_summaries <- function(summaryName="fit_summaries.Rds",inDir="data/main/",outPath="data/proc/summaries/"){
  ## make sure directory exists
  dir.create(outPath,recursive = T,showWarnings = F)
  ff <- list.files(inDir)
  ff <- paste0(inDir,ff,"/",summaryName)
  ff <- ff[file.exists(ff)]
  x <- lapply(ff,readRDS)
  saveRDS(x,paste0(outPath,summaryName))
  return(x)
}

##20 mins 70 cores
apply_landscape_metrics <- function(mainDir="data/main/",cfig_path="config.txt",
                                    lscape_path="landscape.txt",fit_path="sweep_fits",
                                    outPath="data/proc/summaries/landscape_metrics.Rds",cores=70){
  library(parallel)
  
  # Set up parallel cluster
  cl <- makeCluster(cores)
  
  # Source necessary scripts on all workers
  clusterCall(cl, function() {
    source("utils/comparison_functions.R")
    source("utils/ALFA-K.R")
  })
  
  # Export variables to the cluster
  clusterExport(cl, varlist = c("mainDir","cfig_path","lscape_path","fit_path"), envir = environment())
  dirs <- list.files(mainDir)
  
  df <- do.call(rbind,parLapplyLB(cl, dirs, function(di) {
    landscape_path <- paste(mainDir,di,lscape_path,sep="/")
    config_path <- paste(mainDir,di,cfig_path,sep="/")
    fobj <- gen_fitness_object(cfig_path = paste(mainDir,di,cfig_path,sep="/"),
                               lscape_path = paste(mainDir,di,lscape_path,sep="/"))
    
    fit_path_i <- paste(mainDir,di,fit_path,sep="/")
    targets <- list.files(fit_path_i)
    
    df_i <- do.call(rbind,lapply(targets,function(ti){
      tarpath <- paste(fit_path_i,ti,sep="/")
      krigfit <- readRDS(tarpath)
      if("input"%in%names(krigfit)) krigfit <- krigfit$fit
      k <- rbind(do.call(rbind,lapply(rownames(krigfit$xo),s2v)), gen_all_neighbours(rownames(krigfit$xo)))
      
      df <- data.frame(f_est=predict(krigfit$fit,k),
                       f_tru=apply(k, 1, function(xi) getf(xi,fobj)))
      
      res <- landscape_accuracy_metrics(df)
      res$fit_id <- gsub(".Rds","",ti)
      return(res)
    }))
    df_i$sim_id <- di
    return(df_i)
  }))
  saveRDS(df,outPath)
}

get_subdir_combinations <- function(base_dir, subdir_1, subdir_2, only_train_00000 = FALSE) {
  dirs_A <- list.files(base_dir, full.names = TRUE)
  
  do.call(rbind, lapply(dirs_A, function(dir_A) {
    path_1 <- file.path(dir_A, subdir_1)
    path_2 <- file.path(dir_A, subdir_2)
    
    if (!dir.exists(path_1) || !dir.exists(path_2)) return(NULL)
    
    # Find all bottom-level directories in subdir_1
    path_1_dirs <- list.dirs(path_1, full.names = TRUE, recursive = TRUE)
    path_1_dirs <- path_1_dirs[sapply(path_1_dirs, function(d) {
      length(list.dirs(d, recursive = FALSE)) == 0
    })]
    
    # Find all bottom-level directories in subdir_2
    path_2_dirs <- list.dirs(path_2, full.names = TRUE, recursive = TRUE)
    path_2_dirs <- path_2_dirs[sapply(path_2_dirs, function(d) {
      length(list.dirs(d, recursive = FALSE)) == 0
    })]
    
    # Apply special case filter
    if (only_train_00000) {
      if (subdir_1 == "train") {
        path_1_dirs <- path_1_dirs[basename(path_1_dirs) == "00000"]
      }
      if (subdir_2 == "train") {
        path_2_dirs <- path_2_dirs[basename(path_2_dirs) == "00000"]
      }
    }
    
    # If no bottom-level directories, return NULL
    if (length(path_1_dirs) == 0 || length(path_2_dirs) == 0) return(NULL)
    
    # Generate combinations between path_1 and path_2 bottom-level directories
    base_path <- dir_A
    relative_path_1 <- sub(paste0(base_path, "/"), "", path_1_dirs)
    relative_path_2 <- sub(paste0(base_path, "/"), "", path_2_dirs)
    
    expand.grid(
      base_path = base_path,
      path_1 = relative_path_1,
      path_2 = relative_path_2,
      stringsAsFactors = FALSE
    )
  }))
}


compute_population_metrics <- function(metrics=c("angle"), eval_times=seq(2000,3000,100), inDir="data/main/", outPath="data/proc/summaries/train_test_angles.Rds", delta_t = 2000, cores = 70) {
  # Load required libraries
  library(parallel)
  
  # Validate input
  if (!all(metrics %in% c("angle", "wass"))) {
    stop("Metrics must be a subset of c('angle', 'wass')")
  }
  
  # Retrieve train-test paths
  df <- get_subdir_combinations(inDir, subdir_1="train", subdir_2="test_v2", only_train_00000 = TRUE)
  if (nrow(df) == 0) {
    stop("No valid train-test paths found in the specified directory.")
  }
  
  # Set up parallel cluster
  cl <- makeCluster(cores)
  
  # Source necessary scripts on all workers
  clusterCall(cl, function() {
    source("utils/comparison_functions.R")
    source("utils/ALFA-K.R")
    library(transport)
  })
  
  # Export variables to the cluster
  clusterExport(cl, varlist = c("metrics", "eval_times", "delta_t", "df"), envir = environment())
  
  # Parallel processing
  res <- do.call(rbind, parLapplyLB(cl, 1:nrow(df), function(i) {
    tryCatch({
      # Extract paths
      train_path <- paste(df$base_path[i], df$path_1[i], sep = "/")
      test_path <- paste(df$base_path[i], df$path_2[i], sep = "/")
      
      # Process simulations
      x0 <- proc_sim(train_path, times = eval_times)
      x1 <- proc_sim(test_path, times = eval_times - delta_t)
      colnames(x1$x) <- as.numeric(colnames(x1$x))+delta_t
      x_ini <- get_mean(x0,1)
      # Compute metrics
      a <- w <- c()
      if ("angle" %in% metrics) {
        a <- sapply(eval_times, function(t) angle_metric(x0, x1, t = t,x0 = x_ini))
        names(a) <- paste0("a", eval_times)
      }
      if ("wass" %in% metrics) {
        w <- sapply(eval_times, function(t) wasserstein_metric(x0, x1, t = t))
        names(w) <- paste0("w", eval_times)
      }
      return(c(a, w))
    }, error = function(e) {
      # Return a vector of NAs to ensure consistent output
      return(rep(NA, length(metrics) * length(eval_times)))
    })
  }))
  
  # Stop the cluster
  stopCluster(cl)
  
  # Combine results with the data frame
  df <- cbind(df, res)
  saveRDS(df, outPath)
  return(df)
}

## gets all sim output timepoints from pathx that are in time range rangex.
## if rangex is a single number it will find the closest output time and return that.
get_eval_times <- function(pathx,rangex,limit=Inf){
  if(!is.finite(rangex[1])){
    ## this hack after the fact allows us to set e.g. rangex = c(NaN,1200)
    ## then take the final timepoint before 1200 and N=limit-1 timepoints after
    fx <- list.files(pathx)
    fx <- fx[!fx%in%c("log.txt","summary.txt")]
    tx <- as.numeric(gsub(".csv","",fx))
    t0 <- max(tx[tx<rangex[2]])
    tx <- head(tx[tx>=rangex[2]],limit-1)
    return(c(t0,tx))
  }
  if(length(rangex)>2) return(head(rangex,limit))
  fx <- list.files(pathx)
  fx <- fx[!fx%in%c("log.txt","summary.txt")]
  tx <- as.numeric(gsub(".csv","",fx))
  if(length(rangex)==1){
    tmp <- abs(tx-rangex)
    return(head(tx[tmp==min(tmp)],limit))
  }
  rangex <- rangex[order(rangex)]
  return(head(tx[tx>rangex[1]&tx<rangex[2]],limit))
  
}
wasserstein_matrix <- function(path1,path2,range1,range2){
  ## function needs the following to be sourced in the environment to work:
  ## source("utils/comparison_functions.R")
  ## source("utils/ALFA-K.R")
  x1 <- proc_sim(path1,times = get_eval_times(path1,range1))
  x2 <- proc_sim(path2,times=get_eval_times(path2,range2))
  m <- do.call(rbind,lapply(as.numeric(colnames(x1$x)),function(t1){
    sapply(as.numeric(colnames(x2$x)),function(t2){
      wasserstein_distance(test=x1,ref=x2,t=t1,t2=t2)
    })
  }))
  rownames(m) <- colnames(x1$x)
  colnames(m) <- colnames(x2$x)
  return(m)
}

metric_array <- function(path1, path2, range1, range2, metrics = c("wasserstein", "euclidean", "cosine", "jaccard"),limit=Inf,diag_only=F) {
  ## function needs the following to be sourced in the environment to work:
  ## source("utils/comparison_functions.R")
  ## source("utils/ALFA-K.R")
  
  # Process simulations for the given paths and time ranges
  x1 <- proc_sim(path1, times = get_eval_times(path1, range1,limit=limit))
  x2 <- proc_sim(path2, times = get_eval_times(path2, range2,limit=limit))
  
  # Initialize an empty list to store results for each metric
  metric_results <- lapply(metrics, function(metric) {
    do.call(rbind, lapply(as.numeric(colnames(x1$x)), function(t1) {
      sapply(as.numeric(colnames(x2$x)), function(t2) {
        ##N.B the wasserstein function takes a numeric time input then finds the nearest 
        if(diag_only){
          if(which(colnames(x1$x)==t1)!=which(colnames(x2$x)==t2)) return(NaN)
        }
        if (metric == "wasserstein") {
          wasserstein_distance(test = x1, ref = x2, t = t1, t2 = t2)
        } else if (metric == "euclidean") {
          compute_mean_karyotype_distance(x1,x2,t1,t2)
        } else if (metric == "cosine") {
          compute_cosine_similarity(x1,x2,t1,t2)
        } else if (metric == "jaccard") {
          compute_overlap_coefficient(x1,x2,t1,t2)
        } else {
          stop(paste("Unsupported metric:", metric))
        }
      })
    }))
  })
  
  # Combine results into a 3D array with named dimensions
  names(metric_results) <- metrics
  result_array <- abind::abind(metric_results, along = 3)
  dimnames(result_array) <- list(
    rownames = colnames(x1$x),
    colnames = colnames(x2$x),
    metrics = metrics
  )
  
  return(result_array)
}



metric_comps <- function(subdir_1 = "train", subdir_2 = "train", 
                         range_1 = 2000, range_2 = seq(2000, 3000, 200), 
                         metrics = c("wasserstein", "euclidean", "cosine", "jaccard"),
                         limit=Inf,diag_only=F,
                         cores = 70, inDir = "data/main/", only_train_00000 = TRUE, 
                         outPath = "data/proc/summaries/train_train_matrices.Rds") {
  library(parallel)
  
  # Set up parallel cluster
  cl <- makeCluster(cores)
  
  # Source necessary scripts on all workers
  clusterCall(cl, function() {
    source("utils/comparison_functions.R")
    source("utils/ALFA-K.R")
    library(transport)
  })
  
  # Get combinations of directories
  df <- get_subdir_combinations(inDir, subdir_1 = subdir_1, subdir_2 = subdir_2, 
                                only_train_00000 = only_train_00000)
  
  # Export variables to the cluster
  clusterExport(cl, varlist = c("range_1", "range_2", "df", "metric_array", 
                                "get_eval_times", "metrics","limit","diag_only"), envir = environment())
  
  # Parallel processing
  res <- parLapplyLB(cl, 1:nrow(df), function(i) {
    tryCatch({
      # Extract paths
      path1 <- paste(df$base_path[i], df$path_1[i], sep = "/")
      path2 <- paste(df$base_path[i], df$path_2[i], sep = "/")
      
      # Compute metric array
      metric_array(path1, path2, range_1, range_2, metrics = metrics,limit=limit,diag_only = diag_only)
    }, error = function(e) {
      # Return NULL on error to ensure consistent output
      return(NULL)
    })
  })
  
  # Stop the cluster
  stopCluster(cl)
  
  # Prepare output
  output <- list(res = res, df = df)
  saveRDS(output, outPath)
  return(output)
}


compute_wasserstein_distances <- function(subdir_1="train", subdir_2="test", 
                                          compare_only_train_00000 = TRUE,eval_times=c(2000,3000), 
                                          inDir="data/main/",
                                          outPath="data/proc/summaries/wasserstein_distances.Rds", 
                                          delta_t = 2000, cores = 70) {
  #N.B the way this function is setup, always set eval times as the "initial" timepoint & the timepoint of interest
  # Load required libraries
  library(parallel)
  
  # Retrieve train-test paths
  get_subdir_combinations(inDir, subdir_1=subdir_1, subdir_2=subdir_2, compare_only_train_00000 = compare_only_train_00000)
  if (nrow(df) == 0) {
    stop("No valid train-test paths found in the specified directory.")
  }
  
  # Set up parallel cluster
  cl <- makeCluster(cores)
  
  # Source necessary scripts on all workers
  clusterCall(cl, function() {
    source("utils/comparison_functions.R")
    source("utils/ALFA-K.R")
    library(transport)
  })
  
  # Export variables to the cluster
  clusterExport(cl, varlist = c("eval_times", "delta_t", "df"), envir = environment())
  
  # Parallel processing
  res <- data.frame(do.call(rbind, parLapplyLB(cl, 1:nrow(df), function(i) {
    tryCatch({
      # Extract paths
      train_path <- paste(df$base_path[i], df$train_path[i], sep = "/")
      test_path <- paste(df$base_path[i], df$test_path[i], sep = "/")
      
      # Process simulations
      x0 <- proc_sim(train_path, times = eval_times)
      x1 <- proc_sim(test_path, times = eval_times - delta_t)
      colnames(x1$x) <- as.numeric(colnames(x1$x))+delta_t
      
      
      d1 <- wasserstein_distance(x1,t=eval_times[2])
      d2 <- wasserstein_distance(x0,t=eval_times[2])
      d3 <- wasserstein_distance(x1,x0,t=eval_times[2])

      return(c(d1,d2,d3))
    }, error = function(e) {
      # Return a vector of NAs to ensure consistent output
      return(rep(NA, 3))
    })
  })))
  colnames(res) <- c("d_pred","d_ref","d_ref_pred")
  
  # Stop the cluster
  stopCluster(cl)
  
  # Combine results with the data frame
  df <- cbind(df, res)
  df$eval_time <- tail(eval_times,1)
  saveRDS(df, outPath)
  return(df)
}


f <- function(mainDir="data/main/",cores=70,outDir="data/proc/summaries/f2000.Rds"){
  
  cl <- makeCluster(cores)
  clusterCall(cl, function() {
    source("utils/comparison_functions.R")
    source("utils/ALFA-K.R")
    library(transport)
  })
  
  # Export variables to the cluster
  clusterExport(cl, varlist = c("mainDir"), envir = environment())
  ff <- list.files(mainDir)
  
  res <- data.frame(do.call(rbind, parLapplyLB(cl, ff, function(fi) {
    tryCatch({
      inDir <- paste0(mainDir,fi,"/")
      x <- proc_sim(paste0(inDir,"train/00000/"),times=2000)
      xin <- do.call(rbind,lapply(rownames(x$x),s2v))
      fits <- list.files(paste0(inDir,"sweep_fits/"))
      
      preds <- do.call(cbind,lapply(fits,function(fi){
        ai <- readRDS(paste0(inDir,"sweep_fits/",fi))
        fest <- predict(ai$fit,xin)
        c(fest)
      }))
      
      preds <- data.frame(preds)
      fits <- gsub("_00000.Rds","",fits)
      colnames(preds) <- fits
      
      df <- data.frame(N=x$x[,1],f=x$clone.fitness,row.names = NULL)
      preds <- preds[order(df$N,decreasing=T),]
      df <- df[order(df$N,decreasing=T),]
      df$base_dir <- fi
      df <- cbind(df,preds)
      
      return(df)
    }, error = function(e) {
      # Return a vector of NAs to ensure consistent output
      return(NULL)
    })
  })))
  
  
  saveRDS(res,outDir)
  
  
}


aggregate_salehi_preds <- function(inDir="data/salehi/forward_sims/minobs_20/",outPath="data/proc/summaries/salehi_preds_minobs_20.Rds",cores=70,times=seq(0,200,10)){
  library(parallel)
  cl <- makeCluster(cores)
  clusterCall(cl, function() {
    source("utils/ALFA-K.R")
  })
  
  # Export variables to the cluster
  clusterExport(cl, varlist = c("inDir","times"), envir = environment())
  source("utils/ALFA-K.R")
  ff <- list.files(inDir) 
  res <- parLapplyLB(cl, ff, function(fi) {
    x <- tryCatch({
      ## seems that sometimes the predicted population goes extinct. 
      subdir <- paste(inDir,fi,"output",sep="/")
      output_folders <- list.files(subdir)
      x <- lapply(output_folders,function(oi){
        proc_sim(paste(subdir,oi,sep="/"),times=times)
      })
      return(x)
    },error=function(e) return(NULL))
    
    return(x)
  })
  names(res) <- ff
  saveRDS(res,outPath)
}

get_evo_rate <- function(mainDir="data/main/"){
  ff <- list.files(mainDir)
  tarDir <- "train/00000/"
  source("utils/ALFA-K.R")
  get_tt <- function(path){
    tt <- list.files(path)
    tt <- gsub(".csv","",tt)
    tt <- as.numeric(tt[!tt%in%c("log.txt","summary.txt")])
  }
  x <- pbapply::pblapply(ff,function(fi){
    path <- paste(mainDir,fi,tarDir,sep="/")
    tt <- get_tt(path)
    x <- proc_sim(path,tt)
    data.frame(tt=tt,fitness=x$pop.fitness,id=fi)
  })
  saveRDS(x,"data/proc/summaries/train_fitness.Rds")
}

get_xval <- function(mainDir="data/main",subDir="alfak_eq_space_fits",outPath="data/proc/summaries/sweep_xval_eq_space_fits.Rds"){
  setwd("~/projects/ALFA-K")
  ci <- list.files(mainDir)
  fpaths <- paste(mainDir,ci,subDir,sep="/")
  x <- pbapply::pblapply(fpaths,function(f){
    fits <- list.files(f)
    fitnames <- sapply(fits,function(fi){
      tail(unlist(strsplit(fi,split=".Rds")),1)
    })
    fits <- lapply(paste(f,fits,sep="/"),function(fit_path){
      x <- readRDS(fit_path)
      x$xval
    })
    names(fits) <- fitnames
    return(fits)
  })
  names(x) <- ci
  saveRDS(x,outPath)
  return(x)
  
}
