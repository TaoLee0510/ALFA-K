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



compute_population_metrics <- function(metrics=c("angle", "wass"), eval_times=seq(2000,2500,100), inDir="data/main/", outPath="data/proc/summaries/train_test_metrics.Rds", delta_t = 2000, cores = 70) {
  # Load required libraries
  library(parallel)
  
  # Validate input
  if (!all(metrics %in% c("angle", "wass"))) {
    stop("Metrics must be a subset of c('angle', 'wass')")
  }
  
  # Retrieve train-test paths
  get_subdir_combinations(inDir, subdir_1="train", subdir_2="test", only_train_00000 = TRUE)
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
      train_path <- paste(df$base_path[i], df$train_path[i], sep = "/")
      test_path <- paste(df$base_path[i], df$test_path[i], sep = "/")
      
      # Process simulations
      x0 <- proc_sim(train_path, times = eval_times)
      x1 <- proc_sim(test_path, times = eval_times - delta_t)
      colnames(x1$x) <- as.numeric(colnames(x1$x))+delta_t
      
      # Compute metrics
      a <- w <- c()
      if ("angle" %in% metrics) {
        a <- sapply(eval_times, function(t) angle_metric(x0, x1, t = t))
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
get_eval_times <- function(pathx,rangex){
  fx <- list.files(pathx)
  fx <- fx[!fx%in%c("log.txt","summary.txt")]
  tx <- as.numeric(gsub(".csv","",fx))
  if(length(rangex)==1){
    tmp <- abs(tx-rangex)
    return(tx[tmp==min(tmp)])
  }
  rangex <- rangex[order(rangex)]
  return(tx[tx>rangex[1]&tx<rangex[2]])
  
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
  return(m)
}

wasserstein_comps <- function(subdir_1="train",subdir_2="train",range_1 = 2000, 
                              range_2=c(2000,3000), cores = 70,
                              inDir="data/main/", only_train_00000 =T,
                              outPath="data/proc/summaries/train_train_matrices.Rds"){
  
  # Set up parallel cluster
  cl <- makeCluster(cores)
  
  # Source necessary scripts on all workers
  clusterCall(cl, function() {
    source("utils/comparison_functions.R")
    source("utils/ALFA-K.R")
    library(transport)
  })
  
  df <- get_subdir_combinations(inDir, subdir_1=subdir_1, subdir_2=subdir_2, 
                                only_train_00000 = only_train_00000)
  
    
  # Export variables to the cluster
  clusterExport(cl, varlist = c("eval_range", "delta_t", "df","wasserstein_matrix",
                                "get_eval_times"), envir = environment())
  

  
  # Parallel processing
  res <- parLapplyLB(cl, 1:nrow(df), function(i) {
    tryCatch({
      # Extract paths
      path1 <- paste(df$base_path[i], df$train_path[i], sep = "/")
      path2 <- paste(df$base_path[i], df$test_path[i], sep = "/")
      m <- wasserstein_matrix(path1,path2,range_1,range_2)
    }, error = function(e) {
      # Return a vector of NAs to ensure consistent output
      return(NA)
    })
  })
  list(res=res,df=df)
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

