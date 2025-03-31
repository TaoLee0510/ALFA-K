
# Function to resample species counts at each timepoint using multinomial sampling
bootstrap_counts <- function(data) {
  num_species <- nrow(data)
  num_timepoints <- ncol(data)
  boot_data <- matrix(NA, nrow = num_species, ncol = num_timepoints)
  rownames(boot_data) <- rownames(data)
  for (i in 1:num_timepoints) {
    total_counts <- sum(data[, i])
    if (total_counts == 0) {
      boot_data[, i] <- rep(0, num_species)
    } else {
      boot_data[, i] <- rmultinom(1, size = total_counts, prob = data[, i] / total_counts)
    }
  }
  boot_data
}

# Compute finite-difference dx/dt using normalized x
compute_dx_dt <- function(x, timepoints) {
  num_species <- nrow(x)
  dx_dt <- matrix(NA, nrow = num_species, ncol = ncol(x) - 1)
  
  for (i in 1:num_species) {
    dx_dt[i, ] <- diff(x[i, ]) / diff(timepoints)
  }
  dx_dt
}

# Simple helper to avoid numerical overflow
logSumExp <- function(v) {
  m <- max(v)
  m + log(sum(exp(v - m)))
}

find_birth_times <- function(opt_res, time_range,minF=0.001) {
  f_est <- opt_res$f
  x0_est <- opt_res$x0
  num_species <- length(f_est)
  
  birth_times <- rep(NA, num_species)
  
  for (i in seq_len(num_species)) {
    if (f_est[i] <= min(f_est)) next  # Skip least fit clone (never extinct)
    
    birth_fn <- function(t) {
      log_x_t <- log(x0_est[i]) + f_est[i] * t
      denom <- logSumExp(log(x0_est) + f_est * t)
      exp(log_x_t - denom) - minF
    }
    
    # Find root numerically in a wide time range
    root <- try(uniroot(birth_fn, range(time_range), tol = 1e-6)$root, silent = TRUE)
    
    if (!inherits(root, "try-error")) {
      birth_times[i] <- root
    }
  }
  
  birth_times
}

project_forward_log <- function(x0, f, timepoints) {
  # x0 and f are length-K vectors
  # returns a K x length(timepoints) matrix of frequencies
  K <- length(x0)
  out <- matrix(NA, nrow = K, ncol = length(timepoints))
  log_x0 <- log(x0)
  for (i in seq_along(timepoints)) {
    lv <- log_x0 + f * timepoints[i]
    denom <- logSumExp(lv)
    out[, i] <- exp(lv - denom)
  }
  out
}

# Infer x0 by least squares (used as a quick initial guess)
optimize_initial_frequencies <- function(x_obs, f, timepoints) {
  loss_function <- function(log_x0) {
    # x0 must sum to 1
    x0 <- exp(log_x0)
    x0 <- x0 / sum(x0)
    x_pred <- project_forward_log(x0, f, timepoints)
    sum((x_pred - x_obs)^2)
  }
  x_ini <- x_obs[,1] + 1e-6
  x_ini <- x_ini / sum(x_ini)
  opt_result <- optim(log(x_ini), loss_function, method = "BFGS")
  x0_opt <- exp(opt_result$par)
  x0_opt / sum(x0_opt)
}

# Negative log-likelihood under a multinomial model
# Param will be of length (K - 1) + K = 2K - 1
#   first (K-1) entries => free parameters for f (last is -sum)
#   next K entries => log_x0
neg_log_lik <- function(param, counts, timepoints) {
  K <- nrow(counts)
  Tt <- ncol(counts)
  
  # Parse out f
  f_free <- param[1:(K-1)]
  f_full <- c(f_free, -sum(f_free))  # enforce sum(f)=0
  
  # Parse out x0
  log_x0 <- param[K:(2*K - 1)]
  # x_k(t) = exp( log_x0[k] + f[k]*t ) / sum_j exp(...)
  
  nll <- 0
  for (i in seq_len(Tt)) {
    # total count at time i
    lv <- log_x0 + f_full * timepoints[i]
    denom <- logSumExp(lv)
    # Add - sum_k [ count_{k,i} * log( x_k(t_i) ) ]
    # = - sum_k [ c_{k,i} * (lv[k] - denom) ]
    # = - [ sum_k c_{k,i} * lv[k] ] + sum_k [ c_{k,i} ] * denom
    # But we can do it directly:
    for (k in seq_len(K)) {
      if (counts[k, i] > 0) {
        nll <- nll - counts[k, i] * (lv[k] - denom)
      }
    }
  }
  nll
}

# Joint optimization of f and x0 by maximizing the multinomial likelihood
# (equivalently minimizing negative log-likelihood).
joint_optimize <- function(counts, timepoints, f_init, x0_init) {
  K <- length(f_init)
  # param: (K-1) + K = 2K - 1
  # build initial guess
  f_free_init <- f_init[1:(K-1)]  # we ignore the last, which will be -sum
  x0_init_log <- log(x0_init + 1e-12)
  param_init <- c(f_free_init, x0_init_log)
  
  obj_fun <- function(par) neg_log_lik(par, counts, timepoints)
  
  # Optimize with BFGS
  opt <- optim(
    par = param_init,
    fn = obj_fun,
    method = "BFGS",
    control = list(maxit = 200, reltol = 1e-8)
  )
  
  # Parse out final parameter
  f_free_opt <- opt$par[1:(K-1)]
  f_opt <- c(f_free_opt, -sum(f_free_opt))
  
  log_x0_opt <- opt$par[K:(2*K - 1)]
  x0_opt <- exp(log_x0_opt)
  x0_opt <- x0_opt / sum(x0_opt)
  
  list(f = f_opt, x0 = x0_opt)
}

gen_nn_info <- function(fq,pm=0.00005){
  # Generate all single-step neighbors
  nn <- gen_all_neighbours(fq) 
  nn_str <- as.character(apply(nn, 1, paste, collapse = "."))
  
  n_info <- lapply(nn_str, function(ni) {
    # Get parents (frequent clones adjacent to neighbor)
    nj <- as.character(apply(gen_all_neighbours(ni), 1, paste, collapse = "."))
    nj <- nj[nj %in% fq]
    # Compute missegregation probabilities for each parent
    nivec <- s2v(ni)
    pij <- sapply(nj, function(si) {
      si_vec <- as.numeric(unlist(strsplit(si, split = "[.]")))
      prod(sapply(1:length(si_vec), function(k) pij(si_vec[k], nivec[k], pm)))
    })
    
    
    list(ni=ni,nj=nj,pij=pij)
  })
  n_info
}

solve_fitness_bootstrap <- function(data, minobs, nboot = 1000, epsilon = 1e-6,pm=0.00005) {
  # identify rows (karyotypes) with total < minobs
  fq <- rownames(data$x)[rowSums(data$x) > minobs]
  nn <- gen_nn_info(fq,pm)
  fq_vec <- do.call(rbind,lapply(fq,s2v))
  fq_nn <- which(as.matrix(dist(fq_vec))==1)
  timepoints <- as.numeric(colnames(data$x)) * data$dt
  num_species <- length(fq)
  num_timepoints <- ncol(data$x)
  # Matrices to store results
  f_initial_mat <- matrix(NA, nrow = nboot, ncol = num_species)
  f_final_mat   <- matrix(NA, nrow = nboot, ncol = num_species)
  x0_initial_mat <- matrix(NA, nrow = nboot, ncol = num_species)
  x0_final_mat <- matrix(NA, nrow = nboot, ncol = num_species)
  f_nn_mat <- matrix(NA,nrow=nboot,ncol=length(nn))
  
  colnames(f_initial_mat) <- fq
  colnames(f_nn_mat) <- sapply(nn,function(nni) nni$ni)
  colnames(f_final_mat)   <- fq
  colnames(x0_initial_mat) <- fq
  colnames(x0_final_mat) <- fq
  
  # Main bootstrap loop
  for (b in 1:nboot) {
    # 1) Resample counts
    boot_data <- bootstrap_counts(data$x)
    # 2) Convert to frequencies
    x <- apply(boot_data[fq,,drop=F], 2, function(col) col / sum(col))
    
    # 3) dx/dt
    dx_dt <- compute_dx_dt(x, timepoints)
    x_trim <- x[, -1, drop = FALSE]
    
    # 4) Build QP for f
    Q_accum <- matrix(0, nrow = num_species, ncol = num_species)
    r_accum <- rep(0, num_species)
    
    for (t in 1:(num_timepoints - 1)) {
      xt <- x_trim[, t]
      M_t <- diag(xt) - outer(xt, xt)
      Q_accum <- Q_accum + M_t %*% M_t
      r_accum <- r_accum + M_t %*% dx_dt[, t]
    }
    
    Dmat_boot <- 2 * Q_accum + diag(epsilon, num_species)
    dvec_boot <- 2 * r_accum
    A <- matrix(1, nrow = num_species, ncol = 1)
    bvec <- 0
    
    qp_sol <- solve.QP(Dmat_boot, dvec_boot, A, bvec, meq = 1)
    f_qp <- qp_sol$solution
    
    # Infer x0 for these f
    x0_init <- optimize_initial_frequencies(x, f_qp, timepoints)
    
    # Save these "initial" estimates
    f_initial_mat[b, ]  <- f_qp
    x0_initial_mat[b, ] <- x0_init
    
    # 5) Joint optimization of f and x0 to maximize likelihood
    opt_res <- joint_optimize(boot_data[fq,,drop=F], timepoints, f_qp, x0_init)
    birth_times <- find_birth_times(opt_res,time_range=c(-1000,max(timepoints)),minF=0.001)
    peak_times <- timepoints[apply(x,1,which.max)]
    mean_risetime <- mean(peak_times-birth_times,na.rm=T)
    birth_times[is.na(birth_times)] <- peak_times[is.na(birth_times)]-mean_risetime
    
    x0par <- opt_res$x0
    names(x0par) <- fq
    fpar <- opt_res$f
    names(fpar) <- fq
    names(birth_times) <- fq
    
    dfb <- as.matrix(dist(as.numeric(fpar)))
    dfb[upper.tri(dfb)] <- dfb[upper.tri(dfb)]*(-1)
    sdy <- sd(dfb[fq_nn])
    
    f_final_mat[b, ] <- opt_res$f
    x0_final_mat[b,] <- opt_res$x0 
    
    fExp <- function(fc,fp,pij,tt){ ##handle multiple clones
      pij*fp/(fc-fp)*(exp(tt*(fc-fp))-1)
    }
    xfit <- project_forward_log(x0par,fpar,timepoints)
    rownames(xfit) <- fq
    ntot <- colSums(boot_data)
    
    opt_fc <- function(fc,nni,delta_f){
      child <- nni$ni
      xc_est <- colSums(do.call(rbind,lapply(1:length(nni$nj),function(i){
        par <- nni$nj[i]
        tt <- timepoints-birth_times[par]
        ## so this code is sensitive to the exact value of fp and will fail for fp<0.
        ## we dont know how much dividing is actually going on, so makes sense to fix this
        ## value for all parents? That would correspond to a case where all division rates
        ## are equal and fitness advantages are a reduction in death rates? Other alternative
        ## is to add a fixed offset to all values... may need more thought.
        fExp(fc,fpar[par]+delta_f,nni$pij[i],tt)*xfit[par,]
      })))
      xc_est <- pmax(0,pmin(1,xc_est))
      xc_obs <- rep(0,length(timepoints))
      if(child%in%rownames(boot_data)) xc_obs <- boot_data[nni$ni,] 
      -sum(c(dbinom(xc_obs,ntot,prob=xc_est,log=T),dnorm(fc-0.5,mean = 0,sd=sdy,log=T)))
    }
    
    search_interval <- range(fpar)
    interval_range <- diff(range(search_interval))
    search_interval[1] <- search_interval[1] - interval_range
    search_interval[2] <- search_interval[2] + interval_range
    
    delta_f <- 0.5
    
    fc <- sapply(nn,function(nni){
      res <- optimise(opt_fc,interval = search_interval+delta_f,nni=nni,delta_f=delta_f)
      res$minimum-delta_f
    })
    
    f_nn_mat[b,] <- fc
    
    
  
  }
  
  # Return the three matrices (and nothing else)
  list(
    initial_fitness      = f_initial_mat, 
    final_fitness        = f_final_mat,
    initial_frequencies  = x0_initial_mat,
    final_frequencies  = x0_final_mat,
    nn_fitness = f_nn_mat
  )
}

proc_sweep_input <- function(yi,ntp=8,tmax=1200){
  yi$x <- yi$x[, as.numeric(colnames(yi$x)) < 1200]
  yi$x <- yi$x[, (ncol(yi$x) - ntp + 1):ncol(yi$x)]
  is.absent <- rowSums(yi$x)==0
  yi$clone.fitness <- yi$clone.fitness[!is.absent]
  yi$x <- yi$x[!is.absent,]
  names(yi$clone.fitness) <- row.names(yi$x)
  return(yi)
}

test_fit_quality <- function(yi,nboot=30,minobs=20,ntp=4){
  tryCatch({
    data  <- proc_sweep_input(yi$data,ntp=ntp)
    fq_boot <- solve_fitness_bootstrap(data, minobs=minobs, nboot = nboot)
    kfq <- do.call(rbind,lapply(colnames(fq_boot$initial_fitness),s2v))
    ftru <- sapply(1:nrow(kfq),function(i){
      ki <- kfq[i,]
      getf(ki,yi$lscape)
    })
    fpred_ini <- apply(fq_boot$initial_fitness,2,median,na.rm=T)
    fpred_final <- apply(fq_boot$final_fitness,2,median,na.rm=T)
    
    data.frame(nfq=length(ftru),nboot,minobs,ntp,Rsqini=R2R(ftru,fpred_ini),rini = cor(ftru,fpred_ini),
               Rsqfin = R2R(ftru,fpred_final),rfin=cor(ftru,fpred_final))
  },error=function(e) return(NULL))
  
}



wrap_fit <- function(yi,minobs,ntp,nboot){
  data  <- proc_sweep_input(yi$data,ntp=ntp)
  fq_boot <- solve_fitness_bootstrap(data, minobs=minobs, nboot = nboot)
  
  fq_str <- colnames(fq_boot$final_fitness)
  nn_str <- colnames(fq_boot$nn_fitness)
  
  ktrain <- do.call(rbind,lapply(c(fq_str,nn_str),s2v))
  fboot <- cbind(fq_boot$final_fitness,fq_boot$nn_fitness)
  
  ktest <- unique(rbind(k,gen_all_neighbours(nn_str)))
  
  ktest_str <- apply(ktest,1,paste,collapse=".")
  fq_ids <- ktest_str%in%fq_str
  nn_ids <- ktest_str%in%nn_str
  lscape <- yi$lscape
  
  f <- apply(ktest,1,function(ki) getf(ki,lscape))
  

  boot_predictions <- matrix(0, nrow = nrow(ktest), ncol = nboot)
  
  
  for (b in 1:nboot) {
    # Sample a different realization from fboot for each point
    boot_f <- fboot[cbind(sample(1:nrow(fboot), ncol(fboot), replace = TRUE), 1:ncol(fboot))]
    
    # Fit Kriging model on this sampled dataset
    fit_boot <- Krig(ktrain, boot_f, cov.function = "stationary.cov", cov.args = list(Covariance = "Matern", smoothness = 1.5))
    
    # Store predictions
    boot_predictions[, b] <- predict(fit_boot, ktest)
  }
  
  pred_means <- apply(boot_predictions,1,mean)
  pred_medians <- apply(boot_predictions,1,median)
  pred_sd <- apply(boot_predictions,1,sd)
  
  landscape <- data.frame(k=ktest_str,mean=pred_means,median=pred_medians,sd=pred_sd,fq=fq_ids,nn=nn_ids)
  
  summary <- data.frame(r=cor(pred_means,f),
                        rfq=cor(pred_means[fq_ids],f[fq_ids]),
                        rnn=cor(pred_means[nn_ids],f[nn_ids]),
                        rd2=cor(pred_means[!(fq_ids|nn_ids)],f[!(fq_ids|nn_ids)]),
                        R=R2R(f,pred_means),
                        Rfq=R2R(f[fq_ids],pred_means[fq_ids]),
                        Rnn=R2R(f[nn_ids],pred_means[nn_ids]),
                        Rd2=R2R(f[!(fq_ids|nn_ids)],pred_means[!(fq_ids|nn_ids)]))
  
  list(summary=summary,landscape=landscape)
}


library(parallel)

# Set up cluster
num_cores <- 70  # Use one less than the max cores
cl <- makeCluster(num_cores)

# Export necessary objects, including all functions
setwd("~/projects/ALFA-K/")
y <- readRDS("data/proc/sweep_inputs.Rds")
pars <- expand.grid(id=names(y),minobs=c(5,10,20),ntp=c(2,4,8))
clusterExport(cl, ls())  # This exports *everything* in the current environment

# Ensure each worker loads required libraries and sources files
clusterEvalQ(cl, {
  library(quadprog)
  library(nloptr)
 setwd("~/projects/ALFA-K/")
  source("utils/ALFA-K.R")
  source("utils/comparison_functions.R")
})

# Run in parallel
x <- parLapplyLB(cl, 1:nrow(pars), function(i) {
  yi <- y[[pars$id[i]]]
  minobs <- pars$minobs[i]
  ntp <- pars$ntp[i]
  tryCatch(wrap_fit(yi, minobs, ntp, nboot=30), error = function(e) return(NULL))
})

# Stop cluster
stopCluster(cl)

# Collect results
summary <- do.call(rbind, lapply(x, function(xi) xi$summary))
lscapes <- lapply(x, function(xi) xi$landscape)

saveRDS(summary,"data/proc/boot_summary.Rds")
saveRDS(lscapes,"data/proc/boot_landscapes.Rds")

print(0)






