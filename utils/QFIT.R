library(quadprog)
library(nloptr)
library(ggplot2)
library(truncnorm)
source("utils/ALFA-K.R")
setwd("~/projects/008_birthrateLandscape/ALFA-K/")
# Function to resample species counts at each timepoint using multinomial sampling
bootstrap_counts <- function(data) {
  num_species <- nrow(data)
  num_timepoints <- ncol(data)
  boot_data <- matrix(NA, nrow = num_species, ncol = num_timepoints)
  
  for (i in 1:num_timepoints) {
    total_counts <- sum(data[, i])  # Total counts at timepoint i
    if (total_counts == 0) {
      boot_data[, i] <- rep(0, num_species)  # Prevent division by zero
    } else {
      boot_data[, i] <- rmultinom(1, size = total_counts, prob = data[, i] / total_counts)
    }
  }
  return(boot_data)
}

# Function to compute finite-difference dx/dt using normalized x
compute_dx_dt <- function(x, timepoints) {
  num_species <- nrow(x)
  num_timepoints <- ncol(x)
  dx_dt <- matrix(NA, nrow = num_species, ncol = num_timepoints - 1)
  
  for (i in 1:num_species) {
    dx_dt[i, ] <- diff(x[i, ]) / diff(timepoints)
  }
  
  return(dx_dt)
}



solve_fitness_bootstrap <- function(data, minobs, nboot = 1000, epsilon = 1e-6) {
  data$x <- data$x[rowSums(data$x) > minobs, ]
  timepoints <- as.numeric(colnames(data$x)) * data$dt
  data <- data$x
  num_species <- nrow(data)
  num_timepoints <- ncol(data)
  f_boot <- matrix(NA, nrow = nboot, ncol = num_species)
  
  infer_initial_frequencies <- function(x_final, f, total_time) {
    log_S <- log(sum(x_final * exp(-f * total_time)))
    x_initial <- x_final * exp(-f * total_time - log_S)
    return(x_initial)
  }
  
  project_forward_log <- function(x0, f, timepoints) {
    log_x0 <- log(x0)
    projected_x <- matrix(NA, nrow = length(x0), ncol = length(timepoints))
    
    for (t in seq_along(timepoints)) {
      log_x_t <- log_x0 + f * timepoints[t] - log(sum(exp(log_x0 + f * timepoints[t])))
      projected_x[, t] <- exp(log_x_t)
    }
    
    return(projected_x)
  }
  
  optimize_initial_frequencies <- function(x_obs, f, timepoints) {
    loss_function <- function(log_x0) {
      x0 <- exp(log_x0) / sum(exp(log_x0))  # Normalize to sum to 1
      x_pred <- project_forward_log(x0, f, timepoints)
      return(sum((x_pred - x_obs)^2))  # Error across all timepoints
    }
    xini <- x_obs[,1]+1e-6
    xini <- xini/sum(xini)
    log_x0_init <- log(xini)  # Start from first observed frequencies
    opt_result <- optim(log_x0_init, loss_function, method = "BFGS")
    x0_opt <- exp(opt_result$par) / sum(exp(opt_result$par))  # Normalize again
    
    return(x0_opt)
  }
  
  for (b in 1:nboot) {
    # Step 1: Resample species counts at each timepoint
    boot_data <- bootstrap_counts(data)
    
    # Step 2: Compute new relative abundances (frequencies) per timepoint
    x <- apply(boot_data, 2, function(col) col / sum(col))
    
    # Step 3: Compute finite-difference dx/dt using forward differences
    dx_dt <- compute_dx_dt(x, timepoints)
    
    # Step 4: Trim x to match dx/dt dimensions
    x_trim <- x[, -1, drop = FALSE]  # dimensions: num_species x (num_timepoints - 1)
    
    # Step 5: Construct QP matrices without resampling timepoints
    Q_accum <- matrix(0, nrow = num_species, ncol = num_species)
    r_accum <- rep(0, num_species)
    
    for (t in 1:(num_timepoints - 1)) {
      xt <- x_trim[, t]
      M_t <- diag(xt) - outer(xt, xt)
      Q_accum <- Q_accum + M_t %*% M_t
      r_accum <- r_accum + M_t %*% dx_dt[, t]
    }
    
    # Step 6: Solve the QP problem
    Dmat_boot <- 2 * Q_accum + diag(epsilon, num_species)  # Add ridge term
    dvec_boot <- 2 * r_accum
    
    # Step 7: Constraint: sum(f) = 0 (fix gauge)
    A <- matrix(1, nrow = num_species, ncol = 1)  # sum(f) = 0 constraint
    bvec <- 0
    
    # Solve the QP
    qp_sol <- solve.QP(Dmat_boot, dvec_boot, A, bvec, meq = 1)
    fbb <- qp_sol$solution
    x0 <- optimize_initial_frequencies(x,fbb,timepoints)
    f_boot[b, ] <- fbb
  }
  colnames(f_boot) <- rownames(data)
  
  return(list(
    fitness = f_boot,
    infer_initial_frequencies = infer_initial_frequencies,
    optimize_initial_frequencies = optimize_initial_frequencies
  ))
}




##version where we assume start time is the start of our sample data.
boot_nn <- function(f_boot_samples, x, pm = 0.00005) {

  xfq <- x$x[colnames(f_boot_samples),]
  xfq <- apply(xfq,2,function(xi) xi/sum(xi))
  
  fq_vec <- do.call(rbind,lapply(colnames(f_boot_samples),s2v))
  fq_nn <- which(as.matrix(dist(fq_vec))==1)
  
  # Generate all single-step neighbors
  nn <- gen_all_neighbours(colnames(f_boot_samples)) 
  nn_str <- as.character(apply(nn, 1, paste, collapse = "."))
  
  n_info <- lapply(nn_str, function(ni) {
    # Get parents (frequent clones adjacent to neighbor)
    nj <- as.character(apply(gen_all_neighbours(ni), 1, paste, collapse = "."))
    nj <- nj[nj %in% colnames(f_boot_samples)]
    # Compute missegregation probabilities for each parent
    nivec <- s2v(ni)
    pij <- sapply(nj, function(si) {
      si_vec <- as.numeric(unlist(strsplit(si, split = "[.]")))
      prod(sapply(1:length(si_vec), function(k) pij(si_vec[k], nivec[k], pm)))
    })
    
    
    np <- x$x[nj,,drop=F]
    nc <- rep(0,ncol(x$x))
    if(ni%in%rownames(x$x)) nc <- x$x[ni,]
    
    list(np=np,nc=nc,kp=nj,pij=pij)
  })
  names(n_info) <- nn_str
  
  tt <- as.numeric(colnames(x$x))*x$dt
  tt <- tt-min(tt)
  ntot <- colSums(x$x)
  
  fExp <- function(fc,fp,pij){ ##handle multiple clones
    pij*fp/(fc-fp)*(exp(tt*(fc-fp))-1)
  }
  
  logLnn <- function(fc,fp,ni,sdy){
    xc <- colSums(do.call(rbind,lapply(1:length(ni$kp),function(j){
      fExp(fc,fp[j],ni$pij[j])*ni$np[j,]/ntot
    })))
    deltaf <- fc-fp
    negll <- c(-dbinom(ni$nc,size = ntot,prob = xc,log=T), -dnorm(deltaf,mean=0,sd=sdy,log=T))
    negll[!is.finite(negll)] <- 10^9## to prevent warnings when probability is zero
    sum(negll)
    
  }

  boot_nn <- do.call(rbind,pbapply::pblapply(1:nrow(f_boot_samples),function(b){
    fb <- f_boot_samples[b,]
    dfb <- as.matrix(dist(as.numeric(fb)))
    dfb[upper.tri(dfb)] <- dfb[upper.tri(dfb)]*(-1)
    sdy <- sd(dfb[fq_nn])
    optimRange <- range(fb)
    optimRange[1] <- optimRange[1]-sd(fb)
    optimRange[2] <- optimRange[2]+sd(fb)
    
    opt_nn <- sapply(n_info,function(ni){
      fp=fb[ni$kp]
      opt <- optimise(logLnn,interval=optimRange,fp=fp,ni=ni,sdy=sdy)$minimum
    })
    return(opt_nn)
  }))
  
  
  return(boot_nn)
}

proc_sweep_input <- function(yi,ntp=8,tmax=1200){
  yi$x <- yi$x[, as.numeric(colnames(yi$x)) < 1200]
  yi$x <- yi$x[, (ncol(yi$x) - ntp + 1):ncol(yi$x)]
  is.absent <- rowSums(yi$x)==0
  yi$clone.fitness <- yi$clone.fitness[!is.absent]
  yi$x <- yi$x[!is.absent,]
  return(yi)
}


## SO RATHER THAN TRY AND INFER X1 WE ARE JUST GOING TO INTERPOLATE IT FROM THE BOOTSTRAPPED DATA. 
## WE ONLY REALLY NEED IT TO on THE NEAREST NEIGHBOUR ESTIMATION STEP ANYWAY.
## I THINK THAT THE BOOTSTRAPPING OUGHT TO GIVE US THE VARIABILITY WE NEED.
## STILL LOOK FOR WAYS TO AUTOMATE THAT MIN OBS THRESHOLD THOUGH. 
if(FALSE) y <- readRDS("data/proc/sweep_inputs.Rds")

data  <- proc_sweep_input(y[[101]],ntp=8)
fq_boot <- solve_fitness_bootstrap(data, minobs=20, nboot = 100)
#nn_boot <- boot_nn(fq_boot,x)

