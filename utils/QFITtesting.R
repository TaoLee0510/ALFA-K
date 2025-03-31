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



# Hybrid Bootstrap Function (Resample Timepoints + Resample Counts)
solve_fitness_bootstrap <- function(data, timepoints, nboot = 1000, epsilon = 1e-6) {
  num_species <- nrow(data)
  num_timepoints <- ncol(data)
  f_boot <- matrix(NA, nrow = nboot, ncol = num_species)
  x1_boot <- matrix(NA, nrow = nboot, ncol = nrow(data))
  
  for (b in 1:nboot) {
    # Step 1: Resample species counts at each timepoint
    boot_data <- bootstrap_counts(data)
    
    # Step 2: Compute new relative abundances (frequencies) per timepoint
    x <- apply(boot_data, 2, function(col) col / sum(col))
    
    # Step 3: Compute finite-difference dx/dt using forward differences
    dx_dt <- compute_dx_dt(x, timepoints)
    
    # Step 4: Trim x to match dx/dt dimensions
    x_trim <- x[, -1,drop=F]  # dimensions: num_species x (num_timepoints - 1)
    
    # Step 5: Resample timepoints (with replacement)
    indices <- sample(1:(num_timepoints - 1), size = (num_timepoints - 1), replace = TRUE)
    
    # Step 6: Construct QP matrices from resampled time intervals
    Q_accum <- matrix(0, nrow = num_species, ncol = num_species)
    r_accum <- rep(0, num_species)
    
    for (t in indices) {
      xt <- x_trim[, t]
      M_t <- diag(xt) - outer(xt, xt)
      Q_accum <- Q_accum + M_t %*% M_t
      r_accum <- r_accum + M_t %*% dx_dt[, t]
    }
    
    # Step 7: Solve the QP problem
    Dmat_boot <- 2 * Q_accum + diag(epsilon, num_species)  # Add ridge term
    dvec_boot <- 2 * r_accum
    
    # Step 8: Constraint: sum(f) = 0 (fix gauge)
    A <- matrix(1, nrow = num_species, ncol = 1)  # sum(f) = 0 constraint
    bvec <- 0
    
    # Solve the QP
    qp_sol <- solve.QP(Dmat_boot, dvec_boot, A, bvec, meq = 1)
    f_boot[b, ] <- qp_sol$solution
  }
  
  return(f_boot)
}


##aggregated version, seems shit. 
boot_nn <- function(f_boot_samples, x, pm = 0.00005) {
  ftot <- rowSums(x$x)
  ntot <- sum(ftot)
  
  ## precalculate for each nn the influx due to missegregations.
  ## Also  extract its observed frequency:
  
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
    
    np <- ftot[nj]
    nc <- 0
    if(ni%in%names(ftot)) nc <- ftot[ni]
    
    list(np=np,nc=nc,kp=nj,pij=pij)
  })
  names(n_info) <- nn_str
  
  logLnn <- function(fc,fp,pp,np,nc){
    nexpect <- sum(np*pp*fp/(fp-fc))
    dpois(nc,nexpect)
  }
  mean_boot <- mean(f_boot_samples)
  sd_boot <- sd(f_boot_samples)
  boot_nn <- do.call(rbind,pbapply::pblapply(1:nrow(f_boot_samples),function(b){
    fb <- f_boot_samples[b,]
    opt_nn <- sapply(n_info,function(ni){
      fp=fb[ni$kp]
      f <- rtruncnorm(10,b=min(fp),mean=mean_boot,sd=sd_boot)
      p <- sapply(f,function(fc){
        logLnn(fc,fp=fp,pp=ni$pij,np=ni$np,nc=ni$nc)
      })
      p[is.na(p)] <- 0
      if(sum(p)==0) return(NaN)
      sample(f,1,prob=p) 
        
    })
    return(opt_nn)
  }))
  
  
  return(boot_nn)
}


get_emer <- function(tt,xfq,fb,xth=1e-6){
  pop_sigma <- numeric(length(tt))
  for (k in seq_along(tt)) {
    # Weighted average of clone fitness by that clone's frequency at time k
    pop_sigma[k] <- sum(xfq[, k] * fb[rownames(xfq)])
  }
  
  # B) Fit sigma(t) = a + b * t
  lm_sigma <- lm(pop_sigma ~ tt)
  a <- coef(lm_sigma)[1]
  b <- coef(lm_sigma)[2]
  
  cat(sprintf("Estimated linear fit for mean fitness: sigma(t) = %.3f + %.3f * t\n", a, b))
  
  
  ## -- 2) HELPER FUNCTION: Solve for T where x(T) = 1e-6 -- ##
  #
  # We use the replicator ODE solution in which x(t) grows as:
  #    x(t) = x_k * exp( [f_i - (a + b t_k)] (t - t_k)
  #                     - 0.5 * b [ t^2 - t_k^2 ] )
  #
  # Rearrange that to find T when x(T) = 1e-6.
  
  emergence_time <- function(x_k, f_i, t_k, a, b) {
    # x_target:
    x_target <- xth
    
    # We want:  x_target = x_k * exp( ... ), so:
    # log(x_target/x_k) = (f_i - [a + b t_k]) (T - t_k) - 0.5 b (T^2 - t_k^2)
    #
    # Let D = log(x_target) - log(x_k).
    # Let R = f_i - (a + b t_k).
    #
    # => D = R (T - t_k) - 0.5 b [ T^2 - t_k^2 ]
    # Rearrange into standard quadratic form: 0.5 b T^2 + ... = 0
    D <- log(x_target) - log(x_k)
    R <- f_i - (a + b * t_k)
    
    # Move all terms to left:
    #  0.5 b T^2 - R T + [ R t_k + D - 0.5 b t_k^2 ] = 0
    A <- 0.5 * b
    B <- -R
    C <- R * t_k + D - 0.5 * b * (t_k^2)
    
    disc <- B^2 - 4 * A * C
    if (disc < 0) {
      # No real solution. Return NA
      return(NA_real_)
    }
    
    # Two solutions
    sol1 <- (-B + sqrt(disc)) / (2 * A)
    sol2 <- (-B - sqrt(disc)) / (2 * A)
    
    # Usually we want the EARLIER time => min(sol1, sol2)
    # (This can be negative if the clone was at 1e-6 before t=0.)
    T_est <- min(sol1, sol2)
    
    return(T_est)
  }
  
  
  ## -- 3) LOOP OVER CLONES, ESTIMATE EMERGENCE TIME -- ##
  #
  # For each clone, we might pick the *last timepoint* (k = ncol(xfq)) 
  # to define x_k = xfq[i, k] and t_k = tt[k], and solve for T.
  # Or pick whichever timepoint the clone has the largest or best-measured frequency.
  
  res <- data.frame(clone = rownames(xfq),
                    fitness = fb[rownames(xfq)],
                    emergence_time = NA_real_)
  
  for (i in seq_len(nrow(xfq))) {
    clone_id <- rownames(xfq)[i]
    f_i <- fb[clone_id]
    
    # e.g. choose the last timepoint's frequency
    k <- ncol(xfq)
    x_k <- xfq[i, k]*res$emergence_time[i] <- T_est
  }
  
  ## -- 4) Inspect result
  head(res)
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




## SO RATHER THAN TRY AND INFER X1 WE ARE JUST GOING TO INTERPOLATE IT FROM THE BOOTSTRAPPED DATA. 
## WE ONLY REALLY NEED IT TO on THE NEAREST NEIGHBOUR ESTIMATION STEP ANYWAY.
## I THINK THAT THE BOOTSTRAPPING OUGHT TO GIVE US THE VARIABILITY WE NEED.
## STILL LOOK FOR WAYS TO AUTOMATE THAT MIN OBS THRESHOLD THOUGH. 
y <- readRDS("data/proc/sweep_inputs.Rds")

yi <- y[[101]]
data <- yi$x

ntp <- 8
minobs <- 20
data <- data[, as.numeric(colnames(data)) < 1200]
data <- data[, (ncol(data) - ntp + 1):ncol(data)]

f_tru <- yi$clone.fitness[order(rowSums(data), decreasing = TRUE)]
data <- data[order(rowSums(data), decreasing = TRUE), ]

f_tru <- f_tru[rowSums(data)>0]
data <- data[rowSums(data)>0,]

dt <- yi$dt

yi$x <- data
yi$clone.fitness <- f_tru

x <- yi

f_tru <- f_tru[rowSums(data) > minobs]
data <- data[rowSums(data) > minobs, , drop = FALSE]

if (nrow(data) < 2) return(NaN)

timepoints <- as.numeric(colnames(data)) * dt
karyotypes <- row.names(data)

data <- matrix(c(data), ncol = ncol(data), byrow = FALSE)

num_species <- nrow(data)
num_timepoints <- ncol(data)

nboot <- 100
fq_boot <- solve_fitness_bootstrap(data, timepoints, nboot = nboot)+0.5
colnames(fq_boot) <- karyotypes
nn_boot <- boot_nn(fq_boot,x)

library(fields)

fit <- Krig(k,f_est,m=1)
set.panel( 2,2) 
plot(fit)

f_est <- c(colMeans(fq_boot,na.rm = T),colMeans(nn_boot,na.rm=T))
k <- do.call(rbind,lapply(names(f_est),s2v))

tmp <- pbapply::pblapply(1:100,function(i) fit <- Krig(k,f_est,m=1))
pl





se_est <- apply(fq_boot,2,function(i) sd(i)/sqrt(length(i)))

df_pred <- data.frame(tru=f_tru,est=f_est,
                      stderr=apply(fq_boot,2,function(i) sd(i)/sqrt(length(i))),
                      sd = apply(fq_boot,2,sd))

p <- ggplot(df_pred,aes(x=tru,y=est))+
  geom_point()+
  geom_errorbar(aes(ymin=est-stderr,ymax=est+stderr))
p


fest <- c(colMeans(nn_boot,na.rm = T))
fest <- fest[names(fest)%in%rownames(x$x)]
ftru <- x$clone.fitness
names(ftru) <- rownames(x$x)
ftru <- ftru[names(fest)]
plot(ftru,fest)

fest <- c(colMeans(fq_boot,na.rm=T))
fest <- fest[names(fest)%in%rownames(x$x)]
ftru <- x$clone.fitness
names(ftru) <- rownames(x$x)
ftru <- ftru[names(fest)]
plot(ftru,fest)
