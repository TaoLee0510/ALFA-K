##########################################
# Helper functions (global environment)
##########################################

pij <- function(i, j, beta) {
  qij <- 0
  if (abs(i - j) > i) { ## not enough copies for i->j
    return(qij)
  }
  if (j == 0) j <- 2 * i
  s <- seq(abs(i - j), i, by = 2)
  for (z in s) {
    qij <- qij + choose(i, z) * beta^z * (1 - beta)^(i - z) *
      0.5^z * choose(z, (z + i - j) / 2)
  }
  return(qij)
}

s2v <- function(s) as.numeric(unlist(strsplit(s, split = "[.]")))

R2R <- function(obs, pred) {
  obs <- obs - mean(obs)
  pred <- pred - mean(pred)
  1 - sum((pred - obs)^2) / sum((obs - mean(obs))^2)
}

gen_all_neighbours <- function(ids, as.strings = TRUE, remove_nullisomes = TRUE) {
  if (as.strings) 
    ids <- lapply(ids, function(ii) as.numeric(unlist(strsplit(ii, split = "[.]"))))
  nkern <- do.call(rbind, lapply(1:length(ids[[1]]), function(i) {
    x0 <- rep(0, length(ids[[1]]))
    x1 <- x0
    x0[i] <- -1
    x1[i] <- 1
    rbind(x0, x1)
  }))
  n <- do.call(rbind, lapply(ids, function(ii) t(apply(nkern, 1, function(i) i + ii))))
  n <- unique(n)
  nids <- length(ids)
  n <- rbind(do.call(rbind, ids), n)
  n <- unique(n)
  n <- n[-(1:nids), ]
  if (remove_nullisomes) 
    n <- n[apply(n, 1, function(ni) sum(ni < 1) == 0), ]
  n
}

bootstrap_counts <- function(data) {
  num_species <- nrow(data)
  num_timepoints <- ncol(data)
  boot_data <- matrix(NA, nrow = num_species, ncol = num_timepoints)
  rownames(boot_data) <- rownames(data)
  for (i in seq_len(num_timepoints)) {
    total_counts <- sum(data[, i])
    if (total_counts == 0) {
      boot_data[, i] <- rep(0, num_species)
    } else {
      boot_data[, i] <- rmultinom(1, size = total_counts, prob = data[, i] / total_counts)
    }
  }
  boot_data
}

compute_dx_dt <- function(x, timepoints) {
  num_species <- nrow(x)
  dx_dt <- matrix(NA, nrow = num_species, ncol = ncol(x) - 1)
  for (i in seq_len(num_species)) {
    dx_dt[i, ] <- diff(x[i, ]) / diff(timepoints)
  }
  dx_dt
}

logSumExp <- function(v) {
  m <- max(v)
  m + log(sum(exp(v - m)))
}

gen_nn_info <- function(fq, pm = 0.00005) {
  nn <- gen_all_neighbours(fq)
  nn_str <- as.character(apply(nn, 1, paste, collapse = "."))
  n_info <- lapply(nn_str, function(ni) {
    nj <- as.character(apply(gen_all_neighbours(ni), 1, paste, collapse = "."))
    nj <- nj[nj %in% fq]
    nivec <- s2v(ni)
    pij_vals <- sapply(nj, function(si) {
      si_vec <- as.numeric(unlist(strsplit(si, split = "[.]")))
      prod(sapply(1:length(si_vec), function(k) pij(si_vec[k], nivec[k], pm)))
    })
    list(ni = ni, nj = nj, pij = pij_vals)
  })
  n_info
}

neg_log_lik <- function(param, counts, timepoints) {
  K <- nrow(counts)
  Tt <- ncol(counts)
  f_free <- param[1:(K - 1)]
  f_full <- c(f_free, -sum(f_free))
  log_x0 <- param[K:(2 * K - 1)]
  nll <- 0
  for (i in seq_len(Tt)) {
    lv <- log_x0 + f_full * timepoints[i]
    denom <- logSumExp(lv)
    for (k in seq_len(K)) {
      if (counts[k, i] > 0) {
        nll <- nll - counts[k, i] * (lv[k] - denom)
      }
    }
  }
  nll
}

joint_optimize <- function(counts, timepoints, f_init, x0_init) {
  K <- length(f_init)
  f_free_init <- f_init[1:(K - 1)]
  x0_init_log <- log(x0_init + 1e-12)
  param_init <- c(f_free_init, x0_init_log)
  obj_fun <- function(par) neg_log_lik(par, counts, timepoints)
  opt <- optim(par = param_init, fn = obj_fun, method = "BFGS",
               control = list(maxit = 200, reltol = 1e-8))
  f_free_opt <- opt$par[1:(K - 1)]
  f_opt <- c(f_free_opt, -sum(f_free_opt))
  log_x0_opt <- opt$par[K:(2 * K - 1)]
  x0_opt <- exp(log_x0_opt)
  x0_opt <- x0_opt / sum(x0_opt)
  list(f = f_opt, x0 = x0_opt)
}

proc_sweep_input <- function(yi, ntp = 8, tmax = 1200) {
  yi$x <- yi$x[, as.numeric(colnames(yi$x)) < tmax]
  yi$x <- yi$x[, (ncol(yi$x) - ntp + 1):ncol(yi$x)]
  is.absent <- rowSums(yi$x) == 0
  yi$clone.fitness <- yi$clone.fitness[!is.absent]
  yi$x <- yi$x[!is.absent, ]
  names(yi$clone.fitness) <- row.names(yi$x)
  yi
}

project_forward_log <- function(x0, f, timepoints) {
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

optimize_initial_frequencies <- function(x_obs, f, timepoints) {
  loss_function <- function(log_x0) {
    x0 <- exp(log_x0)
    x0 <- x0 / sum(x0)
    x_pred <- project_forward_log(x0, f, timepoints)
    sum((x_pred - x_obs)^2)
  }
  x_ini <- x_obs[, 1] + 1e-6
  x_ini <- x_ini / sum(x_ini)
  opt_result <- optim(log(x_ini), loss_function, method = "BFGS")
  x0_opt <- exp(opt_result$par)
  x0_opt / sum(x0_opt)
}

find_birth_times <- function(opt_res, time_range, minF) {
  f_est <- opt_res$f
  x0_est <- opt_res$x0
  num_species <- length(f_est)
  birth_times <- rep(NA, num_species)
  for (i in seq_len(num_species)) {
    if (f_est[i] <= min(f_est)) next
    birth_fn <- function(t) {
      log_x_t <- log(x0_est[i]) + f_est[i] * t
      denom <- logSumExp(log(x0_est) + f_est * t)
      exp(log_x_t - denom) - minF
    }
    root <- try(uniroot(birth_fn, range(time_range), tol = 1e-6)$root, silent = TRUE)
    if (!inherits(root, "try-error")) {
      birth_times[i] <- root
    }
  }
  birth_times
}

##########################################
# Core functions that use parallel loops
##########################################

solve_fitness_bootstrap <- function(data, minobs, nboot = 1000, epsilon = 1e-6, pm = 0.00005,
                                    n0, nb, passage_times = NULL, cl = NULL) {
  fq <- rownames(data$x)[rowSums(data$x) > minobs]
  nn <- gen_nn_info(fq, pm)
  names(nn) <- sapply(nn, function(nni) nni$ni)
  fq_vec <- do.call(rbind, lapply(fq, s2v))
  fq_nn <- which(as.matrix(dist(fq_vec)) == 1)
  timepoints <- as.numeric(colnames(data$x)) * data$dt
  num_species <- length(fq)
  num_timepoints <- ncol(data$x)
  
  # Define bootstrap_iter with explicit arguments
  bootstrap_iter <- function(b, data, fq, timepoints, num_species, num_timepoints,
                             epsilon, n0, nb, passage_times, nn) {
    boot_data <- bootstrap_counts(data$x)
    x <- apply(boot_data[fq, , drop = FALSE], 2, function(col) col / sum(col))
    dx_dt <- compute_dx_dt(x, timepoints)
    x_trim <- x[, -1, drop = FALSE]
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
    x0_init <- optimize_initial_frequencies(x, f_qp, timepoints)
    opt_res <- joint_optimize(boot_data[fq, , drop = FALSE], timepoints, f_qp, x0_init)
    
    g0 <- log(nb / n0) / diff(timepoints)[1]
    if (!is.null(passage_times))
      g0 <- log(nb / n0) / diff(passage_times * data$dt)[1]
    opt_res$f <- opt_res$f + g0 - sum(opt_res$x0 * opt_res$f)
    birth_times <- find_birth_times(opt_res, time_range = c(-1000, max(timepoints)), minF = 1 / n0)
    peak_times <- timepoints[apply(x, 1, which.max)]
    mean_risetime <- mean(peak_times - birth_times, na.rm = TRUE)
    birth_times[is.na(birth_times)] <- peak_times[is.na(birth_times)] - mean_risetime
    
    x0par <- opt_res$x0
    names(x0par) <- fq
    fpar <- opt_res$f
    names(fpar) <- fq
    names(birth_times) <- fq
    
    dfb <- as.matrix(dist(as.numeric(fpar)))
    dfb[upper.tri(dfb)] <- dfb[upper.tri(dfb)] * (-1)
    f_final <- opt_res$f
    x0_final <- opt_res$x0
    
    fExp <- function(fc, fp, pij_val, tt) {
      pij_val * fp / (fc - fp) * (exp(tt * (fc - fp)) - 1)
    }
    xfit <- project_forward_log(x0par, fpar, timepoints)
    rownames(xfit) <- fq
    ntot <- colSums(boot_data)
    
    opt_fc <- function(fc, nni, prior_mean = NULL, prior_sd = NULL, do_prior = FALSE) {
      child <- nni$ni
      xc_est <- colSums(do.call(rbind, lapply(1:length(nni$nj), function(i) {
        par <- nni$nj[i]
        tt <- timepoints - birth_times[par]
        fExp(fc, fpar[par], nni$pij[i], tt) * xfit[par, ]
      })))
      xc_est <- pmax(0, pmin(1, xc_est))
      xc_obs <- rep(0, length(timepoints))
      if (child %in% rownames(boot_data))
        xc_obs <- boot_data[nni$ni, ]
      res <- dbinom(xc_obs, ntot, prob = xc_est, log = TRUE)
      if (do_prior)
        res <- c(res, dnorm(fc - fpar[nni$nj], mean = prior_mean, sd = prior_sd, log = TRUE))
      res[!is.finite(res)] <- -(10^9)
      -sum(res)
    }
    search_interval <- range(fpar)
    interval_range <- diff(search_interval)
    search_interval[1] <- search_interval[1] - interval_range
    search_interval[2] <- search_interval[2] + interval_range
    nn_present <- names(nn) %in% rownames(boot_data)
    nn_present[nn_present] <- nn_present[nn_present] & sapply(names(nn)[nn_present], function(ni) {
      sum(boot_data[ni, ]) > 0
    })
    fc <- rep(NaN, length(nn))
    if (any(nn_present)) {
      fc[nn_present] <- sapply(nn[nn_present], function(nni) {
        res <- optimise(opt_fc, interval = search_interval, nni = nni, do_prior = FALSE)
        res$minimum
      })
    }
    fc_prior <- unlist(lapply(1:length(fc[nn_present]), function(i) {
      as.numeric(fpar[nn[nn_present][[i]]$nj] - fc[nn_present][i])
    }))
    if (any(!nn_present)) {
      fc[!nn_present] <- sapply(nn[!nn_present], function(nni) {
        res <- optimise(opt_fc, interval = search_interval, nni = nni,
                        prior_mean = mean(fc_prior),
                        prior_sd = sd(fc_prior),
                        do_prior = TRUE)
        res$minimum
      })
    }
    list(f_initial = f_qp,
         f_final = f_final,
         x0_initial = x0_init,
         x0_final = x0_final,
         f_nn = fc)
  }
  
  # Run bootstrap iterations in parallel if cl is provided; otherwise, run serially.
  boot_list <- if (is.null(cl)) {
    lapply(1:nboot, function(b) {
      bootstrap_iter(b, data, fq, timepoints, num_species, num_timepoints,
                     epsilon, n0, nb, passage_times, nn)
    })
  } else {
    parLapplyLB(cl, 1:nboot, bootstrap_iter, data = data, fq = fq, timepoints = timepoints,
                num_species = num_species, num_timepoints = num_timepoints,
                epsilon = epsilon, n0 = n0, nb = nb, passage_times = passage_times, nn = nn)
  }
  
  f_initial_mat <- do.call(rbind, lapply(boot_list, function(x) x$f_initial))
  f_final_mat   <- do.call(rbind, lapply(boot_list, function(x) x$f_final))
  x0_initial_mat <- do.call(rbind, lapply(boot_list, function(x) x$x0_initial))
  x0_final_mat  <- do.call(rbind, lapply(boot_list, function(x) x$x0_final))
  f_nn_mat <- do.call(rbind, lapply(boot_list, function(x) x$f_nn))
  
  colnames(f_initial_mat) <- fq
  colnames(f_nn_mat) <- names(nn)
  colnames(f_final_mat) <- fq
  colnames(x0_final_mat) <- fq
  colnames(x0_initial_mat) <- fq
  
  list(initial_fitness = f_initial_mat, 
       final_fitness = f_final_mat,
       initial_frequencies = x0_initial_mat,
       final_frequencies = x0_final_mat,
       nn_fitness = f_nn_mat)
}

gen_coords <- function(k, cl = NULL) {
  indices <- unlist(lapply(1:22, rep, 2))
  res <- if (is.null(cl)) {
    do.call(rbind, lapply(1:length(k), function(i) {
      k0 <- s2v(k[i])
      nn <- apply(gen_all_neighbours(k[i], remove_nullisomes = FALSE), 1, paste, collapse = ".")
      nc <- k0[indices]
      names(nc) <- nn
      ii <- which(k %in% nn)
      nc <- as.numeric(nc[k[ii]])
      ii <- c(i, ii)
      jj <- rep(i, length(ii))
      nc <- c(sum(k0), nc)
      cbind(ii, jj, nc)
    }))
  } else {
    do.call(rbind, parLapplyLB(cl, 1:length(k), function(i) {
      k0 <- s2v(k[i])
      nn <- apply(gen_all_neighbours(k[i], remove_nullisomes = FALSE), 1, paste, collapse = ".")
      nc <- k0[indices]
      names(nc) <- nn
      ii <- which(k %in% nn)
      nc <- as.numeric(nc[k[ii]])
      ii <- c(i, ii)
      jj <- rep(i, length(ii))
      nc <- c(sum(k0), nc)
      cbind(ii, jj, nc)
    }))
  }
  res
}

fitKrig <- function(fq_boot, nboot, cl = NULL) {
  fboot <- cbind(fq_boot$final_fitness, fq_boot$nn_fitness)
  fq_str <- colnames(fq_boot$final_fitness)
  nn_str <- colnames(fq_boot$nn_fitness)
  ktrain <- do.call(rbind, lapply(c(fq_str, nn_str), s2v))
  ktest <- unique(rbind(ktest = ktrain, gen_all_neighbours(nn_str)))
  ktest_str <- apply(ktest, 1, paste, collapse = ".")
  fq_ids <- ktest_str %in% fq_str
  nn_ids <- ktest_str %in% nn_str
  boot_predictions_list <- if (is.null(cl)) {
    lapply(1:nboot, function(b) {
      boot_f <- fboot[cbind(sample(1:nrow(fboot), ncol(fboot), replace = TRUE), 1:ncol(fboot))]
      fit_boot <- Krig(ktrain, boot_f,
                       cov.function = "stationary.cov",
                       cov.args = list(Covariance = "Matern", smoothness = 1.5))
      predict(fit_boot, ktest)
    })
  } else {
    parLapplyLB(cl, 1:nboot, function(b) {
      boot_f <- fboot[cbind(sample(1:nrow(fboot), ncol(fboot), replace = TRUE), 1:ncol(fboot))]
      fit_boot <- Krig(ktrain, boot_f,
                       cov.function = "stationary.cov",
                       cov.args = list(Covariance = "Matern", smoothness = 1.5))
      predict(fit_boot, ktest)
    })
  }
  boot_predictions <- do.call(cbind, boot_predictions_list)
  pred_means <- apply(boot_predictions, 1, mean)
  pred_medians <- apply(boot_predictions, 1, median)
  pred_sd <- apply(boot_predictions, 1, sd)
  data.frame(k = ktest_str, mean = pred_means, median = pred_medians, sd = pred_sd,
             fq = fq_ids, nn = nn_ids)
}

xval <- function(fq_boot, cl = NULL) {
  fboot <- cbind(fq_boot$final_fitness, fq_boot$nn_fitness)
  fq_str <- colnames(fq_boot$final_fitness)
  nn_str <- colnames(fq_boot$nn_fitness)
  ktrain <- do.call(rbind, lapply(c(fq_str, nn_str), s2v))
  ids <- unlist(lapply(1:length(fq_str), function(i) {
    ki <- c(fq_str[i], as.character(apply(gen_all_neighbours(fq_str[i]), 1, paste, collapse = ".")))
    idi <- rep(i, length(ki))
    names(idi) <- ki
    idi
  }))
  ids <- ids[!duplicated(names(ids))]
  uids <- unique(ids)
  tmp_list <- if (is.null(cl)) {
    lapply(uids, function(id) {
      fboot_shuffled <- apply(fboot, 2, sample)
      fi <- fboot_shuffled[1, ]
      train_k <- ktrain[!(ids == id), ]
      train_f <- fi[!(ids == id)]
      test_k <- ktrain[(ids == id), ]
      test_f <- fi[(ids == id)]
      fit <- Krig(train_k, train_f,
                  cov.function = "stationary.cov",
                  cov.args = list(Covariance = "Matern", smoothness = 1.5))
      est_f <- predict(fit, test_k)
      cbind(test_f, est_f)
    })
  } else {
    parLapplyLB(cl, uids, function(id) {
      fboot_shuffled <- apply(fboot, 2, sample)
      fi <- fboot_shuffled[1, ]
      train_k <- ktrain[!(ids == id), ]
      train_f <- fi[!(ids == id)]
      test_k <- ktrain[(ids == id), ]
      test_f <- fi[(ids == id)]
      fit <- Krig(train_k, train_f,
                  cov.function = "stationary.cov",
                  cov.args = list(Covariance = "Matern", smoothness = 1.5))
      est_f <- predict(fit, test_k)
      cbind(test_f, est_f)
    })
  }
  tmp <- do.call(rbind, tmp_list)
  R2R(tmp[, 1], tmp[, 2])
}

predict_evo <- function(landscape, pred_times, dtest, tm, t, pred_iters, cl = NULL) {
  single_iter <- function(iter) {
    iter_xpred <- matrix(0, nrow = length(pred_times), ncol = nrow(landscape))
    xp <- rep(0, nrow(landscape))
    names(xp) <- landscape$k
    x0 <- rmultinom(1, sum(dtest), dtest / sum(dtest))[, 1]
    x0 <- x0 / sum(x0)
    xp[names(x0)] <- x0
    fb <- rnorm(nrow(landscape), mean = landscape$mean, sd = landscape$sd)
    for (time in t[-1]) {
      mean_fitness <- sum(xp * fb)
      M <- Diagonal(x = fb) - mean_fitness * Diagonal(length(fb)) + tm
      xp <- M %*% xp
      if (time %in% pred_times) {
        idx <- which(time == pred_times)
        iter_xpred[idx, ] <- iter_xpred[idx, ] + as.numeric(xp / sum(xp))
      }
    }
    iter_xpred
  }
  iter_results <- if (is.null(cl)) {
    lapply(1:pred_iters, single_iter)
  } else {
    parLapplyLB(cl, 1:pred_iters, single_iter)
  }
  xpred <- Reduce("+", iter_results)
  for (i in 1:nrow(xpred)) xpred[i, ] <- xpred[i, ] / sum(xpred[i, ])
  rownames(xpred) <- pred_times
  colnames(xpred) <- landscape$k
  xpred
}

tmbuild <- function(p0, coords, dims) {
  xx <- coords[, "nc"] * p0 / 2
  sources <- coords[, "ii"] == coords[, "jj"]
  xx[sources] <- (1 - coords[sources, "nc"] * p0)
  sparseMatrix(i = coords[, "ii"],
               j = coords[, "jj"],
               x = xx,
               dims = dims)
}

##########################################
# Main function: alfak
##########################################

alfak <- function(yi, outdir, passage_times, minobs = 20,
                  nboot = 45,
                  pred_iters = 100,
                  n0 = 1e5,
                  nb = 1e7,
                  pred_times = c(275),
                  pm = 0.00005,
                  num_cores = NULL) {
  
  library(quadprog)
  library(fields)
  library(Matrix)
  library(parallel)
  
  dir.create(outdir, recursive = TRUE)
  
  if (max(rowSums(yi$x)) < minobs)
    stop(paste("no frequent karyotypes detected for minobs", minobs))
  
  # Set up a cluster and export all needed objects robustly
  cl <- NULL
  if (!is.null(num_cores)) {
    cl <- makeCluster(num_cores)
    clusterEvalQ(cl, {
      library(quadprog)
      library(fields)
      library(Matrix)
    })
    # Export all global functions and variables used in parallel workers.
    # (This list is explicit to ensure nothing is missing.)
    clusterExport(cl, varlist = c("pij", "s2v", "R2R", "gen_all_neighbours",
                                  "bootstrap_counts", "compute_dx_dt", "logSumExp",
                                  "gen_nn_info", "neg_log_lik", "joint_optimize",
                                  "optimize_initial_frequencies", "find_birth_times",
                                  "project_forward_log", "solve.QP",
                                  "gen_coords", "fitKrig", "xval", "predict_evo",
                                  "tmbuild"), envir = .GlobalEnv)
    # Also export parameters from this function's environment that will be used
    clusterExport(cl, varlist = c("passage_times"), envir = environment())
  }
  
  fq_boot <- solve_fitness_bootstrap(yi, minobs = minobs, nboot = nboot,
                                     n0 = n0, nb = nb, pm = pm,
                                     passage_times = passage_times, cl = cl)
  saveRDS(fq_boot, file = file.path(outdir, "bootstrap_res.Rds"))
  
  landscape <- fitKrig(fq_boot, nboot, cl = cl)
  saveRDS(landscape, file = file.path(outdir, "landscape.Rds"))
  
  Rxv <- xval(fq_boot, cl = cl)
  saveRDS(Rxv, file = file.path(outdir, "xval.Rds"))
  
  if (length(pred_times) == 0) {
    if (!is.null(cl)) stopCluster(cl)
    return(Rxv)
  }
  
  coords <- gen_coords(landscape$k, cl = cl)
  dims <- rep(nrow(landscape), 2)
  tm <- tmbuild(pm, coords, dims)
  
  dtest <- yi$x[rownames(yi$x) %in% landscape$k, ncol(yi$x)]
  
  t0 <- max(as.numeric(colnames(yi$x)) * yi$dt)
  tend <- max(pred_times)
  t_seq <- t0:tend
  
  xpred <- predict_evo(landscape, pred_times, dtest, tm, t_seq, pred_iters, cl = cl)
  saveRDS(xpred, file = file.path(outdir, "predictions.Rds"))
  
  if (!is.null(cl)) stopCluster(cl)
  return(Rxv)
}




