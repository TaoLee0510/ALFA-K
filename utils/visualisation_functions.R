melt_for_plotting <- function(data_obj, nclones=5, karyotypes=NULL, fit_obj=NULL){
  if(!is.null(fit_obj) && is.null(karyotypes)){
    nclones <- min(nclones, sum(fit_obj$xo$id == "fq"))
  }
  
  dt <- data_obj$dt
  x <- data_obj$x
  
  # If karyotypes are provided, subset by them
  if(!is.null(karyotypes)){
    x <- x[rownames(x) %in% karyotypes, , drop=FALSE]
  } else { 
    # Otherwise, default to selecting top nclones
    x <- x[order(rowSums(x), decreasing=TRUE), , drop=FALSE]
    x <- head(x, nclones)
  }
  
  # Normalize frequencies
  for(i in 1:ncol(x)) x[,i] <- x[,i] / sum(x[,i])
  
  # Reshape data for plotting
  x <- reshape2::melt(x)
  colnames(x) <- c("karyotype", "time", "frequency")
  x$time <- x$time * dt
  
  out <- list(data=x, fit=NULL)
  
  # Handle fit_obj if provided
  if(!is.null(fit_obj)){
    y <- fit_obj$xo
    y <- y[unique(as.character(x$karyotype)), , drop=FALSE]
    
    t <- seq(min(x$time), max(x$time), length.out=100)
    
    z <- y$u0 + y$f_est %*% t(t)
    z <- apply(z, 2, function(zi) exp(zi) / sum(exp(zi)))
    colnames(z) <- t
    rownames(z) <- rownames(y)
    
    z <- reshape2::melt(z)
    colnames(z) <- c("karyotype", "time", "frequency")    
    out$fit <- z
  }
  return(out)
}
