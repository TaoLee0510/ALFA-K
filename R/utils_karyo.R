## helper function convert string karyotype 2.2.2.2 -> numeric
s2v <- function(charvec) as.numeric(unlist(strsplit(charvec,split="[.]")))

## generates all one missegregation neighbours of inputs karyotypes.
gen_all_neighbours <- function(ids, as.strings = TRUE, remove_nullisomes = TRUE) {
  if (as.strings) 
    ids <- lapply(ids, s2v)
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

## returns fitness value using one of the GRF random landscapes as input (tru_lscape)
getf <- function(k,wavelength,tru_lscape){
  Nwaves <- nrow(tru_lscape)
  scalef <- 1/(pi*sqrt(Nwaves))
  d <- apply(tru_lscape,1,function(ci){
    sqrt(sum((k-ci)^2))
  })
  sum(sin(d/wavelength)*scalef)
}