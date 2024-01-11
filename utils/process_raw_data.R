setwd("/mnt/ix1/Shared_Folders_crdlab/HighPloidy_CostBenefits/data/timeSeries_singleCellSeq/TableS1")

fo <- list.files()
has_data <- sapply(fo, function(foi) file.exists(paste0(foi,"/named_mat.csv")))
fpaths <- fo[has_data]


x <- lapply(fpaths, function(fo){
  fi <- paste0(fo,"/named_mat.csv")
  x <- data.table::fread(fi,nrows = 2)
  x <- colnames(x)
  x <- sapply(x, function(xi) unlist(strsplit(xi,split="-"))[2])
  x <- unique(x)
  x <- x[!x=="bin_name"]
  y <- rep(fo,length(x))
  names(y) <- x
  return(y)
})


y <- unlist(x)
z <- readxl::read_xlsx("../41586_2021_3648_MOESM3_ESM.xlsx",sheet = 1)
z <- data.frame(z[,1:6])
z$foldername <- NaN
z_in <- z[z$library_id%in%names(y),]
z_out <- z[!z$library_id%in%names(y),]
z_in$foldername <- y[z_in$library_id]
x <- list(z_in=z_in,z_out=z_out)
saveRDS(x,"~/tmp_unique_ids.Rds")
