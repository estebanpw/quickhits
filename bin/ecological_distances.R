
# Global
setwd("C:/Users/Annette/Documents/Cosas de esteban/metag_vs_metag")


filenames <- list.files("C:/Users/Annette/Documents/Cosas de esteban/metag_vs_metag/quickhits_results_braz", pattern="*.data", full.names=TRUE)

# Table to store the *tables*
tables <- list()
n_files <- length(filenames)

# Read all files
for(i in 1:n_files){
  # Read individually
  tables[[i]] <- scan(filenames[i], skip=1, sep="\n", what = character())  
  
}

# They must have same length
rows <- length(tables[[1]])
cols <- length(strsplit(tables[[1]][1], split = "\t")[[1]])

# Matrix with the column of interest values
data_hits <- array(NA, dim = c(n_files, rows))

for(i in 1:n_files){
  for(j in 1:rows){
    
    cada_campo <- strsplit(tables[[i]][j], split = "\t")[[1]][2]
    data_hits[i,j] <- as.numeric(cada_campo)
  }
}

## Row and col names
all_names <- c()
for(i in 1:n_files){
  name <- strsplit(filenames[i], "/")
  val <- name[[1]][length(name[[1]])]
  all_names[i] <- strsplit(val, "[.]")[[1]][1]
}


# Euclidean Heatmap matrix
euc_dists <- array(NA, dim = c(n_files, n_files))
for(i in 1:n_files){
  for(j in 1:n_files){
    euc_dists[i,j] <- euclid(data_hits[i,], data_hits[j,], rows)
  }
}

# Bray Curtis distance
braycurtis_dists <- array(NA, dim = c(n_files, n_files))
for(i in 1:n_files){
  for(j in 1:n_files){
    braycurtis_dists[i,j] <- braycurtis(data_hits[i,], data_hits[j,], rows)
  }
}

# Plot heatmaps
#heatmap(euc_dists, col= heat.colors(256), labRow = all_names, labCol = all_names)
heatmap(braycurtis_dists, col= heat.colors(256), labRow = all_names, labCol = all_names)


####### Compute euclidean distance
euclid <- function(col_1, col_2, size){
  sum_all <- 0
  for(i in 1:size){
    sum_all <- sum_all + ((col_1[i] - col_2[i])*(col_1[i] - col_2[i]))
  }
  return (sqrt(sum_all))
}

braycurtis <- function(col_1, col_2, size){
  minsum <- 0
  sumsum <- 0
  for(i in 1:size){
    minsum <- minsum + min(col_1[i], col_2[i])
    sumsum <- sumsum + (col_1[i] + col_2[i])
  }
  return (1 - 2 * (minsum/sumsum))
}



