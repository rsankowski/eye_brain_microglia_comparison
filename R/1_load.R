# adjusted from Josip Herman
library(readr)
library(data.table)
library(tidyverse)
library(tools)
library(assertthat)

# Specify directories
file_paths <- file.path("data","counts")
names(file_paths) <- file_paths


# Get complete paths to all files
all_file_paths        <- file.path(file_paths, list.files(file_paths))

names(all_file_paths) <- lapply(strsplit(all_file_paths,split="/"), function(x) { sub(".coutt.csv|.coutt.txt.gz","",x[length(x)]) } )

# Calculate md5sums to check for duplicated
md5sums      <- lapply(all_file_paths, function(x) {md5sum(x)} )

# Check for duplicated data
assert_that(sum(duplicated(md5sums)) == 0)

# Check for duplicated names
assert_that(sum(duplicated(unname(unlist(lapply(strsplit(unlist(files_by_ext),split = "/"),tail,1))))) == 0)

####
#### LOADING
####
# Loading data using lapply
data_list   <- lapply(all_file_paths, function(x) {fread(x, header= T)} )

# Add dataset name prefix to all columns, Merge with remaining gene names
for (d in names(data_list)) { 
  data_list[[d]]$GENEID <- make.unique(gsub("_.*", "", data_list[[d]]$GENEID))
  colnames(data_list[[d]]) <- c("GENEID", paste(d, "_",1:192,sep="" ))
  
}

# Cbind list of data.tables and removing the GENEID column from data.tables
data_list_cbind <- reduce(data_list, full_join, by = "GENEID")
data_list_cbind[is.na(data_list_cbind)]   <- 0

# fuse data
prdata <- as.data.frame(data_list_cbind)
rownames(prdata) <- prdata$GENEID
prdata$GENEID    <- NULL

#remove cells with low expression
cs <- apply(prdata,2,sum)
prdata <- prdata[,cs > 500]

save(prdata, file = file.path("data", "prdata.RData"))
