#################
# Title: Extra Genotype Mass Univariate Analysis
# Author: Andrew DiLernia
# Date: 1/12/2018
# Purpose: Conduct mass univariate analysis for new genotype data
#################

library(data.table)
library(readxl)
library(R.utils)
library(tidyverse)
library(stringr)
library(tidyr)
library(dplyr)
library("foreach")
library("doParallel")

dir <- "/panfs/roc/groups/15/petersla/diler001/Genetics/"
#dir <- "C:/Users/Owner/Google Drive/Genetics_Research/"

path <- dir

setwd(dir)
population <- "CEU"

cleaned <- FALSE

#####Initial downloading of data - only need to do this section once#####

if(file.exists("CEU_Gen_Dat_Check.csv") == FALSE) {
  load("file_names.Rda")
  
  file_names <- c(file_names[which(grepl("CEU", file_names, fixed=TRUE))],
                  file_names[which(grepl("YRI", file_names, fixed=TRUE))])
  
  URL <- "ftp://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/2010-08_phaseII+III/forward/"
  
  URLs = paste0(URL, file_names)
  
  Gen_Get <- function(index, urls = URLs, pop = "CEU")
  {
    files <- file_names[which(grepl(pop, file_names, fixed=TRUE))]
    
    if(file.exists(substr(x = files[[index]], start=0, stop = (nchar(files[[index]]) - 3))) == FALSE) {
      download.file(url = urls[[index]], destfile = files[[index]])
      gunzip(files[[index]], overwrite=TRUE)
    }
    
    fread(substr(x = files[[index]], start=0, stop = (nchar(files[[index]]) - 3)), fill = TRUE)
  }
  
  CEU_Gen_Dat <- map_dfr(.f = Gen_Get, .x = as.list(1:25),
                         pop = "CEU")
  YRI_Gen_Dat <- map_dfr(.f = Gen_Get, .x = as.list(1:25),
                         pop = "YRI")
  
  # Subsetting by columns
  CEU_Gen_Dat <- setDT(CEU_Gen_Dat)[, c("chrom", "pos", "strand", "assembly#", "center",
                                        "protLSID", "assayLSID", "panelLSID", "QCcode") := NULL]
  YRI_Gen_Dat <- setDT(YRI_Gen_Dat)[, c("chrom", "pos", "strand", "assembly#", "center",
                                        "protLSID", "assayLSID", "panelLSID", "QCcode") := NULL]
  
  fwrite(CEU_Gen_Dat, file = "CEU_Gen_Dat_Check.csv")
  fwrite(YRI_Gen_Dat, file = "YRI_Gen_Dat_Check.csv")
}

#####End initial downloading of data#####

#####Loading old SNP genotypes#####

if(file.exists(paste0("snp_matrix.", population, "_Check.csv")) == FALSE) {
  SNP_data <- fread(paste0(path,'snp_matrix.', population, '.csv'))
  
  # Cleans up snp names
  colnames(SNP_data) <- colnames(SNP_data) %>%
    gsub(pattern = " ", replacement = "") %>% tolower()
  
  # Only keeps columns with valid snp names
  SNP_data <- SNP_data %>% as.data.frame()
  SNP_data <- SNP_data[, c(1, which(str_detect(colnames(SNP_data), "rs")))]
  
  colnames(SNP_data)[1] <- "Cell_line"
  
  Cell_lines <- SNP_data$Cell_line %>% gsub(pattern = " ", replacement = "") %>% toupper()
  
  # Adds back in Cell_line column
  SNP_data$Cell_line <- unlist(Cell_lines)
  print(paste0("Number of cell lines: ", length(Cell_lines)))
  
  # Saves original data
  fwrite(SNP_data, file = paste0("snp_matrix.", population, "_Check.csv"))
}

# Obtaining vector of population cell lines
load(paste0("Original_Data/", population, "_IC_data.Rda"))
Cell_lines <- get(paste0(population, "_IC_data"))$Cell_line
remove(list = paste0(population, "_IC_data"))

# Obtaining vector of original snps
SNP_data <- fread(paste0("snp_matrix.", population, "_Check.csv"))
remove("SNP_data")
print(paste0("Number of cell lines: ", length(Cell_lines)))

if(cleaned == FALSE) {
  
  #####Cleaning new genotype data#####
  
  numRows <- nrow(fread(paste0(population, "_Gen_Dat_Check.csv"), select = 1))
  parts <- 100
  chunkSize <- ceiling(numRows/parts)
  
  stopStart <- data.frame(start = seq(1, numRows, by = chunkSize), 
                          stop = c((seq(1, numRows, by = chunkSize)-1)[-1], numRows))

  library(foreach)
  library(doParallel)
  ncores <- 10 # set ncores to something less than or equal to available cores
  cl <- makeCluster(ncores) # creates a cluster with <ncore> cores
  registerDoParallel(cl) # register the cluster
  
  res <- foreach(chunk = 1:nrow(stopStart)) %dopar% {
    library(data.table)
    library(readxl)
    library(tidyverse)
  
  # Loading partition of new genotype data
  Gen_Dat <- fread(paste0(population, "_Gen_Dat_Check.csv"), skip = (chunk-1)*chunkSize, nrows = chunkSize)
  
  # Loading snps to be excluded due to potential errors
  SNPs_to_Exclude <- read_excel(paste0(path, "SNPs_to_Exclude.xlsx"))$rsid %>% unique() %>% 
    gsub(pattern = " ", replacement = "") %>% tolower()
  
  # Cleaning up snp names
  colnames(Gen_Dat)[which(colnames(Gen_Dat) == "rs#")] <- "SNP" 
  Gen_Dat$SNP <- Gen_Dat$SNP %>% gsub(pattern = " ", replacement = "") %>% tolower()
  
  # Cleaning up cell line names
  colnames(Gen_Dat) <- colnames(Gen_Dat) %>% gsub(pattern = " ", replacement = "") %>% toupper()
  
  # Subsetting to only keep necessary cell line columns
  Cell_lines <- Cell_lines %>% unlist() %>% trimws() %>% toupper()
  Gen_Dat <- Gen_Dat[, c(1:2, which(colnames(Gen_Dat) %in% Cell_lines)), with = FALSE]
  
  # Determining which snps are not in error list of snps
  # where errors are from: ftp://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/2010-08_phaseII+III/inconsistencies/
  SNPs_ok <- Gen_Dat$SNP[which(!(Gen_Dat$SNP %in% c(SNPs_to_Exclude)))] %>% unique() 
  
  # Excludes error snps from: ftp://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/2010-08_phaseII+III/inconsistencies/
  Gen_Dat <- as.data.frame(Gen_Dat[SNP %in% SNPs_ok])
  
  #####Recoding SNP values in new data#####
  
  # Creating vector of valid genotypes
  genotypes <- apply(expand.grid(c("A", "T", "G", "C"), c("A", "T", "G", "C")), 
                     1, FUN = paste, collapse = "")
  
  # Setting values not a valid genotype to missing
  temp <- Gen_Dat[, c(1:2)]
  Gen_Dat[Gen_Dat != genotypes[1] & Gen_Dat != genotypes[2] & Gen_Dat != genotypes[3] &
            Gen_Dat != genotypes[4] & Gen_Dat != genotypes[5] & Gen_Dat != genotypes[6] &
            Gen_Dat != genotypes[7] & Gen_Dat != genotypes[8] & Gen_Dat != genotypes[9] & 
            Gen_Dat != genotypes[10] & Gen_Dat != genotypes[11] & Gen_Dat != genotypes[12] &
            Gen_Dat != genotypes[13] & Gen_Dat != genotypes[14] & Gen_Dat != genotypes[15] & 
            Gen_Dat != genotypes[16]] <- NA 
  
  Gen_Dat[, c(1:2)] <- temp
  remove(temp)
  
  # Creating function for converting from string genotype to MAF
  Converter <- function(column, key = Gen_Dat[, 2]) {
    for(r in 1:length(column)) {
      if (!is.na(column[r])){
        # Converts from string genotype to MAF
        column[r] <- str_count(str_sub(column[r], 1, 1), 
                               str_sub(key[r], -1, -1)) + 
          str_count(str_sub(column[r], -1, -1), 
                    str_sub(key[r], -1, -1))
      }
    }
    return(as.integer(column))
  }
  
  # Converting from genotype to MAF using data.table code
  Gen_Dat <- setDT(Gen_Dat)
  Gen_Dat <- Gen_Dat[, (colnames(Gen_Dat)[-c(1:2)]) := lapply(.SD, Converter), 
                     .SDcols = colnames(Gen_Dat)[-c(1:2)]]
  
  # Saving data
  fwrite(Gen_Dat, file = paste0("Partitions/", population, "_Gen_Dat_Clean_Check", p, ".csv"))
  print(paste0("Partition ", p, " complete"))
  }
  stopCluster(cl) #Stop executing in parallel
}

