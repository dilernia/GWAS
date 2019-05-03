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

#dir <- "C:\\Users\\Andrew\\Documents\\Research\\Tobacco Research\\Extra_Gen_Dat\\"
#path <- "C:\\Users\\Andrew\\Documents\\Research\\Tobacco Research\\Data Sets\\"
#dir <- "/home/andrew/andrew.old/Tobacco_Research/"
dir <- "/panfs/roc/groups/15/petersla/diler001/Genetics/"

path <- dir

setwd(dir)
population <- "CEU"

cleaned <- TRUE

#####Initial downloading of data - only need to do this section once#####

#load("file_names.Rda")
#
#file_names <- c(file_names[which(grepl("CEU", file_names, fixed=TRUE))], 
#               file_names[which(grepl("YRI", file_names, fixed=TRUE))])
#
#URL <- "ftp://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/2010-08_phaseII+III/forward/"
#
#URLs = paste0(URL, file_names)
#
#Gen_Get <- function(index, urls = URLs, pop = "CEU")
#{
#  files <- file_names[which(grepl(pop, file_names, fixed=TRUE))]
#  
#  if(file.exists(substr(x = files[[index]], start=0, stop = (nchar(files[[index]]) - 3))) == FALSE) {
#    download.file(url = urls[[index]], destfile = files[[index]])
#    gunzip(files[[index]], overwrite=TRUE)
#  }
#  
#  fread(substr(x = files[[index]], start=0, stop = (nchar(files[[index]]) - 3)), fill = TRUE)
#}
#
#CEU_Gen_Dat <- map_dfr(.f = Gen_Get, .x = as.list(1:25), 
#                       pop = "CEU")
#YRI_Gen_Dat <- map_dfr(.f = Gen_Get, .x = as.list(1:25), 
#                       pop = "YRI")
#
#save(CEU_Gen_Dat, file = "CEU_Gen_Dat.Rda")
#save(YRI_Gen_Dat, file = "YRI_Gen_Dat.Rda")
#
#####End initial downloading of data#####

#####Loading old SNP genotypes#####

#population <- "YRI"
#
#SNP_data <- fread(paste0(path,'snp_matrix.', population, '.csv'))
#
#Cleans up snp names
#colnames(SNP_data) <- colnames(SNP_data) %>% 
#  gsub(pattern = " ", replacement = "") %>% tolower()
#
#Only keeps columns with valid snp names
#SNP_data <- SNP_data %>% as.data.frame()
#SNP_data <- SNP_data[, c(1, which(str_detect(colnames(SNP_data), "rs")))]
#
#SNPs <- colnames(SNP_data)[-c(1)]
#
#save(SNPs, file = paste0("Original_snps_", population, ".Rda"))
#
#colnames(SNP_data)[1] <- "Cell_line"
#
#Cell_lines <- SNP_data$Cell_line %>% gsub(pattern = " ", replacement = "") %>% toupper()
#
#Converting all snp columns to integer type
#SNP_data <- sapply(SNP_data, as.integer)
#
#Adds back in Cell_line column
#SNP_data$Cell_line <- Cell_lines
#
#Saves original data
#save(SNP_data, file = paste0("snp_matrix.", population, ".Rda"))

if(cleaned == FALSE) {

#####Cleaning new genotype data#####

#Loading new genotype data
load(paste0(population, "_Gen_Dat.Rda"))
assign("Gen_Dat", get(paste0(population, "_Gen_Dat")))
assign(paste0(population, "_Gen_Dat"), NULL)

#Subsetting by columns
Gen_Dat <- as.data.frame(Gen_Dat[, c("chrom", "pos", "strand", "assembly#", "center",
                                     "protLSID", "assayLSID", "panelLSID", "QCcode") := NULL])

print(paste0("Current memory used: ", pryr::mem_used()))

#Loading snps to be excluded due to potential errors
SNPs_to_Exclude <- read_excel(paste0(path, "SNPs_to_Exclude.xlsx"))$rsid %>% unique() %>% 
  gsub(pattern = " ", replacement = "") %>% tolower()

#Cleaning up snp names
Gen_Dat <- Gen_Dat %>% dplyr::rename(SNP = `rs#`) 
Gen_Dat$SNP <- Gen_Dat$SNP %>% gsub(pattern = " ", replacement = "") %>% tolower()

#Cleaning up cell line names
colnames(Gen_Dat) <- colnames(Gen_Dat) %>% gsub(pattern = " ", replacement = "") %>% toupper()

#Subsetting to only keep necessary cell line columns
load("All_Cell_Lines.Rda")
Clean_Dat_Cell_IDs <- Clean_Dat_Cell_IDs %>% unlist() %>% trimws() %>% toupper()
Gen_Dat <- Gen_Dat[, c(1:2, which(colnames(Gen_Dat) %in% Clean_Dat_Cell_IDs))]

#Determining which snps are new and not included in old ones or in error list of snps
#SNPs_new <- Gen_Dat$SNP[which(!(Gen_Dat$SNP %in% c(SNPs, SNPs_to_Exclude)))] %>% unique() 

#Only keeps snps which were not used in original analysis and excludes errors from: ftp://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/2010-08_phaseII+III/inconsistencies/
#Gen_Dat <- Gen_Dat[which(!(Gen_Dat$SNP %in% c(SNPs_to_Exclude, CEU_Original_SNPs))), ]

#Excludes error snps from: ftp://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/2010-08_phaseII+III/inconsistencies/
Gen_Dat <- Gen_Dat[which(!(Gen_Dat$SNP %in% SNPs_to_Exclude)), ]

#####Recoding SNP values in new data#####

#Creating vector of valid genotypes
genotypes <- apply(expand.grid(c("A", "T", "G", "C"), c("A", "T", "G", "C")), 
                   1, FUN = paste, collapse = "")

#Setting values not a valid genotype to missing
temp <- Gen_Dat[, c(1:2)]
Gen_Dat[Gen_Dat != genotypes[1] & Gen_Dat != genotypes[2] & Gen_Dat != genotypes[3] &
          Gen_Dat != genotypes[4] & Gen_Dat != genotypes[5] & Gen_Dat != genotypes[6] &
          Gen_Dat != genotypes[7] & Gen_Dat != genotypes[8] & Gen_Dat != genotypes[9] & 
          Gen_Dat != genotypes[10] & Gen_Dat != genotypes[11] & Gen_Dat != genotypes[12] &
          Gen_Dat != genotypes[13] & Gen_Dat != genotypes[14] & Gen_Dat != genotypes[15] & 
          Gen_Dat != genotypes[16]] <- NA 

Gen_Dat[, c(1:2)] <- temp
remove(temp)

#Creating function for converting from string genotype to MAF
Converter <- function(column, key = Gen_Dat[, 2]) {
  for(r in 1:length(column)) {
    if (!is.na(column[r])){
      #Converts from string genotype to MAF
      column[r] <- str_count(str_sub(column[r], 1, 1), 
                                 str_sub(key[r], -1, -1)) + 
                         str_count(str_sub(column[r], -1, -1), 
                                   str_sub(key[r], -1, -1))
    }
  }
  print(paste0("Column complete. ", Sys.time()))
  return(as.integer(column))
}

#Trying with data.table code and function instead
Gen_Dat <- setDT(Gen_Dat)
Gen_Dat <- Gen_Dat[, (colnames(Gen_Dat)[-c(1:2)]) := lapply(.SD, Converter), .SDcols = colnames(Gen_Dat)[-c(1:2)]]

saveRDS(Gen_Dat, file = paste0(population, "_Gen_Dat_Clean.Rda"))

#####Restructuring genotype data for analysis#####
Gen_Dat <- readRDS(paste0(population, "_Gen_Dat_Clean.Rda"))
Cell_line <- colnames(Gen_Dat)[-c(1:2)]
Gen_Dat <- setDT(transpose(Gen_Dat))
colnames(Gen_Dat) <- unlist(Gen_Dat[1])
Gen_Dat <- Gen_Dat[-c(1:2)]
Gen_Dat <- cbind(Cell_line, Gen_Dat)
saveRDS(Gen_Dat, file = paste0(population, "_Gen_Dat_Clean2.Rda"))

}

#####Reading original phenotype data in#####

phenFiles <- paste0(dir, "Original_Data/", population, "_", 
                    c("Cyto", "Repair", "IC", "Mut", "Mutation"), "_data.Rda")
for(f in phenFiles){load(f)}

#Creates list of phenotype data sets
PhenDats <- lapply(X = paste0(population, "_", 
                              c("Cyto", "Repair", "IC", "Mut", "Mutation"), "_data"), 
                   FUN = function(x){eval(parse(text=x))})

#Creating function for cleaning cell line names
PhenDat_Cleaner <- function(dat) {
  dat$Cell_line <- dat$Cell_line %>% gsub(pattern = "\\s", replacement = "") %>% toupper()
  return(dat)
}
PhenDats <- lapply(X = PhenDats, FUN = PhenDat_Cleaner)

#Combining all phenotype datasets into one
PhenData_Combined <- PhenDats %>% 
  Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="Cell_line"), .) %>% 
  dplyr::select(-ends_with(".x")) %>% dplyr::select(-ends_with(".y")) %>% 
  dplyr::select(Cell_line, Ethnicity, Gender, everything())

#Calculating log transformations for IC data, recoding gender variable, log mut diffs
PhenData_Combined <- PhenData_Combined %>% mutate(Avg_IC20 = log(Avg_IC20), Avg_IC50 = log(Avg_IC50), 
                                      Avg_IC80 = log(Avg_IC80), Gender = factor(Gender), 
                                      Mut10v0 = log10((Mut10v0+abs(Mut10v0))/2 + 0.01), 
                                      Mut20v0 = log10((Mut20v0+abs(Mut20v0))/2 + 0.01),
                                      Mut20v10 = log10((Mut20v10+abs(Mut20v10))/2 + 0.01)) %>% 
  dplyr::select(-c(Mut0, Mut10, Mut20))

#####Looping through doing analysis for each partition of data#####
for(p in 1:20) {
#Loading new genotype data
Gen_Dat <- readRDS(paste0(population, "_Gen_Dat_Clean2.Rda"))
Gen_Dat <- setDT(Gen_Dat)

#Separating genotype data into 20 partitions
num_cols <- ncol(Gen_Dat); 
indices <- c(seq(1, num_cols, by = ceiling(num_cols/20)), num_cols)
indices[1] <- 0
Gen_Dat <- Gen_Dat[, (indices[p]+1):indices[p+1]] %>% as.data.frame()

#Converting columns to numeric type
#Gen_Dat <- Gen_Dat[, (colnames(Gen_Dat)[-c(1)]) := lapply(.SD, as.numeric), .SDcols = colnames(Gen_Dat)[-c(1)]]
Gen_Dat[, -c(1)] <- sapply(Gen_Dat[, -c(1)], as.numeric)

#Running mass univariate analysis for given genotype data subset
detectCores()
ncores <- 10 #Set ncores to something less than or equal to available cores
cl <- makeCluster(ncores) # creates a cluster with <ncore> cores
registerDoParallel(cl) # register the cluster

results <- foreach(j = 4:ncol(PhenData_Combined)) %dopar% {
  library(data.table)
  library(tidyverse)
  #Creating function for fitting MLR model
  MLR <- function(gendCol, snpCol, phenCol) {
    temp <- data.frame(g = gendCol, s = snpCol)
    if(temp[complete.cases(temp), c("g")] %>% unique() %>% length() > 1) {
    fit <- lm(phenCol ~ gendCol + snpCol)
    if(dim(summary(fit)$coefficients)[1]>2) {
    tStat <- summary(fit)$coefficients[3,3]
    pval <- summary(fit)$coefficients[3,4]
    return(list(tStat, pval))
    } else{return(list(NA, NA))}
    } else{return(list(NA, NA))}
  }
  
  tempResult <- lapply(FUN = MLR, X = as.data.frame(Gen_Dat)[, -c(1)], gendCol = PhenData_Combined$Gender, 
         phenCol = PhenData_Combined[, j])
  out <- list(tempResult, colnames(PhenData_Combined)[j], Sys.time())
  names(out) <- c("tVal_pVal", "phenotype", "time")
  return(out)
}

stopCluster(cl) #Stop executing in parallel

saveRDS(results, paste0(getwd(), "/Analysis_Results_Test/", population, "_Mass_Uni_Results_", p, ".rds"))
}
