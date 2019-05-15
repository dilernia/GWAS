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
library(foreach)
library(doParallel)
library(qqman)
library(rsnps)

#dir <- "/panfs/roc/groups/15/petersla/diler001/Genetics/"
dir <- "C:/Users/Owner/Google Drive/Genetics_Research/New/"

path <- dir

setwd(dir)
population <- "CEU"

#####Initial downloading of data - only need to do this section once#####

if(file.exists("CEU_Gen_Dat_Check.csv") == FALSE) {
  load("file_names.Rda")
  
  file_names <- c(file_names[which(grepl(population, file_names, fixed=TRUE))])
  
  URL <- "ftp://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/2010-08_phaseII+III/forward/"
  
  URLs = paste0(URL, file_names)
  
  Gen_Get <- function(index, urls = URLs, pop)
  {
    files <- file_names[which(grepl(pop, file_names, fixed=TRUE))]
    
    if(file.exists(substr(x = files[[index]], start=0, stop = (nchar(files[[index]]) - 3))) == FALSE) {
      download.file(url = urls[[index]], destfile = files[[index]])
      gunzip(files[[index]], overwrite=TRUE)
    }
    
    fread(substr(x = files[[index]], start=0, stop = (nchar(files[[index]]) - 3)), fill = TRUE)
  }
  
  Gen_Dat <- rbindlist(lapply(FUN = Gen_Get, 
                              X = 1:length(file_names), pop = population))
  
  # Subsetting by columns
  Gen_Dat <- Gen_Dat[, c("chrom", "pos", "strand", "assembly#", "center",
                         "protLSID", "assayLSID", "panelLSID", "QCcode") := NULL]
  
  fwrite(Gen_Dat, file = paste0(population, "_Gen_Dat_Check.csv"))
}

#####End initial downloading of data#####

# Obtaining vector of population cell lines
load(paste0("Original_Data/", population, "_IC_data.Rda"))
Cell_lines <- get(paste0(population, "_IC_data"))$Cell_line
remove(list = paste0(population, "_IC_data"))

print(paste0("Number of cell lines: ", length(Cell_lines)))

if(file.exists(paste0(population, "Gen_Dat_Clean_check.csv")) == FALSE) {
  
  #####Cleaning new genotype data#####
  
  numRows <- nrow(fread(paste0(population, "_Gen_Dat_Check.csv"), select = 1)) + 1
  parts <- 1000
  chunkSize <- ceiling(numRows/parts)
  lastSize <- numRows - chunkSize*(parts-1)
  
  genColNames <- colnames(fread(paste0(population, "_Gen_Dat_Check.csv"), 
                                 header = T, nrows = 0))

  for(p in 1:parts) {
  
    # Loading partition of new genotype data
    Gen_Dat <- fread(paste0(population, "_Gen_Dat_Check.csv"), 
                     header = (p==1), skip = (p-1)*chunkSize, 
                     nrows = (p==parts)*lastSize + (p!=parts)*chunkSize)
    colnames(Gen_Dat) <- genColNames
    
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
    
    # Creating function for converting from string genotype to MAF
    Converter <- function(column, key = temp[, 2]) {
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
    print(paste0(Sys.time(), ", Partition ", p, " complete"))
  }

# Reading in cleaned partitions and aggregating for analysis
Gen_Dat <- rbindlist(lapply(X = str_subset(list.files("Partitions/"), pattern = population),
  FUN = function(x) {
    print(x); return(fread(paste0("Partitions/", x))[,ALLELES:=NULL])
    }))

Cell_line <- colnames(Gen_Dat)[-c(1)]
SNP <- Gen_Dat[, "SNP", with = F]

# Transposing genotype values (long to wide format)
Gen_Dat <- as.data.table(t(as.matrix(Gen_Dat[, SNP := NULL])))
Gen_Dat <- Gen_Dat[, Cell_lines := Cell_line] 

# Adding back in column names
colnames(Gen_Dat) <- c(unlist(SNP), "Cell_line")

# Rearranging columns
Gen_Dat <- Gen_Dat[, c("Cell_line", unlist(SNP)), with = FALSE]
remove(SNP)

fwrite(Gen_Dat, file = paste0(population, "Gen_Dat_Clean_check.csv"))
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

numCols <- ncol(fread(paste0(population, "Gen_Dat_Clean_check.csv"), 
                      nrows = 0))-1
parts <- 100
chunkSize <- floor(numCols/parts)

startStop <- data.frame(starts = seq(2, numCols+1, by = chunkSize), 
           stops = c(seq(1, numCols+1, by = chunkSize)[-c(1)], numCols))

for(p in 1:nrow(startStop)) {
  # Loading partition of cleaned genotype data
  Gen_Dat <- fread(paste0(population, "Gen_Dat_Clean_check.csv"),
                   select = c(1, startStop[p, 1]:startStop[p, 2]))
  
  # Running mass univariate analysis for given genotype data subset
  detectCores()
  ncores <- 15
  cl <- makeCluster(ncores) # creates a cluster with <ncore> cores
  registerDoParallel(cl) # register the cluster
  
  results <- foreach(j = 5:ncol(PhenData_Combined)) %dopar% {
    library(data.table)
    library(tidyverse)
    
    # Creating function for fitting MLR model
    MLR <- function(gendCol, snpCol, phenCol) {
      temp <- data.frame(gender = gendCol, snp = snpCol, phen = phenCol)
      if(temp[complete.cases(temp), c("gender")] %>% unique() %>% length() > 1) {
        fit <- lm(phenCol ~ gendCol + snpCol)
        if(dim(summary(fit)$coefficients)[1]>2) {
          tStat <- summary(fit)$coefficients[3,3]
          pval <- summary(fit)$coefficients[3,4]
          return(list(tStat, pval))
        } else{return(list(NA, NA))}
      } else{return(list(NA, NA))}
    }
    
    tempResult <- lapply(FUN = MLR, X = as.data.frame(Gen_Dat)[, -c(1)], 
                         gendCol = PhenData_Combined$Gender, 
                         phenCol = PhenData_Combined[, j])
    out <- list(tempResult, colnames(PhenData_Combined)[j], Sys.time())
    names(out) <- c("tVal_pVal", "phenotype", "time")
    return(out)
  }
  
  stopCluster(cl) #Stop executing in parallel
  
  saveRDS(results, paste0("Analysis_Results_Check/", 
                          population, "_Mass_Uni_Results_", p, ".rds"))
}

# Aggregating GWAS Results ------------------------------------------------

readSig <- function(p, population, thresh = 1.01) {
  temp <- readRDS(paste0("Analysis_Results_Check/", population, 
                         "_Mass_Uni_Results_", p, ".rds"))
  
  # Initialize empty data table
  combined <- data.table(
    snp = character(),
    tVal = numeric(),
    pVal = numeric(),
    phenotype = character())
  
  for(i in 1:length(temp)) {
    results <- data.table(snp = names(temp[[i]]$tVal_pVal), 
                          tVal = unlist(temp[[i]]$tVal_pVal)[seq(1, length(unlist(temp[[i]]$tVal_pVal)), by = 2)], 
                          pVal = unlist(temp[[i]]$tVal_pVal)[seq(2, length(unlist(temp[[i]]$tVal_pVal)), by = 2)],
                          phenotype = temp[[i]]$phenotype)
    combined <- rbindlist(list(combined, results))
  }
  print(p)
  return(combined)
}

if(file.exists(paste0(population, "_Combined_Results.csv")) == FALSE) {
Combined_Results <- rbindlist(lapply(FUN = readSig, X = 1:101,
                                     population = population))
fwrite(Combined_Results, paste0(population, "_Combined_Results.csv"))
}

# Creating Manhattan and QQ Plots ------------------------------------------------

# Creating vector of phenotypes
phenos <- c("Avg_X7mG0", "Avg_X7mG6", "Avg_X6mG0", "Avg_X6mG6", "Avg_r7v6.0", "Avg_r7v6.6", "Avg_IC20",     
            "Avg_IC50", "Avg_IC80", "Avg_logmut.0", "Avg_logmut.10", "Avg_logmut.20", "Mut10v0", 
            "Mut20v0", "Mut20v10")

# Preparing genotype data for Manhattan and QQ plots
if(file.exists(paste0(population, "_ManQQ_Check.csv")) == FALSE) {
  load("file_names.Rda")
  
  file_names <- c(file_names[which(grepl(population, file_names, fixed=TRUE))])
  
  URL <- "ftp://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/2010-08_phaseII+III/forward/"
  
  URLs = paste0(URL, file_names)
  
  Gen_Get <- function(index, urls = URLs, pop)
  {
    files <- file_names[which(grepl(pop, file_names, fixed=TRUE))]
    
    if(file.exists(substr(x = files[[index]], start=0, stop = (nchar(files[[index]]) - 3))) == FALSE) {
      download.file(url = urls[[index]], destfile = files[[index]])
      gunzip(files[[index]], overwrite=TRUE)
    }
    
    fread(substr(x = files[[index]], start=0, stop = (nchar(files[[index]]) - 3)), fill = TRUE)
  }
  
  Gen_Dat <- rbindlist(lapply(FUN = Gen_Get, 
                              X = 1:length(file_names), pop = population))
  
  # Subsetting by columns
  Gen_Dat <- Gen_Dat[, c("strand", "assembly#", "center",
                         "protLSID", "assayLSID", "panelLSID", "QCcode") := NULL]
  
  fwrite(Gen_Dat, file = paste0(population, "_ManQQ_Check.csv"))
}

# Reading in GWAS results for specified phenotype
for(p in phenos) {
Combined_Results <- fread(paste0(population, "_Combined_Results.csv"))[phenotype == p][, .(SNP = `snp`, P = pVal)]

# Reading in subset of genotype data
manDat <- fread(paste0(population, "_ManQQ_Check.csv"))[, .(BP = pos, CHR = chrom, SNP = `rs#`)]

# Merging with GWAS results for given phenotype
manDat <- merge(manDat, Combined_Results, by = "SNP")
remove(Combined_Results)

# Removes X, Y, and M chromosome data
manDat$CHR <- manDat$CHR %>% gsub(pattern = "chr", replacement = "") %>% 
  as.integer()
manDat <- manDat[complete.cases(manDat), ]

# Creating Manhattan and QQ plots for each outcome
jpeg(paste0("Plots/New/", population, "_", p, "_Manhat", ".jpeg"))
manhattan(manDat, main = gsub(p, pattern = "Avg_", replacement = ""))
dev.off()
jpeg(paste0("Plots/New/", population, "_", p, "_QQplot", ".jpeg"))
qq(manDat$P,  main = gsub(p, pattern = "Avg_", replacement = ""))
dev.off()
}

# Creating Boxplots -------------------------------------------------------

# Loading and organizing phenotype data

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

# Loading subset of genotype data
impSnps <- c("rs11258248", "rs3781080", "rs6848554", "rs6823445")
impGenes <- c("MCM10", "MCM10", "LOC107986236", "LOC107986236")
boxGen <- fread(paste0(population, "Gen_Dat_Clean_check.csv"), select = c("Cell_line", impSnps))

# Making genotypes a factor
as.data.frame(boxGen)[, impSnps] <- lapply(as.data.frame(boxGen)[, impSnps],
                                           FUN = function(x) {temp <- factor(x);
                                           levels(temp) <- c('aa', 'Aa', 'AA'); return(temp)})

# Merging phenotype and genotype data
boxData <- PhenData_Combined %>% select(Cell_line, Avg_IC20, Avg_IC50) %>% 
  left_join(boxGen, by = c("Cell_line" = "Cell_line"))

# Vector of phenotypes to make boxplots for
phens <- c("Avg_IC20", "Avg_IC50")

for(s in 1:length(impSnps)) {
for(p in phens) {
  snpNum <- s
  pheno <- p
  
# Creating box plots
pLab1 <- gsub(pheno, pattern = "Avg_", replacement = "")
if(str_detect(pheno, pattern = "IC")) {
pLab <- paste0("log(", pLab1, ")")
  } else {
    pLab <- pLab1
}

# Recoding from MAF to AA, Aa, aa format
boxData[, impSnps] <- lapply(boxData[, impSnps], FUN = function(x) {
                               dplyr::recode(as.character(x), 
                                             "0" = "AA", "1" = "Aa", "2" = "aa")})

colVec <- c("#009E73", "#0072B2", "#D55E00")

ggData <- boxData %>% select(c(phens, impSnps[snpNum])) %>% 
  rename("Genotype" = impSnps[snpNum]) %>% filter(!is.na(Genotype)) 
myGG <- ggData %>% ggplot(aes(x = factor(Genotype), y = get(pheno), 
                              fill = factor(Genotype))) + 
  geom_boxplot() + labs(x = "Genotype", y = pLab, 
   title = paste0(pLab1,
   " \n SNP: ", impSnps[snpNum], ", Gene: ", impGenes[snpNum])) +
  scale_fill_manual(values = colVec[which(c("aa", "Aa", "AA") %in% 
                                    sort(unique(ggData$Genotype)))]) +
  ylim(0, 5) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                     legend.position = "none")
ggsave(myGG, device = "pdf",
       file = paste0("Plots/New/", population, "_", p, "_", impSnps[snpNum], "_box.pdf"))
}
}
