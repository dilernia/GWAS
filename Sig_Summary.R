#################
# Title: Summarize significant SNPs
# Author: Andrew DiLernia
# Date: 02/20/2019
# Purpose: Summarize mass univariate results
#################

library(tidyverse)
library(foreach)
library(doParallel)
library(data.table)
library(xlsx)
library(rsnps)
library(qqman)

# SNPs for MCM10 and MGMT were obtained here: https://www.ncbi.nlm.nih.gov/guide/howto/view-all-snps/
Datdir <- "/home/andrew/andrew.old/Gdrive_Backup/Tobacco_Research/New/"

# Reading in GWAS results and only selecting significant snps
readSig <- function(p, population = "CEU", snps = sigSNPs) {
  temp <- readRDS(paste0(Datdir, "Analysis_Results/", population, "_Mass_Uni_Results_", p, ".rds"))
  
  #Initialize empty data table
  combined <- data.table(
    snp = character(),
    tVal = numeric(),
    pVal = numeric(),
    phenotype = character(),
    pop = character())
  
  for(i in 1:length(temp)) {
    results <- data.table(snp = names(temp[[i]]$tVal_pVal), 
                          tVal = unlist(temp[[i]]$tVal_pVal)[seq(1, length(unlist(temp[[i]]$tVal_pVal)), by = 2)], 
                          pVal = unlist(temp[[i]]$tVal_pVal)[seq(2, length(unlist(temp[[i]]$tVal_pVal)), by = 2)],
                          phenotype = temp[[i]]$phenotype,
                          pop = population)
    combined <- rbindlist(list(combined, results[snp %in% snps]))
  }
  print(p)
  return(combined)
}

# Determining significant SNPs
sigSNPs <- bind_rows(readRDS(paste0(Datdir, "Analysis_Results/CEU_Combined_Results_Sig.rds")), 
                     readRDS(paste0(Datdir, "Analysis_Results/YRI_Combined_Results_Sig.rds"))) %>% 
  select(snp) %>% unlist() %>% unique()

inputs <- expand.grid(p = 1:20, population = c("CEU", "YRI"), stringsAsFactors = FALSE)

sigSNPres <- pmap_dfr(.f = readSig, .l = list(p = inputs$p, population = inputs$population))
sigSNPres <- sigSNPres %>% filter(!is.na(tVal) & phenotype != "Viability") %>% 
  mutate(phenotype = gsub(x = phenotype, pattern = "Avg_", replacement = ""))

# Saving significant SNP results
saveRDS(sigSNPres, "sigSNPres.rds")

# Gene specific summaries -------------------------------------------------

# Function for determining MAF's for sample
mafCalc <- function(pop, snps) {
  #Loading new genotype data to obtain chromosome and position data for snps
  Gen_Dat <- readRDS(paste0(Datdir, pop, "_Gen_Dat_Clean2.Rda"))
  
  Gen_Dat <- Gen_Dat[, colnames(Gen_Dat)[which(colnames(Gen_Dat) %in% snps)], with=FALSE]
  
  result <- Gen_Dat %>% sapply(FUN = function(geno){mean(as.numeric(geno), na.rm = TRUE) / 2}) %>% as.data.frame() %>% 
    setDT(keep.rownames = TRUE) %>% setnames(1:2, c("SNP", paste0("MAF_", pop)))
  return(result)
}

# Calculating MAF and saving results for MCM10 gene snps
mcm10snps <- read_csv("MCM10_snps.txt") %>% unlist() %>% str_split(pattern = "\\[") %>% map(1) %>% trimws()
mcm10 <- pmap_dfr(.f = readSig, .l = list(p = inputs$p, population = inputs$population), snps = mcm10snps) %>% filter(!is.na(tVal))
saveRDS(mcm10, "mcm10Res.rds")

mcm10Summary1 <- mafCalc(pop = "CEU", snps = mcm10snps)
saveRDS(mcm10Summary1, "mcm10Summary1.rds")
mcm10Summary2 <- mafCalc(pop = "YRI", snps = mcm10snps)
saveRDS(mcm10Summary2, "mcm10Summary2.rds")
mcm10Summary <- full_join(mcm10Summary1, mcm10Summary2, by = c("SNP")) %>% full_join(mcm10, by = c("SNP" = "snp")) %>% 
  filter(!is.na(tVal) & phenotype != "Viability") %>% mutate(phenotype = gsub(x = phenotype, pattern = "Avg_", replacement = "")) %>% 
  arrange(phenotype, pVal)
saveRDS(mcm10Summary, "mcm10Summary.rds")

# Calculating MAF and saving results for MGMT gene snps
mgmtsnps <- read_csv("MGMT_snps.txt") %>% unlist() %>% str_split(pattern = "\\[") %>% map(1) %>% trimws()
mgmt <- pmap_dfr(.f = readSig, .l = list(p = inputs$p, population = inputs$population), snps = mgmtsnps)  %>% filter(!is.na(tVal))
saveRDS(mgmt, "mgmtRes.rds")

mgmtSummary1 <- mafCalc(pop = "CEU", snps = mgmtsnps)
saveRDS(mgmtSummary1, "mgmtSummary1.rds")
mgmtSummary2 <- mafCalc(pop = "YRI", snps = mgmtsnps)
saveRDS(mgmtSummary2, "mgmtSummary2.rds")
mgmtSummary <- full_join(mgmtSummary1, mgmtSummary2, by = c("SNP")) %>% full_join(mgmt, by = c("SNP" = "snp")) %>% 
  filter(!is.na(tVal) & phenotype != "Viability") %>% mutate(phenotype = gsub(x = phenotype, pattern = "Avg_", replacement = "")) %>% 
  arrange(phenotype, pVal)
saveRDS(mgmtSummary, "mgmtSummary.rds")

# Creating vector of phenotypes
phenotypes <- sigSNPres$phenotype %>% unique()

# Exporting to Excel

# Overall results
for(p in phenotypes) {
write.xlsx(sigSNPres %>% filter(phenotype == p) %>% arrange(snp), paste0(Datdir, p, "_results.xlsx"), row.names = FALSE)
}

# MGMT results
for(p in phenotypes) {
  write.xlsx(mgmtSummary %>% filter(phenotype == p) %>% arrange(SNP), paste0(Datdir, "MGMT_", p, "_results.xlsx"), row.names = FALSE)
}

# MCM10 results
for(p in phenotypes) {
  write.xlsx(mcm10Summary %>% filter(phenotype == p) %>% arrange(SNP), paste0(Datdir, "MCM10_", p, "_results.xlsx"), row.names = FALSE)
}

# Obtaining SNP positions -------------------------------------------------

# Reading in snp vectors for each gene
mcm10snps <- read_csv("MCM10_snps.txt") %>% unlist() %>% str_split(pattern = "\\[") %>% map(1) %>% trimws()
mgmtsnps <- read_csv("MGMT_snps.txt") %>% unlist() %>% str_split(pattern = "\\[") %>% map(1) %>% trimws()

# Obtaining positions of snps for MCM10 gene

# Function for determining BP for snps
bpFind <- function(snps) {
  
  #Loading new genotype data to obtain chromosome and position data for snps
  load("/home/andrew/andrew.old/Gdrive_Backup/Tobacco_Research/CEU_Gen_Dat.Rda")
  CEU_Gen_Dat <- CEU_Gen_Dat[`rs#` %in% snps, .(SNP = `rs#`, BP = pos)]
  CEU_Gen_Dat[, population := "CEU"]
  load("/home/andrew/andrew.old/Gdrive_Backup/Tobacco_Research/YRI_Gen_Dat.Rda")
  YRI_Gen_Dat <- YRI_Gen_Dat[`rs#` %in% snps, .(SNP = `rs#`, BP = pos)]
  YRI_Gen_Dat[, population := "YRI"]
  
  return(rbind(CEU_Gen_Dat, YRI_Gen_Dat))
}

mcm10Positions <- bpFind(mcm10snps)

write.xlsx(mcm10Positions, paste0(Datdir, "mcm10Positions.xlsx"), row.names = FALSE)

mgmtPositions <- bpFind(mgmtsnps)

write.xlsx(mgmtPositions, paste0(Datdir, "mgmtPositions.xlsx"), row.names = FALSE)


# Manhattan and QQ Plots --------------------------------------------------

# Reading in Excel file
sigResults <- read.xlsx(paste0(Datdir, "resultsSummary.xlsx"), 1)
Combined_Results <- readRDS(paste0(pop, "_Combined_Results_Sig.rds"))

outcomes <- Combined_Results$phenotype %>% unique()

# Loading new genotype data to obtain chromosome and position data for snps
load(paste0("/home/andrew/andrew.old/Gdrive_Backup/Tobacco_Research/", pop, 
            "_Gen_Dat.Rda"))
assign("Gen_Dat", get(paste0(pop, "_Gen_Dat")))
assign(paste0(pop, "_Gen_Dat"), NULL)

# Subsetting by columns
Gen_Dat <- as.data.frame(Gen_Dat[, c("strand", "assembly#", "center",
                                     "protLSID", "assayLSID", "panelLSID", "QCcode") := NULL])
Gen_Dat <- Gen_Dat %>% dplyr::select(`rs#`, chrom, pos) %>% dplyr::rename(snps = `rs#`)

for(p in outcomes) {
  outcome <- p
  
  # Filtering to only include given phenotype
  temp <- Combined_Results %>% filter(phenotype == outcome)
  
  # Merging results data to include position and chromosome for each snp
  ManHatDat <- inner_join(temp, Gen_Dat, by = "snps") %>% 
    dplyr::rename(BP = pos, CHR = chrom, P = pVals, SNP = snps)
  
  # Removes X, Y, and M chromosome data
  ManHatDat$CHR <- ManHatDat$CHR %>% gsub(pattern = "chr", replacement = "") %>% as.integer()
  ManHatDat <- ManHatDat[complete.cases(ManHatDat), ]
  
  # Creating and saving Manhattan Plot and QQ-Plot for given outcome
  # Saving manhattan plots
  jpeg(paste0(Plotsdir, population, "_", outcome, "_Manhat", ".jpeg"))
  manhattan(ManHatDat, main = paste0(population, " ", outcome))
  dev.off()
  # Saving qqplots
  jpeg(paste0(Plotsdir, population, "_", outcome, "_QQplot", ".jpeg"))
  qq(ManHatDat$P,  main = paste0(population, " ", outcome))
  dev.off()
}

#colnames(CEU_Gen_Dat)
#"rs#"       "alleles"   "chrom"     "pos"       "strand"    "assembly#" "center"   
# "protLSID"  "assayLSID" "panelLSID" "QCcode"    "NA06984"   "NA06985"
