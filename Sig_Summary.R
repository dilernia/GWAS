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
library(broom)
library(xlsx)
library(rsnps)
library(qqman)

# SNPs for MCM10 and MGMT were obtained here: https://www.ncbi.nlm.nih.gov/guide/howto/view-all-snps/

if(file.exists("sigSNPres.rds") == FALSE) {
  # Reading in GWAS results and only selecting significant snps
  readSig <- function(p, population = "CEU", snps = sigSNPs) {
    temp <- readRDS(paste0("Analysis_Results/", population, "_Mass_Uni_Results_", p, ".rds"))
    
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
  sigSNPs <- bind_rows(readRDS("Analysis_Results/CEU_Combined_Results_Sig.rds"), 
                       readRDS("Analysis_Results/YRI_Combined_Results_Sig.rds")) %>% 
    select(snp) %>% unlist() %>% unique()
  
  inputs <- expand.grid(p = 1:20, population = c("CEU", "YRI"), stringsAsFactors = FALSE)
  
  sigSNPres <- pmap_dfr(.f = readSig, .l = list(p = inputs$p, population = inputs$population))
  sigSNPres <- sigSNPres %>% filter(!is.na(tVal) & phenotype != "Viability") %>% 
    mutate(phenotype = gsub(x = phenotype, pattern = "Avg_", replacement = ""))
  
  # Saving significant SNP results
  saveRDS(sigSNPres, "sigSNPres.rds")
}
sigSNPres <- readRDS("sigSNPres.rds")

# Exploring variation explained by significant SNPs -----------------------

# Determining significant SNP's and phenotypes
snpPheno <- sigSNPres %>% filter(pVal <= 5*10^(-8))
remove("sigSNPres")

# Importing cleaned phenotype data
dir <- "/home/adilernia/andrew.old/Gdrive_Backup/Tobacco_Research/"
phenData <- rbind(readRDS(paste0(dir, "CEU_Pheno_Clean.rds")), 
                  readRDS(paste0(dir, "YRI_Pheno_Clean.rds")))
colnames(phenData) <- colnames(phenData) %>% gsub(pattern = "Avg_", replacement = "")
phenFull <- phenData
phenData <- phenData[, c("Cell_line", "Ethnicity", "Gender", unique(snpPheno$phenotype))]

# Importing cleaned genotype data
if(sum(c(file.exists("sigCEUGenDat.rds"), file.exists("sigYRIGenDat.rds"))) < 2) {
ceuGenDat <- as.data.frame(setDT(readRDS("CEU_Gen_Dat_Clean2.Rda"))[, 
                           c("Cell_line", snpPheno[which(snpPheno$pop == "CEU"), c("snp")]), 
                           with = FALSE])
ceuGenDat <- ceuGenDat[, !duplicated(colnames(ceuGenDat))]
saveRDS(ceuGenDat, "sigCEUGenDat.rds")
yriGenDat <- as.data.frame(setDT(readRDS("YRI_Gen_Dat_Clean2.Rda"))[, 
                           c("Cell_line", snpPheno[which(snpPheno$pop == "YRI"), c("snp")]), 
                           with = FALSE])
yriGenDat <- yriGenDat[, !duplicated(colnames(yriGenDat))]
saveRDS(yriGenDat, "sigYRIGenDat.rds")
}

fullData <- phenData %>% left_join(readRDS("sigCEUGenDat.rds"), by = c("Cell_line" = "Cell_line")) %>% 
  left_join(readRDS("sigYRIGenDat.rds"), by = c("Cell_line" = "Cell_line")) %>% as.data.frame()
fullData[which(!(colnames(fullData) %in% c("Cell_line", "Ethnicity", "Gender")))] <- sapply(fullData[which(!(colnames(fullData) %in% c("Cell_line", "Ethnicity", "Gender")))], 
                                                                                            FUN = function(x) {as.numeric(as.character(x))})

popPhen <- snpPheno[!duplicated(snpPheno[, c("pop", "phenotype")]), c("pop", "phenotype")]

# Reobtaining results for my own sanity
resRerun <- function(x) {
pop <- popPhen[x, "pop"]
pheno <- popPhen[x, "phenotype"]
snps <- snpPheno[which(snpPheno$phenotype == pheno & snpPheno$pop == pop), "snp"]

MLR <- function(gendCol, snpCol, phenCol) {
  temp <- data.frame(g = gendCol, s = snpCol)
    fit <- lm(phenCol ~ gendCol + snpCol)
      tStat <- summary(fit)$coefficients[3,3]
      pval <- summary(fit)$coefficients[3,4]
      return(c(tStat, pval))
}

results <- data.frame(snp = character(), tVal = numeric(), pVal = numeric(), 
                      phenotype = character(), pop = character(), stringsAsFactors = FALSE)

for(i in 1:length(snps)) {
tpVals <- MLR(gendCol = fullData[which(fullData$Ethnicity == pop), "Gender"],
    phenCol = fullData[which(fullData$Ethnicity == pop), pheno],
    snpCol = fullData[which(fullData$Ethnicity == pop), snps[i]])
results[i, ] <- c(snps[i], tpVals[1], tpVals[2], pheno, pop)
}

return(results)
}

testRes <- map_dfr(.x = 1:nrow(popPhen), .f = resRerun) 

# Function for imputing phenotype values with mean
impute <- function(missCol = "IC80", datf = fullData) {
  impVal <- mean(datf[, missCol], na.rm = TRUE)
  datf[which(is.na(datf[, missCol])), missCol] <- impVal
  return(datf[, missCol])
}

impIC80 <- impute(missCol = "IC80", datf = fullData)
implogmut.20 <- impute(missCol = "logmut.20", datf = fullData)

# Function for imputing phenotype values using sample proportions
imputeSNP <- function(missCol = "rs7136711", datf) {
  res <- table(datf[, missCol]) %>% as.data.frame() %>% mutate(Prob = Freq / sum(Freq))
  res <- lapply(res, FUN = function(x){as.numeric(as.character(x))})
  datf[is.na(datf[, missCol]), missCol] <- sample(res$Var1, size = sum(is.na(datf[, missCol])), replace = TRUE, prob = res$Prob)
  datf[, missCol]
}

yriSigSnps <- unique(snpPheno[which(snpPheno$pop == "YRI"), "snp"])
impYRIsnps <- lapply(X = 1:length(yriSigSnps), FUN = function(x) {imputeSNP(missCol = yriSigSnps[x], 
                                               datf = fullData[which(fullData$Ethnicity == "YRI"), yriSigSnps])}) %>% 
  as.data.frame()
colnames(impYRIsnps) <- yriSigSnps

ceuSigSnps <- unique(snpPheno[which(snpPheno$pop == "CEU"), "snp"])
impCEUsnps <- lapply(X = 1:length(ceuSigSnps), FUN = function(x) {imputeSNP(missCol = ceuSigSnps[x], 
                                               datf = fullData[which(fullData$Ethnicity == "CEU"), ceuSigSnps])}) %>% as.data.frame()
colnames(impCEUsnps) <- ceuSigSnps

# Random imputation based on sample proportions
imputeR2 <- function(seed) {
  
  set.seed(seed)

# Imputing genotype values using proportions estimated using modleing above
impData <- fullData
impData[, c("IC80", "logmut.20")] <- cbind(impIC80, implogmut.20)
impData[which(impData$Ethnicity == "CEU"), ceuSigSnps] <- impCEUsnps
impData[which(impData$Ethnicity == "YRI"), yriSigSnps] <- impYRIsnps

# Function for summarize R^2 for each set of significant SNPs
rSqrSummary <- function(p) {
pop <- popPhen[p, "pop"]
pheno <- popPhen[p, "phenotype"]
sigSnps <- snpPheno[which(snpPheno$pop == pop & snpPheno$phenotype == pheno), "snp"]
phenSub <- impData[which(impData$Ethnicity == pop), c(pheno, "Gender", sigSnps)]
if(length(sigSnps) == 1) {
phenSub[, sigSnps] <- as.numeric(phenSub[, sigSnps])
} else {
  phenSub[, sigSnps] <- lapply(phenSub[, sigSnps], FUN = as.numeric)
}

fit <- lm(paste(pheno, " ~ Gender"), data = phenSub)
fitSNPs <- lm(paste(pheno, " ~ ."), data = phenSub)
return(data.frame(Population = pop, Phenotype = pheno, SNPs = paste(sigSnps, collapse = ", "), 
           R2_SNPs = round(summary(fitSNPs)$r.squared - summary(fit)$r.squared, 4),
           R2_Gend = round(summary(fit)$r.squared, 4),
           R2_Full = round(summary(fitSNPs)$r.squared, 4)))
}

aggRes <- map_dfr(.f = rSqrSummary, .x = 1:nrow(popPhen))
return(aggRes)
}

r2Res <- map_dfr(.f = imputeR2, .x = 1:100) %>% group_by(Population, Phenotype, SNPs) %>% 
  summarize(R2_SNPs = mean(R2_SNPs), R2_Gend = mean(R2_Gend), R2_Full = mean(R2_Full))

# Gene specific summaries -------------------------------------------------

# Function for determining MAF's for sample
mafCalc <- function(pop, snps) {
  #Loading new genotype data to obtain chromosome and position data for snps
  Gen_Dat <- readRDS(paste0(pop, "_Gen_Dat_Clean2.Rda"))
  
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
  write.xlsx(sigSNPres %>% filter(phenotype == p) %>% arrange(snp), paste0(p, "_results.xlsx"), row.names = FALSE)
}

# MGMT results
for(p in phenotypes) {
  write.xlsx(mgmtSummary %>% filter(phenotype == p) %>% arrange(SNP), paste0("MGMT_", p, "_results.xlsx"), row.names = FALSE)
}

# MCM10 results
for(p in phenotypes) {
  write.xlsx(mcm10Summary %>% filter(phenotype == p) %>% arrange(SNP), paste0("MCM10_", p, "_results.xlsx"), row.names = FALSE)
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

write.xlsx(mcm10Positions, paste0("mcm10Positions.xlsx"), row.names = FALSE)

mgmtPositions <- bpFind(mgmtsnps)

write.xlsx(mgmtPositions, paste0("mgmtPositions.xlsx"), row.names = FALSE)


# Manhattan and QQ Plots --------------------------------------------------

# Reading in Excel file
sigResults <- read.xlsx(paste0("resultsSummary.xlsx"), 1)
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

# Creating Dot Plots ------------------------------------------------------

phenVar <- "IC50raw"
phenLabs <- list(`r7v6.6` = bquote('7-mG/'*~O^6*'-mG at 6 h'),
                 IC50 = bquote(~IC[50]*" ("*mu*"M NMUr)"))

# Transforming IC measures back to original units
phenFull <- phenFull %>% mutate(IC20raw = exp(IC20), IC50raw = exp(IC50), IC80raw = exp(IC80))

meanData <- phenFull %>% group_by(Ethnicity) %>% summarize(Mean = mean(get(phenVar)))

phenFull %>% ggplot(aes(x = get(phenVar), fill = Ethnicity)) + geom_dotplot(stackdir = "center") +
  facet_grid(rows = vars(Ethnicity)) + theme_bw() + theme(legend.position = "none") + 
  scale_y_continuous(NULL, breaks = NULL) + scale_fill_manual(values=c("royalblue1", "green4")) +
  labs(x = bquote(~IC[50]*" ("*mu*"M NMUr)")) + 
  geom_vline(data = meanData, aes(xintercept = Mean), linetype = "dashed")

ggsave(paste0("dotPlot_", phenVar, ".png"), device = "png",
       width = 77*2, height = 51.5*2, units = "mm")

# Effect of Age or Gender -------------------------------------------------

gRes <- map_dfr(.x = colnames(phenFull)[which(!(colnames(phenFull) %in% c("Cell_line", "Ethnicity", "Gender")))], 
        .f = function(phen){
fit <- lm(get(phen) ~ Gender, data = phenFull)
res <- broom::tidy(fit)[2, ]
res$term <- phen
return(res)
})

gRes <- gRes %>% filter(!(term %in% c("IC20raw", "IC50raw", "IC80raw", "turnover"))) %>% 
  rename(Phenotype = term, Effect = estimate, SE = std.error) %>% 
  select(Phenotype, Effect, SE, `p.value`) %>% mutate(Phenotype = gsub(gsub(gsub(Phenotype, pattern = "IC20", replacement = "log.IC20"),
                                                  pattern = "IC50", replacement = "log.IC50"),
                                                  pattern = "IC80", replacement = "log.IC80"))

# Saving to csv
write.csv(gRes, "genderEffect.csv")
