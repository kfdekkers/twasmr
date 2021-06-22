# set seed and options

set.seed(1)
options(stringsAsFactors = F)

# load libraries

.libPaths(c(.libPaths(), "/home/kfdekkers/researchdrive/kdekkers/Scratch/library/"))

library(edgeR)
library(SNPlocs.Hsapiens.dbSNP142.GRCh37)
library(BBMRIomics)

cores <- 16

# set working directory

setwd("/home/kfdekkers/researchdrive/kdekkers/Scratch/MRTWAS/")

# source functions
source("Scripts/functions.R")

# load data

bbmri.data(rnaSeqData_ReadCounts_BIOS_cleaned)
load("predcells.rdata")
gwas <- read.delim("gwas.txt")
gwas$SNP <- gsub("chr", "", gwas$SNP_hg19)

# get dosages

ranges <- GRanges(seqnames = seqnames(counts), ranges = ranges(counts))
counts <- counts[, which(!is.na(counts$imputation_id))]
snps <- unique(gwas$SNP)
chr <- gsub("[:].*$", "", snps)
loc <- as.numeric(gsub("^.*[:]", "", snps))
dosages <- GRanges(seqnames = chr, ranges=IRanges(start = loc, end = loc))
dosages <- getGenotypes(imputation_id = counts$imputation_id, snps = dosages, type = "HRC", geno = "DS", BASE = VM_BASE_DATA, BPPARAM = MulticoreParam(cores))

# get reference and alternative allele info

folder <- "/home/kfdekkers/researchdrive/RSC_BIOS/RP3_data/HRC_Imputation/CODAM/results/unzipped/"
files <- list.files(folder, recursive = T, pattern = "info")
files <- files[!grepl("md5", files)]
files <- unlist(sapply(unique(paste0("chr", chr, ".")), function(x) files[grepl(x, files, fixed = T)]))

info <- bplapply(files, function(file) {
  
  data <- read.delim(paste0(folder, file))
  data <- data[, c("SNP", "REF.0.", "ALT.1.")]
  data[which(data$SNP %in% snps), ]
  
}, BPPARAM = MulticoreParam(cores))

info <- do.call(rbind, info)

save(dosages, ps, file="Data/data2.rdata")

# make effect allele the reference allele

gwas$alt <- tolower(info$ALT.1.[match(gwas$SNP, info$SNP)])
gwas$ref <- tolower(info$REF.0.[match(gwas$SNP, info$SNP)])
dosages <- dosages[match(gwas$SNP, rownames(dosages)), ]
rownames(dosages) <- gwas$rsid
colnames(dosages) <- counts$uuid[match(colnames(dosages), counts$imputation_id)]
dosages <- t(dosages)

for (i in 1:nrow(gwas)) {
  
  if (gwas$A1[i] == gwas$alt[i]) {
    
    dosages[, i] <- 2 - dosages[, i]
    gwas$alt[i] <- gwas$A2[i]
    gwas$ref[i] <- gwas$A1[i]
    
  }
  
}

# check

table(paste0(gwas$A1, gwas$A2) == paste0(gwas$ref, gwas$alt))

# get and clean variables

lipids <- c("tg", "ldl", "hdl")
ss <- colData(counts)
ldl <- ss$TotChol - ss$HDLchol - (ss$Triglycerides / 5)
var <- na.omit(data.frame(ss[c("Triglycerides", "HDLchol")], ldl))
var <- var[, c(1, 3, 2)]
colnames(var) <- lipids

wbcc <- wbcc[match(rownames(ss), rownames(wbcc)), ]

cells <- lapply(unique(ss$biobank_id), function(biobank) {
  
  if (biobank %in% c("CODAM", "PAN")) res <- wbcc[which(ss$biobank_id == biobank), ]
  if (biobank %in% c("LLS")) res <- data.frame(ss[which(ss$biobank_id == biobank), c("Neut_Perc", "Lymph_Perc", "Mono_Perc", "Eos_Perc", "Baso_Perc", "WBC")], wbcc[which(ss$biobank_id == biobank), "RBC", drop = F])
  if (biobank %in% c("RS")) res <- data.frame(ss[which(ss$biobank_id == biobank), c("Lymph_Perc", "Mono_Perc", "WBC",  "RBC")], wbcc[which(ss$biobank_id == biobank), c("Eos_Perc", "Baso_Perc", "Neut_Perc")])
  if (biobank %in% c("LL", "NTR")) res <- ss[which(ss$biobank_id == biobank), ]
  res[, c("Neut_Perc", "Lymph_Perc", "Mono_Perc", "Eos_Perc", "WBC", "RBC")]
  
})

cells <- do.call(rbind, cells)

cov <- na.omit(data.frame(ss[, c("Sampling_Age", "Sex", "biobank_id")], cells))
cov$Sex <- as.numeric(as.factor(cov$Sex))
smoking <- ss$Smoking
medication <- ss$LipidMed
medication <- as.character(medication)
medication[which(!medication == "no")] <- "yes"
bmi <- ss$Weight / (ss$Height / 100) ^ 2
sens <- data.frame(ss[, colnames(wbcc)], smoking, medication, bmi)
sens$smoking <- as.numeric(as.factor(sens$smoking))
sens$medication <- as.numeric(as.factor(sens$medication))

# filter and normalize transcripts

counts <- counts[which(!seqnames(counts) %in% c("X", "Y")), ]
counts <- counts[which(rowSums(assays(counts)$data > 0) > 0.8 * ncol(counts)), ]

# match everything

id <- rownames(var)[which(rownames(var) %in% rownames(cov) & rownames(var) %in% rownames(dosages) & rownames(var) %in% rownames(sens) & rownames(var) %in% colnames(counts))]
var <- var[match(id, rownames(var)), ]
cov <- cov[match(id, rownames(cov)), ]
dosages <- dosages[match(id, rownames(dosages)), ]
sens <- sens[match(id, rownames(sens)), ]
counts <- counts[, match(id, colnames(counts))]

# check if all rownames and colnames match

all(sapply(list(rownames(var), rownames(cov), rownames(dosages), rownames(sens)), FUN = identical, colnames(counts)))
all(colnames(dosages) == gwas$rsid)

# RIN transform function

RIN <- function(x) qnorm((rank(x, "keep") - 0.5) / sum(!is.na(x)))

# RIN transform or scale variables

var <- apply(var, 2, RIN)
dosages <- apply(dosages, 2, scale)
rownames(dosages) <- rownames(var)
counts <- DGEList(counts = assays(counts)$data)
counts <- calcNormFactors(counts)
counts <- cpm(counts)
counts <- t(counts)
counts <- apply(counts, 2, RIN)

# create lipid polygenic scores

ps.fun <- function(dosages, weights) {
  
  scale(colSums(t(dosages) * weights))
  
}

lipids.ps <- lapply(lipids, function(lipid) ps.fun(dosages[, which(gwas$Lipid == lipid)], gwas$beta[which(gwas$Lipid == lipid)]))
lipids.ps <- do.call(cbind, lipids.ps)
colnames(ps) <- lipids

# get bmi and blood pressure dosages

bmi.gwas <- read.delim("bmi.tsv")
bp.gwas <- read.delim("bp.tsv")
bmi.bp <- rbind(bmi.gwas, bp.gwas)

snps <- snpsById(SNPlocs.Hsapiens.dbSNP142.GRCh37, bmi.bp$SNPS, ifnotfound = "drop")
seqlevels(snps) <- gsub("ch", "", seqlevels(snps))
bmi.bp.dosages <- getGenotypes(imputation_id = ss$imputation_id, snps = snps, type = "HRC", geno = "DS", BASE = VM_BASE_DATA, BPPARAM = MulticoreParam(cores))

snps <- as.data.frame(snps)
snps$SNP <- paste(snps$seqnames, snps$pos, sep = ":")
folder <- "/home/kfdekkers/researchdrive/RSC_BIOS/RP3_data/HRC_Imputation/CODAM/results/unzipped/"
files <- list.files(folder, recursive = T, pattern = "info")
files <- unlist(sapply(paste0("chr", unique(snps$seqnames)), function(x) files[grep(x, files)]))

info <- bplapply(files, function(file) {
  
  data <- read.delim(paste0(folder, file))
  data <- data[, c("SNP", "REF.0.", "ALT.1.")]
  data$file <- file
  data[which(data$SNP %in% snps$SNP), ]
  
}, BPPARAM = MulticoreParam(cores))

info <- do.call(rbind, info)

# make effect allele the reference allele

bmi.bp$SNP <- snps$SNP[match(bmi.bp$SNPS, snps$RefSNP_id)]
bmi.bp$alt <- info$ALT.1.[match(bmi.bp$SNP, info$SNP)]
bmi.bp$ref <- info$REF.0.[match(bmi.bp$SNP, info$SNP)]
bmi.bp$A1 <- substr(bmi.bp$STRONGEST.SNP.RISK.ALLELE, nchar(bmi.bp$STRONGEST.SNP.RISK.ALLELE), nchar(bmi.bp$STRONGEST.SNP.RISK.ALLELE))
bmi.bp <- bmi.bp[which(!is.na(bmi.bp$alt)), ]
bmi.bp.dosages <- bmi.bp.dosages[match(bmi.bp$SNP, rownames(bmi.bp.dosages)), ]
rownames(bmi.bp.dosages) <- bmi.bp$SNPS
colnames(bmi.bp.dosages) <- ss$uuid[match(colnames(bmi.bp.dosages), ss$imputation_id)]
bmi.bp.dosages <- t(bmi.bp.dosages)

for (i in 1:nrow(bmi.bp)) {
  
  if (bmi.bp$A1[i] == bmi.bp$alt[i]) {
    
    bmi.bp.dosages[, i] <- 2 - bmi.bp.dosages[, i]
    bmi.bp$ref[i] <- tolower(bmi.bp$A1[i])
    
  }
  
}

# create polygenic scores for bmi and bp

bmi.bp.ps <- lapply(unique(bmi.bp$MAPPED_TRAIT), function(x) {
  
  ps.fun(bmi.bp.dosages[, which(bmi.bp$MAPPED_TRAIT == x)], bmi.bp$OR.or.BETA[which(bmi.bp$MAPPED_TRAIT == x)])
  
})

names(bmi.bp.ps) <- unique(bmi.bp$MAPPED_TRAIT)
bmi.bp.ps <- do.call(cbind, bmi.bp.ps)
bmi.bp.ps <- bmi.bp.ps[match(rownames(var), rownames(bmi.bp.ps)), ]
bmi.bp.ps <- bmi.bp.ps[, c(1:3)]
colnames(bmi.bp.ps) <- c("BMI", "SBP", "DBP")

ps <- cbind(bmi.bp.ps, lipid.ps)

save(var, cov, dosages, counts, sens, ranges, ss, lipid.ps, bmi.bp.ps, bmi.bp.dosages, ps, file = "Data/data.rdata")
