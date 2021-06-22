# set seed and options

set.seed(1)
options(stringsAsFactors = F)

# load libraries

.libPaths(c(.libPaths(), "/etc/miniconda/lib/R/library/", "/home/kfdekkers/researchdrive/kdekkers/Scratch/library/"))
library(BiocGenerics, lib.loc = "/home/kfdekkers/researchdrive/Scratch/library/")
library(cate)
library(bacon)
library(BiocParallel)
library(rio)
library(BBMRIomics)

# set working directory

setwd("/home/kfdekkers/researchdrive/kdekkers/Scratch/MRTWAS/")

# source functions

source("Scripts/functions.R")

# lipids

lipids <- c("tg", "ldl", "hdl")

# perform twas analysis

perform.twas <- function(lipid, sens_var = NA, factors = 5) {
  
  if (!is.na(sens_var)) {
    
    cov <- data.frame(cov, sens[, sens_var])
    cov <- na.omit(cov)
    cov$biobank_id <- factor(cov$biobank_id)
    var <- var[match(rownames(cov), rownames(var)), ]
    counts <- counts[match(rownames(cov), rownames(counts)), ]
    
  }
  
  twas <- bplapply(unique(cov$biobank_id), function(biobank) {
    
    var <- var[which(cov$biobank_id == biobank), lipid]
    counts <- counts[which(cov$biobank_id == biobank), ]
    cov <- cov[which(cov$biobank_id == biobank), ]
    cov <- cov[, vapply(cov, function(x) length(unique(x)) > 1, logical(1L))]
    cov <- cov[, which(!colnames(cov) == "biobank_id")]
    z <- factor.fun(var, cov, counts, factors)
    cov <- cbind(cov, z)
    twas <- twas.fun.lm(var, cov, counts)
    twas_bc <- bacon.fun(twas)
    twas <- cbind(twas, twas_bc)
    twas
    
  }, BPPARAM=MulticoreParam(length(unique(cov$biobank_id))))
  
  names(twas) <- unique(cov$biobank_id)
  
  twas <- do.call(cbind, twas)
  es <- twas[, grepl("es_bc", colnames(twas))]
  se <- twas[, grepl("se_bc", colnames(twas))]
  twas_meta <- meta.fun(es, se)
  twas_meta <- cbind(twas, twas_meta)
  twas_meta$padj <- p.adjust(twas_meta$p, method = "bonf", n = sum(!is.na(twas_meta$p)))
  twas_meta
  
}

twas <- bplapply(lipids,  perform.twas, BPPARAM = MulticoreParam(length(lipids)))
names(twas) <- lipids
save(twas, file = "Data/twas.rdata")

twas.smoking <- bplapply(lipids,  function(lipid) perform.twas(lipid, "smoking"), BPPARAM = MulticoreParam(length(lipids)))
names(twas.smoking) <- lipids

twas.medication <- bplapply(lipids,  function(lipid) perform.twas(lipid, "medication"), BPPARAM = MulticoreParam(length(lipids)))
names(twas.medication) <- lipids

save(twas.smoking, twas.medication, file = "Data/twas.sens.rdata")

# perform factor analysis

perform.factor <- function(factors=5) {
  
  z <- bplapply(unique(cov$biobank_id), function(biobank) {
    
    counts <- counts[which(cov$biobank_id == biobank), ]
    cov <- cov[which(cov$biobank_id == biobank), ]
    cov <- cov[, which(!colnames(cov) == "biobank_id")]
    factor.fun(cov[, 1], cov[, -1], counts, factors)
    
  }, BPPARAM=MulticoreParam(length(unique(cov$biobank_id))))
  
  z <- do.call(rbind, z)
  
}

z <- perform.factor()
save(z, file="Data/z.rdata")

# perform polygenic score analysis

perform.ps.str <- function(lipid) {
  
  ps <- lapply(lipids, function(x) {
    
    var <- var[, x]
    dosages <- dosages[, which(gwas$Lipid == x)]
    gwas <- gwas[which(gwas$Lipid == x), ]
    ps.fun(dosages, gwas$beta)
    
  })
  ps <- do.call(cbind, ps)
  colnames(ps) <- lipids
  var <- var[, lipid]
  ps2 <- ps[, which(!colnames(ps) == lipid)]
  ps <- ps[, lipid]
  fit <- lm(ps ~ var + cov$biobank_id)
  res <- data.frame(round(1 - (anova(fit)[nrow(anova(fit)), 2] / (anova(fit)[nrow(anova(fit)), 2] + anova(fit)["var", 2])), 6), round(abs(coef(summary(fit))["var", 3])^2, 6), n=length(var))
  colnames(res) <- c("rsq", "f", "n")
  
  str <- lapply(unique(cov$biobank_id), function(biobank) {
    
    var <- var[which(cov$biobank_id == biobank)]
    ps <- ps[which(cov$biobank_id == biobank)]
    ps2 <- ps2[which(cov$biobank_id == biobank), ]
    bmi.bp.ps <- bmi.bp.ps[which(cov$biobank_id == biobank), ]
    cov <- cov[which(cov$biobank_id == biobank), ]
    cov <- cov[, which(!colnames(cov) == "biobank_id")]
    
    apply(cbind(var, cov, ps2, bmi.bp.ps), 2, function(x) {
      
      coef(summary(lm(ps ~ x)))[2, 1:2]
      
    })
    
    
  })
  str <- do.call(cbind, str)
  
  str_meta <- lapply(unique(colnames(str)), function(name) {
    
    es <- str[1, colnames(str) == name]
    se <- str[2, colnames(str) == name]
    meta.fun(t(data.frame(es)), t(data.frame(se)))
    
  })
  
  names(str_meta) <- unique(colnames(str))
  str_meta <- do.call(cbind, str_meta)
  str_meta1 <- str_meta[, 1:4]
  colnames(str_meta1) <- c("es", "se", "p", "qp")
  str_meta2 <- str_meta[, 5:ncol(str_meta)]
  p <- str_meta2[, grepl(".p", colnames(str_meta2), fixed = T)]
  names(p) <- paste0(unique(colnames(str))[2:length(unique(colnames(str)))], "_p")
  data.frame(str_meta1, p, res)
  
}

ps.str <- bplapply(lipids, perform.ps.str, BPPARAM = MulticoreParam(length(lipids)))
names(ps.str) <- lipids
save(ps.str, file = "Data/ps.str.rdata")

# perform fieller analysis

perform.fieller <- function(lipid, rsids = NA, rsids2 = NA, gi = NA) {
  
  ps.str <- ps.str[[lipid]]
  twas <- twas[[lipid]]
  
  twas.ps <- bplapply(unique(cov$biobank_id), function(biobank) {
    
    var <- var[which(cov$biobank_id == biobank), lipid]
    counts <- counts[which(cov$biobank_id == biobank), ]
    dosages <- dosages[which(cov$biobank_id == biobank), which(gwas$Lipid == lipid)]
    gwas <- gwas[which(gwas$Lipid == lipid), ]
    z <- z[which(cov$biobank_id == biobank), ]
    ps <- ps.fun(dosages, gwas$beta)
    if (!is.na(rsids)) ps <- ps.fun(dosages[, which(!colnames(dosages) %in% rsids)], gwas$beta[which(!gwas$rsid %in% rsids)])
    if (!is.na(gi)) gi <- gi[which(cov$biobank_id == biobank), ]
    cov <- cov[which(cov$biobank_id == biobank), ]
    cov <- cov[, which(!colnames(cov) == "biobank_id")]
    if (!is.na(rsids2)) cov <- cbind(cov, dosages[, which(colnames(dosages) == rsids2)])
    if (!is.na(gi)) cov <- cbind(cov, gi)
    res <- twas.fun.lm(ps, cbind(cov, z), counts)
    bacon.fun(res)
    
  }, BPPARAM=MulticoreParam(length(unique(cov$biobank_id))))
  
  names(twas.ps) <- unique(cov$biobank_id)
  twas.ps <- do.call(cbind, twas.ps)
  es <- twas.ps[, grepl("es_bc", colnames(twas.ps))]
  se <- twas.ps[, grepl("se_bc", colnames(twas.ps))]
  twas.ps <- meta.fun(es, se)
  twas.ps <- twas.ps[which(rownames(twas.ps) %in% rownames(twas[twas$padj < 0.05, ])), ]
  
  res <- bplapply(1:nrow(twas.ps), function(i) {
    
    fieller.fun(as.numeric(twas.ps$es[i]), as.numeric(twas.ps$se[i]), as.numeric(ps.str$n), as.numeric(ps.str$es), as.numeric(ps.str$se), as.numeric(ps.str$n))
  
  }, BPPARAM=MulticoreParam(length(unique(cov$biobank_id))))
  
  res <- as.data.frame(do.call(rbind, res))
  rownames(res) <- rownames(twas.ps)
  res
  
}

fieller <- bplapply(lipids,  perform.fieller, BPPARAM=MulticoreParam(length(lipids)))
names(fieller) <- lipids

# perform correction for snps within 1 mb

perform.fieller.cor <- function(lipid) {
  
  fieller <- fieller[[lipid]]
  gwas <- gwas[which(gwas$Lipid == lipid), ]
  fieller <- fieller[which(fieller$p < 0.05), ]
  ranges <- ranges[which(names(ranges) %in% rownames(fieller))]
  snps <- GRanges(seqnames = gsub(":.*$", "", gwas$SNP), ranges = IRanges(start = as.numeric(gsub("^.*:", "", gwas$SNP)), end = as.numeric(gsub("^.*:", "", gwas$SNP)), names = gwas$rsid))
  hits <- findOverlaps(ranges, snps, maxgap = 1000000, ignore.strand = T)
  ranges <- ranges[queryHits(hits)]
  ranges$rsid <- names(snps[subjectHits(hits)])
  res <- lapply(1:length(ranges), function(i) {
    
    fieller <- perform.fieller(lipid, rsids2 = ranges$rsid[i])
    fieller <- fieller[which(rownames(fieller) == names(ranges)[i]), ]
    
  })
  
  do.call(rbind, res)
  
}

fieller.cor <- bplapply(lipids,  perform.fieller.cor, BPPARAM=MulticoreParam(length(lipids)))
names(fieller.cor) <- lipids
save(fieller.cor, file = "Data/fieller.cor.rdata")

# correction for snps within 1mb

perform.fieller.cor2 <- function(lipid) {
  
  fieller <- fieller[[lipid]]
  gwas <- gwas[which(gwas$Lipid == lipid), ]
  fieller <- fieller[which(fieller$p < 0.05), ]
  ranges <- ranges[names(ranges) %in% rownames(fieller)]
  snps <- GRanges(seqnames = gsub(":.*$", "", gwas$SNP), ranges = IRanges(start=as.numeric(gsub("^.*:", "", gwas$SNP)), end = as.numeric(gsub("^.*:", "", gwas$SNP)), names = gwas$rsid))
  hits <- findOverlaps(ranges, snps, maxgap = 1000000, ignore.strand = T)
  ranges <- ranges[queryHits(hits)]
  ranges$rsid <- names(snps[subjectHits(hits)])
  
  res <- lapply(1:length(ranges), function(i) {
    
    fieller <- perform.fieller(lipid, rsids=ranges$rsid[i])
    fieller <- fieller[which(rownames(fieller) == names(ranges)[i]), ]
    
  })
  
  do.call(rbind, res)
  
}

fieller.cor2 <- bplapply(lipids,  perform.fieller.cor2, BPPARAM = MulticoreParam(length(lipids)))
names(fieller.cor2) <- lipids

# perform twas per snp

perform.twas.snp <- function() {
  
  twas <- bplapply(unique(cov$biobank_id), function(biobank) {
    
    dosages <- dosages[which(cov$biobank_id == biobank), ]
    counts <- counts[which(cov$biobank_id == biobank), ]
    z <- z[which(cov$biobank_id == biobank),]
    cov <- cov[which(cov$biobank_id == biobank), ]
    cov <- cov[, which(!colnames(cov) == "biobank_id")]
    
    res <- bplapply(1:ncol(dosages), function(i) {
      
      res <- twas.fun.lm(dosages[, i], cbind(cov, z), counts)
      bacon.fun(res)
      
    }, BPPARAM = MulticoreParam(5))
    
    do.call(rbind, res)
    
  }, BPPARAM=MulticoreParam(length(unique(cov$biobank_id))))
  
  names(twas) <- unique(cov$biobank_id)
  twas <- do.call(cbind, twas)
  es <- twas[, grepl("es_bc", colnames(twas))]
  se <- twas[, grepl("se_bc", colnames(twas))]
  twas_meta <- meta.fun(es, se)
  twas_meta$rsid <- rep(colnames(dosages), each=ncol(counts))
  twas_meta$transcript <- rep(colnames(counts), ncol(dosages))
  rownames(twas_meta) <- NULL
  twas_meta
  
}

twas.snp <- perform.twas.snp()
save(twas.snp, file="Data/twas.snps.rdata")

# perform qstat analysis

perform.qstat.fun <- function(lipid) {
  
  twas.snp <- twas.snp[rep(gwas$Lipid == lipid, each=length(unique(twas.snp$transcript))), ]
  gwas <- gwas[which(gwas$Lipid == lipid), ]
  twas <- twas[[lipid]]
  twas.snp <- twas.snp[which(twas.snp$transcript %in% rownames(twas[which(twas$padj < 0.05), ])), ]
  
  res <- lapply(unique(twas.snp$transcript), function(transcript) {
    
    twas.snp <- twas.snp[which(twas.snp$transcript == transcript), ]
    gwas <- gwas[match(twas.snp$rsid, gwas$rsid), ]
    
    r_input <- data.frame(twas.snp$rsid, gwas$beta, twas.snp$es, gwas$se, twas.snp$se)
    colnames(r_input) <- c("SNP", "beta.exposure", "beta.outcome", "se.exposure", "se.outcome")
    
    es_egger = summary(lm(r_input$beta.outcome~r_input$beta.exposure, weights=r_input$se.outcome^-2))$coef[2,1]
    se_egger = summary(lm(r_input$beta.outcome~r_input$beta.exposure, weights=r_input$se.outcome^-2))$coef[2,2]/
      min(summary(lm(r_input$beta.outcome~r_input$beta.exposure, weights=r_input$se.outcome^-2))$sigma, 1)
    p_egger <- 2*pt(-abs(es_egger/se_egger), nrow(r_input)-2)
    
    es_pleiotropy = summary(lm(r_input$beta.outcome~r_input$beta.exposure, weights=r_input$se.outcome^-2))$coef[1,1]
    se_pleiotropy = summary(lm(r_input$beta.outcome~r_input$beta.exposure, weights=r_input$se.outcome^-2))$coef[1,2]/
      min(summary(lm(r_input$beta.outcome~r_input$beta.exposure, weights=r_input$se.outcome^-2))$sigma, 1)
    p_pleiotropy <- 2*pt(-abs(es_pleiotropy/se_pleiotropy), nrow(r_input)-2)
    
    snps <- NULL
    qp <- NULL
    Total_Q_chi <- 0
    
    while (Total_Q_chi < 0.05) {
      
      r_input <- r_input[!r_input$SNP %in% snps, ]
      Ratios <-r_input[,3]/r_input[,2]
      F <- r_input[,2]^2/r_input[,4]^2
      mf <- mean(F)
      W <-((r_input[,2]^2)/(r_input[,5]^2))
      Wj <-sqrt(W)
      BetaWj <-Ratios*Wj
      IVW.Model <-lm(BetaWj~-1+Wj)
      EstimatesIVW <-summary(lm(IVW.Model))
      IVW.Slope<-EstimatesIVW$coefficients[1]
      IVW.SE<-EstimatesIVW$coefficients[2]
      IVW_CI<-confint(IVW.Model)
      DF<-length(r_input[,1])-1
      Qj<-W*(Ratios-IVW.Slope)^2
      Total_Q<-sum(Qj)
      Total_Q_chi<-pchisq(Total_Q,length(r_input[,2])-1,lower.tail = F)
      W<- ((r_input[,5]^2+(IVW.Slope^2*r_input[,4]^2))/r_input[,2]^2)^-1
      Wj<-sqrt(W)
      BetaWj<-Ratios*Wj
      IVW.Model<-lm(BetaWj~-1+Wj)
      EstimatesIVW<-summary(lm(BetaWj~-1+Wj))
      IVW.Slope<-EstimatesIVW$coefficients[1]
      IVW.SE<-EstimatesIVW$coefficients[2]
      IVW_CI<-confint(IVW.Model)
      Qj<-W*(Ratios-IVW.Slope)^2
      Total_Q<-sum(Qj)
      Total_Q_chi<-pchisq(Total_Q,length(r_input[,2])-1,lower.tail = F)
      qp <- c(qp, Total_Q_chi)
      if (Total_Q_chi < 0.05) snps <- c(snps, as.character(r_input$SNP)[which.max(Qj)])
      
    }
    
    data.frame(qp=paste(qp, collapse=","), rsid=paste(snps, collapse=","), es_egger=es_egger, se_egger=se_egger, p_egger=p_egger, p_pleiotropy=p_pleiotropy)
    
  })
  
  res <- do.call(rbind, res)
  res[which(res == "")] <- NA
  rownames(res) <- unique(twas.snp$transcript)
  res
  
}

qstat <- bplapply(lipids, qstat.fun, BPPARAM = MulticoreParam(length(lipids)))
names(qstat) <- lipids

# perform fieller adjusted for qstat

perform.fieller.qstat <- function(lipid) {
  
  qstat <- qstat[[lipid]]
  qstat <- na.omit(qstat)
  
  res <- lapply(1:nrow(qstat), function(i) {
    
    rsids <- unlist(strsplit(as.character(qstat$rsid[i]), split = ","))
    fieller <- perform.fieller(lipid, rsids)
    fieller <- fieller[which(rownames(fieller) == rownames(qstat)[i]), ]
    
  })
  
  do.call(rbind, res)
  
}

fieller.qstat <- bplapply(lipids,  perform.fieller.qstat, BPPARAM = MulticoreParam(length(lipids)))
names(fieller.qstat) <- lipids

fieller.comb <- lapply(lipids, function(lipid) {
  
  fieller <- fieller[[lipid]]
  qstat <- qstat[[lipid]]
  fieller.qstat <- fieller.qstat[[lipid]]
  colnames(fieller.qstat) <- paste0(colnames(fieller.qstat), ".cor")
  cbind(fieller, qstat[match(rownames(fieller), rownames(qstat)), ], fieller.qstat[match(rownames(fieller), rownames(fieller.qstat)), ])

  })

names(fieller.comb) <- lipids
save(fieller.comb, file="Data/fieller.rdata")

# perform fieller adjusted for ps

perform.fieller.ps <- function(lipid, gi) {
  
  qstat <- qstat[[lipid]]
  qstat <- na.omit(qstat)
  gi <- ps[, gi, drop = F]
  
  res <- lapply(1:nrow(qstat), function(i) {
    
    rsids <- unlist(strsplit(as.character(qstat$rsid[i]), split=","))
    fieller <- perform.fieller(lipid, rsids, gi=gi)
    fieller[which(rownames(fieller) == rownames(qstat)[i]), ]
    
  })
  
  res <- do.call(rbind, res)
  colnames(res) <- paste0(colnames(res), ".cor")
  fieller <- perform.fieller(lipid, gi = gi)
  cbind(fieller, qstat[match(rownames(fieller), rownames(qstat)), ], res[match(rownames(fieller), rownames(res)), ])
  
}

fieller.ps <- lapply(lipids, function(lipid) bplapply(c("BMI", "SBP", "DBP", "tg", "ldl", "hdl"), function(gi) perform.fieller.ps(lipid, gi), BPPARAM=MulticoreParam(3)))
names(fieller.ps) <- lipids

fieller.ps <- lapply(fieller.qstat2, function(x) {
  
  names(x) <- c("BMI", "SBP", "DBP", "tg", "ldl", "hdl")
  x
  
})

save(fieller.ps, file="Data/fieller.cor.ps.rdata")

eqtls <- import("2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt")
freq <- import("2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt")
eqtls <- eqtls[which(eqtls$BonferroniP < 0.05), ]
eqtls <- eqtls[which(eqtls$Gene %in% unlist(lapply(twas, function(x) rownames(x[which(x$padj < 0.05), ])))), ]
eqtls <- eqtls[order(eqtls$Pvalue), ]
eqtls$beta <- eqtls$Zscore / sqrt(2 * freq$AlleleB_all[match(eqtls$SNP, freq$SNP)] * (1 - freq$AlleleB_all[match(eqtls$SNP, freq$SNP)]) * (eqtls$NrSamples + eqtls$Zscore ^ 2))
eqtls$se <- 1 / sqrt(2 * freq$AlleleB_all[match(eqtls$SNP, freq$SNP)] * (1 - freq$AlleleB_all[match(eqtls$SNP, freq$SNP)]) * (eqtls$NrSamples + eqtls$Zscore ^ 2))
eqtls <- eqtls[which(!duplicated(eqtls$Gene)), ]

eqtls.dosages <- GRanges(seqnames = eqtls$SNPChr, ranges = IRanges(start = eqtls$SNPPos, end = eqtls$SNPPos))
eqtls.dosages <- getGenotypes(imputation_id = ss$imputation_id[which(ss$uuid %in% rownames(counts))], snps = eqtls.dosages, type = "HRC", geno = "DS", BASE = VM_BASE_DATA, BPPARAM = MulticoreParam(15))
colnames(eqtls.dosages) <- ss$uuid[match(colnames(eqtls.dosages), ss$imputation_id)]

eqtls$SNP2 <- paste(eqtls$SNPChr, eqtls$SNPPos, sep = ":")
eqtls.dosages <- eqtls.dosages[, match(rownames(counts), colnames(eqtls.dosages))]
eqtls.dosages <- eqtls.dosages[na.omit(match(eqtls$SNP2, rownames(eqtls.dosages))), ]
eqtls <- eqtls[na.omit(match(rownames(eqtls.dosages), eqtls$SNP2)), ]

# perform fieller analysis for reverse directions

perform.fieller.reverse <- function(lipid) {
  
  twas <- twas[[lipid]]
  var <- var[, lipid]
  eqtls.dosages <- eqtls.dosages[which(eqtls$Gene %in% rownames(twas[which(twas$padj < 0.05), ])), ]
  eqtls <- eqtls[which(eqtls$Gene %in% rownames(twas[which(twas$padj < 0.05), ])), ]

  metrics <- lapply(1:nrow(eqtls), function(i) {
    
    eqtl <- eqtls.dosages[i, ]
    fit <- lm(counts[, which(colnames(counts) == eqtls$Gene[i])] ~ eqtl + cov$biobank_id)
    data.frame(round(1 - (anova(fit)[nrow(anova(fit)), 2] / (anova(fit)[nrow(anova(fit)), 2] + anova(fit)["eqtl", 2])), 6), round(abs(coef(summary(fit))["eqtl", 3])^2, 6), n=length(eqtl))
  
  })
  metrics <- do.call(rbind, metrics)
  colnames(metrics) <- c("rsq", "f", "n")
  rownames(metrics) <- eqtls$Gene

  fieller <- bplapply(unique(cov$biobank_id), function(biobank) {
    
    var <- var[which(cov$biobank_id == biobank)]
    counts <- counts[which(cov$biobank_id == biobank), ]
    eqtls.dosages <- eqtls.dosages[, which(cov$biobank_id == biobank)]
    z <- z[which(cov$biobank_id == biobank), ]
    cov <- cov[which(cov$biobank_id == biobank), ]
    cov <- cov[, which(!colnames(cov) == "biobank_id")]
    
    res <- lapply(1:nrow(eqtls), function(i) {
      
      fit1 <- twas.fun.lm(eqtls.dosages[i, ], cbind(cov, z), counts)
      fit1 <- bacon.fun(fit1)[, 1:2]
      fit1 <- fit1[which(rownames(fit1) == eqtls$Gene[i]), ]
      fit2 <- coef(summary(lm(eqtls.dosages[i, ] ~ var)))[2, 1:2]
      res <- data.frame(fit1[1], fit1[2], fit2[1], fit2[2])
      colnames(res) <- paste(rep(c("fit1", "fit2"), each = 2), c("es", "se"), sep="_")
      res
      
    })
    do.call(rbind, res)
    
  }, BPPARAM=MulticoreParam(length(unique(cov$biobank_id))))
  
  names(fieller) <- unique(cov$biobank_id)
  fieller <- do.call(cbind, fieller)
  rownames(fieller) <- eqtls$Gene
  
  es <- fieller[, grepl("fit1_es", colnames(fieller))]
  se <- fieller[, grepl("fit1_se", colnames(fieller))]
  fieller_meta1 <- meta.fun(es, se)
  
  es <- fieller[, grepl("fit2_es", colnames(fieller))]
  se <- fieller[, grepl("fit2_se", colnames(fieller))]
  fieller_meta2 <- meta.fun(es, se)

  fieller <- bplapply(1:nrow(fieller_meta1), function(i) {
    
    fieller.fun(as.numeric(fieller_meta2[i, 1]), as.numeric(fieller_meta2[i, 2]), as.numeric(metrics$n[i]), as.numeric(fieller_meta1[i, 1]), as.numeric(fieller_meta1[i, 2]), as.numeric(metrics$n[i]))
  
  }, BPPARAM = MulticoreParam(5))
  
  fieller <- do.call(rbind, fieller)
  rownames(fieller) <- rownames(fieller_meta1)
  cbind(fieller, metrics)
  
}

fieller.reverse <- bplapply(lipids,  perform.fieller.reverse, BPPARAM = MulticoreParam(length(lipids)))
names(fieller.reverse) <- lipids
save(fieller.reverse, file = "Data/fieller.reverse.rdata")
