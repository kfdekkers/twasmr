# set seed and options

set.seed(1)
options(stringsAsFactors = F)

# load libraries

.libPaths(c(.libPaths(), "/etc/miniconda/lib/R/library/", "/home/kfdekkers/researchdrive/kdekkers/Scratch/library/"))
library(scales)
library(EnsDb.Hsapiens.v86)
library(rio)
library(reshape2)
library(ggplot2)
library(ggforce)
library(ggrepel)
library(GGally)
library(ggpubr)
library(enrichR)

# load 

setwd("/home/kfdekkers/researchdrive/kdekkers/Scratch/MRTWAS/")

# load data

load("Data/data.rdata")
load("Data/twas.rdata")
load("Data/twas.sens.rdata")
load("Data/ps.str.rdata")
load("Data/fieller.rdata")
load("Data/fieller.cor.rdata")
load("Data/fieller.cor.ps.rdata")
load("Data/fieller.reverse.rdata")
lipids <- c("tg", "hdl", "ldl")

# TWAS table

table1 <- lapply(lipids, function(lipid) {
  
  twas <- twas[[lipid]]
  twas$Lipid <- lipid
  twas$se <- abs(twas$se)
  twas$ensembl <- rownames(twas)
  symbol <- select(EnsDb.Hsapiens.v86, keys = twas$ensembl, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
  twas$symbol <- sapply(twas$ensembl, function(ensembl) {
    
    symbol <- symbol$SYMBOL[which(symbol$GENEID == ensembl)]
    paste0(symbol, collapse = ",")
    
  })
  twas[order(twas$p), c("Lipid", "ensembl", "symbol", "es", "se", "p", "padj")]
  
})

table1 <- do.call(rbind, table1)
table1$Lipid <- factor(table1$Lipid)
levels(table1$Lipid) <- c("HDL-C", "LDL-C", "TG")
colnames(table1) <- c("Lipid", "EnsemblID", "Gene", "Estimate", "SE", "P.value", "adj.P.value")
export(table1[which(table1$adj.P.value < 0.05), ], "Results/TableS1.tsv")

sig <- table1[which(table1$adj.P.value < 0.05), ]
table(sig$Lipid)
min(sig$P.value)
lapply(unique(sig$Lipid), function(lipid) {
  
  sum(sig$EnsemblID[sig$Lipid == lipid] %in% sig$EnsemblID[!sig$Lipid == lipid])
  
})
max(sig$P.value)
nrow(sig)

# MR table

table2 <- lapply(lipids, function(lipid) {
  
  fieller <- fieller.comb[[lipid]]
  fieller1 <- fieller[is.na(fieller$rsid), c("es", "se", "cil", "ciu", "p", "rsid", "qp")]
  fieller2 <- fieller[!is.na(fieller$rsid), c("es.cor", "se.cor", "cil.cor", "ciu.cor", "p.cor", "rsid", "qp")]
  colnames(fieller2) <- c("es", "se", "cil", "ciu", "p", "rsid", "qp")
  fieller <- rbind(fieller1, fieller2)
  fieller$se <- abs(fieller$se)
  fieller$Lipid <- lipid
  fieller$ensembl <- rownames(fieller)
  symbol <- select(EnsDb.Hsapiens.v86, keys = fieller$ensembl, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
  fieller$symbol <- sapply(fieller$ensembl, function(ensembl) {
    
    symbol <- symbol$SYMBOL[symbol$GENEID == ensembl]
    paste0(symbol, collapse = ",")
    
  })
  fieller$padj <- p.adjust(fieller$p, method = "BH", n = sum(!is.na(fieller$p)))
  fieller$qp <- sapply(fieller$qp, function(x) {
    x <- unlist(strsplit(x, split = ","))
    x[length(x)]
  })
  fieller[order(fieller$p), c("Lipid", "ensembl", "symbol", "es", "cil", "ciu", "p", "padj", "rsid", "qp")]
  
})

table2 <- do.call(rbind, table2)
table2$Lipid <- factor(table2$Lipid)
levels(table2$Lipid) <- c("HDL-C", "LDL-C", "TG")
colnames(table2) <- c("Lipid", "EnsemblID", "Gene", "Estimate", "CI.lower", "CI.upper", "P.value", "adj.P.value", "Removed.rsID", "Q.P.value")

# Table S3

tables3 <- lapply(ps.str, function(x) x[, c("es", "se", "p", "rsq", "f")])
tables3 <- do.call(rbind, tables3)
tables3 <- tables3[c(1, 3, 2), ]
tables3 <- data.frame(Lipid = c("TG", "HDL-C", "LDL-C"), tables3)
colnames(tables3) <- c("Lipid", "Estimate", "SE", "P-value", "R-squared", "F-statistic")
export(tables3, "Results/TableS3.tsv")

tables4 <- lapply(ps.str, function(x) x[, c("Sampling_Age_p", "Sex_p", "Neut_Perc_p", "Lymph_Perc_p", "Mono_Perc_p", "Eos_Perc_p", "WBC_p", "RBC_p")])
tables4 <- do.call(rbind, tables4)
tables4 <- tables4[c(1, 3, 2), ]
tables4 <- data.frame(Lipid = c("TG", "HDL-C", "LDL-C"), tables4)
colnames(tables4) <- c("Lipid", "Age", "Sex", "Neutrophils", "Lymphocytes", "Monocytes", "Eosinophils", "WBC", "RBC")
export(tables4, "Results/TableS4.tsv")

power.fun <- function(N, alpha, byx, bOLS, R2xz, varx, vary, epower) {
  
  threschi <- qchisq(1 - alpha, 1) 
  f.value <- 1 + N * R2xz / (1 - R2xz)
  con <- (bOLS - byx) * varx 
  vey <- vary - byx * varx * (2 * bOLS - byx)
  
  if (vey < 0) {
    
    data.frame(Error = "Error: Invalid input. The provided parameters result in a negative estimate for variance of the error term in the two-stage least squares model.")
    
  } else {
    
    if (is.na(epower)) {
      
      b2sls <- byx + con / (N * R2xz)
      v2sls <- vey / (N * R2xz * varx)
      NCP <- b2sls^2 / v2sls
      # 2-sided test
      power <- 1 - pchisq(threschi, 1, NCP)
      data.frame(Parameter = c("Power", "NCP", "F-statistic"), Value = c(power, NCP, f.value), Description = c("", "Non-Centrality-Parameter", "The strength of the instrument"))    
      
    } else {
      
      # Calculation of sample size given power
      z1 <- qnorm(1 - alpha / 2)
      z2 <- qnorm(epower)
      Z  <- (z1 + z2)^2
      # Solve quadratic equation in N
      a <- (byx * R2xz)^2
      b <- R2xz * (2 * byx * con - Z * vey / varx)
      c <- con^2
      N1 <- ceiling((-b + sqrt(b^2 - 4 * a * c)) / (2 * a))
      data.frame(Parameter = "Sample Size", Value = N1)
      
    }
  }
}

perform.power <- function(lipid) {
  
  lipid2 <- "TG"
  if (lipid == "hdl") lipid2 <- "HDL-C"
  if (lipid == "ldl") lipid2 <- "LDL-C"
  
  table1 <- table1[which(table1$Lipid == lipid2), ]
  table2 <- table2[which(table2$Lipid == lipid2), ]
  table1 <- table1[which(table1$adj.P.value < 0.05), ]
  
  table1 <- table1[match(table2$EnsemblID, table1$EnsemblID), ]
  counts <- counts[, match(table2$EnsemblID, colnames(counts))]
  ps.str <- ps.str[[lipid]]
  
  power <- sapply(1:nrow(table2), function(i) {
    round(power.fun(ps.str$n, 0.05, table2$Estimate[i], table1$Estimate[i], ps.str$rsq, var(var[, lipid], na.rm = T), var(counts[, i], na.rm = T), epower = NA)[1, 2], 6)
  })
  n <- sapply(1:nrow(table2), function(i) {
    round(power.fun(NA, 0.05, table2$Estimate[i], table1$Estimate[i], ps.str$rsq, var(var[, lipid]), var(counts[, i]), epower = 0.8)[1, 2], 6)
  })
  table2$Power <- power
  table2$N.power <- n
  table2
  
}

table2 <- do.call(rbind, lapply(lipids, function(lipid) perform.power(lipid)))

sig <- table2[which(table2$adj.P.value < 0.05), ]
table(sig$Lipid)

# Figure S2

data <- table2
data$Identified <- sapply(data$N.power, function(n) sum(data$N.power <= n) / length(data$N.power))
data <- data[which(!round(data$Identified, 2) == 1), ]
ggplot(data, aes(x = N.power, y = Identified)) + geom_density(stat = "identity", fill = "grey", color = "grey") + scale_x_log10(expand = c(0, 0), breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) + scale_y_continuous(expand = expansion(mult = c(0, .01)), labels = scales::percent, breaks = c(0, 0.25, 0.5, 0.75, 1)) + theme_minimal() + theme(axis.line = element_line()) + labs(x = "Study size needed", y = "Effects identifyable") + geom_segment(aes(x = 3229, xend = 3229, y = 0.3, yend = 0.1), arrow = arrow(length = unit(0.2, "cm"), type = "closed", angle = 15)) + expand_limits(x = 1)
ggsave("Results/FigureS2.png", height = 80, width = 174, units = "mm")
ggsave("Results/FigureS2.pdf", height = 80, width = 174, units = "mm")

remove.opposite <- data.frame(table1, table2[match(paste(table1$Lipid, table1$EnsemblID), paste(table2$Lipid, table2$EnsemblID)), ])
remove.opposite <- remove.opposite[which(remove.opposite$adj.P.value.1 < 0.05 & sign(remove.opposite$Estimate) != sign(remove.opposite$Estimate.1)), ]
tables5 <- table2[which(!paste0(table2$Lipid, table2$EnsemblID) %in% remove & table2$adj.P.value < 0.05), ]
table(tables5$Lipid, is.na(tables5$rsID))

table2 <- table2[which(!paste0(table2$Lipid, table2$EnsemblID) %in% remove), ]

# Table S6

table3 <- lapply(lipids, function(lipid) {
  
  fieller <- fieller.reverse[[lipid]]
  fieller <- fieller[, c("es", "se", "p")]
  fieller$se <- abs(fieller$se)
  fieller$Lipid <- lipid
  fieller$ensembl <- rownames(fieller)
  symbol <- select(EnsDb.Hsapiens.v86, keys = fieller$ensembl, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
  fieller$symbol <- sapply(fieller$ensembl, function(ensembl) {
    
    symbol <- symbol$SYMBOL[symbol$GENEID == ensembl]
    paste0(symbol, collapse = ",")
    
  })
  fieller$padj <- p.adjust(fieller$p, method = "BH", n = sum(!is.na(fieller$p)))
  fieller[order(fieller$p), c("Lipid", "ensembl", "symbol", "es", "se", "p", "padj")]
  
})

table3 <- do.call(rbind, table3)
table3$Lipid <- factor(table3$Lipid)
levels(table3$Lipid) <- c("HDL-C", "LDL-C", "TG")
colnames(table3) <- c("Lipid", "EnsemblID", "Gene", "Estimate", "SE", "P.value", "adj.P.value")
table3 <- table3[paste0(table3$Lipid, table3$EnsemblID) %in% paste0(tables5$Lipid, tables5$EnsemblID), ]
table3$adj.P.value <- unlist(lapply(unique(table3$Lipid), function(lipid) {
  
  table3 <- table3[table3$Lipid == lipid, ]
  p.adjust(table3$P.value, method = "BH", n = sum(!is.na(table3$P.value)))
  
}))
export(table3, "Results/TableS6.tsv")
table(table3$adj.P.value < 0.05)

# Figure 1A

data <- table1
data$Lipid <- as.character(data$Lipid)
data$Lipid[data$adj.P.value >= 0.05] <- "none"
data$Lipid <- factor(data$Lipid, levels = c("TG", "HDL-C", "LDL-C", "none"))
data <- rbind(data[which(data$Lipid == "none"), ], data[which(data$Lipid == "TG"), ], data[which(data$Lipid == "HDL-C"), ], data[which(data$Lipid == "LDL-C"), ])
plot1a <- ggplot(data, aes(x = Estimate, y = -log10(P.value), color = Lipid)) + geom_point(size = 1) + 
  labs(x = "Effect size", y = "-log10(P-value)") + xlim(-0.5, 0.5) +
  scale_y_continuous(expand = c(0, 0), breaks = c(0, 25, 50, 75, 100), limits = c(0, 110)) + 
  geom_hline(yintercept = min(-log10(data$P.value[data$Lipid != "none"])), linetype = 5) + theme_minimal() +
  theme(text = element_text(size = 10), axis.line.x = element_line()) + geom_vline(xintercept = 0) +
  scale_color_manual(values = c(muted("blue"), muted("red"), muted("yellow", 80, 80), "lightgrey"), breaks = c("TG", "HDL-C", "LDL-C"))

# Figure 1B

data <- table1[table1$adj.P.value < 0.05, ]
data$Lipid <- factor(data$Lipid, levels = c("TG", "HDL-C", "LDL-C"))
table <- table(data$EnsemblID, data$Lipid)
table <- apply(table, 1, function(x) paste(names(x[which(x > 0)]), collapse = ", "))
data$color <- table[match(data$EnsemblID, names(table))]
data$color[data$color == "TG"] <- "TG only"
data$color[data$color == "HDL-C"] <- "HDL-C only"
data$color[data$color == "LDL-C"] <- "LDL-C only"
data$color <- factor(data$color, levels = c("TG only", "HDL-C only", "LDL-C only", "TG, HDL-C", "TG, LDL-C", "HDL-C, LDL-C", "TG, HDL-C, LDL-C"))
plot1b <- ggplot(data, aes(Lipid, fill = color)) + geom_bar() + theme_minimal() + 
  theme(axis.line.x = element_line(), text = element_text(size = 10), panel.grid.major.x = element_blank()) + scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = c(muted("blue"), muted("red"), muted("yellow", 80, 80), muted("purple"), muted("green"), muted("orange", 60, 75), "tan4")) + 
  labs(fill = "Overlap", y = "Count")
ggsave("Results/Figure1.pdf", ggarrange(plot1a, plot1b, labels = c("A", "B"), nrow = 1, align = "v", widths = c(1.1, 1)), height = 80, width = 174, units = "mm")
ggsave("Results/Figure1.png", ggarrange(plot1a, plot1b, labels = c("A", "B"), nrow = 1, align = "v", widths = c(1.1, 1)), height = 80, width = 174, units = "mm", dpi = 600)

# Figure 2A

data <- data.frame(table1, table2[match(paste(table1$Lipid, table1$Ensembl), paste(table2$Lipid, table2$Ensembl)), ])
data$Lipid <- as.character(data$Lipid)
data$Lipid[which(data$adj.P.value.1 >= 0.05)] <- "none"
data$Lipid[which(paste0(data$Lipid, data$EnsemblID) %in% remove)] <- "none"
data$Lipid <- factor(data$Lipid, levels = c("TG", "HDL-C", "LDL-C", "none"))
data$Gene[which(data$Gene == "")] <- data$EnsemblID[which(data$Gene == "")]
data$Gene[which(data$adj.P.value.1 >= 0.05)] <- ""
data <- data[order(data$P.value.1), ]
data$Gene[which(data$Lipid != "TG" | !data$Gene %in% c("ABCA1", "SQLE", "HPGDS", "CYP11A1", "ACSL6", "SLC27A2", "SREBF2", "ABCG1", "TNFRSF21", "PLD3", "IL4", "IL1RL1", "HDC", "HRH4", "FCER1A", "MS4A2", "NTRK1", "PTGER3"))] <- ""
# data <- do.call(rbind, data)
data <- rbind(data[which(data$Lipid == "none"), ], data[which(data$Lipid == "TG"), ], data[which(data$Lipid == "HDL-C"), ], data[which(data$Lipid == "LDL-C"), ])
plot2a <- ggplot(data, aes(x = Estimate, y = Estimate.1, color = Lipid, label = Gene)) + geom_point(size = 0.5) + 
  labs(x = "Association effect size", y = "Causal effect size") + xlim(-0.6, 0.6) + ylim(-0.6, 0.6) +
  coord_equal() + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + theme_minimal() + 
  theme(text = element_text(size = 10)) + geom_abline(slope = 1, intercept = 0, linetype = 5) +
  scale_color_manual(values = c(muted("blue"), muted("red"), muted("yellow", 80, 80), "lightgrey"), breaks = c("TG", "HDL-C")) + geom_text_repel(size = 1.5, show.legend = F)
plot2a <- ggplot(data, aes(x = Estimate, y = Estimate.1, color = Lipid, label = Gene)) + geom_point(size = 1) + 
  labs(x = "Association effect size", y = "Causal effect size") + xlim(-0.6, 0.6) + ylim(-0.6, 0.6) +
  coord_equal() + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + theme_minimal() + 
  theme(text = element_text(size = 10)) + geom_abline(slope = 1, intercept = 0, linetype = 5) +
  scale_color_manual(values = c(muted("blue"), muted("red"), muted("yellow", 80, 80), "lightgrey"), breaks = c("TG", "HDL-C"))

# Figure 2B

data <- table2[which(table2$adj.P.value < 0.05 & !paste0(table2$Lipid, table2$EnsemblID) %in% remove), ]
data$Lipid <- factor(data$Lipid, levels = c("TG", "HDL-C", "LDL-C"))
table <- table(data$EnsemblID, data$Lipid)
table <- apply(table, 1, function(x) paste(names(x[which(x > 0)]), collapse = ", "))
data$color <- as.character(table[match(data$EnsemblID, names(table))])
data$color[which(data$color == "TG")] <- "TG only"
data$color <- factor(data$color, levels = c("TG only", "TG, HDL-C"))
plot2b <- ggplot(data, aes(Lipid, fill = color)) + geom_bar() + theme_minimal() + 
  theme(axis.line.x = element_line(), text = element_text(size = 10), panel.grid.major.x = element_blank()) + scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = c(muted("blue"), muted("purple"))) + 
  labs(fill = "Overlap", y = "Count")
ggsave("Results/Figure2.pdf", ggarrange(plot2a, plot2b, labels = c("A", "B"), nrow = 1, align = "v", widths = c(1.5, 1)), height = 80, width = 174, units = "mm")
ggsave("Results/Figure2.png", ggarrange(plot2a, plot2b, labels = c("A", "B"), nrow = 1, align = "v", widths = c(1.5, 1)), height = 80, width = 174, units = "mm", dpi = 600)

# Figure S1

lipids2 <- c("TG", "LDL-C", "HDL-C")
base <- lapply(1:3, function(i) {
  twas <- twas[[i]][c("es", "padj")]
  twas$Lipid <- lipids2[i]
  twas
})
base <- do.call(rbind, base)

smoking <- lapply(1:3, function(i) {
  twas.smoking <- twas.smoking[[i]][c("es", "padj")]
  twas.smoking$Lipid <- lipids2[i]
  twas.smoking
})
smoking <- do.call(rbind, smoking)

medication <- lapply(1:3, function(i) {
  twas.medication <- twas.medication[[i]][c("es", "padj")]
  twas.medication$Lipid <- lipids2[i]
  twas.medication
})
medication <- do.call(rbind, medication)

colnames(base) <- paste("Base", colnames(base), sep = ".")
colnames(smoking) <- paste("Smoking", colnames(smoking), sep = ".")
colnames(medication) <- paste("Medication", colnames(medication), sep = ".")

data <- cbind(base, smoking, medication)
data <- melt(data[, c("Base.padj", "Base.Lipid", grep("es", colnames(data), value = T))], id = c("Base.Lipid", "Base.es", "Base.padj"))
levels(data$variable) <- c("Smoking behavior", "Lipid-lowering medication use")
data$Base.Lipid[data$Base.padj >= 0.05] <- "none"
data$Base.Lipid <- factor(data$Base.Lipid, levels = c("TG", "HDL-C", "LDL-C", "none"))
data <- rbind(data[data$Base.Lipid == "none", ], data[data$Base.Lipid == "TG", ], data[data$Base.Lipid == "HDL-C", ], data[data$Base.Lipid == "LDL-C", ])
ggplot(data, aes(x = Base.es, y = value, color = Base.Lipid)) + geom_point() + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  geom_abline(linetype = "dashed") + labs(x = "Base model effect size", y = "Adjusted model effect size", color = "Lipid") + xlim(-0.5, 0.5) +
  theme_minimal() + theme(text = element_text(size = 10)) + facet_wrap(~variable, scales = "free") +  ylim(-0.5, 0.5) +
  scale_color_manual(values = c(muted("blue"), muted("red"), muted("yellow", 80, 80), "lightgrey"), breaks = c("TG", "HDL-C", "LDL-C"))
ggsave("Results/FigureS1.pdf", height = 80, width = 174, units = "mm")
ggsave("Results/FigureS1.png", height = 80, width = 174, units = "mm", dpi = 600)

# Figure S3

none <- lapply(lipids, function(lipid) {
  
  fieller <- fieller.comb[[lipid]]
  fieller$Lipid <- lipid
  fieller$ensembl <- rownames(fieller)
  fieller[order(fieller$p), c("Lipid", "ensembl", "es", "cil", "ciu", "p")]
  
})
none <- do.call(rbind, none)
none$se <- abs(none$se)
none$Lipid <- factor(none$Lipid)
levels(none$Lipid) <- c("HDL-C", "LDL-C", "TG")
colnames(none) <- c("Lipid", "EnsemblID", "Estimate", "CI.lower", "CI.upper", "P.value")
none$Adjustment <- "None"

cor <- lapply(lipids, function(lipid) {
  
  fieller <- fieller.cor[[lipid]]
  fieller$Lipid <- lipid
  fieller$ensembl <- rownames(fieller)
  fieller[order(fieller$p), c("Lipid", "ensembl", "es", "cil", "ciu", "p")]
  
})
cor <- do.call(rbind, cor)
cor$Lipid <- factor(cor$Lipid)
levels(cor$Lipid) <- c("HDL-C", "LDL-C", "TG")
colnames(cor) <- c("Lipid", "EnsemblID", "Estimate", "CI.lower", "CI.upper", "P.value")
cor$Adjustment <- "Local variant"

q <- table2[, c("Lipid", "EnsemblID", "Estimate", "CI.lower", "CI.upper", "P.value")]
q$Adjustment <- "Q-stat"

data <- rbind(none, q, cor)
data <- data[paste(data$Lipid, data$EnsemblID) %in% paste(cor$Lipid, cor$EnsemblID), ]

remove.local <- data[which(data$P.value >= 0.05 & data$Adjustment != "Q-stat"), ]

data$Lipid <- factor(data$Lipid, levels = c("TG", "HDL-C", "LDL-C"))
data$Adjustment <- factor(data$Adjustment, levels = rev(c("None", "Q-stat", "Local variant")))
symbol <- select(EnsDb.Hsapiens.v86, keys = data$EnsemblID, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
data$Gene <- sapply(data$EnsemblID, function(ensembl) {
  
  symbol <- symbol$SYMBOL[symbol$GENEID == ensembl]
  paste0(symbol, collapse = ",")
  
})

ggplot(data, aes(x = Gene, y = Estimate, color = Adjustment, ymin = CI.lower, ymax = CI.upper)) + geom_pointrange(position = position_dodge(width = 0.5), size = 0.3) + geom_hline(yintercept = 0) +
  theme_minimal() + theme(axis.line.x = element_line(), text = element_text(size = 10), panel.spacing = unit(1, "lines"), panel.grid.major.y = element_blank()) + facet_col(~ Lipid, scales = "free", space = "free") + coord_flip() + labs(color = "MR method", x = "") + guides(color = guide_legend(reverse = T)) + scale_y_continuous(expand = c(0, 0))
ggsave("Results/FigureS3.pdf", height = 250, width = 174, units = "mm")
ggsave("Results/FigureS3.png", height = 250, width = 174, units = "mm", dpi = 600)

# Figure S4

data <- lapply(c("tg", "hdl", "ldl"), function(lipid) lapply(c("tg", "hdl", "ldl", "BMI", "SBP"), function(adj) {
  
  fieller <- fieller.qstat2[[lipid]][[adj]]
  fieller1 <- fieller[is.na(fieller$rsid), c("es", "cil", "ciu", "p", "rsid", "qp")]
  fieller2 <- fieller[!is.na(fieller$rsid), c("es.cor", "cil.cor", "ciu.cor", "p.cor", "rsid", "qp")]
  colnames(fieller2) <- c("es", "cil", "ciu", "p", "rsid", "qp")
  fieller <- rbind(fieller1, fieller2)
  fieller$Lipid <- lipid
  fieller$Adjustment <- adj
  fieller$ensembl <- rownames(fieller)
  symbol <- select(EnsDb.Hsapiens.v86, keys = fieller$ensembl, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
  fieller$symbol <- sapply(fieller$ensembl, function(ensembl) {
    
    symbol <- symbol$SYMBOL[symbol$GENEID == ensembl]
    paste0(symbol, collapse = ",")
    
  })
  fieller$padj <- p.adjust(fieller$p, method = "BH", n = sum(!is.na(fieller$p)))
  fieller$qp <- sapply(fieller$qp, function(x) {
    x <- unlist(strsplit(x, split = ","))
    x[length(x)]
  })
  fieller[order(fieller$p), c("Adjustment", "Lipid", "ensembl", "symbol", "es", "cil", "ciu", "p", "padj", "rsid", "qp")]
  
}))
data <- do.call(rbind, do.call(rbind, data))
data$Lipid <- factor(data$Lipid)
levels(data$Lipid) <- c("HDL-C", "LDL-C", "TG")
data$Adjustment <- factor(data$Adjustment)
levels(data$Adjustment) <- c("BMI", "HDL-C", "LDL-C", "Blood pressure", "TG")
colnames(data) <- c("Adjustment", "Lipid", "EnsemblID", "Gene", "Estimate", "CI.lower", "CI.upper", "P.value", "adj.P.value", "rsID", "Q.stat")

egger <- lapply(lipids, function(lipid) {
  
  fieller <- fieller.comb[[lipid]]
  fieller[which(fieller$p_pleiotropy >= 0.05), ] <- rep(NA, ncol(fieller))
  fieller <- fieller[, c("es_egger", "se_egger", "p_egger")]
  colnames(fieller) <- c("es", "se", "p")
  fieller$Lipid <- lipid
  fieller$cil <- fieller$es - qnorm(0.975) * fieller$se
  fieller$ciu <- fieller$es + qnorm(0.975) * fieller$se
  fieller$ensembl <- rownames(fieller)
  symbol <- select(EnsDb.Hsapiens.v86, keys = fieller$ensembl, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
  fieller$symbol <- sapply(fieller$ensembl, function(ensembl) {
    
    symbol <- symbol$SYMBOL[symbol$GENEID == ensembl]
    paste0(symbol, collapse = ",")
    
  })
  fieller[order(fieller$p), c("Lipid", "ensembl", "symbol", "es", "cil", "ciu", "p")]
  
})
egger <- do.call(rbind, egger)
egger$Lipid <- factor(egger$Lipid)
levels(egger$Lipid) <- c("HDL-C", "LDL-C", "TG")
colnames(egger) <- c("Lipid", "EnsemblID", "Gene", "Estimate", "CI.lower", "CI.upper", "P.value")
egger$Adjustment <- "Egger"

data <- rbind(data[, c("Lipid", "EnsemblID", "Gene", "Estimate", "CI.lower", "CI.upper", "P.value", "Adjustment")], egger)
data <- data[which(paste(data$Lipid, data$EnsemblID) %in% paste(table2$Lipid, table2$EnsemblID)[which(table2$adj.P.value < 0.05)]), ]
data[which(as.character(data$Adjustment) == as.character(data$Lipid)), ] <- rep(NA, ncol(data))
data <- data[which(!is.na(data$Lipid)), ]

remove.multivariable <- data[which(data$P.value >= 0.05), ]

data$Lipid <- factor(data$Lipid, levels = c("TG", "HDL-C", "LDL-C"))
data$Adjustment <- factor(data$Adjustment, levels = rev(c("TG", "HDL-C", "LDL-C", "BMI", "Blood pressure", "Egger")))
symbol <- select(EnsDb.Hsapiens.v86, keys = data$EnsemblID, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
data$Gene <- sapply(data$EnsemblID, function(ensembl) {
  
  symbol <- symbol$SYMBOL[which(symbol$GENEID == ensembl)]
  paste0(symbol, collapse = ",")
  
})
data$Gene[which(data$Gene == "")] <- data$EnsemblID[which(data$Gene == "")]

ggplot(data, aes(x = Gene, y = Estimate, color = Adjustment, ymin = CI.lower, ymax = CI.upper)) + geom_pointrange(position = position_dodge(width = 0.7), size = 0.3) + geom_hline(yintercept = 0) +
  theme_minimal() + theme(axis.line.x = element_line(), text = element_text(size = 10), panel.spacing = unit(1, "lines"), panel.grid.major.y = element_blank()) + facet_col(~ Lipid, scales = "free", space = "free") + coord_flip() + labs(color = "Adjustment", x = "") + guides(color = guide_legend(reverse = T)) + scale_y_continuous(expand = c(0, 0))
ggsave("Results/FigureS4.pdf", height = 900, width = 174, units = "mm")
ggsave("Results/FigureS4.png", height = 900, width = 174, units = "mm", dpi = 600)

table2$Pleiotropy <- NA
table2$Pleiotropy[match(paste(remove.multivariable$Lipid, remove.multivariable$EnsemblID), paste(table2$Lipid, table2$EnsemblID))] <- as.character(remove.multivariable$Adjustment)
table2$Pleiotropy[which(table2$Pleiotropy == "TG")] <- "TG GIV"
export(table2[which(table2$adj.P.value < 0.05), ], "Results/TableS5.tsv")

remove <- c(paste(remove.opposite$Lipid, remove.opposite$EnsemblID), paste(remove.local$Lipid, remove.local$EnsemblID), paste(remove.multivariable$Lipid, remove.multivariable$EnsemblID))
data <- table2[which(!paste(table2$Lipid, table2$EnsemblID) %in% remove), ]
data <- data[which(data$adj.P.value < 0.05), ]

# Figure 3A

genes <- data
data <- read.delim("GSE107011_Processed_data_TPM.txt")
X <- gsub("[.].*$", "", as.character(data$X))
X2 <- genes$Gene[match(X, genes$EnsemblID)]
celltypes <- sub("_", ":", colnames(data)[2:ncol(data)])
celltypes <- sub("^.*:", "", celltypes)
data <- lapply(unique(celltypes), function(celltype) {
  data <- data[, grepl(celltype, colnames(data))]
  rowMeans(data)
})
data <- as.data.frame(do.call(cbind, data))
colnames(data) <- unique(celltypes)
data <- cbind(X = X, X2 = X2, data)

data <- data[match(genes$EnsemblID[which(genes$Lipid == "TG")], data$X), ]
data <- melt(data, id = c("X", "X2"))

medians <- lapply(unique(data$variable), function(x) {
  
  data <- data[which(data$variable == x), ]
  median(data$value, na.rm = T)
  
})
data$variable <- factor(data$variable, levels = unique(data$variable)[order(unlist(medians))])

lips <- c("ABCA1", "SQLE", "HPGDS", "CYP11A1", "ACSL6", "SLC27A2", "SREBF2", "ABCG1", "PLD3")
asthma <- c("IL4", "IL1RL1", "HDC", "HRH4", "FCER1A", "MS4A2", "CCR3", "HPGDS", "CYP11A1", "PTGER3", "NTRK1")

lips <- data[data$X2 %in% lips, ]
asthma <- data[data$X2 %in% asthma, ]

data$Type <- "All"
lips$Type <- "Lipid metabolism"
asthma$Type <- "Allergy"

data <- rbind(data, lips, asthma)

plot1 <- ggplot(data, aes(x = variable, y = log1p(value), fill = Type)) + geom_boxplot(outlier.shape = NA, coef = 0) + theme_minimal() + 
  theme(axis.line = element_line(), axis.ticks = element_line(), text = element_text(size = 10), axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Cell type", y = "logTPM", fill = "Genes") + scale_y_continuous(expand = c(0, 0))

# Figure 3B

enrich <- enrichr(genes$Gene[which(genes$Lipid == "TG")], c("BioPlanet_2019", "WikiPathways_2019_Human", "KEGG_2019_Human", "Elsevier_Pathway_Collection", "BioCarta_2015", "Reactome_2016", "HumanCyc_2016", "NCI-Nature_2016", "Panther_2016", "MSigDB_Hallmark_2020"))

enrich <- lapply(1:length(enrich), function(i) {
  
  res <- enrich[[i]]
  res$Database <- names(enrich)[i]
  res
  
})
enrich <- do.call(rbind, enrich)
enrich <- enrich[which(enrich$Adjusted.P.value < 0.05), ]
enrich <- enrich[order(enrich$P.value), ]
enrich <- enrich[, !grepl("Old", colnames(enrich))]
export(enrich, "Results/TableS8.tsv")
enrich$Term <- gsub(" Homo.*$", "", enrich$Term)
enrich$Term <- gsub(" WP.*$", "", enrich$Term)
enrich$Database <- gsub("_", " ", enrich$Database)
enrich$Term <- paste0(enrich$Term, "\n(", gsub(";", ", ", enrich$Genes), ")")
enrich$Term[12] <- paste0(" ", enrich$Term[12]) # fix double term
enrich$Term[49] <- paste0(" ", enrich$Term[49]) # fix double term

enrich$Term <- factor(enrich$Term, levels = enrich$Term[order(enrich$P.value, decreasing = T)])
plot2 <- ggplot(enrich[1:15, ], aes(x = -log10(P.value), y = Term, fill = Database)) + geom_col() + theme_minimal() + labs(y = "", x = "-log10(P-value)") + theme(text = element_text(size = 9), axis.line = element_line(), axis.ticks = element_line(), panel.grid = element_blank()) + scale_x_continuous(expand = c(0, 0))

ggsave("Results/Figure3.pdf", ggarrange(plot2, plot1, labels = c("A", "B"), nrow = 2, heights = c(1.5, 1)), height = 200, width = 174, units = "mm")
ggsave("Results/Figure3.png", ggarrange(plot2, plot1, labels = c("A", "B"), nrow = 2, heights = c(1.5, 1)), height = 200, width = 174, units = "mm", dpi = 600)

# Table S9

fieller.fun <- function(betaA, seA, nA, betaB, seB, nB) {
  
  Q <- betaA / betaB
  n <- nA + nB
  g <- (qt(0.975, (n - 2)) * seB / betaB) ^ 2
  
  if (g >= 1) {
    
    res <- rep(NA, 4)
    
  } else {
    
    seQ <- Q / (1 - g) * sqrt((1 - g) * seA ^ 2/betaA ^ 2 + seB ^ 2 / betaB ^ 2)
    ci <- Q/(1 - g) + c(-1, 1) * qt(0.975, (n - 2)) * seQ
    t <- Q / (1 - g) / seQ
    p <- 2 * pt(-abs(t), df = n - 2)
    res <- c(Q, ci, p)
    
  }
  names(res) <- c("es", "cil", "ciu", "p")
  res
  
}

# import data

eqtlgen <- import("2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt")
freq <- import("2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt")

ige <- import("total_log_ige.01-Mar-2010_15-31.imputed.RATIO.gt.0.3.MAF.gt.0.01.csv")
adult <- import("ADULT1_ADULT2_ONSET_ASTHMA.20180716.allchr.assoc.GC.txt")
child <- import("CHILD_ONSET_ASTHMA.20180501.allchr.assoc.GC.txt")
allergy <- import("SHARE-without23andMe.LDSCORE-GC.SE-META.v0.txt")

# clean data

genes <- genes[which(genes$Lipid == "TG"), ]

eqtlgen <- eqtlgen[which(eqtlgen$Gene %in% genes$EnsemblID), ]
eqtlgen <- eqtlgen[order(eqtlgen$Pvalue), ]
eqtlgen$beta <- eqtlgen$Zscore / sqrt(2 * freq$AlleleB_all[match(eqtlgen$SNP, freq$SNP)] * (1-freq$AlleleB_all[match(eqtlgen$SNP, freq$SNP)]) * (eqtlgen$NrSamples + eqtlgen$Zscore^2))
eqtlgen$se <- 1 / sqrt(2 * freq$AlleleB_all[match(eqtlgen$SNP, freq$SNP)] * (1-freq$AlleleB_all[match(eqtlgen$SNP, freq$SNP)]) * (eqtlgen$NrSamples + eqtlgen$Zscore^2))

# genes on ige

exposure <- eqtlgen[which(eqtlgen$SNP %in% ige$SNP), ]
exposure <- exposure[order(exposure$Pvalue), ]
exposure <- exposure[which(!duplicated(exposure$Gene)), ]

outcome <- ige[na.omit(match(exposure$SNP, ige$SNP)), ]
outcome$beta[which(tolower(outcome$MAJOR) == tolower(exposure$AssessedAllele))] <- -outcome$beta[which(tolower(outcome$MAJOR) == tolower(exposure$AssessedAllele))]

ige_fieller_eqtlgen <- lapply(1:nrow(exposure), function(i) {
  
  fieller.fun(outcome$beta[i], outcome$se[i], outcome$N[i], exposure$beta[i], exposure$se[i], exposure$NrSamples[i])
  
})
ige_fieller_eqtlgen <- do.call(rbind, ige_fieller_eqtlgen)
ige_fieller_eqtlgen <- data.frame(EnsemblID = exposure$Gene, Gene = genes$Gene[match(exposure$Gene, genes$EnsemblID)], rsID = outcome$SNP, GWAS = "IgE", Estimate = ige_fieller_eqtlgen[, 1], CI.Lower = ige_fieller_eqtlgen[, 2], CI.Upper = ige_fieller_eqtlgen[, 3], P.value = ige_fieller_eqtlgen[, 4], adj.P.value = p.adjust(ige_fieller_eqtlgen[, 4], method = "BH", n = nrow(exposure)))

# genes on child asthma

exposure <- eqtlgen[which(eqtlgen$SNP %in% child$SNP), ]
exposure <- exposure[order(exposure$Pvalue), ]
exposure <- exposure[which(!duplicated(exposure$Gene)), ]

outcome <- child[na.omit(match(exposure$SNP, child$SNP)), ]
outcome$BETA <- log(outcome$BETA)
outcome$BETA[which(tolower(outcome$ALLELE0) == tolower(exposure$AssessedAllele))] <- -outcome$BETA[which(tolower(outcome$ALLELE0) == tolower(exposure$AssessedAllele))]

child_fieller_eqtlgen <- lapply(1:nrow(exposure), function(i) {
  
  fieller.fun(outcome$BETA[i], outcome$SE[i], outcome$N[i], exposure$beta[i], exposure$se[i], exposure$NrSamples[i])
  
})
child_fieller_eqtlgen <- do.call(rbind, child_fieller_eqtlgen)
child_fieller_eqtlgen <- data.frame(EnsemblID = exposure$Gene, Gene = genes$Gene[match(exposure$Gene, genes$EnsemblID)], rsID = outcome$SNP, GWAS = "Childhood-onset asthma", Estimate = child_fieller_eqtlgen[, 1], CI.Lower = child_fieller_eqtlgen[, 2], CI.Upper = child_fieller_eqtlgen[, 3], P.value = child_fieller_eqtlgen[, 4], adj.P.value = p.adjust(child_fieller_eqtlgen[, 4], method = "BH", n = nrow(exposure)))

# genes on adult asthma

exposure <- eqtlgen[which(eqtlgen$SNP %in% adult$SNP), ]
exposure <- exposure[order(exposure$Pvalue), ]
exposure <- exposure[which(!duplicated(exposure$Gene)), ]

outcome <- adult[na.omit(match(exposure$SNP, adult$SNP)), ]
outcome$BETA <- log(outcome$BETA)
outcome$BETA[which(tolower(outcome$ALLELE0) == tolower(exposure$AssessedAllele))] <- -outcome$BETA[which(tolower(outcome$ALLELE0) == tolower(exposure$AssessedAllele))]

adult_fieller_eqtlgen <- lapply(1:nrow(exposure), function(i) {
  
  fieller.fun(outcome$BETA[i], outcome$SE[i], outcome$N[i], exposure$beta[i], exposure$se[i], exposure$NrSamples[i])
  
})
adult_fieller_eqtlgen <- do.call(rbind, adult_fieller_eqtlgen)
adult_fieller_eqtlgen <- data.frame(EnsemblID = exposure$Gene, Gene = genes$Gene[match(exposure$Gene, genes$EnsemblID)], rsID = outcome$SNP, GWAS = "Adult-onset asthma", Estimate = adult_fieller_eqtlgen[, 1], CI.Lower = adult_fieller_eqtlgen[, 2], CI.Upper = adult_fieller_eqtlgen[, 3], P.value = adult_fieller_eqtlgen[, 4], adj.P.value = p.adjust(adult_fieller_eqtlgen[, 4], method = "BH", n = nrow(exposure)))

fieller_eqtlgen <- rbind(child_fieller_eqtlgen, adult_fieller_eqtlgen, ige_fieller_eqtlgen)
fieller_eqtlgen <- fieller_eqtlgen[fieller_eqtlgen$adj.P.value < 0.05, ]

# genes on allergy

exposure <- eqtlgen[which(eqtlgen$SNP %in% allergy$RS_ID), ]
exposure <- exposure[order(exposure$Pvalue), ]
exposure <- exposure[which(!duplicated(exposure$Gene)), ]

outcome <- allergy[na.omit(match(exposure$SNP, allergy$RS_ID)), ]
outcome$BETA[which(tolower(outcome$OTHER_ALLELE) == tolower(exposure$AssessedAllele))] <- -outcome$BETA[which(tolower(outcome$OTHER_ALLELE) == tolower(exposure$AssessedAllele))]

allergy_fieller_eqtlgen <- lapply(1:nrow(exposure), function(i) {
  
  fieller.fun(outcome$BETA[i], outcome$SE[i], outcome$N[i], exposure$beta[i], exposure$se[i], exposure$NrSamples[i])
  
})
allergy_fieller_eqtlgen <- do.call(rbind, allergy_fieller_eqtlgen)
allergy_fieller_eqtlgen <- data.frame(EnsemblID = exposure$Gene, Gene = genes$Gene[match(exposure$Gene, genes$EnsemblID)], rsID = outcome$RS_ID, GWAS = "Allergy", Estimate = allergy_fieller_eqtlgen[, 1], CI.Lower = allergy_fieller_eqtlgen[, 2], CI.Upper = allergy_fieller_eqtlgen[, 3], P.value = allergy_fieller_eqtlgen[, 4], adj.P.value = p.adjust(allergy_fieller_eqtlgen[, 4], method = "BH", n = nrow(exposure)))

fieller_eqtlgen <- rbind(allergy_fieller_eqtlgen, child_fieller_eqtlgen, adult_fieller_eqtlgen, ige_fieller_eqtlgen)
fieller_eqtlgen <- fieller_eqtlgen[which(fieller_eqtlgen$adj.P.value < 0.05), ]

export(fieller_eqtlgen[order(fieller_eqtlgen$GWAS, fieller_eqtlgen$P.value), ], "Results/TableS9.tsv")

# Table S7

bartel_tg <- import("pgen.1005274.s021.xlsx", which = 6)
bartel_tg <- bartel_tg$mRNA[which(p.adjust(bartel_tg$`p-value`, method = "bonf") < 0.05)]
inouye_tg <- c("HDC", "CPA3", "SLC45A3", "GATA2", "MS4A2", "FCER1A", "SNORD13", "SPRYD5", "C1ORF186", "HIST1H4C", "MS4A3", "C21ORF7", "C1ORF150", "SNORD46", "FBXW8", "RTN2", "MYLIP", "DAPK1", "ABCA1", "RPS27", "ENPP3", "COX7C", "TMEM140", "HSH2D", "RGS18", "TRIM48", "SSR2", "MRPL40", "HS.132563", "CD52", "CEBPD", "SDPR")

bartel_hdl <- import("pgen.1005274.s021.xlsx", which = 2)
bartel_hdl <- bartel_hdl$mRNA[which(p.adjust(bartel_hdl$`p-value`, method = "bonf") < 0.05)]
inouye_hdl <- c("HDC", "SL45A3", "FCER1A")

remove <- c(paste(remove.opposite$Lipid, remove.opposite$EnsemblID), paste(remove.local$Lipid, remove.local$EnsemblID), paste(remove.multivariable$Lipid, remove.multivariable$EnsemblID))
data <- table2[which(!paste(table2$Lipid, table2$EnsemblID) %in% remove), ]
data <- data[which(data$adj.P.value < 0.05), ]

tg <- data.frame(Lipid = "TG", Gene = data$Gene[which(data$Lipid == "TG")], Bartel = data$Gene[which(data$Lipid == "TG")] %in% bartel_tg, Inouye = data$Gene[which(data$Lipid == "TG")] %in% inouye_tg)
hdl <- data.frame(Lipid = "HDL-C", Gene = data$Gene[which(data$Lipid == "HDL-C")], Bartel = data$Gene[which(data$Lipid == "HDL-C")] %in% bartel_hdl, Inouye = data$Gene[which(data$Lipid == "HDL-C")] %in% inouye_hdl)
tables7 <- rbind(tg, hdl)
export(tables7, "Results/TableS7.tsv")

table(tg$Bartel | tg$Inouye)
table(hdl$Bartel | hdl$Inouye)
