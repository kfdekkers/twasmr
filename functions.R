# factor fun

factor.fun <- function(var, cov, counts, factors = 5) {
  
  design <- data.frame(var, cov)
  design <- model.matrix(~ ., design)
  fit <- cate.fit(design[, 2, drop = F], design[, -2], counts, r = factors, calibrate = F)
  z <- fit$Z
  colnames(z) <- paste0("factor", 1:factors)
  z
  
}

# twas fun (lm)

twas.fun.lm <- function(var, cov, counts) {
  
  design <- cbind(var, cov)
  design <- model.matrix(~ ., design)
  X <- as.matrix(design[, 2, drop = F])
  Z <- as.matrix(design[, -2])
  Y <- as.matrix(counts)
  
  k <- ncol(Z)
  n <- nrow(Y)
  
  U1 <- crossprod(Z, Y)
  U2 <- solve(crossprod(Z), U1)
  ytr <- Y - Z %*% U2
  U3 <- crossprod(Z, X)
  U4 <- solve(crossprod(Z), U3)
  Xtr <- X - Z %*% U4
  Xtr2 <- colSums(Xtr**2)
  
  b <- crossprod(ytr, Xtr)
  b <- b / matrix(Xtr2, nr = nrow(b), nc = ncol(b), byrow = T)
  term1 <- colSums(ytr^2)
  term2 <- matrix(Xtr2, nc = ncol(b), nr = nrow(b), byrow = T) * (b ** 2)
  sig <- (term1 - term2) / (n-k-2)
  err <- sqrt(sig * matrix(1 / Xtr2, nc = ncol(sig), nr = nrow(sig), byrow = T))
  p <- 2 * pt(abs(b/err), n-ncol(Z)-1, lower.tail = F)
  res <- data.frame(es = b, se = err, p = p)
  res
  
}

# bacon fun

bacon.fun <- function(twas) {
  
  bc <- bacon(effectsizes = as.matrix(twas$es), standarderrors = as.matrix(twas$se))
  res <- data.frame(es_bc = es(bc), se_bc = se(bc))
  rownames(res) <- rownames(twas)
  res
  
}

# meta fun

meta.fun <- function(es, se) {
  
  W <- 1 / se ^ 2
  V <- 1 / rowSums(W)
  TS <- rowSums(es * W) * V
  Z <- TS / sqrt(V)
  Q <- rowSums((es - TS) ^ 2 * W)
  Qc <- qchisq(0.95, df = ncol(es) - 1)
  Qp <- pchisq(Q, df = ncol(es) - 1, lower.tail = F)
  data.frame(es = TS, se = sqrt(V), cil = TS - qnorm(0.975) * sqrt(V), ciu = TS + qnorm(0.975) * sqrt(V), p = 2 * pnorm(-abs(Z)), qp = Qp)
  
}

# polygenic score fun

ps.fun <- function(dosages, weights) {
  
  scale(colSums(t(dosages) * weights))
  
}

# fieller fun

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
    res <- c(Q, p, ci)
    
  }
  names(res) <- c("es", "p", "cil", "ciu")
  res
  
}

