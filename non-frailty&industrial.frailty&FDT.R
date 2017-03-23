library(data.table)
##
## Intensity Model Without Frailty
##
intensity.model <- function(data, varibles, defaults, fdt=FALSE, sort.by, bin=5){
  ## Extract Data
  varibles <- as.matrix(data[,varibles])
  defaults <- data[,defaults]
  sort.by <- data[,sort.by]
  ## Log-Likelihood Function:
  loglike <- function(varibles, par, defaults){
    -sum(dpois(defaults, exp(par[1] + varibles %*% par[2:length(par)]), log = T))
  }
  k <- ncol(varibles)+1
  result <- optim(par = rep(0, times = k), fn = loglike,
                      defaults = defaults, varibles = varibles,
                      #control = list(maxit = 50000, reltol = 1e-10),
                      hessian = TRUE, method = "BFGS")
  ## Result Collection
  if(class(try(solve(result$hessian),silent=T))=="matrix"){
    SE <- sqrt(diag(solve(result$hessian)))
    t.value <- result$par/SE
    bic <- (result$value)*2 + k*(log(length(defaults))-log(2*pi))
    par <- rbind(result$par, t.value)
    colnames(par) <- c("Intercept", colnames(varibles))
    rownames(par) <- c("par", "t.value")
    lnlik <- -result$value
  }
  else {
    SE <- NA
    t.value <- NA
    bic <- (result$value)*2 + k*(log(length(defaults))-log(2*pi))
    par <- rbind(result$par, t.value)
    colnames(par) <- c("Intercept", colnames(varibles))
    rownames(par) <- c("par", "t.value")
    lnlik <- -result$value
    print("WARNING: The Hessian matrix isn't invertible.")
  }
  res <- list(PAR = par,
              Log.Likelihood = lnlik,
              BIC = bic)
  if(fdt){
    intensity <- exp(res$PAR[1, 1] + varibles %*% res$PAR[1, 2:ncol(res$PAR)])
    res$FDT <- lapply(bin, FDT, indicator=defaults, intensity=intensity, sort.by=sort.by)
  }
  res
}

##
## Intensity Model With Industrial Frailty
##
industrial.frailty <- function(data, varibles, defaults, industry, tolerance=0.001, fdt=FALSE, sort.by, bin=5){
  ## Extract Data
  par <- c(as.numeric(intensity.model(data, varibles, defaults)$PAR[1,]), 1)
  varibles <- as.matrix(data[, varibles])
  defaults <- data[, defaults]
  industry <- data[, industry]
  sort.by <- data[,sort.by]
  ## Log-Likelihood Function:
  loglike <- function(dat, par, Di, Ai, Ci){
    theta <- par[length(par)]
    dat[, lambda:=exp(par[1]+varibles%*%par[2:(length(par)-1)])]
    -(sum(((1/theta)-1+Di)*(digamma(Ai)-log(Ci))-(Ai/Ci)/theta)
      -length(unique(dat[,industry]))*(log(gamma(1/theta))+(1/theta)*log(theta))
      +sum(dat[,defaults]*log(dat[,lambda]))-sum((Ai/Ci)*dat[, sum(lambda), by=industry][,V1]))
  }
  dat <- data.table(cbind.data.frame(varibles, defaults, industry, sort.by))
  i = 0
  while(i < length(par)){
    theta <- par[length(par)]
    dat[, lambda:=exp(par[1]+varibles%*%par[2:(length(par)-1)])]
    Di <- dat[, sum(defaults), by=industry][,V1]
    Ai <- (1/theta) + Di
    Ci <- (1/theta) + dat[, sum(lambda), by=industry][,V1]
    result <- optim(par = par, fn = loglike,
                    dat = dat, Di = Di, Ai = Ai, Ci = Ci,
                    hessian = TRUE, method = "BFGS")
    par0 <- par
    par  <- result$par
    i    <- sum(abs(par-par0) < tolerance)
  }
  ## Result Collection
  k <- ncol(varibles)+1
  if(class(try(solve(result$hessian),silent=T))=="matrix"){
    SE <- sqrt(diag(solve(result$hessian)))
    t.value <- result$par/SE
    bic <- (result$value)*2 + k*(log(length(defaults))-log(2*pi))
    par <- rbind(result$par, t.value)
    colnames(par) <- c("Intercept", colnames(varibles), "theta")
    rownames(par) <- c("par", "t.value")
    lnlik <- -result$value
  }
  else{
    SE <- NA
    t.value <- NA
    bic <- (result$value)*2 + k*(log(length(defaults))-log(2*pi))
    par <- rbind(result$par, t.value)
    colnames(par) <- c("Intercept", colnames(varibles), "theta")
    rownames(par) <- c("par", "t.value")
    lnlik <- -result$value
    print("WARNING: The Hessian matrix isn't invertible.")
  }
  res <- list(PAR = par,
              LnLikelihood = lnlik,
              BIC = bic)
  if(fdt){
    frailty <- cbind(dat[, unique(industry), by=industry][,V1], Ai/Ci)
    colnames(frailty) <- c("industry", "industrial.frailty")
    frailty <- data.table(frailty, key="industry")
    setkey(dat, "industry")
    dat[frailty, intensity:=lambda * industrial.frailty]
    res$FDT <- lapply(bin, FDT, indicator=dat[,defaults], intensity=dat[,intensity], sort.by=dat[,sort.by])
  }
  res
}

##
##  Fisher Dispersion Test
##
FDT <- function(indicator, intensity, sort.by, bin){
  data <- cbind.data.frame(indicator, intensity, sort.by)
  data <- data[order(data[,3]),]
  Z <- cumsum(data[,1])
  C <- cumsum(data[,2])
  k <- ceiling(sum(data[,2])/bin)
  a <- seq(from=bin, to=bin*k, by=bin)
  g <- NULL
  for(i in 1:length(a)){
    g[i] <- which.min(abs(C-a[i]))
  }
  if(identical(g[length(g)-1],g[length(g)])){g <- g[1:(length(g)-1)]}
  C <- C[g]
  Z <- Z[g]
  W <- sum(((c(Z[1], diff(Z))-c(C[1], diff(C)))^2)/c(C[1], diff(C)))
  df <- k-1
  p.value <- 1-pchisq(W, df)
  res <- round(c(df, W, p.value), 3)
  names(res) <- c("df", "statistics", "p.value")
  res
}

