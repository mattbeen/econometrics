setwd('/Users/matthew/Documents/Documents - Hong’s MacBook Air/ECON7350 Applied Econometrics')

library(forecast)
library(dplyr)
library(zoo)
library(aTSA)
library(ARDL)

source("ardl_irfs_ci.R")


# create a function to estimate a range of ADF regression specifications
# in levels along with the AICs and BICs
ADF_estimate_lev <- function(y,  p_max = 9){
  TT <- length(y)
  ADF_est <- list()
  ic <- matrix(nrow = 3 * (1 + p_max), ncol = 5)
  colnames(ic) <- c("const", "trend", "p", "aic", "bic")
  i <- 0
  for (const in 0:1)
  {
    for (p in 0:p_max)
    {
      i <- i + 1
      try(silent = T, expr =
            {
              ADF_est[[i]] <- Arima(diff(y), xreg = y[-TT],
                                    order = c(p, 0, 0),
                                    include.mean = as.logical(const),
                                    include.drift = F)
              ic[i,] <- c(const, 0, p, ADF_est[[i]]$aic,
                          ADF_est[[i]]$bic)
            })
    }
    
    if (const)
    {
      # only add a specification with trend if there is a
      # constant (i.e., exclude no constant with trend)
      for (p in 0:p_max)
      {
        i <- i + 1
        try(silent = T, expr =
              {
                ADF_est[[i]] <- Arima(diff(y), xreg = y[-TT],
                                      order = c(p, 0, 0),
                                      include.mean = as.logical(const),
                                      include.drift = T)
                ic[i,] <- c(const, 1, p, ADF_est[[i]]$aic,
                            ADF_est[[i]]$bic)
              })
      }
    }
  }
  
  ic_aic <- ic[order(ic[,4]),][1:10,]
  ic_bic <- ic[order(ic[,5]),][1:10,]
  
  return(list(ADF_est = ADF_est, ic = ic,
              ic_aic = ic_aic, ic_bic = ic_bic))
}

# create a function to estimate a range of ADF regression specifications
# in differences along with the AICs and BICs
ADF_estimate_diff <- function(y, p_max = 9){
  TT <- length(diff(y))
  ADF_est_diff <- list()
  ic_diff <- matrix(nrow = 3 * (1 + p_max), ncol = 5)
  colnames(ic_diff) <- c("const", "trend", "p", "aic", "bic")
  i <- 0
  for (const in 0:1)
  {
    for (p in 0:p_max)
    {
      i <- i + 1
      try(silent = T, expr =
            {
              ADF_est_diff[[i]] <- Arima(diff(diff(y)),
                                         xreg = diff(y)[-TT],
                                         order = c(p, 0, 0),
                                         include.mean = as.logical(const),
                                         include.drift = F)
              ic_diff[i,] <- c(const, 0, p, ADF_est_diff[[i]]$aic,
                               ADF_est_diff[[i]]$bic)
            })
    }
    
    if (const)
    {
      # only add a specification with trend if there is a
      # constant (i.e., exclude no constant with trend)
      for (p in 0:p_max)
      {
        i <- i + 1
        try(silent = T, expr =
              {
                ADF_est_diff[[i]] <- Arima(diff(diff(y)),
                                           xreg = diff(y)[-TT],
                                           order = c(p, 0, 0),
                                           include.mean = as.logical(const),
                                           include.drift = T)
                ic_diff[i,] <- c(const, 1, p, ADF_est_diff[[i]]$aic,
                                 ADF_est_diff[[i]]$bic)
              })
      }
    }
  }
  
  ic_aic_diff <- ic_diff[order(ic_diff[,4]),][1:10,]
  ic_bic_diff <- ic_diff[order(ic_diff[,5]),][1:10,]
  
  return(list(ADF_est_diff = ADF_est_diff,
              ic_diff = ic_diff,
              ic_aic_diff = ic_aic_diff,
              ic_bic_diff = ic_bic_diff))  
}

# Read data 
data <- read.delim("report1.csv", header = TRUE,  sep = ",")

# Extract Variables
dates <- as.yearqtr(data$quarter)
cpi <- data$cpi_inflation
employ <- data$unemployment_rate
cash_rate <- data$cash_rate

# Plot
plot(data$cpi, type = "l", xlab = "Time (Quarters)",
     main = "CPI")
plot(data$unemployment_rate, type = "l", xlab = "Time (Quarters)",
     main = "Unemployment Rate")
plot(data$cash_rate, type = "l", xlab = "Time (Quarters)",
     main = "Cash Target Rate")

# ADF CPI
CPI_ADF_lev <- ADF_estimate_lev(cpi, p_max = 15)
print(CPI_ADF_lev$ic_aic)
print(CPI_ADF_lev$ic_bic)

cpi_adq_set <- as.matrix(arrange(as.data.frame(
  rbind(CPI_ADF_lev$ic_aic[c(2, 3, 6),],
        CPI_ADF_lev$ic_bic[c(1, 3, 4),])),
  const, trend, p))
cpi_adq_idx <- match(data.frame(t(cpi_adq_set[, 1:3])),
                     data.frame(t(CPI_ADF_lev$ic[, 1:3])))

for (i in 1:length(cpi_adq_idx)){
  checkresiduals(CPI_ADF_lev$ADF_est[[cpi_adq_idx[i]]])
}

cpi_adq_set <- as.matrix(arrange(as.data.frame(
  rbind(CPI_ADF_lev$ic_aic[c(2, 6),],
        CPI_ADF_lev$ic_bic[c(1, 3),])),
  const, trend, p))
cpi_adq_idx <- match(data.frame(t(cpi_adq_set[, 1:3])),
                     data.frame(t(CPI_ADF_lev$ic[, 1:3])))

cpi_adq_set
adf.test(cpi, nlag = 15)

# Unemployment rate
unem_ADF_lev <- ADF_estimate_lev(employ, p_max = 15)
print(unem_ADF_lev$ic_aic)
print(unem_ADF_lev$ic_bic)

unem_adq_set <- as.matrix(arrange(as.data.frame(
  rbind(unem_ADF_lev$ic_aic[c(1, 4, 8),],
        unem_ADF_lev$ic_bic[c(2, 6, 7),])),
  const, trend, p))
unem_adq_idx <- match(data.frame(t(unem_adq_set[, 1:3])),
                     data.frame(t(unem_ADF_lev$ic[, 1:3])))

for (i in 1:length(unem_adq_idx)){
  checkresiduals(unem_ADF_lev$ADF_est[[unem_adq_idx[i]]])
}

unem_adq_set
adf.test(employ, nlag = 10)

# diff unem
unem_ADF_diff <- ADF_estimate_diff(employ, p_max = 15)
print(unem_ADF_diff$ic_aic_diff)
print(unem_ADF_diff$ic_bic_diff)

unem_adq_set_diff <- as.matrix(arrange(as.data.frame(
  unem_ADF_diff$ic_bic_diff[c(1, 2, 3, 5),]),
  const, trend, p))
unem_adq_idx_diff <- match(data.frame(
  t(unem_adq_set_diff[, 1:3])),
  data.frame(
    t(unem_ADF_diff$ic_diff[, 1:3])))

for (i in 1:length(unem_adq_idx_diff))
{
  checkresiduals(
    unem_ADF_diff$ADF_est_diff[[unem_adq_idx_diff[i]]])
}

unem_adq_set_diff
adf.test(diff(employ))

# Cash rate
cash_ADF_lev <- ADF_estimate_lev(cash_rate, p_max = 15)
print(cash_ADF_lev$ic_aic)
print(cash_ADF_lev$ic_bic)

cash_adq_set <- as.matrix(arrange(as.data.frame(
  rbind(cash_ADF_lev$ic_bic[c(1, 2, 3,5, 8),])),
  const, trend, p))
cash_adq_idx <- match(data.frame(t(cash_adq_set[, 1:3])),
                      data.frame(t(cash_ADF_lev$ic[, 1:3])))

for (i in 1:length(cash_adq_idx)){
  checkresiduals(cash_ADF_lev$ADF_est[[cash_adq_idx[i]]])
}

cash_adq_set
adf.test(cash_rate, nlag = 10)

# diff cash rate
cash_ADF_diff <- ADF_estimate_diff(cash_rate, p_max = 15)
print(cash_ADF_diff$ic_aic_diff)
print(cash_ADF_diff$ic_bic_diff)

cash_adq_set_diff <- as.matrix(arrange(as.data.frame(
  cash_ADF_diff$ic_bic_diff[c(1, 2, 3, 4),]),
  const, trend, p))
cash_adq_idx_diff <- match(data.frame(
  t(cash_adq_set_diff[, 1:3])),
  data.frame(
    t(cash_ADF_diff$ic_diff[, 1:3])))

for (i in 1:length(cash_adq_idx_diff))
{
  checkresiduals(
    cash_ADF_diff$ADF_est_diff[[cash_adq_idx_diff[i]]])
}

cash_adq_set_diff
adf.test(diff(cash_rate))

# A regression in R is implemented using the lm function.
eg_reg <- lm( cpi_inflation ~ unemployment_rate + cash_rate,data)
summary(eg_reg)
eg_res <- eg_reg$residuals
eg_res

# same approach as in Q1 but with eg_res instead of data
egr_ADF_lev <- ADF_estimate_lev(eg_res, p_max = 15)
print(egr_ADF_lev$ic_aic)
print(egr_ADF_lev$ic_bic)

# We will only consider specifications without a constant or trend since we are
# focusing on residuals.
# AIC: no constant + no trend and lags 10-13;
# BIC: no constant + no trend and lags 0-2,4-6;

rbind(egr_ADF_lev$ic_aic[c(1, 2, 3, 5),],
      egr_ADF_lev$ic_bic[c(1, 2, 3, 4, 5,10),])

egr_adq_set <- as.matrix(arrange(as.data.frame(
  rbind(egr_ADF_lev$ic_aic[c(1, 2, 3, 5),],
        egr_ADF_lev$ic_bic[c(1, 2, 3, 4, 5,10),])),
  const, trend, p))
egr_adq_idx <- match(data.frame(t(egr_adq_set[, 1:3])),
                     data.frame(t(egr_ADF_lev$ic[, 1:3])))

for (i in 1:length(egr_adq_idx))
{
  checkresiduals(egr_ADF_lev$ADF_est[[egr_adq_idx[i]]])
}
# All residuals look OK.
print(egr_adq_set)
# Now, we can use the function coint.test form the aTSA package to implement
# the E-G test for the specs in the adequate set.
eg_test <- matrix(nrow = 1, ncol = 4)
colnames(eg_test) <- rep("", 4)
rownames(eg_test) <- c("No const, no trend")
for (l in 1:4){
  eg_l <- coint.test(cpi, cbind(employ,cash_rate), nlag = l, output = F)
  eg_test[, l] <- eg_l[1, 3]
  colnames(eg_test)[l] <- paste("Lag", l)
}
print(eg_test)


# ARDL
ardl_est <- list()
ic_ardl <- matrix(nrow = 5 * 6 * 6, ncol = 5)
colnames(ic_ardl) <- c("p", "l", "c","aic", "bic")

# cycle through all the possible variants
ptm <- proc.time() # time the search
i <- 0;
for (p in 1:5)
{
  for (l in 0:5)
      {
    for (c in 0:5) {
        i <- i + 1
        ardl_est[[i]] <- ardl(cpi_inflation ~ unemployment_rate + cash_rate,
                              data, order = c(p, l, c))
        ic_ardl[i,] <- c(p, l, c,AIC(ardl_est[[i]]),
                         BIC(ardl_est[[i]]))
        }
      }
}
print(proc.time() - ptm) # display how long it took

ic_aic_ardl <- ic_ardl[order(ic_ardl[,4]),][1:10,]
ic_bic_ardl <- ic_ardl[order(ic_ardl[,5]),][1:10,]

# find the intersection of AIC and BIC preferred sets
ic_int_ardl <- intersect(as.data.frame(ic_aic_ardl),
                         as.data.frame(ic_bic_ardl))

# find the union of AIC and BIC preferred sets
ic_union_ardl <- union(as.data.frame(ic_aic_ardl),
                       as.data.frame(ic_bic_ardl))

adq_set_ardl <- as.matrix(arrange(ic_int_ardl, p, l, c))
adq_idx_ardl <- match(data.frame(t(adq_set_ardl[, 1:3])),
                      data.frame(t(ic_ardl[, 1:3])))

# Residual anaylsis
for (i in 1:length(adq_idx_ardl))
{
  order <- adq_set_ardl[i,1:3]
  acf(as.numeric(ardl_est[[adq_idx_ardl[i]]]$residuals),
      lag.max = 10, xlim = c(1, 10), xaxp = c(1, 10, 1),
      ylim = c(-0.25, 0.25), yaxp = c(-0.15, 0.15, 2),
      main = paste("Residuals ACF for ARDL(",
                   order[1], ", ",
                   order[2], ", ",
                   order[3], ")", sep = ""))
}


# Equilibrium Relationship
for (i in 1:length(adq_idx_ardl))
{
  order <- adq_set_ardl[i,1:3]
  ecm_sr <- recm(ardl_est[[adq_idx_ardl[i]]], case = 2)
  ecm_lrm <- multipliers(ardl_est[[adq_idx_ardl[i]]])
  print(
    paste("ECM of ARDL(",
          order[1], ", ",
          order[2], ", ",
          order[3], ")", sep = "")
  )
  print(summary(ecm_sr))
  
  print(
    paste("LRM of ARDL(",
          order[1], ", ",
          order[2], ", ",
          order[3], ")", sep = "")
  )
  print((ecm_lrm))
  
}



# DY
z68 = qnorm(1 - (1 - .68) / 2)
z95 = qnorm(1 - (1 - .95) / 2)
j <- 1 # select responses to productivity
shock_name <- "CPI Shock"
y_min <- Inf
y_max <- -Inf
irfs_ci <- list()
lrms_ci <- array(dim = c(3, 5, length(adq_idx_ardl)),
                 dimnames = list(c("Constant", "Unemployment Rate","Cash Rate"),
                                 c("Estimate", "lb95", "lb68",
                                   "ub68", "ub95"),
                                 1:length(adq_idx_ardl)))
adj_ci <- matrix(nrow = 5, ncol = length(adq_idx_ardl),
                 dimnames = list(c("Estimate", "lb95", "lb68",
                                   "ub68", "ub95"),
                                 1:length(adq_idx_ardl)))

for (i in 1:length(adq_idx_ardl))
{
  # compute irfs and CIs
  irfs_ci_i <- ardl_irfs_ci(ardl_est[[adq_idx_ardl[i]]], conf = 0.68)
  irfs_ci[[i]] <- cbind(irfs_ci_i$lb[, j],
                        irfs_ci_i$md[, j],
                        irfs_ci_i$ub[, j])
  y_min <- min(y_min, irfs_ci_i$lb[, j])
  y_max <- max(y_max, irfs_ci_i$ub[, j])
  
  # compute lrms and CIs
  lrms_ci_i <- multipliers(ardl_est[[adq_idx_ardl[i]]])[, 1:3]
  lrms_ci[, , i] <- cbind(lrms_ci_i[, 2],
                          lrms_ci_i[, 2] - z95 * lrms_ci_i[, 3],
                          lrms_ci_i[, 2] - z68 * lrms_ci_i[, 3],
                          lrms_ci_i[, 2] + z68 * lrms_ci_i[, 3],
                          lrms_ci_i[, 2] + z95 * lrms_ci_i[, 3])
  
  # compute speed of adjustment and CIs
  adj_ci_i <- summary(recm(ardl_est[[adq_idx_ardl[i]]], case = 2))
  adj_ci[, i] <- c(adj_ci_i$coefficients["ect", 1],
                   adj_ci_i$coefficients["ect", 1]
                   - z95 * adj_ci_i$coefficients["ect", 2],
                   adj_ci_i$coefficients["ect", 1]
                   - z68 * adj_ci_i$coefficients["ect", 2],
                   adj_ci_i$coefficients["ect", 1]
                   + z68 * adj_ci_i$coefficients["ect", 2],
                   adj_ci_i$coefficients["ect", 1]
                   + z95 * adj_ci_i$coefficients["ect", 2])
}

for (i in 1:length(adq_idx_ardl))
{
  order <- adq_set_ardl[i,1:3]
  plot(0:40, irfs_ci[[i]][, 2], type = "l", ylim = c(y_min, y_max),
       ylab = "Impulse Response", xlab = "Horizon",
       main = paste("ARDL(", order[1], ", ", order[2], ", ",
                    order[3],"): Cumulative IRFs to ", shock_name, sep = ""))
  lines(0:40, irfs_ci[[i]][, 1], type = "l", col = "blue")
  lines(0:40, irfs_ci[[i]][, 3], type = "l", col = "blue")
  lines(0:40, numeric(41), type = "l", col = "red")
}

for (j in 1:3)
{
  assign(paste0("lrm_ci_", dimnames(lrms_ci)[[1]][j]), lrms_ci[j,,])
}
