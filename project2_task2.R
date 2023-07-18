setwd('/Users/matthew/Documents/Documents - Hong’s MacBook Air/ECON7350 Applied Econometrics')
library(dplyr)
library(zoo)
library(vars)
library(urca)

# Read data 
data <- read.delim("report1.csv", header = TRUE,  sep = ",")
data_new <- read.delim("report2.csv", header = TRUE,  sep = ",")

dates <- as.yearqtr(data$quarter)
cpi <- data$cpi_inflation
gdp <- data$real_gdp
unemploy <- data$unemployment_rate
cash <- data$cash_rate
x <- cbind(cpi, unemploy,cash)

new_cpi <- data_new$cpi_inflation
new_dates <- as.yearqtr(data_new$quarter)
new_gdp <- data_new$real_gdp
new_unemploy <- data_new$unemployment_rate
new_cash <- data_new$cash_rate
new_x <- cbind(new_cpi, new_unemploy,new_cash)

# VAR model
VAR_est <- list()
ic_var <- matrix(nrow = 20,  ncol = 3)
colnames(ic_var) <- c("p", "aic", "bic")
for (p in 1:20){
  VAR_est[[p]] <- VAR(x, p)
  ic_var[p,] <- c(p, AIC(VAR_est[[p]]),
                  BIC(VAR_est[[p]]))
}

ic_aic_var <- ic_var[order(ic_var[,2]),]
ic_bic_var <- ic_var[order(ic_var[,3]),]

ic_aic_var
ic_bic_var

# Adequate set
adq_set_var <- as.matrix(ic_var[5:7,])
adq_idx_var <- c(5:7)

# Residual
nmods <- length(adq_idx_var)
for (i in 1:nmods){
  p <- adq_idx_var[i]
  print(paste0("Checking VAR(", p, ")"))
  print(serial.test(VAR_est[[p]], lags.bg = 1,
                    type = "BG"))
}

# Adjusted Adequate set
adq_set_var <- as.matrix(rbind(ic_var[6,],ic_var[7,]))
adq_idx_var <- c(6,7)

# Cholesky Decomposition
# To investigate the orders and see the responses are sensitive to ordering
orders <- perms(1:3)
vnames <- c("cpi", "unemploy","cash")
for (i in 1:3)
{
  for (j in 1:3)
  {
    for (k in 1:nrow(orders))
    {
      title_i_j_k <- paste0("Response of ",
                            vnames[i],
                            " to a shock in ",
                            vnames[j],
                            "; x = (",
                            paste0(vnames[orders[k,]],
                                   collapse = ", "),
                            ")'")
      
      irf_i_j_k <- irf(VAR(x[,orders[k,]], 3),
                       n.ahead = 40,
                       response = vnames[i],
                       impulse = vnames[j],
                       boot = TRUE)
      
      plot(irf_i_j_k, main = title_i_j_k)
      cat("\r", title_i_j_k, "  ", sep = "")
    }
  }
}

nmods <- length(adq_idx_var)
for (i in 1:nmods){
  p <- adq_idx_var[i]
  print(paste0("VAR(", p, ")"))
  print(summary(ca.jo(x, type = "trace", K = p)))
}

FEVD_est <- fevd(VAR_est[[6]], n.ahead = 40)
plot(FEVD_est, mar = c(2,1,2,1),
     oma = c(0,1,0,1))

n <- 3
VECM_est <- list()
ic_vecm <- matrix(nrow = 4 * (1 + n),  ncol = 4)
colnames(ic_vecm) <- c("p", "r", "aic", "bic")
i <- 0
for (p in 2:4)
{
  for (r in 0:n)
  {
    i <- i + 1
    if (r == n)
    {
      VECM_est[[i]] <- VAR(x, p)
    }
    else if (r == 0)
    {
      VECM_est[[i]] <- VAR(diff(x), p - 1)
    }
    else
    {
      VECM_est[[i]] <- vec2var(ca.jo(x, K = p), r)
    }
    ic_vecm[i,] <- c(p, r, AIC(VECM_est[[i]]),
                     BIC(VECM_est[[i]]))
  }
}

ic_aic_vecm <- ic_vecm[order(ic_vecm[,3]),][1:5,]
ic_bic_vecm <- ic_vecm[order(ic_vecm[,4]),][1:5,]

ic_aic_vecm
ic_bic_vecm

ic_int_vecm <- intersect(as.data.frame(ic_aic_vecm),
                         as.data.frame(ic_bic_vecm))

ic_int_vecm

adq_set_vecm <- as.matrix(arrange(as.data.frame(
  ic_int_vecm), p, r))
adq_idx_vecm <- match(data.frame(t(adq_set_vecm[, 1:2])),
                      data.frame(t(ic_vecm[, 1:2])))

# do a final check of the residuals
nmods <- length(adq_idx_vecm)
for (i in 1:nmods)
{
  p <- adq_set_vecm[i, 1]
  r <- adq_set_vecm[i, 2]
  print(paste0("Checking VECM(", p, ", ", r, ")"))
  print(serial.test(VECM_est[[adq_idx_vecm[i]]],
                    type = "BG"))
}

p <- 6
for (r in 1:3)
{
  if (r < n)
  {
    vec_pr <- ca.jo(x, type = "trace", spec = "transitory")
    beta_pr <- cajorls(vec_pr, r)$beta
    
    bpvals_pr <- beta_pr
    bpvals_pr[,] <- NA
    
    for (i in (r + 1):n)
    {
      for (j in 1:r)
      {
        H <- beta_pr
        H[i, j] <- 0
        bpvals_pr[i, j] <- blrtest(vec_pr, H, r)@pval[1]
      }
    }
    
    print(paste0("Results for VECM(", p, ", ", r, "):"))
    print(beta_pr)
    print(bpvals_pr)
  }
}

vnames <- c("cpi","unemploy","cash")
nmods <- length(adq_idx_vecm)
for (i in 1:3)
{
  for (j in 1:3)
  {
    for (imod in 1:nmods)
    {
      p <- adq_set_vecm[imod, 1]
      r <- adq_set_vecm[imod, 2]
      title_i_j <- paste0("VECM(", p, ", ", r,
                          "): Response of ", vnames[i],
                          " to a shock in ", vnames[j])
      
      irf_i_j <- irf(VECM_est[[adq_idx_vecm[imod]]],
                     n.ahead = 40,
                     response = vnames[i],
                     impulse = vnames[j],
                     boot = TRUE)
      
      plot(irf_i_j, main = title_i_j)
      cat("\r", title_i_j, "  ", sep = "")
    }
  }
}
# new model
data_new <- read.delim("report2.csv", header = TRUE,  sep = ",")
new_cpi <- data_new$cpi_inflation
new_dates <- as.yearqtr(data_new$quarter)
new_unemploy <- data_new$unemployment_rate
new_cash <- data_new$cash_rate
new_x <- cbind(new_cpi,new_unemploy,new_cash)

# VECM
n <- 3
VECM_est_long <- list()
ic_vecm_long <- matrix(nrow = 4 * (1 + n),  ncol = 4)
colnames(ic_vecm_long) <- c("p", "r", "aic", "bic")
i <- 0
for (p in 5:7)
{
  for (r in 0:n)
  {
    i <- i + 1
    if (r == n)
    {
      VECM_est_long[[i]] <- VAR(new_x, p)
    }
    else if (r == 0)
    {
      VECM_est_long[[i]] <- VAR(diff(new_x), p - 1)
    }
    else
    {
      VECM_est_long[[i]] <- vec2var(ca.jo(new_x, K = p), r)
    }
    ic_vecm_long[i,] <- c(p, r, AIC(VECM_est_long[[i]]),
                          BIC(VECM_est_long[[i]]))
  }
}

ic_aic_vecm_long <- ic_vecm_long[order(ic_vecm_long[,3]),][1:10,]
ic_bic_vecm_long <- ic_vecm_long[order(ic_vecm_long[,4]),][1:10,]
ic_int_vecm_long <- intersect(as.data.frame(ic_aic_vecm_long),
                              as.data.frame(ic_bic_vecm_long))

adq_set_vecm_long <- as.matrix(arrange(as.data.frame(
  ic_int_vecm_long), p, r))
adq_idx_vecm_long <- match(data.frame(t(adq_set_vecm_long[, 1:2])),
                           data.frame(t(ic_vecm_long[, 1:2])))

# IRF
vnames <- c("new_cpi","new_unemploy","new_cash")
nmods <- length(adq_idx_vecm_long)
for (i in 1:3)
{
  for (j in 1:3)
  {
    for (imod in 1:nmods)
    {
      p <- adq_set_vecm_long[imod, 1]
      r <- adq_set_vecm_long[imod, 2]
      title_i_j <- paste0("VECM(", p, ", ", r,
                          "): Response of ", vnames[i],
                          " to a shock in ", vnames[j])
      
      irf_i_j <- irf(VECM_est_long[[adq_idx_vecm_long[imod]]],
                     n.ahead = 40,
                     response = vnames[i],
                     impulse = vnames[j],
                     boot = TRUE)
      
      plot(irf_i_j, main = title_i_j)
    }
  }
}

