setwd('/Users/matthew/Documents/Documents - Hong’s MacBook Air/ECON7350 Applied Econometrics')
library(forecast)
library(dplyr)
library(zoo)
library(aTSA)

# Read data 
data <- read.delim("report1.csv", header = TRUE,  sep = ",")

# Extract Variables
dates <- as.yearqtr(data$quarter)
cpi <- data$cpi_inflation

# Plot Inflation along time
plot(dates, cpi, type = "l", xlab = "Time (Quarters)",
     main = "Australia CPI Inflation")

# ADF Regression Model
TT <- length(cpi)
ADF_est <- list()
ic <- matrix( nrow = 30, ncol = 5 )
colnames(ic) <- c("cons", "trend", "p", "aic", "bic")
i <- 0
for (const in 0:1)
{
  for (p in 0:9)
  {
    i <- i + 1
    ADF_est[[i]] <- Arima(diff(cpi), xreg = cpi[-TT], # exclude last observation + add variable
                          order = c(p, 0, 0),
                          include.mean = as.logical(const), # include constant
                          include.drift = F)
    ic[i,] <- c(const, 0, p, ADF_est[[i]]$aic,
                ADF_est[[i]]$bic)
  }
  
  if (const)
  {
    # only add a specification with trend if there is a
    # constant (i.e., exclude no constant with trend)
    for (p in 0:9)
    {
      i <- i + 1
      ADF_est[[i]] <- Arima(diff(cpi), xreg = cpi[-TT],
                            order = c(p, 0, 0),
                            include.mean = as.logical(const),
                            include.drift = T) # whether include trend term
      ic[i,] <- c(const, 1, p, ADF_est[[i]]$aic,
                  ADF_est[[i]]$bic)
    }
  }
}

ic_aic <- ic[order(ic[,4]),][1:10,]
ic_bic <- ic[order(ic[,5]),][1:10,]

ic_aic
ic_bic

ic_int <- intersect(as.data.frame(ic_aic),as.data.frame(ic_bic))

ic_int

adq_set <- as.matrix(arrange(as.data.frame(ic_int),
                                   const, trend, p))
adq_idx <- match(data.frame(t(adq_set[, 1:5])),
                       data.frame(t(ic[, 1:5])))

for (i in 1:length(adq_idx)){
  checkresiduals(ADF_est[[adq_idx[i]]])
}

adq_set
adf.test(cpi,nlag = 10)

TT <- length(diff(cpi))
ADF_est_diff <- list()
ic_diff <- matrix( nrow = 30, ncol = 5 )
colnames(ic_diff) <- c("cons", "trend", "p", "aic", "bic")
i <- 0
for (const in 0:1)
{
  for (p in 0:9)
  {
    i <- i + 1
    ADF_est_diff[[i]] <- Arima(diff(diff(cpi)),
                               xreg = diff(cpi)[-TT],
                               order = c(p, 0, 0),
                               include.mean = as.logical(const),
                               include.drift = F)
    ic_diff[i,] <- c(const, 0, p, ADF_est_diff[[i]]$aic,
                     ADF_est_diff[[i]]$bic)
  }
  
  if (const)
  {
    # only add a specification with trend if there is a
    # constant (i.e., exclude no constant with trend)
    for (p in 0:9)
    {
      i <- i + 1
      ADF_est_diff[[i]] <- Arima(diff(diff(cpi)),
                                 xreg = diff(cpi)[-TT],
                                 order = c(p, 0, 0),
                                 include.mean = as.logical(const),
                                 include.drift = T)
      ic_diff[i,] <- c(const, 1, p, ADF_est_diff[[i]]$aic,
                       ADF_est_diff[[i]]$bic)
    }
  }
}

ic_aic_diff <- ic_diff[order(ic_diff[,4]),][1:10,]
ic_bic_diff <- ic_diff[order(ic_diff[,5]),][1:10,]


ic_aic_diff
ic_bic_diff

ic_int_diff <- intersect(as.data.frame(ic_aic_diff),as.data.frame(ic_bic_diff))

ic_int_diff

adq_set_diff <- as.matrix(arrange(as.data.frame(ic_int_diff),
                             const, trend, p))

adq_idx_diff <- match(data.frame(t(adq_set_diff[, 1:3])),
                      data.frame(t(ic_diff[, 1:3])))

for (i in 1:length(adq_idx_diff)){
  checkresiduals(ADF_est_diff[[adq_idx_diff[i]]])
}


adq_set_diff
adf.test(diff(cpi),nlag = 10)

# build an adequate set from ARIMAs with a possible linear trend
TT <- length(cpi)
ARIMA_est <- list()
ic_arima <- matrix( nrow = 2 * 3 * 4 * 4, ncol = 7 )
colnames(ic_arima) <- c("d", "cons", "trend", "p", "q", "aic", "bic")
i <- 0
for (d in 0:1)
{
  for (const in 0:1)
  {
    for (p in 0:3)
    {
      for (q in 0:3)
      {
        i <- i + 1
        d1 <- as.logical(d)
        c1 <- as.logical(const)
        
        try(silent = T, expr =
              {
                ARIMA_est[[i]] <- Arima(cpi, order = c(p, d, q),
                                        include.constant = c1)
                
                ic_arima[i,] <- c(d, const, 0, p, q,
                                  ARIMA_est[[i]]$aic,
                                  ARIMA_est[[i]]$bic)
              })
        
        if (const)
        {
          # only add a specification with trend if there is a
          # constant (i.e., exclude no constant with trend)
          i <- i + 1
          
          if (d1)
          {
            x <- c(0,cumsum(1:(TT - 1)))
          }
          else
          {
            x <- NULL
          }
          
          try(silent = T, expr =
                {
                  ARIMA_est[[i]] <- Arima(cpi, order = c(p, d, q),
                                          xreg = x,
                                          include.constant = c1,
                                          include.drift = T)
                  
                  ic_arima[i,] <- c(d, const, 1, p, q,
                                    ARIMA_est[[i]]$aic,
                                    ARIMA_est[[i]]$bic)
                })
        }
      }
    }
  }
}
ic_aic_arima <- ic_arima[order(ic_arima[,6]),][1:10,]
ic_bic_arima <- ic_arima[order(ic_arima[,7]),][1:10,]

# find the intersection of AIC and BIC preferred sets
ic_int_arima <- intersect(as.data.frame(ic_aic_arima),
                          as.data.frame(ic_bic_arima))

adq_set_arima <- as.matrix(arrange(as.data.frame(ic_int_arima),
  d, const, trend, p))
adq_idx_arima <- match(data.frame(t(adq_set_arima[, 1:5])),
                       data.frame(t(ic_arima[, 1:5])))

# Check the residuals for specifications in the adequate set.
for (i in 1:length(adq_idx_arima))
{
  checkresiduals(ARIMA_est[[adq_idx_arima[i]]])
}

# White noise residuals are rejected at 99% confidence interval for ARIMA(2,1,3)
# Remove ARIMA(2,1,3) from the adequate set

adq_set_arima <- adq_set_arima[1:8,]
adq_idx_arima <- adq_idx_arima[1:8]

# Do the forecasting
pre2021 <- as.Date(dates) <= as.Date("2021-12-31")
date_fcst <- append(data$quarter,c("2022Q1","2022Q2","2022Q3","2022Q4","2023Q1","2023Q2","2023Q3","2023Q4"))
dates_fcst <- as.yearqtr(date_fcst)
hrz <- 8
xticks <- c(TT - 3 * hrz + c(1, 2 * hrz, 3 * hrz ,4 * hrz))
# xticks <- c(dates[129],dates[182],dates[TT])
# actual_cpi <- as.matrix(cpi[!pre2021])
fcst_cpi <- vector(mode = "list", length(adq_set_arima))
for (i in 1:nrow(adq_set_arima))
{
  model_p_d_q <- adq_set_arima[i, c(4, 1, 5)]
  fcst_cpi[[i]] <- forecast::forecast(
    Arima(cpi[pre2021], model_p_d_q,
          include.constant = adq_set_arima[i,2]),
    h = hrz, level = c(68, 95))
  
  title_p_d_q <- paste("ARIMA(",
                       as.character(model_p_d_q[1]), ", ",
                       as.character(model_p_d_q[2]), ", ",
                       as.character(model_p_d_q[3]), ")",
                       sep = "")
  
  plot(fcst_cpi[[i]], include = hrz * 2, ylab = colnames(cpi),
       main = title_p_d_q, xaxt = "n",
       ylim = c(-1, 10))
  #lines(sum(pre2021) + 1:hrz, actual_cpi)
  lines(1:TT, rep(0, TT), col = "red")
  axis(1, at = xticks, labels = dates_fcst[xticks])
}
