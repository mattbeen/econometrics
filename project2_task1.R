setwd('/Users/matthew/Documents/Documents - Hong’s MacBook Air/ECON7350 Applied Econometrics')
library(zoo)
library(vars)
library(pracma)

# Read data 
data <- read.delim("report1.csv", header = TRUE,  sep = ",")
data_new <- read.delim("report2.csv", header = TRUE,  sep = ",")

dates <- as.yearqtr(data$quarter)
cpi <- data$cpi_inflation
gdp <- data$real_gdp
unemploy <- data$unemployment_rate
cash <- data$cash_rate
x <- cbind(cpi, gdp, unemploy,cash)

new_cpi <- data_new$cpi_inflation
new_dates <- as.yearqtr(data_new$quarter)
new_gdp <- data_new$real_gdp
new_unemploy <- data_new$unemployment_rate
new_cash <- data_new$cash_rate
new_x <- cbind(new_cpi, new_gdp, new_unemploy,new_cash)


# Generate Graph
plot(dates, cpi, type = "l", xlab = "Time (Quarters)", ylab = "CPI Inflation")
plot(dates, gdp, type = "l", xlab = "Time (Quarters)", ylab = "Real GDP")
plot(dates, unemploy, type = "l", xlab = "Time (Quarters)", ylab = "Unemployment Rate")
plot(dates, cash, type = "l", xlab = "Time (Quarters)", ylab = "Cash Rate")

# We will consider VARs with p= 1,...,20.
# Note: the VAR command does not allow estimating a VAR(0),
# so we will not worry about this specification.
VAR_est <- list()
ic_var <- matrix(nrow = 20,  ncol = 3)
colnames(ic_var) <- c("p", "aic", "bic")
for (p in 1:20)
{
  VAR_est[[p]] <- VAR(x, p)
  ic_var[p,] <- c(p, AIC(VAR_est[[p]]),
                  BIC(VAR_est[[p]]))
}

ic_aic_var <- ic_var[order(ic_var[,2]),]
ic_bic_var <- ic_var[order(ic_var[,3]),]

ic_aic_var
ic_bic_var

# AIC, BIC agree on p = 2,3,4,5
adq_set_var <- as.matrix(ic_var[2:5,])
adq_idx_var <- c(2:5)

# Check the residuals: the vars package provides a few
# different tests through the serial.test function;
# we will use the LM test by setting type = "BG", but
# other tests are just as valid.
nmods <- length(adq_idx_var)
for (i in 1:nmods)
{
  p <- adq_idx_var[i]
  print(paste0("Checking VAR(", p, ")"))
  print(serial.test(VAR_est[[p]], type = "BG"))
}
# Autocorrelations appear to be quite high for three models p=2,3,5.
# Only model with p = 4 is rejected 
# This is a concern, so we proceed by checking all 20 VAR
# specifications.
for (p in 1:20)
{
  print(paste0("Checking VAR(", p, ")"))
  print(serial.test(VAR_est[[p]], type = "BG"))
}
# White noise residuals are rejected only model with p = 4,6
# We will proceed with p = 4, 6 as the adequate set.
adq_set_var <- as.matrix(rbind(ic_var[4,],ic_var[6,]))
adq_idx_var <- c(4,6)

# intercept and slope coefficients
nmods <- length(adq_idx_var)
for (i in 1:nmods)
{
  p <- adq_idx_var[i]
  print(paste0("VAR(", p, ") has ",
               3 * (1 + 3 * p),
               " coefficients."))
}
# The vars package provides a handy function "roots" to
# ascertain the stability of the estimated VARs.
# Unfortunately, it does not provide confidence intervals
# for the estimated roots.
# Also, be careful not to use the function "stability" for
# this purpose -- that function searches for potential
# structural breaks.
nmods <- length(adq_idx_var)
for (i in 1:nmods)
{
  p <- adq_idx_var[i]
  print(paste0("VAR(", p,
               "): Maximum absolute eigenvalue is ",
               max(vars::roots(VAR_est[[p]]))))
}

# Use the predict function to generate forecasts; plot will
# then automatically produce 95% predictive intervals
# (forecast uncertainty only, not estimation uncertainty).
hrz <- 8
VAR_fcst <- list()
xlim <- c(length(dates) - 3 * hrz,
          length(dates) + hrz)
ylim <- c(-1,max(cpi) + 1)
nmods <- length(adq_idx_var)
for (i in 1:nmods)
{
  p <- adq_idx_var[i]
  VAR_fcst[[i]] <- predict(VAR_est[[p]], n.ahead = hrz)
  title_p <- paste("Forecast for CPI Inflation of VAR(",
                   as.character(p), ")",
                   sep = "")
  plot(VAR_fcst[[i]], names = "cpi",
       xlim = xlim, ylim = ylim,
       main = title_p,
       xlab = "Horizon",
       ylab = "CPI")
}

# All period data
VAR_est <- list()
ic_var <- matrix(nrow = 20,  ncol = 3)
colnames(ic_var) <- c("p", "aic", "bic")
for (p in 1:20)
{
  VAR_est[[p]] <- VAR(new_x, p)
  ic_var[p,] <- c(p, AIC(VAR_est[[p]]),
                  BIC(VAR_est[[p]]))
}

ic_aic_var <- ic_var[order(ic_var[,2]),]
ic_bic_var <- ic_var[order(ic_var[,3]),]

ic_aic_var
ic_bic_var

# Adequate set from AIC, BIC
adq_set_var <- as.matrix(ic_var[2:5,])
adq_idx_var <- c(2:5)

# Check White Noise Residuals
nmods <- length(adq_idx_var)
for (i in 1:nmods)
{
  p <- adq_idx_var[i]
  print(paste0("Checking VAR(", p, ")"))
  print(serial.test(VAR_est[[p]], type = "BG"))
}

# Autocorrelations appear to be quite high for all models.
for (p in 1:20)
{
  print(paste0("Checking VAR(", p, ")"))
  print(serial.test(VAR_est[[p]], type = "BG"))
}
# only model with p = 6 is can be rejected with 99% confidence interval
adq_set_var <- as.matrix(ic_var[6,])
adq_idx_var <- c(6)

# use the model to forecast 2023 CPI inflation
hrz = 4
VAR_fcst <- list()
xlim <- c(length(new_dates) - 7 * hrz,
          length(new_dates) + hrz)
ylim <- c(-1,max(new_cpi) + 1)
for (i in 1:nmods)
{
  p <- adq_idx_var[i]
  VAR_fcst[[i]] <- predict(VAR_est[[p]],
                           n.ahead = hrz)
  plot(VAR_fcst[[i]], names = "new_cpi",
       xlim = xlim, ylim = ylim,
       main = "Forecast for CPI Inflation",
       xlab = "Horizon",
       ylab = "RRP")
}