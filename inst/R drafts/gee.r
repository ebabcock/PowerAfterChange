library(geepack)

# Ensure correct factor ordering (descending = event of interest is "1")
b2$exceed_bin <- factor(b2$exceed_bin, levels = c(1, 0))  # descending order

# Convert class variables to factors
b2$site   <- factor(b2$site)
b2$season <- factor(b2$season)
b2$wyr    <- factor(b2$wyr)
b2$post   <- factor(b2$post)

# Fit GEE model with AR(1) correlation structure
gee_model <- geeglm(
  exceed_bin ~ post + wyr + season,
  data    = b2,
  family  = binomial(link = "logit"),
  id      = site,
  corstr  = "ar1"
)

# View results (equivalent to GEEEmpPEst / empirical parameter estimates)
gee_estimates <- summary(gee_model)
print(gee_estimates)

# Extract coefficient table as a data frame
gee_estimates_df <- as.data.frame(coef(summary(gee_model)))
print(gee_estimates_df)