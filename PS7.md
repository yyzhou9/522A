

# Q1 

## (a)

To show:

$$
\hat{\beta}_{FE} = \beta + \left( \sum_{i=1}^{N} \sum_{t=1}^{T} (x_{it} - \bar{x}_i)(x_{it} - \bar{x}_i)' \right)^{-1} \sum_{i=1}^{N} \sum_{t=1}^{T} (x_{it} - \bar{x}_i)(u_{it} - \bar{u}_i).
$$

Since $E[\hat{\beta}_{FE}] = \beta$, and demeaning removes $\alpha_i$ and $z_{i}'\theta$,

$$
\hat{\beta}_{FE} = \beta + \left( \sum_{i=1}^{N} \sum_{t=1}^{T} (x_{it} - \bar{x}_i)(x_{it} - \bar{x}_i)' \right)^{-1} \sum_{i=1}^{N} \sum_{t=1}^{T} (x_{it} - \bar{x}_i)u_{it}.
$$

Demeaning leads to $\bar{u}_i$ terms averaging to zero.

## (b)

To show:

$$
\sum_{i=1}^{N} \sum_{t=1}^{T} (x_{it} - \bar{x}_i)(u_{it} - \bar{u}_i) = \sum_{i=1}^{N} \sum_{t=1}^{T} (x_{it} - \bar{x}_i)u_{it}.
$$

Expand and simplify the left-hand side:

$$
\sum_{i=1}^{N} \sum_{t=1}^{T} (x_{it}u_{it} - x_{it}\bar{u}_i - \bar{x}_iu_{it} + \bar{x}_i\bar{u}_i).
$$

Terms with $\bar{u}_i$ and $\bar{x}_i$ cancel out, leaving:

$$
\sum_{i=1}^{N} \sum_{t=1}^{T} (x_{it} - \bar{x}_i)(u_{it} - \bar{u}_i) = \sum_{i=1}^{N} \sum_{t=1}^{T} (x_{it} - \bar{x}_i)u_{it}.
$$





# Q2

## (a)

```R
# Load necessary library
library(plm)

# Set seed for reproducibility
set.seed(123)

# Define the number of individuals and time periods
N <- 1000
T <- 2

# Generate the data for alpha_i, u_it, x_it, z_i
alpha_i <- rnorm(N, mean = 0, sd = 1)
u1_it <- rnorm(N * T, mean = 0, sd = 1)
u2_it <- rnorm(N * T, mean = 0, sd = 1)
x1_it <- alpha_i + u1_it
x2_it <- alpha_i + u2_it
z_i <- rnorm(N, mean = 0, sd = 1)

# Initialize y_it matrix
y_it <- matrix(NA, nrow = N, ncol = T)

# Generate y_it data based on the model
for (i in 1:N) {
  for (t in 1:T) {
    y_it[i, t] <- 1 + x1_it[i] + x2_it[i] + z_i[i] + alpha_i[i] + ifelse(t==1, u1_it[i], u2_it[i])
  }
}

# Creating a dataframe for analysis
data <- data.frame(
  id = rep(1:N, each = T),
  time = rep(1:T, times = N),
  y_it = as.vector(y_it),
  x1_it = rep(x1_it, each = T),
  x2_it = rep(x2_it, each = T),
  z_i = rep(z_i, each = T)
)

# View the first few rows of the dataframe
head(data)

```




## (b)

```R
# Numerical demonstration of near orthogonality using LLN
# For simplicity, let's demonstrate it for a single i

i <- 1 # Choosing one individual for demonstration
x_diff <- data$x1_it[data$id == i & data$time == 2] - data$x1_it[data$id == i & data$time == 1]
v_diff <- u2_it[i] - u1_it[i]

# Calculate the sample covariance as a numerical demonstration
covariance <- cov(x_diff, v_diff)
covariance

```




## (c)

```R
# OLS estimation for t = 1
ols_model_t1 <- lm(y_it ~ x1_it + x2_it + z_i, data = data[data$time == 1, ])
summary(ols_model_t1)

```




## (d)

Given the fixed effects model:
$$ y_{it} = \alpha y_{it-1} + x'_{it}\beta + z'_{i}\theta + c_i + \epsilon_{it}, $$
the fixed effect estimator $\hat{\beta}_{FE}$ is:
$$ \hat{\beta}_{FE} = \left(\sum (x_{it} - \overline{x}_{i})' (x_{it} - \overline{x}_{i})\right)^{-1} \left(\sum (x_{it} - \overline{x}_{i})' (y_{it} - \overline{y}_{i})\right). $$


## (e)

Time-invariant variables are eliminated in the fixed effect model due to time-demeaning:
$$ \tilde{y}_{it} - \tilde{y}_{it-1} = (x_{it} - x_{it-1})'\beta + (\epsilon_{it} - \epsilon_{it-1}), $$
therefore, $z_i$ and the constant term cannot be estimated.


## (f)

```R
# First stage: OLS estimation of x1_it and x2_it
first_stage <- plm(y_it ~ x1_it + x2_it, data = pdata, model = "within", index = c("id", "time"))

# Obtain residuals from the first stage
pdata$residuals <- residuals(first_stage)

# Second stage: OLS estimation for z_i using the residuals
second_stage <- lm(residuals ~ z_i, data = pdata)
summary(second_stage)

```




## (g)

First stage residuals are used in the two-stage estimation process:
$$ \tilde{\epsilon}_{it} = y_{it} - x'_{it}\hat{\beta}_{FE}. $$
Assuming instruments are uncorrelated with $\tilde{\epsilon}_{it}$, the IV estimator of $\theta$ is:
$$ \hat{\theta}_{IV} = (Z'Z)^{-1}Z'\tilde{\epsilon}. $$
Under the central limit theorem, $\hat{\theta}_{IV}$ is asymptotically normal, provided the instruments are not weakly correlated with the endogenous regressors.

## (h)

```r
# Regenerate data using the new definition of x2it
x2_it_new <- alpha_i + 0.1 * u2_it

# Update the dataframe
data$x2_it <- rep(x2_it_new, each = T)

# Re-estimate OLS for t = 1 with the new x2_it
ols_model_t1_new <- lm(y_it ~ x1_it + x2_it + z_i, data = data[data$time == 1, ])
summary(ols_model_t1_new)

```



# Q3

## (a)

```R
# Assuming the environment already has the 'pdata' dataframe from previous context.
# Estimate the coefficients by OLS of the differenced model
diff_y <- diff(pdata$y_it)  # Difference of y_it: y_i2 - y_i1
diff_x1 <- pdata$x1_it[2:2*N] - pdata$x1_it[1:2*N-1]  # Difference of x1_it: x1_i2 - x1_i1
diff_x2 <- pdata$x2_it[2:2*N] - pdata$x2_it[1:2*N-1]  # Difference of x2_it: x2_i2 - x2_i1

# Running the OLS regression on the differences
ols_diff_model <- lm(diff_y ~ diff_x1 + diff_x2)

# Output the summary of the model to check the coefficients
summary(ols_diff_model)

```




## (b)

For $T = 2$, the OLS estimator using first differences is equivalent to the fixed effects estimator because both approaches effectively remove the individual-specific effects $c_i$. The first-differenced estimator for the coefficient of $x_{it}$ is:

$$ \Delta y_{it} = \Delta x'_{it}\beta + \Delta \epsilon_{it}, $$

where $\Delta$ denotes the first difference operator. Since the fixed effects are differenced out, the estimated coefficients from the first-differenced OLS regression will numerically be the same as the fixed effects estimator.


## (c)

The fixed effect estimator is consistent under the given scenario because the assumptions necessary for consistency, such as the absence of serial correlation in $\epsilon_{it}$ and the exogeneity of $x_{it}$, hold true. By differencing the data, the fixed effects $c_i$ are eliminated, which would otherwise bias the OLS estimates.


## (d)

```R
# Assuming that 'data' is a dataframe with the original data
# and 'alpha_i', 'u_it', 'z_i' have been generated previously

# Generate x2it differently
data$x2it <- alpha_i + 2 * data$u_it + rnorm(N, mean = 0, sd = 1)

# Create the panel data frame
pdata_new_x2 <- pdata.frame(id = rep(1:N, each = T), 
                            time = rep(1:T, times = N), 
                            y_it = data$y_it, 
                            x1_it = data$x1_it, 
                            x2_it = data$x2it, 
                            z_i = data$z_i)

# Estimate the fixed effects model with the new x2it
fe_model_new_x2 <- plm(y_it ~ x1_it + x2_it, data = pdata_new_x2, model = "within", index = c("id", "time"))

# View the summary of the fixed effects model
summary(fe_model_new_x2)

```


## (e)

Even if $x_{2it}$ is generated in a different manner, as long as the new generation process does not induce correlation between $x_{2it}$ and the error term $\epsilon_{it}$, and there is no serial correlation in $\epsilon_{it}$, the fixed effect estimator remains consistent. The within transformation used in the fixed effects estimator continues to remove the individual-specific effects $c_i$.


## (f)

```R
# Assuming the initial setup is the same and 'data' is a dataframe with original data

# Generate sample for t = 1, ..., 20
T_long <- 20
data_long <- expand.grid(id = 1:N, time = 1:T_long)
data_long$alpha_i <- rep(alpha_i, times = T_long)
data_long$u_it <- rnorm(N * T_long, mean = 0, sd = 1)
data_long$x1_it <- data_long$alpha_i + data_long$u_it
data_long$x2it <- data_long$alpha_i + 0.5 * data_long$u_it

# Assuming y_it is generated the same way as in the previous example, just expanded for T_long
data_long$y_it <- 1 + data_long$x1_it + data_long$x2it + rep(z_i, times = T_long) + data_long$alpha_i + data_long$u_it

# Convert to pdata.frame for plm usage
pdata_long <- pdata.frame(data_long, index = c("id", "time"))

# Estimate the fixed effects model with the extended data
fe_model_long <- plm(y_it ~ x1_it + x2_it, data = pdata_long, model = "within")

# View the summary of the fixed effects model
summary(fe_model_long)

```


## (g)

The fixed effect estimator's consistency is preserved as long as the assumptions of no serial correlation and strict exogeneity are maintained. When data is generated across multiple time periods ($t = 1, \ldots, 20$), these conditions need to be checked to ensure they still hold. The consistency of the fixed effect estimator is contingent upon these conditions.

# Q4 

## (a)

```R
# Generate data for the dynamic panel model

set.seed(123) # For reproducibility

N <- 1000 # Number of individuals
T <- 20   # Number of time periods

# Pre-allocate space for our variables
y_it <- matrix(NA, nrow = N, ncol = T)
x1_it <- matrix(NA, nrow = N, ncol = T)
x2_it <- matrix(NA, nrow = N, ncol = T)
alpha_i <- rnorm(N, mean = 0, sd = 1)
z_i <- rnorm(N, mean = 0, sd = 1)

# Initial values
y_i0 <- rnorm(N, mean = 0, sd = 1)
u1_i1 <- rnorm(N, mean = 0, sd = 1)
u2_i1 <- rnorm(N, mean = 0, sd = 1)

# Generate data based on the given model
for (i in 1:N) {
  x1_it[i, 1] <- alpha_i[i] + u1_i1[i]
  x2_it[i, 1] <- alpha_i[i] + u2_i1[i]
  for (t in 1:T) {
    v_it <- rnorm(1, mean = 0, sd = 1)
    if (t == 1) {
      y_it[i, t] <- 1 + 0.6*y_i0[i] + x1_it[i, t] + x2_it[i, t] + z_i[i] + alpha_i[i] + v_it
    } else {
      y_it[i, t] <- 1 + 0.6*y_it[i, t-1] + x1_it[i, t-1] + x2_it[i, t-1] + z_i[i] + alpha_i[i] + v_it
      x1_it[i, t] <- alpha_i[i] + rnorm(1, mean = 0, sd = 1)
      x2_it[i, t] <- alpha_i[i] + rnorm(1, mean = 0, sd = 1)
    }
  }
}

# Convert the data to a plm-compatible format
library(plm)
pdata <- pdata.frame(data.frame(id = rep(1:N, each = T), time = rep(1:T, times = N),
                                y_it = as.vector(t(y_it)), x1_it = as.vector(t(x1_it)),
                                x2_it = as.vector(t(x2_it)), z_i = rep(z_i, each = T)),
                     index = c("id", "time"))

```




## (b)

Fixed effects estimator using data for $t = 1$ and $t = 2$ captures the within-individual variation by differencing out the individual-specific effects $\alpha_i$. When applying this to a dynamic panel model, the lagged dependent variable introduces bias known as "Nickell bias." The bias decreases as $T$ increases but is particularly pronounced when $T$ is small.


## (c)

The estimators from fixed effects are not consistent for all coefficients in the dynamic panel model when $T$ is small because of the endogeneity of the lagged dependent variable. The fixed effects estimator for $\beta$ will be biased due to the correlation between the lagged $y_{it}$ and the fixed effect $c_i$.


## (d)




## (e)

The first estimator in (b) consistently estimates the coefficient of the lagged dependent variable only when $T$ is sufficiently large to mitigate the "Nickell bias". For $T = 2$, the estimator does not consistently estimate the parameter due to the bias introduced by the within transformation of the lagged dependent variable.


## (f)

```R
# Using the panel data 'pdata' from part (a)
# Generate data with alternative method for t = 3, 4, ..., 20
set.seed(123) # Ensure reproducibility

# Generate variables using the specified dynamic panel model structure
for (i in 1:N) {
  for (t in 3:20) {
    pdata$x1_it[pdata$id == i & pdata$time == t] <- alpha_i[i] + rnorm(1, mean = 0, sd = 1)
    pdata$x2_it[pdata$id == i & pdata$time == t] <- alpha_i[i] + rnorm(1, mean = 0, sd = 1)
    pdata$y_it[pdata$id == i & pdata$time == t] <- 1 + 0.6 * pdata$y_it[pdata$id == i & pdata$time == t - 1] +
      pdata$x1_it[pdata$id == i & pdata$time == t] +
      pdata$x2_it[pdata$id == i & pdata$time == t] + z_i[i] + alpha_i[i] + rnorm(1, mean = 0, sd = 1)
  }
}

# Now pdata has the updated values for x1_it and x2_it for all t

```




## (g)

```R
# Estimate the coefficients using IV/2SLS
library(AER)  # For ivreg() function

# Define the model formula
iv_formula <- y_it ~ lag(y_it, 1) + x1_it + x2_it | lag(x1_it, 1) + lag(x2_it, 1) + lag(y_it, 2)

# Run the IV regression on pdata
iv_model <- ivreg(iv_formula, data = pdata)

# Output the summary of the model
summary(iv_model)

```




## (h)

```R
# Generate data as in (f) but with alternative instruments for t = 3, 4, ..., 20
# We use the same pdata from part (a) and continue from there

# Define the new instruments based on the problem statement
pdata$inst1 <- lag(pdata$x1_it, 1) - lag(pdata$x1_it, 2)
pdata$inst2 <- lag(pdata$x2_it, 1) - lag(pdata$x2_it, 2)

# Define the model formula for IV/2SLS using the new instruments
iv_formula_new <- y_it ~ lag(y_it, 1) + x1_it + x2_it | inst1 + inst2

# Run the IV regression on pdata
iv_model_new <- ivreg(iv_formula_new, data = pdata)

# Output the summary of the model
summary(iv_model_new)

```




## (i)

```R
# Assuming the environment has 'pdata' from part (a) with the appropriate structure

# Estimate coefficients using IV for dynamic panel data
library(AER)  # Loading the AER package for ivreg()

# Using the lagged differences of x1it and x2it as instruments
iv_formula_i <- y_it ~ lag(y_it, 1) + x1_it + x2_it |
                 lag(x1_it, 1) - lag(x1_it, 2) + lag(x2_it, 1) - lag(x2_it, 2)

# Run the IV regression for all t = 1, ..., 20
iv_model_i <- ivreg(iv_formula_i, data = pdata)

# Output the summary of the model
summary(iv_model_i)

```




## (j)

When the coefficient on the lagged dependent variable is increased to 0.9, the first IV estimator tends to have a larger standard error because the endogeneity problem becomes more pronounced. The larger coefficient implies a stronger correlation between the lagged dependent variable and the error term, increasing the potential for bias and inconsistency in the OLS estimator, and thus the IV estimator will exhibit a larger standard error to reflect this greater potential variability.


## (k)

With the change in the model where $v_{it} = 0.5v_{it-1} + e_{it}$, the estimators considered in (f) and (g) would still be consistent under the assumption that the $e_{it}$ is uncorrelated with past values of the dependent variable and other regressors. The lagged dependent variable and the predetermined regressors in the model no longer serve as valid instruments if they are correlated with $v_{it}$ due to the introduced serial correlation. This necessitates the search for external instruments that are uncorrelated with the composite error term $v_{it}$.


## (l)

```R
# Generate new data as described in part (k) of the problem
set.seed(123) # Set seed for reproducibility

N <- 1000     # Number of individuals
T <- 20       # Number of time periods

# Generating initial v_i0 and e_it as per the new model structure
v_i0 <- rnorm(N, mean = 0, sd = 1)
e_it <- rnorm(N * T, mean = 0, sd = 1)

# Generate v_it based on the new process
v_it <- matrix(NA, nrow = N, ncol = T)
v_it[, 1] <- 0.5 * v_i0 + e_it[1:N]

for (t in 2:T) {
  v_it[, t] <- 0.5 * v_it[, t - 1] + e_it[(t - 1) * N + 1:t * N]
}

# Update pdata with the new v_it
pdata$v_it <- as.vector(t(v_it))

# Recalculate y_it using the new v_it
for (i in 1:N) {
  for (t in 1:T) {
    pdata$y_it[pdata$id == i & pdata$time == t] <- 1 + 0.6 * pdata$y_it[pdata$id == i & pdata$time == t - 1] +
      pdata$x1_it[pdata$id == i & pdata$time == t] +
      pdata$x2_it[pdata$id == i & pdata$time == t] + pdata$z_i[pdata$id == i] +
      alpha_i[i] + pdata$v_it[pdata$id == i & pdata$time == t]
  }
}

# Re-estimate the fixed effects model using the new data
fe_model_l <- plm(y_it ~ lag(y_it, 1) + x1_it + x2_it, data = pdata, model = "within")

# Output the summary of the model
summary(fe_model_l)

```



# Q5 

## (a)

To show that fixed effect estimator of $\beta$ is consistent and asymptotically normal:
Given model $y_{it} = x'_{it}\beta + z'_i\theta + c_i + \epsilon_{it}$ with assumptions:

1. $E(\epsilon_{it} | x_{i1}, \ldots, x_{iT}, z_i, c_i) = 0$
2. $Var(\epsilon_{it} | x_{i1}, \ldots, x_{iT}, z_i, c_i) = \sigma^2 < \infty$
3. $Cov(\epsilon_{is}, \epsilon_{it} | x_{i1}, \ldots, x_{iT}, z_i, c_i) = 0$ for $s \neq t$

Consistency follows from $E(\hat{\beta}_{FE} - \beta) = 0$ and Law of Large Numbers.

Asymptotic normality follows from $\sqrt{n}(\hat{\beta}_{FE} - \beta) \rightarrow N(0, \sigma^2(Q'Q)^{-1})$ by Central Limit Theorem where $Q = \frac{1}{\sqrt{n}}\sum_{i=1}^{n}x_{it}\epsilon_{it}$.


## (b)

Efficient GMM estimator $\hat{\beta}_{GMM}$ is obtained by minimizing the GMM criterion:

$\hat{\beta}_{GMM} = \arg \min_{\beta} (\sum_{i=1}^{n} x_{it}\epsilon_{it})'W(\sum_{i=1}^{n} x_{it}\epsilon_{it})$

where $W$ is an optimal weight matrix derived from the variance of the errors.


## (c)

For paired observations $(\epsilon_{i1t}, \epsilon_{i2t})$ with correlation and odd $i$:

Even with correlation, fixed effect estimator remains consistent if cross-sectional independence is maintained.


## (d)

For efficient GMM estimator with paired observations:

Use a block-diagonal weight matrix $W$ to construct an estimator that accounts for the within-pair correlation.

$\hat{\beta}_{GMM} = \arg \min_{\beta} (\sum_{i=1}^{n/2} x_{it}\epsilon_{it})'W(\sum_{i=1}^{n/2} x_{it}\epsilon_{it})$

where the weight matrix $W$ compensates for the correlated errors within pairs.




# Q6 

## (a)

Consider the panel data model $y_{it} = \alpha y_{it-1} + x'_{it}\beta + z'_i\theta + c_i + \epsilon_{it}$.

Instrumental variables (IV) are $x_{it-1} - x_{it-2}$ and $x_{it} - x_{it-1}$.

Differenced model is $y_{it} - y_{it-1} = \alpha(y_{it-1} - y_{it-2}) + (x'_{it} - x'_{it-1})\beta + \epsilon_{it} - \epsilon_{it-1}$.

Assuming $E(\epsilon_{it}|x_{i1}, \ldots, x_{iT}, z_i, c_i) = 0$ and $Var(\epsilon_{it}|x_{i1}, \ldots, x_{iT}, z_i, c_i) = \sigma^2$ with $Cov(\epsilon_{is}, \epsilon_{it}|x_{i1}, \ldots, x_{iT}, z_i, c_i) = 0$ for $s \neq t$,

The estimator $\hat{\beta}$ is consistent because:

$E[(x_{it-1} - x_{it-2})(\epsilon_{it} - \epsilon_{it-1})] = E[(x_{it} - x_{it-1})(\epsilon_{it} - \epsilon_{it-1})] = 0$.

Asymptotic normality follows from the central limit theorem since $\sqrt{n}(\hat{\beta} - \beta) \xrightarrow{d} N(0, \sigma^2(Var[X])^{-1})$.


## (b)

Efficient GMM estimator $\hat{\beta}_{GMM}$ for the differenced model:

$\hat{\beta}_{GMM} = \arg \min_{\beta} (\sum_{i=1}^{n} (x_{it-1} - x_{it-2})'(\epsilon_{it} - \epsilon_{it-1}))'W(\sum_{i=1}^{n} (x_{it-1} - x_{it-2})'(\epsilon_{it} - \epsilon_{it-1}))$

$W$ is a weight matrix that is a function of instruments and error terms.


## (c)

With paired observations $(\epsilon_{i1t}, \epsilon_{i2t})$ and correlated errors for odd $i$, IV estimator:

The presence of correlation between $\epsilon_{i1t}$ and $\epsilon_{i2t}$ for odd $i$ does not affect the consistency of the IV estimator if $E[\epsilon_{i1t} \epsilon_{i2t}|X, Z] = 0$ is maintained.

The estimator remains asymptotically normal as pairs are treated as individual observations, and instruments are valid.


## (d)

Efficient GMM estimator in the context of paired observations with correlated errors:

Adjust $W$ in the GMM criterion to reflect the correlation structure of the errors.

$\hat{\beta}_{GMM} = \arg \min_{\beta} (\sum_{i=1}^{n/2} (x_{i1t-1} - x_{i1t-2})'(\epsilon_{i1t} - \epsilon_{i1t-1}))'W(\sum_{i=1}^{n/2} (x_{i1t-1} - x_{i1t-2})'(\epsilon_{i1t} - \epsilon_{i1t-1}))$

$W$ accounts for the pairing of observations and error correlation within each pair.



