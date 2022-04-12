## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----load_package-------------------------------------------------------------
library(stratallo)

## ----example_dopt_pop---------------------------------------------------------
# Define example population
N <- c(3000, 4000, 5000, 2000) # Strata sizes.
S <- c(48, 79, 76, 17) # Standard deviations of a study variable in strata.
a <- N * S

## ----example_dopt_M-----------------------------------------------------------
M <- c(100, 90, 70, 80) # Upper bounds constraints imposed on the sample sizes in strata.
all(M <= N)
n <- 190 # Total sample size.
n < sum(M)

# Solution to Problem 1.
opt <- dopt(n = n, a = a, M = M)
opt
sum(opt) # Equals to n.
all(opt <= M) # Does not violate upper bounds constraints.
# Variance of the pi-estimator that corresponds to a given optimal allocation.
var_tst_si(opt, N, S)

## ----example_dopt_m-----------------------------------------------------------

m <- c(50, 120, 1, 1) # Lower bounds constraints imposed on the sample sizes in strata.
n > sum(m)

# Solution to Problem 2.
opt <- dopt(n = n, a = a, m = m)
opt
sum(opt) # Equals to n.
all(opt >= m) # Does not violate lower bounds constraints.
# Variance of the pi-estimator that corresponds to a given optimal allocation.
var_tst_si(opt, N, S)

## ----example_dopt_Neyman------------------------------------------------------

# Tschuprov-Neyman allocation (no inequality constraints).
opt <- dopt(n = n, a = a)
opt
sum(opt) # Equals to n.
# Variance of the pi-estimator that corresponds to a given optimal allocation.
var_tst_si(opt, N, S)

## ----example_nopt-------------------------------------------------------------

a <- c(3000, 4000, 5000, 2000)
b <- 70000
M <- c(100, 90, 70, 80)
D <- 1e6 # Variance constraint.

opt <- nopt(D, a, b, M)
sum(opt)

