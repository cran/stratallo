## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----load_package-------------------------------------------------------------
library(stratallo)

## ----pop----------------------------------------------------------------------
N <- c(3000, 4000, 5000, 2000) # Strata sizes.
S <- c(48, 79, 76, 16) # Standard deviations of a study variable in strata.
a <- N * S
n <- 190 # Total sample size.

## ----opt_Neyman---------------------------------------------------------------
opt <- opt(n = n, a = a)
opt
sum(opt) == n
# Variance of the stratified estimator that corresponds to optimum allocation.
var_st_tsi(opt, N, S)

## ----opt_M--------------------------------------------------------------------
M <- c(100, 90, 70, 80) # Upper bounds imposed on the sample sizes in strata.
all(M <= N)
n <= sum(M)

# Solution to Problem 1.
opt <- opt(n = n, a = a, M = M)
opt
sum(opt) == n
all(opt <= M) # Does not violate upper-bounds constraints.
# Variance of the stratified estimator that corresponds to optimum allocation.
var_st_tsi(opt, N, S)

## ----opt_m--------------------------------------------------------------------
m <- c(50, 120, 1, 2) # Lower bounds imposed on the sample sizes in strata.
n >= sum(m)

# Solution to Problem 2.
opt <- opt(n = n, a = a, m = m)
opt
sum(opt) == n
all(opt >= m) # Does not violate lower-bounds constraints.
# Variance of the stratified estimator that corresponds to optimum allocation.
var_st_tsi(opt, N, S)

## ----opt_box------------------------------------------------------------------
m <- c(100, 90, 500, 50) # Lower bounds imposed on sample sizes in strata.
M <- c(300, 400, 800, 90) # Upper bounds imposed on sample sizes in strata.
n <- 1284
n >= sum(m) && n <= sum(M)

# Optimum allocation under box constraints.
opt <- opt(n = n, a = a, m = m, M = M)
opt
sum(opt) == n
all(opt >= m & opt <= M) # Does not violate any lower or upper bounds constraints.
# Variance of the stratified estimator that corresponds to optimum allocation.
var_st_tsi(opt, N, S)

## ----optcost------------------------------------------------------------------
a <- c(3000, 4000, 5000, 2000)
a0 <- 70000
unit_costs <- c(0.5, 0.6, 0.6, 0.3) # c_h, h = 1,...4.
M <- c(100, 90, 70, 80)
V <- 1e6 # Variance constraint.
V >= sum(a^2 / M) - a0

opt <- optcost(V = V, a = a, a0 = a0, M = M, unit_costs = unit_costs)
opt
sum(a^2 / opt) - a0 == V
all(opt <= M)

## ----rounding-----------------------------------------------------------------
m <- c(100, 90, 500, 50)
M <- c(300, 400, 800, 90)
n <- 1284

# Optimum, non-integer allocation under box constraints.
opt <- opt(n = n, a = a, m = m, M = M)
opt

opt_int <- round_oric(opt)
opt_int

## ----finit_prec1--------------------------------------------------------------
N <- c(3000, 4000, 5000, 2000)
S <- c(48, 79, 76, 17)
a <- N * S
n <- 190
opt <- opt(n = n, a = a) # which after simplification is (n / sum(a)) * a
opt

## ----finit_prec2--------------------------------------------------------------
sum(opt) == n

## ----finit_prec3--------------------------------------------------------------
options(digits = 22)
sum(opt)

sum((n / sum(a)) * a) == n # mathematically, it should be TRUE!

