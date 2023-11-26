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
A <- N * S
n <- 190 # Total sample size.

## ----opt_Neyman---------------------------------------------------------------
xopt <- opt(n = n, A = A)
xopt
sum(xopt) == n
# Variance of the st. estimator that corresponds to the optimum allocation.
var_st_tsi(xopt, N, S)

## ----opt_M--------------------------------------------------------------------
M <- c(100, 90, 70, 80) # Upper bounds imposed on the sample sizes in strata.
all(M <= N)
n <= sum(M)

# Solution to Problem 1.
xopt <- opt(n = n, A = A, M = M)
xopt
sum(xopt) == n
all(xopt <= M) # Does not violate upper-bounds constraints.
# Variance of the st. estimator that corresponds to the optimum allocation.
var_st_tsi(xopt, N, S)

## ----opt_m--------------------------------------------------------------------
m <- c(50, 120, 1, 2) # Lower bounds imposed on the sample sizes in strata.
n >= sum(m)

# Solution to Problem 2.
xopt <- opt(n = n, A = A, m = m)
xopt
sum(xopt) == n
all(xopt >= m) # Does not violate lower-bounds constraints.
# Variance of the st. estimator that corresponds to the optimum allocation.
var_st_tsi(xopt, N, S)

## ----opt_box------------------------------------------------------------------
m <- c(100, 90, 500, 50) # Lower bounds imposed on sample sizes in strata.
M <- c(300, 400, 800, 90) # Upper bounds imposed on sample sizes in strata.
n <- 1284
n >= sum(m) && n <= sum(M)

# Optimum allocation under box constraints.
xopt <- opt(n = n, A = A, m = m, M = M)
xopt
sum(xopt) == n
all(xopt >= m & xopt <= M) # Does not violate any lower or upper bounds constraints.
# Variance of the st. estimator that corresponds to the optimum allocation.
var_st_tsi(xopt, N, S)

## ----optcost------------------------------------------------------------------
A <- c(3000, 4000, 5000, 2000)
A0 <- 70000
unit_costs <- c(0.5, 0.6, 0.6, 0.3) # c_h, h = 1,...4.
M <- c(100, 90, 70, 80)
V <- 1e6 # Variance constraint.
V >= sum(A^2 / M) - A0

xopt <- optcost(V = V, A = A, A0 = A0, M = M, unit_costs = unit_costs)
xopt
sum(A^2 / xopt) - A0 == V
all(xopt <= M)

## ----rounding-----------------------------------------------------------------
m <- c(100, 90, 500, 50)
M <- c(300, 400, 800, 90)
n <- 1284

# Optimum, non-integer allocation under box constraints.
xopt <- opt(n = n, A = A, m = m, M = M)
xopt

xopt_int <- round_oric(xopt)
xopt_int

## ----finit_prec1--------------------------------------------------------------
N <- c(3000, 4000, 5000, 2000)
S <- c(48, 79, 76, 17)
a <- N * S
n <- 190
xopt <- opt(n = n, A = A) # which after simplification is (n / sum(a)) * a
xopt

## ----finit_prec2--------------------------------------------------------------
sum(xopt) == n

## ----finit_prec3--------------------------------------------------------------
options(digits = 22)
sum(xopt)

sum((n / sum(A)) * A) == n # mathematically, it should be TRUE!

