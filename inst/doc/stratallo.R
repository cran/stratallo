## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----load_package, echo=FALSE, message=FALSE, warning=FALSE-------------------
library(stratallo)

## ----datasets-----------------------------------------------------------------
data(package = "stratallo")

## ----pop----------------------------------------------------------------------
N <- c(3000, 4000, 5000, 2000) # Strata sizes.
S <- c(48, 179, 176, 16) # Standard deviations of a study variable in strata.
A <- N * S
n <- 190 # Total sample size.

## ----opt_Neyman---------------------------------------------------------------
x <- opt(n = n, A = A)
x

# Variance of the st. estimator that corresponds to the optimum allocation.
var_stsi(x, N, S)

# Round non-integer allocation.

x_int <- round(x)
x_int
sum(x_int)

x_int_oric <- round_oric(x)
x_int_oric
sum(x_int_oric)
x_int_oric <= N
var_stsi(x_int_oric, N, S)

## ----opt_M--------------------------------------------------------------------
M <- c(100, 90, 70, 80) # Upper bounds.
all(M <= N)
n <= sum(M)

x <- opt(n = n, A = A, M = M)
x

# Variance of the st. estimator that corresponds to the optimum allocation.
var_stsi(x, N, S)

## ----opt_box------------------------------------------------------------------
m <- c(100, 90, 500, 50) # Lower bounds.
M <- c(300, 400, 800, 90) # Upper bounds.
n <- 1284
n >= sum(m) && n <= sum(M)

x <- opt(n = n, A = A, m = m, M = M)
x

var_stsi(x, N, S)

## ----optcost------------------------------------------------------------------
A <- c(3000, 4000, 5000, 2000)
A0 <- 70000
unit_costs <- c(0.5, 0.6, 0.6, 0.3) # c_h, h = 1,...,4.
M <- c(100, 90, 70, 80)
V <- 1e6 # Variance constraint.
V >= sum(A^2 / M) - A0

optcost(V = V, A = A, A0 = A0, M = M, unit_costs = unit_costs)

## ----dopt---------------------------------------------------------------------
# Three domains with 2, 2, and 3 strata, respectively.
H_counts <- c(2, 2, 3)
N <- c(140, 110, 135, 190, 200, 40, 70)
S <- c(180, 20, 5, 4, 35, 9, 40)
total <- c(2, 3, 5)
kappa <- c(0.5, 0.2, 0.3)
n <- 828

dopt(n, H_counts, N, S, total, kappa)

## ----fprec_pop----------------------------------------------------------------
N <- 101:104 # strata sizes
S <- 1001:1004 # standard deviations in strata
A <- N * S
n <- 409L # total sample size

## ----fprec_x------------------------------------------------------------------
x <- opt(n = n, A = A)
x

## ----fprec_sumx---------------------------------------------------------------
sum(x) == n
sum((n / sum(A)) * A) == n

## ----finit_sumx22-------------------------------------------------------------
options(digits = 22)
sum(x)

