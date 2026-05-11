ok  <- 0L
fail <- 0L

check <- function(label, expr) {
  tryCatch({
    expr
    message("  [PASS] ", label)
    ok  <<- ok + 1L
  }, error = function(e) {
    message("  [FAIL] ", label, " — ", conditionMessage(e))
    fail <<- fail + 1L
  })
}

message("── R package tests ──────────────────────────────────")

# Load checks
check("library(tidyverse)",  library(tidyverse))
check("library(magrittr)",   library(magrittr))
check("library(splines)",    library(splines))
check("library(fixest)",     library(fixest))
check("library(glmmTMB)",    library(glmmTMB))
check("library(MuMIn)",      library(MuMIn))
check("library(INLA)",       library(INLA))

# Functional checks
check("bs() produces a matrix", {
  x <- seq(0, 1, length.out = 50)
  mat <- bs(x, df = 4)
  stopifnot(is.matrix(mat), ncol(mat) == 4)
})

check("fixest::fepois fits on toy data", {
  set.seed(42)
  df <- data.frame(
    y      = rpois(200, lambda = 5),
    x      = rnorm(200),
    group  = rep(letters[1:10], 20)
  )
  fit <- fepois(y ~ x | group, data = df)
  stopifnot(!is.null(coef(fit)))
})

check("glmmTMB fits on toy data", {
  set.seed(42)
  df <- data.frame(
    y     = rpois(100, lambda = 3),
    x     = rnorm(100),
    group = rep(1:10, 10)
  )
  fit <- glmmTMB(y ~ x + (1 | group), data = df, family = poisson)
  stopifnot(!is.null(fixef(fit)))
})

check("INLA version readable", {
  inla.version()  # prints to stdout; succeeding without error is the check
})

message("─────────────────────────────────────────────────────")
message(sprintf("Results: %d passed, %d failed", ok, fail))

if (fail > 0L) quit(status = 1L)
