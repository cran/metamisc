### Some stuff necessary for testing
set.seed(8092017)
n <- 100
n.cov <- 3
td <- data.frame(matrix(rbinom(n * (n.cov + 1), 1, .5), ncol = n.cov + 1, nrow = n))
td.ig <- td + 1 # For inverse gaussian and Gamma.

### To be included:
# one-stage
# predFUN.

test_that("The predict functions predict accurately.", {
  m.bi <- glm(td, family = binomial)
  expect_true(is.function(pm <- metamisc:::getPredictMethod(m.bi, two.stage = TRUE)))
  expect_true(all(unlist(pm(m.bi, td, coef(m.bi))) == unlist(m.bi$fitted.values))) # == intentionally ignores names.

  m.lm <- lm(td)
  expect_true(is.function(pm <- metamisc:::getPredictMethod(m.lm, two.stage = TRUE)))
  expect_true(all.equal(unlist(pm(m.lm, td, coef(m.lm))) ,as.matrix(unlist(m.lm$fitted.values)),
                        use.names = F, check.attributes = F)) # also intentionally ignores names.

  m.no <- glm(td)
  expect_true(is.function(pm <- metamisc:::getPredictMethod(m.no, two.stage = TRUE)))
  expect_true(all(unlist(pm(m.no, td, coef(m.no))) == unlist(m.no$fitted.values))) # == intentionally ignores names.

  m.gm <- glm(td.ig, family = Gamma)
  expect_true(is.function(pm <- metamisc:::getPredictMethod(m.gm, two.stage = TRUE)))
  expect_true(all.equal(unlist(pm(m.gm, td.ig, coef(m.gm))) ,as.matrix(unlist(m.gm$fitted.values)),
                        use.names = F, check.attributes = F)) # also intentionally ignores names.

  m.ig <- glm(td.ig, family = inverse.gaussian)
  expect_true(is.function(pm <- metamisc:::getPredictMethod(m.ig, two.stage = TRUE)))
  expect_true(all(unlist(pm(m.ig, td.ig, coef(m.ig))) == unlist(m.ig$fitted.values))) # == intentionally ignores names.

  m.po <- glm(td, family = poisson)
  expect_true(is.function(pm <- metamisc:::getPredictMethod(m.po, two.stage = TRUE)))
  expect_true(all(unlist(pm(m.po, td, coef(m.po))) == unlist(m.po$fitted.values))) # == intentionally ignores names.

  m.q <- glm(td, family = quasi)
  expect_true(is.function(pm <- metamisc:::getPredictMethod(m.q, two.stage = TRUE)))
  expect_true(all.equal(unlist(pm(m.q, td, coef(m.q))) ,as.matrix(unlist(m.q$fitted.values)),
                        use.names = F, check.attributes = F)) # also intentionally ignores names.

  m.qb <- glm(td, family = quasibinomial)
  expect_true(is.function(pm <- metamisc:::getPredictMethod(m.qb, two.stage = TRUE)))
  expect_true(all(unlist(pm(m.qb, td, coef(m.qb))) == unlist(m.qb$fitted.values))) # == intentionally ignores names.

  m.qp <- glm(td, family = quasipoisson)
  expect_true(is.function(pm <- metamisc:::getPredictMethod(m.qp, two.stage = TRUE)))
  expect_true(all.equal(unlist(pm(m.qp, td, coef(m.qp))) ,as.matrix(unlist(m.qp$fitted.values)),
                        use.names = F, check.attributes = F)) # also intentionally ignores names.
})

test_that("metapred produces a model.", {
  expect_true(is.list(mp <- metapred(data = td, strata = "X4", family = binomial))) # binomial
  expect_true(inherits(mp, "metapred"))
  expect_true(is.list(mp$stepwise))
  expect_true(is.list(mp$FUN))
  expect_true(is.call(mp$call))

  expect_true(is.list(mp <- metapred(data = td, strata = "X4", family = binomial(link = "log")))) # binomial, loglink
  expect_true(inherits(mp, "metapred"))
  expect_true(is.list(mp$stepwise))
  expect_true(is.list(mp$FUN))
  expect_true(is.call(mp$call))

  expect_true(is.list(mp <- metapred(data = td, strata = "X4"))) # gaussian
  expect_true(inherits(mp, "metapred"))
  expect_true(is.list(mp$stepwise))
  expect_true(is.list(mp$FUN))
  expect_true(is.call(mp$call))

  expect_true(is.list(mp <- metapred(data = td.ig, strata = "X4", family = Gamma))) # Gamma
  expect_true(inherits(mp, "metapred"))
  expect_true(is.list(mp$stepwise))
  expect_true(is.list(mp$FUN))
  expect_true(is.call(mp$call))

  expect_true(is.list(mp <- metapred(data = td.ig, strata = "X4", family = inverse.gaussian))) # inverse.gaussian
  expect_true(inherits(mp, "metapred"))
  expect_true(is.list(mp$stepwise))
  expect_true(is.list(mp$FUN))
  expect_true(is.call(mp$call))

  expect_true(is.list(mp <- metapred(data = td, strata = "X4", family = poisson))) # poisson
  expect_true(inherits(mp, "metapred"))
  expect_true(is.list(mp$stepwise))
  expect_true(is.list(mp$FUN))
  expect_true(is.call(mp$call))

  expect_true(is.list(mp <- metapred(data = td, strata = "X4", family = quasi))) # quasi
  expect_true(inherits(mp, "metapred"))
  expect_true(is.list(mp$stepwise))
  expect_true(is.list(mp$FUN))
  expect_true(is.call(mp$call))

  expect_true(is.list(mp <- metapred(data = td, strata = "X4", family = quasibinomial))) # quasibinomial
  expect_true(inherits(mp, "metapred"))
  expect_true(is.list(mp$stepwise))
  expect_true(is.list(mp$FUN))
  expect_true(is.call(mp$call))

  expect_true(is.list(mp <- metapred(data = td, strata = "X4", family = quasipoisson))) # quasipoisson
  expect_true(inherits(mp, "metapred"))
  expect_true(is.list(mp$stepwise))
  expect_true(is.list(mp$FUN))
  expect_true(is.call(mp$call))
})

test_that("metapred's stepwise is WAD.", {
  expect_true(is.list(mp <- metamisc:::metapred(data = td, strata = "X4" ))) 
  expect_true(inherits(mp, "metapred"))
  expect_true(is.list(mp$stepwise))
  expect_true(is.list(mp$FUN))
  expect_true(is.call(mp$call))
  expect_length(mp$coefficients, 2) # One is selected due to random fluctuation.
  
  set.seed(324234)
  td.none <- data.frame(matrix(rbinom(n * (n.cov + 1), 1, .5), ncol = n.cov + 1, nrow = n))
  expect_true(is.list(mp <- metamisc:::metapred(data = td.none, strata = "X4" ))) 
  expect_true(inherits(mp, "metapred"))
  expect_true(is.list(mp$stepwise))
  expect_true(is.list(mp$FUN))
  expect_true(is.call(mp$call))
  expect_length(mp$coefficients, 1) # None is selected, because the data is pure noise.
  
  td.all <- data.frame(matrix(rbinom(n * (n.cov + 1), 1, .5), ncol = n.cov + 1, nrow = n))
  td.all[ , 1] <- rowSums(td.all)
  expect_true(is.list(mp <- metamisc:::metapred(data = td.all, strata = "X4" ))) 
  expect_true(inherits(mp, "metapred"))
  expect_true(is.list(mp$stepwise))
  expect_true(is.list(mp$FUN))
  expect_true(is.call(mp$call))
  expect_length(mp$coefficients, 3) # All are selected, as predictors are good predictors.
  
  expect_true(is.list(mp <- metamisc:::metapred(data = td.none, strata = "X4", stepwise = FALSE ))) 
  expect_true(inherits(mp, "metapred"))
  expect_true(is.list(mp$stepwise))
  expect_true(is.list(mp$FUN))
  expect_true(is.call(mp$call))
  expect_length(mp$coefficients, 3)  # All noise predictors are selected, because stepwise = F.
})

test_that("coef.metapred gets the coefficients", {
  mp <- metapred(data = td, strata = "X4", family = binomial)
  expect_true(is.numeric(coef(mp)))

  gl.b <- glm(formula = X1 ~ X2 + X3, data = td, family  = binomial)
  expect_equal(length(coef(mp)), n.cov - 1)
})

# This one can be a little annoying
# test_that("print.metapred prints a metapred object", {
#   mp <- metapred(data = td, strata = "X4", family = binomial)
#   cat("\n")
#   print(mp)
# })

test_that("metapred.predict predicts.", {
  mp <- metapred(data = td, strata = "X4", family = binomial)
  p <- predict(mp, newdata = td)
  expect_true(is.numeric(p))
  expect_true(all(p <= 1))
  expect_true(all(p >= 0))
})

test_that("metapred.family and $family get the family.", {
  expect_true(is.list(mp <- metapred(data = td, strata = "X4", family = binomial(link = "log"))))
  gl <- glm(formula = X1 ~ X2 + X3, data = td, family = binomial(link = "log"))
  expect_equal(gl$family, mp$family)
  expect_equal(family(gl), family(mp))
})

test_that("metapred.formula and $formula get the formula (test-dependent).", {# formula of glm is specific for this data set!
  expect_true(is.list(mp <- metapred(data = td, strata = "X4", family = binomial(link = "log"))))
  gl <- glm(formula = X1 ~ X3, data = td, family = binomial(link = "log"))
  expect_equal(formula(gl), formula(mp))
  expect_equal(gl$formula, mp$formula) # note that glm(...)$formula is unpredictable.
})

