# =====================================================================
# Analysis of moments approximation
#                                                        Eduardo Junior
#                                                    edujrrib@gmail.com
#                                                            2017-08-01
# =====================================================================

#-----------------------------------------------------------------------
# Load package and functions
source("functions.R")
source("lattice-panels.R")

#-----------------------------------------------------------------------
# Convergence of Z(lambda, nu) constant
computeZ <- function(lambda, nu, maxit = 1e4, tol = 1e-5) {
    z <- vector("numeric", maxit)
    j = 1
    z[j] <- exp(j * log(lambda) - nu * lfactorial(j))
    while (abs(z[j] - 0) > tol && j <= maxit) {
        j = j + 1
        z[j] <- exp(j * log(lambda) - nu * lfactorial(j))
    }
    if (lambda == 1 & nu == 0) z <- Inf
    return(sum(z, na.rm = TRUE) + 1)
}

grid <- expand.grid(
    lambda = c(0.5, 1, 5, 10, 30, 50),
    nu = seq(0, 1, length.out = 11))
grid$z <- apply(grid, 1, function(par) {
    computeZ(lambda = par[1], nu = par[2], maxit = 1e4)
})

xt <- xtabs(z ~ nu + lambda, data = grid)
xt

#-----------------------------------------------------------------------
# Study the approximation

# Parameters
aux <- expand.grid(
    mu = seq(2, 30, length.out = 50),
    phi = seq(log(0.3), log(2.5), length.out = 50))

# Compute approximately and numerically moments
moments <- mapply(FUN = calc_moments,
                  mu = aux$mu,
                  phi = aux$phi,
                  MoreArgs = list(sumto = 300),
                  SIMPLIFY = FALSE)
grid <- cbind(aux, t(do.call(cbind, moments)))
grid <- transform(grid, va = mu / exp(phi))

#-------------------------------------------
# Errors in approximations for E[Y] and V[Y]
grid <- transform(grid,
                  emu = (mu - mean)^2,
                  eva = (va - var)^2)

myreg <- colorRampPalette(c("gray90",  "gray20"))
xy1 <- levelplot(emu ~ phi + mu, data = grid,
                 aspect = "fill",
                 col.regions = myreg,
                 xlab = expression(phi),
                 ylab = expression(mu),
                 sub = "(a)",
                 colorkey = list(space = "top"),
                 par.settings = ps2,
                 panel = function(x, y, z, ...) {
                     panel.levelplot(x, y, z, ...)
                     panel.curve(10 - ( exp(x) - 1)/(2 * exp(x)),
                                 lty = 2)
                     panel.abline(v = 0, lty = 2)
                 })

xy2 <- levelplot(eva ~ phi + mu, data = grid,
                 aspect = "fill",
                 col.regions = myreg,
                 xlab = expression(phi),
                 ylab = expression(mu),
                 sub = "(b)",
                 colorkey = list(space = "top"),
                 par.settings = ps2,
                 panel = function(x, y, z, ...) {
                     panel.levelplot(x, y, z, ...)
                     panel.curve((10 - ( exp(x) - 1)/
                                  (2 * exp(x)))/exp(x), lty = 2)
                     panel.abline(v = 0, lty = 2)
                 })

print(xy1, split = c(1, 1, 2, 1), more = TRUE)
print(xy2, split = c(2, 1, 2, 1), more = FALSE)

#-----------------------------------------------------------------------
# COM-Poisson probability mass function (mean parametrization)
gridpar <- expand.grid(mu = c(2, 8, 15), phi = c(-1, 0, 1))
y <- 0:30
py <- mapply(FUN = dcmp,
             mu = gridpar$mu,
             phi = gridpar$phi,
             MoreArgs = list(y = y, sumto = 100),
             SIMPLIFY = FALSE)
gridpar <- cbind(gridpar[rep(1:nrow(gridpar), each = length(y)), ],
              y = y,
              py = unlist(py))

# COM-Poisson p.m.f. to different combination betwenn phi and mu
leg_phi <- parse(
    text = paste("phi == \"",
                 formatC(unique(gridpar$phi), 1, format = "f"),
                 "\""))
barchart(py ~ y | factor(mu),
         groups = factor(phi),
         data = gridpar,
         horizontal = FALSE,
         layout = c(NA, 1),
         as.table = TRUE,
         axis = axis.grid,
         origin = 0,
         xlim = extendrange(y, f = 0.01),
         border = "transparent",
         scales = list(x = list(at = pretty(y))),
         ylab = expression(P(Y==y)),
         xlab = expression(y),
         auto.key = list(
             columns = 3,
             rectangles = FALSE,
             lines = TRUE,
             text = leg_phi
         ),
         strip = strip.custom(
             strip.names = TRUE,
             var.name = expression(mu == ""),
             sep = ""))

#-----------------------------------------------------------------------
# Study indexex on the reparametrized COM-Poisson distribution
psaux <- modifyList(
    ps0, list(superpose.line = list(
                  col = myreg(length(unique(grid$phi)))),
              layout.widths = list(
                  left.padding = -0.5)
              )
)

#-------------------------------------------
# Mean and variance relationship
xy1 <- xyplot(var ~ mean | "Mean-Variance",
              groups = phi,
              data = grid,
              type = "l",
              lwd = 2,
              axis = axis.grid,
              xlab = expression(E(Y)),
              ylab = expression(var(Y)),
              sub = "(a)",
              legend = list(
                  top = list(
                      fun = draw.colorkey,
                      args = list(
                          key = list(
                              space = "top",
                              col = myreg(length(unique(grid$phi))),
                              at = unique(grid$phi),
                              draw = FALSE)))),
              par.settings = psaux,
              panel = function(x, y, ...) {
                  panel.xyplot(x, y, ...)
                  panel.curve(1*x, min(x), max(x), lty = 2)
              },
              page = function(...) {
                  grid.text(expression(phi),
                            just = "bottom",
                            x = unit(0.15, "npc"),
                            y = unit(0.83, "npc"))
              })

#-------------------------------------------
# Dispersion index (DI)
xy2 <- xyplot(var / mean ~ mu | "Dispersion index",
              groups = phi,
              data = grid,
              type = "l",
              lwd = 2,
              axis = axis.grid,
              xlab = expression(mu),
              ylab = "",
              sub = "(b)",
              legend = list(
                  top = list(
                      fun = draw.colorkey,
                      args = list(
                          key = list(
                              space = "top",
                              col = myreg(length(unique(grid$phi))),
                              at = unique(grid$phi),
                              draw = FALSE)))),
              par.settings = psaux,
              panel = function(x, y, ...) {
                  panel.xyplot(x, y, ...)
                  panel.curve(1 + 0*x, min(x), max(x), lty = 2)
              },
              page = function(...) {
                  grid.text(expression(phi),
                            just = "bottom",
                            x = unit(0.15, "npc"),
                            y = unit(0.83, "npc"))
              })

#-------------------------------------------
# Zero-inflation index
prob0 <- unlist(mapply(FUN = dcmp,
                       mu = grid$mu,
                       phi = grid$phi,
                       MoreArgs = list(y = 0, sumto = 300),
                       SIMPLIFY = FALSE))
grid <- transform(grid, zero_index = 1 + log(prob0) / mean)

xy3 <- xyplot(zero_index ~ mu | "Zero-inflation index",
              groups = phi,
              data = grid,
              type = "l",
              lwd = 2,
              axis = axis.grid,
              xlab = expression(mu),
              ylab = "",
              sub = "(c)",
              legend = list(
                  top = list(
                      fun = draw.colorkey,
                      args = list(
                          key = list(
                              space = "top",
                              col = myreg(length(unique(grid$phi))),
                              at = unique(grid$phi),
                              draw = FALSE)))),
              par.settings = psaux,
              panel = function(x, y, ...) {
                  panel.xyplot(x, y, ...)
                  panel.curve(0*x, min(x), max(x), lty = 2)
              },
              page = function(...) {
                  grid.text(expression(phi),
                            just = "bottom",
                            x = unit(0.15, "npc"),
                            y = unit(0.83, "npc"))
              })

#-------------------------------------------
# Heavy tail-index
hi_fun <- Vectorize(FUN = function(x, mu, phi, sumto) {
    probs <- dcmp(y = c(x, x + 1), mu = mu, phi = phi, sumto = sumto)
    probs[2] / probs[1]
}, vectorize.args = "x")

x <- 40:150
his <- mapply(FUN = hi_fun,
              mu = grid$mu,
              phi = grid$phi,
              MoreArgs = list(x = x, sumto = 300),
              SIMPLIFY = FALSE)

aux1 <- lapply(his, function(k) data.frame("x" = x, "heavy_index" = k))
aux1 <- do.call(rbind, aux1)
aux2 <- grid[c("mu", "phi")][rep(1:nrow(grid), each = length(x)), ]
heavy_data <- cbind(aux1, aux2)

# Choose one specific mu of unique(grid$mu)
choosemu <- 18
xy4 <- xyplot(heavy_index ~ x | "Heavy-tail index",
              groups = phi,
              # data = heavy_data,
              data = subset(heavy_data, mu == choosemu),
              type = "l",
              lwd = 2,
              axis = axis.grid,
              xlab = expression(x),
              ylab = "",
              sub = "(d)",
              legend = list(
                  top = list(
                      fun = draw.colorkey,
                      args = list(
                          key = list(
                              space = "top",
                              col = myreg(length(unique(grid$phi))),
                              at = unique(grid$phi),
                              draw = FALSE)))),
              par.settings = psaux,
              panel = function(x, y, ...) {
                  panel.xyplot(x, y, ...)
                  panel.curve(choosemu / (x + 1), from = min(x),
                              to = max(x), lty = 2)
              },
              page = function(...) {
                  grid.text(expression(phi),
                            just = "bottom",
                            x = unit(0.15, "npc"),
                            y = unit(0.83, "npc"))
              })

print(xy1, split = c(1, 1, 4, 1), more = TRUE)
print(xy2, split = c(2, 1, 4, 1), more = TRUE)
print(xy3, split = c(3, 1, 4, 1), more = TRUE)
print(xy4, split = c(4, 1, 4, 1), more = FALSE)

#-----------------------------------------------------------------------
# Likelihood surfaces

# Random generation
rcmp <- function(n, mu, phi, upper = 500L, sumto = 300L) {
    x <- 0:upper
    fx <- dcmp(x, mu, phi, sumto = sumto)
    sample(x, size = n, replace = TRUE, prob = fx)
}


# Moments function
moments_compois <- Vectorize(
    function(mu, phi, upper = 500L, sumto = 300L) {
        x <- 0:upper
        fx <- dcmp(x, mu, phi, sumto = sumto)
        me <- sum(x * fx)
        va <- sum(x^2 * fx) - me^2
        return(c("Expectation" = me, "Variance" = va))
    }, vectorize.args = c("mu", "phi"))

# System equation to find parameters
library(rootSolve)
system_equation <- function(param, Expec, DI, trace = FALSE) {
    par1 <- param[1]
    par2 <- param[2]
    if (trace) print(param)
    moments <- moments_compois(par1, par2)
    eq1 <- moments[1] -  Expec
    eq2 <- moments[2] / moments[1] - DI
    return(c(eq1, eq2))
}

# Dispersion indexes
dis <- c(0.5, 1, 3)
names(dis) <- dis

# Find parameters for fixed expectation and different DI's
pars <- lapply(dis, function(di) {
    multiroot(system_equation, start = c(5, 0),
              Expec = 5, DI = di)$root
})

# Simulate 1000 observations of the DW distribution
set.seed(97380)
das <- lapply(pars, function(obj) {
    rcmp(1000L, mu = obj[1], phi = obj[2])
})

# Verify simulations
do.call("cbind", lapply(das, function(x) {
    c("Mean" = mean(x), "Var" = var(x), "DI" = var(x)/mean(x))
}))

# Fit models
index <- seq_along(dis)
names(index) <- dis

# Fit on original parametrization
models_orpar <- lapply(index, function(i) {
    lambda <- mu2lambda(pars[[i]][1], pars[[i]][2])[[1]]
    start <- c("phi" = pars[[i]][2],
               "(Intercept)" = log(lambda))
    y <- das[[i]]
    fitcm(y ~ 1, sumto = 100L, model = "CP", start = start)
})

# Fit on propose reparametrization
models_repar <- lapply(index, function(i) {
    start <- c("phi2" = pars[[i]][2],
               "(Intercept)" = log(pars[[i]][1]))
    y <- das[[i]]
    fitcm(y ~ 1, sumto = 100L, model = "CP2", start = start)
})

# Same maximum loglikelihood
cbind("Original" = vapply(models_orpar, logLik, numeric(1)),
      "Propose"  = vapply(models_repar, logLik, numeric(1)))

#-------------------------------------------
# Deviance contours

# Original parametrization
devs_orpar <- purrr::map_df(index, function(i) {
    m0 <- models_orpar[[i]]
    co <- unname(coef(m0))
    # Wald intervals
    interval <- confint(m0, level = 0.999, method = "quad")
    aux <- lapply(as.data.frame(t(interval)),
                  function(x) seq(x[1], x[2], length.out = 50))
    grid <- do.call("expand.grid", list(aux, KEEP.OUT.ATTRS = FALSE))
    # Compute loglikelihood and deviance for grid
    grid$loglik <- apply(grid, 1, function(par) {
        -llcmp(params = par, y = m0@data$y, X = m0@data$X, sumto = 100L)
    })
    grid$deviance <- -2 * (grid$loglik - logLik(m0))
    grid$lambda <- exp(grid[, 2])
    transform(grid,
              "lambda" = exp(grid[, 2]),
              "reallambda" = mu2lambda(pars[[i]][1], pars[[i]][2])[[1]],
              "realphi" = pars[[i]][2],
              "fitlambda" = exp(co[2]),
              "fitphi" = co[1])
}, .id = "di")

# Propose parametrization
devs_repar <- purrr::map_df(index, function(i) {
    m0 <- models_repar[[i]]
    co <- unname(coef(m0))
    # Wald intervals
    interval <- confint(m0, level = 0.999, method = "quad")
    aux <- lapply(as.data.frame(t(interval)),
                  function(x) seq(x[1], x[2], length.out = 50))
    grid <- do.call("expand.grid", list(aux, KEEP.OUT.ATTRS = FALSE))
    # Compute loglikelihood and deviance for grid
    grid$loglik <- apply(grid, 1, function(par) {
        -llcmp2(params = par, y = m0@data$y, X = m0@data$X, sumto = 100L)
    })
    grid$deviance <- -2 * (grid$loglik - logLik(m0))
    transform(grid,
              "mu" = exp(grid[, 2]),
              "realmu" = pars[[i]][1],
              "realphi" = pars[[i]][2],
              "fitmu" = exp(co[2]),
              "fitphi" = co[1])
}, .id = "di")

#-------------------------------------------
# Plots deviance surfaces
niveis <- c(0.9, 0.95, 0.99)
cortes <- qchisq(niveis, df = 2)

# Original parametrization
xy1 <- useOuterStrips(
    levelplot(
        deviance ~ lambda + phi | di + "Original parametrization",
        data = devs_orpar,
        scales = list(y = list(rot = 90, relation = "free"),
                      x = "free"),
        cuts = 30,
        xlab = expression(lambda),
        ylab = expression(phi),
        colorkey = list(space = "right"),
        panel = function(x, y, z, at, region, ...,
                         subscripts = subscripts) {
            ##
            flambda <- devs_orpar$fitlambda[subscripts]
            fphi <- devs_orpar$fitphi[subscripts]
            ##
            tlambda <- devs_orpar$reallambda[subscripts]
            tphi <- devs_orpar$realphi[subscripts]
            ##
            panel.levelplot(x, y, z, at = at, region = TRUE,
                            ..., subscripts = subscripts)
            panel.contourplot(x, y, z, ..., at = cortes,
                              contour = TRUE, region = FALSE,
                              subscripts = subscripts)
            panel.abline(v = flambda, h = fphi, lty = 2)
            panel.points(x = tlambda, y = tphi,
                         lty = 2, pch = 19)
        }),
    strip = strip.custom(
        strip.names = TRUE,
        sep = " = ",
        var.name = "DI"
    )
)

# Propose parametrization
xy2 <- useOuterStrips(
    levelplot(
        deviance ~ mu + phi2 | di + "Proposed parametrization",
        data = devs_repar,
        scales = list(y = list(rot = 90, relation = "free"),
                      x = "free"),
        cuts = 30,
        xlab = expression(mu),
        ylab = expression(phi),
        colorkey = list(space = "right"),
        panel = function(x, y, z, at, region, ...,
                         subscripts = subscripts) {
            ##
            fmu <- devs_repar$fitmu[subscripts]
            fphi <- devs_repar$fitphi[subscripts]
            ##
            tmu <- devs_repar$realmu[subscripts]
            tphi <- devs_repar$realphi[subscripts]
            ##
            panel.levelplot(x, y, z, at = at, region = TRUE,
                            ..., subscripts = subscripts)
            panel.contourplot(x, y, z, ..., at = cortes,
                              contour = TRUE, region = FALSE,
                              subscripts = subscripts)
            panel.abline(v = fmu, h = fphi, lty = 2)
            panel.points(x = tmu, y = tphi,
                         lty = 2, pch = 19)
        }),
    strip = strip.custom(
        strip.names = TRUE,
        sep = " = ",
        var.name = "DI"
    )
)

# Organize plots
print(xy1, position = c(0.0, 0.49, 1.0, 1.0), more = TRUE)
print(xy2, position = c(0.0, 0.0, 0.995, 0.51), more = FALSE)
