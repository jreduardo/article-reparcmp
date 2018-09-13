# =====================================================================
# Analysis of moments approximation
#                                                        Eduardo Junior
#                                                    edujrrib@gmail.com
#                                                            2017-08-01
# =====================================================================

#-----------------------------------------------------------------------
# Load package and functions
source("helper01_general-functions.R")
source("helper02_lattice-panels.R")

library(dplyr)

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
    mu = seq(3, 30, length.out = 55),
    phi = seq(-1.5, 1.8, length.out = 50),
    KEEP.OUT.ATTRS = FALSE)

grid <- with(aux, t(calc_moments(mu, phi))) %>%
    as.data.frame() %>%
    cbind(aux, .) %>%
    transform(va = mu / exp(phi)) %>%
    transform(emu = (mu - mean),
              eva = (va - var))

max(abs(grid$emu))
max(abs(grid$eva))

#-------------------------------------------
# Errors in approximations for E[Y] and V[Y]

myreg <- colorRampPalette(c("gray20", "gray95"))
xy1 <-
    levelplot(emu ~ phi + mu, ,
              aspect = "fill",
              xlab = expression(phi),
              ylab = expression(mu),
              sub = "(a)",
              at = seq(-0.2, 0.07, length.out = 18),
              par.settings = ps2,
              colorkey = list(space = "top"),
              col.regions = myreg,
              panel = function(x, y, z, ...) {
                  panel.levelplot(x, y, z, ...)
                  panel.curve(10 - (exp(x) - 1)/(2 * exp(x)),
                              lty = 2)
                  panel.abline(v = 0, lty = 2)
              },
              data = grid)

xy2 <-
    levelplot(eva ~ phi + mu, ,
              aspect = "fill",
              xlab = expression(phi),
              ylab = expression(mu),
              sub = "(b)",
              at = seq(-0.15, 10, length.out = 18),
              par.settings = ps2,
              colorkey = list(space = "top"),
              # col.regions = myreg,
              panel = function(x, y, z, ...) {
                  panel.levelplot(x, y, z, ...)
                  panel.curve(10 - (exp(x) - 1)/(2 * exp(x)),
                              lty = 2)
                  panel.abline(v = 0, lty = 2)
              },
              data = grid)

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
grid <- transform(grid,
                  dispersion_index = var / mean)
xy2 <- xyplot(dispersion_index ~ mu | "Dispersion index",
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
choosemu <- 25
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
