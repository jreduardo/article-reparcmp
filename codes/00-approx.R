# =====================================================================
# Analysis of moments approximation
#                                                        Eduardo Junior
#                                                    edujrrib@gmail.com
#                                                            2017-08-01
# =====================================================================

#-----------------------------------------------------------------------
# Load package and functions
source("config.R")
source("functions.R")

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

#-------------------------------------------
# Mean and variance relationship
aux <- expand.grid(
    mu = seq(2, 30, length.out = 50),
    phi = seq(log(0.3), log(2.5), length.out = 50))

moments <- mapply(FUN = calc_moments,
                  mu = aux$mu,
                  phi = aux$phi,
                  MoreArgs = list(sumto = 300),
                  SIMPLIFY = FALSE)
grid <- cbind(aux, t(do.call(cbind, moments)))
grid <- transform(grid, va = mu / exp(phi))

col <- brewer.pal(n = 8, name = "RdBu")
myreg <- colorRampPalette(colors = col)
xy1 <- xyplot(var ~ mean,
              groups = phi,
              data = grid,
              type = "l",
              lwd = 2,
              axis = axis.grid,
              xlab = expression(E(X)),
              ylab = expression(V(X)),
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
              par.settings = modifyList(
                  ps2, list(superpose.line = list(
                                col = myreg(length(unique(grid$phi))))
                            )
              ),
              panel = function(x, y, ...) {
                  panel.xyplot(x, y, ...)
                  panel.curve(1*x, min(x), max(x), lty = 2)
              },
              page = function(...) {
                  grid.text(expression(log(nu)),
                            just = "bottom",
                            x = unit(0.10, "npc"),
                            y = unit(0.83, "npc"))
              })

#-------------------------------------------
# Errors in approximations for E[Y] and V[Y]
grid <- transform(grid,
                  emu = (mu - mean)^2,
                  eva = (va - var)^2)

myreg <- colorRampPalette(c("gray90",  "gray20"))(100)
xy2 <- levelplot(emu ~ phi + mu, data = grid,
                 aspect = "fill",
                 col.regions = myreg,
                 xlab = expression(phi),
                 ylab = expression(mu),
                 sub = "(b)",
                 colorkey = list(space = "top"),
                 par.settings = ps2,
                 panel = function(x, y, z, ...) {
                     panel.levelplot(x, y, z, ...)
                     panel.curve(10 - ( exp(x) - 1)/(2 * exp(x)),
                                 lty = 2)
                     panel.abline(v = 0, lty = 2)
                 })

xy3 <- levelplot(eva ~ phi + mu, data = grid,
                 aspect = "fill",
                 col.regions = myreg,
                 xlab = expression(phi),
                 ylab = expression(mu),
                 sub = "(c)",
                 colorkey = list(space = "top"),
                 par.settings = ps2,
                 panel = function(x, y, z, ...) {
                     panel.levelplot(x, y, z, ...)
                     panel.curve((10 - ( exp(x) - 1)/
                                  (2 * exp(x)))/exp(x), lty = 2)
                     panel.abline(v = 0, lty = 2)
                 })

print(xy1, split = c(1, 1, 3, 1), more = TRUE)
print(xy2, split = c(2, 1, 3, 1), more = TRUE)
print(xy3, split = c(3, 1, 3, 1), more = FALSE)

#-----------------------------------------------------------------------
# COM-Poisson probability mass function (mean parametrization)
grid <- expand.grid(mu = c(2, 8, 15), phi = c(-1, 0, 1))
y <- 0:30
py <- mapply(FUN = dcmp,
             mu = grid$mu,
             phi = grid$phi,
             MoreArgs = list(y = y, sumto = 100),
             SIMPLIFY = FALSE)
grid <- cbind(grid[rep(1:nrow(grid), each = length(y)), ],
              y = y,
              py = unlist(py))

# COM-Poisson p.m.f. to different combination betwenn phi and mu
leg_phi <- parse(
    text = paste("phi == \"",
                 formatC(unique(grid$phi), 1, format = "f"),
                 "\""))
barchart(py ~ y | factor(mu),
         groups = factor(phi),
         data = grid,
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
