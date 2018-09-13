#' % Reparametrization of COM-Poisson Regression Models with Applications in the Analysis of Experimental Count Data
#' % [Eduardo Elias Ribeiro Junior](http://www.leg.ufpr.br/~eduardojr) <br> [Walmes Marques Zeviani](http://www.leg.ufpr.br/~walmes) <br> [Wagner Hugo Bonat](http://www.leg.ufpr.br/~wagner) <br> [Clarice Garcia Borges Demétrio](http://www4.esalq.usp.br/pesquisa/node/25) <br> [John Hinde](http://www.nuigalway.ie/maths/jh/)
#' % August 31, 2017

#+ include=FALSE
library(knitr)
options(width = 110, digits = 4)
opts_chunk$set(
    warning = FALSE,
    message = FALSE,
    echo = TRUE,
    cache = TRUE,
    fig.width = 9,
    fig.height = 5,
    fig.align = "center",
    dev.args = list(family = "Palatino")
    )

#' # 00_explore-model #
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

#+ fig.height=4, fig.width=9
print(xy1, split = c(1, 1, 4, 1), more = TRUE)
print(xy2, split = c(2, 1, 4, 1), more = TRUE)
print(xy3, split = c(3, 1, 4, 1), more = TRUE)
print(xy4, split = c(4, 1, 4, 1), more = FALSE)

#' # 01_simulation-study #
#=======================================================================
# Simulation study about estimators
#                                                        Eduardo Junior
#                                                    edujrrib@gmail.com
#                                                            2017-08-04
#=======================================================================

##----------------------------------------------------------------------
## Load package and functions
source("helper01_general-functions.R")
source("helper02_lattice-panels.R")

library(plyr)
library(tidyr)
library(bbmle)

#-------------------------------------------
# Random generation for the COM-Poisson distribution
rcmp <- Vectorize(FUN = function(n, mu, phi) {
    par <- mu2lambda(mu, phi)
    with(par, compoisson::rcom(n = n, lambda, nu))
}, vectorize.args = "mu")

#-----------------------------------------------------------------------
# Simulation study
B <- 1000

beta <- c("b0" = 2, "b1" = 0.5, "b21" = 0.8, "b22" = -0.8)
phis <- c(0, -1.6, -1, 1.8)
names(phis) <- sprintf("phi=%s", phis)

sizes <- c(50, 100, 300, 1000)
names(sizes) <- sprintf("n=%s", sizes)

#-------------------------------------------
# Justify the choices of the parameters simulation
x1 <- seq(0, 1, length.out = 100)
x2 <- rep(letters[1:3], length.out = 100)
X <- model.matrix(~ x1 + x2)
mu <- exp(X %*% beta)

daexplain <- ldply(lapply(phis, function(phi) {
    va <- compute_variance(mu, phi, sumto = 1000)
    data.frame(mu = mu, va = va, di = va / mu, x1 = x1, x2 = x2)
}), .id = "phi")

fl <- parse(text = gsub("=", "==", levels(daexplain$phi)))
xy1 <- xyplot(mu ~ x1, groups = x2,
              type = c("g", "l"),
              lwd = 2,
              xlab = "Continous variable",
              ylab = "Counts average",
              auto.key = list(
                  column = 3,
                  points = FALSE,
                  lines = TRUE,
                  title = "Categorical variable",
                  cex.title = 1.1
              ))
xy2 <- xyplot(di ~ x1 | phi, groups = x2,
              type = c("g", "l"),
              lwd = 2,
              # scales = "free",
              layout = c(2, 2),
              xlab = "Continous variable",
              ylab = "Dispersion index",
              data = daexplain,
              strip = strip.custom(factor.levels = fl))

print(xy1, split = c(1, 1, 2, 1), more = TRUE)
print(xy2, split = c(2, 1, 2, 1), more = FALSE)

#-------------------------------------------
# Simulation
# CAUTION LONG TIME CONSUMING 12H, LOAD THE RESULTS OF COMPRESSED FILE
#+ eval=FALSE
library(parallel)
results <- mclapply(sizes, function(n) {
    x1 <- seq(0, 1, length.out = n)
    x2 <- rep(letters[1:3], length.out = n)
    X <- model.matrix(~ x1 + x2)
    mu <- exp(X %*% beta)
    lapply(phis, function(phi) {
        lapply(1:B, function(i) {
            y <- rcmp(1, mu, phi)
            ti <- system.time(
                fit <- try(fitcm(y ~ x1 + x2, start = c(phi, beta),
                                 model = "CP2", sumto = 200))
            )
            if (! class(fit) == "try-error") {
                ll <- -fit@min
                vc <- fit@vcov
                co <- fit@coef
                ci <- confint(fit, method = "quad")
            } else {
                ll <- NA
                vc <- matrix(NA, ncol = length(beta) + 1,
                             nrow = length(beta) + 1)
                co <- rep(NA, length(beta) + 1)
                ci <- matrix(NA, nrow = length(beta) + 1, ncol = 2)
            }
            colnames(ci) <- c("lwr", "upr")
            out <- list("coef" = co,
                        "cint" = ci,
                        "loglik" = ll,
                        "vcov" = vc,
                        "ctime" = ti)
            return(out)
        }) # Close replicate
    }) # Close phi
}, mc.cores = 4) # Close sizes
saveRDS(results, "simulation.rds", compress = FALSE)

#-------------------------------------------
# Load the results
results <- readRDS("simulation.rds")
results_time <- readRDS("simulation2.rds")

#-----------------------------------------------------------------------
# Compute bias
std <- lapply(results["n=50"], function(x) {
    ind <- names(x); names(ind) <- ind
    lapply(ind, function(y) {
        bhat <- t(vapply(x[[y]], "[[", double(5), "coef"))
        real <- matrix(c(phis[y], beta), byrow = TRUE,
                       nrow = B, ncol = 5)
        matrix(apply(bhat - real, 2, sd, na.rm = TRUE),
               byrow = TRUE, nrow = B, ncol = 5)
    })
})

aux <- ldply(lapply(results, function(x) {
    ind <- names(x); names(ind) <- ind
    out <- lapply(ind, function(y) {
        bhat <- t(vapply(x[[y]], "[[", double(5), "coef"))
        real <- matrix(c(phis[y], beta), byrow = TRUE,
                       nrow = B, ncol = 5)
        (bhat - real) / std[[1]][[y]]
        # (t(x[[y]]) - real) / std
    })
    ldply(out, .id = "phi")
}), .id = "n")
str(aux)

# Where is NA's?
aux$na <- apply(aux, 1, function(x) as.integer(sum(is.na(x)) != 0))
xtabs(na ~ n + phi, data = aux)

# Organize the results
aux <- na.omit(aux)
bias <- gather(aux, param, bias, phi2:b22, factor_key = TRUE)
bias$phi <- ordered(bias$phi, c("phi=-1.6", "phi=-1", "phi=0",
                                "phi=1.8"))

#-------------------------------------------
# Distributions of bias
fl <- parse(text = gsub("=", "==", levels(bias$phi)))
yl <- parse(text = c("hat(phi)",
                     paste0("hat(beta)[", c(0, 1, 21, 22), "]")))

bwplot(param ~ bias | phi, groups = n, data = bias,
       box.width = 0.15,
       pch = "|",
       grid = TRUE,
       as.table = TRUE,
       layout = c(4, 1),
       xlab = "Standardized Bias",
       fill = colorRampPalette(c("gray80",  "gray20"))(4),
       strip = strip.custom(factor.levels = fl),
       panel = panel.superpose,
       scales = list(y = list(labels = yl)),
       whisker.width = 0.4,
       auto.key = list(
           column = 4,
           rectangles = TRUE,
           points = FALSE,
           title = "Sample size",
           cex.title = 1
       ),
       par.settings = list(
           plot.symbol = list(pch = 20, cex = 0.5)
       ),
       panel.groups = function(x, y, group.number, ...) {
           my.panel.bwplot(x, y + (group.number - 2.5) / 5, ...)
           panel.abline(v = 0, lty = 2)
       })

#-------------------------------------------
# Confidence interval for bias
ci <- function(x) {
    ci <- mean(x) + c(-1, 0, 1) * qnorm(0.975) * sd(x)
    names(ci) <- c("lwr", "fit", "upr")
    ci
}
cidata <- aggregate(bias ~ n + phi + param, data = bias, ci)

key <- list(
    type = "o",
    divide = 1,
    columns = 4,
    title = "Sample size",
    cex.title = 1.1,
    lines = list(pch = c(21:24), cex = 0.8),
            text = list(names(sizes))
)

fl <- parse(text = gsub("=", "==", levels(bias$phi)))
yl <- parse(text = c("hat(phi)",
                     paste0("hat(beta)[", c(0, 1, 21, 22), "]")))

xyplot(param ~ bias | phi, groups = n, data = bias,
       panel = panel.superpose,
       layout = c(4, 1),
       key = key,
       strip = strip.custom(factor.levels = fl),
       scales = list(y = list(labels = yl)),
       ylab = "",
       xlab = "Standardized Bias",
       jitter.y = TRUE,
       factor = 0.1,
       gap = 0.3,
       alpha = 0.6,
       col = "gray80",
       cex = 0.7,
       panel.groups = function(x, y, group.number,
                               subscripts, gap, ...){
           noise <- centfac(factor(levels(bias$n)), space = gap)
           noise <- sort(noise)
           panel.xyplot(x, y + noise[group.number], ...)
       }) +
    as.layer(
        segplot(param ~ bias[, "lwr"] + bias[, "upr"] | phi,
                centers = bias[, "fit"],
                data = cidata,
                draw = FALSE,
                horizontal = TRUE,
                layout = c(4, 1),
                groups = n, gap = 0.3,
                key = key,
                pch = 21:24,
                cex = 0.7,
                lwd = 2,
                panel = function(...) {
                    panel.groups.segplot(...)
                    panel.abline(v = 0, lty = 2)
                    panel.abline(h = 1:5, col = "lightgray", lty = 2)
                })
    )

# Alternative graph (figure with all point is very large in size memory)
bwplot(param ~ bias | phi, groups = n, data = bias,
       pch = "|",
       grid = TRUE,
       as.table = TRUE,
       layout = c(4, 1),
       xlab = "Standardized Bias",
       fill = "gray80",
       col = "red",
       strip = strip.custom(factor.levels = fl),
       panel = panel.superpose,
       scales = list(y = list(labels = yl)),
       box.width = 0.05,
       whisker.width = 0.4,
       key = key,
       par.settings = list(
           plot.symbol = list(pch = 20, cex = 0.5, col = "gray70"),
           box.rectangle = list(col = "gray70"),
           box.umbrella=list(col = "gray70", lwd = 2)
       ),
       panel.groups = function(x, y, group.number, ...) {
           my.panel.bwplot(x, (y - 0.05) + (group.number - 2.5) / 5, ...)
           panel.abline(v = 0, lty = 2)
       }) +
    segplot(param ~ bias[, "lwr"] + bias[, "upr"] | phi,
            centers = bias[, "fit"],
            data = cidata,
            draw = FALSE,
            horizontal = TRUE,
            layout = c(4, 1),
            groups = n, gap = 0.3,
            key = key,
            pch = 21:24,
            cex = 0.7,
            lwd = 2,
            panel = function(...) {
                panel.groups.segplot(...)
                panel.abline(v = 0, lty = 2)
                panel.abline(h = 1:5, col = "lightgray", lty = 2)
            })

#-----------------------------------------------------------------------
# Compute coverage rate
aux <- ldply(lapply(results, function(x) {
    ind <- names(x); names(ind) <- ind
    out <- lapply(ind, function(y) {
        cint <- lapply(x[[y]], function(z) {
            cint_center <- z[["cint"]] - c(phis[y], beta)
            apply(cint_center, 1, function(w) sum(sign(w)))
        })
        do.call("rbind", cint)
    })
    ldply(out, .id = "phi")
}), .id = "n")

str(aux)

# Organize the results
aux <- na.omit(aux)
coverage <- gather(aux, param, coverage, phi2:b22, factor_key = TRUE)
coverage$phi <- ordered(coverage$phi, c("phi=-1.6", "phi=-1", "phi=0",
                                        "phi=1.8"))
coverage$n <- as.numeric(gsub("n=([0-9]+)", "\\1", coverage$n))

covdata <- aggregate(coverage ~ n + phi + param,
                     data = coverage,
                     function(x) sum(x == 0) / length(x))

fl <- parse(text = gsub("=", "==", levels(bias$phi)))
yl <- parse(text = c("hat(phi)",
                     paste0("hat(beta)[", c(0, 1, 21, 22), "]")))
xyplot(coverage ~ n | param,
       type = c("g", "p", "l"),
       groups = phi,
       data = covdata,
       pch = 19,
       lwd = 1.5,
       layout = c(NA, 1),
       xlab = "Sample size",
       ylab = "Coverage rate",
       par.settings = list(
           superpose.symbol = list(pch = 19),
           layout.heights = list(strip = 1.5)
       ),
       auto.key = list(
           column = 4,
           lines = TRUE,
           text = fl,
           title = "Dispersion levels",
           cex.title = 1.1
       ),
       strip = strip.custom(
           factor.levels = yl
       ),
       panel = function(x, y, ...) {
           panel.xyplot(x, y, ...)
           panel.abline(h = 0.95, lty = 2)
       })

#-----------------------------------------------------------------------
# Get and show the confidence intervals
aux <- ldply(lapply(results, function(x) {
    ind <- names(x); names(ind) <- ind
    out <- lapply(ind, function(y) {
        da <- lapply(x[[y]], function(z) {
            cint <- z[["cint"]]
            da <- data.frame(cint)
            da$param <- rownames(cint)
            da$real <- c(phis[y], beta)
            da$ind <- with(da, as.integer(lwr < real & upr > real))
            da
        })
        names(da) <- seq(da)
        ldply(da, .id = "i")
    })
    ldply(out, .id = "phi")
}), .id = "n")

# Organize the results
cint <- na.omit(aux)
cint$scenario <- ordered(
    with(cint, paste(n, phi, sep = "\n")),
    unique(with(cint, paste(n, phi, sep = "\n"))))

fl <- do.call(
    c, lapply(strsplit(gsub("=", "==", levels(cint$scenario)),
                       split = "\n"),
              function(x) {
                  # parse(text = sprintf("atop(%s, %s)", x[1], x[2]))
                  parse(text = sprintf("phantom()[%s]^{%s}", x[1], x[2]))
              }))
yl <- parse(text = c(paste0("beta[", c(0, 1, 21, 22), "]"), "phi"))

useOuterStrips(
    segplot(i ~ lwr + upr | param + scenario,
            data = cint,
            draw = FALSE,
            scales = list(relation = "free", draw = FALSE),
            horizontal = FALSE,
            lattice.options = list(
                layout.widths = list(strip.left = list(x = 2.5)),
                layout.heights = list(strip = list(x = 1.2))
            ),
            panel = function(x, y, z, subscripts, ...) {
                col <- ifelse(cint$ind, "gray60", "black")
                panel.segplot(x, y, z, subscripts,
                              col = col[subscripts],
                              ...)
                panel.abline(h = cint$real[subscripts], lty = 2)
            }),
    strip.left = strip.custom(
        factor.levels = fl,
        par.strip.text = list(cex = 1.5)
    ),
    strip = strip.custom(
       factor.levels =  yl
    )
)

#-----------------------------------------------------------------------
# Obtain the covariance between dispersion and regression parameters
aux <- ldply(
    lapply(results, function(x) {
    ind <- names(x); names(ind) <- ind
    out <- lapply(ind, function(y) {
        vcovs <- lapply(x[[y]], "[[", "vcov")
        index <- vapply(lapply(vcovs, is.na), sum, integer(1))
        vcovs <- vcovs[!index]
        corrs <- lapply(vcovs, cov2cor)
        covs <- do.call(rbind, lapply(vcovs, "[", 1, 2:5))
        cors <- do.call(rbind, lapply(corrs, "[", 1, 2:5))
        ldply(list("cov" = covs, "cor" = cors), .id = "scale")
    })
    ldply(out, .id = "phi")
}), .id = "n")
str(aux)

# Organize the results
dacov <- gather(aux, param, value, b0:b22, factor_key = TRUE)
dacov$phi <- ordered(dacov$phi,
                     c("phi=-1.6", "phi=-1", "phi=0", "phi=1.8"))
str(dacov)

# Distributions of covariance/correlations
fl <- parse(text = gsub("=", "==", levels(dacov$phi)))
yl <- parse(text = paste0("hat(beta)[", c(0, 1, 21, 22), "]"))

# Boxplots
useOuterStrips(
    bwplot(param ~ value | phi + scale,
           groups = n,
           data = dacov,
           box.width = 0.12,
           pch = "|",
           grid = TRUE,
           as.table = TRUE,
           layout = c(4, 2),
           xlab = "Observed values in the simulations",
           fill = colorRampPalette(c("gray80",  "gray20"))(4),
           panel = panel.superpose,
           scales = list(y = list(labels = yl), x = "free"),
           whisker.width = 0.4,
           auto.key = list(
               column = 4,
               rectangles = TRUE,
               points = FALSE,
               title = "Sample size",
               cex.title = 1
           ),
           par.settings = list(
               plot.symbol = list(pch = 19, cex = 0.1)
           ),
           panel.groups = function(x, y, group.number, ...) {
               my.panel.bwplot(x, y + (group.number - 2.5) / 5, ...)
               panel.abline(v = 0, lty = 2)
           }),
    strip = strip.custom(factor.levels = fl),
    strip.left = strip.custom(
        factor.levels = c("Covariance", "Correlation"))
)

#-------------------------------------------
# Times to fit
#+ eval=FALSE
library(purrr)
datime <-
    map_dfr(results_time, .id = "size", function(datan) {
        map_dfr(datan, .id = "phi", function(datap) {
            map_dfr(datap, function(sims) {
                tibble::tibble(time = sims[["ctime"]]["elapsed"])
            })
        })
    })
saveRDS(datime, "simultimes.rds")

#-------------------------------------------
datime <- readRDS("simultimes.rds")
datime %>%
    mutate(size = factor(sizes[size], levels = sizes),
           phi = factor(phi, levels = names(phis)[c(2, 3, 1, 4)])) %>%
    na.omit() %>%
    bwplot(time ~ size | phi,
           # groups = phi,
           layout = c(4, 1),
           scales = list(y = "free"),
           ylab = "Time (seconds)",
           xlab = "Sample size",
           axis = axis.grid,
           strip = strip.custom(
               factor.levels = parse(
                   text = gsub("=", "==", levels(.$phi)))
           ),
           data = .)

#-----------------------------------------------------------------------
# Likelihood surfaces (simulation with simple no covariate case)

# Settings of the simulation
fixed_mu <- 5
fixed_size <- 1000L

# Simulate 1000 observations of the DW distribution
set.seed(97380)
das <- lapply(phis, function(phi) {
    rcmp(fixed_size, mu = fixed_mu, phi = phi)
})

# Compute approximately expected dispersion indexes
dis <- vapply(phis, compute_variance, double(1),
              mu = fixed_mu, sumto = 500L) / 10
dis

# Verify mean, variance and dispersion indexes in the simulations
ldply(lapply(das, function(x) {
    c("Mean" = mean(x), "Var" = var(x), "DI" = var(x)/mean(x))
}), .id = "phi") %>%
    mutate(out, "Expected DI" = dis)

# Fit models
index <- seq_along(phis)
names(index) <- phis

# Fit on original parametrization
models_orpar <- lapply(index, function(i) {
    lambda <- mu2lambda(fixed_mu, phis[i])[[1]]
    start <- c("phi" = phis[i], "(Intercept)" = log(lambda))
    y <- das[[i]]
    fitcm(y ~ 1, sumto = 300L, model = "CP", start = start)
})

# Fit on propose reparametrization
models_repar <- lapply(index, function(i) {
    start <- c("phi2" = phis[i], "(Intercept)" = log(fixed_mu))
    y <- das[[i]]
    fitcm(y ~ 1, sumto = 300L, model = "CP2", start = start)
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
              "reallambda" = mu2lambda(fixed_mu, phis[i])[[1]],
              "realphi" = phis[i],
              "fitlambda" = exp(co[2]),
              "fitphi" = co[1])
})

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
              "realmu" = fixed_mu,
              "realphi" = phis[i],
              "fitmu" = exp(co[2]),
              "fitphi" = co[1])
})

# #-------------------------------------------
# # Plots deviance surfaces
# output <- readRDS("orthogonality.rds")
# devs_orpar <- output$devs_orpar
# devs_repar <- output$devs_repar

niveis <- c(0.9, 0.95, 0.99)
cortes <- qchisq(niveis, df = 2)
fl <- parse(text = paste("phi==", sort(phis)))

# Original parametrization
devs_orpar$fphi <- ordered(devs_orpar$realphi)
xy1 <- useOuterStrips(
    levelplot(
        deviance ~ lambda + phi | fphi + "Original parametrization",
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
    strip = strip.custom(factor.levels = fl)
)

# Propose parametrization
devs_repar$fphi <- ordered(devs_repar$realphi)
xy2 <- useOuterStrips(
    levelplot(
        deviance ~ mu + phi2 | fphi + "Proposed parametrization",
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
    strip = strip.custom(factor.levels = fl)
)

# Organize plots
print(xy1, position = c(0.0, 0.49, 1.0, 1.0), more = TRUE)
print(xy2, position = c(0.0, 0.0, 0.995, 0.51), more = FALSE)

#' # 91_defoliation-analysis #
# =====================================================================
#  Analysis of Artificial Defoliation in Cotton Plants
#                                                        Eduardo Junior
#                                                    edujrrib@gmail.com
#                                                            2017-08-05
# =====================================================================

#-----------------------------------------------------------------------
# Load package and functions
source("helper01_general-functions.R")
source("helper02_lattice-panels.R")

library(bbmle)
library(multcomp)
library(dplyr)

# Colors for legends
cols <- trellis.par.get("superpose.line")$col

#-----------------------------------------------------------------------
# Load data
# data(cottonBolls, package = "cmpreg")
cottonBolls <- read.table("../data/cottonBolls.txt",
                          header = TRUE, sep = "\t")
str(cottonBolls)

#-----------------------------------------------------------------------
# Exploratory analysis

# Scatter plot with a beewarm taste
xy1 <- xyplot(ncap ~ des | est,
              data = cottonBolls,
              layout = c(2, 3),
              as.table = TRUE,
              grid = TRUE,
              type = c("p", "smooth"),
              xlab = "Artificial defoliation level",
              ylab = "Number of bolls produced",
              spread = 0.05,
              panel = panel.beeswarm)

# Sample variance vs sample mean (evidence in favor of the
# underdispersion).
mv <- cottonBolls %>%
    group_by(est, des) %>%
    summarise(mu = mean(ncap), va = var(ncap))
xlim <- ylim <- extendrange(c(mv$mu, mv$va), f = 0.05)

xy2 <- xyplot(va ~ mu,
              data = mv,
              grid = TRUE,
              type = c("p", "r"),
              xlim = xlim,
              ylim = ylim,
              xlab = expression("Sample"~"mean"~(bar(y))),
              ylab = expression("Sample"~"variance"~(s^2)),
              panel = function(...) {
                  panel.xyplot(...)
                  panel.abline(a = 0, b = 1, lty = 2)
              })

print(xy1, split = c(1, 1, 2, 1), more = TRUE)
print(xy2, split = c(2, 1, 2, 1), more = FALSE)

#-----------------------------------------------------------------------
# Fit models
mnames <- c("PO", "C1", "C2", "QP")

# Predictor, following Zeviani et al. (2014)
form1 <- ncap ~ est:(des + I(des^2))

m1PO <- glm(form1, data = cottonBolls, family = poisson)
system.time(
    m1C1 <- fitcm(form1, data = cottonBolls, model = "CP", sumto = 50)
)
system.time(
    m1C2 <- fitcm(form1, data = cottonBolls, model = "CP2", sumto = 50)
)
m1QP <- glm(form1, data = cottonBolls, family = quasipoisson)

models.ncap <- list(m1PO, m1C1, m1C2, m1QP)
names(models.ncap) <- mnames

# Numbers of calls to loglik and numerical gradient
models.ncap$C1@details$counts
models.ncap$C2@details$counts

#-----------------------------------------------------------------------
# LRT and profile extra parameter

# LRT between Poisson and COM-Poisson (test: phi == 0)
getAnova(m1PO, m1C2)

profs.ncap <- lapply(list(c(m1C1, "phi"), c(m1C2, "phi2")),
                     function(x) myprofile(x[[1]], x[[2]]))
profs.ncap <- do.call("rbind", profs.ncap)

# Plot profile extra parameter
snames <- parse(text = c("'COM-Poisson'~(phi)",
                         "'COM-Poisson'[mu]~(phi)"))
xyprofile(profs.ncap, namestrip = snames,
          ylim = extendrange(c(0, 3)),
          xlab = "Precision parameter")

#-----------------------------------------------------------------------
# Goodness of fit measures and estimate parameters

# GoF measures
measures.ncap <- sapply(models.ncap, function(x)
    c("LogLik" = logLik(x), "AIC" = AIC(x), "BIC" = BIC(x)))
measures.ncap

# Get the estimates
est <- lapply(models.ncap, FUN = function(x) getCoefs(x))
est.ncap <- do.call(cbind, est)
est.ncap

#-----------------------------------------------------------------------
# Prediction

# Data for prediction
pred <- with(cottonBolls,
             expand.grid(
                 est = levels(est),
                 des = seq(min(des), max(des), l = 20)
             ))
qn <- qnorm(0.975) * c(fit = 0, lwr = -1, upr = 1)

# Design matrix for prediction
X <- model.matrix(update(form1, NULL~.), pred)

# Considering Poisson
aux <- exp(confint(
    glht(m1PO, linfct = X), calpha = univariate_calpha())$confint)
colnames(aux) <- c("fit", "lwr", "upr")
aux <- data.frame(modelo = "Poisson", aux)
predPO.ncap <- cbind(pred, aux)

# Considering COM-Poisson
aux <- predictcm(m1C1, newdata = X)
aux <- data.frame(modelo = "COM-Poisson", aux)
predC1.ncap <- cbind(pred, aux)

# Considering COM-Poisson (mean parametrization)
aux <- predictcm(m1C2, newdata = X)
aux <- data.frame(modelo = "COM-Poisson2", aux)
predC2.ncap <- cbind(pred, aux)

# Considering Quasi-Poisson
aux <- exp(confint(
    glht(m1QP, linfct = X), calpha = univariate_calpha())$confint)
colnames(aux) <- c("fit", "lwr", "upr")
aux <- data.frame(modelo = "Quasi-Poisson", aux)
predQP.ncap <- cbind(pred, aux)

# Representing the confidence intervals
pred.ncap <- rbind(predPO.ncap, predC1.ncap, predC2.ncap, predQP.ncap)

# Legend
key <- list(columns = 2,
            lines = list(col = 1, lty = rev(1:4)),
            text = list(parse(
                text = c("'Poisson'", "'COM-Poisson'",
                         "'COM-Poisson'[mu]", "'Quasi-Poisson'"))
                ))

# Graph
update(xy1, layout = c(NA, 1), type = c("g", "p"),
       alpha = 0.6, key = key) +
    as.layer(
        xyplot(fit ~ des | est,
               auto.key = TRUE,
               data = pred.ncap,
               groups = modelo,
               type = "l",
               layout = c(NA, 1),
               as.table = TRUE,
               col = 1,
               ly = pred.ncap$lwr,
               uy = pred.ncap$ upr,
               cty = "bands",
               fill = "gray80",
               alpha = 0.1,
               panel = panel.superpose,
               panel.groups = panel.cbH,
               prepanel = prepanel.cbH,
               lty = rev(1:4))
    )

#-----------------------------------------------------------------------
# Correlation between estimates
corr.ncap <- purrr::map_df(models.ncap[c("C1", "C2")],
                           function(x) cov2cor(vcov(x))[1, -1])
round(corr.ncap, 5)
#' # 92_soybean-analysis #
# =====================================================================
# Analysis of number of grains in Soyabean under Soil Moisture and
# Potassium Doses
#                                                        Eduardo Junior
#                                                    edujrrib@gmail.com
#                                                            2017-08-05
# =====================================================================

#-----------------------------------------------------------------------
# Load package and functions
source("helper01_general-functions.R")
source("helper02_lattice-panels.R")

library(plyr)
library(tidyr)
library(bbmle)

# Colors for legends
cols <- trellis.par.get("superpose.line")$col

#-----------------------------------------------------------------------
# Load data
# data(soyaBeans, package = "cmpreg")
soyaBeans <- read.table("../data/soyaBeans.txt",
                        header = TRUE, sep = "\t")
soyaBeans$umid <- as.factor(soyaBeans$umid)
soyaBeans <- soyaBeans[-74, ] ## Incorrect observation
soyaBeans <- transform(soyaBeans, K = K / 100)

#-----------------------------------------------------------------------
# Exploratory analysis

# Scatter plot
xy1 <- xyplot(ngra ~ K | umid,
              data = soyaBeans,
              xlab = "Potassium fertilization level",
              ylab = "Number of grains per pot",
              type = c("p", "g", "smooth"),
              as.table =  TRUE,
              layout = c(2, 2),
              strip = strip.custom(
                  strip.names = TRUE, var.name = "humidity",
                  factor.levels = paste0(levels(soyaBeans$umid), "%")))

# Sample variance vs sample mean (evidence in favor of the
# overdispersion).
mv <- soyaBeans %>%
    group_by(K, umid) %>%
    summarise(mu = mean(ngra), va = var(ngra))
xlim <- ylim <- extendrange(c(mv$mu, mv$va), f = 0.05)

xy2 <- xyplot(va ~ mu,
              data = mv,
              type = c("p", "r", "g"),
              xlim = xlim,
              ylim = ylim,
              xlab = expression("Sample"~"mean"~(bar(y))),
              ylab = expression("Sample"~"variance"~(s^2)),
              panel = function(...) {
                  panel.xyplot(...)
                  panel.abline(a = 0, b = 1, lty = 2)
              })

print(xy1, split = c(1, 1, 2, 1), more = TRUE)
print(xy2, split = c(2, 1, 2, 1), more = FALSE)

#-----------------------------------------------------------------------
# Fit models
mnames <- c("PO", "C1", "C2", "QP")

# Predictor
form2 <-  ngra ~ bloc + umid * K + I(K^2)

m2PO <- glm(form2, data = soyaBeans, family = poisson)
system.time(
    m2C1 <- fitcm(form2, data = soyaBeans, model = "CP", sumto = 700)
)
system.time(
    m2C2 <- fitcm(form2, data = soyaBeans, model = "CP2", sumto = 700)
)
m2QP <- glm(form2, data = soyaBeans, family = quasipoisson)

models.ngra <- list(m2PO, m2C1, m2C2, m2QP)
names(models.ngra) <- mnames

## Numbers of calls to loglik and numerical gradient
models.ngra$C1@details$counts
models.ngra$C2@details$counts

#-----------------------------------------------------------------------
# LRT and profile extra parameter

# LRT between Poisson and COM-Poisson (test: phi == 0)
getAnova(m2PO, m2C2)

profs.ngra <- lapply(list(c(m2C1, "phi"), c(m2C2, "phi2")),
                     function(x) myprofile(x[[1]], x[[2]]))
profs.ngra <- do.call("rbind", profs.ngra)

# Plot profile extra parameter
snames <- parse(text = c("'COM-Poisson'~(phi)",
                         "'COM-Poisson'[mu]~(phi)"))
xyprofile(profs.ngra, namestrip = snames,
          ylim = extendrange(c(0, 3)),
          xlab = "Precision parameter",
          scales.x = NULL)

#-----------------------------------------------------------------------
# Goodness of fit measures and estimate parameters

# GoF measures
measures.ngra <- sapply(models.ngra, function(x)
    c("LogLik" = logLik(x), "AIC" = AIC(x), "BIC" = BIC(x)))
measures.ngra

# Get the estimates
est <- lapply(models.ngra, FUN = function(x) getCoefs(x))
est.ngra <- do.call(cbind, est)
est.ngra

#-----------------------------------------------------------------------
# Prediction

# Data for prediction
pred <- with(soyaBeans,
             expand.grid(
                 bloc = factor(levels(bloc)[1], levels = levels(bloc)),
                 umid = levels(umid),
                 K = seq(min(K), max(K), l = 20)
             ))
qn <- qnorm(0.975) * c(fit = 0, lwr = -1, upr = 1)

# Design matrix for prediction
X <- model.matrix(update(form2, NULL ~ .), pred)
bl <- attr(X, "assign") == 1
X[, bl] <- X[, bl] + 1/(sum(bl) + 1)

# Considering Poisson
aux <- exp(confint(
    glht(m2PO, linfct = X), calpha = univariate_calpha())$confint)
colnames(aux) <- c("fit", "lwr", "upr")
aux <- data.frame(modelo = "Poisson", aux)
predPO.ngra <- cbind(pred, aux)

# Considering COM-Poisson
aux <- predictcm(m2C1, newdata = X)
aux <- data.frame(modelo = "COM-Poisson", aux)
predC1.ngra <- cbind(pred, aux)

# Considering COM-Poisson (mean parametrization)
aux <- predictcm(m2C2, newdata = X)
aux <- data.frame(modelo = "COM-Poisson2", aux)
predC2.ngra <- cbind(pred, aux)

# Considering Quasi-Poisson
aux <- exp(confint(
    glht(m2QP, linfct = X), calpha = univariate_calpha())$confint)
colnames(aux) <- c("fit", "lwr", "upr")
aux <- data.frame(modelo = "Quasi-Poisson", aux)
predQP.ngra <- cbind(pred, aux)

# Representing the confidence intervals
pred.ngra <- rbind(predPO.ngra, predC1.ngra, predC2.ngra, predQP.ngra)

# Legend
key <- list(columns = 2,
            lines = list(col = 1, lty = rev(1:4)),
            text = list(parse(
                text = c("'Poisson'", "'COM-Poisson'",
                         "'COM-Poisson'[mu]", "'Quasi-Poisson'"))
                )
            )

# Graph
update(xy1, layout = c(NA, 1), type = c("p", "g"),
       alpha = 0.6, key = key) +
    as.layer(
        xyplot(fit ~ K | umid,
               data = pred.ngra,
               groups = modelo,
               type = "l",
               col = 1,
               ly = pred.ngra$lwr,
               uy = pred.ngra$upr,
               cty = "bands",
               fill = "gray80",
               alpha = 0.1,
               panel = panel.superpose,
               panel.groups = panel.cbH,
               prepanel = cmpreg::prepanel.cbH,
               lty = rev(1:4))
    )

#-----------------------------------------------------------------------
# Correlation between estimates
corr.ngra <- purrr::map_df(models.ngra[c("C1", "C2")],
                           function(x) cov2cor(vcov(x))[1, -1])
round(corr.ngra, 5)

#' # 93_nitrofen-analysis #
# =====================================================================
# Analysis of number of grains in Soyabean under Soil Moisture and
# Potassium Doses
#                                                        Eduardo Junior
#                                                    edujrrib@gmail.com
#                                                            2017-08-05
# =====================================================================

#-----------------------------------------------------------------------
# Load package and functions
source("helper01_general-functions.R")
source("helper02_lattice-panels.R")

library(plyr)
library(tidyr)
library(bbmle)

# Colors for legends
cols <- trellis.par.get("superpose.line")$col

#-----------------------------------------------------------------------
# Load data
# data(nitrofen, package = "boot")
# data(Paula, package = "labestData")
# nitrofen <- PaulaEx4.6.20
nitrofen <- read.table("../data/nitrofen.txt",
                       header = TRUE, sep = "\t")
nitrofen <- transform(nitrofen, dose = dose / 100)

#-----------------------------------------------------------------------
# Exploratory analysis

# Scatter plot
xy1 <- xyplot(novos ~ dose,
              data = nitrofen,
              xlab = "Nitrofen concentration level",
              ylab = "Number of live offspring",
              type = c("p", "g", "smooth"),
              spread = 0.05,
              panel = panel.beeswarm)

# Sample variance vs sample mean (evidence in favor of the
# no relationship between mean and variance).
mv <- nitrofen %>%
    group_by(dose) %>%
    summarise(mu = mean(novos), va = var(novos))
xlim <- ylim <- extendrange(c(mv$mu, mv$va), f = 0.05)

xy2 <- xyplot(va ~ mu,
              data = mv,
              type = c("p", "r", "g"),
              xlim = xlim,
              ylim = ylim,
              xlab = expression("Sample"~"mean"~(bar(y))),
              ylab = expression("Sample"~"variance"~(s^2)),
              panel = function(...) {
                  panel.xyplot(...)
                  panel.abline(a = 0, b = 1, lty = 2)
              })

print(xy1, split = c(1, 1, 2, 1), more = TRUE)
print(xy2, split = c(2, 1, 2, 1), more = FALSE)

#-----------------------------------------------------------------------
# Fit models
mnames <- c("PO", "C1", "C2", "QP")

# Predictors
form31 <-  novos ~ dose
form32 <-  novos ~ dose + I(dose^2)
form33 <-  novos ~ dose + I(dose^2) + I(dose^3)

predictors <- list("pred1" = form31, "pred2" = form32, "pred3" = form33)
fmodels.ovos <- lapply(predictors, function(form) {
    PO <- glm(form, data = nitrofen, family = poisson)
    C1 <- fitcm(form, data = nitrofen, model = "CP", sumto = 100)
    C2 <- fitcm(form, data = nitrofen, model = "CP2", sumto = 100)
    QP <- glm(form, data = nitrofen, family = quasipoisson)
    list("PO" = PO, "C1" = C1, "C2" = C2, "QP" = QP)
})

#-----------------------------------------------------------------------
# LRT for nested models

# Poisson
auxPO <- lapply(fmodels.ovos, function(x) x$PO)
do.call("getAnova", auxPO)

# COM-Poisson standard
auxC1 <- lapply(fmodels.ovos, function(x) x$C1)
do.call("getAnova", auxC1)

# COM-Poisson mean-parameterized
auxC2 <- lapply(fmodels.ovos, function(x) x$C2)
do.call("getAnova", auxC2)

# Quasi-Poisson
auxQP <- lapply(fmodels.ovos, function(x) x$QP)
do.call("getAnova", auxQP)

#--------------------------------------------
# Separe the choose models
form3 <- form33
m3PO <- fmodels.ovos$pred3$PO
m3C1 <- fmodels.ovos$pred3$C1
m3C2 <- fmodels.ovos$pred3$C2
m3QP <- fmodels.ovos$pred3$QP

models.ovos <- list(m3PO, m3C1, m3C2, m3QP)
names(models.ovos) <- mnames

# Numbers of calls to loglik and numerical gradient
models.ovos$C1@details$counts
models.ovos$C2@details$counts

#-----------------------------------------------------------------------
# LRT and profile extra parameter

# LRT between Poisson and COM-Poisson (test: phi == 0)
getAnova(m3PO, m3C2)

profs.novos <- lapply(list(c(m3C1, "phi"), c(m3C2, "phi2")),
                     function(x) myprofile(x[[1]], x[[2]]))
profs.novos <- do.call("rbind", profs.novos)

# Plot profile extra parameter
snames <- parse(text = c("'COM-Poisson'~(phi)",
                         "'COM-Poisson'[mu]~(phi)"))
xyprofile(profs.novos, namestrip = snames,
          ylim = extendrange(c(0, 3)),
          xlab = "Precision parameter",
          scales.x = NULL) +
    layer(panel.abline(v = 0, col = 2, lty = 2))

#-----------------------------------------------------------------------
# Goodness of fit measures and estimate parameters

# GoF measures
measures.ovos <- sapply(models.ovos, function(x)
    c("LogLik" = logLik(x), "AIC" = AIC(x), "BIC" = BIC(x)))
measures.ovos

# Get the estimates
est <- lapply(models.ovos, FUN = function(x) getCoefs(x))
est.ovos <- do.call(cbind, est)
est.ovos

#-----------------------------------------------------------------------
# Prediction

# Data for prediction
pred <- with(nitrofen,
             data.frame("dose" = seq(min(dose), max(dose),
                                     length.out = 100)))
qn <- qnorm(0.975) * c(fit = 0, lwr = -1, upr = 1)

# Design matrix for prediction
X <- model.matrix(update(form3, NULL ~ .), pred)

# Considering Poisson
aux <- exp(confint(
    glht(m3PO, linfct = X), calpha = univariate_calpha())$confint)
colnames(aux) <- c("fit", "lwr", "upr")
aux <- data.frame(modelo = "Poisson", aux)
predPO.novos <- cbind(pred, aux)

# Considering COM-Poisson
aux <- predictcm(m3C1, newdata = X)
aux <- data.frame(modelo = "COM-Poisson", aux)
predC1.novos <- cbind(pred, aux)

# Considering COM-Poisson (mean parametrization)
aux <- predictcm(m3C2, newdata = X)
aux <- data.frame(modelo = "COM-Poisson2", aux)
predC2.novos <- cbind(pred, aux)

# Considering Quasi-Poisson
aux <- exp(confint(
    glht(m3QP, linfct = X), calpha = univariate_calpha())$confint)
colnames(aux) <- c("fit", "lwr", "upr")
aux <- data.frame(modelo = "Quasi-Poisson", aux)
predQP.novos <- cbind(pred, aux)

# Representing the confidence intervals
pred.novos <- rbind(predPO.novos, predC1.novos, predC2.novos, predQP.novos)
ord <- order(pred.novos$dose, pred.novos$modelo)
pred.novos <- pred.novos[ord, ]

# Legend
key <- list(columns = 2,
            lines = list(col = 1, lty = rev(1:4)),
            text = list(parse(
                text = c("'Poisson'", "'COM-Poisson'",
                         "'COM-Poisson'[mu]", "'Quasi-Poisson'"))
                )
            )

# Graph
update(xy1, type = c("g", "p"), alpha = 0.6, key = key) +
    as.layer(
        xyplot(fit ~ dose,
               auto.key = TRUE,
               data = pred.novos,
               groups = modelo,
               type = "l",
               col = 1,
               ly = pred.novos$lwr,
               uy = pred.novos$upr,
               cty = "bands",
               fill = "gray80",
               alpha = 0.1,
               panel = panel.superpose,
               panel.groups = panel.cbH,
               prepanel = prepanel.cbH,
               lty = rev(1:4))
    )

#-----------------------------------------------------------------------
# Correlation between estimates
corr.ovos <- purrr::map_df(models.ovos[c("C1", "C2")],
                           function(x) cov2cor(vcov(x))[1, -1])
round(corr.ovos, 5)

#' # 94_timestofit-study #
# =====================================================================
# Comparasion of the computational times
#                                                        Eduardo Junior
#                                                    edujrrib@gmail.com
#                                                            2017-08-05
# =====================================================================

#-----------------------------------------------------------------------
# Load package and functions
source("helper01_general-functions.R")
source("helper02_lattice-panels.R")

library(bbmle)
library(microbenchmark)

# Colors for legends
cols <- trellis.par.get("superpose.line")$col

#-----------------------------------------------------------------------
# Load datas

#-------------------------------------------
# Section 6.1
data1 <- read.table("../data/cottonBolls.txt",
                    header = TRUE, sep = "\t")
form1 <- ncap ~ est:(des + I(des^2))

#-------------------------------------------
# Section 6.2
data2 <- read.table("../data/soyaBeans.txt",
                        header = TRUE, sep = "\t")
data2$umid <- as.factor(data2$umid)
data2 <- data2[-74, ] ## Incorrect observation
data2 <- transform(data2, K = K / 100)
form2 <- ngra ~ bloc + umid * K + I(K^2)

#-------------------------------------------
# Section 6.3
data3 <- read.table("../data/nitrofen.txt",
                       header = TRUE, sep = "\t")
data3 <- transform(data3, dose = dose / 100)
form3 <- novos ~ dose + I(dose^2) + I(dose^3)

#-----------------------------------------------------------------------
# Compare computation times by 50 repetitions

# Times to fit case studies
#+ eval=FALSE
bench1 <- microbenchmark(
    "CMP  " = fitcm(form1, data = data1, model = "CP" , sumto = 50),
    "CMPmu" = fitcm(form1, data = data1, model = "CP2", sumto = 50),
    times = 50)
bench2 <- microbenchmark(
    "CMP  " = fitcm(form2, data = data2, model = "CP" , sumto = 700),
    "CMPmu" = fitcm(form2, data = data2, model = "CP2", sumto = 700),
    times = 50)
bench3 <- microbenchmark(
    "CMP  " = fitcm(form3, data = data3, model = "CP" , sumto = 100),
    "CMPmu" = fitcm(form3, data = data3, model = "CP2", sumto = 100),
    times = 50)

#-----------------------------------------------------------------------
# Organize results
#+ eval=FALSE
benchs <- list("Cotton (under)" = bench1,
               "Soybean (over)" = bench2,
               "Nitrofen (equi)" = bench3)
benchs <- purrr::map_dfr(benchs, identity, .id = "case")

class(benchs) <- "data.frame"
saveRDS(benchs, "comparetimes.rds")

#-------------------------------------------
benchs <- readRDS("comparetimes.rds")

xlabs <- expression("COM-Poisson", "COM-Poisson"[mu])
bwplot(time ~ expr | case,
       ylab = "Time (seconds)",
       scales = list(
           y = list(relation = "free"),
           x = list(at = 1:2, labels = xlabs)
       ),
       data = benchs)

#' # helper01_general-functions #
# =====================================================================
# General functions
#                                                        Eduardo Junior
#                                                    edujrrib@gmail.com
#                                                            2017-08-01
# =====================================================================

#=======================================================================
# Tranformations
lambda2mu <- function(lambda, nu) {
    phi <- log(nu)
    mu <- lambda^(1/nu) - (nu - 1)/(2 * nu)
    list("mu" = mu, "phi" = phi)
}

mu2lambda <- function(mu, phi) {
    nu <- exp(phi)
    lambda <- (mu + (nu-1)/(2*nu))^nu
    list("lambda" = lambda, "nu" = nu)
}

#=======================================================================
# Negative log-likelihood functions

# COM-Poisson Model (Original parameterization)
llcmp <- function (params, y, X, offset = NULL, sumto = NULL) {
    if (!is.null(offset)) {
        stop("This model doesn't support offset yet.")
    }
    phi <- params[1]
    nu <- exp(phi)
    Xb <- X %*% params[-1]
    i <- 0:sumto
    zs <- sapply(Xb, function(loglambda) {
        sum(exp(i * loglambda - nu * lfactorial(i)))
    })
    Z <- sum(log(zs))
    ll <- sum(y * Xb - nu * lfactorial(y)) - Z
    return(-ll)
}

# COM-Poisson Model (mean-parameterization)
llcmp2 <- function (params, y, X, offset = NULL, sumto = NULL) {
    if (!is.null(offset)) {
        stop("This model doesn't support offset yet.")
    }
    phi2 <- params[1]
    nu <- exp(phi2)
    #-------------------------------------------
    mu <- exp(X %*% params[-1])
    Xb <- nu * log(mu + (nu - 1)/(2 * nu))
    #-------------------------------------------
    i <- 0:sumto
    zs <- sapply(Xb, function(loglambda) {
        sum(exp(i * loglambda - nu * lfactorial(i)))
    })
    Z <- sum(log(zs))
    ll <- sum(y * Xb - nu * lfactorial(y)) - Z
    return(-ll)
}

#=======================================================================
# Probability mass function (mean parametrization)
dcmp <- function(y, mu, phi, sumto = 500L) {
    nu <- exp(phi)
    k0 <- nu * log(mu + (nu - 1)/(2 * nu))
    j <- 1:sumto
    Z <- 1 + sum(exp(j * k0 - nu * lfactorial(j)))
    l <- y * k0 - nu * lfactorial(y) - log(Z)
    return(exp(l))
}

#=======================================================================
# Compute momentes mean and variance
compute_mean <- function(mu, phi, sumto = NULL, upto = 500) {
    y <- 0:upto
    names(mu) <- NULL
    names(phi) <- NULL
    pars <- data.frame(mu = mu, phi = phi)
    apply(pars, 1L, function(p) {
        py <- dcmp(y, p[1L], p[2L])
        sum(y * py)
    })
}

compute_variance <- function(mu, phi, sumto = NULL, upto = 500) {
    y <- 0:upto
    names(mu) <- NULL
    names(phi) <- NULL
    pars <- data.frame(mu = mu, phi = phi)
    apply(pars, 1L, function(p) {
        py <- dcmp(y, p[1L], p[2L])
        ey <- sum(y * py)
        sum((y - ey)^2 * py)
    })
}

calc_moments <- function(mu, phi, sumto = 500, upto = 500) {
    y <- 0:upto
    names(mu) <- NULL
    names(phi) <- NULL
    pars <- data.frame(mu = mu, phi = phi)
    apply(pars, 1L, function(p) {
        py <- dcmp(y, p[1L], p[2L], sumto = sumto)
        ey <- sum(y * py)
        e2 <- sum(y^2 * py)
        vy <- e2 - ey^2
        c("mean" = ey, "var" = vy)
    })
}

#=======================================================================
# Framework for fit count models (using bbmle::mle2 function)
fitcm <- function(formula, data, start = NULL,
                  model = c("CP", "CP2"), ...) {
    dots <- list(...)
    model <- match.arg(model)
    if (missing(data)) {
        data <- environment(formula)
    }
    #-------------------------------------------
    frame <- model.frame(formula, data)
    terms <- attr(frame, "terms")
    y <- model.response(frame)
    X <- model.matrix(terms, frame)
    off <- model.offset(frame)
    if (is.null(start)) {
        m0 <- glm.fit(x = X, y = y, offset = off, family = poisson())
        start <- c(0, m0$coefficients)
    }
    #-------------------------------------------
    if (model == "CP") {
        names(start)[1] <- "phi"
        bbmle::parnames(llcmp) <- names(start)
        sumto <- dots$sumto
        dots$sumto <- NULL
        args <- append(
            list(minuslog = llcmp, start = start,
                 data = list(y = y, X = X, offset = off,
                             sumto = sumto, model = model),
                 vecpar = TRUE), dots)
        fit <- suppressWarnings(
            do.call(what = bbmle::mle2, args = args))
    }
    if (model == "CP2") {
        names(start)[1] <- "phi2"
        bbmle::parnames(llcmp2) <- names(start)
        sumto <- dots$sumto
        dots$sumto <- NULL
        args <- append(
            list(minuslog = llcmp2, start = start,
                 data = list(y = y, X = X, offset = off,
                             sumto = sumto, model = model),
                 vecpar = TRUE), dots)
        fit <- suppressWarnings(
            do.call(what = bbmle::mle2, args = args))
    }
    slot(fit, "call.orig") <- match.call()
    return(fit)
}

#=======================================================================
# Analysis of Deviance Table to glm and mle2 objetcs
getAnova <- function(..., print = TRUE) {
    model.list <- list(...)
    test <- "Chisq"
    cls <- sapply(model.list, function(x) class(x)[1])
    if (all(cls == "glm")) {
        fam <- sapply(model.list, function(x) x$family$family)
        if (all(fam == "quasipoisson")) {
            test <- "F"
        }
    }
    nps <- sapply(model.list, function(x) attr(logLik(x), "df"))
    aic <- sapply(model.list, function(x) AIC(x))
    if (test == "Chisq") {
        lls <- sapply(model.list, function(x) logLik(x))
        cst <- 2 * diff(lls)
        pvs <- pchisq(cst, df = diff(nps), lower.tail = FALSE)
        tab <- cbind(
            "np" = nps,
            "ll" = lls,
            "AIC" = aic,
            "2*dif" = c(NA, cst),
            "dif np" = c(NA, diff(nps)),
            "Pr(>Chisq)" = c(NA, pvs))
    }
    if (test == "F") {
        ## ver pág. 187 expresssão 7.10 em Modern Applied, Ripley
        sigma <- sapply(model.list, function(x) summary(x)$dispersion)
        df.sigma <- sapply(model.list, function(x) x$df.residual)
        dev <- sapply(model.list, function(x) deviance(x))
        dfs <- c(NA, -diff(dev))
        dnp <- c(NA, diff(nps))
        F <- c(NA, -diff(dev))/(dnp * sigma)
        pvs <- pf(F, df1 = dnp, df2 = df.sigma, lower.tail = FALSE)
        tab <- cbind(
            "np" = nps + 1,
            "dev" = dev,
            "AIC" = aic,
            "F" = F,
            "dif np" = dnp,
            "Pr(>F)" = pvs)
    }
    rownames(tab) <- 1:length(model.list)
    if (print) printCoefmat(tab, na.print = "", cs.ind = 1)
    invisible(tab)
}

# Get the estimated coefficients of the parameters
getCoefs <- function(model, digits = 3, rownames = NULL) {
    #-------------------------------------------
    glue <- function(est, d = digits) {
        est
        # apply(est, 1, FUN = function(x)
        #     paste0(format(x[1], digits = d, nsmall = d), " (",
        #            format(x[2], digits = d, nsmall = d), ")"))
    }
    #-------------------------------------------
    cls <- class(model)[1]
    est <- switch(
        cls,
        "glm" = {
            fam <- model$family$family
            if (!fam %in% c("poisson", "quasipoisson")) {
                stop("Família de distribuições não contemplada")
            }
            switch(
                fam,
                "poisson" = {
                    glue(rbind(c(NA, NA),
                               summary(model)$coef[, c(1, 3)]))
                },
                "quasipoisson" = {
                    glue(rbind(c(summary(model)$dispersion, NA),
                          summary(model)$coef[, c(1, 3)]))
                }
            )
        },
        "mle2" = {
            glue(bbmle::summary(model)@coef[, c(1, 3)])
        }
    )
    if (!is.null(rownames)) {
        rownames(est) <- rownames
    }
    colnames(est) <- c("Est", "Est/SE")
    return(est)
}

#=======================================================================
# Predict Values

# Get the interval values of linear predictors (for Delta method with
# cholesky decomposition)
cholV_eta <- function(V, X, b, qn) {
    eta <- c(X %*% b)
    if (length(qn) > 1) {
        U <- chol(V)
        stderr <- sqrt(
            apply(X %*% t(U), MARGIN = 1, FUN = function(x) sum(x^2))
        )
        eta <- sweep(x = outer(stderr, qn, FUN = "*"),
                     MARGIN = 1,
                     STATS = eta,
                     FUN = "+")
    }
    return(eta)
}

# Get the means (inverse link function)
calc_mean <- function(eta, dispersion, model = c("CP", "CP2"),
                      tol = 1e-5, offset = NULL, ...) {
    if (is.null(offset)) {
        offset <- 1L
    }
    model <- match.arg(model)
    names(eta) <- NULL
    names(dispersion) <- NULL
    pars <- data.frame(eta = eta, dispersion = dispersion)
    ##-------------------------------------------
    if (model == "CP") {
        ##-------------------------------------------
        pars <- transform(pars, lambda = exp(eta),
                          nu = exp(dispersion))
        sumto <- list(...)$sumto
        ##-------------------------------------------
        approxmu <- with(pars, lambda^(1/nu) - (nu - 1)/(2 * nu))
        sigma <- with(pars, (1/nu) * approxmu)
        sigma <- ifelse(sigma < 0, 0, sigma)
        ymax <- ceiling(max(approxmu + 5 * sqrt(sigma)))
        etamax <- max(pars$eta)
        dispersionmin <- min(pars$dispersion)
        pmax <- exp(
            -llcmp(params = c(dispersionmin, etamax), y = ymax, X = 1,
                   sumto = sumto))
        while (pmax > tol) {
            ymax <- ymax + 1
            pmax <- exp(
                -llcmp(params = c(dispersionmin, etamax),
                        y = ymax, X = 1, sumto = sumto))
        }
        y <- 1:ymax
        estmean <- mapply(
            eta = pars$eta,
            nu = pars$nu,
            MoreArgs = list(y = y),
            FUN = function(eta, nu, y) {
                i <- 0:sumto
                zs <- sapply(eta, function(eta) {
                    sum(exp(i * eta - nu * lfactorial(i)))
                })
                py <- exp(y * eta - nu * lfactorial(y) - sum(log(zs)))
                sum(y * py)
            },
            SIMPLIFY = TRUE)
    }
    if (model == "CP2") {
        estmean <- exp(eta)
    }
    names(estmean) <- NULL
    return(c(estmean))
}

# Calculate the confidence interval of predict values
predictcm <- function(object, newdata, offset = NULL, level = 0.95) {
    #-------------------------------------------
    if (missing(newdata)) {
        X <- object@data$X
    } else {
        X <- newdata
    }
    #-------------------------------------------
    qn <- -qnorm((1 - level[1])/2)
    qn <- qn * c(fit = 0, lwr = -1, upr = 1)
    #--------------------------------------------
    V <- vcov(object)
    Vc <- V[-1, -1] - V[-1, 1] %*% solve(V[1, 1]) %*% V[1, -1]
    # Vc <- V[-1, -1]
    eta <- cholV_eta(Vc, X, b = coef(object)[-1], qn = qn)
    #-------------------------------------------
    args <- list(dispersion = coef(object)[1],
                 model = object@data$model,
                 offset = offset)
    if (object@data$model == "CP") {
        args$sumto <- object@data$sumto
    }
    pred <- apply(as.matrix(eta), MARGIN = 2,
                  FUN = function(x) {
                      do.call(what = calc_mean,
                              args = append(list(eta = x), args))
                  })
    colnames(pred) <- names(qn)
    return(pred)
}

# To profile log-likelihood
myprofile <- function(model, which) {
    sel <- c("param", "z", "focal")
    prof <- profile(model, which = which)
    as.data.frame(prof)[, sel]
}

#' # helper02_lattice-panels #
# =====================================================================
# Auxiliary funtions for lattice graphics
#                                                        Eduardo Junior
#                                                    edujrrib@gmail.com
#                                                            2017-08-01
# =====================================================================

#-----------------------------------------------------------------------
# Lattice configuration
library(lattice)
library(latticeExtra)
library(grid)

## http://www.magesblog.com/2013/04/how-to-change-alpha-value-of-colours-in.html
add.alpha <- function(col, alpha = 1){
    apply(sapply(col, col2rgb)/255, 2,
          function(x) rgb(x[1], x[2], x[3], alpha = alpha))
}

## Define Colors
mycol <- colorRampPalette(c("gray80", "gray10"))(4)
myreg <- colorRampPalette(c("gray90", "gray20"))(100)

## Trellis graphical style.
ps1 <- list(
    superpose.symbol = list(
        col = mycol, pch = 1,
        fill = add.alpha(mycol, alpha = 0.4)),
    box.rectangle = list(col = 1, fill = c("gray70")),
    box.umbrella = list(col = 1, lty = 1),
    box.dot = list(pch = "|"),
    dot.symbol = list(col = 1, pch = 19),
    dot.line = list(col = "gray50", lty = 3),
    plot.symbol = list(col = 1),
    plot.line = list(col = 1),
    plot.polygon = list(col = "gray95"),
    superpose.line = list(col = mycol, lty = 1),
    superpose.polygon = list(col = mycol),
    strip.background = list(col = c("gray90", "gray70")),
    regions = list(col = myreg),
    par.sub.text = list(
        font = 1, just = "left", cex = 0.9,
        x = grid::unit(10, "mm"))
    )

# Space economic settings
ac <- list(pad1 = 0.5, pad2 = 0.5, tck = 0.5)
ps2 <- list(
    layout.widths = list(
        left.padding = 0.25,
        right.padding = -1,
        ylab.axis.padding = 0),
    layout.heights = list(
        bottom.padding = 0.25,
        top.padding = 0,
        axis.xlab.padding = 0,
        xlab.top = 0),
    axis.components = list(
        bottom = ac, top = ac,
        left = ac, right = ac)
)

trellis.par.set(ps0 <- lattice:::updateList(ps1, ps2), warn = FALSE)

#=======================================================================
# Lattice function for graphics

#-------------------------------------------
# Plot profile in lattice
xyprofile <- function(prof, conf = c(0.9, 0.95, 0.99),
                      namestrip = NULL, subset = 4,
                      scales.x = "free", ...) {
    #-------------------------------------------
    conf <- conf[order(conf, decreasing = TRUE)]
    da <- subset(prof, abs(z) <= subset)
    #-------------------------------------------
    fl <- levels(da$param)
    if (!is.null(namestrip)) {
        fl <- namestrip
    }
    xyplot(abs(z) ~ focal | param,
           data = da,
           layout = c(NA, 1),
           # ylab = expression(abs(z)~~(sqrt(~Delta~"deviance"))),
           ylab = expression(sqrt(2~(L(hat(phi))-L(phi)))),
           scales = list(x = scales.x),
           type = c("l", "g"),
           strip = strip.custom(
               factor.levels = fl),
           panel = function(x, y, subscripts, ...) {
               conf <- c(0.9, 0.95, 0.99)
               hl <- sqrt(qchisq(conf, 1))
               #-------------------------------------------
               mle <- x[y == 0]
               xl <- x[x < mle]; yl <- y[x < mle]
               xr <- x[x > mle]; yr <- y[x > mle]
               #-------------------------------------------
               funleft <- approxfun(x = yl, y = xl)
               funright <- approxfun(x = yr, y = xr)
               vxl <- funleft(hl)
               vxr <- funright(hl)
               vz <- c(hl, hl)
               ##-------------------------------------------
               panel.xyplot(x, y, ...)
               panel.arrows(c(vxl, vxr), 0, c(vxl, vxr), vz,
                            code = 1, length = 0.1, lty = 2,
                            col = "gray40")
               panel.segments(vxl, vz, vxr, vz, lty = 2,
                              col = "gray40")
               panel.abline(h = 0, v = mle, lty = 3)
               panel.text(x = rep(mle, 2), y = vz + 0.1,
                          labels = paste(conf*100, "%"),
                          col = "gray20")
           }, ...)
}

#-------------------------------------------
# Panel for envelope bands
panel.cbH <- function(x, y, ly, uy,
                      subscripts, cty,
                      col.line = plot.line$col,
                      lwd = plot.line$lwd,
                      desloc = NULL,
                      fill = 1, alpha = 0.1, length = 0.05, ...) {
    plot.line <- trellis.par.get("plot.line")
    if (is.null(desloc)) {
        desloc <- rep(0, length(uy))
    }
    y <- as.numeric(y)
    x <- as.numeric(x)
    or <- order(x)
    ly <- as.numeric(ly[subscripts])
    uy <- as.numeric(uy[subscripts])
    xo <- x[or]
    yo <- y[or]
    lyo <- ly[or]
    uyo <- uy[or]
    desl <- desloc[subscripts]
    if (cty == "bands") {
        panel.polygon(c(xo, rev(xo)), c(lyo, rev(uyo)), col = fill,
                      alpha = alpha, border = NA)
        panel.lines(xo, lyo, lty = list(...)$lty, lwd = 0.5,
                    col = col.line)
        panel.lines(xo, uyo, lty = list(...)$lty, lwd = 0.5,
                    col = col.line)
    }
    if (cty == "bars") {
        panel.arrows(xo + desl, lyo, xo + desl, uyo, length = length,
                     code = 3, angle = 90, col = col.line, lwd = lwd)
    }
    panel.xyplot(x + desl, y, subscripts = subscripts,
                 col.line = col.line, lwd = lwd, ...)
}

prepanel.cbH <- function(y, ly, uy, subscripts) {
    ly <- as.numeric(ly[subscripts])
    uy <- as.numeric(uy[subscripts])
    y <- as.numeric(y[subscripts])
    list(ylim = range(y, uy, ly, finite = TRUE))
}

#-------------------------------------------
# Panel for boxplot with `whisker.width` argument
my.panel.bwplot <- function(x, y, box.ratio = 1, box.width = box.ratio/(1 +
    box.ratio), horizontal = TRUE, pch = box.dot$pch, col = box.dot$col,
    alpha = box.dot$alpha, cex = box.dot$cex, font = box.dot$font,
    fontfamily = box.dot$fontfamily, fontface = box.dot$fontface,
    fill = box.rectangle$fill, varwidth = FALSE, notch = FALSE,
    notch.frac = 0.5, ..., levels.fos = if (horizontal) sort(unique(y)) else sort(unique(x)),
    stats = boxplot.stats, coef = 1.5, do.out = TRUE,
    identifier = "bwplot", whisker.width) {
    if (all(is.na(x) | is.na(y)))
        return()
    x <- as.numeric(x)
    y <- as.numeric(y)
    box.dot <- trellis.par.get("box.dot")
    box.rectangle <- trellis.par.get("box.rectangle")
    box.umbrella <- trellis.par.get("box.umbrella")
    plot.symbol <- trellis.par.get("plot.symbol")
    fontsize.points <- trellis.par.get("fontsize")$points
    if (!notch)
        notch.frac <- 0
    if (horizontal) {
        blist <- tapply(x, factor(y, levels = levels.fos), stats,
            coef = coef, do.out = do.out)
        blist.stats <- t(sapply(blist, "[[", "stats"))
        blist.out <- lapply(blist, "[[", "out")
        blist.height <- box.width
        if (varwidth) {
            maxn <- max(table(y))
            blist.n <- sapply(blist, "[[", "n")
            blist.height <- sqrt(blist.n/maxn) * blist.height
        }
        blist.conf <- if (notch)
            t(sapply(blist, "[[", "conf")) else blist.stats[, c(2, 4), drop = FALSE]
        xbnd <- cbind(blist.stats[, 3], blist.conf[, 2], blist.stats[,
            4], blist.stats[, 4], blist.conf[, 2], blist.stats[,
            3], blist.conf[, 1], blist.stats[, 2], blist.stats[,
            2], blist.conf[, 1], blist.stats[, 3])
        ytop <- levels.fos + blist.height/2
        ybot <- levels.fos - blist.height/2
        ybnd <- cbind(ytop - notch.frac * blist.height/2, ytop,
            ytop, ybot, ybot, ybot + notch.frac * blist.height/2,
            ybot, ybot, ytop, ytop, ytop - notch.frac * blist.height/2)
        xs <- cbind(xbnd, NA_real_)
        ys <- cbind(ybnd, NA_real_)
        panel.polygon(t(xs), t(ys), lwd = box.rectangle$lwd,
            lty = box.rectangle$lty, col = fill, alpha = box.rectangle$alpha,
            border = box.rectangle$col, identifier = paste(identifier,
                "box", sep = "."))
        panel.segments(
            c(blist.stats[, 2], blist.stats[, 4]),
            rep(levels.fos, 2),
            c(blist.stats[, 1], blist.stats[, 5]),
            rep(levels.fos, 2),
            col = box.umbrella$col,
            alpha = box.umbrella$alpha, lwd = box.umbrella$lwd,
            lty = box.umbrella$lty, identifier = paste(identifier,
                "whisker", sep = "."))
        panel.segments(
            c(blist.stats[, 1], blist.stats[, 5]),
            levels.fos - blist.height * whisker.width/2,
            c(blist.stats[, 1], blist.stats[, 5]),
            levels.fos + blist.height * whisker.width/2,
            col = box.umbrella$col, alpha = box.umbrella$alpha,
            lwd = box.umbrella$lwd, lty = box.umbrella$lty,
            identifier = paste(identifier,
                "cap", sep = "."))
        if (all(pch == "|")) {
            mult <- if (notch)
                1 - notch.frac else 1
            panel.segments(blist.stats[, 3], levels.fos - mult *
                blist.height/2, blist.stats[, 3], levels.fos +
                mult * blist.height/2, lwd = box.rectangle$lwd,
                lty = box.rectangle$lty, col = box.rectangle$col,
                alpha = alpha, identifier = paste(identifier,
                  "dot", sep = "."))
        } else {
            panel.points(x = blist.stats[, 3], y = levels.fos,
                pch = pch, col = col, alpha = alpha, cex = cex,
                fontfamily = fontfamily, fontface = lattice:::chooseFace(fontface,
                  font), fontsize = fontsize.points, identifier = paste(identifier,
                  "dot", sep = "."))
        }
        panel.points(x = unlist(blist.out), y = rep(levels.fos,
            sapply(blist.out, length)), pch = plot.symbol$pch,
            col = plot.symbol$col, alpha = plot.symbol$alpha,
            cex = plot.symbol$cex, fontfamily = plot.symbol$fontfamily,
            fontface = lattice:::chooseFace(plot.symbol$fontface, plot.symbol$font),
            fontsize = fontsize.points, identifier = paste(identifier,
                "outlier", sep = "."))
    } else {
        blist <- tapply(y, factor(x, levels = levels.fos), stats,
            coef = coef, do.out = do.out)
        blist.stats <- t(sapply(blist, "[[", "stats"))
        blist.out <- lapply(blist, "[[", "out")
        blist.height <- box.width
        if (varwidth) {
            maxn <- max(table(x))
            blist.n <- sapply(blist, "[[", "n")
            blist.height <- sqrt(blist.n/maxn) * blist.height
        }
        blist.conf <- if (notch)
            sapply(blist, "[[", "conf") else t(blist.stats[, c(2, 4), drop = FALSE])
        ybnd <- cbind(blist.stats[, 3], blist.conf[2, ], blist.stats[,
            4], blist.stats[, 4], blist.conf[2, ], blist.stats[,
            3], blist.conf[1, ], blist.stats[, 2], blist.stats[,
            2], blist.conf[1, ], blist.stats[, 3])
        xleft <- levels.fos - blist.height/2
        xright <- levels.fos + blist.height/2
        xbnd <- cbind(xleft + notch.frac * blist.height/2, xleft,
            xleft, xright, xright, xright - notch.frac * blist.height/2,
            xright, xright, xleft, xleft, xleft + notch.frac *
                blist.height/2)
        xs <- cbind(xbnd, NA_real_)
        ys <- cbind(ybnd, NA_real_)
        panel.polygon(t(xs), t(ys), lwd = box.rectangle$lwd,
            lty = box.rectangle$lty, col = fill, alpha = box.rectangle$alpha,
            border = box.rectangle$col, identifier = paste(identifier,
                "box", sep = "."))
        panel.segments(rep(levels.fos, 2), c(blist.stats[, 2],
            blist.stats[, 4]), rep(levels.fos, 2), c(blist.stats[,
            1], blist.stats[, 5]), col = box.umbrella$col, alpha = box.umbrella$alpha,
            lwd = box.umbrella$lwd, lty = box.umbrella$lty, identifier = paste(identifier,
                "whisker", sep = "."))
        panel.segments(levels.fos - blist.height/2, c(blist.stats[,
            1], blist.stats[, 5]), levels.fos + blist.height/2,
            c(blist.stats[, 1], blist.stats[, 5]), col = box.umbrella$col,
            alpha = box.umbrella$alpha, lwd = box.umbrella$lwd,
            lty = box.umbrella$lty, identifier = paste(identifier,
                "cap", sep = "."))
        if (all(pch == "|")) {
            mult <- if (notch)
                1 - notch.frac else 1
            panel.segments(levels.fos - mult * blist.height/2,
                blist.stats[, 3], levels.fos + mult * blist.height/2,
                blist.stats[, 3], lwd = box.rectangle$lwd, lty = box.rectangle$lty,
                col = box.rectangle$col, alpha = alpha, identifier = paste(identifier,
                  "dot", sep = "."))
        } else {
            panel.points(x = levels.fos, y = blist.stats[, 3],
                pch = pch, col = col, alpha = alpha, cex = cex,
                fontfamily = fontfamily, fontface = lattice:::chooseFace(fontface,
                  font), fontsize = fontsize.points, identifier = paste(identifier,
                  "dot", sep = "."))
        }
        panel.points(x = rep(levels.fos, sapply(blist.out, length)),
            y = unlist(blist.out), pch = plot.symbol$pch, col = plot.symbol$col,
            alpha = plot.symbol$alpha, cex = plot.symbol$cex,
            fontfamily = plot.symbol$fontfamily, fontface = lattice:::chooseFace(plot.symbol$fontface,
                plot.symbol$font), fontsize = fontsize.points,
            identifier = paste(identifier, "outlier", sep = "."))
    }
}

#-------------------------------------------
# Panel for segplots divide by groups
centfac <- function (group, space = NULL) {
    stopifnot(is.factor(group))
    if (is.null(space)) {
        space <- 0.5/nlevels(group)
    }
    d <- 2 * ((as.integer(group) - 1)/(nlevels(group) - 1)) -
        1
    return(space * d)
}

panel.groups.segplot <- function (x, y, z, centers, groups, gap = NULL,
                                  data, subscripts,  ...) {
    if (!missing(data)) {
        data <- eval(data, envir = parent.frame())
        groups <- data[, deparse(substitute(groups))]
    }
    stopifnot(is.factor(groups))
    stopifnot(length(groups) == length(z))
    if (is.null(gap)) {
        gap <- 0.5/nlevels(groups)
    }
    d <- 2 * ((as.numeric(groups) - 1)/(nlevels(groups) - 1)) -
        1
    z <- as.numeric(z) + gap * d
    panel.segplot(x, y, z, centers = centers,
                  subscripts = subscripts, ...)
}

panel.beeswarm <- function (x, y, subscripts, spread, ...) {
    xx <- x
    yy <- y
    aux <- by(cbind(yy, xx, subscripts),
              INDICES = xx,
              FUN = function(i) {
                  or <- order(i[, 1])
                  ys <- i[or, 1]
                  yt <- table(ys)
                  dv <- sapply(unlist(yt), FUN = function(j) {
                      seq(from = 1, to = j, length.out = j) - (j + 1)/2
                  })
        if (!is.list(dv)) {
            dv <- as.list(dv)
        }
        xs <- i[or, 2] + spread * do.call(c, dv)
        cbind(x = xs, y = ys, subscripts = subscripts[or])
    })
    aux <- do.call(rbind, aux)
    panel.xyplot(aux[, 1], aux[, 2], subscripts = aux[, 3], ...)
}
