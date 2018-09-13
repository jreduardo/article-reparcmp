
#=======================================================================
# General settings
#=======================================================================

## ----setup
library(knitr)
library(xtable)
options(digits = 3, OutDec = ".",
        xtable.caption.placement = "top",
        xtable.booktabs = TRUE,
        xtable.sanitize.text.function = identity,
        xtable.size = "small",
        xtable.math.style.negative = TRUE)
opts_chunk$set(
    warning = FALSE,
    message = FALSE,
    echo = FALSE,
    results = "hide",
    fig.width = 7,
    fig.height = 5,
    out.width = "1\\textwidth",
    fig.align = "center",
    fig.pos = "!htp",
    # dev = "tikz"
    dev.args = list(family = "Palatino")
    )

divide_matrix <- function(mat, divide) {
    mat1 <- mat[, 1:divide]
    mat2 <- mat[, (divide + 1):ncol(mat)]
    nc1 <- ncol(mat1)
    nc2 <- ncol(mat2)
    cdiff <- abs(ncol(mat1) - ncol(mat2))
    if (cdiff) {
        complete <- matrix(NA, nrow = nrow(mat), ncol = cdiff)
        colnames(complete) <- rep(" ", cdiff)
        if (nc1 < nc2) {
            mat1 <- cbind(mat1, complete)
        } else {
            mat2 <- cbind(mat2, complete)
        }
    }
    return(list(mat1, mat2))
}

source("./codes/lattice-panels.R")
source("./codes/functions.R")

# Colors for legends
cols <- trellis.par.get("superpose.line")$col

# Useful packages
library(bbmle)
library(methods)
library(multcomp)
library(plyr)
library(tidyr)
library(dplyr)

#=======================================================================
# Convergence of Z
#=======================================================================

## ----convergenceZ
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

## ----results-Z
caption <- paste(
    "Values for $Z(\\lambda, \\nu)$ constant (numerically computed)",
    "for values of $\\lambda$ (0.5 to 50) and $\\phi$ (0 to 1)")

# Remove excessive scientific notation
xt2 <- gsub("e\\+00", "    ", format(xt))
xt2 <- gsub("Inf", "divergent$^{**}$", format(xt2))
xt2[1, -1] <- "divergent$^{*\\,\\,\\,}$"

print(xtable(xt2, digits = -2,
             caption = caption,
             align = "C|CCCCCC",
             # align = "c|cccccc",
             label = "tab:convergenceZ"),
      include.colnames = FALSE,
      tabular.environment = "tabularx",
      only.contents = TRUE,
      width = "\\textwidth",
      add.to.row = list(
          pos = list(0, 0),
          command = c(
              "& \\multicolumn{6}{c}{$\\bm{\\lambda}$} \\\\\n",
              "$\\bm{\\nu}$ & 0.5 & 1 & 5 & 10 & 30 & 50 \\\\\n")
      ))

#=======================================================================
# Study the approximation
#=======================================================================

## ----data-approx
# Mean and variance relationship
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

maemu <- max(abs(grid$emu))
maeva <- max(abs(grid$eva))

## ----approx-plot
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

# COM-Poisson probabilities
## ----pmf-cmp
parg <- expand.grid(mu = c(2, 8, 15), phi = c(-0.7, 0, 0.7))
y <- 0:30
py <- mapply(FUN = dcmp,
             mu = parg$mu,
             phi = parg$phi,
             MoreArgs = list(y = y, sumto = 100),
             SIMPLIFY = FALSE)
parg <- cbind(parg[rep(1:nrow(parg), each = length(y)), ],
              y = y, py = unlist(py))

leg_phi <- parse(
    text = paste("phi == \"",
                 formatC(unique(parg$phi), 1, format = "f"),
                 "\""))
barchart(py ~ y | factor(mu),
         groups = factor(phi),
         data = parg,
         horizontal = FALSE,
         layout = c(NA, 1),
         as.table = TRUE,
         axis = axis.grid,
         origin = 0,
         xlim = extendrange(y, f = 0.01),
         border = "transparent",
         # scales = list(x = list(at = pretty(y))),
         ylab = expression(P(Y==y)),
         xlab = expression(y),
         par.settings = list(
             superpose.polygon = list(
                 col = c("gray40", "gray60", "gray10")),
             superpose.line = list(
                 col = c("gray40", "gray60", "gray10"),
                 lwd = 2)
         ),
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

#-------------------------------------------
# Indexes

## ----compute-indexes
# Dispersion index
grid <- transform(grid, dispersion_index = var / mean)

# Zero-inflation index
prob0 <- mapply(FUN = dcmp,
                mu = grid$mu,
                phi = grid$phi,
                MoreArgs = list(y = 0, sumto = 300),
                SIMPLIFY = FALSE)
prob0 <- unlist(prob0)
grid <- transform(grid, zero_index = 1 + log(prob0) / mean)

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

# Study indexes on the reparametrized COM-Poisson distribution
## ----indexes-plot
mygpar <- gpar(cex = 1.2)
psaux <- modifyList(
    ps0, list(superpose.line = list(
                  col = myreg(length(unique(grid$phi)))),
              layout.widths = list(
                  left.padding = -0.5)
              )
)

# Mean and variance relationship
xy1 <- xyplot(var ~ mu | "Mean-Variance",
              groups = phi,
              data = grid,
              type = "l",
              lwd = 2,
              axis = axis.grid,
              xlab = expression(mu),
              ylab = "",
              # ylab = expression(Var(Y)),
              scales = list(y = list(rot = 90)),
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
                            gp = mygpar,
                            x = unit(0.020, "npc"),
                            y = unit(0.875, "npc"))
              })

# Dispersion index (DI)
xy2 <- xyplot(var / mean ~ mu | "Dispersion index",
              groups = phi,
              data = grid,
              type = "l",
              lwd = 2,
              axis = axis.grid,
              xlab = expression(mu),
              ylab = "",
              scales = list(y = list(rot = 90)),
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
                            gp = mygpar,
                            x = unit(0.020, "npc"),
                            y = unit(0.875, "npc"))
              })

# Zero-inflation index
xy3 <- xyplot(zero_index ~ mu | "Zero-inflation index",
              groups = phi,
              data = grid,
              type = "l",
              lwd = 2,
              axis = axis.grid,
              xlab = expression(mu),
              ylab = "",
              scales = list(y = list(rot = 90)),
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
                            gp = mygpar,
                            x = unit(0.020, "npc"),
                            y = unit(0.875, "npc"))
              })

# Heavy tail-index
# Choose one specific mu of unique(grid$mu)
choosemu <- 25
xy4 <- xyplot(heavy_index ~ x | "Heavy-tail index",
              groups = phi,
              # data = heavy_data,
              data = subset(heavy_data, mu == choosemu),
              type = "l",
              lwd = 2,
              axis = axis.grid,
              xlab = expression(y),
              ylab = "",
              scales = list(y = list(rot = 90)),
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
                            gp = mygpar,
                            x = unit(0.020, "npc"),
                            y = unit(0.875, "npc"))
              })

print(xy1, split = c(1, 1, 4, 1), more = TRUE)
print(xy2, split = c(2, 1, 4, 1), more = TRUE)
print(xy3, split = c(3, 1, 4, 1), more = TRUE)
print(xy4, split = c(4, 1, 4, 1), more = FALSE)

midi <- min(grid$dispersion_index)
madi <- max(grid$dispersion_index)

#=======================================================================
# Simulation study
#=======================================================================

# Configuration
## ----load-simulation
B <- 1000
beta <- c("b0" = 2, "b1" = 0.5, "b21" = 0.8, "b22" = -0.8)
phis <- c(0, -1.6, -1, 1.8)
names(phis) <- sprintf("phi=%s", phis)

sizes <- c(50, 100, 300, 1000)
names(sizes) <- sprintf("n=%s", sizes)

# Load results
results <- readRDS("./codes/simulation.rds")

## ----justpars

# Justify the choices of the parameters simulation
x1 <- seq(0, 1, length.out = 100)
x2 <- rep(letters[1:3], length.out = 100)
X <- model.matrix(~ x1 + x2)
mu <- exp(X %*% beta)

daexplain <- ldply(lapply(phis, function(phi) {
    va <- compute_variance(mu, phi, sumto = 300)
    data.frame(mu = mu, va = va, di = va / mu, x1 = x1, x2 = x2)
}), .id = "phi")

labx2 <- expression("Level of categorical variable"~(x[2]))
labx1 <- expression("Values of continous variable"~(x[1]))

fl <- parse(text = gsub("=", "==", levels(daexplain$phi)))
xy1 <- xyplot(mu ~ x1, groups = x2,
              type = c("g", "l"),
              lwd = 2,
              xlab = labx1,
              ylab = "Average count",
              auto.key = list(
                  column = 3,
                  points = FALSE,
                  lines = TRUE,
                  title = labx2,
                  cex.title = 1.1
              ))
xy2 <- xyplot(di ~ x1 | phi, groups = x2,
              type = c("g", "l"),
              lwd = 2,
              # scales = "free",
              layout = c(2, 2),
              xlab = labx1,
              ylab = "Dispersion index",
              data = daexplain,
              strip = strip.custom(factor.levels = fl))

print(xy1, split = c(1, 1, 2, 1), more = TRUE)
print(xy2, split = c(2, 1, 2, 1), more = FALSE)

# Compute standardized bias
## ----bias-data
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
        # (bhat - real)                 # raw bias
        (bhat - real) / std[[1]][[y]] # standardized bias
    })
    ldply(out, .id = "phi")
}), .id = "n")

# Organize the results
aux <- na.omit(aux)
bias <- gather(aux, param, bias, phi2:b22, factor_key = TRUE)
bias$phi <- ordered(bias$phi, c("phi=-1.6", "phi=-1", "phi=0",
                                "phi=1.8"))

# Distributions of bias and average with confidence intervals
## ----bias-plot
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

# Compute coverage rate
## ----coverage-data
aux <- ldply(lapply(results, function(x) {
    ind <- names(x); names(ind) <- ind
    out <- lapply(ind, function(y) {
        ind <- lapply(x[[y]], function(z) {
            real <- c(phis[[y]], beta)
            cint <- z[["cint"]]
            ind <- as.integer(cint[, 1] < real & cint[, 2] > real)
            names(ind) <- c("phi2", names(beta))
            ind
        })
        do.call("rbind", ind)
    })
    ldply(out, .id = "phi")
}), .id = "n")

# Organize the results
aux <- na.omit(aux)
coverage <- gather(aux, param, coverage, phi2:b22, factor_key = TRUE)
coverage$phi <- ordered(coverage$phi, c("phi=-1.6", "phi=-1", "phi=0",
                                        "phi=1.8"))
coverage$n <- as.numeric(gsub("n=([0-9]+)", "\\1", coverage$n))
covdata <- aggregate(coverage ~ n + phi + param, data = coverage,
                     function(x) sum(x == 1) / length(x))

## ----coverage-plot
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
           layout.heights = list(strip = 1.5)
       ),
       key = list(
           column = 4,
           type = "o",
           divide = 1,
           text = list(fl),
           lines = list(pch = 19, col = cols[1:4])
       ),
       strip = strip.custom(
           factor.levels = yl
       ),
       panel = function(x, y, ...) {
           panel.xyplot(x, y, ...)
           panel.abline(h = 0.95, lty = 2)
       })

# Obtain the covariance between dispersion and regression parameters
## ----ortho-plot
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
bwplot(param ~ value | phi,
       groups = n,
       data = subset(dacov, scale == "cor"),
       box.width = 0.12,
       pch = "|",
       grid = TRUE,
       as.table = TRUE,
       layout = c(4, 1),
       xlab = "Empirical correlations with dispersion parameter estimator",
       fill = colorRampPalette(c("gray80",  "gray20"))(4),
       panel = panel.superpose,
       scales = list(y = list(labels = yl), x = "free"),
       whisker.width = 0.4,
       strip = strip.custom(factor.levels = fl),
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
       })

## ----ortho-surf

# Plots deviance surfaces
output <- readRDS("codes/orthogonality.rds")
devs_orpar <- output$devs_orpar
devs_repar <- output$devs_repar

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
        colorkey = FALSE,
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
    strip = strip.custom(factor.levels = fl),
    strip.left = strip.custom(par.strip.text=list(cex = 0.85))
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
        colorkey = FALSE,
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
    strip = strip.custom(factor.levels = fl),
    strip.left = strip.custom(par.strip.text=list(cex = 0.85))
)

# Organize plots
print(xy1, position = c(0.0, 0.48, 1.0, 1.0), more = TRUE)
print(xy2, position = c(0.0, 0.0, 1.0, 0.52), more = FALSE)

## ----fit-cotton

# Load data
# data(cottonBolls, package = "cmpreg")
cottonBolls <- read.table("./data/cottonBolls.txt",
                          header = TRUE, sep = "\t")
cottonBolls$est <- ordered(
    cottonBolls$est,
    c("vegetative", "flower bud", "blossom", "boll", "boll open")
)

# Fit models
mnames <- c("PO", "C1", "C2", "QP")

# Predictor, following Zeviani et al. (2014)
form1 <- ncap ~ est:(des + I(des^2))

m1PO <- glm(form1, data = cottonBolls, family = poisson)
time11 <- system.time(
    m1C1 <- fitcm(form1, data = cottonBolls, model = "CP", sumto = 50)
)
time12 <- system.time(
    m1C2 <- fitcm(form1, data = cottonBolls, model = "CP2", sumto = 50)
)
m1QP <- glm(form1, data = cottonBolls, family = quasipoisson)

models.ncap <- list(m1PO, m1C1, m1C2, m1QP)
names(models.ncap) <- mnames

# Numbers of calls to loglik and numerical gradient
c11 <- models.ncap$C1@details$counts
c12 <- models.ncap$C2@details$counts

# LRT between Poisson and COM-Poisson (test: phi == 0)
lrt.ncap <- getAnova(m1PO, m1C2)

# Goodness of fit measures and estimate parameters
## ----results-cotton

# GoF measures
measures.ncap <- lapply(models.ncap, function(x)
    c("LogLik" = logLik(x), "AIC" = AIC(x), "BIC" = BIC(x)))

# Get the estimates
co1 <- coef(m1C2)
est <- lapply(models.ncap, FUN = function(x) getCoefs(x))
est.ncap <- do.call(cbind, est)

# Organize in table
pnames <- c("\\phi\\,,\\,\\sigma", "\\beta_0",
            paste0("\\beta_{1", 1:5, "}"),
            paste0("\\beta_{2", 1:5, "}"))
rownames(est.ncap) <- paste0("$", pnames, "$")
meds <- apply(do.call(cbind, measures.ncap), 1, function(x) {
    x <- formatC(x, getOption("digits"), format = "f")
    x <- gsub("NA", "-", x)
    paste(paste0("\\multicolumn{2}{c}{$", x, "$}"),
          collapse = " & ")
})
text_gof <- paste(paste(names(meds), "&", meds),
                  collapse = "\\\\\n ")

append_gof <- list(
    pos = list(nrow(est.ncap)),
    command = paste("\\specialrule{0.01em}{0.3em}{0.3em} \n",
                    text_gof, "\\\\\n",
                    "\\bottomrule"))
print.xtable(xtable(est.ncap, digits = 4,
                    label = "tab:coef-cotton"),
             hline.after = 0,
             only.contents = TRUE,
             add.to.row = append_gof)

## ----pred-cotton

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
key <- list(columns = 4,
            cex = 0.9,
            lines = list(col = 1, lty = rev(1:4)),
            text = list(parse(
                text = c("'Poisson'", "'COM-Poisson'",
                         "'COM-Poisson'[mu]", "'Quasi-Poisson'"))
                ))

# Graph
xyplot(ncap ~ des | est,
       data = cottonBolls,
       layout = c(NA, 1),
       as.table = TRUE,
       grid = TRUE,
       type = "p",
       xlab = "Artificial defoliation level",
       ylab = "Number of bolls produced",
       spread = 0.05,
       key = key,
       alpha = 0.6,
       panel = panel.beeswarm) +
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

# Correlation between estimates
## ----corr-cotton
corr.ncap <- do.call("rbind",
                     lapply(models.ncap[c("C1", "C2")],
                            function(x) cov2cor(vcov(x))[1, -1]))

# Organize on table
rownames(corr.ncap) <- paste0("COM-Poisson", c("", "$_\\mu$"))
colnames(corr.ncap) <- gsub("beta", "hat{\\\\beta}", pnames[-1])
corrncap_aux <- divide_matrix(corr.ncap, divide = 6)

## ----cres-cotton
print(xtable(corrncap_aux[[1]], digits = 4),
      only.contents = TRUE,
      hline.after = FALSE,
      sanitize.colnames.function = function(x) sprintf("$%s$", x))
cat("\\midrule", sep = "\n")

print(xtable(corrncap_aux[[2]], digits = 4),
      only.contents = TRUE,
      sanitize.colnames.function = function(x) sprintf("$%s$", x))

## ----fit-soy

# Load data
# data(soyaBeans, package = "cmpreg")
soyBeans <- read.table("./data/soyaBeans.txt",
                        header = TRUE, sep = "\t")
soyBeans$umid <- as.factor(soyBeans$umid)
soyBeans <- soyBeans[-74, ] # Incorrect observation
soyBeans <- transform(soyBeans, K = K / 100)

# Fit models

form2 <-  ngra ~ bloc + umid * K + I(K^2)

m2PO <- glm(form2, data = soyBeans, family = poisson)
time21 <- system.time(
    m2C1 <- fitcm(form2, data = soyBeans, model = "CP", sumto = 700)
)
time22 <- system.time(
    m2C2 <- fitcm(form2, data = soyBeans, model = "CP2", sumto = 700)
)
m2QP <- glm(form2, data = soyBeans, family = quasipoisson)

models.ngra <- list(m2PO, m2C1, m2C2, m2QP)
names(models.ngra) <- mnames

co2 <- coef(m2C2)

# Numbers of calls to loglik and numerical gradient
c21 <- models.ngra$C1@details$counts
c22 <- models.ngra$C2@details$counts

# # Profile extra parameter
# profs.ngra <- lapply(list(c(m2C1, "phi"), c(m2C2, "phi2")),
#                      function(x) myprofile(x[[1]], x[[2]]))
# profs.ngra <- do.call("rbind", profs.ngra)

# LRT between Poisson and COM-Poisson (test: phi == 0)
lrt.ngra <- getAnova(m2PO, m2C2)

## ----desc-soy

# Exploratory analysis

# Scatter plot
xy1 <- xyplot(ngra ~ K | umid,
              data = soyBeans,
              xlab = "Potassium fertilization level",
              ylab = "Number of grains per pot",
              type = c("p", "g", "smooth"),
              as.table =  TRUE,
              layout = c(2, 2),
              strip = strip.custom(
                  strip.names = TRUE, var.name = "moisture",
                  factor.levels = paste0(levels(soyBeans$umid), "%")))

# Sample variance vs sample mean (evidence in favor of the
# overdispersion).
mv <- soyBeans %>%
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

# Goodness of fit measures and estimate parameters
## ----results-soya

# GoF measures
measures.ngra <- lapply(models.ngra, function(x)
    c("LogLik" = logLik(x), "AIC" = AIC(x), "BIC" = BIC(x)))

# Get the estimates
co2 <- coef(m2C2)
est <- lapply(models.ngra, FUN = function(x) getCoefs(x))
est.ngra <- do.call(cbind, est)

# Organize in table
pnames <- c("\\phi\\,,\\,\\sigma", "\\beta_0",
            paste0("\\gamma_{", 1:4, "}"),
            paste0("\\tau_{", 1:2, "}"),
            "\\beta_1", "\\beta_2",
            paste0("\\beta_{3", 1:2, "}"))

rownames(est.ngra) <- paste0("$", pnames, "$")
meds <- apply(do.call(cbind, measures.ngra), 1, function(x) {
    x <- formatC(x, getOption("digits"), format = "f")
    x <- gsub("NA", "-", x)
    paste(paste0("\\multicolumn{2}{c}{$", x, "$}"),
          collapse = " & ")
})
text_gof <- paste(paste(names(meds), "&", meds),
                  collapse = "\\\\\n ")

append_gof <- list(
    pos = list(nrow(est.ngra)),
    command = paste("\\specialrule{0.01em}{0.3em}{0.3em} \n",
                    text_gof, "\\\\\n",
                    "\\bottomrule"))
print.xtable(xtable(est.ngra, digits = 4,
                    label = "tab:coef-soy"),
             hline.after = 0,
             only.contents = TRUE,
             add.to.row = append_gof)

## Correlation between estimates
## ----corr-soy

corr.ngra <- do.call("rbind",
                     lapply(models.ngra[c("C1", "C2")],
                            function(x) cov2cor(vcov(x))[1, -1]))

## Organize on table
rownames(corr.ngra) <- paste0("COM-Poisson", c("", "$_\\mu$"))
colnames(corr.ngra) <- gsub("(beta|tau|gamma)", "hat{\\\\\\1}",
                            pnames[-1])
corrngra_aux <- divide_matrix(corr.ncap, divide = 7)

## ----cres-soy

print(xtable(corrngra_aux[[1]], digits = 4),
      only.contents = TRUE,
      hline.after = FALSE,
      sanitize.colnames.function = function(x) sprintf("$%s$", x))
cat("\\midrule", sep = "\n")

print(xtable(corrngra_aux[[2]], digits = 4),
      only.contents = TRUE,
      sanitize.colnames.function = function(x) sprintf("$%s$", x))

## ----pred-soy

# Prediction

# Data for prediction
pred <- with(soyBeans,
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
key <- list(columns = 4,
            cex = 0.9,
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


## ----fit-ovos

# Load data
# data(nitrofen, package = "boot")
# data(Paula, package = "labestData")
# nitrofen <- PaulaEx4.6.20
nitrofen <- read.table("./data/nitrofen.txt",
                       header = TRUE, sep = "\t")
nitrofen <- transform(nitrofen, dose = dose / 100)

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

## ----anova-ovos1
auxPO <- lapply(fmodels.ovos, function(x) x$PO)
tab <- do.call("getAnova", c(print = FALSE, auxPO))
tab <- cbind(tab, NA)
rownames(tab) <- paste("Preditor", rownames(tab))
digits <- c(1, 0, 3, 3, 3, 0, -2, 3)
print(xtable(tab, digits = digits),
      include.colnames = FALSE,
      hline.after = NULL,
      only.contents = TRUE)


## ----anova-ovos2
auxC1 <- lapply(fmodels.ovos, function(x) x$C1)
tab <- do.call("getAnova", c(print = FALSE, auxC1))
tab <- cbind(tab, sapply(auxC1, function(x) coef(x)[1]))
rownames(tab) <- paste("Preditor", rownames(tab))
print(xtable(tab, digits = digits),
      include.colnames = FALSE,
      hline.after = NULL,
      only.contents = TRUE)


## ----anova-ovos3
auxC2 <- lapply(fmodels.ovos, function(x) x$C2)
tab <- do.call("getAnova", c(print = FALSE, auxC2))
tab <- cbind(tab, sapply(auxC2, function(x) coef(x)[1]))
rownames(tab) <- paste("Preditor", rownames(tab))
print(xtable(tab, digits = digits),
      include.colnames = FALSE,
      hline.after = NULL,
      only.contents = TRUE)

## ----anova-ovos4
auxQP <- lapply(fmodels.ovos, function(x) x$QP)
tab <- do.call("getAnova", c(print = FALSE, auxQP))
tab <- cbind(tab, sapply(auxQP, function(x) summary(x)$dispersion))
rownames(tab) <- paste("Preditor", rownames(tab))
print(xtable(tab, digits = digits),
      include.colnames = FALSE,
      hline.after = NULL,
      only.contents = TRUE)

# Goodness of fit measures and estimate parameters
## ----coef-ovos

# GoF measures
measures.ovos <- lapply(models.ovos, function(x)
    c("LogLik" = logLik(x), "AIC" = AIC(x), "BIC" = BIC(x)))

# Get the estimates
co2 <- coef(m2C2)
est <- lapply(models.ovos, FUN = function(x) getCoefs(x))
est.ovos <- do.call(cbind, est)[-1, ]

# Organize in table
pnames <- paste0("\\beta_{", 0:3, "}")
rownames(est.ovos) <- paste0("$", pnames, "$")
print.xtable(xtable(est.ovos, digits = 4),
             hline.after = c(0, nrow(est.ovos)),
             only.contents = TRUE)

# Prediction
## ----pred-ovos

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
xyplot(novos ~ dose,
       data = nitrofen,
       xlab = "Nitrofen concentration level",
       ylab = "Number of live offspring",
       grid = TRUE,
       alpha = 0.6,
       key = key,
       spread = 0.05,
       panel = panel.beeswarm) +
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

# Correlation between estimates
## ----corr-ovos

corr.ovos <- do.call("rbind",
                     lapply(models.ovos[c("C1", "C2")],
                            function(x) cov2cor(vcov(x))[1, -1]))

# Organize on table
rownames(corr.ovos) <- paste0("COM-Poisson", c("", "$_\\mu$"))
colnames(corr.ovos) <- gsub("(beta)", "hat{\\\\\\1}", pnames)

caption <- paste("Empirical correlations between $\\hat{\\phi}$ and",
                 "$\\hat{\\bm{\\beta}}$ for the two parametrizations",
                 "of COM-Poisson model fit to equidispersed data.")
print(xtable(corr.ovos,
             align = c("lrrrr"),
             caption = caption,
             digits = 4,
             label = "tab:corr-ovos"),
      # size = "small",
      sanitize.rownames.function = identity,
      sanitize.colnames.function = function(x) sprintf("$%s$", x),
      # width = "\\textwidth,"
      tabular.environment = "tabular")
