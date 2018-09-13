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
phis <- c(-1.6, -1, 0, 1.8)
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
              as.table = TRUE,
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

# Save into compress file
out <- list("models_orpar" = models_orpar,
            "models_repar" = models_repar,
            "devs_orpar" = devs_orpar,
            "devs_repar" = devs_repar)
saveRDS(out, "orthogonality.rds", compress = FALSE)
