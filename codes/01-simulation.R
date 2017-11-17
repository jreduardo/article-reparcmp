#=======================================================================
# Simulation study about estimators
#                                                        Eduardo Junior
#                                                    edujrrib@gmail.com
#                                                            2017-08-04
#=======================================================================

##----------------------------------------------------------------------
## Load package and functions
source("lattice-panels.R")
source("functions.R")
library(plyr)
library(tidyr)
library(bbmle)

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
# Random generation for the COM-Poisson distribution
rcmp <- Vectorize(FUN = function(n, mu, phi) {
    par <- mu2lambda(mu, phi)
    with(par, compoisson::rcom(n = n, lambda, nu))
}, vectorize.args = "mu")

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
            fit <- try(fitcm(y ~ x1 + x2, start = c(phi, beta),
                             model = "CP2", sumto = 200))
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
            out <- list("coef" = co, "cint" = ci,
                        "loglik" = ll, "vcov" = vc)
            return(out)
        }) # Close replicate
    }) # Close phi
}, mc.cores = 3) # Close sizes
str(results)

saveRDS(results, "simulation.rds", compress = FALSE)

#-------------------------------------------
# Load the results
results <- readRDS("simulation.rds")

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

covdata <- aggregate(coverage ~ n + phi + param, data = coverage,
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
