#-----------------------------------------------------------------------
# Load package and functions

## ---- setup
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
    fig.pos = "H",
    # dev = "tikz"
    dev.args = list(family = "Palatino")
    )

## ---- load-packages
source("../codes/lattice-panels.R")
source("../codes/functions.R")
library(purrr)
library(dplyr)
library(tidyr)
cols <- trellis.par.get("superpose.line")$col

## ---- counts-ranging
rcmp <- Vectorize(FUN = function(n, mu, phi) {
    par <- mu2lambda(mu, phi)
    with(par, compoisson::rcom(n = n, lambda, nu))
}, vectorize.args = "mu")


phis <- c(-1.6, -1, 0, 1.8)
names(phis) <- sprintf("phi=%s", phis)
fl <- parse(text = gsub("=", "==", names(phis)))
da <- map_dfr(phis, function(p) {
    tibble(phi = factor(p, levels = phis),
           y   = c(rcmp(1000, 27, p)))
})

histogram(~y | phi,
          type="density",
          layout = c(4, 1),
          axis = axis.grid,
          breaks = seq(0, max(da$y) + 1, 1),
          strip = strip.custom(factor.levels = fl),
          panel = function(x, ...) {
              panel.histogram(x, ...)
          },
          data = da)


## ----approx-var
compute_apvar <- function(mu, phi) {
    mu * exp(-phi) + (exp(phi) - 1) / (2 * exp(2 * phi))
}

mu <- seq(1, 30, length.out = 500)
mus <- c(2, 8, 20)
names(mus) <- sprintf("mu=%s", mus)

phi <- seq(-3.5, 5, length.out = 500)
phis <- c(-2, 0, 2)
names(phis) <- sprintf("phi=%s", phis)

fa <- parse(text = gsub("=", "==", names(mus)))
xy1 <-
    map_dfr(mus, function(m) {
    tibble(mu = factor(m, levels = mus),
           va = compute_apvar(m, phi = phi),
           phi = phi)}) %>%
    xyplot(va ~ phi,
           groups = mu,
           type = "l",
           xlab = expression(phi),
           ylab = expression(aVar(Y[i])),
           sub = "(a)",
           lwd = 3,
           grid = TRUE,
           key = list(
               corner = c(0.9, 0.9),
               lines = list(col = cols[1:3], lwd = 2),
               text = list(fa)
           ),
           panel = function(...) {
               panel.xyplot(...)
               panel.abline(h = 0, lty = 2)
           },
           data = .)

fb <- parse(text = gsub("=", "==", names(phis)))
xy2 <-
    map_dfr(phis, function(p) {
        tibble(phi = factor(p, levels = phis),
               va = compute_apvar(mu, phi = p),
               mu = mu)}) %>%
    xyplot(va ~ mu,
           groups = phi,
           type = "l",
           xlab = expression(mu),
           ylab = expression(aVar(Y[i])),
           sub = "(b)",
           lwd = 3,
           grid = TRUE,
           key = list(
               corner = c(0.1, 0.9),
               lines = list(col = cols[1:3], lwd = 2),
               text = list(fb)
           ),
           panel = function(...) {
               panel.xyplot(...)
               panel.abline(h = 0, lty = 2)
           },
           data = .)

plot(xy1, split = c(1, 1, 2, 1), more = TRUE)
plot(xy2, split = c(2, 1, 2, 1), more = FALSE)

## ---- explore-constrain
compute_mu <- function(lambda, nu) lambda^(1/nu) - (nu-1)/(2*nu)
grid <- crossing(lambda = seq(0.001, 0.1, length.out = 100),
                 nu = seq(0.01, 7, length.out = 100)) %>%
    mutate(mu = apply(., 1L,
                      function(x) compute_mu(x[1L], x[2L])),
           rmu = apply(., 1L,
                       function(x) compoisson::com.mean(x[1L], x[2L])),
           id = mu > 0)

xy1 <-
    levelplot(id ~ nu * lambda,
              xlab = expression(nu),
              ylab = expression(lambda),
              sub = "(a)",
              col.regions = c("gray90", "gray50"),
              colorkey = FALSE,
              key = list(
                  space = "top",
                  columns = 2,
                  rectangle = list(col = c("gray90", "gray50")),
                  text = list(expression(mu<0, mu>0))
              ),
              panel = function(...) {
                  panel.levelplot(...)
                  panel.curve(((x-1)/(2*x))^x, lwd = 1.2)
                  panel.abline(v = 1, lty = 2, lwd = 1.2)
              },
              data = grid)
xy2 <-
    levelplot(rmu ~ nu * lambda,
              xlab = expression(nu),
              ylab = expression(lambda),
              sub = "(b)",
              colorkey = list(space = "top"),
              panel = function(...) {
                  panel.levelplot(...)
                  panel.curve(((x-1)/(2*x))^x, lwd = 1.2)
                  panel.abline(v = 1, lty = 2, lwd = 1.2)
              },
              data = grid)

print(xy1, position = c(0, 0, 0.5, 0.94), more = TRUE)
print(xy2, position = c(0.5, 0, 1, 1), more = FALSE)


## ---- contraints-case-studies

cottonBolls <- read.table("../data/cottonBolls.txt",
                          header = TRUE, sep = "\t")
form1 <- ncap ~ est:(des + I(des^2))
m1C2 <- fitcm(form1, data = cottonBolls, model = "CP2", sumto = 50)

soyaBeans <- read.table("../data/soyaBeans.txt",
                        header = TRUE, sep = "\t")[-74, ] %>%
    mutate(umid = as.factor(umid), K = K / 100)
form2 <-  ngra ~ bloc + umid * K + I(K^2)
m2C2 <- fitcm(form2, data = soyaBeans, model = "CP2", sumto = 700)

nitrofen <- read.table("../data/nitrofen.txt",
                       header = TRUE, sep = "\t")
nitrofen <- transform(nitrofen, dose = dose / 100)
form3 <-  novos ~ dose + I(dose^2) + I(dose^3)
m3C2 <- fitcm(form3, data = nitrofen, model = "CP2", sumto = 100)

com = c(paste("\\toprule\n",
              "& & \\multicolumn{3}{c}{$\\hat{\\lambda}_i$}\\\\\n",
              "\\cmidrule(lr){3-5}\n"),
        "\\bottomrule\n")
list(m1C2, m2C2, m3C2) %>%
    set_names(paste(c("Cotton", "Soybean", "Nitrofen"),
                    "experiment")) %>%
    map_dfr(.id = "case", .f = function(m) {
        eta <- unique(c(m@data$X %*% m@coef[-1]))
        tibble(mu = exp(eta), phi = m@coef[1])
    }) %>%
    group_by(case, phi) %>%
    summarise_at(vars(mu), funs(min, median, max)) %>%
    ungroup() %>%
    set_names("Case study", "$\\hat{\\phi}$",
              "Minimun", "Median", "Maximum") %>%
    xtable() %>%
    print.xtable(
        include.rownames = FALSE,
        hline.after = FALSE,
        add.to.row = list(
            as.list(c(-1, nrow(.))),
            com)
    )

## ---- computation-times
benchs <- readRDS("../codes/comparetimes.rds") %>%
    mutate(case = forcats::fct_inorder(case))
xlabs <- expression("COMPo"(lambda[i], phi),
                    "COMPo"(mu[i], phi))
bwplot(time/10e6 ~ expr | case,
       ylab = "Time (seconds)",
       scales = list(
           y = list(relation = "free"),
           x = list(at = 1:2, labels = xlabs)
       ),
       data = benchs)

## ---- no-include
