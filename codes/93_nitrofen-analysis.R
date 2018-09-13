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
