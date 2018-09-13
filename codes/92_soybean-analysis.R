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
