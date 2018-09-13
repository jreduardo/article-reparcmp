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

# Times to the first case study
bench1 <- microbenchmark(
    "CMP  " = fitcm(form1, data = data1, model = "CP" , sumto = 50),
    "CMPmu" = fitcm(form1, data = data1, model = "CP2", sumto = 50),
    times = 50)

#-------------------------------------------
# Unit: seconds
#   expr      min       lq     mean   median       uq      max neval cld
#  CMP   2.134262 2.144675 2.167662 2.157547 2.187935 2.326526    50   b
#  CMPmu 1.533308 1.549595 1.568905 1.559458 1.581462 1.627994    50  a
#-------------------------------------------

# Times to the second case study
bench2 <- microbenchmark(
    "CMP  " = fitcm(form2, data = data2, model = "CP" , sumto = 700),
    "CMPmu" = fitcm(form2, data = data2, model = "CP2", sumto = 700),
    times = 50)

#-------------------------------------------
# Unit: seconds
#   expr       min        lq      mean    median        uq       max neval cld
#  CMP   14.240183 14.331817 14.413927 14.392338 14.436526 15.027403    50   b
#  CMPmu  6.778149  6.840705  6.878897  6.856285  6.911408  7.047864    50  a
#-------------------------------------------

# Times to the second case study
bench3 <- microbenchmark(
    "CMP  " = fitcm(form3, data = data3, model = "CP" , sumto = 100),
    "CMPmu" = fitcm(form3, data = data3, model = "CP2", sumto = 100),
    times = 50)

#-------------------------------------------
# Unit: seconds
#   expr      min       lq     mean   median       uq      max neval cld
#  CMP   1.125297 1.143321 1.189389 1.160909 1.237385 1.357786    50   b
#  CMPmu 1.034631 1.041504 1.091136 1.071659 1.128976 1.258067    50  a
#-------------------------------------------

#-----------------------------------------------------------------------
# Organize results
benchs <- list("Cotton (under)" = bench1,
               "Soybean (over)" = bench2,
               "Nitrofen (equi)" = bench3)
benchs <- purrr::map_dfr(benchs, identity, .id = "case")

class(benchs) <- "data.frame"
saveRDS(benchs, "comparetimes.rds")

xlabs <- expression("COM-Poisson", "COM-Poisson"[mu])
bwplot(time ~ expr | case,
       ylab = "Time (seconds)",
       scales = list(
           y = list(relation = "free"),
           x = list(at = 1:2, labels = xlabs)
       ),
       data = benchs)
