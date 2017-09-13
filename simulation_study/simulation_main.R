library(tikzDevice)
library(glmnet)
library(lars)
library(elasticnet)
library(ncvreg)
library(parcor)
library(MASS)
library(mvtnorm)
library(knitr)
library(ellipse)
library(pdc)
library(knitr)

# Simulation 1:
tikz('simulation_study/figs/simu1_part1.tex', width=6,height=5)
simu1_part1 <- run_simulation(nb_simu = 500, true_beta = c(1, rep(0,58), 1),
                              n = 80, n_star = 800, sigma = 1, cor_X = 0.5,
                              xlim = c(1.035, 1.095), ylim = c(0.02, 0.1),
                              xlab = "PE$_{y}$",
                              ylab = "MSE$_{\\beta}$",
                              title = "Simulation 1: $n < p$")

dev.off()
save.image("data/simu1_part1.RData")
rm(list = ls())


# Simulation 2:
tikz('simulation_study/figs/simu2_part1.tex', width=6,height=5)
simu1_part2 <- run_simulation(nb_simu = 500, true_beta = rep(c(0.3,0),5),
                              n = 100, n_star = 1000, sigma = 1, cor_X = 0.75,
                              xlim = c(1.065,1.17), ylim = c(0.18,0.46),
                              xlab = "PE$_{y}$",
                              ylab = "MSE$_{\\beta}$",
                              title = "Simulation 2: $n < p$")

dev.off()
save.image("data/simu2_part1.RData")
rm(list = ls())


# Simulation 3:
tikz('simulation_study/figs/simu3_part1.tex', width=6,height=5)
simu1_part2 <- run_simulation(nb_simu = 500,
                              true_beta = c(rep(c(2,0,1),2), rep(0,16), rep(0.1,6),
                                            rep(0,16), rep(c(2,0,1), 2)),
                              n = 100, n_star = 1000, sigma = 2, cor_X = 0.5,
                              xlim = c(4.73,5.4), ylim = c(0.79,1.8),
                              xlab = "PE$_{y}$",
                              ylab = "MSE$_{\\beta}$",
                              title = "Simulation 3: $n < p$")

dev.off()
save.image("data/simu3_part1.RData")
rm(list = ls())


