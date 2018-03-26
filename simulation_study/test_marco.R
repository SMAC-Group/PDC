simu1_part1 <- run_simulation(nb_simu = 10,
                              true_beta = c(1, rep(0,58), 1),
                              n = 80, n_star = 800, sigma = 1,
                              cor_X = 0.5,
                              adaptive = FALSE)

make_main_graph(simu1_part1)

plot(simu1_part1, xlim = c(1.035, 1.09))
     ,
     ylim = c(2.35, 5.2), xlab = expression(PE[y]),
     ylab = expression(MSE[beta]))
plot(simu1_part1, add_correct = FALSE, add_selected = FALSE)
simu1_part1$summary_table
