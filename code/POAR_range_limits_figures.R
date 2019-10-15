source('all_vr_bayes.R')


# Posterior predictive checks ---------------------------------------------


ppc_dens_overlay(data_all$y_s, y_s_sim)
ppc_dens_overlay(data_all$y_g, y_g_sim)+xlim(0, 100)
ppc_dens_overlay(data_all$y_f, y_f_sim)
ppc_dens_overlay(data_all$y_p, y_p_sim)+xlim(0, 100)
ppc_dens_overlay(data_all$y_v, y_v_sim) ## maybe need beta-binomial?
ppc_dens_overlay(data_all$y_m, y_m_sim) ## maybe need beta-binomial?


# Core vital rates --------------------------------------------------------
# estimate sex- and longitude-specific vital rates for small, medium, and large plants


# Seed viability ----------------------------------------------------------
viab_pars <- rstan::summary(fit_full, pars=c("v0","a_v"))[[1]][,"mean"]

plot(data_all$SR_v,(data_all$y_v / data_all$tot_seeds_v),
     cex = 3 * (data_all$tot_seeds_v / max(data_all$tot_seeds_v)))
lines(seq(0,1,0.01),
      viab_pars["v0"] * (1 - seq(0,1,0.01) ^ viab_pars["a_v"]))