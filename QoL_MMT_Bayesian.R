library(tidybayes)
library(tidyverse)
library(magrittr)
library(haven)
library(brms)
library(bayesplot)
library(cowplot)
library(gmodels)
library(psych)
library(table1)
#===============================
df <- readRDS(file = "QoL_MMT_Vietnam.rds")

#----- Descriptive analysis
#===============================================================================
# Set up template for tables
render.cat <- function(x) {
    c("", sapply(stats.default(x), 
                 function(y) with(y, sprintf("%d (%0.1f)", FREQ, PCT))))
}

#----- Table 1: Participant’s characteristics
table1(~ gender_f + edu_f + marital + occu  + dose_gr + sideEff01_f + urine_f + 
           hiv_f + hbv_f + hcv_f  + age + duration + dose, data = df, 
       overall = "Total", render.categorical = render.cat)

#----- Table 2: Quality of life and family functioning among MMT patients
table1(~ c1 + c2 + c3 + c4 + c5 + apgar + phys + psyco + social + envir | apgar_cat, 
       data = df, overall = "Total", render.categorical = render.cat)


#----- Figure 1: Correlation between family functioning and four aspects of quality of life among MMT patients

df %>% select(phys, psyco, social, envir, apgar) %>% gather(-apgar, key = "scale", value = "value") %>% 
    mutate(scale = recode(scale, "phys" = "Physical health",
                               "psyco" = "Psychological",
                               "social" = "Social relationships",
                               "envir" = "Environment")) %>%
    ggplot(aes(x = apgar, y = value)) +
    geom_jitter(size = 2, alpha = 0.5, show.legend = F, color = "gray10") + 
    #scale_color_viridis(option = "C", direction = -1) +
    expand_limits(ylim = c(0,100)) +
    geom_smooth(method = "lm", fill = "grey50", color = "black") +
    facet_wrap(~scale, nrow = 2) + 
    labs(x = "Family functioning score",
         y = "Quality of life scores") + 
    theme_bw() +
    theme(
        axis.title.y = element_text(size = 14, color = "black",  face = "bold"),
        axis.title.x = element_text(size = 14, color = "black",  face = "bold"),
        axis.text.y = element_text(size = 13, color = "black"),
        axis.text.x = element_text(size = 13, color = "black"),
        strip.text = element_text(size = 13, face = "bold"),
        strip.background  = element_rect(color = NULL, fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    ) ->> f1

f1

# save figure to *tiff format
tiff("Figure1.tiff", units="in", width= 8, height=7, res=500, compression = "lzw")
f1
dev.off()

# Calculate correlation coefficients
df %>% select(phys, psyco, social, envir, apgar) %>% psych::pairs.panels(stars = T)

#----- Bayesian model
#===============================================================================
# Set 8 cores parallel on CPU
options(mc.cores = 8)

# Define a loss function (rmse, mae) for Posterior predictive checks
rmse <- function(y, yrep) {
    yrep_mean <- colMeans(yrep)
    sqrt(mean((yrep_mean - y)^2))
}

mae <- function(y, yrep) {
    yrep_mean <- colMeans(yrep)
    mean(abs(yrep_mean - y))
}


# Model Details

# HRQoLi ~ N(μi, σ2)
# μi ~ α + b1*APGARi + b2*X2i + … + bk*Xki , i = 1, …, n.
# Vague prior: b1 ~ N(0, 100^2)
# Equivocal prior: b1 ~ N(0, 1.02^2)
# Optimistic prior: b1 ~ N(2, 1.02^2)
# Pessimistic prior: b1 ~ N(-2, 1.02^2)

#--- set up prior for 4 Scenarios
# brms uses SD to present normal distribution (i.e., N(mu, sigma))

prior_vague = c(set_prior(prior = "normal(0, 100)", class = "b", coef = "apgar"))
prior_equal = c(set_prior(prior = "normal(0, 1.020408)", class = "b", coef = "apgar"))
prior_opti = c(set_prior(prior = "normal(2, 1.020408)", class = "b", coef = "apgar"))
prior_pessi = c(set_prior(prior = "normal(-2, 1.020408)", class = "b", coef = "apgar"))

#====== Scale 1: Physical health
#===============================
# Get non-informative prior for remainded coefficient
get_prior(phys ~ apgar + age + edu + marital + gender,
          data = df)

# fit models
phys_prior1 = brm(phys ~ apgar + age + edu + marital + gender,
                data = df, 
                seed = 123,
                iter = 8000,
                warmup = 1000,
                chains = 4,
                thin = 2,
                prior = prior_vague)

phys_prior2 = brm(phys ~ apgar + age + edu + marital + gender,
                  data = df, 
                  seed = 123,
                  iter = 8000,
                  warmup = 1000,
                  chains = 4,
                  thin = 2,
                  prior = prior_equal)

phys_prior3 = brm(phys ~ apgar + age + edu + marital + gender,
                  data = df, 
                  seed = 123,
                  iter = 8000,
                  warmup = 1000,
                  chains = 4,
                  thin = 2,
                  prior = prior_opti)

phys_prior4 = brm(phys ~ apgar + age + edu + marital + gender,
                  data = df, 
                  seed = 123,
                  iter = 8000,
                  warmup = 1000,
                  chains = 4,
                  thin = 2,
                  prior = prior_pessi)

# Graphical posterior predictive checks
color_scheme_set("viridisB")

# Set theme for figures
pp_theme <- function(...) {
    theme(
        plot.subtitle = element_text(face = "italic", size = 9, hjust = 0.5),
        plot.title = element_text(face = "bold", hjust = 0.5)
    )
}

ppc_phys1 <- pp_check(phys_prior1, type='stat', stat='mean') +
    labs(title = "A - Vague Prior",
         x = "Score") + pp_theme()
ppc_phys2 <- pp_check(phys_prior2, type='stat', stat='mean') +
    labs(title = "B - Equivocal Prior ",
         x = "Score") + pp_theme()
ppc_phys3 <- pp_check(phys_prior3, type='stat', stat='mean') +
    labs(title = "C - Optimistic Prior",
         x = "Score") + pp_theme()
ppc_phys4 <- pp_check(phys_prior4, type='stat', stat='mean') +
    labs(title = "D - Pessimistic Prior",
         x = "Score") + pp_theme()


tiff("Supf5.tiff", units="in", width= 12, height=8, res=600, compression = "lzw")
gridExtra::grid.arrange(ppc_phys1, ppc_phys3, ppc_phys2, ppc_phys4)
dev.off()

# Calculate loss
for (i in list(phys_prior1, phys_prior2, phys_prior3, phys_prior4)) {
    # perform k-fold cross validation
    kf <- kfold(i, save_fits = TRUE, chains = 1)
    # predict responses and evaluate the loss (rmse, mae)
    kfp <- kfold_predict(kf)
    rmse_val <- rmse(y = kfp$y, yrep = kfp$yrep)
    mae_val <- mae(y = kfp$y, yrep = kfp$yrep)
    print(rmse_val)
    print(mae_val)
}


# Summary results
plot(phys_prior1, combo = c("dens", "trace"), theme = theme_bw())
summary(phys_prior1)
phys_prior1 %>% gather_draws(b_apgar) %>% median_hdi()

plot(phys_prior2, combo = c("dens", "trace"), theme = theme_bw())
summary(phys_prior2)
phys_prior2 %>% gather_draws(b_apgar) %>% median_hdi()

plot(phys_prior3, combo = c("dens", "trace"), theme = theme_bw())
summary(phys_prior3)
phys_prior3 %>% gather_draws(b_apgar) %>% median_hdi()

plot(phys_prior4, combo = c("dens", "trace"), theme = theme_bw())
summary(phys_prior4)
phys_prior4 %>% gather_draws(b_apgar) %>% median_hdi()

# Posterior 
# Hypothesis testing
# Scenario 1: Vague prior
hypothesis(phys_prior1, "apgar > 0") %>% .$hypothesis
hypothesis(phys_prior1, "apgar > 1") %>% .$hypothesis
hypothesis(phys_prior1, "apgar > 2") %>% .$hypothesis

# Scenario 2: Equivocal prior
hypothesis(phys_prior2, "apgar > 0") %>% .$hypothesis
hypothesis(phys_prior2, "apgar > 1") %>% .$hypothesis
hypothesis(phys_prior2, "apgar > 2") %>% .$hypothesis

# Scenario 3: Optimistic prior
hypothesis(phys_prior3, "apgar > 0") %>% .$hypothesis
hypothesis(phys_prior3, "apgar > 1") %>% .$hypothesis
hypothesis(phys_prior3, "apgar > 2") %>% .$hypothesis

# Scenario 4: Pessimistic prior
hypothesis(phys_prior4, "apgar > 0") %>% .$hypothesis
hypothesis(phys_prior4, "apgar > 1") %>% .$hypothesis
hypothesis(phys_prior4, "apgar > 2") %>% .$hypothesis

# Save MCMC for visualization (Sup figures)
phys_MCMC1 = phys_prior1 %>% gather_draws(b_apgar)
phys_MCMC2 = phys_prior2 %>% gather_draws(b_apgar)
phys_MCMC3 = phys_prior3 %>% gather_draws(b_apgar) 
phys_MCMC4 = phys_prior4 %>% gather_draws(b_apgar) 

#===============================
#====== Scale 2: Psychological
#===============================
# Get non-informative prior for remainded coefficient
get_prior(psyco ~ apgar + age + edu + marital + gender,
          data = df)

# fit model
psyco_prior1 = brm(psyco ~ apgar + age + edu + marital + gender,
                  data = df, 
                  seed = 123,
                  iter = 8000,
                  warmup = 1000,
                  chains = 4,
                  thin = 2,
                  prior = prior_vague)

psyco_prior2 = brm(psyco ~ apgar + age + edu + marital + gender,
                  data = df, 
                  seed = 123,
                  iter = 8000,
                  warmup = 1000,
                  chains = 4,
                  thin = 2,
                  prior = prior_equal)

psyco_prior3 = brm(psyco ~ apgar + age + edu + marital + gender,
                  data = df, 
                  seed = 123,
                  iter = 8000,
                  warmup = 1000,
                  chains = 4,
                  thin = 2,
                  prior = prior_opti)

psyco_prior4 = brm(psyco ~ apgar + age + edu + marital + gender,
                  data = df, 
                  seed = 123,
                  iter = 8000,
                  warmup = 1000,
                  chains = 4,
                  thin = 2,
                  prior = prior_pessi)

# Graphical posterior predictive checks
color_scheme_set("brightblue")
ppc_psyco1 <- pp_check(psyco_prior1, type='stat', stat='mean') +
    labs(title = "A - Vague Prior",
         x = "Score") + pp_theme()
ppc_psyco2 <- pp_check(psyco_prior2, type='stat', stat='mean') +
    labs(title = "B - Equivocal Prior ",
         x = "Score") + pp_theme()
ppc_psyco3 <- pp_check(psyco_prior3, type='stat', stat='mean') +
    labs(title = "C - Optimistic Prior",
         x = "Score") + pp_theme()
ppc_psyco4 <- pp_check(psyco_prior4, type='stat', stat='mean') +
    labs(title = "D - Pessimistic Prior",
         x = "Score") + pp_theme()


tiff("Supf6.tiff", units="in", width= 12, height=8, res=600, compression = "lzw")
gridExtra::grid.arrange(ppc_psyco, ppc_psyco3, ppc_psyco2, ppc_psyco4)
dev.off()

# Calculate loss
for (i in list(psyco_prior1, psyco_prior2, psyco_prior3, psyco_prior4)) {
    # perform k-fold cross validation
    kf <- kfold(i, save_fits = TRUE, chains = 1)
    # predict responses and evaluate the loss (rmse, mae)
    kfp <- kfold_predict(kf)
    rmse_val <- rmse(y = kfp$y, yrep = kfp$yrep)
    mae_val <- mae(y = kfp$y, yrep = kfp$yrep)
    print(rmse_val)
    print(mae_val)
}

# summary data
plot(psyco_prior1, combo = c("dens", "trace"), theme = theme_bw())
summary(psyco_prior1)
psyco_prior1 %>% gather_draws(b_apgar) %>% median_hdi()

plot(psyco_prior2, combo = c("dens", "trace"), theme = theme_bw())
summary(psyco_prior2)
psyco_prior2 %>% gather_draws(b_apgar) %>% median_hdi()

plot(psyco_prior3, combo = c("dens", "trace"), theme = theme_bw())
summary(psyco_prior3)
psyco_prior3 %>% gather_draws(b_apgar) %>% median_hdi()

plot(psyco_prior4, combo = c("dens", "trace"), theme = theme_bw())
summary(psyco_prior4)
psyco_prior4 %>% gather_draws(b_apgar) %>% median_hdi()

# Posterior 
# Hypothesis testing
# Scenario 1: Vague prior
hypothesis(psyco_prior1, "apgar > 0") %>% .$hypothesis
hypothesis(psyco_prior1, "apgar > 1") %>% .$hypothesis
hypothesis(psyco_prior1, "apgar > 2") %>% .$hypothesis

# Scenario 2: Equivocal prior
hypothesis(psyco_prior2, "apgar > 0") %>% .$hypothesis
hypothesis(psyco_prior2, "apgar > 1") %>% .$hypothesis
hypothesis(psyco_prior2, "apgar > 2") %>% .$hypothesis

# Scenario 3: Optimistic prior
hypothesis(psyco_prior3, "apgar > 0") %>% .$hypothesis
hypothesis(psyco_prior3, "apgar > 1") %>% .$hypothesis
hypothesis(psyco_prior3, "apgar > 2") %>% .$hypothesis

# Scenario 4: Pessimistic prior
hypothesis(psyco_prior4, "apgar > 0") %>% .$hypothesis
hypothesis(psyco_prior4, "apgar > 1") %>% .$hypothesis
hypothesis(psyco_prior4, "apgar > 2") %>% .$hypothesis


# Save MCMC for visualization (Sup figures)
psyco_MCMC1 = psyco_prior1 %>% gather_draws(b_apgar)
psyco_MCMC2 = psyco_prior2 %>% gather_draws(b_apgar)
psyco_MCMC3 = psyco_prior3 %>% gather_draws(b_apgar) 
psyco_MCMC4 = psyco_prior4 %>% gather_draws(b_apgar) 

#===================================
#====== Scale 3: Social relationships
#===================================
# Get non-informative prior for remainded coefficient
get_prior(social ~ apgar + age + edu + marital + gender,
          data = df)

# fit model
social_prior1 = brm(social ~ apgar + age + edu + marital + gender,
                   data = df, 
                   seed = 123,
                   iter = 8000,
                   warmup = 1000,
                   chains = 4,
                   thin = 2,
                   prior = prior_vague)

social_prior2 = brm(social ~ apgar + age + edu + marital + gender,
                   data = df, 
                   seed = 123,
                   iter = 8000,
                   warmup = 1000,
                   chains = 4,
                   thin = 2,
                   prior = prior_equal)

social_prior3 = brm(social ~ apgar + age + edu + marital + gender,
                   data = df, 
                   seed = 123,
                   iter = 8000,
                   warmup = 1000,
                   chains = 4,
                   thin = 2,
                   prior = prior_opti)

social_prior4 = brm(social ~ apgar + age + edu + marital + gender,
                   data = df, 
                   seed = 123,
                   iter = 8000,
                   warmup = 1000,
                   chains = 4,
                   thin = 2,
                   prior = prior_pessi)

# Graphical posterior predictive checks
color_scheme_set("mix-purple-yellow")
ppc_social1 <- pp_check(social_prior1, type='stat', stat='mean') +
    labs(title = "A - Vague Prior",
         x = "Score") + pp_theme()
ppc_social2 <- pp_check(social_prior2, type='stat', stat='mean') +
    labs(title = "B - Equivocal Prior ",
         x = "Score") + pp_theme()
ppc_social3 <- pp_check(social_prior3, type='stat', stat='mean') +
    labs(title = "C - Optimistic Prior",
         x = "Score") + pp_theme()
ppc_social4 <- pp_check(social_prior4, type='stat', stat='mean') +
    labs(title = "D - Pessimistic Prior",
         x = "Score") + pp_theme()


tiff("Supf7.tiff", units="in", width= 12, height=8, res=600, compression = "lzw")
gridExtra::grid.arrange(ppc_social1, ppc_social3, ppc_social2, ppc_social4)
dev.off()

# Calculate loss
for (i in list(social_prior1, social_prior2, social_prior3, social_prior4)) {
    # perform k-fold cross validation
    kf <- kfold(i, save_fits = TRUE, chains = 1)
    # predict responses and evaluate the loss (rmse, mae)
    kfp <- kfold_predict(kf)
    rmse_val <- rmse(y = kfp$y, yrep = kfp$yrep)
    mae_val <- mae(y = kfp$y, yrep = kfp$yrep)
    print(rmse_val)
    print(mae_val)
}

# Summary data
plot(social_prior1, combo = c("dens", "trace"), theme = theme_bw())
summary(social_prior1)
social_prior1 %>% gather_draws(b_apgar) %>% median_hdi()

plot(social_prior2, combo = c("dens", "trace"), theme = theme_bw())
summary(social_prior2)
social_prior2 %>% gather_draws(b_apgar) %>% median_hdi()

plot(social_prior3, combo = c("dens", "trace"), theme = theme_bw())
summary(social_prior3)
social_prior3 %>% gather_draws(b_apgar) %>% median_hdi()

plot(social_prior4, combo = c("dens", "trace"), theme = theme_bw())
summary(social_prior4)
social_prior4 %>% gather_draws(b_apgar) %>% median_hdi()

# Posterior 
# Hypothesis testing
# Scenario 1: Vague prior
hypothesis(social_prior1, "apgar > 0") %>% .$hypothesis
hypothesis(social_prior1, "apgar > 1") %>% .$hypothesis
hypothesis(social_prior1, "apgar > 2") %>% .$hypothesis

# Scenario 2: Equivocal prior
hypothesis(social_prior2, "apgar > 0") %>% .$hypothesis
hypothesis(social_prior2, "apgar > 1") %>% .$hypothesis
hypothesis(social_prior2, "apgar > 2") %>% .$hypothesis

# Scenario 3: Optimistic prior
hypothesis(social_prior3, "apgar > 0") %>% .$hypothesis
hypothesis(social_prior3, "apgar > 1") %>% .$hypothesis
hypothesis(social_prior3, "apgar > 2") %>% .$hypothesis

# Scenario 4: Pessimistic prior
hypothesis(social_prior4, "apgar > 0") %>% .$hypothesis
hypothesis(social_prior4, "apgar > 1") %>% .$hypothesis
hypothesis(social_prior4, "apgar > 2") %>% .$hypothesis

# Save MCMC for visualization (Sup figures)
social_MCMC1 = social_prior1 %>% gather_draws(b_apgar)
social_MCMC2 = social_prior2 %>% gather_draws(b_apgar)
social_MCMC3 = social_prior3 %>% gather_draws(b_apgar) 
social_MCMC4 = social_prior4 %>% gather_draws(b_apgar) 


#===========================
#====== Scale 4: Environment
#===========================
# Get non-informative prior for remainded coefficient
get_prior(envir ~ apgar + age + edu + marital + gender,
          data = df)

# fit model
envir_prior1 = brm(envir ~ apgar + age + edu + marital + gender,
                   data = df, 
                   seed = 123,
                   iter = 8000,
                   warmup = 1000,
                   chains = 4,
                   thin = 2,
                   prior = prior_vague)

envir_prior2 = brm(envir ~ apgar + age + edu + marital + gender,
                   data = df, 
                   seed = 123,
                   iter = 8000,
                   warmup = 1000,
                   chains = 4,
                   thin = 2,
                   prior = prior_equal)

envir_prior3 = brm(envir ~ apgar + age + edu + marital + gender,
                   data = df, 
                   seed = 123,
                   iter = 8000,
                   warmup = 1000,
                   chains = 4,
                   thin = 2,
                   prior = prior_opti)

envir_prior4 = brm(envir ~ apgar + age + edu + marital + gender,
                   data = df, 
                   seed = 123,
                   iter = 8000,
                   warmup = 1000,
                   chains = 4,
                   thin = 2,
                   prior = prior_pessi)

# Graphical posterior predictive checks
color_scheme_set("teal")
ppc_envir1 <- pp_check(envir_prior1, type='stat', stat='mean') +
    labs(title = "A - Vague Prior",
         x = "Score") + pp_theme()
ppc_envir2 <- pp_check(envir_prior2, type='stat', stat='mean') +
    labs(title = "B - Equivocal Prior ",
         x = "Score") + pp_theme()
ppc_envir3 <- pp_check(envir_prior3, type='stat', stat='mean') +
    labs(title = "C - Optimistic Prior",
         x = "Score") + pp_theme()
ppc_envir4 <- pp_check(envir_prior4, type='stat', stat='mean') +
    labs(title = "D - Pessimistic Prior",
         x = "Score") + pp_theme()


tiff("Supf8.tiff", units="in", width= 12, height=8, res=600, compression = "lzw")
gridExtra::grid.arrange(ppc_envir1, ppc_envir3, ppc_envir2, ppc_envir4)
dev.off()

# Calculate loss
for (i in list(envir_prior1, envir_prior2, envir_prior3, envir_prior4)) {
    # perform k-fold cross validation
    kf <- kfold(i, save_fits = TRUE, chains = 1)
    # predict responses and evaluate the loss (rmse, mae)
    kfp <- kfold_predict(kf)
    rmse_val <- rmse(y = kfp$y, yrep = kfp$yrep)
    mae_val <- mae(y = kfp$y, yrep = kfp$yrep)
    print(rmse_val)
    print(mae_val)
}

# Summary data

plot(envir_prior1, combo = c("dens", "trace"), theme = theme_bw())
summary(envir_prior1)
envir_prior1 %>% gather_draws(b_apgar) %>% median_hdi()

plot(envir_prior2, combo = c("dens", "trace"), theme = theme_bw())
summary(envir_prior2)
envir_prior2 %>% gather_draws(b_apgar) %>% median_hdi()

plot(envir_prior3, combo = c("dens", "trace"), theme = theme_bw())
summary(envir_prior3)
envir_prior3 %>% gather_draws(b_apgar) %>% median_hdi()

plot(envir_prior4, combo = c("dens", "trace"), theme = theme_bw())
summary(envir_prior4)
envir_prior4 %>% gather_draws(b_apgar) %>% median_hdi()

# Posterior 
# Hypothesis testing
# Scenario 1: Vague prior
hypothesis(envir_prior1, "apgar > 0") %>% .$hypothesis
hypothesis(envir_prior1, "apgar > 1") %>% .$hypothesis
hypothesis(envir_prior1, "apgar > 2") %>% .$hypothesis

# Scenario 2: Equivocal prior
hypothesis(envir_prior2, "apgar > 0") %>% .$hypothesis
hypothesis(envir_prior2, "apgar > 1") %>% .$hypothesis
hypothesis(envir_prior2, "apgar > 2") %>% .$hypothesis

# Scenario 3: Optimistic prior
hypothesis(envir_prior3, "apgar > 0") %>% .$hypothesis
hypothesis(envir_prior3, "apgar > 1") %>% .$hypothesis
hypothesis(envir_prior3, "apgar > 2") %>% .$hypothesis

# Scenario 4: Pessimistic prior
hypothesis(envir_prior4, "apgar > 0") %>% .$hypothesis
hypothesis(envir_prior4, "apgar > 1") %>% .$hypothesis
hypothesis(envir_prior4, "apgar > 2") %>% .$hypothesis

# Save MCMC for visualization (Sup figures)
envir_MCMC1 = envir_prior1 %>% gather_draws(b_apgar)
envir_MCMC2 = envir_prior2 %>% gather_draws(b_apgar)
envir_MCMC3 = envir_prior3 %>% gather_draws(b_apgar) 
envir_MCMC4 = envir_prior4 %>% gather_draws(b_apgar) 

#----- Data visualization for supplement figures (1, 2, 3, and 4)
# add theme for visualization
mytheme <- function(...) {
    theme(
        axis.line.y = element_blank(),
        panel.background = element_rect(fill = "white"),
        plot.subtitle = element_text(face = "italic", size = 9, hjust = 0.5),
        plot.title = element_text(face = "bold", hjust = 0.5)
    )
}
#===============================================================================
# Supplemental Fig. 1: Prior and posterior distributions for the relationship between 
# family functioning and physical health aspect
#===============================

Supf1a <- ggplot(data = data.frame(x = c(-6, 6)), aes(x)) +
    stat_function(fun = dnorm, n = 182,
                  args = list(mean = 0, sd = 100),
                  color="black",linetype=2) +
    stat_function(fun = dnorm,
                  n = 182,
                  args = list(mean = mean(phys_MCMC1$.value),
                              sd = sd(phys_MCMC1$.value)),color="black") +
    geom_area(stat = "function", fun = dnorm,args = list(mean = mean(phys_MCMC1$.value),
                                                         sd = sd(phys_MCMC1$.value)),
              fill = "#969696", xlim = c(0,6),alpha=0.5) +
    geom_area(stat = "function", fun = dnorm,args = list(mean = mean(phys_MCMC1$.value),
                                                         sd = sd(phys_MCMC1$.value)),
              fill = "#d9d9d9", xlim = c(-6,0),alpha=0.5) +
    geom_vline(xintercept = 0, linetype = 3)+
    scale_y_continuous(breaks = NULL)+
    scale_x_continuous(limits = c(-6,6),
                       labels = c(seq(-6,6, by = 1)),
                       breaks = c(seq(-6,6, by = 1))) +
    mytheme() + 
    labs(x="Coefficient", 
         y = NULL, 
         title = "A - Vague Prior",
         subtitle = expression('Prior: Coef ~ N(0, 100'^2*')'))


Supf1b <- ggplot(data = data.frame(x = c(-6, 6)), aes(x)) +
    stat_function(fun = dnorm, n = 182,
                  args = list(mean = 0, sd = 1.020408),
                  color="black",linetype=2) +
    stat_function(fun = dnorm,
                  n = 182,
                  args = list(mean = mean(phys_MCMC2$.value),
                              sd = sd(phys_MCMC2$.value)),color="black") +
    geom_area(stat = "function", fun = dnorm,args = list(mean = mean(phys_MCMC2$.value),
                                                         sd = sd(phys_MCMC2$.value)),
              fill = "#969696", xlim = c(0,6),alpha=0.5) +
    geom_area(stat = "function", fun = dnorm,args = list(mean = mean(phys_MCMC2$.value),
                                                         sd = sd(phys_MCMC2$.value)),
              fill = "#d9d9d9", xlim = c(-6,0),alpha=0.5) +
    geom_vline(xintercept = 0, linetype = 3)+
    scale_y_continuous(breaks = NULL)+
    scale_x_continuous(limits = c(-6,6),
                       labels = c(seq(-6,6, by = 1)),
                       breaks = c(seq(-6,6, by = 1))) +
    mytheme() + 
    labs(x="Coefficient", 
         y = NULL, 
         title = "B - Equivocal Prior",
         subtitle = expression('Prior: Coef ~ N(0, 1.02'^2*')'))


Supf1c <- ggplot(data = data.frame(x = c(-6, 6)), aes(x)) +
    stat_function(fun = dnorm, n = 182,
                  args = list(mean = 2, sd = 1.020408),
                  color="black",linetype=2) +
    stat_function(fun = dnorm,
                  n = 182,
                  args = list(mean = mean(phys_MCMC3$.value),
                              sd = sd(phys_MCMC3$.value)),color="black") +
    geom_area(stat = "function", fun = dnorm,args = list(mean = mean(phys_MCMC3$.value),
                                                         sd = sd(phys_MCMC3$.value)),
              fill = "#969696", xlim = c(0,6),alpha=0.5) +
    geom_area(stat = "function", fun = dnorm,args = list(mean = mean(phys_MCMC3$.value),
                                                         sd = sd(phys_MCMC3$.value)),
              fill = "#d9d9d9", xlim = c(-6,0),alpha=0.5) +
    geom_vline(xintercept = 0, linetype = 3)+
    scale_y_continuous(breaks = NULL)+
    scale_x_continuous(limits = c(-6,6),
                       labels = c(seq(-6,6, by = 1)),
                       breaks = c(seq(-6,6, by = 1))) +
    mytheme() + 
    labs(x="Coefficient", 
         y = NULL, 
         title = "C - Optimistic Prior",
         subtitle = expression('Prior: Coef ~ N(2, 1.02'^2*')'))

Supf1d <- ggplot(data = data.frame(x = c(-6, 6)), aes(x)) +
    stat_function(fun = dnorm, n = 182,
                  args = list(mean = -2, sd = 1.020408),
                  color="black",linetype=2) +
    stat_function(fun = dnorm,
                  n = 182,
                  args = list(mean = mean(phys_MCMC4$.value),
                              sd = sd(phys_MCMC4$.value)),color="black") +
    geom_area(stat = "function", fun = dnorm,args = list(mean = mean(phys_MCMC4$.value),
                                                         sd = sd(phys_MCMC4$.value)),
              fill = "#969696", xlim = c(0,6),alpha=0.5) +
    geom_area(stat = "function", fun = dnorm,args = list(mean = mean(phys_MCMC4$.value),
                                                         sd = sd(phys_MCMC4$.value)),
              fill = "#d9d9d9", xlim = c(-6,0),alpha=0.5) +
    geom_vline(xintercept = 0, linetype = 3)+
    scale_y_continuous(breaks = NULL)+
    scale_x_continuous(limits = c(-6,6),
                       labels = c(seq(-6,6, by = 1)),
                       breaks = c(seq(-6,6, by = 1))) +
    mytheme() + 
    labs(x="Coefficient", 
         y = NULL, 
         title = "D - Pessimistic Prior",
         subtitle = expression('Prior: Coef ~ N(-2, 1.02'^2*')'))

gridExtra::grid.arrange(Supf1a, Supf1c, Supf1b, Supf1d, ncol = 2)

tiff("Supf1_Phys.tiff", units="in", width= 12, height=7, res=600, compression = "lzw")
gridExtra::grid.arrange(Supf1a, Supf1c, Supf1b, Supf1d, ncol = 2)
dev.off()

#===============================
# Supplemental Fig. 2: Prior and posterior distributions for the relationship between 
#family functioning and psychological health aspect

Supf2a <- ggplot(data = data.frame(x = c(-6, 6)), aes(x)) +
    stat_function(fun = dnorm, n = 182,
                  args = list(mean = 0, sd = 100),
                  color="black",linetype=2) +
    stat_function(fun = dnorm,
                  n = 182,
                  args = list(mean = mean(psyco_MCMC1$.value),
                              sd = sd(psyco_MCMC1$.value)),color="black") +
    geom_area(stat = "function", fun = dnorm,args = list(mean = mean(psyco_MCMC1$.value),
                                                         sd = sd(psyco_MCMC1$.value)),
              fill = "#969696", xlim = c(0,6),alpha=0.5) +
    geom_area(stat = "function", fun = dnorm,args = list(mean = mean(psyco_MCMC1$.value),
                                                         sd = sd(psyco_MCMC1$.value)),
              fill = "#d9d9d9", xlim = c(-6,0),alpha=0.5) +
    geom_vline(xintercept = 0, linetype = 3)+
    scale_y_continuous(breaks = NULL)+
    scale_x_continuous(limits = c(-6,6),
                       labels = c(seq(-6,6, by = 1)),
                       breaks = c(seq(-6,6, by = 1))) +
    mytheme() + 
    labs(x="Coefficient", 
         y = NULL, 
         title = "A - Vague Prior",
         subtitle = expression('Prior: Coef ~ N(0, 100'^2*')'))


Supf2b <- ggplot(data = data.frame(x = c(-6, 6)), aes(x)) +
    stat_function(fun = dnorm, n = 182,
                  args = list(mean = 0, sd = 1.020408),
                  color="black",linetype=2) +
    stat_function(fun = dnorm,
                  n = 182,
                  args = list(mean = mean(psyco_MCMC2$.value),
                              sd = sd(psyco_MCMC2$.value)),color="black") +
    geom_area(stat = "function", fun = dnorm,args = list(mean = mean(psyco_MCMC2$.value),
                                                         sd = sd(psyco_MCMC2$.value)),
              fill = "#969696", xlim = c(0,6),alpha=0.5) +
    geom_area(stat = "function", fun = dnorm,args = list(mean = mean(psyco_MCMC2$.value),
                                                         sd = sd(psyco_MCMC2$.value)),
              fill = "#d9d9d9", xlim = c(-6,0),alpha=0.5) +
    geom_vline(xintercept = 0, linetype = 3)+
    scale_y_continuous(breaks = NULL)+
    scale_x_continuous(limits = c(-6,6),
                       labels = c(seq(-6,6, by = 1)),
                       breaks = c(seq(-6,6, by = 1))) +
    mytheme() + 
    labs(x="Coefficient", 
         y = NULL, 
         title = "B - Equivocal Prior",
         subtitle = expression('Prior: Coef ~ N(0, 1.02'^2*')'))


Supf2c <- ggplot(data = data.frame(x = c(-6, 6)), aes(x)) +
    stat_function(fun = dnorm, n = 182,
                  args = list(mean = 2, sd = 1.020408),
                  color="black",linetype=2) +
    stat_function(fun = dnorm,
                  n = 182,
                  args = list(mean = mean(psyco_MCMC3$.value),
                              sd = sd(psyco_MCMC3$.value)),color="black") +
    geom_area(stat = "function", fun = dnorm,args = list(mean = mean(psyco_MCMC3$.value),
                                                         sd = sd(psyco_MCMC3$.value)),
              fill = "#969696", xlim = c(0,6),alpha=0.5) +
    geom_area(stat = "function", fun = dnorm,args = list(mean = mean(psyco_MCMC3$.value),
                                                         sd = sd(psyco_MCMC3$.value)),
              fill = "#d9d9d9", xlim = c(-6,0),alpha=0.5) +
    geom_vline(xintercept = 0, linetype = 3)+
    scale_y_continuous(breaks = NULL)+
    scale_x_continuous(limits = c(-6,6),
                       labels = c(seq(-6,6, by = 1)),
                       breaks = c(seq(-6,6, by = 1))) +
    mytheme() + 
    labs(x="Coefficient", 
         y = NULL, 
         title = "C - Optimistic Prior",
         subtitle = expression('Prior: Coef ~ N(2, 1.02'^2*')'))

Supf2d <- ggplot(data = data.frame(x = c(-6, 6)), aes(x)) +
    stat_function(fun = dnorm, n = 182,
                  args = list(mean = -2, sd = 1.020408),
                  color="black",linetype=2) +
    stat_function(fun = dnorm,
                  n = 182,
                  args = list(mean = mean(psyco_MCMC4$.value),
                              sd = sd(psyco_MCMC4$.value)),color="black") +
    geom_area(stat = "function", fun = dnorm,args = list(mean = mean(psyco_MCMC4$.value),
                                                         sd = sd(psyco_MCMC4$.value)),
              fill = "#969696", xlim = c(0,6),alpha=0.5) +
    geom_area(stat = "function", fun = dnorm,args = list(mean = mean(psyco_MCMC4$.value),
                                                         sd = sd(psyco_MCMC4$.value)),
              fill = "#d9d9d9", xlim = c(-6,0),alpha=0.5) +
    geom_vline(xintercept = 0, linetype = 3)+
    scale_y_continuous(breaks = NULL)+
    scale_x_continuous(limits = c(-6,6),
                       labels = c(seq(-6,6, by = 1)),
                       breaks = c(seq(-6,6, by = 1))) +
    mytheme() + 
    labs(x="Coefficient", 
         y = NULL, 
         title = "D - Pessimistic Prior",
         subtitle = expression('Prior: Coef ~ N(-2, 1.02'^2*')'))

gridExtra::grid.arrange(Supf2a, Supf2c, Supf2b, Supf2d, ncol = 2)

tiff("Supf2_Psyco.tiff", units="in", width= 12, height=7, res=600, compression = "lzw")
gridExtra::grid.arrange(Supf2a, Supf2c, Supf2b, Supf2d, ncol = 2)
dev.off()

#===============================
# Supplemental Fig. 3: Prior and posterior distributions for the relationship between 
# family functioning and social relationships aspect

Supf3a <- ggplot(data = data.frame(x = c(-6, 6)), aes(x)) +
    stat_function(fun = dnorm, n = 182,
                  args = list(mean = 0, sd = 100),
                  color="black",linetype=2) +
    stat_function(fun = dnorm,
                  n = 182,
                  args = list(mean = mean(social_MCMC1$.value),
                              sd = sd(social_MCMC1$.value)),color="black") +
    geom_area(stat = "function", fun = dnorm,args = list(mean = mean(social_MCMC1$.value),
                                                         sd = sd(social_MCMC1$.value)),
              fill = "#969696", xlim = c(0,6),alpha=0.5) +
    geom_area(stat = "function", fun = dnorm,args = list(mean = mean(social_MCMC1$.value),
                                                         sd = sd(social_MCMC1$.value)),
              fill = "#d9d9d9", xlim = c(-6,0),alpha=0.5) +
    geom_vline(xintercept = 0, linetype = 3)+
    scale_y_continuous(breaks = NULL)+
    scale_x_continuous(limits = c(-6,6),
                       labels = c(seq(-6,6, by = 1)),
                       breaks = c(seq(-6,6, by = 1))) +
    mytheme() + 
    labs(x="Coefficient", 
         y = NULL, 
         title = "A - Vague Prior",
         subtitle = expression('Prior: Coef ~ N(0, 100'^2*')'))


Supf3b <- ggplot(data = data.frame(x = c(-6, 6)), aes(x)) +
    stat_function(fun = dnorm, n = 182,
                  args = list(mean = 0, sd = 1.020408),
                  color="black",linetype=2) +
    stat_function(fun = dnorm,
                  n = 182,
                  args = list(mean = mean(social_MCMC2$.value),
                              sd = sd(social_MCMC2$.value)),color="black") +
    geom_area(stat = "function", fun = dnorm,args = list(mean = mean(social_MCMC2$.value),
                                                         sd = sd(social_MCMC2$.value)),
              fill = "#969696", xlim = c(0,6),alpha=0.5) +
    geom_area(stat = "function", fun = dnorm,args = list(mean = mean(social_MCMC2$.value),
                                                         sd = sd(social_MCMC2$.value)),
              fill = "#d9d9d9", xlim = c(-6,0),alpha=0.5) +
    geom_vline(xintercept = 0, linetype = 3)+
    scale_y_continuous(breaks = NULL)+
    scale_x_continuous(limits = c(-6,6),
                       labels = c(seq(-6,6, by = 1)),
                       breaks = c(seq(-6,6, by = 1))) +
    mytheme() + 
    labs(x="Coefficient", 
         y = NULL, 
         title = "B - Equivocal Prior",
         subtitle = expression('Prior: Coef ~ N(0, 1.02'^2*')'))


Supf3c <- ggplot(data = data.frame(x = c(-6, 6)), aes(x)) +
    stat_function(fun = dnorm, n = 182,
                  args = list(mean = 2, sd = 1.020408),
                  color="black",linetype=2) +
    stat_function(fun = dnorm,
                  n = 182,
                  args = list(mean = mean(social_MCMC3$.value),
                              sd = sd(social_MCMC3$.value)),color="black") +
    geom_area(stat = "function", fun = dnorm,args = list(mean = mean(social_MCMC3$.value),
                                                         sd = sd(social_MCMC3$.value)),
              fill = "#969696", xlim = c(0,6),alpha=0.5) +
    geom_area(stat = "function", fun = dnorm,args = list(mean = mean(social_MCMC3$.value),
                                                         sd = sd(social_MCMC3$.value)),
              fill = "#d9d9d9", xlim = c(-6,0),alpha=0.5) +
    geom_vline(xintercept = 0, linetype = 3)+
    scale_y_continuous(breaks = NULL)+
    scale_x_continuous(limits = c(-6,6),
                       labels = c(seq(-6,6, by = 1)),
                       breaks = c(seq(-6,6, by = 1))) +
    mytheme() + 
    labs(x="Coefficient", 
         y = NULL, 
         title = "C - Optimistic Prior",
         subtitle = expression('Prior: Coef ~ N(2, 1.02'^2*')'))

Supf3d <- ggplot(data = data.frame(x = c(-6, 6)), aes(x)) +
    stat_function(fun = dnorm, n = 182,
                  args = list(mean = -2, sd = 1.020408),
                  color="black",linetype=2) +
    stat_function(fun = dnorm,
                  n = 182,
                  args = list(mean = mean(social_MCMC4$.value),
                              sd = sd(social_MCMC4$.value)),color="black") +
    geom_area(stat = "function", fun = dnorm,args = list(mean = mean(social_MCMC4$.value),
                                                         sd = sd(social_MCMC4$.value)),
              fill = "#969696", xlim = c(0,6),alpha=0.5) +
    geom_area(stat = "function", fun = dnorm,args = list(mean = mean(social_MCMC4$.value),
                                                         sd = sd(social_MCMC4$.value)),
              fill = "#d9d9d9", xlim = c(-6,0),alpha=0.5) +
    geom_vline(xintercept = 0, linetype = 3)+
    scale_y_continuous(breaks = NULL)+
    scale_x_continuous(limits = c(-6,6),
                       labels = c(seq(-6,6, by = 1)),
                       breaks = c(seq(-6,6, by = 1))) +
    mytheme() + 
    labs(x="Coefficient", 
         y = NULL, 
         title = "D - Pessimistic Prior",
         subtitle = expression('Prior: Coef ~ N(-2, 1.02'^2*')'))

gridExtra::grid.arrange(Supf3a, Supf3c, Supf3b, Supf3d, ncol = 2)

tiff("Supf3_Social.tiff", units="in", width= 12, height=7, res=600, compression = "lzw")
gridExtra::grid.arrange(Supf3a, Supf3c, Supf3b, Supf3d, ncol = 2)
dev.off()

#===============================
# Supplemental Fig. 4: Prior and posterior distributions for the relationship between 
# family functioning and environment aspect

Supf4a <- ggplot(data = data.frame(x = c(-6, 6)), aes(x)) +
    stat_function(fun = dnorm, n = 182,
                  args = list(mean = 0, sd = 100),
                  color="black",linetype=2) +
    stat_function(fun = dnorm,
                  n = 182,
                  args = list(mean = mean(envir_MCMC1$.value),
                              sd = sd(envir_MCMC1$.value)),color="black") +
    geom_area(stat = "function", fun = dnorm,args = list(mean = mean(envir_MCMC1$.value),
                                                         sd = sd(envir_MCMC1$.value)),
              fill = "#969696", xlim = c(0,6),alpha=0.5) +
    geom_area(stat = "function", fun = dnorm,args = list(mean = mean(envir_MCMC1$.value),
                                                         sd = sd(envir_MCMC1$.value)),
              fill = "#d9d9d9", xlim = c(-6,0),alpha=0.5) +
    geom_vline(xintercept = 0, linetype = 3)+
    scale_y_continuous(breaks = NULL)+
    scale_x_continuous(limits = c(-6,6),
                       labels = c(seq(-6,6, by = 1)),
                       breaks = c(seq(-6,6, by = 1))) +
    mytheme() + 
    labs(x="Coefficient", 
         y = NULL, 
         title = "A - Vague Prior",
         subtitle = expression('Prior: Coef ~ N(0, 100'^2*')'))


Supf4b <- ggplot(data = data.frame(x = c(-6, 6)), aes(x)) +
    stat_function(fun = dnorm, n = 182,
                  args = list(mean = 0, sd = 1.020408),
                  color="black",linetype=2) +
    stat_function(fun = dnorm,
                  n = 182,
                  args = list(mean = mean(envir_MCMC2$.value),
                              sd = sd(envir_MCMC2$.value)),color="black") +
    geom_area(stat = "function", fun = dnorm,args = list(mean = mean(envir_MCMC2$.value),
                                                         sd = sd(envir_MCMC2$.value)),
              fill = "#969696", xlim = c(0,6),alpha=0.5) +
    geom_area(stat = "function", fun = dnorm,args = list(mean = mean(envir_MCMC2$.value),
                                                         sd = sd(envir_MCMC2$.value)),
              fill = "#d9d9d9", xlim = c(-6,0),alpha=0.5) +
    geom_vline(xintercept = 0, linetype = 3)+
    scale_y_continuous(breaks = NULL)+
    scale_x_continuous(limits = c(-6,6),
                       labels = c(seq(-6,6, by = 1)),
                       breaks = c(seq(-6,6, by = 1))) +
    mytheme() + 
    labs(x="Coefficient", 
         y = NULL, 
         title = "B - Equivocal Prior",
         subtitle = expression('Prior: Coef ~ N(0, 1.02'^2*')'))


Supf4c <- ggplot(data = data.frame(x = c(-6, 6)), aes(x)) +
    stat_function(fun = dnorm, n = 182,
                  args = list(mean = 2, sd = 1.020408),
                  color="black",linetype=2) +
    stat_function(fun = dnorm,
                  n = 182,
                  args = list(mean = mean(envir_MCMC3$.value),
                              sd = sd(envir_MCMC3$.value)),color="black") +
    geom_area(stat = "function", fun = dnorm,args = list(mean = mean(envir_MCMC3$.value),
                                                         sd = sd(envir_MCMC3$.value)),
              fill = "#969696", xlim = c(0,6),alpha=0.5) +
    geom_area(stat = "function", fun = dnorm,args = list(mean = mean(envir_MCMC3$.value),
                                                         sd = sd(envir_MCMC3$.value)),
              fill = "#d9d9d9", xlim = c(-6,0),alpha=0.5) +
    geom_vline(xintercept = 0, linetype = 3)+
    scale_y_continuous(breaks = NULL)+
    scale_x_continuous(limits = c(-6,6),
                       labels = c(seq(-6,6, by = 1)),
                       breaks = c(seq(-6,6, by = 1))) +
    mytheme() + 
    labs(x="Coefficient", 
         y = NULL, 
         title = "C - Optimistic Prior",
         subtitle = expression('Prior: Coef ~ N(2, 1.02'^2*')'))

Supf4d <- ggplot(data = data.frame(x = c(-6, 6)), aes(x)) +
    stat_function(fun = dnorm, n = 182,
                  args = list(mean = -2, sd = 1.020408),
                  color="black",linetype=2) +
    stat_function(fun = dnorm,
                  n = 182,
                  args = list(mean = mean(envir_MCMC4$.value),
                              sd = sd(envir_MCMC4$.value)),color="black") +
    geom_area(stat = "function", fun = dnorm,args = list(mean = mean(envir_MCMC4$.value),
                                                         sd = sd(envir_MCMC4$.value)),
              fill = "#969696", xlim = c(0,6),alpha=0.5) +
    geom_area(stat = "function", fun = dnorm,args = list(mean = mean(envir_MCMC4$.value),
                                                         sd = sd(envir_MCMC4$.value)),
              fill = "#d9d9d9", xlim = c(-6,0),alpha=0.5) +
    geom_vline(xintercept = 0, linetype = 3)+
    scale_y_continuous(breaks = NULL)+
    scale_x_continuous(limits = c(-6,6),
                       labels = c(seq(-6,6, by = 1)),
                       breaks = c(seq(-6,6, by = 1))) +
    mytheme() + 
    labs(x="Coefficient", 
         y = NULL, 
         title = "D - Pessimistic Prior",
         subtitle = expression('Prior: Coef ~ N(-2, 1.02'^2*')'))

gridExtra::grid.arrange(Supf4a, Supf4c, Supf4b, Supf4d, ncol = 2)

tiff("Supf4_Envir.tiff", units="in", width= 12, height=7, res=600, compression = "lzw")
gridExtra::grid.arrange(Supf4a, Supf4c, Supf4b, Supf4d, ncol = 2)
dev.off()

#===============================================================================
