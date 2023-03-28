#-------------------------------------------------------
# A model fit for the gapminder data using the year as a 
# categorical feature ----------------------------------
#-------------------------------------------------------
# Exemplifying:
# Package loading  ----------------------------------
library(magrittr)
library(ggplot2)
library(tidymodels)
library(firatheme)
devtools::load_all(".")
load("data/gapminder_recent_g20.RData")

# Dataset split  ------------------------------------
set.seed(2022)

# Account for random things that happened in each year
df_real     <- gapminder_recent_g20 %>% 
  select(year, country, lifeExp, year0, decade0, continent, gdpPercap) |> 
  set_names(c('X1', 'country', 'y', "X2", "X3", "continent", "X4"))
#View(df_real)
df_real$X4 <- log(df_real$X4)

# df_real |> 
#   group_by(continent) |> 
#   summarise(n = n_distinct(country))
conts <- c("Americas", "Europe", "Asia")
df_real <- df_real |> filter(continent %in% conts)

# turn groups into numeric 
num_group <- data.frame(continent = unique(df_real$continent)) |> 
  mutate(num_continent = 1:n())
groups_table <- df_real |> 
  group_by(continent, country) |> 
  slice(1) |> 
  select(continent, country) |> 
  #mutate(num_continent = as.numeric(continent)) |> 
  group_by(continent) |> 
  mutate(num_country = 1:n()) |> 
  inner_join(num_group, by = 'continent') |> 
  mutate(num_country = paste0(num_continent, num_country))

df_real <- df_real |> 
  inner_join(groups_table, by = c("continent", "country"))
#-------------------------------------------------------------------------
#-------------------------------------------------------------------------

run_gapminder <- function(x){
  to_remove   <- sample(years, 15)
  train       <- df_real |> filter(!(X1 %in% to_remove))
  test        <- df_real |> filter(X1 %in% to_remove)
  num_trees   <- 6
  
  # Model parameters -----------------------------------
  group_variable <-  "num_continent"
  subgroup_variable <-  "num_country"
  formula        <- y ~ X1 + X2
  
  # Running the model ----------------------------------
  hb_model <- nhebart(formula,
                     data           = train,
                     group_variable = group_variable,
                     subgroup_variable = subgroup_variable,
                     num_trees = num_trees,
                     priors = list(
                       alpha = 0.95, # Prior control list
                       beta = 2,
                       nu = 2,
                       lambda = 0.1,
                       tau_mu = 16 * num_trees,
                       shape_sigma_phi = 0.5,
                       scale_sigma_phi = 1,
                       sample_sigma_phi = TRUE,
                       shape_sigma_gamma = 0.5,
                       scale_sigma_gamma = 1
                     ), 
                     inits = list(tau         = 1,
                                  tau_phi     = 3, 
                                  sigma_phi   = 3,
                                  sigma_gamma = 2),
                     MCMC = list(iter = 500, # Number of iterations
                                 burn = 100, # Size of burn in
                                 thin = 1,
                                 sigma_phi_sd   = 2,
                                 sigma_gamma_sd = 2)
  )
  
  pp <- predict_nhebart(newX = test, 
                        new_groups = test$num_continent, 
                        new_subgroups = test$num_country,
                        hebart_posterior = hb_model, 
                        type = "mean")
  
  rmse_mhebart <-  sqrt(mean((pp - test$y)^2))
  test$preds <- pp
  
  lme_ss <- lme4::lmer(y ~ X1 + X2 + (1|continent) + (1|continent/country), train)
  pplme <- predict(lme_ss, test)
  rmse_lmer <- sqrt(mean((pplme - test$y)^2)) # 3.991818
  
  return(list(train = train, test = test, 
              hb_model = hb_model, rmse_mhebart = rmse_mhebart, 
              pred_hebart = pp, 
              lme_model = lme_ss, rmse_lme = rmse_lmer, 
              pred_lme = pplme
  )) 
}

runs <- tibble(
  all = map(1:10, run_gapminder)
)

saveRDS(runs, "models_gapminder.rds")
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
# runs <- models_gapminder

# -------------------------------------------------------------
new_lme <- function(train, test){
  
  lme_ss <- me_ss <- lme4::lmer(y ~ X1 + X2 + (X1 | continent) + (X1 |continent:country), train)
  pplme <- predict(lme_ss, test)
  rmse_lmer <- sqrt(mean((pplme - test$y)^2)) # 3.991818
  
  return(list(lme_model_new = lme_ss, 
              rmse_lme_new = rmse_lmer, 
              pred_lme_new = pplme
  )) 
}

runs2 <- runs |> 
  mutate(train = map(all, "train"), 
         test = map(all, "test"), 
         new_lme = map2(train, test, new_lme)
  ) 

runs3 <-  runs2 |> 
  mutate(
    rmse_lme = map_dbl(new_lme,"rmse_lme_new"),
    rmse_nhebart = map_dbl(all,"rmse_mhebart"),
    pred_lme = map(new_lme, "pred_lme_new"),
    pred_nhebart = map(all, "pred_hebart"),
    id = 1:n(), 
    test = map(all, "test")
  )
runs3$rmse_lme
runs3$rmse_nhebart

#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
runs3 <-  runs3  |> 
  unnest(c(pred_lme, pred_nhebart, test)) |> 
  dplyr::select(y, pred_nhebart, pred_lme, rmse_lme, 
                rmse_nhebart, country, continent, X1, X2) |> 
  mutate(sd_lme = sd(pred_lme)/sqrt(n()), 
         sd_nhebart = sd(pred_nhebart)/sqrt(n())) |> 
  group_by(X1, X2, country, continent) |> 
  reframe(
    y = mean(y), 
    pred_lme_mean = mean(pred_lme), 
    pred_nhebart_mean = mean(pred_nhebart), 
    low_lme = pred_lme_mean - 1.95 * sd_lme,
    upp_lme = pred_lme_mean + 1.95 * sd_lme,
    low_nhebart = pred_nhebart_mean - 1.95 * sd_nhebart,
    upp_nhebart = pred_nhebart_mean + 1.95 * sd_nhebart
  ) |> 
  group_by(X1, X2, country, continent) |> 
  slice(1)


# Plots -------------------------------------------------------
avg_rmse_lme <- mean(runs2$rmse_lme)
upp_lme <- avg_rmse_lme + 1.96 * sd(runs2$rmse_lme)/sqrt(10)
low_lme <- avg_rmse_lme - 1.96 * sd(runs2$rmse_lme)/sqrt(10)


avg_rmse <- mean(runs2$rmse_nhebart)
upp <- avg_rmse + 1.96 * sd(runs2$rmse_nhebart)/sqrt(10)
low <- avg_rmse - 1.96 * sd(runs2$rmse_nhebart)/sqrt(10)


BottleRocket2 = c("#FAD510", "#CB2314", "#0E86D4",
                           "#1E1E1E", "#18A558")
                           
selected_countries <- c(
  "South Africa", "Russia", 
  "Turkey", "Indonesia", "Brazil",
  "Argentina", "France", "United States"
)


runs3 |>  
  filter(country %in% selected_countries) |> 
  ggplot(aes(x = X1, y = pred_nhebart_mean)) +
  facet_wrap(~country, ncol = 2, scales = 'free_y') +
  geom_ribbon(aes(ymin=low_nhebart, ymax=upp_nhebart),
              fill = BottleRocket2[3], alpha = 0.5) + 
  geom_ribbon(aes(ymin=low_lme, ymax=upp_lme),
              fill = BottleRocket2[2], alpha = 0.3) + 
  geom_line(aes(colour = BottleRocket2[2]), size = 0.7) +
  geom_line(aes(y = pred_lme_mean), 
            colour = BottleRocket2[2], size = 0.7) +
  geom_line(aes(colour = BottleRocket2[3]), size = 0.7) +
  geom_line(colour = BottleRocket2[3], size = 0.7) +
  geom_point(aes(x = X1, y = y, colour =  'black'), size = 0.25) + 
  geom_point(aes(x = X1, y = y), 
             colour =  'black', size = 0.25) + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 7)) +
  labs(y = "Life expectancy (in years)", 
       x = 'Passage of time', 
       title = 
         
         paste0("Average predictions per country for Nested-RE HEBART and LME \n", 
                "Nested-RE HEBART Test RMSE: ", 
                paste0(round(avg_rmse, 2), " [", round(low, 2), ",", round(upp, 2), "]"), "\n", 
                "LME Test RMSE: ", 
                paste0(round(avg_rmse_lme, 2), " [", round(low_lme, 2), ",", round(upp_lme, 2), "]"), "\n"
         )) + 
  theme_linedraw(12) +
  scale_colour_manual(
    name="Source:",
    values = c("black", BottleRocket2[3], BottleRocket2[2]), 
    labels = c("Data", "Nested-RE HEBART", "LME"), 
    
    guide = guide_legend(override.aes = list(
      size = c(2, 2, 2)))) +
  theme(panel.spacing.x = unit(0.5, "lines"), 
        legend.position = "bottom")

ggsave(file = "examples/images/predictions_plot_gapminder_nhebart.png",
       width = 8, height = 8)
# ------------------------------------------------------------------------------
# This should be predictions per country
# ------------------------------------------------------------------------------

pred_continent <- function(train, test, hb_model){
  pp <- predict_nhebart(newX = test, 
                        new_groups = test$num_continent, 
                        new_subgroups = NULL,
                        hebart_posterior = hb_model, 
                        type = "mean")
  
  rmse_nhebart <-  sqrt(mean((pp - test$y)^2))
  
  suppressWarnings(
  lme_ss <- lme4::lmer(y ~ X1 + X2 + (X1 |continent), train)
  )
  pplme <- predict(lme_ss, test)
  rmse_lmer <- sqrt(mean((pplme - test$y)^2)) # 3.991818
  
  return(list(
    pred_nhebart = pp, pred_lme = pplme, 
    rmse_nhebart = rmse_nhebart, rmse_lme = rmse_lmer
  )) 
  
} 
# ------------------------------------------------------------------------------
runs_cont <-  runs |> 
  mutate(
    id = 1:n(), 
    test = map(all, "test"),
    train = map(all, "train"),
    hb_model = map(all, "hb_model"),
    continent = pmap(list(train, test, hb_model), pred_continent)
  )

runs2 <-  runs_cont |> 
  mutate(
    rmse_lme = map_dbl(continent,"rmse_lme"),
    rmse_nhebart = map_dbl(continent,"rmse_nhebart"),
    pred_lme = map(continent, "pred_lme"),
    pred_nhebart = map(continent, "pred_nhebart")
  )

# ------------------------------------------------------------------------------

runs4 <-  runs2  |> 
  select(-continent) |> 
  unnest(c(pred_lme, pred_nhebart, test)) |> 
  dplyr::select(y, pred_nhebart, pred_lme, rmse_lme, 
                rmse_nhebart, X1, X2, continent) |> 
  mutate(sd_lme = sd(pred_lme)/sqrt(n()), 
         sd_nhebart = sd(pred_nhebart)/sqrt(n())) |> 
  group_by(X1, X2, continent) |> 
  reframe(
    y = mean(y), 
    pred_lme_mean = mean(pred_lme), 
    pred_nhebart_mean = mean(pred_nhebart), 
    low_lme = pred_lme_mean - 1.95 * sd_lme,
    upp_lme = pred_lme_mean + 1.95 * sd_lme,
    low_mhebart = pred_nhebart_mean - 1.95 * sd_nhebart,
    upp_mhebart = pred_nhebart_mean + 1.95 * sd_nhebart
  ) |> 
  group_by(X1, X2, continent) |> 
  slice(1)


avg_rmse_lme <- mean(runs2$rmse_lme)
upp_lme <- avg_rmse_lme + 1.96 * sd(runs2$rmse_lme)/sqrt(10)
low_lme <- avg_rmse_lme - 1.96 * sd(runs2$rmse_lme)/sqrt(10)


avg_rmse <- mean(runs2$rmse_nhebart)
upp <- avg_rmse + 1.96 * sd(runs2$rmse_nhebart)/sqrt(10)
low <- avg_rmse - 1.96 * sd(runs2$rmse_nhebart)/sqrt(10)


BottleRocket2 = c("#FAD510", "#CB2314", "#0E86D4",
                           "#1E1E1E", "#18A558")
                          

runs4 |>  
  ggplot(aes(x = X1, y = pred_nhebart_mean)) +
  facet_wrap(~continent, ncol = 3, scales = 'free_y') +
  geom_ribbon(aes(ymin=low_mhebart, ymax=upp_mhebart),
              fill = BottleRocket2[3], alpha = 0.5) + 
  geom_ribbon(aes(ymin=low_lme, ymax=upp_lme),
              fill = BottleRocket2[2], alpha = 0.3) + 
  geom_line(aes(colour = BottleRocket2[2]), size = 0.7) +
  geom_line(aes(y = pred_lme_mean), 
            colour = BottleRocket2[2], size = 0.7) +
  geom_line(aes(colour = BottleRocket2[3]), size = 0.7) +
  geom_line(colour = BottleRocket2[3], size = 0.7) +
  geom_point(aes(x = X1, y = y, colour =  'black'), size = 0.25) + 
  geom_point(aes(x = X1, y = y), 
             colour =  'black', size = 0.25) + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 7)) +
  labs(y = "Life expectancy (in years)", 
       x = 'Passage of time', 
       title = 
         
         paste0("Average predictions per continent for Nested-RE HEBART and LME \n", 
                "Nested-RE HEBART Test RMSE: ", 
                paste0(round(avg_rmse, 2), " [", round(low, 2), ",", round(upp, 2), "]"), "\n", 
                "LME Test RMSE: ", 
                paste0(round(avg_rmse_lme, 2), " [", round(low_lme, 2), ",", round(upp_lme, 2), "]"), "\n"
         )) + 
  theme_linedraw(12) +
  scale_colour_manual(
    name="Source:",
    values = c("black", BottleRocket2[3], BottleRocket2[2]), 
    labels = c("Data", "Nested-RE HEBART", "LME"), 
    
    guide = guide_legend(override.aes = list(
      size = c(2, 2, 2)))) +
  theme(panel.spacing.x = unit(0.5, "lines"), 
        legend.position = "bottom")

ggsave(file = "examples/images/predictions_plot_gapminder_continent_nhebart.png",
       width = 8, height = 4)

