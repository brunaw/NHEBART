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
set.seed(2023)
years <- unique(df_real$X1)
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
                      MCMC = list(iter = 20, # Number of iterations
                                  burn = 10, # Size of burn in
                                  thin = 1,
                                  sigma_phi_sd   = 2,
                                  sigma_gamma_sd = 2)
  )
  
pp <- predict_nhebart(newX = test, 
                      new_groups = test$num_continent, 
                      new_subgroups = test$num_country,
                      hebart_posterior = hb_model, 
                      type = "mean")

rmse_nhebart <-  sqrt(mean((pp - test$y)^2))
test$preds <- pp

lme_ss <- lme4::lmer(y ~ X1 + X2 + (X1|continent) + (X1|continent:country), train)
pplme <- predict(lme_ss, test)
rmse_lmer <- sqrt(mean((pplme - test$y)^2)) 

c(rmse_lmer, rmse_nhebart)

test$preds_lme <- pplme



selected_countries <- c(
  "South Africa", "Russia", 
  "Turkey", "Indonesia", "Brazil",
  "Argentina", "France", "United States"
)



BottleRocket2 = c("#FAD510", "#CB2314", "#0E86D4",
                           "#1E1E1E", "#18A558")
                           
test |>  
  filter(country %in% selected_countries) |> 
  ggplot(aes(x = X1, y = preds)) +
  facet_wrap(~country, ncol = 2, scales = 'free_y') +
  geom_line(aes(colour = BottleRocket2[3]), size = 0.7) +
  geom_line(colour = BottleRocket2[3], size = 0.7) +
  geom_line(aes(y = preds_lme, colour = BottleRocket2[2]), size = 0.7) +
  geom_line(aes(y = preds_lme), 
            colour = BottleRocket2[2], size = 0.7) +
  geom_point(aes(x = X1, y = y), colour =  'black', size = 0.25) + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 7)) +
  labs(y = "Life expectancy (in years)", 
       x = 'Passage of time') + 
  theme_linedraw(12) +
  scale_colour_manual(
    name="Source:",
    values = c(BottleRocket2[3], BottleRocket2[2]), 
    labels = c("Nested-RE HEBART", "LME"), 
    
    guide = guide_legend(override.aes = list(
      size = c(2, 2)))) +
  theme(panel.spacing.x = unit(0.5, "lines"), 
        legend.position = "bottom")
#----------------------------------------------------------------------------
