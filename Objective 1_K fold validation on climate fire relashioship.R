rm(list=ls(all=T))
graphics.off()
## > 0. Packages and working directory
libs <- c('data.table','tibble', 'dplyr',"tidyr", 'readr','stringr','imputeTS','dplR','tidyverse',"lubridate","hrbrthemes",
          'car','nlme','lme4','ggplot2',"sjPlot", 'ggpubr', 'ggeffects',"gridExtra", "lmerTest","grid","cowplot")
invisible(lapply(libs, library, character.only = T))


setwd("C:/Users/Mensah/Documents/Mensah 2023-Research/fire plot 1/Analaysis") #set wd

data <- read.csv("fireplotdata.csv")
head(data)
str(data)
names(data)


library(dplyr)
library(spaMM)
library(sf)
library(Metrics)


data_sf <- st_as_sf(data, coords = c("Longitude", "Latitude"), crs = 4326)# Converting data to sf object

# Creating folds with cv_spatial()
set.seed(123)
sb <- cv_spatial(
  x = data_sf,
  k = 5,
  selection = "random",
)


folds_list <- sb$folds_list  # list of train/test indices
results_nb <- list()
fitted_nb_models <- list()  # initializing model storage

for (i in seq_along(folds_list)) {
  cat("Fitting fold", i, "\n")
  
  train_idx <- folds_list[[i]][[1]]
  test_idx  <- folds_list[[i]][[2]]
  
  train_data <- data_sf[train_idx, ]
  test_data  <- data_sf[test_idx, ]
  
  train_data$X <- st_coordinates(train_data)[,1]
  train_data$Y <- st_coordinates(train_data)[,2]
  test_data$X  <- st_coordinates(test_data)[,1]
  test_data$Y  <- st_coordinates(test_data)[,2]
  
  train_df <- st_drop_geometry(train_data)
  test_df  <- st_drop_geometry(test_data)
  
  ## Fitting Spatial NB
  nb_model <- fitme(
    Fire_frequency ~ scale(rainfall_seasonality)+scale(rainfall)+scale(temperature) + Matern(1 | X + Y),
    data = train_df,
    family = "negbin",
    control = list(fix_predVar = TRUE)
  )
  
  # Saving model
  fitted_nb_models[[i]] <- nb_model
  
  # Predicting & storing
  nb_preds <- predict(nb_model, newdata = test_df, type = "response")
  results_nb[[i]] <- data.frame(
    observed = test_df$Fire_frequency,
    predicted = nb_preds,
    fold = i,
    model = "Spatial NB"
  )
}

for (i in seq_along(fitted_nb_models)) {
     cat("\n==== Fold", i, "====\n")
     print(summary(fitted_nb_models[[i]]))
   }

library(ggplot2)
library(dplyr)

# Combining results from all folds into one data.frame
all_results_nb <- bind_rows(results_nb)

# Calculating residuals
all_results_nb <- all_results_nb %>%
  mutate(residual = observed - predicted)

# Plot residuals vs predicted values
ggplot(all_results_nb, aes(x = predicted, y = residual)) +
  geom_point(size=2, color="#00AFBB")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "Residuals vs Predicted Fire Frequency (Spatial NB model)",
    x = "Predicted Fire Frequency",
    y = "Residuals (Observed - Predicted)"
  ) ->a;a




library(dplyr)

# List to store tidy summaries per fold
coef_summaries <- vector("list", length = length(fitted_nb_models))

for (i in seq_along(fitted_nb_models)) {
  cat("\n==== Fold", i, "====\n")
  
  model_summary <- summary(fitted_nb_models[[i]])
  print(model_summary)
  
  # Extracting fixed effects table (beta_table) if it exists
  if (!is.null(model_summary$beta_table)) {
    coef_df <- as.data.frame(model_summary$beta_table)
    coef_df$term <- rownames(coef_df)
    coef_df$fold <- i
    
    # Reordering columns for readability while keeping t-value
    coef_df <- coef_df %>% select(fold, term, Estimate, `Cond. SE`, `t-value`)
    
    coef_summaries[[i]] <- coef_df
  } else {
    warning(paste("No beta_table found for fold", i))
  }
}

# Combining all folds coefficients into one data.frame
all_coefs <- bind_rows(coef_summaries)

# Summarizing across folds
summary_table <- all_coefs %>%
  group_by(term) %>%
  summarise(
    mean_estimate = mean(Estimate),
    min_estimate = min(Estimate),
    max_estimate = max(Estimate),
    pooled_se = sqrt(mean(`Cond. SE`^2)),
    lower_95 = mean_estimate - 1.96 * pooled_se,
    upper_95 = mean_estimate + 1.96 * pooled_se,
    mean_t_value = mean(`t-value`),
    approx_p_value = 2 * (1 - pnorm(abs(mean_t_value))),
    folds = n(),
    .groups = "drop"
  ) %>%
  arrange(term)

print(summary_table)


library(ggplot2)

library(ggplot2)
library(dplyr)

all_coefs_no_intercept <- all_coefs %>% filter(term != "(Intercept)")

ggplot(all_coefs_no_intercept, aes(x = term, y = Estimate, fill = term)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_jitter(width = 0.1, color = "black", size = 1.5, alpha = 0.8) +
  labs(
    title = "",
    x = "Predictor",
    y = "Estimate"
  ) +coord_flip()+
  theme(legend.position = "none")


# Summarizing stats for each predictor term (excluding Intercept)
summary_stats <- all_coefs_no_intercept %>%
  group_by(term) %>%
  summarise(
    mean_est = mean(Estimate),
    lower = mean_est - 1.96 * sd(Estimate),
    upper = mean_est + 1.96 * sd(Estimate),
    .groups = "drop"
  )


all_coefs_clean <- all_coefs_no_intercept %>%
  mutate(term = recode(term,
                       "scale(rainfall_seasonality)" = "Rainfall seasonality",
                       "scale(rainfall)" = "Rainfall",
                       "scale(temperature)" = "Temperature"
  ))

summary_stats_clean <- summary_stats %>%
  mutate(term = recode(term,
                       "scale(rainfall_seasonality)" = "Rainfall seasonality",
                       "scale(rainfall)" = "Rainfall",
                       "scale(temperature)" = "Temperature"
  ))

# Plot
ggplot(all_coefs_clean, aes(x = term, y = Estimate, fill = term)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_jitter(width = 0.1, color = "black", size = 1.5, alpha = 0.7) +
  geom_point(data = summary_stats_clean, aes(x = term, y = mean_est), color = "red", size = 3, inherit.aes = FALSE) +
  geom_errorbar(
    data = summary_stats_clean,
    aes(x = term, ymin = lower, ymax = upper),
    width = 0.2,
    color = "red",
    inherit.aes = FALSE
  ) + geom_hline(yintercept = 0, linetype = "dashed", color = "gray30")+
  labs(
    title = "",
    x = "",
    y = "Estimate"
  ) +coord_flip()+theme_bw()+
  theme(legend.position = "none")->c
c
