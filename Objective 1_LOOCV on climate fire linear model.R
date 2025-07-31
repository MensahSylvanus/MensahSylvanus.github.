
data <- read.csv("fireplotdata.csv")
head(data)
str(data)
names(data)


library(nlme)
library(dplyr)
library(ggplot2)

# Scaling predictors
data <- data %>%
  mutate(
    rainfall_seasonality_scaled = scale(rainfall_seasonality)[,1],
    rainfall_scaled = scale(rainfall)[,1],
    temperature_scaled = scale(temperature)[,1]
  )

n <- nrow(data)
results_list <- vector("list", length = n)
normality_results <- data.frame(fold = integer(), shapiro_p = numeric())

for (i in seq_len(n)) {
  cat("Fitting LOOCV for observation", i, "\n")
  
  train_data <- data[-i, ]
  test_data  <- data[i, , drop = FALSE]
  
  gls_model <- gls(
    log(Fire_frequency + 1) ~ rainfall_seasonality_scaled+rainfall_scaled+temperature_scaled,
    correlation = corExp(form = ~ Longitude + Latitude),
    data = train_data,
    method = "ML"
  )
  
  # Residuals on training set (log scale)
  train_pred <- predict(gls_model, newdata = train_data)
  train_resid <- log(train_data$Fire_frequency + 1) - train_pred
  
  # Shapiro-Wilk normality test on residuals for each fold
  shapiro_test <- shapiro.test(train_resid)
  normality_results <- rbind(normality_results, data.frame(fold = i, shapiro_p = shapiro_test$p.value))
  
  pred_log <- predict(gls_model, newdata = test_data)
  
  results_list[[i]] <- data.frame(
    observed = log(test_data$Fire_frequency + 1),
    predicted_raw = pred_log,
    obs_index = i
  )
}

all_loocv <- bind_rows(results_list)

# Calculating residuals on log scale
all_loocv <- all_loocv %>%
  mutate(residual = observed - predicted_raw)

# Computing Residual Standard Error (RSE)
RSE <- sd(all_loocv$residual)
CF <- exp((RSE^2)/2)
cat("Residual Standard Error (RSE):", RSE, "\n")
cat("Correction Factor (CF):", CF, "\n")

# Back-transforming predictions with correction factor applied
all_loocv <- all_loocv %>%
  mutate(
    predicted = (exp(predicted_raw) - 1) * CF,
    observed_orig = exp(observed) - 1,
    residual_orig = observed_orig - predicted
  )

# Plot residuals vs corrected predicted values on original scale
ggplot(all_loocv, aes(x = predicted, y = residual_orig)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_point(size=2, color="#00AFBB")+
 # geom_smooth(method="lm", size=1,linetype=1, se=T,color = "#FC4E07")+
  labs(
    title = "Residuals vs Predicted Fire Frequency (Spatial GLS model)",
    x = "Predicted Fire Frequency (corrected)",
    y = "Residuals (Observed - Predicted)"
  ) ->b;b

# Shapiro-Wilk p-values per fold
print(normality_results)

cat("Folds with residuals passing normality test (p > 0.05):", 
    sum(normality_results$shapiro_p > 0.05), "out of", n, "\n")





library(broom)   # for tidy model summaries

n <- nrow(data)
results_list <- vector("list", length = n)
normality_results <- data.frame(fold = integer(), shapiro_p = numeric())

# Creating list to store model summaries per fold
model_summaries <- vector("list", length = n)

for (i in seq_len(n)) {
  cat("Fitting LOOCV for observation", i, "\n")
  
  train_data <- data[-i, ]
  test_data  <- data[i, , drop = FALSE]
  
  gls_model <- gls(
    log(Fire_frequency + 1) ~ rainfall_seasonality_scaled+rainfall_scaled+temperature_scaled,
    correlation = corExp(form = ~ Longitude + Latitude),
    data = train_data,
    method = "ML"
  )
  
  # Extracting coefficients, standard errors, t-values and p-values:
  coef_table <- summary(gls_model)$tTable
  coef_df <- as.data.frame(coef_table)
  coef_df$term <- rownames(coef_df)
  rownames(coef_df) <- NULL
  coef_df$fold <- i
  
  # Adding to model_summaries list
  model_summaries[[i]] <- coef_df[, c("fold", "term", "Value", "Std.Error", "t-value", "p-value")]
  
  # Computing R² (pseudo R² by squared correlation between observed and fitted)
  fitted_vals <- predict(gls_model)
  r2 <- cor(log(train_data$Fire_frequency + 1), fitted_vals)^2
  model_summaries[[i]]$R2 <- r2
  
  # Storing residuals for normality test
  train_resid <- log(train_data$Fire_frequency + 1) - fitted_vals
  shapiro_test <- shapiro.test(train_resid)
  normality_results <- rbind(normality_results, data.frame(fold = i, shapiro_p = shapiro_test$p.value))
  
  pred_log <- predict(gls_model, newdata = test_data)
  
  results_list[[i]] <- data.frame(
    observed = log(test_data$Fire_frequency + 1),
    predicted_raw = pred_log,
    obs_index = i
  )
}

all_loocv <- bind_rows(results_list)
all_summaries <- bind_rows(model_summaries)

# Summaries
print(all_summaries)

# Summary stats for R2 across folds
cat("Mean R2 across folds:", mean(unique(all_summaries$R2)), "\n")


library(dplyr)
library(tidyr)

# Summarizing coefficients across folds
coef_summary <- all_summaries %>%
  group_by(term) %>%
  summarise(
    mean_estimate = mean(Value),
    min_estimate = min(Value),
    max_estimate = max(Value),
    sd_estimate = sd(Value),
    mean_se = mean(Std.Error),
    n = n(),
    lower_ci = mean_estimate - 1.96 * mean_se,
    upper_ci = mean_estimate + 1.96 * mean_se,
    mean_t = mean(`t-value`),
    mean_p = mean(`p-value`),
    .groups = "drop"
  )

print(coef_summary)

library(dplyr)
library(ggplot2)

# Renaming terms
all_summaries_clean <- all_summaries %>%
  filter(term != "(Intercept)") %>%
  mutate(term = recode(term,
                       "rainfall_seasonality_scaled" = "Rainfall Seasonality",
                       "rainfall_scaled" = "Rainfall",
                       "temperature_scaled" = "Temperature"
  ))

# Calculating mean estimate and 95% CI per term across LOOCV folds
summary_stats_clean <- all_summaries_clean %>%
  group_by(term) %>%
  summarise(
    mean_est = mean(Value),
    se_est = sd(Value) / sqrt(n()),
    lower = mean_est - 1.96 * se_est,
    upper = mean_est + 1.96 * se_est,
    .groups = "drop"
  )

# Plot coefficient distributions across LOOCV folds
ggplot(all_summaries_clean, aes(x = term, y = Value, fill = term)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_jitter(width = 0.1, color = "black", size = 1.5, alpha = 0.7) +
  geom_point(
    data = summary_stats_clean,
    aes(x = term, y = mean_est),
    color = "red",
    size = 3,
    inherit.aes = FALSE
  ) +
  geom_errorbar(
    data = summary_stats_clean,
    aes(x = term, ymin = lower, ymax = upper),
    width = 0.2,
    color = "red",
    inherit.aes = FALSE
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
  labs(
    title = "",
    x = "",
    y = "Estimate"
  ) +coord_flip()+theme_bw()+
  theme(legend.position = "none")->d
d
