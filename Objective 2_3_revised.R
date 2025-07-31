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

----------------------------------------------------------------------
### Obj. 2. Effects of fire frequency on tree diversity metrics
### Obj. 3. Effects of fire frequency on stand attributes 
----------------------------------------------------------------------

  ### Species richness

#checking for overdispersion
Fire_sp_mod0<-glm(Species_richness~Fire_frequency, family=poisson, data)
summary(Fire_sp_mod0)
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model, type = "pearson")
  sum(rp^2) / rdf
}
overdisp_fun(Fire_sp_mod0)### [1] 2.990721 showing overdispersion
library(performance)
check_overdispersion(Fire_sp_mod0)
# Overdispersion test

#dispersion ratio =   2.991
#Pearson's Chi-Squared = 239.258
#                p-value = < 0.001

#Overdispersion detected.

## Dealing with overdispersion
Fire_sp_mod01<-MASS::glm.nb(Species_richness~Fire_frequency, data)
summary(Fire_sp_mod01)
overdisp_fun(Fire_sp_mod01)### [1] 0.9535338 showing no overdispersion
library(performance)
check_overdispersion(Fire_sp_mod01)
# Overdispersion test

#dispersion ratio =  0.954
#Pearson's Chi-Squared = 76.283
#                p-value =  0.597

#No overdispersion detected.


## Testing model with random effects
Fire_sp_mod02<-glmer.nb(Species_richness~Fire_frequency+(1|Plot_size)+(1|Vegetation), data)
Fire_sp_mod02_reduced <- glmer.nb(Species_richness ~ Fire_frequency+(1|Plot_size), data = data)
Fire_sp_mod03_reduced <- glmer.nb(Species_richness ~ Fire_frequency+(1|Vegetation), data = data)
summary(Fire_sp_mod02)
anova(Fire_sp_mod02, Fire_sp_mod02_reduced)
anova(Fire_sp_mod02, Fire_sp_mod03_reduced)
anova(Fire_sp_mod02_reduced, Fire_sp_mod03_reduced)

## Testing spatial autocorrelation

Fire_sp_spatial <- fitme(
  Species_richness ~ Fire_frequency + (1|Plot_size)+ Matern(1|Longitude+Latitude),data = data,family = negbin())
summary(Fire_sp_spatial)

Fire_sp_nospatial <- fitme(
  Species_richness ~ Fire_frequency + (1|Plot_size),  data = data,  family = negbin())
summary(Fire_sp_nospatial)

logLik(Fire_sp_spatial)
logLik(Fire_sp_nospatial)

AIC(Fire_sp_spatial)
AIC(Fire_sp_nospatial )

### The spatial model (cAIC = 460.72) fits better than the non-spatial one (465.03)


### For the remaining diversity metrics, see "Testing for random effects" file; random effects were dropped in the subsequent analyses.
### Retained models for cross validation:

mod1 <- fitme(Species_richness ~ Fire_frequency + Matern(1|Longitude+Latitude),data = data,family = negbin())
mod2 <- gls(responses ~ Fire_frequency, correlation = corExp(form = ~ Longitude + Latitude),
  data = train_data,method = "ML")

  
### Here we performed a spatial cross-validation loop to fit two types of models depending on the response variable:
### For Species_richness: we used spaMM::fitme() with a spatial negative binomial model.
### For all other responses: we used nlme::gls() with a spatial correlation structure.
 
 
 library(spaMM)
 library(nlme)
 library(sf)
 library(dplyr)
 library(purrr)

responses <- c("Species_richness","FDis", "sesMPD", "sesMNTD", "cvDBH", "skwDBH", 
               "Tree_density", "LargeTrees_density", "RemTrees_density", "HLorey")### 
predictor <- "Fire_frequency"

# Converting data to sf object
data_sf <- st_as_sf(data, coords = c("Longitude", "Latitude"), crs = 4326)

# Create folds with cv_spatial()
#set.seed(123)
#sb <- cv_spatial(x = data_sf,k = 5, selection = "random")# not to run because already implemented in objective 1:k fold cross validation

folds_list <- sb$folds_list 

# Storing results
cv_results <- list()
fitted_models <- list()
coef_results <- list()
rmse_results <- list()


for (response in responses) {
  cat("\n==== Fitting models for:", response, "====\n")
  
  response_results <- list()
  model_list <- list()
  coef_list <- list()
  rmse_list <- list()
  
  for (i in seq_along(folds_list)) {
    cat("  Fold", i, "\n")
    
    train_idx <- folds_list[[i]][[1]]
    test_idx  <- folds_list[[i]][[2]]
    
    train_data <- data_sf[train_idx, ]
    test_data  <- data_sf[test_idx, ]
    
    train_data$Longitude <- st_coordinates(train_data)[, 1]
    train_data$Latitude  <- st_coordinates(train_data)[, 2]
    test_data$Longitude <- st_coordinates(test_data)[, 1]
    test_data$Latitude  <- st_coordinates(test_data)[, 2]
    
    train_df <- st_drop_geometry(train_data)
    test_df  <- st_drop_geometry(test_data)
    
    # Fitting model
    if (response == "Species_richness") {
      formula_nb <- as.formula(paste(response, "~", predictor, "+ Matern(1 | Longitude + Latitude)"))
      model <- fitme(formula_nb, data = train_df, family = negbin())
      preds <- predict(model, newdata = test_df, type = "response")
      
      # Extracting fixed effects
      beta_table <- as.data.frame(summary(model)$beta_table)
      beta_table$term <- rownames(beta_table)
      beta_table$response <- response
      beta_table$fold <- i
      # Filter out intercept
      beta_table <- beta_table %>% dplyr::filter(term != "(Intercept)")
      
      coef_list[[i]] <- beta_table
      
    } else {
      vars_to_check <- c(response, predictor, "Longitude", "Latitude")
      train_df_clean <- train_df %>%
        dplyr::select(all_of(vars_to_check)) %>%
        drop_na()
      
      if (nrow(train_df_clean) < 5) {
        warning(paste("Skipping fold", i, "for response", response, "due to too few observations after NA removal."))
        next
      }
      
      formula_gls <- as.formula(paste(response, "~", "scale(",predictor,")"))
      model <- gls(
        formula_gls,
        data = train_df_clean,
        correlation = corExp(form = ~ Longitude + Latitude),
        method = "ML"
      )
      preds <- predict(model, newdata = test_df)
      
      # Extracting coefficients
      model_summary <- summary(model)
      coef_table <- as.data.frame(model_summary$tTable)
      coef_table$term <- rownames(coef_table)
      coef_table$response <- response
      coef_table$fold <- i
      # Filtering out intercept
      coef_table <- coef_table %>% dplyr::filter(term != "(Intercept)")
      
      coef_list[[i]] <- coef_table
    }
    
    # Computing residuals and RMSE
    obs_vals <- test_df[[response]]
    rmse_val <- sqrt(mean((obs_vals - preds)^2, na.rm = TRUE))
    
    rmse_list[[i]] <- data.frame(
      response = response,
      fold = i,
      RMSE = rmse_val
    )
    
    # Storing predictions
    response_results[[i]] <- data.frame(
      observed = obs_vals,
      predicted = preds,
      fold = i,
      model = ifelse(response == "Species_richness", "fitme", "gls"),
      response = response
    )
    
    model_list[[i]] <- model
  }
  
  # Combining and storing per-response
  cv_results[[response]] <- bind_rows(response_results)
  fitted_models[[response]] <- model_list
  coef_results[[response]] <- bind_rows(coef_list)
  rmse_results[[response]] <- bind_rows(rmse_list)
}


all_cv_results <- bind_rows(cv_results);all_cv_results
all_coef_list  <- bind_rows(coef_results);all_coef_list
all_rmse_list  <- bind_rows(rmse_results);all_rmse_list


### Residual Plotting, Coefficient Summary, and RMSE Calculation
library(ggplot2)
library(dplyr)

all_cv_results <- all_cv_results %>%
  mutate(residual = observed - predicted)

# Plot residuals vs predicted per response
ggplot(all_cv_results, aes(x = predicted, y = residual)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  facet_wrap(~response, scales = "free") +
  labs(
    title = "Residuals vs Predicted Values by Response",
    x = "Predicted",
    y = "Residual"
  ) 



#####

all_coef_list <- all_coef_list %>%
  mutate(
    coef_estimate = ifelse(response == "Species_richness", Estimate, Value),
    coef_se = ifelse(response == "Species_richness", `Cond. SE`, Std.Error));all_coef_list


plot_df <- all_coef_list %>%
  filter(term == "scale(Fire_frequency)") %>%
  mutate(
    lower = coef_estimate - coef_se,
    upper = coef_estimate + coef_se
  );plot_df


plot_df$response <- factor(plot_df$response,
                           levels = c("Species_richness", "FDis", "sesMPD", "sesMNTD",
                                      "cvDBH", "skwDBH", "Tree_density", "LargeTrees_density", "RemTrees_density", "HLorey"),
                           labels = c("Species richness", "FDis", "sesMPD", "sesMNTD",
                                      "cvDBH", "skwDBH", "Tree density", "Large tree density", "Remaining tree density","Lorey's height")
                          )

plot_df <- plot_df %>%
  mutate(x_jitter = as.numeric(factor(response)) + runif(n(), -0.05, 0.05))

ggplot(plot_df, aes(x = reorder(response, coef_estimate), y = coef_estimate)) +
  geom_violin(aes(x = response, y = coef_estimate, fill = response),
              alpha = 0.8, trim = FALSE, color = NA) +
  
  geom_point(aes(x = x_jitter), color = "black", size = 1.5, alpha = 0.7) +  # Replaces geom_jitter
  geom_errorbar(aes(x = x_jitter, ymin = lower, ymax = upper), 
                width = 0.1, color = "grey50") +  # Aligned with points
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
  
  theme(
    axis.text.x = element_text(angle = 40, hjust = 1)
  )+labs( x = NULL,    y = "Slope estimate",
          title = ""
  ) +
  theme(legend.position = "none")



################# Re-running models on overall dataset

head(data_sf)

# Initializing result storage
model_list <- list()
coef_list <- list()

for (response in responses) {
  cat("\n==== Fitting full model for:", response, "====\n")
  
  data_sf$Longitude <- st_coordinates(data_sf)[, 1]
  data_sf$Latitude  <- st_coordinates(data_sf)[, 2]
  data_df <- st_drop_geometry(data_sf)
  
  if (response == "Species_richness") {
    formula_nb <- as.formula(paste(response, "~", predictor, "+ Matern(1 | Longitude + Latitude)"))
    model <- fitme(formula_nb, data = data_df, family = negbin())
    
    beta_table <- as.data.frame(summary(model)$beta_table)
    beta_table$term <- rownames(beta_table)
    beta_table$response <- response
    coef_list[[response]] <- beta_table
  } else {
    vars_to_check <- c(response, predictor, "Longitude", "Latitude")
    data_clean <- data_df %>%
      dplyr::select(all_of(vars_to_check)) %>%
      drop_na()
    
    if (nrow(data_clean) < 5) {
      warning(paste("Skipping", response, "due to too few observations."))
      next
    }
    
    formula_gls <- as.formula(paste(response, "~", predictor))#paste("scale(",response,")", "~", predictor)
    model <- gls(
      formula_gls,
      data = data_clean,
      correlation = corExp(form = ~ Longitude + Latitude),
      method = "ML"
    )
    
    model_summary <- summary(model)
    coef_table <- as.data.frame(model_summary$tTable)
    coef_table$term <- rownames(coef_table)
    coef_table$response <- response
    coef_list[[response]] <- coef_table
  }
}

# Combining and inspecting coefficients
all_coef_main <- bind_rows(coef_list)
print(all_coef_main)
all_coef_main2 <- all_coef_main %>%
  mutate(
    coef_estimate = ifelse(response == "Species_richness", Estimate, Value),
    coef_se = ifelse(response == "Species_richness", `Cond. SE`, Std.Error));all_coef_main2



# Plotting main model coefficients with error bars
ggplot(all_coef_main2 %>% filter(term != "(Intercept)"), 
       aes(x = reorder(response, coef_estimate), y = coef_estimate)) +
  geom_point(size = 1.5, color = "#FC4E07") +
  geom_errorbar(aes(ymin = coef_estimate - 1.96*coef_se, ymax = coef_estimate + 1.96*coef_se), 
                width = 0.3, color = "#00AFBB",size=0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
  coord_flip() +
  labs(
    x = NULL,
    y = "Coefficient Estimate (± SE)",
    title = ""
  ) 
