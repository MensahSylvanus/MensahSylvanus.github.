library(lme4)
library(dplyr)

# Define responses and predictor
responses <- c("Species_richness", "FDis", "sesMPD", "sesMNTD", "cvDBH", "skwDBH", 
               "Tree_density", "LargeTrees_density", "RemTrees_density", "HLorey")
predictor <- "Fire_frequency"

# Store results
anova_results <- list()

for (resp in responses) {
  cat("\n\n--- Testing random effects for:", resp, "---\n")
  
  # Build model formulas
  full_formula      <- as.formula(paste(resp, "~", predictor, "+ (1|Plot_size) + (1|Vegetation)"))
  reduced_plot      <- as.formula(paste(resp, "~", predictor, "+ (1|Vegetation)"))
  reduced_veg       <- as.formula(paste(resp, "~", predictor, "+ (1|Plot_size)"))
  
  # Fit models
  full_mod          <- lmer(full_formula, data = data, REML = FALSE)
  mod_no_plot       <- lmer(reduced_plot, data = data, REML = FALSE)
  mod_no_veg        <- lmer(reduced_veg, data = data, REML = FALSE)
  
  # Model with only Vegetation
  veg_only_formula  <- as.formula(paste(resp, "~", predictor, "+ (1|Vegetation)"))
  plot_only_formula <- as.formula(paste(resp, "~", predictor, "+ (1|Plot_size)"))
  
  mod_veg_only      <- lmer(veg_only_formula, data = data, REML = FALSE)
  mod_plot_only     <- lmer(plot_only_formula, data = data, REML = FALSE)
  
  # Compare models
  comp_no_plot  <- anova(full_mod, mod_no_plot)
  comp_no_veg   <- anova(full_mod, mod_no_veg)
  comp_plot_vs_veg <- anova(mod_plot_only, mod_veg_only)
  
  # Store results
  anova_results[[resp]] <- list(
    no_plot         = comp_no_plot,
    no_veg          = comp_no_veg,
    plot_vs_veg     = comp_plot_vs_veg,
    full_model_fit  = summary(full_mod)
  )
  
  # Print summaries
  cat("\n??? Comparing Full vs No Plot_size:\n")
  print(comp_no_plot)
  
  cat("\n??? Comparing Full vs No Vegetation:\n")
  print(comp_no_veg)
  
  cat("\n??? Comparing Plot_size-only vs Vegetation-only:\n")
  print(comp_plot_vs_veg)
}




# Initialize an empty data frame
summary_stats <- data.frame()

for (resp in names(anova_results)) {
  res <- anova_results[[resp]]
  
  # Extract LRT statistics and p-values
  comp1 <- res$no_plot
  comp2 <- res$no_veg
  comp3 <- res$plot_vs_veg
  
  # Create a row for each comparison
  df_comp <- data.frame(
    response = rep(resp, 3),
    comparison = c("Full vs No Plot_size", "Full vs No Vegetation", "Plot_size vs Vegetation"),
    df = c(comp1$Df[2], comp2$Df[2], comp3$Df[2]),
    AIC_1 = c(comp1$AIC[1], comp2$AIC[1], comp3$AIC[1]),
    AIC_2 = c(comp1$AIC[2], comp2$AIC[2], comp3$AIC[2]),
    LRT_stat = c(comp1$Chisq[2], comp2$Chisq[2], comp3$Chisq[2]),
    p_value = c(comp1$`Pr(>Chisq)`[2], comp2$`Pr(>Chisq)`[2], comp3$`Pr(>Chisq)`[2])
  )
  
  summary_stats <- rbind(summary_stats, df_comp)
}

# Round for readability
summary_stats <- summary_stats %>%
  mutate(
    AIC_diff = AIC_2 - AIC_1,
    p_value = signif(p_value, 3),
    LRT_stat = round(LRT_stat, 2)
  )

# Print summary
print(summary_stats)
