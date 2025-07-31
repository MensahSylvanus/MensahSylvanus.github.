
----------------------------------------------------------------------------------------------
#### testing the effects of fire frequency on aboveground carbon (AGC)
----------------------------------------------------------------------------------------------

folds_list <- sb$folds_list 

# Storing results

results_gls <- list()
fitted_gls_models <- list()
normality_results <- data.frame()
r2_results <- data.frame()


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
  
  ## Fitting Spatial GLS
  gls_model <- gls(
    AGC ~ Fire_frequency,
    correlation = corExp(form = ~ Longitude + Latitude),
    data = train_df,
    method = "ML"
  )
  
  # Saving model
  fitted_gls_models[[i]] <- gls_model
  
  # Predicting
  gls_pred <- predict(gls_model, newdata = test_df)
  observed = test_df$AGC
  test_resid = observed - gls_pred
  
  #storing
  results_gls[[i]] <- data.frame(
    observed = observed,
    predicted = gls_pred,
    test_resid = test_resid,
    fold = i,
    model = "SpatialGLS")
    
    # Shapiro-Wilk normality test on residuals for this fold
    shapiro_test <- shapiro.test(test_resid)
    normality_results <- rbind(normality_results, data.frame(fold = i, shapiro_p = shapiro_test$p.value))
                            
    # Computing R² as squared Pearson correlation
     r2 <- cor(observed, gls_pred, use = "complete.obs")^2
     r2_results <- rbind(r2_results, data.frame(fold = i, R2 = r2))
    
}

for (i in seq_along(fitted_gls_models)) {
  cat("\n==== Fold", i, "====\n")
  print(summary(fitted_gls_models[[i]]))
}

library(ggplot2)
library(dplyr)


# Combining results from all folds into one data.frame
all_results_gls <- bind_rows(results_gls);all_results_gls

# Calculating residuals
all_results_gls <- all_results_gls %>%
  mutate(residual = observed - predicted)

# Plotting residuals vs predicted values
ggplot(all_results_gls, aes(x = predicted, y = residual)) +
  geom_point(size=1.5, color="#00AFBB")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "Residuals vs Predicted AGC (Spatial GLS)",
    x = "Predicted AGC",
    y = "Residuals (Observed - Predicted)"
  ) ->a;a

all_normality_results <- dplyr::bind_rows(normality_results);all_normality_results



# Slope coefficients and pvalues

# Initializing list to hold coefficient summaries
gls_coef_list <- list()

for (i in seq_along(fitted_gls_models)) {
  model <- fitted_gls_models[[i]]
  summary_model <- summary(model)
  coef_table <- as.data.frame(summary_model$tTable)
  
  coef_table$term <- rownames(coef_table)
  coef_table$fold <- i
  
  gls_coef_list[[i]] <- coef_table
}

# Combining all into one data frame
gls_coef_df <- dplyr::bind_rows(gls_coef_list) %>%
  dplyr::select(fold, term, Estimate=Value, Std.Error, t.value = `t-value`, p.value = `p-value`);gls_coef_df



#R2
r2_results
meanR2<-mean(r2_results$R2);meanR2#> [1] 0.0847434


##### Decision tree

data <- read.csv("fireplotdata.csv")
#theme_set(theme_bw())
library(sjPlot)
set_theme(theme_sjplot2(base_size = 12),axis.title.color = "black",axis.textcolor.x = "black",
          axis.tickslen = 0, # hides tick marks
          axis.title.size = 0.9,
          axis.textsize = 0.9,
          legend.size = 1.1,
          legend.title.size = 1.1,plot.margin=unit(c(0,0,0.1,0), "cm"),
          geom.label.size = 3.5)
a<- ggplot(data,  aes(x=Fire_frequency, y=AGC))+geom_point(size=2, color="#00AFBB")+
  geom_smooth(method="lm", size=1,linetype=1, se=T,color = "#FC4E07")+
 #   geom_quantile(quantiles = 0.99, size=1, colour = "black",linetype=1)+
  labs(x="Fire frequency", y= "AGC (Mg/ha)"); a


library(ggparty)
data$`Fire frequency` <- data$Fire_frequency
tr_tree <- lmtree(AGC~ `Fire frequency`, data, caseweights = T)

treeplot1<-ggparty(tr_tree , terminal_space = 0.4, horizontal = F) +
  geom_edge() +
  geom_node_splitvar() +
  geom_edge_label() +
  geom_node_plot(
    gglist = list(geom_boxplot(aes(y=AGC, x=""),position=position_dodge(), color = "#FC4E07",fill="lightgrey",outlier.shape = NA),
                  geom_jitter(aes(y = AGC, x = ""), width = 0.1, size = 1.5, color = "#00AFBB", alpha = 0.8),
                  theme_bw()+theme(plot.margin = margin(t = -40, r = 10, b = -40, l = 0)),
                  ylab("AGC (Mg/ha)"), xlab("")),
    ids = "terminal", # not necessary since default
    shared_axis_labels = T);treeplot1


plot_grid(a,treeplot1, ncol=2, labels =c("a","b"), align="hv")

AGC_fire_mod0<-gls(AGC~Fire_frequency, data)
AGC_fire_mod_ac<-update(AGC_fire_mod0,correlation = corExp(form = ~ Longitude+Latitude))
summary(AGC_fire_mod_ac)

extract_model_stats <- function(model) {
  r2 <- mean(r2_results$R2)
  p <- summary(model)$tTable["Fire_frequency", "p-value"]
  label <- paste0("R² = ", round(r2 * 100, 2), "%\np = ", format.pval(p, digits = 1, eps = .001))
  return(label)
}

label_text <- extract_model_stats(AGC_fire_mod_ac);label_text

grid.text(label_text, x=0.40, y=0.90, gp=gpar(cex = 0.8, font = 4, col = "darkblue"))# Fig 2



#--------------------------------------------
### Exploring correlation between  variables
#-------------------------------------------

data <- read.csv("fireplotdata.csv")
names(data)
selected<-c("Fire_frequency","rainfall" ,"rainfall_seasonality" ,"temperature" ,"Total_N" ,"clay_silt" ,  "Species_richness", "FDis" ,"sesMNTD", "sesMPD","cvDBH", "skwDBH",
            "Tree_density","LargeTrees_density","RemTrees_density", "HLorey" )
selected<-data.frame(data[selected])
names(selected)<-c("Fire_freq","rainfall" ,"rainfall_seas" ,"temp" , "Total_N" ,"clay_silt" ,  "Spe_rich", "FDis" ,"sesMNTD", "sesMPD","cvDBH", "skwDBH",
                   "Tree_density","LT_density","RemT_density", "HLorey" )
library(GGally)
ggcorr(selected,high = "#3B9AB2",
       mid = "#EEEEEE",
       low = "#F21A00",
       label = TRUE,label_round = 2,legend.position = "",
       label_size = 5,
       color = "black")### Fig S6

dev.off()




#################################################################################################################################################
#########   Based on correlation matrix between the predictors , we excluded overall tree density and density of                      
#########   95% remaining trees because they were highly correlated (r=1), and were also strongly correlated with species richness (r >0.8).
#########   In addition, because sesMPD and sesMNTD were also strongly correlated, we excluded one variable (sesMPD).                           
#########   For climate, we used rainfall seasonality only because it stood as the most important variable driving fire frequency               
#########   With the remaining candidate variables, we performed a multiple model to determine 
########     how AGC related with fire, climate, soil (physical and chemical properties), species diversity 
########     (phylogenetic, taxonomic and functional trait attributes) and stand structural attributes
################################################################################################################################################

###  Performing the MuMIn with AGC as the response and RVI

data <- read.csv("fireplotdata.csv")
nrow(data)
data<-na.omit(data);nrow(data);# removing plots (4) with low species richness, where sesMNTD could not be computed

all_variables<-names(data[c(6:22)])
# Standardizing predictors
data <- data %>%
  mutate(across(all_of(all_variables), ~scale(.)[,1]))
names(data)

AGC_mod<-lme(AGC~Fire_frequency+rainfall_seasonality+clay_silt+Total_N+Species_richness+FDis+sesMNTD+cvDBH+skwDBH+LargeTrees_density,
             random= list (~1|Vegetation, ~1|Plot_size), data)
AGC_mod.ac<-update(AGC_mod,correlation = corExp(form = ~ Longitude+Latitude))
AGC_mod.resid<-as.data.frame(AGC_mod$residuals)
AGC_mod.ac.resid<-as.data.frame(AGC_mod.ac$residuals)

Longlatdist <- as.matrix(dist(cbind(data$Longitude, data$Latitude)))
Longlatdist.inv <- 1/Longlatdist
diag(Longlatdist.inv) <- 0

Moran.I(AGC_mod.resid$fixed, Longlatdist.inv)
Moran.I(AGC_mod.ac.resid$fixed, Longlatdist.inv)
shapiro.test(AGC_mod.resid$fixed)
shapiro.test(AGC_mod.ac.resid$fixed)
anova(AGC_mod,AGC_mod.ac)# Table S4
summary(AGC_mod) # Table S5
#summary(AGC_mod.ac)
cbind(car::vif(AGC_mod))
r.squaredGLMM(AGC_mod)

library(MuMIn)
AGC_mod_mumin <- MuMIn::dredge(AGC_mod,beta="sd",trace = T, evaluate = TRUE, rank = "AICc")

summary(model.avg(AGC_mod_mumin, fit = TRUE))
final.mod<-model.avg(AGC_mod_mumin, fit = TRUE, delta<2)
summary(final.mod)  # Table S6 
final.mod<-get.models(AGC_mod_mumin, 3)[[1]]
summary(final.mod)
r.squaredGLMM(final.mod)
cbind(car::vif(final.mod))


### Plotting RVI (Figure 4)
cbind(sw(AGC_mod_mumin))
db.RVI <- data.frame(cbind(sw(AGC_mod_mumin)))

db.RVI$Variables<-row.names(db.RVI)
db.RVI$Model <-  c("Structure",'Diversity','Diversity',"Climate","Disturbance",'Diversity','Diversity',"Soil","Soil",'Diversity')
colnames(db.RVI)<-c("RVI","Variables","Model")

db.RVI$Variables <- factor(db.RVI$Variables,levels = c
                           ("rainfall_seasonality","Fire_frequency",
                             "clay_silt","Total_N",
                             "Species_richness","sesMNTD", 'FDis',"skwDBH","cvDBH",
                             "LargeTrees_density"), ordered=F)
db.RVI$Model <- factor(db.RVI$Model,levels = 
                         c("Climate","Disturbance","Soil", 'Diversity',"Structure"))

db.RVI$rel.sum<-db.RVI$RVI/sum(db.RVI$RVI)

RVI_plot<-ggplot(data = db.RVI, aes(x=Variables, y=rel.sum)) + 
  geom_bar(width=0.8,stat="identity",position=position_dodge(width=0.9),aes(fill=Variables))+
  xlab("")+ylab("Relative sum of weights")+ 
  #geom_text(aes(label=round(RVI,digits = 2)), hjust=-0.2, angle=90)+
  theme(legend.title = element_blank())+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  theme(legend.key.size = unit(0.55, 'cm'),legend.text= element_text(size = 9));RVI_plot





#----------------------------------------------------------------------------------------------
#  #### testing the effects of fire frequency, climate, soil diversity and structures  on aboveground carbon (AGC)
#----------------------------------------------------------------------------------------------

folds_list <- sb$folds_list 

# Storing results

results_gls <- list()
fitted_gls_models <- list()
normality_results <- data.frame()
r2_results <- data.frame()


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
  
  predictor_vars <- c("Fire_frequency", "rainfall_seasonality", "clay_silt", "Total_N", 
                      "sesMNTD", "FDis", "Species_richness", "cvDBH", "skwDBH", 
                      "LargeTrees_density", "X", "Y", "AGC")
  train_df_clean <- train_df %>%
    dplyr::select(all_of(predictor_vars)) %>%
    drop_na()
  
  test_df_clean <- test_df %>%
    dplyr::select(all_of(predictor_vars)) %>%
    drop_na()
  
  ## Fitting Spatial GLS
  gls_model <- gls(
    AGC ~ sesMNTD+Species_richness+LargeTrees_density,
    correlation = corExp(form = ~ X + Y),
    data = train_df_clean,
    method = "ML"
  )
  
  # Saving model
  fitted_gls_models[[i]] <- gls_model
  
  # Predicting
  gls_pred <- predict(gls_model, newdata = test_df_clean)
  observed = test_df_clean$AGC
  test_resid = observed - gls_pred
  
  #storing
  results_gls[[i]] <- data.frame(
    observed = observed,
    predicted = gls_pred,
    test_resid = test_resid,
    fold = i,
    model = "SpatialGLS")
  
  # Shapiro-Wilk normality test on residuals for this fold
  shapiro_test <- shapiro.test(test_resid)
  normality_results <- rbind(normality_results, data.frame(fold = i, shapiro_p = shapiro_test$p.value))
  
  # Computing R² as squared Pearson correlation
  r2 <- cor(observed, gls_pred, use = "complete.obs")^2
  r2_results <- rbind(r2_results, data.frame(fold = i, R2 = r2))
  
}

for (i in seq_along(fitted_gls_models)) {
  cat("\n==== Fold", i, "====\n")
  print(summary(fitted_gls_models[[i]]))
}

library(ggplot2)
library(dplyr)


# Combining results from all folds into one data.frame
all_results_gls <- bind_rows(results_gls);all_results_gls

# Calculating residuals
all_results_gls <- all_results_gls %>%
  mutate(residual = observed - predicted)

# Plotting residuals vs predicted values
ggplot(all_results_gls, aes(x = predicted, y = residual)) +
  geom_point(size=1.5, color="#00AFBB")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "Residuals vs Predicted AGC (Spatial GLS)",
    x = "Predicted AGC",
    y = "Residuals (Observed - Predicted)"
  ) ->a;a

all_normality_results <- dplyr::bind_rows(normality_results);all_normality_results



# Slope coefficients and pvalues

# Initializing list to hold coefficient summaries
gls_coef_list <- list()

for (i in seq_along(fitted_gls_models)) {
  model <- fitted_gls_models[[i]]
  summary_model <- summary(model)
  coef_table <- as.data.frame(summary_model$tTable)
  
  coef_table$term <- rownames(coef_table)
  coef_table$fold <- i
  
  gls_coef_list[[i]] <- coef_table
}

# Combining all into one data frame
gls_coef_df <- dplyr::bind_rows(gls_coef_list) %>%
  dplyr::select(fold, term, Estimate=Value, Std.Error, t.value = `t-value`, p.value = `p-value`);gls_coef_df



#R2
r2_results
meanR2<-mean(r2_results$R2);meanR2#> 
  
  
  
  
#----------------------------------------------------------------------------------------------
######## structural equation model (SEM)
#----------------------------------------------------------------------------------------------



#### Piecewise SEM

data <- read.csv("fireplotdata.csv")
data<-na.omit(data)# removing plots (4) with low species richness, where sesMNTD could not be computed
data$LargeTrees_densty<-data$LargeTrees_density*data$Plot_size# convertng LargeTrees_density back to absolute density (to remove plot size effect) because plot size used as random later 

var_to_scale<-c("LargeTrees_densty", "rainfall_seasonality","sesMNTD","AGC")

# Standardizing predictors
data <- data %>%mutate(across(all_of(var_to_scale), ~scale(.)[,1]))

m0<-glm.nb(Fire_frequency~rainfall_seasonality, data)
m1<-gls(LargeTrees_densty~Fire_frequency+rainfall_seasonality,
        correlation = corExp(form = ~ Longitude + Latitude),
        data = data,
        method = "ML")
m2<-glm.nb(Species_richness~Fire_frequency+rainfall_seasonality, data)
m3<-gls(sesMNTD~Fire_frequency+rainfall_seasonality,
        correlation = corExp(form = ~ Longitude + Latitude),
        data = data,
        method = "ML")
m4<-gls(AGC ~ Fire_frequency+rainfall_seasonality+LargeTrees_densty+Species_richness+sesMNTD,
  correlation = corExp(form = ~ Longitude + Latitude),
  data = data,
  method = "ML")

library(piecewiseSEM)

sem0 <- psem(m0,m1,m2,m3,m4)
summary(sem0,conserve = T, .progressBar = F,data=data)
sem <- update(sem0, sesMNTD %~~% Species_richness)
summary(sem,conserve = T, .progressBar = F)
sem_summary <- summary(sem, conserve = TRUE, .progressBar = FALSE)

# Extracting path coefficients table
coefs <- sem_summary$coefficients

dir.create("sem_outputs")
# Saving to CSV
write.csv(coefs, "sem_outputs/psem_coefficients.csv", row.names = FALSE)



### Lavaan SEM

data <- read.csv("fireplotdata.csv")
data$LargeTrees_densty<-data$LargeTrees_density*data$Plot_size# convertng LargeTrees_density back to absolute density (to remove plot size effect) because plot size used as random later 
data$Fire_frequency_t<-log(1+data$Fire_frequency)
data$Species_richness_t<-log(1+data$Species_richness)
var_to_scale<-c("LargeTrees_densty", "Species_richness_t","Fire_frequency_t", "rainfall_seasonality","sesMNTD","AGC")

# Standardizing predictors
data <- data %>%mutate(across(all_of(var_to_scale), ~scale(.)[,1]))

fire.sem <- '
Fire_frequency_t~a*rainfall_seasonality
LargeTrees_densty~b*Fire_frequency_t+c*rainfall_seasonality
Species_richness_t~d*Fire_frequency_t+e*rainfall_seasonality
sesMNTD~f*Fire_frequency_t+g*rainfall_seasonality
AGC ~ h*Fire_frequency_t+i*rainfall_seasonality+j*LargeTrees_densty+k*Species_richness_t+l*sesMNTD
Species_richness_t ~~ sesMNTD
# indirect effect of fire on AGC
ind.fire.on.AGC:=(b*j)
tot.fire.on.AGC:=(b*j)+h

# indirect effect of rain_seas. on Large trees
ind.rain.on.LargeTrees_densty:=(a*b)
tot.rain.on.LargeTrees_densty:=(a*b)+c

# indirect effect of rain_seas. on AGC
ind.rain.on.AGC via fire:=(a*b*j)
tot.rain.on.AGC via LTD:=c*j
'

library(lavaan)
fit.fire.sem <- sem(fire.sem, data=data,check.gradient = FALSE) ## Fitting the model with sem functions
summary(fit.fire.sem, modindices=TRUE)
summary(fit.fire.sem, rsq=TRUE, standardized=TRUE, fit.measures=TRUE)
fitMeasures(fit.fire.sem, c("chisq", "df", "pvalue", "cfi", "rmsea","gfi","srmr", "AIC"))
lsem<-standardizedSolution(fit.fire.sem, type="std.all")
lsem
# Saving to CSV
write.csv(lsem, "sem_outputs/lsem_coefficients.csv", row.names = FALSE)

