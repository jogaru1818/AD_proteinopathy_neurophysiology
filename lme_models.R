
### Script wrote by Jonathan Gallego Rudolf - PhD at McGill University

### Run and plot LME models---- 

### Script used to calculate LME models to test for the association between MEG and amyloid
### across the whole-brain, and its interaction with temporal tau accumulation. 
### 1. Get data and load required libraries and functions
### 2. Run LME models
### 3. Calculate model residuals to generate plots
### 4. Main effect plots
### 5. Interaction plots
### 6. Interaction plots for WB Ab*Tau models



### 1. Get data and load required libraries and functions ----

### Load data and required packages
require(ggplot2)  ### Plot visualization 
require(nlme)     ### Run lme models
library(sjPlot)
library(ggpubr)
require(MuMIn)
require(grid)

### Load function to separate residual data into low and high tau for plotting 
### the interactions (based on entorhinal Tau+ 2SD YNG threshold 1.25; based 
### on tau meta-ROI Tau+ 2SD YNG threshold 1.25)
source("/Users/jogar/Documents/PREVENT_AD/Multimodal_data/lme_analysis/scripts/f_separate_tau_data.R")


### Load data
# data=read.csv("C:/Users/jogar/Documents/data_PAD.csv",header=T)

### Verify that categorical variables are defined as factors
data$sex = as.factor(data$sex); data$amyloid_status = as.factor(data$amyloid_status)






### 2. Run LME models----

### Create a function to run LME models and store the results
run_lme_model <- function(outcome, predictors, covariates, filename, path, show_sum,data){
  formula <- as.formula(paste(outcome, " ~ ",predictors," + ", covariates))
  lme_model <- lme(formula, data = data, random = ~1 | ids/rois, control = lmeControl(opt = "optim"), na.action = na.exclude)
  summary_lme_model <- summary(lme_model)
  if(show_sum=="yes"){
    print(summary_lme_model)
  }
  capture.output(summary_lme_model, file = paste(path,filename,".txt"))
}


### Define variables of interest
meg_variables=c("delta","theta","alpha","beta")
pet_variables=c("Ab_data_newr","Ab_data_newr * Tau_ent_newr","Ab_data_newr * Tau_meta_roi_newr",
                "Ab_data_newr * Tau_meta_roi_newr_no_ec","Ab_data_newr * Tau_data_newr")

### Set values to run the model
##### MODIFY OUTCOME, PREDICTORS, COVARIATES AND FILENAME TO RUN DIFFERENT MODELS
outcome=meg_variables[3]
predictors=pet_variables[2]
covariates="age + sex + edu + hipp_vol + trials"; ### ALTERNATIVELY: covariates="age + sex + edu + hipp_vol + trials + offset + exponent"
filename="model_test"
path="C:/Users/jogar/Documents/PREVENT_AD/Multimodal_data/Figures/"
show_sum="yes"

### Run LME
run_lme_model(outcome, predictors, covariates, filename, path, show_sum,data)




### 3. Calculate model residuals to generate plots----

### Define the outcome variables
outcomes <- c("delta", "alpha", "Ab_data_newr", "Tau_data_newr", "Tau_ent_newr", "Tau_meta_roi_newr")

### Create an empty data frame to store residuals
resid_data <- data.frame(); resid_data <- data.frame(data$ids); names(resid_data)[names(resid_data) == 'data.ids'] <- 'ids'

### Loop through the outcome variables to calculate model residuals
for (outcome in outcomes) {
  formula <- as.formula(paste(outcome, "~ age + sex + edu + hipp_vol + trials"))
  resids <- residuals(lme(formula, data=data, random=~1|ids/rois, control=lmeControl(opt="optim"), na.action=na.exclude))
  resid_data[outcome] <- data.frame(resids)
}

### Add extra variables
resid_data$high_tau_ent = data$high_tau_ent;  resid_data$high_tau_meta_roi = data$high_tau_meta_roi;
resid_data$amyloid_status = data$amyloid_status; resid_data$amyloid_tau_status = data$amyloid_tau_status; 
resid_data$tau_ent = data$Tau_ent_newr; resid_data$tau_meta_roi = data$Tau_meta_roi_newr;




### 4. Main effect plots----

### Set path to save figures
path="C:/Users/jogar/Documents/PREVENT_AD/Multimodal_data/Figures/"

main_effects_plot <- function(data, x_var, y_var, group_var, color_var, fill_var, path, save) {
  plot <- ggplot(data = data, aes_string(x = x_var, y = y_var)) +
    geom_smooth(aes_string(group = group_var, color = color_var, fill = fill_var), method = lm, se = FALSE, linetype = "solid", size = 0.2, show.legend = FALSE) +
    scale_color_gradient(low = "black", high = "#F8766D") +
    geom_smooth(color = 'black', fill = 'gray', method = lm, se = TRUE, linetype = "solid", size = 1, show.legend = FALSE) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_rect(fill = "transparent", color = NA), plot.background = element_rect(fill = "transparent", color = NA))
  print(plot)
  if(save=="yes"){
    ggsave(paste0(y_var, "_", x_var, "_main_", color_var, ".png"), plot = plot, units = "in", width = 4, height = 4, dpi = 600, bg = "white", path = path)
  }
}

# Example usage:
main_effects_plot(resid_data, "Ab_data_newr", "delta", "ids", "tau_ent", "tau_ent", path,"no")
main_effects_plot(resid_data, "Ab_data_newr", "delta", "ids", "tau_meta_roi", "tau_meta_roi",path,"no")
main_effects_plot(resid_data, "Ab_data_newr", "alpha", "ids", "tau_ent", "tau_ent",path,"no")
main_effects_plot(resid_data, "Ab_data_newr", "alpha", "ids", "tau_meta_roi", "tau_meta_roi",path,"no")




### 5. Interaction plots----

interaction_plot <- function(data_low_tau, data_high_tau, x_var, y_var,tau_roi,path,save){
  if(tau_roi=="tau_ent"){f_separate_tau_data(resid_data,"entorhinal",1)}
  if(tau_roi=="tau_meta_roi"){f_separate_tau_data(resid_data,"meta_roi",1)}
  if(tau_roi=="wb"){data_low_tau=data_low_tau; data_high_tau=data_high_tau;}
  plot <- ggplot() +
    geom_smooth(data = data_low_tau, aes_string(x = x_var, y = y_var), color = 'black', fill = 'black', method = lm, se = TRUE, linetype = "solid", size = 1.5, show.legend = FALSE) +
    geom_smooth(data = data_high_tau, aes_string(x = x_var, y = y_var), color = '#F8766D', fill = '#F8766D', method = lm, se = TRUE, linetype = "solid", size = 1.5, show.legend = FALSE) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_rect(fill = "transparent", color = NA), plot.background = element_rect(fill = "transparent", color = NA))
  print(plot)
  if(save=="yes"){
    ggsave(paste0(y_var,"_x_",x_var,"_",tau_roi,".png"), plot = plot, units = "in", width = 4, height = 4, dpi = 600, bg = "white", path = path)
  }
}

interaction_plot(data_low_tau, data_high_tau, "Ab_data_newr", "delta", "tau_ent",path, "no")
interaction_plot(data_low_tau, data_high_tau, "Ab_data_newr", "alpha", "tau_ent",path, "no")
interaction_plot(data_low_tau, data_high_tau, "Ab_data_newr", "delta", "tau_meta_roi",path, "no")
interaction_plot(data_low_tau, data_high_tau, "Ab_data_newr", "alpha", "tau_meta_roi",path, "no")