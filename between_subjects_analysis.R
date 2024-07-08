### Script wrote by Jonathan Gallego Rudolf - PhD at McGill University

### Between-subjects analysis for comparing whole brain MEG power and longitudinal cognitive slopes
### across the three PET biomarker-defined groups: (Ab-/Tau-, Ab+/Tau-, Ab+/Tau+

### Perform between-subjects analysis using ANCOVA:

### Run ANCOVAs for each MEG frequency band, then calculate the FDR corrected p values and 
### perform posthoc comparisons to assess significant differences in MEG power across PET groups. 
### Generate boxplots to show the statistical comparisons. 
### Perform between subjects comparison of cognitive slopes (RBANS) based on 
### a) The Ab/Tau group classification and b) The continuous Ab and Tau values
### Repeat analyses removing outliers (based on Cook's distance) to assess robustness of the associations

### 1. Get data and load libraries and functions
### 2. Calculate Ab+ Tau+ thresholds from YNG
### 3. Run ANCOVA models to look at MEG differences between groups
### 4. ANCOVAs for comparing cognition slopes across PET groups




### 1. Get data and load libraries and functions----

require(ggplot2)
library(ggseg)
library(dplyr)
library(tidyr)
library(multcomp)
library(car)
library(see)
library(ggpubr)
library(ggsignif)


library(ggsegExtra)
require(grid)
require(FSA) 
library(tidyverse)
library(rstatix)
library(broom)
library(lmPerm)


##### LOAD data_meg_pet_cogn spreadsheet

### Verify that categorical variables are defined as factors
data_meg_pet_cogn$amyloid_tau_status=as.factor(data_meg_pet_cogn$amyloid_tau_status)
data_meg_pet_cogn$amyloid_status=as.factor(data_meg_pet_cogn$amyloid_status)
data_meg_pet_cogn$high_tau_ent=as.factor(data_meg_pet_cogn$high_tau_ent)


### Remove the Ab-/Tau+ subject
data_meg_pet_cogn2=data_meg_pet_cogn; data_meg_pet_cogn2=data_meg_pet_cogn2[complete.cases(data_meg_pet_cogn2$amyloid_tau_status), ]








### 2. Calculate Ab+ Tau+ thresholds from YNG----

### Calculate Ab+ Tau+ thresholds taking 2SD above the mean of 11 YNG individuals

### Import YNG subjects data
yng_data_ab=read.csv("C:/Users/jogar/Documents/PREVENT_AD/Multimodal_data/PET/newr_PET/YNG_NAV_BL_newr.csv",header=T)
yng_data_tau=read.csv("C:/Users/jogar/Documents/PREVENT_AD/Multimodal_data/PET/newr_PET/YNG_TAU_BL_newr.csv",header=T)

### Calculate Ab threshold
yng_ab_idx=yng_data_ab[,70]; yng_mean_ab_idx=mean(yng_ab_idx); yng_sd_ab_idx=sd(yng_ab_idx);
ab_thres_yng=yng_mean_ab_idx+(3*yng_sd_ab_idx)

### Calculate entorhinal tau threshold
yng_tau_ent=yng_data_tau[,c(6,40)]; yng_tau_ent=rowMeans(yng_tau_ent); yng_mean_tau_ent=mean(yng_tau_ent); yng_sd_tau_ent=sd(yng_tau_ent);
tau_ent_thres_yng=yng_mean_tau_ent+(2*yng_sd_tau_ent)

### Calculate meta-ROI tau threshold
yng_tau_meta_roi=yng_data_tau[,c(6,7,9,13,15,16,40,41,43,47,49,50)]; yng_tau_meta_roi=rowMeans(yng_tau_meta_roi); yng_mean_tau_meta_roi=mean(yng_tau_meta_roi); yng_sd_tau_meta_roi=sd(yng_tau_meta_roi);
tau_meta_roi_thres_yng=yng_mean_tau_meta_roi+(2*yng_sd_tau_meta_roi)







### 3. Run ANCOVA models to look at MEG differences between groups----

### Personalize theme for plots
My_Theme = theme(plot.title = element_text(hjust = 0.5, size = 16, family = "TT Arial"), axis.text.x = element_text(family = "TT Arial"), 
                 axis.text.y = element_text(family = "TT Arial"), axis.title.x = element_text(family = "TT Arial"), 
                 axis.title.y = element_text(family = "TT Arial"), legend.text=element_text(family = "TT Arial"),
                 legend.title =element_text(family = "TT Arial"))


### Create function to run ANCOVA models and generate plots
run_ANCOVA <- function(meg_power, predictors, data, out, save, theme, band, plot_params) {
  
  ### Run model either with or without removing the outliers
  formula <- as.formula(paste(meg_power, "~", paste(predictors, collapse = "+")))
  model <- aov(formula, data = data)
  if(out=="yes"){
    mod=model; cooksD <- cooks.distance(mod)
    influential <- cooksD[(cooksD > (3 * mean(cooksD, na.rm = TRUE)))]; influential;
    names_of_influential <- names(influential);
    outliers <- data[names_of_influential,];
    data_meg_pet_cogn_wo_out <- suppressMessages(data %>% anti_join(outliers));
    formula <- as.formula(paste(meg_power, "~", paste(predictors, collapse = "+")))
    model=aov(formula, data=data_meg_pet_cogn_wo_out) 
  }
  summary_model <- summary(model); print(summary_model)
  
  ### Check model assumptions
  shap=shapiro.test(residuals(model)); print(shap)
  if(out=="no"){formula <- as.formula(paste(meg_power, "~", predictors[1])); lev=leveneTest(formula, data=data); print(lev)} 
  if(out=="yes"){formula <- as.formula(paste(meg_power, "~", predictors[1])); lev=leveneTest(formula, data=data_meg_pet_cogn_wo_out); print(lev)} 
  if(save=="yes" && out=="no"){
    capture.output(summary_model, file = paste(meg_power,"_groups_ancova.txt", sep=""))
    capture.output(shap, file = paste(meg_power,"_ancova_shapiro.txt", sep=""))
    capture.output(lev, file = paste(meg_power,"_ancova_levene.txt", sep=""))
  }
  if(save=="yes" && out=="yes"){
    capture.output(summary_model, file = paste(meg_power,"_groups_ancova_wo.txt", sep=""))
    capture.output(shap, file = paste(meg_power,"_ancova_shapiro_wo.txt", sep=""))
    capture.output(lev, file = paste(meg_power,"_ancova_levene_wo.txt", sep=""))
  }
  
  ### Extract p values
  p_value_main <- summary_model[[1]][["Pr(>F)"]][1]
  posthoc <- TukeyHSD(model, "amyloid_tau_status")
  posthoc_p <- c(posthoc$amyloid_tau_status[10], posthoc$amyloid_tau_status[11], posthoc$amyloid_tau_status[12])
  p1=substring(as.character(posthoc_p[1]),1,5); if(posthoc_p[1]<0.06){p1=paste("p adj =", p1)}; if(posthoc_p[1]>0.06){p1="ns"};
  p2=substring(as.character(posthoc_p[2]),1,5); if(posthoc_p[2]<0.06){p2=paste("p adj =", p2)}; if(posthoc_p[2]>0.06){p2="ns"};
  p3=substring(as.character(posthoc_p[3]),1,5); if(posthoc_p[3]<0.06){p3=paste("p adj =", p3)}; if(posthoc_p[3]>0.06){p3="ns"};
  
  
  ### Generate plot
  if(out=="no"){data_plot=data}
  if(out=="yes"){data_plot=data_meg_pet_cogn_wo_out}
  if(band=="delta"){y="mean_meg_delta"}; if(band=="theta"){y="mean_meg_theta"}; if(band=="alpha"){y="mean_meg_alpha"}; if(band=="beta"){y="mean_meg_beta"};
  p <- ggplot(data_plot, aes(x=amyloid_tau_status, y=.data[[y]], fill=amyloid_tau_status)) + geom_violin(trim=FALSE) + scale_fill_manual(values = c("#00B0F6","#00BA38","#F8766D")) +
    geom_boxplot(width=0.3, outlier.shape=NA) + geom_jitter(shape=16, position=position_jitter(0.1), size=1) +
    labs(y=paste("MEG",  band, "power", sep=" "), x = "Group") + My_Theme + theme_modern() + theme(legend.position="none") + scale_y_continuous(limits = c(plot_params[1], plot_params[2]), breaks = seq(plot_params[3], plot_params[4], by = plot_params[5]), labels = scales::number_format(accuracy = 0.01)) + 
    geom_signif(comparisons = list(c("Ab-/Tau-", "Ab+/Tau-")), annotations = p1, y_position = plot_params[6], vjust = plot_params[9], textsize = 3) +
    geom_signif(comparisons = list(c("Ab-/Tau-", "Ab+/Tau+")), annotations = p2, y_position = plot_params[7], vjust = plot_params[9], textsize = 3) +
    geom_signif(comparisons = list(c("Ab+/Tau-", "Ab+/Tau+")), annotations = p3, y_position = plot_params[8], vjust = plot_params[9], textsize = 3) 
  plot(p)
  if(save=="yes" && out=="no"){ggsave(paste(band, "_amyloid_tau_status_violin.png", sep=""),plot = p,units="in", width=4, height=4, dpi=600,bg="white",path=path)}
  if(save=="yes" && out=="yes"){ggsave(paste(band, "_amyloid_tau_status_violin_wo.png", sep=""),plot = p,units="in", width=4, height=4, dpi=600,bg="white",path=path)}
  return(c(p_value_main = p_value_main, posthoc_p = posthoc_p))
}


### Run ANCOVA models
### Define predictors
predictors=c("amyloid_tau_status", "age", "sex", "education", "hippocampal_volume", "MEG_trials")

### Set parameters for each plot
plot_params_delta=c(-0.1, 0.9, 0, 0.8, 0.2, 0.75,  0.85, 0.8, -0.05)
plot_params_theta=c(0, 0.5, 0.05, 0.45, 0.1, 0.425, 0.475, 0.45,  -0.05)
plot_params_alpha=c(-0.1, 0.9, 0, 0.8, 0.2, 0.75, 0.85, 0.8, -0.05)
plot_params_beta=c(-0.05, 0.45, 0, 0.4, 0.1, 0.375, 0.425, 0.4, -0.05)

### Run models
delta_results <- suppressWarnings(run_ANCOVA("mean_meg_delta", predictors, data_meg_pet_cogn2, "no", "no", My_Theme, "delta", plot_params_delta))
theta_results <- suppressWarnings(run_ANCOVA("mean_meg_theta", predictors, data_meg_pet_cogn2, "no", "no", My_Theme, "theta", plot_params_theta))
alpha_results <- suppressWarnings(run_ANCOVA("mean_meg_alpha", predictors, data_meg_pet_cogn2, "no", "no", My_Theme, "alpha", plot_params_alpha))
beta_results <- suppressWarnings(run_ANCOVA("mean_meg_beta", predictors, data_meg_pet_cogn2, "no", "no", My_Theme, "beta", plot_params_beta))

### Gather p-values and apply multiple comparison correction
p_vals=c(delta_results[1],theta_results[1],alpha_results[1],beta_results[1]); p_vals=p.adjust(p_vals,method = "fdr")
post_hoc_p_vals=cbind(c(delta_results[2],delta_results[3],delta_results[4]), c(theta_results[2],theta_results[3],theta_results[4]),
                      c(alpha_results[2],alpha_results[3],alpha_results[4]), c(beta_results[2],beta_results[3],beta_results[4]))
colnames(post_hoc_p_vals)=c("delta", "theta", "alpha", "beta")


### Run models without outliers
delta_results_wo <- suppressWarnings(run_ANCOVA("mean_meg_delta", predictors, data_meg_pet_cogn2, "yes", "no", My_Theme, "delta", plot_params_delta))
theta_results_wo <- suppressWarnings(run_ANCOVA("mean_meg_theta", predictors, data_meg_pet_cogn2, "yes", "no", My_Theme, "theta", plot_params_theta))
alpha_results_wo <- suppressWarnings(run_ANCOVA("mean_meg_alpha", predictors, data_meg_pet_cogn2, "yes", "no", My_Theme, "alpha", plot_params_alpha))
beta_results_wo <- suppressWarnings(run_ANCOVA("mean_meg_beta", predictors, data_meg_pet_cogn2, "yes", "no", My_Theme, "beta", plot_params_beta))

### Gather p-values and apply multiple comparison correction (without outliers)
p_vals_wo=c(delta_results_wo[1],theta_results_wo[1],alpha_results_wo[1],beta_results_wo[1]); p_vals_wo=p.adjust(p_vals_wo,method = "fdr")
post_hoc_p_vals_wo=cbind(c(delta_results_wo[2],delta_results_wo[3],delta_results_wo[4]), c(theta_results_wo[2],theta_results_wo[3],theta_results_wo[4]),
                         c(alpha_results_wo[2],alpha_results_wo[3],alpha_results_wo[4]), c(beta_results_wo[2],beta_results_wo[3],beta_results_wo[4]))
colnames(post_hoc_p_vals_wo)=c("delta", "theta", "alpha", "beta")








### 4. ANCOVAs for comparing cognition slopes across PET groups ----

run_ANCOVA_cognition <- function(cognition, predictors, data, out, save, theme, plot_params) {
  
  ### Run model either with or without removing the outliers
  formula <- as.formula(paste(cognition, "~", paste(predictors, collapse = "+")))
  model <- aov(formula, data = data)
  if(out=="yes"){
    mod=model; cooksD <- cooks.distance(mod)
    influential <- cooksD[(cooksD > (3 * mean(cooksD, na.rm = TRUE)))]; influential;
    names_of_influential <- names(influential);
    outliers <- data[names_of_influential,];
    data_meg_pet_cogn_wo_out <- suppressMessages(data %>% anti_join(outliers));
    formula <- as.formula(paste(cognition, "~", paste(predictors, collapse = "+")))
    model=aov(formula, data=data_meg_pet_cogn_wo_out) 
  }
  summary_model <- summary(model); print(summary_model)
  
  ### Check model assumptions
  shap=shapiro.test(residuals(model)); print(shap)
  if(out=="no"){formula <- as.formula(paste(cognition, "~", predictors[1])); lev=leveneTest(formula, data=data); print(lev)} 
  if(out=="yes"){formula <- as.formula(paste(cognition, "~", predictors[1])); lev=leveneTest(formula, data=data_meg_pet_cogn_wo_out); print(lev)} 
  if(save=="yes" && out=="no"){
    capture.output(summary_model, file = paste(cognition,"_groups_ancova.txt", sep=""))
    capture.output(shap, file = paste(cognition,"_ancova_shapiro.txt", sep=""))
    capture.output(lev, file = paste(cognition,"_ancova_levene.txt", sep=""))
  }
  if(save=="yes" && out=="yes"){
    capture.output(summary_model, file = paste(cognition,"_groups_ancova_wo.txt", sep=""))
    capture.output(shap, file = paste(cognition,"_ancova_shapiro_wo.txt", sep=""))
    capture.output(lev, file = paste(cognition,"_ancova_levene_wo.txt", sep=""))
  }
  
  ### Extract p values
  p_value_main <- summary_model[[1]][["Pr(>F)"]][1]
  posthoc <- TukeyHSD(model, "amyloid_tau_status")
  posthoc_p <- c(posthoc$amyloid_tau_status[10], posthoc$amyloid_tau_status[11], posthoc$amyloid_tau_status[12])
  p1=substring(as.character(posthoc_p[1]),1,6); if(posthoc_p[1]<0.06){p1=paste("p adj =", p1)}; if(posthoc_p[1]>0.06){p1="ns"};
  p2=substring(as.character(posthoc_p[2]),1,6); if(posthoc_p[2]<0.06){p2=paste("p adj =", p2)}; if(posthoc_p[2]>0.06){p2="ns"};
  p3=substring(as.character(posthoc_p[3]),1,6); if(posthoc_p[3]<0.06){p3=paste("p adj =", p3)}; if(posthoc_p[3]>0.06){p3="ns"};
  
  
  ### Generate plot
  if(out=="no"){data_plot=data}
  if(out=="yes"){data_plot=data_meg_pet_cogn_wo_out}
  if(cognition=="attention"){y="attention"}; if(cognition=="immediate_memory"){y="immediate_memory"}; if(cognition=="delayed_memory"){y="delayed_memory"}; 
  if(cognition=="language"){y="language"}; if(cognition=="visuospatial_construction"){y="visuospatial_construction"}; if(cognition=="total"){y="total"}; if(cognition=="att_mem_score"){y="att_mem_score"};
  p <- ggplot(data_plot, aes(x=amyloid_tau_status, y=.data[[y]], fill=amyloid_tau_status)) + geom_violin(trim=FALSE) + scale_fill_manual(values = c("#00B0F6","#00BA38","#F8766D")) +
    geom_boxplot(width=0.3, outlier.shape=NA) + geom_jitter(shape=16, position=position_jitter(0.1), size=1) +
    labs(y=paste("MEG",  cognition, "power", sep=" "), x = "Group") + My_Theme + theme_modern() + theme(legend.position="none") + scale_y_continuous(limits = c(plot_params[1], plot_params[2]), breaks = seq(plot_params[3], plot_params[4], by = plot_params[5]), labels = scales::number_format(accuracy = 0.01)) + 
    geom_signif(comparisons = list(c("Ab-/Tau-", "Ab+/Tau-")), annotations = p1, y_position = plot_params[6], vjust = plot_params[9], textsize = 3) +
    geom_signif(comparisons = list(c("Ab-/Tau-", "Ab+/Tau+")), annotations = p2, y_position = plot_params[7], vjust = plot_params[9], textsize = 3) +
    geom_signif(comparisons = list(c("Ab+/Tau-", "Ab+/Tau+")), annotations = p3, y_position = plot_params[8], vjust = plot_params[9], textsize = 3) 
  plot(p)
  if(save=="yes" && out=="no"){ggsave(paste(cognition, "_amyloid_tau_status_violin.png", sep=""),plot = p,units="in", width=4, height=4, dpi=600,bg="white",path=path)}
  if(save=="yes" && out=="yes"){ggsave(paste(cognition, "_amyloid_tau_status_violin_wo.png", sep=""),plot = p,units="in", width=4, height=4, dpi=600,bg="white",path=path)}
  return(c(p_value_main = p_value_main, posthoc_p = posthoc_p))
}



### Run models

### Define predictors
predictors=c("amyloid_tau_status", "age", "sex", "education", "hippocampal_volume", "MEG_trials", "days_to_MEG")

### Set parameters for each plot
plot_params_attention=c(-16, 14, -14, 10, 6, 10,  13, 11.5, -0.05)
plot_params_imm_mem=c(-14, 17, -12, 12, 6, 13, 16, 14.5,  -0.05)
plot_params_del_mem=c(-14, 17, -12, 12, 6, 13, 16, 14.5, -0.05)
# plot_params_language=c(-14, 17, -12, 12, 6, 13, 16, 14.5, -0.05)
# plot_params_visuo=c(-14, 17, -12, 12, 6, 13, 16, 14.5, -0.05)
plot_params_total=c(-12, 12, -12, 8, 4, 8, 11, 9.5, -0.05)
plot_params_att_mem=c(-8, 13, -8, 12, 4, 9, 12, 10.5, -0.05) 

### Set NA values to subject with one cognitive timepoint 
data_meg_pet_cogn3=data_meg_pet_cogn2[complete.cases(data_meg_pet_cogn2$total), ]


### Run cognition ANCOVA models
attention_results <- suppressWarnings(run_ANCOVA_cognition("attention", predictors, data_meg_pet_cogn3, "no", "no", My_Theme, plot_params_attention))
imm_mem_results <- suppressWarnings(run_ANCOVA_cognition("immediate_memory", predictors, data_meg_pet_cogn3, "no", "no", My_Theme, plot_params_imm_mem))
del_mem_results <- suppressWarnings(run_ANCOVA_cognition("delayed_memory", predictors, data_meg_pet_cogn3, "no", "no", My_Theme, plot_params_del_mem))
# language_results <- suppressWarnings(run_ANCOVA_cognition("language", predictors, data_meg_pet_cogn3, "no", "no", My_Theme, plot_params_language))
# visuo_results <- suppressWarnings(run_ANCOVA_cognition("visuospatial_construction", predictors, data_meg_pet_cogn3, "no", "no", My_Theme, plot_params_visuo))
total_results <- suppressWarnings(run_ANCOVA_cognition("total", predictors, data_meg_pet_cogn3, "no", "no", My_Theme, plot_params_total))
att_mem_results <- suppressWarnings(run_ANCOVA_cognition("att_mem_score", predictors, data_meg_pet_cogn3, "no", "no", My_Theme, plot_params_att_mem))

### Gather p-values and apply multiple comparison correction (only ran for 4 tests - attention, imm_mem, del_mem and att_mem_score)
p_vals=c(attention_results[1],imm_mem_results[1],del_mem_results[1],att_mem_results[1]); p_vals=p.adjust(p_vals,method = "fdr")
post_hoc_p_vals=cbind(c(attention_results[2],attention_results[3],attention_results[4]), c(imm_mem_results[2],imm_mem_results[3],imm_mem_results[4]),
                      c(del_mem_results[2],del_mem_results[3],del_mem_results[4]), c(att_mem_results[2],att_mem_results[3],att_mem_results[4]))
colnames(post_hoc_p_vals)=c("attention", "immediate_memory", "delayed_memory", "att_mem_score")



### Run cognition ANCOVA models without outliers
attention_results_wo <- suppressWarnings(run_ANCOVA_cognition("attention", predictors, data_meg_pet_cogn3, "yes", "no", My_Theme, plot_params_attention))
imm_mem_results_wo <- suppressWarnings(run_ANCOVA_cognition("immediate_memory", predictors, data_meg_pet_cogn3, "yes", "no", My_Theme, plot_params_imm_mem))
del_mem_results_wo <- suppressWarnings(run_ANCOVA_cognition("delayed_memory", predictors, data_meg_pet_cogn3, "yes", "no", My_Theme, plot_params_del_mem))
# language_results_wo <- suppressWarnings(run_ANCOVA_cognition("language", predictors, data_meg_pet_cogn3, "yes", "no", My_Theme, plot_params_language))
# visuo_results_wo <- suppressWarnings(run_ANCOVA_cognition("visuospatial_construction", predictors, data_meg_pet_cogn3, "yes", "no", My_Theme, plot_params_visuo))
total_results_wo <- suppressWarnings(run_ANCOVA_cognition("total", predictors, data_meg_pet_cogn3, "yes", "no", My_Theme, plot_params_total))
att_mem_results_wo <- suppressWarnings(run_ANCOVA_cognition("att_mem_score", predictors, data_meg_pet_cogn3, "yes", "no", My_Theme, plot_params_att_mem))

### Gather p-values and apply multiple comparison correction (only ran for 4 tests - attention, imm_mem, del_mem and att_mem_score)
p_vals_wo=c(attention_results_wo[1],imm_mem_results_wo[1],del_mem_results_wo[1],att_mem_results_wo[1]); p_vals=p.adjust(p_vals,method = "fdr")
post_hoc_p_vals_wo=cbind(c(attention_results_wo[2],attention_results_wo[3],attention_results_wo[4]), c(imm_mem_results_wo[2],imm_mem_results_wo[3],imm_mem_results_wo[4]),
                         c(del_mem_results_wo[2],del_mem_results_wo[3],del_mem_results_wo[4]), c(att_mem_results_wo[2],att_mem_results_wo[3],att_mem_results_wo[4]))
colnames(post_hoc_p_vals_wo)=c("attention", "immediate_memory", "delayed_memory", "att_mem_score")



