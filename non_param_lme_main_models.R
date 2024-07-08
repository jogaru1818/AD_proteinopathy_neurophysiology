
### Script wrote by Jonathan Gallego Rudolf - PhD at McGill University

### Non-parametric hypotheisis testing on main LME models---- 

### Script used to run non-parametric asseessment of the robustness for the main effect
### of Ab and the Ab x Tau interaction. The script takes the Ab data and shuffles the
### order of the regions for each individual, while keeping their MEG and 

### 1. Create a function to run LME models and store the results
### 2. Shuffle the data and re run LME models
### 3. Save data from null LME models
### 4. Compute p values




### 1. Create a function to run LME models and store the results----
run_lme_model_perm <- function(outcome, predictors, out_coeff, data) {
  formula <- as.formula(paste(outcome, " ~ ",predictors," + age + sex + edu + hipp_vol + trials"))
  lme_model <- lme(formula, data = data, random = ~1 | ids/rois, control = lmeControl(opt = "optim"), na.action = na.exclude)
  summary_lme_model <- summary(lme_model)
  if(out_coeff=="main_eff"){
    t_ab <- coef(summary_lme_model)[2, "t-value"];p_ab <- coef(summary_lme_model)[2, "p-value"]
    return(data.frame(t_ab, p_ab))
  }
  if(out_coeff=="int"){
    t_ab <- coef(summary_lme_model)[2, "t-value"];p_ab <- coef(summary_lme_model)[2, "p-value"]
    t_tau <- coef(summary_lme_model)[3, "t-value"];p_tau <- coef(summary_lme_model)[3, "p-value"]
    t_ab_tau <- coef(summary_lme_model)[9, "t-value"];p_ab_tau <- coef(summary_lme_model)[9, "p-value"]
    return(data.frame(t_ab, p_ab, t_tau, p_tau, t_ab_tau, p_ab_tau))
  }   
  if(out_coeff=="int_aper"){
    t_ab <- coef(summary_lme_model)[2, "t-value"];p_ab <- coef(summary_lme_model)[2, "p-value"]
    t_tau <- coef(summary_lme_model)[3, "t-value"];p_tau <- coef(summary_lme_model)[3, "p-value"]
    t_ab_tau <- coef(summary_lme_model)[11, "t-value"];p_ab_tau <- coef(summary_lme_model)[11, "p-value"]
    return(data.frame(t_ab, p_ab, t_tau, p_tau, t_ab_tau, p_ab_tau))
  }
}






### 2. Shuffle the data and re run LME models---- 

### Initialize variables for all models
stats_table_delta_ab=stats_table_theta_ab=stats_table_alpha_ab=stats_table_beta_ab=data.frame()
stats_table_delta_ab_etau=stats_table_theta_ab_etau=stats_table_alpha_ab_etau=stats_table_beta_ab_etau=data.frame()
stats_table_delta_ab_tau_mr=stats_table_theta_ab_tau_mr=stats_table_alpha_ab_tau_mr=stats_table_beta_ab_tau_mr=data.frame()
stats_table_delta_ab_tau_wb=stats_table_theta_ab_tau_wb=stats_table_alpha_ab_tau_wb=stats_table_beta_ab_tau_wb=data.frame()
stats_table_delta_ab_tau_mr_no_ec=stats_table_theta_ab_tau_mr_no_ec=stats_table_alpha_ab_tau_mr_no_ec=stats_table_beta_ab_tau_mr_no_ec=data.frame()

### Initialize variables for all models (acc aperiodics)
stats_table_delta_ab_aper=stats_table_theta_ab_aper=stats_table_alpha_ab_aper=stats_table_beta_ab_aper=data.frame()
stats_table_delta_ab_etau_aper=stats_table_theta_ab_etau_aper=stats_table_alpha_ab_etau_aper=stats_table_beta_ab_etau_aper=data.frame()
stats_table_delta_ab_tau_mr_aper=stats_table_theta_ab_tau_mr_aper=stats_table_alpha_ab_tau_mr_aper=stats_table_beta_ab_tau_mr_aper=data.frame()
stats_table_delta_ab_tau_wb_aper=stats_table_theta_ab_tau_wb_aper=stats_table_alpha_ab_tau_wb_aper=stats_table_beta_ab_tau_wb_aper=data.frame()
stats_table_delta_ab_tau_mr_no_ec_aper=stats_table_theta_ab_tau_mr_no_ec_aper=stats_table_alpha_ab_tau_mr_no_ec_aper=stats_table_beta_ab_tau_mr_no_ec_aper=data.frame()


### Start iterations (select n)
for (j in 1:1000){
  print(j)
  ### Shuffle Ab and tau labels
  a=sample(1:68, 68, replace=F)
  Ab_rand=vector(); Tau_rand=vector()
  for (i in seq(1,nrow(data),68)){
    sub_ab=data$Ab_data_newr[seq(i,i+67)];
    sub_rand_ab=sub_ab[a]
    Ab_rand=cbind(Ab_rand,sub_rand_ab)
    sub_tau=data$Tau_data_newr[seq(i,i+67)];
    sub_rand_tau=sub_tau[a]
    Tau_rand=cbind(Tau_rand,sub_rand_tau)
  }
  
  ### Add shuffled variables to data frame
  Ab_rand=c(Ab_rand); Tau_rand=c(Tau_rand); data=cbind(data,Ab_rand,Tau_rand)
  

  ### Run null LME models MEG power ~ Amyloid
  stats_table_delta_ab <- rbind(stats_table_delta_ab, run_lme_model_perm("delta", "Ab_rand", "main_eff",data))
  stats_table_theta_ab <- rbind(stats_table_theta_ab, run_lme_model_perm("theta", "Ab_rand", "main_eff",data))
  stats_table_alpha_ab <- rbind(stats_table_alpha_ab, run_lme_model_perm("alpha", "Ab_rand", "main_eff",data))
  stats_table_beta_ab <- rbind(stats_table_beta_ab, run_lme_model_perm("beta", "Ab_rand", "main_eff",data))

  ### Run null LME models MEG power ~ Amyloid + aperiodics
  stats_table_delta_ab_aper <- rbind(stats_table_delta_ab_aper, run_lme_model_perm("delta", "Ab_rand + offset + exponent", "main_eff",data))
  stats_table_theta_ab_aper <- rbind(stats_table_theta_ab_aper, run_lme_model_perm("theta", "Ab_rand + offset + exponent", "main_eff",data))
  stats_table_alpha_ab_aper <- rbind(stats_table_alpha_ab_aper, run_lme_model_perm("alpha", "Ab_rand + offset + exponent", "main_eff",data))
  stats_table_beta_ab_aper <- rbind(stats_table_beta_ab_aper, run_lme_model_perm("beta", "Ab_rand + offset + exponent", "main_eff",data))

  ### Run null LME models MEG power ~ Amyloid * Tau entorhinal
  stats_table_delta_ab_etau <- rbind(stats_table_delta_ab_etau, run_lme_model_perm("delta", "Ab_rand * Tau_ent_newr", "int",data))
  stats_table_theta_ab_etau <- rbind(stats_table_theta_ab_etau, run_lme_model_perm("theta", "Ab_rand * Tau_ent_newr", "int",data))
  stats_table_alpha_ab_etau <- rbind(stats_table_alpha_ab_etau, run_lme_model_perm("alpha", "Ab_rand * Tau_ent_newr", "int",data))
  stats_table_beta_ab_etau <- rbind(stats_table_beta_ab_etau, run_lme_model_perm("beta", "Ab_rand * Tau_ent_newr", "int",data))

  ### Run null LME models MEG power - Amyloid * Tau entorhinal + aperiodics
  stats_table_delta_ab_etau_aper <- rbind(stats_table_delta_ab_etau_aper, run_lme_model_perm("delta", "Ab_rand * Tau_ent_newr + offset + exponent", "int_aper",data))
  stats_table_theta_ab_etau_aper <- rbind(stats_table_theta_ab_etau_aper, run_lme_model_perm("theta", "Ab_rand * Tau_ent_newr + offset + exponent", "int_aper",data))
  stats_table_alpha_ab_etau_aper <- rbind(stats_table_alpha_ab_etau_aper, run_lme_model_perm("alpha", "Ab_rand * Tau_ent_newr + offset + exponent", "int_aper",data))
  stats_table_beta_ab_etau_aper <- rbind(stats_table_beta_ab_etau_aper, run_lme_model_perm("beta", "Ab_rand * Tau_ent_newr + offset + exponent", "int_aper",data))

  
  ### Run null LME models MEG power ~ Amyloid * Tau meta-ROI
  stats_table_delta_ab_tau_mr <- rbind(stats_table_delta_ab_tau_mr, run_lme_model_perm("delta", "Ab_rand * Tau_meta_roi_newr", "int",data))
  stats_table_theta_ab_tau_mr <- rbind(stats_table_theta_ab_tau_mr, run_lme_model_perm("theta", "Ab_rand * Tau_meta_roi_newr", "int",data))
  stats_table_alpha_ab_tau_mr <- rbind(stats_table_alpha_ab_tau_mr, run_lme_model_perm("alpha", "Ab_rand * Tau_meta_roi_newr", "int",data))
  stats_table_beta_ab_tau_mr <- rbind(stats_table_beta_ab_tau_mr, run_lme_model_perm("beta", "Ab_rand * Tau_meta_roi_newr", "int",data))

  ### Run null LME models MEG power - Amyloid * Tau meta-ROI + aperiodics
  stats_table_delta_ab_tau_mr_aper <- rbind(stats_table_delta_ab_tau_mr_aper, run_lme_model_perm("delta", "Ab_rand * Tau_meta_roi_newr + offset + exponent", "int_aper",data))
  stats_table_theta_ab_tau_mr_aper <- rbind(stats_table_theta_ab_tau_mr_aper, run_lme_model_perm("theta", "Ab_rand * Tau_meta_roi_newr + offset + exponent", "int_aper",data))
  stats_table_alpha_ab_tau_mr_aper <- rbind(stats_table_alpha_ab_tau_mr_aper, run_lme_model_perm("alpha", "Ab_rand * Tau_meta_roi_newr + offset + exponent", "int_aper",data))
  stats_table_beta_ab_tau_mr_aper <- rbind(stats_table_beta_ab_tau_mr_aper, run_lme_model_perm("beta", "Ab_rand * Tau_meta_roi_newr + offset + exponent", "int_aper",data))

  
  ### Run null LME models MEG power ~ Amyloid * Tau WB
  stats_table_delta_ab_tau_wb <- rbind(stats_table_delta_ab_tau_wb, run_lme_model_perm("delta", "Ab_rand * Tau_rand", "int",data))
  stats_table_theta_ab_tau_wb <- rbind(stats_table_theta_ab_tau_wb, run_lme_model_perm("theta", "Ab_rand * Tau_rand", "int",data))
  stats_table_alpha_ab_tau_wb <- rbind(stats_table_alpha_ab_tau_wb, run_lme_model_perm("alpha", "Ab_rand * Tau_rand", "int",data))
  stats_table_beta_ab_tau_wb <- rbind(stats_table_beta_ab_tau_wb, run_lme_model_perm("beta", "Ab_rand * Tau_rand", "int",data))

  ### Run null LME models MEG power - Amyloid * Tau WB + aperiodics
  stats_table_delta_ab_tau_wb <- rbind(stats_table_delta_ab_tau_wb, run_lme_model_perm("delta", "Ab_rand * Tau_rand + offset + exponent", "int_aper",data))
  stats_table_theta_ab_tau_wb <- rbind(stats_table_theta_ab_tau_wb, run_lme_model_perm("theta", "Ab_rand * Tau_rand + offset + exponent", "int_aper",data))
  stats_table_alpha_ab_tau_wb <- rbind(stats_table_alpha_ab_tau_wb, run_lme_model_perm("alpha", "Ab_rand * Tau_rand + offset + exponent", "int_aper",data))
  stats_table_beta_ab_tau_wb <- rbind(stats_table_beta_ab_tau_wb, run_lme_model_perm("beta", "Ab_rand * Tau_rand + offset + exponent", "int_aper",data))

  ### Run null LME models MEG power ~ Amyloid * Tau meta-ROI no entorhinal
  stats_table_delta_ab_tau_mr_no_ec <- rbind(stats_table_delta_ab_tau_mr_no_ec, run_lme_model_perm("delta", "Ab_rand * Tau_meta_roi_newr_no_ec", "int",data))
  stats_table_theta_ab_tau_mr_no_ec <- rbind(stats_table_theta_ab_tau_mr_no_ec, run_lme_model_perm("theta", "Ab_rand * Tau_meta_roi_newr_no_ec", "int",data))
  stats_table_alpha_ab_tau_mr_no_ec <- rbind(stats_table_alpha_ab_tau_mr_no_ec, run_lme_model_perm("alpha", "Ab_rand * Tau_meta_roi_newr_no_ec", "int",data))
  stats_table_beta_ab_tau_mr_no_ec <- rbind(stats_table_beta_ab_tau_mr_no_ec, run_lme_model_perm("beta", "Ab_rand * Tau_meta_roi_newr_no_ec", "int",data))

  ### Run null LME models MEG power - Amyloid * Tau meta-ROI no entorhinal + aperiodics
  stats_table_delta_ab_tau_mr_no_ec_aper <- rbind(stats_table_delta_ab_tau_mr_no_ec_aper, run_lme_model_perm("delta", "Ab_rand * Tau_meta_roi_newr_no_ec + offset + exponent", "int_aper",data))
  stats_table_theta_ab_tau_mr_no_ec_aper <- rbind(stats_table_theta_ab_tau_mr_no_ec_aper, run_lme_model_perm("theta", "Ab_rand * Tau_meta_roi_newr_no_ec + offset + exponent", "int_aper",data))
  stats_table_alpha_ab_tau_mr_no_ec_aper <- rbind(stats_table_alpha_ab_tau_mr_no_ec_aper, run_lme_model_perm("alpha", "Ab_rand * Tau_meta_roi_newr_no_ec + offset + exponent", "int_aper",data))
  stats_table_beta_ab_tau_mr_no_ec_aper <- rbind(stats_table_beta_ab_tau_mr_no_ec_aper, run_lme_model_perm("beta", "Ab_rand * Tau_meta_roi_newr_no_ec + offset + exponent", "int_aper",data))

  ### Remove shuffled variables before starting next iteration
  data=subset(data, select = -c(Ab_rand,Tau_rand))
}






### 3. Save data from null LME models----
write.csv(stats_table_delta_ab, paste(path, "/delta_ab_null_data_lme", sep = ""))
write.csv(stats_table_theta_ab, paste(path, "/theta_ab_null_data_lme", sep = ""))
write.csv(stats_table_alpha_ab, paste(path, "/alpha_ab_null_data_lme", sep = ""))
write.csv(stats_table_beta_ab, paste(path, "/beta_ab_null_data_lme", sep = ""))

### Save data from null LME models acc aperiodic
write.csv(stats_table_delta_ab_aper, paste(path, "/delta_ab_aper_null_data_lme", sep = ""))
write.csv(stats_table_theta_ab_aper, paste(path, "/theta_ab_aper_null_data_lme", sep = ""))
write.csv(stats_table_alpha_ab_aper, paste(path, "/alpha_ab_aper_null_data_lme", sep = ""))
write.csv(stats_table_beta_ab_aper, paste(path, "/beta_ab_aper_null_data_lme", sep = ""))


### Save data from null LME models
write.csv(stats_table_delta_ab_etau, paste(path, "/delta_ab_etau_null_data_lme", sep = ""))
write.csv(stats_table_theta_ab_etau, paste(path, "/theta_ab_etau_null_data_lme", sep = ""))
write.csv(stats_table_alpha_ab_etau, paste(path, "/alpha_ab_etau_null_data_lme", sep = ""))
write.csv(stats_table_beta_ab_etau, paste(path, "/beta_ab_etau_null_data_lme", sep = ""))

### Save data from null LME models acc aperiodic
write.csv(stats_table_delta_ab_etau_aper, paste(path, "/delta_ab_etau_aper_null_data_lme", sep = ""))
write.csv(stats_table_theta_ab_etau_aper, paste(path, "/theta_ab_etau_aper_null_data_lme", sep = ""))
write.csv(stats_table_alpha_ab_etau_aper, paste(path, "/alpha_ab_etau_aper_null_data_lme", sep = ""))
write.csv(stats_table_beta_ab_etau_aper, paste(path, "/beta_ab_etau_aper_null_data_lme", sep = ""))


### Save data from null LME models
write.csv(stats_table_delta_ab_tau_mr, paste(path, "/delta_ab_tau_mr_null_data_lme", sep = ""))
write.csv(stats_table_theta_ab_tau_mr, paste(path, "/theta_ab_tau_mr_null_data_lme", sep = ""))
write.csv(stats_table_alpha_ab_tau_mr, paste(path, "/alpha_ab_tau_mr_null_data_lme", sep = ""))
write.csv(stats_table_beta_ab_tau_mr, paste(path, "/beta_ab_tau_mr_null_data_lme", sep = ""))

### Save data from null LME models acc aperiodic
write.csv(stats_table_delta_ab_tau_mr_aper, paste(path, "/delta_ab_tau_mr_aper_null_data_lme", sep = ""))
write.csv(stats_table_theta_ab_tau_mr_aper, paste(path, "/theta_ab_tau_mr_aper_null_data_lme", sep = ""))
write.csv(stats_table_alpha_ab_tau_mr_aper, paste(path, "/alpha_ab_tau_mr_aper_null_data_lme", sep = ""))
write.csv(stats_table_beta_ab_tau_mr_aper, paste(path, "/beta_ab_tau_mr_aper_null_data_lme", sep = ""))



#### Save data from null LME models
write.csv(stats_table_delta_ab_tau_wb, paste(path, "/delta_ab_tau_wb_abtau_same_null_data_lme", sep = ""))
write.csv(stats_table_theta_ab_tau_wb, paste(path, "/theta_ab_tau_wb_abtau_same_null_data_lme", sep = ""))
write.csv(stats_table_alpha_ab_tau_wb, paste(path, "/alpha_ab_tau_wb_abtau_same_null_data_lme", sep = ""))
write.csv(stats_table_beta_ab_tau_wb, paste(path, "/beta_ab_tau_wb_abtau_same_null_data_lme", sep = ""))

#### Save data from null LME models acc aperiodic
write.csv(stats_table_delta_ab_tau_wb_aper, paste(path, "/delta_ab_tau_wb_abtau_same_aper_null_data_lme", sep = ""))
write.csv(stats_table_theta_ab_tau_wb_aper, paste(path, "/theta_ab_tau_wb_abtau_same_aper_null_data_lme", sep = ""))
write.csv(stats_table_alpha_ab_tau_wb_aper, paste(path, "/alpha_ab_tau_wb_abtau_same_aper_null_data_lme", sep = ""))
write.csv(stats_table_beta_ab_tau_wb_aper, paste(path, "/beta_ab_tau_wb_abtau_same_aper_null_data_lme", sep = ""))



### Save data from null LME models
write.csv(stats_table_delta_ab_tau_mr_no_ec, paste(path, "/delta_ab_tau_mr_no_ec_null_data_lme", sep = ""))
write.csv(stats_table_theta_ab_tau_mr_no_ec, paste(path, "/theta_ab_tau_mr_no_ec_null_data_lme", sep = ""))
write.csv(stats_table_alpha_ab_tau_mr_no_ec, paste(path, "/alpha_ab_tau_mr_no_ec_null_data_lme", sep = ""))
write.csv(stats_table_beta_ab_tau_mr_no_ec, paste(path, "/beta_ab_tau_mr_no_ec_null_data_lme", sep = ""))

### Save data from null LME models acc aperiodic
write.csv(stats_table_delta_ab_tau_mr_no_ec_aper, paste(path, "/delta_ab_tau_mr_no_ec_aper_null_data_lme", sep = ""))
write.csv(stats_table_theta_ab_tau_mr_no_ec_aper, paste(path, "/theta_ab_tau_mr_no_ec_aper_null_data_lme", sep = ""))
write.csv(stats_table_alpha_ab_tau_mr_no_ec_aper, paste(path, "/alpha_ab_tau_mr_no_ec_aper_null_data_lme", sep = ""))
write.csv(stats_table_beta_ab_tau_mr_no_ec_aper, paste(path, "/beta_ab_tau_mr_no_ec_aper_null_data_lme", sep = ""))








### 4. Compute p values----
compute_p_val <- function(data, stats_table, outcome, predictor,aper,effect) {
  if(effect=="main"){cutoff=quantile(stats_table$t_ab,probs = seq(0,1,2.5/100)); low=cutoff[2]; high=cutoff[40]}
  if(effect=="int"){cutoff=quantile(stats_table$t_ab_tau,probs = seq(0,1,2.5/100)); low=cutoff[2]; high=cutoff[40]}
  if(aper=="no"){covariates="age + sex + edu + hipp_vol + trials"}; 
  if(aper=="yes"){covariates="age + sex + edu + hipp_vol + trials + offset + exponent"};
  formula <- as.formula(paste(outcome, " ~ ",predictor," + ", covariates))
  lme_model <- lme(formula, data = data, random = ~1 | ids/rois, control = lmeControl(opt = "optim"), na.action = na.exclude)
  if(effect=="main"){real_t_val=coef(summary(lme_model))[2, "t-value"]}; 
  if(effect=="int" && aper=="no"){real_t_val=coef(summary(lme_model))[9, "t-value"]}; 
  if(effect=="int" && aper=="yes"){real_t_val=coef(summary(lme_model))[11, "t-value"]};
  if(effect=="main"){
  if(outcome=="delta" || outcome=="theta"){p_val_np=sum(real_t_val>stats_table$t_ab); p_val_np=p_val_np/1000}
  if(outcome=="alpha" || outcome=="beta"){p_val_np=sum(real_t_val<stats_table$t_ab); p_val_np=p_val_np/1000}
    }
  if(effect=="int"){
    if(outcome=="delta" || outcome=="theta"){p_val_np=sum(real_t_val<stats_table$t_ab_tau); p_val_np=p_val_np/1000}
    if(outcome=="alpha" || outcome=="beta"){p_val_np=sum(real_t_val>stats_table$t_ab_tau); p_val_np=p_val_np/1000}
  }
  print(p_val_np)
  assign('p_val_np', p_val_np, envir=.GlobalEnv)
  #assign(paste(outcome,"_t_ab",sep=""), real_t_val, envir=.GlobalEnv)
  #return(data.frame(p_val_np))
  }


  ### Define variables of interest
  meg_variables=c("delta","theta","alpha","beta")
  pet_variables=c("Ab_data_newr","Ab_data_newr * Tau_ent_newr","Ab_data_newr * Tau_meta_roi_newr",
                  "Ab_data_newr * Tau_meta_roi_newr_no_ec","Ab_data_newr * Tau_data_newr")
  
  
  ### Ab main models
  ### Load data
  stats_table_delta_ab=read.csv('C:/Users/jogar/Documents/PREVENT_AD/Multimodal_data/lme_analysis/Main_models/newr_PET/null_data_perm/delta_ab_null_data_lme', sep=",", header = TRUE)
  stats_table_theta_ab=read.csv('C:/Users/jogar/Documents/PREVENT_AD/Multimodal_data/lme_analysis/Main_models/newr_PET/null_data_perm/theta_ab_null_data_lme', sep=",", header = TRUE)
  stats_table_alpha_ab=read.csv('C:/Users/jogar/Documents/PREVENT_AD/Multimodal_data/lme_analysis/Main_models/newr_PET/null_data_perm/alpha_ab_null_data_lme', sep=",", header = TRUE)
  stats_table_beta_ab=read.csv('C:/Users/jogar/Documents/PREVENT_AD/Multimodal_data/lme_analysis/Main_models/newr_PET/null_data_perm/beta_ab_null_data_lme', sep=",", header = TRUE)
  
  ### Get p values
  p_val_delta_ab=compute_p_val(data,stats_table_delta_ab,meg_variables[1],pet_variables[1],"no","main")
  p_val_theta_ab=compute_p_val(data,stats_table_theta_ab,meg_variables[2],pet_variables[1],"no","main")
  p_val_alpha_ab=compute_p_val(data,stats_table_alpha_ab,meg_variables[3],pet_variables[1],"no","main")
  p_val_beta_ab=compute_p_val(data,stats_table_beta_ab,meg_variables[4],pet_variables[1],"no","main")
  
  ### Save summary
  p_vals=c(p_val_delta_ab,p_val_theta_ab,p_val_alpha_ab,p_val_beta_ab); p_vals_fdr=p.adjust(p_vals,method = "fdr")
  model=lme(delta~Ab_data_newr+age+sex+edu+hipp_vol+trials,data=data,random=~1|ids/rois,control = lmeControl(opt = "optim")); real_t_val=coef(summary(model))[2, "t-value"]; 
  summary_stats_ab=data.frame(p_val=c(p_val_delta_ab, p_val_theta_ab, p_val_alpha_ab, p_val_beta_ab),
                              p_val_fdr=p_vals_fdr)
  print(summary_stats_ab)

  
  
  ### Ab main models + aperiodics
  ### Load data
  stats_table_delta_ab_aper=read.csv('C:/Users/jogar/Documents/PREVENT_AD/Multimodal_data/lme_analysis/Main_models/newr_PET/null_data_perm/delta_ab_aper_null_data_lme', sep=",", header = TRUE)
  stats_table_theta_ab_aper=read.csv('C:/Users/jogar/Documents/PREVENT_AD/Multimodal_data/lme_analysis/Main_models/newr_PET/null_data_perm/theta_ab_aper_null_data_lme', sep=",", header = TRUE)
  stats_table_alpha_ab_aper=read.csv('C:/Users/jogar/Documents/PREVENT_AD/Multimodal_data/lme_analysis/Main_models/newr_PET/null_data_perm/alpha_ab_aper_null_data_lme', sep=",", header = TRUE)
  stats_table_beta_ab_aper=read.csv('C:/Users/jogar/Documents/PREVENT_AD/Multimodal_data/lme_analysis/Main_models/newr_PET/null_data_perm/beta_ab_aper_null_data_lme', sep=",", header = TRUE)
  
  ### Get p values
  p_val_delta_ab_aper=compute_p_val(data,stats_table_delta_ab_aper,meg_variables[1],pet_variables[1],"yes","main")
  p_val_theta_ab_aper=compute_p_val(data,stats_table_theta_ab_aper,meg_variables[2],pet_variables[1],"yes","main")
  p_val_alpha_ab_aper=compute_p_val(data,stats_table_alpha_ab_aper,meg_variables[3],pet_variables[1],"yes","main")
  p_val_beta_ab_aper=compute_p_val(data,stats_table_beta_ab_aper,meg_variables[4],pet_variables[1],"yes","main")
  
  ### Save summary
  p_vals=c(p_val_delta_ab_aper,p_val_theta_ab_aper,p_val_alpha_ab_aper,p_val_beta_ab_aper); p_vals_fdr=p.adjust(p_vals,method = "fdr")
  model=lme(delta~Ab_data_newr+age+sex+edu+hipp_vol+trials,data=data,random=~1|ids/rois,control = lmeControl(opt = "optim")); real_t_val=coef(summary(model))[2, "t-value"]; 
  summary_stats_ab_aper=data.frame(p_val=c(p_val_delta_ab_aper, p_val_theta_ab_aper, p_val_alpha_ab_aper, p_val_beta_ab_aper),
                              p_val_fdr=p_vals_fdr)
  print(summary_stats_ab_aper)
  
  
  
  
  
  
  ### Ab x etau models
  ### Load data
  stats_table_delta_ab_etau=read.csv('C:/Users/jogar/Documents/PREVENT_AD/Multimodal_data/lme_analysis/Main_models/newr_PET/null_data_perm/delta_ab_etau_null_data_lme', sep=",", header = TRUE)
  stats_table_theta_ab_etau=read.csv('C:/Users/jogar/Documents/PREVENT_AD/Multimodal_data/lme_analysis/Main_models/newr_PET/null_data_perm/theta_ab_etau_null_data_lme', sep=",", header = TRUE)
  stats_table_alpha_ab_etau=read.csv('C:/Users/jogar/Documents/PREVENT_AD/Multimodal_data/lme_analysis/Main_models/newr_PET/null_data_perm/alpha_ab_etau_null_data_lme', sep=",", header = TRUE)
  stats_table_beta_ab_etau=read.csv('C:/Users/jogar/Documents/PREVENT_AD/Multimodal_data/lme_analysis/Main_models/newr_PET/null_data_perm/beta_ab_etau_null_data_lme', sep=",", header = TRUE)
  
  ### Get p values
  p_val_delta_ab_etau=compute_p_val(data,stats_table_delta_ab_etau,meg_variables[1],pet_variables[2],"no","int")
  p_val_theta_ab_etau=compute_p_val(data,stats_table_theta_ab_etau,meg_variables[2],pet_variables[2],"no","int")
  p_val_alpha_ab_etau=compute_p_val(data,stats_table_alpha_ab_etau,meg_variables[3],pet_variables[2],"no","int")
  p_val_beta_ab_etau=compute_p_val(data,stats_table_beta_ab_etau,meg_variables[4],pet_variables[2],"no","int")
  
  ### Save summary
  p_vals=c(p_val_delta_ab_etau,p_val_theta_ab_etau,p_val_alpha_ab_etau,p_val_beta_ab_etau); p_vals_fdr=p.adjust(p_vals,method = "fdr")
  model=lme(delta~Ab_data_newr+age+sex+edu+hipp_vol+trials,data=data,random=~1|ids/rois,control = lmeControl(opt = "optim")); real_t_val=coef(summary(model))[2, "t-value"]; 
  summary_stats_ab_etau=data.frame(p_val=c(p_val_delta_ab_etau, p_val_theta_ab_etau, p_val_alpha_ab_etau, p_val_beta_ab_etau),
                              p_val_fdr=p_vals_fdr)
  print(summary_stats_ab_etau)
  
  
  
  ### Ab x etau models + aperiodics
  ### Load data
  stats_table_delta_ab_etau_aper=read.csv('C:/Users/jogar/Documents/PREVENT_AD/Multimodal_data/lme_analysis/Main_models/newr_PET/null_data_perm/delta_ab_etau_aper_null_data_lme', sep=",", header = TRUE)
  stats_table_theta_ab_etau_aper=read.csv('C:/Users/jogar/Documents/PREVENT_AD/Multimodal_data/lme_analysis/Main_models/newr_PET/null_data_perm/theta_ab_etau_aper_null_data_lme', sep=",", header = TRUE)
  stats_table_alpha_ab_etau_aper=read.csv('C:/Users/jogar/Documents/PREVENT_AD/Multimodal_data/lme_analysis/Main_models/newr_PET/null_data_perm/alpha_ab_etau_aper_null_data_lme', sep=",", header = TRUE)
  stats_table_beta_ab_etau_aper=read.csv('C:/Users/jogar/Documents/PREVENT_AD/Multimodal_data/lme_analysis/Main_models/newr_PET/null_data_perm/beta_ab_etau_aper_null_data_lme', sep=",", header = TRUE)
  
  ### Get p values
  p_val_delta_ab_etau_aper=compute_p_val(data,stats_table_delta_ab_etau_aper,meg_variables[1],pet_variables[2],"yes","int")
  p_val_theta_ab_etau_aper=compute_p_val(data,stats_table_theta_ab_etau_aper,meg_variables[2],pet_variables[2],"yes","int")
  p_val_alpha_ab_etau_aper=compute_p_val(data,stats_table_alpha_ab_etau_aper,meg_variables[3],pet_variables[2],"yes","int")
  p_val_beta_ab_etau_aper=compute_p_val(data,stats_table_beta_ab_etau_aper,meg_variables[4],pet_variables[2],"yes","int")
  
  ### Save summary
  p_vals=c(p_val_delta_ab_etau_aper,p_val_theta_ab_etau_aper,p_val_alpha_ab_etau_aper,p_val_beta_ab_etau_aper); p_vals_fdr=p.adjust(p_vals,method = "fdr")
  model=lme(delta~Ab_data_newr+age+sex+edu+hipp_vol+trials,data=data,random=~1|ids/rois,control = lmeControl(opt = "optim")); real_t_val=coef(summary(model))[2, "t-value"]; 
  summary_stats_ab_etau_aper=data.frame(p_val=c(p_val_delta_ab_etau_aper, p_val_theta_ab_etau_aper, p_val_alpha_ab_etau_aper, p_val_beta_ab_etau_aper),
                                   p_val_fdr=p_vals_fdr)
  print(summary_stats_ab_etau_aper)
  
  
  
  
  ### Ab x tau meta-ROI models
  ### Load data
  stats_table_delta_ab_tau_mr=read.csv('C:/Users/jogar/Documents/PREVENT_AD/Multimodal_data/lme_analysis/Main_models/newr_PET/null_data_perm/delta_ab_tau_mr_null_data_lme', sep=",", header = TRUE)
  stats_table_theta_ab_tau_mr=read.csv('C:/Users/jogar/Documents/PREVENT_AD/Multimodal_data/lme_analysis/Main_models/newr_PET/null_data_perm/theta_ab_tau_mr_null_data_lme', sep=",", header = TRUE)
  stats_table_alpha_ab_tau_mr=read.csv('C:/Users/jogar/Documents/PREVENT_AD/Multimodal_data/lme_analysis/Main_models/newr_PET/null_data_perm/alpha_ab_tau_mr_null_data_lme', sep=",", header = TRUE)
  stats_table_beta_ab_tau_mr=read.csv('C:/Users/jogar/Documents/PREVENT_AD/Multimodal_data/lme_analysis/Main_models/newr_PET/null_data_perm/beta_ab_tau_mr_null_data_lme', sep=",", header = TRUE)
  
  ### Get p values
  p_val_delta_ab_tau_mr=compute_p_val(data,stats_table_delta_ab_tau_mr,meg_variables[1],pet_variables[3],"no","int")
  p_val_theta_ab_tau_mr=compute_p_val(data,stats_table_theta_ab_tau_mr,meg_variables[2],pet_variables[3],"no","int")
  p_val_alpha_ab_tau_mr=compute_p_val(data,stats_table_alpha_ab_tau_mr,meg_variables[3],pet_variables[3],"no","int")
  p_val_beta_ab_tau_mr=compute_p_val(data,stats_table_beta_ab_tau_mr,meg_variables[4],pet_variables[3],"no","int")
  
  ### Save summary
  p_vals=c(p_val_delta_ab_tau_mr,p_val_theta_ab_tau_mr,p_val_alpha_ab_tau_mr,p_val_beta_ab_tau_mr); p_vals_fdr=p.adjust(p_vals,method = "fdr")
  model=lme(delta~Ab_data_newr+age+sex+edu+hipp_vol+trials,data=data,random=~1|ids/rois,control = lmeControl(opt = "optim")); real_t_val=coef(summary(model))[2, "t-value"]; 
  summary_stats_ab_tau_mr=data.frame(p_val=c(p_val_delta_ab_tau_mr, p_val_theta_ab_tau_mr, p_val_alpha_ab_tau_mr, p_val_beta_ab_tau_mr),
                                   p_val_fdr=p_vals_fdr)
  print(summary_stats_ab_tau_mr)
  
  
  
  ### Ab x tau mr models + aperiodics
  ### Load data
  stats_table_delta_ab_tau_mr_aper=read.csv('C:/Users/jogar/Documents/PREVENT_AD/Multimodal_data/lme_analysis/Main_models/newr_PET/null_data_perm/delta_ab_tau_mr_aper_null_data_lme', sep=",", header = TRUE)
  stats_table_theta_ab_tau_mr_aper=read.csv('C:/Users/jogar/Documents/PREVENT_AD/Multimodal_data/lme_analysis/Main_models/newr_PET/null_data_perm/theta_ab_tau_mr_aper_null_data_lme', sep=",", header = TRUE)
  stats_table_alpha_ab_tau_mr_aper=read.csv('C:/Users/jogar/Documents/PREVENT_AD/Multimodal_data/lme_analysis/Main_models/newr_PET/null_data_perm/alpha_ab_tau_mr_aper_null_data_lme', sep=",", header = TRUE)
  stats_table_beta_ab_tau_mr_aper=read.csv('C:/Users/jogar/Documents/PREVENT_AD/Multimodal_data/lme_analysis/Main_models/newr_PET/null_data_perm/beta_ab_tau_mr_aper_null_data_lme', sep=",", header = TRUE)
  
  ### Get p values
  p_val_delta_ab_tau_mr_aper=compute_p_val(data,stats_table_delta_ab_tau_mr_aper,meg_variables[1],pet_variables[3],"yes","int")
  p_val_theta_ab_tau_mr_aper=compute_p_val(data,stats_table_theta_ab_tau_mr_aper,meg_variables[2],pet_variables[3],"yes","int")
  p_val_alpha_ab_tau_mr_aper=compute_p_val(data,stats_table_alpha_ab_tau_mr_aper,meg_variables[3],pet_variables[3],"yes","int")
  p_val_beta_ab_tau_mr_aper=compute_p_val(data,stats_table_beta_ab_tau_mr_aper,meg_variables[4],pet_variables[3],"yes","int")
  
  ### Save summary
  p_vals=c(p_val_delta_ab_tau_mr_aper,p_val_theta_ab_tau_mr_aper,p_val_alpha_ab_tau_mr_aper,p_val_beta_ab_tau_mr_aper); p_vals_fdr=p.adjust(p_vals,method = "fdr")
  model=lme(delta~Ab_data_newr+age+sex+edu+hipp_vol+trials,data=data,random=~1|ids/rois,control = lmeControl(opt = "optim")); real_t_val=coef(summary(model))[2, "t-value"]; 
  summary_stats_ab_tau_mr_aper=data.frame(p_val=c(p_val_delta_ab_tau_mr_aper, p_val_theta_ab_tau_mr_aper, p_val_alpha_ab_tau_mr_aper, p_val_beta_ab_tau_mr_aper),
                                        p_val_fdr=p_vals_fdr)
  print(summary_stats_ab_tau_mr_aper)
  
  
  
  ### Ab x tau meta-ROI no entorhinal models
  ### Load data
  stats_table_delta_ab_tau_mr_no_ec=read.csv('C:/Users/jogar/Documents/PREVENT_AD/Multimodal_data/lme_analysis/Main_models/newr_PET/null_data_perm/delta_ab_tau_mr_no_ec_null_data_lme', sep=",", header = TRUE)
  stats_table_theta_ab_tau_mr_no_ec=read.csv('C:/Users/jogar/Documents/PREVENT_AD/Multimodal_data/lme_analysis/Main_models/newr_PET/null_data_perm/theta_ab_tau_mr_no_ec_null_data_lme', sep=",", header = TRUE)
  stats_table_alpha_ab_tau_mr_no_ec=read.csv('C:/Users/jogar/Documents/PREVENT_AD/Multimodal_data/lme_analysis/Main_models/newr_PET/null_data_perm/alpha_ab_tau_mr_no_ec_null_data_lme', sep=",", header = TRUE)
  stats_table_beta_ab_tau_mr_no_ec=read.csv('C:/Users/jogar/Documents/PREVENT_AD/Multimodal_data/lme_analysis/Main_models/newr_PET/null_data_perm/beta_ab_tau_mr_no_ec_null_data_lme', sep=",", header = TRUE)
  
  ### Get p values
  p_val_delta_ab_tau_mr_no_ec=compute_p_val(data,stats_table_delta_ab_tau_mr_no_ec,meg_variables[1],pet_variables[4],"no","int")
  p_val_theta_ab_tau_mr_no_ec=compute_p_val(data,stats_table_theta_ab_tau_mr_no_ec,meg_variables[2],pet_variables[4],"no","int")
  p_val_alpha_ab_tau_mr_no_ec=compute_p_val(data,stats_table_alpha_ab_tau_mr_no_ec,meg_variables[3],pet_variables[4],"no","int")
  p_val_beta_ab_tau_mr_no_ec=compute_p_val(data,stats_table_beta_ab_tau_mr_no_ec,meg_variables[4],pet_variables[4],"no","int")
  
  ### Save summary
  p_vals=c(p_val_delta_ab_tau_mr_no_ec,p_val_theta_ab_tau_mr_no_ec,p_val_alpha_ab_tau_mr_no_ec,p_val_beta_ab_tau_mr_no_ec); p_vals_fdr=p.adjust(p_vals,method = "fdr")
  model=lme(delta~Ab_data_newr+age+sex+edu+hipp_vol+trials,data=data,random=~1|ids/rois,control = lmeControl(opt = "optim")); real_t_val=coef(summary(model))[2, "t-value"]; 
  summary_stats_ab_tau_mr_no_ec=data.frame(p_val=c(p_val_delta_ab_tau_mr_no_ec, p_val_theta_ab_tau_mr_no_ec, p_val_alpha_ab_tau_mr_no_ec, p_val_beta_ab_tau_mr_no_ec),
                                     p_val_fdr=p_vals_fdr)
  print(summary_stats_ab_tau_mr_no_ec)
  
  
  
  ### Ab x tau mr models + aperiodics
  ### Load data
  stats_table_delta_ab_tau_mr_no_ec_aper=read.csv('C:/Users/jogar/Documents/PREVENT_AD/Multimodal_data/lme_analysis/Main_models/newr_PET/null_data_perm/delta_ab_tau_mr_no_ec_aper_null_data_lme', sep=",", header = TRUE)
  stats_table_theta_ab_tau_mr_no_ec_aper=read.csv('C:/Users/jogar/Documents/PREVENT_AD/Multimodal_data/lme_analysis/Main_models/newr_PET/null_data_perm/theta_ab_tau_mr_no_ec_aper_null_data_lme', sep=",", header = TRUE)
  stats_table_alpha_ab_tau_mr_no_ec_aper=read.csv('C:/Users/jogar/Documents/PREVENT_AD/Multimodal_data/lme_analysis/Main_models/newr_PET/null_data_perm/alpha_ab_tau_mr_no_ec_aper_null_data_lme', sep=",", header = TRUE)
  stats_table_beta_ab_tau_mr_no_ec_aper=read.csv('C:/Users/jogar/Documents/PREVENT_AD/Multimodal_data/lme_analysis/Main_models/newr_PET/null_data_perm/beta_ab_tau_mr_no_ec_aper_null_data_lme', sep=",", header = TRUE)
  
  ### Get p values
  p_val_delta_ab_tau_mr_no_ec_aper=compute_p_val(data,stats_table_delta_ab_tau_mr_no_ec_aper,meg_variables[1],pet_variables[4],"yes","int")
  p_val_theta_ab_tau_mr_no_ec_aper=compute_p_val(data,stats_table_theta_ab_tau_mr_no_ec_aper,meg_variables[2],pet_variables[4],"yes","int")
  p_val_alpha_ab_tau_mr_no_ec_aper=compute_p_val(data,stats_table_alpha_ab_tau_mr_no_ec_aper,meg_variables[3],pet_variables[4],"yes","int")
  p_val_beta_ab_tau_mr_no_ec_aper=compute_p_val(data,stats_table_beta_ab_tau_mr_no_ec_aper,meg_variables[4],pet_variables[4],"yes","int")
  
  ### Save summary
  p_vals=c(p_val_delta_ab_tau_mr_no_ec_aper,p_val_theta_ab_tau_mr_no_ec_aper,p_val_alpha_ab_tau_mr_no_ec_aper,p_val_beta_ab_tau_mr_no_ec_aper); p_vals_fdr=p.adjust(p_vals,method = "fdr")
  model=lme(delta~Ab_data_newr+age+sex+edu+hipp_vol+trials,data=data,random=~1|ids/rois,control = lmeControl(opt = "optim")); real_t_val=coef(summary(model))[2, "t-value"]; 
  summary_stats_ab_tau_mr_no_ec_aper=data.frame(p_val=c(p_val_delta_ab_tau_mr_no_ec_aper, p_val_theta_ab_tau_mr_no_ec_aper, p_val_alpha_ab_tau_mr_no_ec_aper, p_val_beta_ab_tau_mr_no_ec_aper),
                                          p_val_fdr=p_vals_fdr)
  print(summary_stats_ab_tau_mr_no_ec_aper)
  
  
  
  
  ### Ab x tau WB models
  ### Load data
  stats_table_delta_ab_tau_wb=read.csv('C:/Users/jogar/Documents/PREVENT_AD/Multimodal_data/lme_analysis/Main_models/newr_PET/null_data_perm/delta_ab_tau_wb_null_data_lme', sep=",", header = TRUE)
  stats_table_theta_ab_tau_wb=read.csv('C:/Users/jogar/Documents/PREVENT_AD/Multimodal_data/lme_analysis/Main_models/newr_PET/null_data_perm/theta_ab_tau_wb_null_data_lme', sep=",", header = TRUE)
  stats_table_alpha_ab_tau_wb=read.csv('C:/Users/jogar/Documents/PREVENT_AD/Multimodal_data/lme_analysis/Main_models/newr_PET/null_data_perm/alpha_ab_tau_wb_null_data_lme', sep=",", header = TRUE)
  stats_table_beta_ab_tau_wb=read.csv('C:/Users/jogar/Documents/PREVENT_AD/Multimodal_data/lme_analysis/Main_models/newr_PET/null_data_perm/beta_ab_tau_wb_null_data_lme', sep=",", header = TRUE)
  
  ### Get p values
  p_val_delta_ab_tau_wb=compute_p_val(data,stats_table_delta_ab_tau_wb,meg_variables[1],pet_variables[5],"no","int")
  p_val_theta_ab_tau_wb=compute_p_val(data,stats_table_theta_ab_tau_wb,meg_variables[2],pet_variables[5],"no","int")
  p_val_alpha_ab_tau_wb=compute_p_val(data,stats_table_alpha_ab_tau_wb,meg_variables[3],pet_variables[5],"no","int")
  p_val_beta_ab_tau_wb=compute_p_val(data,stats_table_beta_ab_tau_wb,meg_variables[4],pet_variables[5],"no","int")
  
  ### Save summary
  p_vals=c(p_val_delta_ab_tau_wb,p_val_theta_ab_tau_wb,p_val_alpha_ab_tau_wb,p_val_beta_ab_tau_wb); p_vals_fdr=p.adjust(p_vals,method = "fdr")
  summary_stats_ab_tau_wb=data.frame(p_val=c(p_val_delta_ab_tau_wb, p_val_theta_ab_tau_wb, p_val_alpha_ab_tau_wb, p_val_beta_ab_tau_wb),
                                           p_val_fdr=p_vals_fdr)
  print(summary_stats_ab_tau_wb)
  
  
  
  ### Ab x tau WB
  ### Load data
  stats_table_delta_ab_tau_wb_aper=read.csv('C:/Users/jogar/Documents/PREVENT_AD/Multimodal_data/lme_analysis/Main_models/newr_PET/null_data_perm/delta_ab_tau_wb_aper_null_data_lme', sep=",", header = TRUE)
  stats_table_theta_ab_tau_wb_aper=read.csv('C:/Users/jogar/Documents/PREVENT_AD/Multimodal_data/lme_analysis/Main_models/newr_PET/null_data_perm/theta_ab_tau_wb_aper_null_data_lme', sep=",", header = TRUE)
  stats_table_alpha_ab_tau_wb_aper=read.csv('C:/Users/jogar/Documents/PREVENT_AD/Multimodal_data/lme_analysis/Main_models/newr_PET/null_data_perm/alpha_ab_tau_wb_aper_null_data_lme', sep=",", header = TRUE)
  stats_table_beta_ab_tau_wb_aper=read.csv('C:/Users/jogar/Documents/PREVENT_AD/Multimodal_data/lme_analysis/Main_models/newr_PET/null_data_perm/beta_ab_tau_wb_aper_null_data_lme', sep=",", header = TRUE)
  
  ### Get p values
  p_val_delta_ab_tau_wb_aper=compute_p_val(data,stats_table_delta_ab_tau_wb_aper,meg_variables[1],pet_variables[5],"yes","int")
  p_val_theta_ab_tau_wb_aper=compute_p_val(data,stats_table_theta_ab_tau_wb_aper,meg_variables[2],pet_variables[5],"yes","int")
  p_val_alpha_ab_tau_wb_aper=compute_p_val(data,stats_table_alpha_ab_tau_wb_aper,meg_variables[3],pet_variables[5],"yes","int")
  p_val_beta_ab_tau_wb_aper=compute_p_val(data,stats_table_beta_ab_tau_wb_aper,meg_variables[4],pet_variables[5],"yes","int")
  
  ### Save summary
  p_vals=c(p_val_delta_ab_tau_wb_aper,p_val_theta_ab_tau_wb_aper,p_val_alpha_ab_tau_wb_aper,p_val_beta_ab_tau_wb_aper); p_vals_fdr=p.adjust(p_vals,method = "fdr")
  model=lme(delta~Ab_data_newr+age+sex+edu+hipp_vol+trials,data=data,random=~1|ids/rois,control = lmeControl(opt = "optim")); real_t_val=coef(summary(model))[2, "t-value"]; 
  summary_stats_ab_tau_wb_aper=data.frame(p_val=c(p_val_delta_ab_tau_wb_aper, p_val_theta_ab_tau_wb_aper, p_val_alpha_ab_tau_wb_aper, p_val_beta_ab_tau_wb_aper),
                                                p_val_fdr=p_vals_fdr)
  print(summary_stats_ab_tau_wb_aper)
  
  
  