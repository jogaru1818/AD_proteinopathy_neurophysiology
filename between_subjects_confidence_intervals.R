### Script wrote by Jonathan Gallego Rudolf - PhD at McGill University

### Script used to calculate confidence intervals for the biomarker-defined PET groups 
### using bootstrapping resampling to complement the between subject-analysis


### Confidence intervals are computed for the Ab/Tau PET groups in the
### 1. MEG spectral power between-subjetcts analysis
### 2. Longitudinal cognitve score slopes between-subjects analysis 




### Prepare data
### The input is the data_meg_pet_cogn2 dataset generated in the between_subjects_analysis.R script

### Set path
path="C:/Users/jogar/Documents/"





### 1. MEG spectral power between-subjects analysis----

### MEG power
### Calculate real mean for each subgroup
delta_real_means=aggregate(data_meg_pet_cogn2$mean_meg_delta, list(data_meg_pet_cogn2$amyloid_tau_status), FUN=mean) 
theta_real_means=aggregate(data_meg_pet_cogn2$mean_meg_theta, list(data_meg_pet_cogn2$amyloid_tau_status), FUN=mean) 
alpha_real_means=aggregate(data_meg_pet_cogn2$mean_meg_alpha, list(data_meg_pet_cogn2$amyloid_tau_status), FUN=mean) 
beta_real_means=aggregate(data_meg_pet_cogn2$mean_meg_beta, list(data_meg_pet_cogn2$amyloid_tau_status), FUN=mean) 


### Save real MEG power means
real_means_MEG_power=data.frame(cbind(delta_real_means$x,theta_real_means$x,alpha_real_means$x,beta_real_means$x))
colnames(real_means_MEG_power)=c("delta","theta","alpha","beta")
rownames(real_means_MEG_power)=c("Ab-/Tau-","Ab+/Tau-","Ab+/Tau+")
# write.csv(real_means_MEG_power,paste(path,"real_means_MEG_power.csv", sep=""))


### Separate datasets based on amyloid-tau status
x=which(data_meg_pet_cogn2$amyloid_tau_status=="Ab-/Tau-"); data_group_1=data_meg_pet_cogn2[x,]
x=which(data_meg_pet_cogn2$amyloid_tau_status=="Ab+/Tau-"); data_group_2=data_meg_pet_cogn2[x,]
x=which(data_meg_pet_cogn2$amyloid_tau_status=="Ab+/Tau+"); data_group_3=data_meg_pet_cogn2[x,]



### Use bootstrapping resampling to obtain condifidence intervals for the mean of each group
### Initialize variable
confidence_intervals_data=vector()

### Iterations
for (i in 1:10000){ # select number of iterations for the bootstrapping resampling
  print(i)
  ### Define that the imblanced proportion between the group's sample size will be preserved for the resampling
  rand_1=sample(1:50, replace=T); rand_2=sample(1:42, replace=T); rand_3=sample(1:11, replace=T)
  
  ### Grab data from each group for this iteration
  data_group_1_boot_samp=data_group_1[rand_1,]; data_group_2_boot_samp=data_group_2[rand_2,]; data_group_3_boot_samp=data_group_3[rand_3,]

  
  ### MEG power
  ### Calculate mean for each group for each bootstrapping iteration
  delta_boot_mean_group_1=mean(data_group_1_boot_samp$mean_meg_delta); theta_boot_mean_group_1=mean(data_group_1_boot_samp$mean_meg_theta)
  alpha_boot_mean_group_1=mean(data_group_1_boot_samp$mean_meg_alpha); beta_boot_mean_group_1=mean(data_group_1_boot_samp$mean_meg_beta)
  
  delta_boot_mean_group_2=mean(data_group_2_boot_samp$mean_meg_delta); theta_boot_mean_group_2=mean(data_group_2_boot_samp$mean_meg_theta)
  alpha_boot_mean_group_2=mean(data_group_2_boot_samp$mean_meg_alpha); beta_boot_mean_group_2=mean(data_group_2_boot_samp$mean_meg_beta)
  
  delta_boot_mean_group_3=mean(data_group_3_boot_samp$mean_meg_delta); theta_boot_mean_group_3=mean(data_group_3_boot_samp$mean_meg_theta)
  alpha_boot_mean_group_3=mean(data_group_3_boot_samp$mean_meg_alpha); beta_boot_mean_group_3=mean(data_group_3_boot_samp$mean_meg_beta)
  
  
  ### Create dataset with the mean values of each PET group for each band obtained from each  
  ### iteration to generate the mean distribution from which the CI will be estimated
  confidence_intervals_vec=cbind(boot_num=i,delta_boot_mean_group_1,delta_boot_mean_group_2,delta_boot_mean_group_3,
                                 theta_boot_mean_group_1,theta_boot_mean_group_2,theta_boot_mean_group_3,
                                 alpha_boot_mean_group_1,alpha_boot_mean_group_2,alpha_boot_mean_group_3,
                                 beta_boot_mean_group_1,beta_boot_mean_group_2,beta_boot_mean_group_3)
  confidence_intervals_data=rbind(confidence_intervals_data,confidence_intervals_vec)
}

### Transform to datasets and save resampled data
confidence_intervals_data=as.data.frame(confidence_intervals_data)
write.csv(confidence_intervals_data,paste(path,"ci_data_meg_power_imbalanced.csv", sep=""))



### Calculate confidence intervals for the mean of each group

### delta Ab-/Tau-
cutoff_1=quantile(confidence_intervals_data$delta_boot_mean_group_1,probs = seq(0,1,2.5/100))
delta_low_1=cutoff_1[2]; delta_high_1=cutoff_1[40]; delta_boot_mean_1=mean(confidence_intervals_data$delta_boot_mean_group_1)

### delta Ab+/Tau-
cutoff_2=quantile(confidence_intervals_data$delta_boot_mean_group_2,probs = seq(0,1,2.5/100))
delta_low_2=cutoff_2[2]; delta_high_2=cutoff_2[40]; delta_boot_mean_2=mean(confidence_intervals_data$delta_boot_mean_group_2)

### delta Ab+/Tau+
cutoff_3=quantile(confidence_intervals_data$delta_boot_mean_group_3,probs = seq(0,1,2.5/100))
delta_low_3=cutoff_3[2]; delta_high_3=cutoff_3[40]; delta_boot_mean_3=mean(confidence_intervals_data$delta_boot_mean_group_3)


### Immediate memory Ab-/Tau-
cutoff_1=quantile(confidence_intervals_data$theta_boot_mean_group_1,probs = seq(0,1,2.5/100))
theta_low_1=cutoff_1[2]; theta_high_1=cutoff_1[40]; theta_boot_mean_1=mean(confidence_intervals_data$theta_boot_mean_group_1)

### Immediate memory Ab+/Tau-
cutoff_2=quantile(confidence_intervals_data$theta_boot_mean_group_2,probs = seq(0,1,2.5/100))
theta_low_2=cutoff_2[2]; theta_high_2=cutoff_2[40]; theta_boot_mean_2=mean(confidence_intervals_data$theta_boot_mean_group_2)

### Immediate memory Ab+/Tau+
cutoff_3=quantile(confidence_intervals_data$theta_boot_mean_group_3,probs = seq(0,1,2.5/100))
theta_low_3=cutoff_3[2]; theta_high_3=cutoff_3[40]; theta_boot_mean_3=mean(confidence_intervals_data$theta_boot_mean_group_3)


### Delayed memory Ab-/Tau-
cutoff_1=quantile(confidence_intervals_data$alpha_boot_mean_group_1,probs = seq(0,1,2.5/100))
alpha_low_1=cutoff_1[2]; alpha_high_1=cutoff_1[40]; alpha_boot_mean_1=mean(confidence_intervals_data$alpha_boot_mean_group_1)

### Delayed memory Ab+/Tau-
cutoff_2=quantile(confidence_intervals_data$alpha_boot_mean_group_2,probs = seq(0,1,2.5/100))
alpha_low_2=cutoff_2[2]; alpha_high_2=cutoff_2[40]; alpha_boot_mean_2=mean(confidence_intervals_data$alpha_boot_mean_group_2)

### Delayed memory Ab+/Tau+
cutoff_3=quantile(confidence_intervals_data$alpha_boot_mean_group_3,probs = seq(0,1,2.5/100))
alpha_low_3=cutoff_3[2]; alpha_high_3=cutoff_3[40]; alpha_boot_mean_3=mean(confidence_intervals_data$alpha_boot_mean_group_3)


### Att mem socre Ab-/Tau-
cutoff_1=quantile(confidence_intervals_data$beta_boot_mean_group_1,probs = seq(0,1,2.5/100))
beta_low_1=cutoff_1[2]; beta_high_1=cutoff_1[40]; beta_boot_mean_1=mean(confidence_intervals_data$beta_boot_mean_group_1)

### Att mem socre Ab+/Tau-
cutoff_2=quantile(confidence_intervals_data$beta_boot_mean_group_2,probs = seq(0,1,2.5/100))
beta_low_2=cutoff_2[2]; beta_high_2=cutoff_2[40]; beta_boot_mean_2=mean(confidence_intervals_data$beta_boot_mean_group_2)

### Att mem socre Ab+/Tau+
cutoff_3=quantile(confidence_intervals_data$beta_boot_mean_group_3,probs = seq(0,1,2.5/100))
beta_low_3=cutoff_3[2]; beta_high_3=cutoff_3[40]; beta_boot_mean_3=mean(confidence_intervals_data$beta_boot_mean_group_3)



### Save as dataframe
meg_power_ci=data.frame(rbind(c(delta_low_1,delta_high_1,delta_boot_mean_1),c(delta_low_2,delta_high_2,delta_boot_mean_2),c(delta_low_3,delta_high_3,delta_boot_mean_3),
                              c(theta_low_1,theta_high_1,theta_boot_mean_1),c(theta_low_2,theta_high_2,theta_boot_mean_2),c(theta_low_3,theta_high_3,theta_boot_mean_3),
                              c(alpha_low_1,alpha_high_1,alpha_boot_mean_1),c(alpha_low_2,alpha_high_2,alpha_boot_mean_2),c(alpha_low_3,alpha_high_3,alpha_boot_mean_3),
                              c(beta_low_1,beta_high_1,beta_boot_mean_1),c(beta_low_2,beta_high_2,beta_boot_mean_2),c(beta_low_3,beta_high_3,beta_boot_mean_3)))

colnames(meg_power_ci)=c("low_CI", "high_CI", "boot_mean")
rownames(meg_power_ci)=c("delta_Ab-/Tau-", "delta_Ab+/Tau-", "delta_Ab+/Tau+",
                         "theta_Ab-/Tau-", "theta_Ab+/Tau-", "theta_Ab+/Tau+",
                         "alpha_Ab-/Tau-", "alpha_Ab+/Tau-", "alpha_Ab+/Tau+",
                         "beta_Ab-/Tau-", "beta_Ab+/Tau-", "beta_Ab+/Tau+")

write.csv(meg_power_ci,paste(path,"ci_meg_power.csv", sep=""))









### 2. Longitudinal cognitve score slopes between-subjects analysis----

### Calculate real mean for each subgroup
attention_real_means=aggregate(data_meg_pet_cogn3$attention, list(data_meg_pet_cogn3$amyloid_tau_status), FUN=mean) 
imm_mem_real_means=aggregate(data_meg_pet_cogn3$immediate_memory, list(data_meg_pet_cogn3$amyloid_tau_status), FUN=mean) 
del_mem_real_means=aggregate(data_meg_pet_cogn3$delayed_memory, list(data_meg_pet_cogn3$amyloid_tau_status), FUN=mean) 
att_mem_score_real_means=aggregate(data_meg_pet_cogn3$att_mem_score, list(data_meg_pet_cogn3$amyloid_tau_status), FUN=mean) 


### Save real cognition means
real_means_cogn=data.frame(cbind(attention_real_means$x,imm_mem_real_means$x,del_mem_real_means$x,att_mem_score_real_means$x))
colnames(real_means_cogn)=c("attention","imm_mem","del_mem","att_mem_score")
rownames(real_means_cogn)=c("Ab-/Tau-","Ab+/Tau-","Ab+/Tau+")
# write.csv(real_means_cogn,paste(path,"real_means_cogn.csv", sep=""))


### Separate datasets based on amyloid-tau status
x=which(data_meg_pet_cogn3$amyloid_tau_status=="Ab-/Tau-"); data_group_1_cogn=data_meg_pet_cogn3[x,]
x=which(data_meg_pet_cogn3$amyloid_tau_status=="Ab+/Tau-"); data_group_2_cogn=data_meg_pet_cogn3[x,]
x=which(data_meg_pet_cogn3$amyloid_tau_status=="Ab+/Tau+"); data_group_3_cogn=data_meg_pet_cogn3[x,]


### Initialize variable
confidence_intervals_cognition_data=vector()


### Iterations
for (i in 1:10000){ # select number of iterations for the bootstrapping resampling
  print(i)
  ### Define that the imblanced proportion between the group's sample size will be preserved for the resampling
  rand_1=sample(1:49, replace=T); rand_2=sample(1:42, replace=T); rand_3=sample(1:11, replace=T)
  ### Grab data from each group for this iteration
  data_group_1_boot_samp=data_group_1_cogn[rand_1,]; data_group_2_boot_samp=data_group_2_cogn[rand_2,]; data_group_3_boot_samp=data_group_3_cogn[rand_3,]
  
  ### Calculate mean for each group for each bootstrapping iteration
  attention_boot_mean_group_1=mean(data_group_1_boot_samp$attention); imm_mem_boot_mean_group_1=mean(data_group_1_boot_samp$immediate_memory)
  del_mem_boot_mean_group_1=mean(data_group_1_boot_samp$delayed_memory); att_mem_score_boot_mean_group_1=mean(data_group_1_boot_samp$att_mem_score)
  
  attention_boot_mean_group_2=mean(data_group_2_boot_samp$attention); imm_mem_boot_mean_group_2=mean(data_group_2_boot_samp$immediate_memory)
  del_mem_boot_mean_group_2=mean(data_group_2_boot_samp$delayed_memory); att_mem_score_boot_mean_group_2=mean(data_group_2_boot_samp$att_mem_score)
  
  attention_boot_mean_group_3=mean(data_group_3_boot_samp$attention); imm_mem_boot_mean_group_3=mean(data_group_3_boot_samp$immediate_memory)
  del_mem_boot_mean_group_3=mean(data_group_3_boot_samp$delayed_memory); att_mem_score_boot_mean_group_3=mean(data_group_3_boot_samp$att_mem_score)
  
  ### Create dataset with the mean values of each PET group for each band obtained from each  
  ### iteration to generate the mean distribution from which the CI will be estimated
  confidence_intervals_cognition_vec=cbind(boot_num=i,attention_boot_mean_group_1,attention_boot_mean_group_2,attention_boot_mean_group_3,
                                           imm_mem_boot_mean_group_1,imm_mem_boot_mean_group_2,imm_mem_boot_mean_group_3,
                                           del_mem_boot_mean_group_1,del_mem_boot_mean_group_2,del_mem_boot_mean_group_3,
                                           att_mem_score_boot_mean_group_1,att_mem_score_boot_mean_group_2,att_mem_score_boot_mean_group_3)
  confidence_intervals_cognition_data=rbind(confidence_intervals_cognition_data,confidence_intervals_cognition_vec)
}


### Transform to datasets and save resampled data
confidence_intervals_cognition_data=as.data.frame(confidence_intervals_cognition_data)
write.csv(confidence_intervals_cognition_data,paste(path,"ci_data_cognition_imbalanced.csv", sep=""))


### Calculate confidence intervals for the mean of each group

### Attention Ab-/Tau-
cutoff_1=quantile(confidence_intervals_cognition_data$attention_boot_mean_group_1,probs = seq(0,1,2.5/100))
att_low_1=cutoff_1[2]; att_high_1=cutoff_1[40]; att_boot_mean_1=mean(confidence_intervals_cognition_data$attention_boot_mean_group_1)

### Attention Ab+/Tau-
cutoff_2=quantile(confidence_intervals_cognition_data$attention_boot_mean_group_2,probs = seq(0,1,2.5/100))
att_low_2=cutoff_2[2]; att_high_2=cutoff_2[40]; att_boot_mean_2=mean(confidence_intervals_cognition_data$attention_boot_mean_group_2)

### Attention Ab+/Tau+
cutoff_3=quantile(confidence_intervals_cognition_data$attention_boot_mean_group_3,probs = seq(0,1,2.5/100))
att_low_3=cutoff_3[2]; att_high_3=cutoff_3[40]; att_boot_mean_3=mean(confidence_intervals_cognition_data$attention_boot_mean_group_3)


### Immediate memory Ab-/Tau-
cutoff_1=quantile(confidence_intervals_cognition_data$imm_mem_boot_mean_group_1,probs = seq(0,1,2.5/100))
imm_mem_low_1=cutoff_1[2]; imm_mem_high_1=cutoff_1[40]; imm_mem_boot_mean_1=mean(confidence_intervals_cognition_data$imm_mem_boot_mean_group_1)

### Immediate memory Ab+/Tau-
cutoff_2=quantile(confidence_intervals_cognition_data$imm_mem_boot_mean_group_2,probs = seq(0,1,2.5/100))
imm_mem_low_2=cutoff_2[2]; imm_mem_high_2=cutoff_2[40]; imm_mem_boot_mean_2=mean(confidence_intervals_cognition_data$imm_mem_boot_mean_group_2)

### Immediate memory Ab+/Tau+
cutoff_3=quantile(confidence_intervals_cognition_data$imm_mem_boot_mean_group_3,probs = seq(0,1,2.5/100))
imm_mem_low_3=cutoff_3[2]; imm_mem_high_3=cutoff_3[40]; imm_mem_boot_mean_3=mean(confidence_intervals_cognition_data$imm_mem_boot_mean_group_3)


### Delayed memory Ab-/Tau-
cutoff_1=quantile(confidence_intervals_cognition_data$del_mem_boot_mean_group_1,probs = seq(0,1,2.5/100))
del_mem_low_1=cutoff_1[2]; del_mem_high_1=cutoff_1[40]; del_mem_boot_mean_1=mean(confidence_intervals_cognition_data$del_mem_boot_mean_group_1)

### Delayed memory Ab+/Tau-
cutoff_2=quantile(confidence_intervals_cognition_data$del_mem_boot_mean_group_2,probs = seq(0,1,2.5/100))
del_mem_low_2=cutoff_2[2]; del_mem_high_2=cutoff_2[40]; del_mem_boot_mean_2=mean(confidence_intervals_cognition_data$del_mem_boot_mean_group_2)

### Delayed memory Ab+/Tau+
cutoff_3=quantile(confidence_intervals_cognition_data$del_mem_boot_mean_group_3,probs = seq(0,1,2.5/100))
del_mem_low_3=cutoff_3[2]; del_mem_high_3=cutoff_3[40]; del_mem_boot_mean_3=mean(confidence_intervals_cognition_data$del_mem_boot_mean_group_3)


### Att mem socre Ab-/Tau-
cutoff_1=quantile(confidence_intervals_cognition_data$att_mem_score_boot_mean_group_1,probs = seq(0,1,2.5/100))
att_mem_score_low_1=cutoff_1[2]; att_mem_score_high_1=cutoff_1[40]; att_mem_score_boot_mean_1=mean(confidence_intervals_cognition_data$att_mem_score_boot_mean_group_1)

### Att mem socre Ab+/Tau-
cutoff_2=quantile(confidence_intervals_cognition_data$att_mem_score_boot_mean_group_2,probs = seq(0,1,2.5/100))
att_mem_score_low_2=cutoff_2[2]; att_mem_score_high_2=cutoff_2[40]; att_mem_score_boot_mean_2=mean(confidence_intervals_cognition_data$att_mem_score_boot_mean_group_2)

### Att mem socre Ab+/Tau+
cutoff_3=quantile(confidence_intervals_cognition_data$att_mem_score_boot_mean_group_3,probs = seq(0,1,2.5/100))
att_mem_score_low_3=cutoff_3[2]; att_mem_score_high_3=cutoff_3[40]; att_mem_score_boot_mean_3=mean(confidence_intervals_cognition_data$att_mem_score_boot_mean_group_3)


### Save as dataframe
cognition_ci=data.frame(rbind(c(att_low_1,att_high_1,att_boot_mean_1),c(att_low_2,att_high_2,att_boot_mean_2),c(att_low_3,att_high_3,att_boot_mean_3),
                              c(imm_mem_low_1,imm_mem_high_1,imm_mem_boot_mean_1),c(imm_mem_low_2,imm_mem_high_2,imm_mem_boot_mean_2),c(imm_mem_low_3,imm_mem_high_3,imm_mem_boot_mean_3),
                              c(del_mem_low_1,del_mem_high_1,del_mem_boot_mean_1),c(del_mem_low_2,del_mem_high_2,del_mem_boot_mean_2),c(del_mem_low_3,del_mem_high_3,del_mem_boot_mean_3),
                              c(att_mem_score_low_1,att_mem_score_high_1,att_mem_score_boot_mean_1),c(att_mem_score_low_2,att_mem_score_high_2,att_mem_score_boot_mean_2),c(att_mem_score_low_3,att_mem_score_high_3,att_mem_score_boot_mean_3)))

colnames(cognition_ci)=c("low_CI", "high_CI", "boot_mean")
rownames(cognition_ci)=c("att_Ab-/Tau-", "att_Ab+/Tau-", "att_Ab+/Tau+",
                         "imm_mem_Ab-/Tau-", "imm_mem_Ab+/Tau-", "imm_mem_Ab+/Tau+",
                         "del_mem_Ab-/Tau-", "del_mem_Ab+/Tau-", "del_mem_Ab+/Tau+",
                         "att_mem_score_Ab-/Tau-", "att_mem_score_Ab+/Tau-", "att_mem_score_Ab+/Tau+")

write.csv(cognition_ci,paste(path,"ci_cognition.csv", sep=""))




