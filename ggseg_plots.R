
### Script wrote by Jonathan Gallego Rudolf - PhD at McGill University

### Script used to generate plots of data overlayed on the brain (parcellated using 
### the Desikan Killiany atalas) by means of ggseg toolbox. Input data as in the format
### of the glm spreadsheet. Script calculates mean values per each region of interest 
### in order to plot group level maps. Also to create plots for individual 
### subjects.

### 1. Get data and load libraries
### 2. Calculate average across subjects for each ROI
### 3. Generate ggseg plots
### 4. Extra plots Ab and Tau ROIs




### 1. Get data and load libraries----

### Load data and required packages
require(ggplot2)
library(dplyr)
library(tidyr)
library(ggseg)
library(ggseg3d)
require(nlme)
require(MuMIn)
require(grid)
library(viridis) # Parula colorscale


### Load data
data = read.csv(file.choose(),header=T)

### Set path to save data 
path="C:/Users/jogar/Documents/PREVENT_AD/Multimodal_data/Figures"

### Separate datasets based on amyloid-tau status
x=which(data$amyloid_tau_status=="Ab-/Tau-"); data_group_1=data[x,]
x=which(data$amyloid_tau_status=="Ab+/Tau-"); data_group_2=data[x,]
x=which(data$amyloid_tau_status=="Ab+/Tau+"); data_group_3=data[x,]





### 2. Calculate average across subjects for each ROI----

### Initialize variables

### PET
ab_avg=data[FALSE,]; tau_avg=data[FALSE,]; ab_idx_rois=data[FALSE,]; tau_ent=data[FALSE,]; tau_meta_roi_rois=data[FALSE,]

### MEG power
delta_avg=data[FALSE,]; theta_avg=data[FALSE,]; alpha_avg=data[FALSE,]; beta_avg=data[FALSE,]; 

### FOOOF parameters - Aperiodic
exponent_avg <- data[FALSE,]; offset_avg <- data[FALSE,]
exponent_std <- data[FALSE,]; offset_std <- data[FALSE,]


### Generate function to average data per ROI o calculate MEG and PET estimates to plot in ggseg
average_rois_ggseg <- function(data, group) {
  ### Calculate average across individuals
  for (i in 1:68){
    idx=seq(i,nrow(data),68)
    # PET data
    ab=data$Ab_data_newr[idx]; ab=mean(ab); ab_avg=rbind(ab_avg,ab)
    tau=data$Tau_data_newr[idx]; tau=mean(tau); tau_avg=rbind(tau_avg,tau)
    ab_idx_rois=c(0,1,1,0,0,0,1,0,1,0,1,0,1,1,0,0,1,1,1,0,0,1,0,1,1,1,1,1,1,1,1,0,0,0); ab_idx_rois=rep(ab_idx_rois,2); ab_idx_rois=as.factor(ab_idx_rois)
    tau_meta_roi_rois=rep(0,34); tau_meta_roi_rois[c(5,6,8,12,14,15)]=1; tau_meta_roi_rois=rep(tau_meta_roi_rois,2); tau_meta_roi_rois=as.factor(tau_meta_roi_rois)
    tau_ent=rep(0,34); tau_ent[c(5)]=1; tau_ent=rep(tau_ent,2); tau_ent=as.factor(tau_ent)
  
    # MEG data
    delta=data$delta[idx]; delta=mean(delta); delta_avg=rbind(delta_avg,delta)
    theta=data$theta[idx]; theta=mean(theta); theta_avg=rbind(theta_avg,theta)
    alpha=data$alpha[idx]; alpha=mean(alpha); alpha_avg=rbind(alpha_avg,alpha)
    beta=data$beta[idx]; beta=mean(beta); beta_avg=rbind(beta_avg,beta)
  
    # Spec params
    exponent=data$exponent[idx]; exponent_m=mean(exponent); exponent_avg=rbind(exponent_avg,exponent_m)
    exponent_s=sd(exponent); exponent_std=rbind(exponent_std,exponent_s); 
    offset=data$offset[idx]; offset_m=mean(offset); offset_avg=rbind(offset_avg,offset_m); 
    offset_s=sd(offset); offset_std=rbind(offset_std,offset_s)
  }

  ### Set colnames PET data
  colnames(ab_avg) <- c("Ab_SUVR"); colnames(tau_avg) <- c("Tau_SUVR")
  ab_idx_rois=data.frame(ab_idx_rois); colnames(ab_idx_rois)=c("ab_idx_rois") 
  tau_meta_roi_rois=data.frame(tau_meta_roi_rois); colnames(tau_meta_roi_rois)=c("tau_meta_roi_rois")
  tau_ent=data.frame(tau_ent); colnames(tau_ent)=c("tau_ent")

  ### Set colnames MEG data
  colnames(delta_avg) <- c("Delta"); colnames(theta_avg) <- c("Theta"); 
  colnames(alpha_avg) <- c("Alpha"); colnames(beta_avg) <- c("Beta");

  ### Set colnames for spec param data
  colnames(exponent_avg) <- c("Exponent"); colnames(offset_avg) <- c("Offset")

  ### Create dataset 
  data_ggseg <- tibble(
  region = rep(c("bankssts","caudal anterior cingulate","caudal middle frontal","cuneus",
                 "entorhinal","fusiform","inferior parietal","inferior temporal","isthmus cingulate",
                 "lateral occipital","lateral orbitofrontal","lingual","medial orbitofrontal",
                 "middle temporal","parahippocampal","paracentral","pars opercularis","pars orbitalis",
                 "pars triangularis","pericalcarine","postcentral","posterior cingulate","precentral",
                 "precuneus","rostral anterior cingulate","rostral middle frontal","superior frontal",
                 "superior parietal","superior temporal","supramarginal","frontal pole","temporal pole",
                 "transverse temporal", "insula"), 2), 
  Ab_SUVR=ab_avg$Ab_SUVR, Tau_SUVR=tau_avg$Tau_SUVR,Delta=delta_avg$Delta,Theta=theta_avg$Theta,
  Alpha=alpha_avg$Alpha, Beta=beta_avg$Beta,Exponent=exponent_avg$Exponent, Offset=offset_avg$Offset,
  Tau_mr=tau_meta_roi_rois$tau_meta_roi_rois,Ab_idx=ab_idx_rois$ab_idx_rois, Tau_ent=tau_ent$tau_ent,
  hemi = c(rep("left", 34), rep("right", 34)))
  
  if(group=="all"){assign('data_ggseg_all', data_ggseg, envir=.GlobalEnv)}; 
  if(group=="g1"){assign('data_ggseg_g1', data_ggseg, envir=.GlobalEnv)};
  if(group=="g2"){assign('data_ggseg_g2', data_ggseg, envir=.GlobalEnv)}; 
  if(group=="g3"){assign('data_ggseg_g3', data_ggseg, envir=.GlobalEnv)};
}


average_rois_ggseg(data,"all")
average_rois_ggseg(data_group_1,"g1")
average_rois_ggseg(data_group_2,"g2")
average_rois_ggseg(data_group_3,"g3")





### 3. Generate ggseg plots----

### Theme for plots
My_Theme = theme(plot.title = element_text(hjust = 0.5, size = 16, family = "TT Arial"), axis.text.x = element_text(family = "TT Arial"), 
                 axis.text.y = element_text(family = "TT Arial"), axis.title.x = element_text(family = "TT Arial"), 
                 axis.title.y = element_text(family = "TT Arial"), legend.text=element_text(family = "TT Arial"),
                 legend.title =element_text(family = "TT Arial"))


# Function to generate PET data plots
generate_ggseg_plots <- function(data, variable, plot_limits, plot_breaks, title, output_filename, path,labels,decimals) {
  if(decimals=="no"){suvr_labels_format <- scales::number_format(accuracy = 1)}
  if(decimals=="yes"){suvr_labels_format <- scales::number_format(accuracy = 0.01)}
  
  ## Plot with no labels
  if(labels=="no"){
    p <- ggplot(data) + geom_brain(atlas = dk, mapping=aes(fill=get(variable)), color="black", position = position_brain(hemi ~ side)) + 
      scale_fill_gradient(low="yellow", high="red",limits = plot_limits, breaks = plot_breaks, labels = suvr_labels_format, na.value = "grey90") + 
      theme_void() + theme(legend.position = ""); plot(p)
    ggsave(paste(output_filename, ".png", sep = ""), plot = p, units = "in", width = 6, height = 6, dpi = 600, bg = "white", path = path)
    return(p)
  }
  
  if(labels=="yes"){
    p <- ggseg(.data = data, colour = "black", position = "stacked", mapping = aes(fill = get(variable))) +
      scale_fill_gradient(low = "yellow", high = "red", limits = plot_limits, breaks = plot_breaks, labels = suvr_labels_format, na.value = "grey90") + My_Theme + theme_void() + theme(legend.title = element_blank()); plot(p)
    ggsave(paste(output_filename, ".png", sep = ""), plot = p, units = "in", width = 6, height = 6, dpi = 600, bg = "white", path = path)
    return(p)
  }
}


### PET plots
generate_ggseg_plots(data_ggseg, "Ab_SUVR", c(0.8, 1.82), seq(0.8, 2, by = 0.4), "Ab_SUVR", "Ab_SUVR",path,"yes","yes")
generate_ggseg_plots(data_ggseg, "Tau_SUVR", c(0.8, 1.4), seq(0.8, 1.4, by = 0.2), "Tau_SUVR", "Tau_SUVR",path,"yes","yes")

### MEG plots
generate_ggseg_plots(data_ggseg, "Delta", c(0.2, 0.5), seq(0.2, 0.5, by = 0.1), "Delta", "Delta",path,"yes","yes")
generate_ggseg_plots(data_ggseg, "Theta", c(0.15, 0.3), seq(0.15, 0.3, by = 0.05), "Theta", "Theta",path,"yes","yes")
generate_ggseg_plots(data_ggseg, "Alpha", c(0.146, 0.45), seq(0.15, 0.45, by = 0.1), "Alpha", "Alpha",path,"yes","yes")
generate_ggseg_plots(data_ggseg, "Beta", c(0.05, 0.2), seq(0.05, 0.2, by = 0.05), "Beta", "Beta",path,"yes","yes")

### Aperiodics
generate_ggseg_plots(data_ggseg, "Exponent", c(0.6, 1.01), seq(0.6, 1, by = 0.1), "Exponent", "Exponent",path,"yes","yes")
generate_ggseg_plots(data_ggseg, "Offset", c(-1.5, -0.8), seq(-1.4, -0.8, by = 0.2), "Offset", "Offset",path,"yes","yes")






### 4. Extra plots Ab and Tau ROIs----

### Amyloid global idx rois
col_idx=which(colnames(data_ggseg)=="Ab_idx")
data_ggseg[, col_idx][data_ggseg[, col_idx] == 0] <- NA
p <- ggplot(data_ggseg) + geom_brain(atlas = dk, mapping=aes(fill=Ab_idx), color="black", position = position_brain(hemi ~ side)) +
  scale_fill_manual(values=c("#fad473"), na.value = "grey92") + theme_void() +  theme(legend.position = "")
plot(p)
ggsave("Ab_idx_rois.png",plot = p,units="in", width=6, height=6, dpi=600,bg="white",path=path)

### Tau ent roi
col_idx=which(colnames(data_ggseg)=="Tau_ent")
data_ggseg[, col_idx][data_ggseg[, col_idx] == 0] <- NA
p <- ggplot(data_ggseg) +
  geom_brain(atlas = dk, mapping=aes(fill=Tau_ent), color="black", position = position_brain(hemi ~ side)) + 
  scale_fill_manual(values=c("#fad473"), na.value = "grey90") +  theme_void() +  theme(legend.position = "")
plot(p)
ggsave("Tau_ent_rois.png",plot = p,units="in", width=6, height=6, dpi=600,bg="white",path=path)


### Tau meta-ROI rois
col_idx=which(colnames(data_ggseg)=="Tau_mr")
data_ggseg[, col_idx][data_ggseg[, col_idx] == 0] <- NA
p <- ggplot(data_ggseg) +
  geom_brain(atlas = dk, mapping=aes(fill=Tau_mr), color="black", position = position_brain(hemi ~ side)) + 
  scale_fill_manual(values=c("#fad473"), na.value = "grey94") + theme_void() +  theme(legend.position = "")
plot(p)
ggsave("Tau_mr_rois.png",plot = p,units="in", width=6, height=6, dpi=600,bg="white",path=path)




