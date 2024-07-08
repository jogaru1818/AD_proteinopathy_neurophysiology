### Script wrote by Jonathan Gallego Rudolf - PhD at McGill University

# Separate data based on tau threshold (either for entorhinal or tau meta ROI)

# input data in glm spreadsheet format
# select if entorhinal tau or tau meta ROI was used
# Select method to calculate threshold:
# m=1 use threshold derived from Melissa's paper
# m=2 use threshold to separate top quartile from rest of the data
# m=3 use threshold to separate data into tertiles
# m=4 use threshold to separate top and low quartile from rest of the data (marginal effects)
# m=5 use threshold to separate data into quintiles
# m=6 use threshold to separate data into quartiles


f_separate_tau_data <- function(data,tau_roi,m) {
  x <- data[1:nrow(data),]
  if (tau_roi=="entorhinal"){
    col_idx=which(colnames(data)=="high_tau_ent")
  }
  if (tau_roi=="meta_roi"){
    col_idx=which(colnames(data)=="high_tau_meta_roi")
  }
  y <- data[1:nrow(data),col_idx]
  mat <-x[FALSE,]
  if(m==1||m==2){
  ### Low tau
  for (i in 1:nrow(data)){
      if(y[i]==0){
        z <- x[i,]
      }
    if(y[i]==0){
      mat <- (rbind(mat, z))
    }
  }
  data_low_tau=mat
  assign('data_low_tau', data_low_tau, envir=.GlobalEnv)
  #### High tau
  mat2 <-x[FALSE,]
  for (i in 1:nrow(data)){
    if(y[i]==1){
      z <- x[i,]
    }
    if(y[i]==1){
      mat2 <- (rbind(mat2, z))
    }
  }
  data_high_tau=mat2
  assign('data_high_tau', data_high_tau, envir=.GlobalEnv)
  }
  if(m==3||m==4){
    ### Tau level 1
    for (i in 1:nrow(data)){
      if(y[i]==1){
        z <- x[i,]
      }
      if(y[i]==1){
        mat <- (rbind(mat, z))
      }
    }
    data_low_tau=mat
    assign('data_low_tau', data_low_tau, envir=.GlobalEnv)
    
    mat <-x[FALSE,]
    ### Tau level 2
    for (i in 1:nrow(data)){
      if(y[i]==2){
        z <- x[i,]
      }
      if(y[i]==2){
        mat <- (rbind(mat, z))
      }
    }
    data_medium_tau=mat
    assign('data_medium_tau', data_medium_tau, envir=.GlobalEnv)
    
    
    mat <-x[FALSE,]
    ### Tau level 3
    for (i in 1:nrow(data)){
      if(y[i]==3){
        z <- x[i,]
      }
      if(y[i]==3){
        mat <- (rbind(mat, z))
      }
    }
    data_high_tau=mat
    assign('data_high_tau', data_high_tau, envir=.GlobalEnv)
  }
  if(m==5){
    ### Tau level 1
    for (i in 1:nrow(data)){
      if(y[i]==1){
        z <- x[i,]
      }
      if(y[i]==1){
        mat <- (rbind(mat, z))
      }
    }
    data_tau_1=mat
    assign('data_tau_cut_1', data_tau_1, envir=.GlobalEnv)
    
    mat <-x[FALSE,]
    ### Tau level 2
    for (i in 1:nrow(data)){
      if(y[i]==2){
        z <- x[i,]
      }
      if(y[i]==2){
        mat <- (rbind(mat, z))
      }
    }
    data_tau_2=mat
    assign('data_tau_cut_2', data_tau_2, envir=.GlobalEnv)
    
    
    mat <-x[FALSE,]
    ### Tau level 3
    for (i in 1:nrow(data)){
      if(y[i]==3){
        z <- x[i,]
      }
      if(y[i]==3){
        mat <- (rbind(mat, z))
      }
    }
    data_tau_3=mat
    assign('data_tau_cut_3', data_tau_3, envir=.GlobalEnv)
    
    mat <-x[FALSE,]
    ### Tau level 4
    for (i in 1:nrow(data)){
      if(y[i]==4){
        z <- x[i,]
      }
      if(y[i]==4){
        mat <- (rbind(mat, z))
      }
    }
    data_tau_4=mat
    assign('data_tau_cut_4', data_tau_4, envir=.GlobalEnv)
    
    mat <-x[FALSE,]
    ### Tau level 5
    for (i in 1:nrow(data)){
      if(y[i]==5){
        z <- x[i,]
      }
      if(y[i]==5){
        mat <- (rbind(mat, z))
      }
    }
    data_tau_5=mat
    assign('data_tau_cut_5', data_tau_5, envir=.GlobalEnv)
  }
  
  if(m==6){
    ### Tau level 1
    for (i in 1:nrow(data)){
      if(y[i]==1){
        z <- x[i,]
      }
      if(y[i]==1){
        mat <- (rbind(mat, z))
      }
    }
    data_tau_1=mat
    assign('data_tau_cut_1', data_tau_1, envir=.GlobalEnv)
    
    mat <-x[FALSE,]
    ### Tau level 2
    for (i in 1:nrow(data)){
      if(y[i]==2){
        z <- x[i,]
      }
      if(y[i]==2){
        mat <- (rbind(mat, z))
      }
    }
    data_tau_2=mat
    assign('data_tau_cut_2', data_tau_2, envir=.GlobalEnv)
    
    
    mat <-x[FALSE,]
    ### Tau level 3
    for (i in 1:nrow(data)){
      if(y[i]==3){
        z <- x[i,]
      }
      if(y[i]==3){
        mat <- (rbind(mat, z))
      }
    }
    data_tau_3=mat
    assign('data_tau_cut_3', data_tau_3, envir=.GlobalEnv)
    
    mat <-x[FALSE,]
    ### Tau level 4
    for (i in 1:nrow(data)){
      if(y[i]==4){
        z <- x[i,]
      }
      if(y[i]==4){
        mat <- (rbind(mat, z))
      }
    }
    data_tau_4=mat
    assign('data_tau_cut_4', data_tau_4, envir=.GlobalEnv)

  }
}