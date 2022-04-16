####Useful Libraries####

# BiocManager::install("ropls")
# install.packages("PMA")
# BiocManager::install("Cardinal")
# devtools::install_github("HerveAbdi/data4PCCAR") # install data4PCCAR
# install.packages("sparsepca")
# devtools::install_github("derekbeaton/gsvd")
# devtools::install_github("derekbeaton/gpls",subdir = "Package")
# remotes::install_github("coolbutuseless/ggpattern")
# install.packages("geometry")
# install.packages("pals")#color palette for 19 symptoms
# install.packages("rstatix")

library(pls)
library(ggplot2)
library(reshape2)
library(permute)
library(doMC)
library(cowplot)
library(candisc)
library(mixOmics)
library(ropls)
library(TExPosition)
library(PMA)
library(Cardinal)
library(data4PCCAR)
library(sparsepca)
library(GPLS)
library(dplyr)
library(ggpattern)
library(tidyverse)
library(caret)
library(geometry)
library(psych)
library(rstatix)
library(corrplot)
library(pals)

####Useful Functions####
switch_names <- function(dataCol,Rem,Add,separ,side){
  if (side=="L"){
    #Switches anything that is to the left of Rem (including Rem) with Add, separated by 'separ'
    tmp<-substr(dataCol, regexpr(Rem, dataCol)+nchar(Rem), regexpr(Rem, dataCol)+50)  
    switched<-paste(Add, tmp, sep=separ)
    return(switched)
  } else if (side =="R"){
    #Switches anything that is to the right of Rem (including Rem) with Add, separated by 'separ'
    tmp<-substr(dataCol, regexpr(Rem, dataCol)-50, regexpr(Rem, dataCol)-1)
    switched<-paste(tmp, Add, sep=separ)
    return(switched)
  } else {
    "Must specify R or L"
  }
}
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
thresh_conn <- function(data,lbl_col,sub_cols,metrs,crit){
  #data=df containing all data. Should contain the following cols: mean conn  sub metric
  #lbls=col within the df containing the labels
  #crit=decimal showing across how many participants the connection should be present
  # lbls_tmp<-data.frame(unique(conn_data_tmp$conn))#for troubleshoot
  # colnames(lbls_tmp)<-"lbls"
  lbls<-data.frame(unique(lbl_col))
  colnames(lbls)<-"lbls"
  
  #Since all weights are calculated from the same bundles, we don't need to threshold for every metric separately, we just need to check once
  # metrs_tmp<-as.data.frame(metrics_tmp)#for troubleshoot
  # colnames(metrs_tmp)<-"metric"
  metrs<-as.data.frame(metrs)
  colnames(metrs)<-"metric"
  
  # sub_cols_tmp<-data.frame(unique(conn_data_tmp$sub))#for troubleshoot
  
  subs<-data.frame(unique(sub_cols))
  # threshd_conn<-str_conn_feats_QCpass#for troubleshoot
  threshd_conn<-data
  for (i in 1:dim(lbls)[1]){
    tmp<-subset(data,conn==lbls[i,])
    
    for (j in 1:dim(metrs)[1]){
      # # tmp<-subset(str_conn_feats,conn==lbls[i,] & metric==metrs[j,])#for troubleshoot
      # tmp<-subset(data,conn==lbls[i,] & metric==metrs[j,])
      
      #Load all rows for a feature, to get the row names for later use
      # tmp<-subset(threshd_conn,conn==lbls_tmp[i,])#for troubleshoot
      # message(paste("Conn",lbls[i,],sep=" "))
      
      
      # tmp_unique<-subset(str_conn_feats_QCpass,conn==lbls_tmp[i,] & metric==metrs[1,])#for troubleshoot
      # tmp_unique<-subset(data,conn==lbls[i,] & metric==metrs[1,])
      # tmp_unique<-subset(str_conn_feats_QCpass,conn=="10_53" & metric=="fa")#for troubleshoot
      # tmp_unique_tmp<-subset(conn_data_tmp,conn==lbls_tmp[i,] & metric==metrs_tmp[j,])
      
      tmp_unique<-subset(data,conn==lbls[i,] & metric==metrs[j,]) #Pick out all connections for a given metric
      tmp_unique[is.na(tmp_unique$mean)==1,c("mean")]<-0#Within a given set of connections for a given metric, replace NAs by 0 (it makes it easier to do the ratio calculation)
      # dim(tmp)
      rat<-1-(sum(tmp_unique$mean==0)/dim(subs)[1])#Should now work even without the na.rm
      # rat<-1-(sum(tmp_unique$mean==0,na.rm=T)/dim(subs)[1])#Here is the issue, I'm doing na.rm, which just skips NaNs. That's why I ended up having NaNs later on. A threshold of 100 means I should remove these conns altogether.
      # threshd_conn[rownames(tmp),c("mean")]<-if(rat<1,yes = 0,no = str_conn_feats[rownames(tmp),c("mean")])#for troubleshoot
      
      # threshd_conn[rownames(tmp),c("mean")]<-ifelse(rat<crit,yes = 0,no = data[rownames(tmp),c("mean")])
      if (rat>=crit){
        message(paste("Conn",lbls[i,],"metric",metrs[j,],"meets criterion of",crit,"(",rat,")","keeping conn",sep=" "))
      } else if (rat<crit){
        message(paste("Conn",lbls[i,],"metric",metrs[j,],"does not meet criterion of",crit,"(",rat,")","setting to 0",sep=" "))
        # threshd_conn[rownames(tmp),c("mean")]<-0#CHANGING TO TURN TO NA, TO AVOID ISSUES WITH KEEPING 0s WHEN USING LESS THAN 1 thresh
        # threshd_conn[rownames(tmp),c("mean")]<-NA
        threshd_conn[rownames(tmp_unique),c("mean")]<-NA #Here put NAs for the entire set of connections for a given metric
        
        
      }
    }
  }
  message("WARNING: Keep in mind that only conns that don't meet criterion are put to NA, missing conns within conns that met criterion (i.e.: when t<1.00) are listed as 0 and should be imputed")
  return(threshd_conn)
}
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
reg_out_site <- function(data,lbl_col,sub_cols,metr_col,nuisance_variables){
  #data=df containing all data. Should contain the following cols: mean conn  sub metric
  #lbls=col within the df containing the labels
  #crit=decimal showing across how many participants the connection should be present
  # lbls<-data.frame(unique(t_str_conn_feats_trim_site$conn))#for troubleshoot
  
  # nuisance_vars_tmp = c("sex","ehi_y_ss_scoreb","mri_info_deviceserialnumber","pubertal_stage") # For troubleshoot
  # residuals_tmp<-reg_out_site(t_conn_data_trim_site_tmp,lbl_col = t_conn_data_trim_site_tmp$conn,sub_cols = t_conn_data_trim_site_tmp$sub,
  #                         metr_col = t_conn_data_trim_site_tmp$metric,nuisance_variables =t_conn_data_trim_site_tmp[,c(nuisance_vars_tmp)])# For troubleshoot
  
  
  lbls<-data.frame(unique(lbl_col))
  colnames(lbls)<-"lbls"
  
  # data_tmp<-t_conn_data_trim_site_tmp#for troubleshoot
  # lbls_tmp<-data.frame(unique(t_conn_data_trim_site_tmp$conn))#for troubleshoot
  # colnames(lbls_tmp)<-"lbls"#for troubleshoot
  
  
  metrs<-data.frame(unique(metr_col))
  colnames(metrs)<-"metric"
  
  # metrs_tmp<-data.frame(unique(t_conn_data_trim_site_tmp$metric)) #For troubleshoot
  # colnames(metrs_tmp)<-"metric" #For troubleshoot
  
  # subs<-data.frame(unique(t_str_conn_feats_trim_site$sub))#for troubleshoot
  subs<-data.frame(unique(sub_cols))
  tot_feats=dim(lbls)[1]*dim(metrs)[1]
  cnt=0
  data$residuals<-NA
  
  # subs_tmp<-data.frame(unique(t_conn_data_trim_site_tmp$sub))#for troubleshoot
  # tot_feats_tmp=dim(lbls_tmp)[1]*dim(metrs_tmp)[1]#for troubleshoot
  # cnt_tmp=0#for troubleshoot
  # data_tmp$residuals<-NA#for troubleshoot
  
  for (i in 1:dim(lbls)[1]){
    for (j in 1:dim(metrs)[1]){
      cnt=cnt+1
      message(paste("Conn",lbls[i,],"metric",metrs[j,],"[",cnt,"/",tot_feats,"]",sep=" "))
      # tmp<-subset(t_str_conn_feats_trim_site,conn==lbls[i,] & metric==metrs[j,])#for troubleshoot
      tmp<-subset(data,conn==lbls[i,] & metric==metrs[j,])
      
      tmp_tmp<-subset(data_tmp,conn==lbls_tmp[i,] & metric==metrs_tmp[j,])#For troubleshoot
      if (dim(tmp)[1]!=0){
        
        # #example of the command to run this function
        # reg_out_site(t_conn_data_trim_site,lbl_col = t_conn_data_trim_site$conn,sub_cols = t_conn_data_trim_site$sub,
        #              metr_col = t_conn_data_trim_site$metric,nuisance_variables = nuisance_dat[,c(nuisance_vars)])
        
        
        lmod<-lm(mean ~ .,data = tmp[,c("mean",nuisance_variables)]) #Select the data to be just the mean and nuisance vars, and then regress out nuisance vars (the '.' in the formula means "use all other vars in data")
        
        # nuisance_variables_tmp<-nuisance_vars_tmp
        # lmod<-lm(mean ~ .,data = tmp_tmp[,c("mean",nuisance_variables_tmp )]) #For troubleshoot
        # summary(lmod)
        data[rownames(tmp),c("residuals")]<-residuals(lmod)
      }
    }
  }
  return(data)
}
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
scale_by_metric_conn <- function(data,metr_col){
  # metrs_tmp<-data.frame(unique(t_str_conn_feats_trim_site$metric))#for troubleshoot
  metrs<-data.frame(unique(metr_col))
  colnames(metrs)<-"metric"
  
  data$scaled_mean<-NA
  
  for (j in 1:dim(metrs)[1]){
    tmp<-subset(data,metric==metrs[j,])
    # tmptmp<-subset(t_str_conn_feats_trim_site,metric==metrs_tmp[j,]) #for troubleshoot
    tmp<-data.frame(tmp)
    # tmptmp<-data.frame(tmptmp) #for troubleshoot
    conns<-data.frame(unique(tmp$conn))
    colnames(conns)<-"conn"
    # conns_tmp<-data.frame(unique(tmptmp$conn)) # for troubleshoot
    for (i in 1:dim(conns)[1]){
      message(paste("Conn",conns[i,],"metric",metrs[j,],sep=" "))
      tmp_conn<-subset(tmp,conn==conns[i,],c("residuals"))
      # colnames(tmp_conn)<-c("sub","mean","conn","metric","keep","mri_info_deviceserialnumber","residuals","scaled_mean")
      # tmp_conn_tmp<-subset(tmptmp,conn==conns_tmp[i,],c("residuals")) #for troubleshoot
      # dim(tmp_conn_tmp)
      tmp_conn<-data.frame(tmp_conn)
      # tmp_conn_tmp<-data.frame(tmp_conn_tmp) # for troubleshoot
      data[rownames(tmp_conn),c("scaled_mean")]<-scale(as.vector(tmp_conn$residuals))
      # t_str_conn_feats_trim_site[rownames(tmp_conn_tmp),c("scaled_mean")]<-scale(as.vector(tmp_conn_tmp$residuals)) #for troubleshoot
    }
  }
  return(data$scaled_mean)
}
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
reshape_conn_mat <- function(data,lbl_names,mtrName){
  # data<-ld_fa_lf1
  print("This function requires a column with the connection numbers named 'conn', which it splits up into separate ROIs. It also needs lbls, a data frame that contains lbls, num, and full_name")
  sym_mat<-matrix(data=NA,nrow=dim(lbl_names)[1],ncol=dim(lbl_names)[1])
  data$lbl1<-NA
  data$lbl2<-NA
  data[,c("lbl1")]<-substr(data$conn,1,regexpr("_",data$conn)-1)
  data[,c("lbl2")]<-substr(data$conn,regexpr("_",data$conn)+1,regexpr("_",data$conn)+5)
  for (i in 1:dim(lbl_names)[1]){
    for (j in 1:dim(lbl_names)[1]){
      if (i==j){
        sym_mat[i,j]<-0
      } else if (i<j){
        tmp<-data[data$lbl1==lbl_names[i,c("num")]&data$lbl2==lbl_names[j,c("num")],colnames(data)==mtrName]
        # tmp<-data[data$lbl1==lbls[i,c("num")]&data$lbl2==lbls[j,c("num")],colnames(data)=="Loadings"]#For troubleshoot
        if (length(tmp)==0){
          sym_mat[i,j]<-0
          sym_mat[j,i]<-0
        } else {
          sym_mat[i,j]<-tmp
          sym_mat[j,i]<-tmp
        }
        
        
      }
      
      
    }
  }
  return(sym_mat)
}
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
apply_pls <- function(X,Y,ncomps,scaleflag,permflag,outflag){
  # data:the mydata set up such that xx is predictors, yy is outcome
  # ncomps: number of components to fit
  # scaleflag: option to decide if we want to scale the input data (both X and Y) or not
  # permflag: option to determine if we want to obtain the full output from tepPLS (permflag=F) or individual outputs (permflag=T).
  # The second option is useful for permutation testing, which uses this function.
  if (scaleflag==T){
    pls_model<-tepPLS(X,Y,k = ncomps,scale1 = TRUE,scale2 = TRUE,graphs = FALSE)
  } else if (scaleflag==F){
    pls_model<-tepPLS(X,Y,k = ncomps,scale1 = FALSE,scale2 = FALSE,graphs = FALSE)
  }
  #Here the function should behave differently if it's used for permutation or to get the actual scores
  if (permflag==F){
    return(pls_model)
  }else if (permflag==T){
    message("Specify output flag: 1-Eigenvectors, 2-PercentCovarianceExplained, 3-XSaliences, 4-YSaliences, 5-Xscores, 6-Yscores")
    if (outflag==1){
      return(pls_model$TExPosition.Data$eigs)
    } else if (outflag==2){
      return(pls_model$TExPosition.Data$t)
    } else if (outflag==3){
      return(pls_model$TExPosition.Data$lx)
    } else if (outflag==4){
      return(pls_model$TExPosition.Data$ly)
    } else if (outflag==5){
      return(pls_model$TExPosition.Data$pdq$p)
    } else if (outflag==6){
      return(pls_model$TExPosition.Data$pdq$q)
    }
    
  }
  
}
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
pre_proc <- function(conn_data,symp_data,nuisance_dat,metrics,nuisance_vars,nuisance_vars_symps,thresh,tr_ratio,out_dir,resume_thresh,thresh_file,resume_regout,regout_file,resume_scale,scale_file){
  # conn_data:the connectivity data after QC removal
  # symp_data:the symptoms data
  # site_dat: site data
  # thresh: threshold 
  # tr_ratio: training/testing group ratio
  # This function saves the outputs of some of the more time-consuming steps so that, if one needs to rerun it, it can resume from one of those steps.
  
  metrics_df<-as.data.frame(metrics)#Turning metrics into a data frame (will be used further down)
  colnames(metrics_df)<-"metric"
  
  conn_data<-conn_data[conn_data$metric %in% metrics,]#Restricting the conn_data to the metrics the user inputed
  conn_data$metric<-factor(conn_data$metric, levels=metrics)
  
  message("We are starting with ",dim(data.frame(unique(conn_data$sub)))[1]," subs, processing metrics ",metrics,sep=" ")
  
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@THRESHOLDING@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  if (resume_thresh==FALSE){
    # t_str_conn_feats<-thresh_conn(str_conn_feats_QCpass,str_conn_feats_QCpass$conn,str_conn_feats_QCpass$sub,str_conn_feats_QCpass$metric,1)#For troubleshoot
    message("Thresholding conns")
    t_conn_data<-thresh_conn(conn_data,conn_data$conn,conn_data$sub,metrics,thresh)
    
    message("Flagging thresholded conns")
    t_conn_data$keep<-0
    t_conn_data[is.na(t_conn_data$mean)==0,c("keep")]<-1#HERE I FLAG CONNS TO REMOVE. I'M FLAGGING THEM BY LOOKING AT THE ONES THAT HAVE NANS. REMEMBER, THOSE THAT ARE MISSING BUT FROM CONNS THAT PASSED THRESHOLD, THEY ARE 0, NOT NANS.
    
    message(paste("After thresholding with t=",thresh,", ",sum(t_conn_data$keep),"/",dim(t_conn_data)[1]," features were kept"))
    
    write.table(t_conn_data,paste(out_dir,"/conn_features_thresh_",thresh,sep=""),quote = F) ##REMEMBER: CHANGE THE OUTPUT NAME HERE ONCE YOU ARE SURE THE SCRIPT WORKS##
    
    t_conn_data_trim<-t_conn_data[t_conn_data$keep==1,]#Trimming first, since we don't want to use missing data to regress out site. HERE IS WHERE I'M REMOVING ROWS BELONGING TO CONNS THAT DID NOT PASS THRESHOLD.
    message(paste("There are ",dim(t_conn_data)[1]-dim(t_conn_data_trim)[1]," features that were removed after thresholding"))
    
  } else if (resume_thresh==TRUE){
    message("Resuming from threshold step")
    t_conn_data<-read.table(thresh_file, sep=' ', header=T)
    t_conn_data<-data.frame(t_conn_data)
    
    t_conn_data_trim<-t_conn_data[t_conn_data$keep==1,]#Trimming first, since we don't want to use missing data to regress out site. 
    message(paste("There are ",dim(t_conn_data)[1]-dim(t_conn_data_trim)[1]," features that were removed after thresholding"))
    
    # 
    # message("Resuming from threshold step")
    # t_conn_data<-read.table('/Users/Guido/Desktop/Projects/PCHStudy/Data/conn_features_thresh_0.9_newesttest.txt', sep=' ', header=T)
    # t_conn_data<-data.frame(t_conn_data)
    # 
    # t_conn_data_trim<-t_conn_data[t_conn_data$keep==1,]#Trimming first, since we don't want to use missing data to regress out site. 
    # message(paste("There are ",dim(t_conn_data)[1]-dim(t_conn_data_trim)[1]," features that were removed after thresholding"))
    # 
    # 
  }
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@SPLITTING DATASETS@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  message("Splitting dataset into training and testing sets")
  # colnames(t_conn_data_trim_site_tmp)<-c("subjectkey","mean","conn","metric","keep","mri_info_deviceserialnumber","residuals")
  subs<-data.frame(unique(t_conn_data_trim$sub))
  # dim(subs)
  colnames(subs)<-"subjectkey"
  smp_size=floor(tr_ratio*nrow(subs))
  set.seed(123)
  train_idx<-sample(seq_len(nrow(subs)),size=smp_size)
  subs$grp<-2
  subs[train_idx,c("grp")]<-1
  
  t_conn_data_trim<-merge(t_conn_data_trim,subs,by.x = "sub",by.y="subjectkey")
  
  ##HERE IS WHERE THE DATA IS ACTUALLY SPLIT IN TWO SETS##
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@SPLIT CONN DATA: TRAINING SET@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  t_conn_data_trim_train<-t_conn_data_trim[t_conn_data_trim$grp==1,]
  message(paste("There are ",dim(as.data.frame(unique(t_conn_data_trim_train$sub)))[1],"subs in the train set, ",
                dim(as.data.frame(unique(t_conn_data_trim$sub)))[1]-dim(as.data.frame(unique(t_conn_data_trim_train$sub)))[1],"in the test set"))
  
  #Changing data structure#
  message("Changing data structure - long to wide")
  t_conn_data_trim_train_halfwide<-reshape(t_conn_data_trim_train[,c("sub","mean","conn","metric")],idvar = c("sub","metric"),timevar="conn",direction="wide")
  
  t_conn_data_trim_train_wide<-reshape(t_conn_data_trim_train_halfwide,idvar = "sub",timevar="metric",direction="wide")
  message(paste("Data in wide has the following dimensions:",dim(t_conn_data_trim_train_wide)[1],dim(t_conn_data_trim_train_wide)[2],sep=" "))
  
  tmp_full_miss<-t_conn_data_trim_train_wide[,colSums(is.na(t_conn_data_trim_train_wide))==
                                               dim(t_conn_data_trim_train_wide)[1]]#Here I'm calculating how many NAs are in a column. If the number of NAs in a column == the dim[2] of the table,
  # it means the entire column is filled with NAs. RECALL THAT I REMOVED SUBTHRESHOLD CONNS ALREADY. I'M HAVING TO DO THIS AGAIN BECAUSE SOME METRICS HAD CONNS THAT SURVIVED THRESHOLD WHEREAS SOME DIDN'T.
  # WHEN I RESHAPE TO WIDE, IT'S EXPECTING EVERYTHING TO BE REPEATED (E.G.: IF I HAVE METRIC=FA, CONN=10_11, IT'S EXPECTING ME TO HAVE METRIC=NUFO, CONN=10_11. IF IT'S MISSING, IT WILL CREATE IT IN THE WIDE TABLE AND FILL IT WITH NANS.
  # I NEED TO FLAG THESE AND REMOVE THEM.
  
  message("Due to lack of some conns in specific metrics but not others, ",dim(tmp_full_miss)[2])
  message("Conns that were introduced by error and now removed: ",colnames(tmp_full_miss))
  t_conn_data_trim_train_wide<-t_conn_data_trim_train_wide[,!colSums(is.na(t_conn_data_trim_train_wide))==dim(t_conn_data_trim_train_wide)[1]] #SO HERE I'M KEEPING THE COLUMNS THAT WEREN'T INTRODUCED BY MISTAKE DURING THE RESHAPING.
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  
  
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@IMPUTING CONN DATA: TRAINING SET@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  if (thresh<1){
    message("Imputing missing conn data - training set")
    
    t_conn_data_trim_train_imput<-t_conn_data_trim_train_wide[,-1]#Needs to be the full wide table. First remove the subjects.
    
    for (i in 1:dim(t_conn_data_trim_train_imput)[2]){
      idx_miss<-rownames(t_conn_data_trim_train_imput[is.na(t_conn_data_trim_train_imput[,i])==1,])
      idx_miss<-data.frame(idx_miss)
      idx_comp<-rownames(t_conn_data_trim_train_imput[is.na(t_conn_data_trim_train_imput[,i])==0,])
      idx_comp<-data.frame(idx_comp)
      if (dim(idx_miss)[1]!=0 && dim(idx_miss)[1]<dim(idx_comp)[1]){
        imp_size=dim(as.data.frame(idx_miss))[1]
        set.seed(i)
        idx_rep<-sample(as.vector(idx_comp$idx_comp),size=imp_size)
        t_conn_data_trim_train_imput[rownames(t_conn_data_trim_train_imput)==as.numeric(as.character(idx_miss$idx_miss)),i]<-t_conn_data_trim_train_imput[idx_rep,i]
      }
    }
  } else if (thresh==1){
    t_conn_data_trim_train_imput<-t_conn_data_trim_train_wide
  }
  
  
  t_conn_data_trim_train_imput$subjectkey<-t_conn_data_trim_train_wide$sub#Add back subjectkey to perform the merge of the nuisance vars
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#@@@@@@@@@@@@@
  
  
  ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@SPLIT SYMPTOM DATA: TRAINING SET#@@@@@@@@@@@@@@@@@@@@@@@
  message("Splitting symptoms into train/test sets")
  symptoms_train<-symp_data[as.character(symp_data$subjectkey) %in% t_conn_data_trim_train_imput$subjectkey,] # filter symptom data by looking at subs that are in training set
  
  ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@IMPUTING MISSING SYMPTOM DATA: TRAINING SET#@@@@@@@@@@@@@@@@@@@@@@@
  message("Imputing missing symp data - training set")
  for (i in 1:dim(symptoms_train)[2]){
    idx_miss<-rownames(symptoms_train[is.na(symptoms_train[,i])==1,])
    if (length(idx_miss)!=0){
      imp_size=dim(as.data.frame(idx_miss))[1]
      # set.seed(123) # Should we set a seed?
      idx_comp<-rownames(symptoms_train[is.na(symptoms_train[,i])==0,])
      set.seed(i)
      idx_rep<-sample(idx_comp,size=imp_size)
      symptoms_train[idx_miss,i]<-symptoms_train[idx_rep,i]
    }
  }
  
  ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#@@@@@@@@@@@@@@@@@@@@@@@
  
  ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@IMPUTING MISSING NUISANCE DATA: TRAINING SET#@@@@@@@@@@@@@@@@@@@@@@@
  nuisance_dat_train<-nuisance_dat[as.character(nuisance_dat$subjectkey) %in% t_conn_data_trim_train_imput$subjectkey,] # filter symptom data by looking at subs that are in training set
  message("Imputing missing nuisance covariates - training set")
  for (i in 1:dim(nuisance_dat_train)[2]){
    idx_miss<-rownames(nuisance_dat_train[is.na(nuisance_dat_train[,i])==1,])
    if (length(idx_miss)!=0){
      imp_size=dim(as.data.frame(idx_miss))[1]
      # set.seed(123) # Should we set a seed?
      idx_comp<-rownames(nuisance_dat_train[is.na(nuisance_dat_train[,i])==0,])
      set.seed(i)
      idx_rep<-sample(idx_comp,size=imp_size)
      nuisance_dat_train[idx_miss,i]<-nuisance_dat_train[idx_rep,i]
    }
  }
  
  ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
  ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@REGRESSING OUT NUISANCE VARS: TRAINING SET#@@@@@@@@@@@@@@@@@@@@@@@
  if (resume_regout==FALSE){
    message("Regressing out nuisance vars - step 1: merging nuisance vars")
    
    t_conn_data_trim_train_nuisance<-merge(t_conn_data_trim_train_imput,nuisance_dat_train)#By merging like this I'm just keeping the nuisance vars that are in train set
    
    message("There are",dim(data.frame(unique(t_conn_data_trim_train_nuisance$subjectkey)))[1]," subs (",
            dim(data.frame(unique(t_conn_data_trim_train_imput$subjectkey)))[1]-dim(data.frame(unique(t_conn_data_trim_train_nuisance$subjectkey)))[1],"lost) after incorporating nuisance info")
    
    ##DO THE SAME FOR SYMPTOMS##
    symptoms_train_nuisance<-merge(symptoms_train,nuisance_dat_train)#By merging like this I'm just keeping the nuisance vars that are in train set
    
    message("There are",dim(data.frame(unique(symptoms_train_nuisance$subjectkey)))[1]," subs (",
            dim(data.frame(unique(symptoms_train$subjectkey)))[1]-dim(data.frame(unique(symptoms_train_nuisance$subjectkey)))[1],"lost) after incorporating nuisance info")
    
    
    
    ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@REGRESSING OUT NUISANCE VARS: TRAINING SET#@@@@@@@@@@@@@@@@@@@@@@@
    message("Regressing out site - step 2: regress out site")
    message("Brain")
    t_conn_data_trim_train_regout<-matrix(data = NA,nrow = dim(t_conn_data_trim_train_nuisance)[1],ncol = dim(t_conn_data_trim_train_nuisance[,-which(colnames(t_conn_data_trim_train_nuisance) %in%
                                                                                                                                                        c("subjectkey","interview_age","t_since_inj",nuisance_vars))])[2])
    t_conn_data_trim_train_regout<-as.data.frame(t_conn_data_trim_train_regout)
    dim(t_conn_data_trim_train_regout)
    colnames(t_conn_data_trim_train_regout)<-colnames(t_conn_data_trim_train_nuisance[,-which(colnames(t_conn_data_trim_train_nuisance) %in% c("subjectkey","interview_age","t_since_inj",nuisance_vars))])
    for (i in 1:dim(t_conn_data_trim_train_regout)[2]){
      #Loop through the data frame but only the connections
      message("Regressing out ", nuisance_vars," for ",colnames(t_conn_data_trim_train_regout)[i]," [",i,"/",dim(t_conn_data_trim_train_regout)[2],"]")
      in_dat<-t_conn_data_trim_train_nuisance[,c(colnames(t_conn_data_trim_train_regout)[i],nuisance_vars)]
      colnames(in_dat)<-c("mean",nuisance_vars)
      # Iteratively change your input data for the model. Since the col names in the output table have been preset, you can use those names to select the input conn. I have to do this because you can't input variable names
      #in lm apparently. But I don't want to write out the nuisance vars instead we want to change them. So by doing it this way, we can use variable nuisance vars. The '.' in the formula means: use all other columns.
      lmod<-lm(mean ~ .,data = in_dat)#
      t_conn_data_trim_train_regout[,i]<-residuals(lmod)
    }
    
    t_conn_data_trim_train_regout$subjectkey<-t_conn_data_trim_train_nuisance$subjectkey#Add back subjectkey
    t_conn_data_trim_train_regout<-as.data.frame(t_conn_data_trim_train_regout)
    
    # write.table(t_conn_data_trim_train_regout,paste(out_dir,"conn_features_thresh_",thresh,"_regd_train.txt",sep=""),quote = F)
    
    message("Symptoms")
    symptoms_train_regout<-matrix(data = NA,nrow = dim(symptoms_train_nuisance)[1],ncol = dim(symptoms_train_nuisance[,-which(colnames(symptoms_train_nuisance) %in%
                                                                                                                                c("subjectkey","interview_age","t_since_inj","mri_info_deviceserialnumber",nuisance_vars_symps))])[2])
    symptoms_train_regout<-as.data.frame(symptoms_train_regout)
    dim(symptoms_train_regout)
    #NOTE: notice that I'm regressing out different things. For symptoms, I'm not regressing out mri site. So that's why I have nuisance_vars_symps
    colnames(symptoms_train_regout)<-colnames(symptoms_train_nuisance[,-which(colnames(symptoms_train_nuisance) %in% c("subjectkey","interview_age","t_since_inj","mri_info_deviceserialnumber",nuisance_vars_symps))])
    for (i in 1:dim(symptoms_train_regout)[2]){
      #Loop through the data frame but only the connections
      message("Regressing out ", nuisance_vars_symps," for ",colnames(symptoms_train_regout)[i]," [",i,"/",dim(symptoms_train_regout)[2],"]")
      in_dat<-symptoms_train_nuisance[,c(colnames(symptoms_train_regout)[i],nuisance_vars_symps)]
      colnames(in_dat)<-c("mean",nuisance_vars_symps)
      # Iteratively change your input data for the model. Since the col names in the output table have been preset, you can use those names to select the input conn. I have to do this because you can't input variable names
      #in lm apparently. But I don't want to write out the nuisance vars instead we want to change them. So by doing it this way, we can use variable nuisance vars. The '.' in the formula means: use all other columns.
      lmod<-lm(mean ~ .,data = in_dat)#
      symptoms_train_regout[,i]<-residuals(lmod)
    }
    
    symptoms_train_regout$subjectkey<-symptoms_train_nuisance$subjectkey#Add back subjectkey
    symptoms_train_regout<-as.data.frame(symptoms_train_regout)
    
    write.table(symptoms_train_regout,paste(out_dir,"symp_features_thresh_",thresh,"_regd_train.txt",sep=""),quote = F)
    
    
    
    
  } else if (resume_regout==TRUE){
    message("Resuming from regressing out site step")
    
    t_conn_data_trim_train_regout<-read.table(regout_file, sep=' ', header=T)
    t_conn_data_trim_train_regout<-data.frame(t_conn_data_trim_train_regout)
    
    #NOTE: these lines below are wrong since I'm importing the same file as above. It's just a placeholder, most likely I will be removing this option of resuming here.
    symptoms_train_regout<-read.table(regout_file, sep=' ', header=T)
    symptoms_train_regout<-data.frame(t_conn_data_trim_train_regout)
    
  }
  ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
  ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ADDITIONAL DATA TRANSFORMATIONS: TRAINING SET#@@@@@@@@@@@@@@@@@@@@@@@
  #Check if there are any NANs
  #And for the moment, set the entire column to 0, as another way to trim, note that there might be mistakes so this trim needs to be double checked
  message("Checking that NaNs have been removed")
  tmp_wide_off<-t_conn_data_trim_train_regout
  for (i in 1:dim(tmp_wide_off)[2]){
    if (sum(is.na(tmp_wide_off[,i]))>0){
      sub<-tmp_wide_off[is.na(tmp_wide_off[,i])==T,c("subjectkey")]
      col<-colnames(t_conn_data_trim_regout)[i]
      message("Found missing value in col ",col," for sub ",sub," remember to threshold again!")
      # tmp_wide[is.na(tmp_wide[,i])==T,i]<-0
      # tmp_wide[,i]<-NA
    }
  }
  if (sum(is.na(tmp_wide_off))==0){
    message("You got rid of all NaNs, you may proceed ...")
  }
  
  message("Quick sanity check: Are symptoms and connectivity features the same dims?")
  message("Conn feats: ",dim(t_conn_data_trim_train_regout)[1],"Symps:",dim(symptoms_train_regout)[1])
  ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#@@@@@@@@@@@@@@@@@@@@@@@
  
  
  ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@SETTING UP DF OUTPUT: TRAINING SET#@@@@@@@@@@@@@@@@@@@@@@@
  message("Setting up DF output - Training set")
  ##I"M CLOSE. THE THING IS, IF I DON'T SPLIT THE METRICS HERE I'll have to do it outside the function and that will add code
  conn_feat_names<-as.data.frame(colnames(t_conn_data_trim_train_regout))
  colnames(conn_feat_names)<-"names"
  split_conn_feats_train<-vector(mode = "list", length = dim(metrics_df)[1])
  for(i in 1:dim(metrics_df)[1]){
    mtr_find<-apply(conn_feat_names,1,function(x){regexpr(metrics_df[i,],x)})
    split_conn_feats_train[[i]]<-t_conn_data_trim_train_regout[,mtr_find>0]
  }
  
  names(split_conn_feats_train)<-metrics
  ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#@@@@@@@@@@@@@@@@@@@@@@@
  
  
  
  ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@SPLIT CONN DATA: TESTING SET#@@@@@@@@@@@@@@@@@@@@@@@
  t_conn_data_trim_test<-t_conn_data_trim[t_conn_data_trim$grp==2,]
  
  #Changing data structure#
  message("Changing data structure - long to wide") 
  #NOTE: had to proceed with this two-step reshape because I had multiple repeating units
  t_conn_data_trim_test_halfwide<-reshape(t_conn_data_trim_test[,c("sub","mean","conn","metric")],idvar = c("sub","metric"),timevar="conn",direction="wide")
  
  t_conn_data_trim_test_wide<-reshape(t_conn_data_trim_test_halfwide,idvar = "sub",timevar="metric",direction="wide")
  
  message(paste("Data in wide has the following dimensions:",dim(t_conn_data_trim_test_wide)[1],dim(t_conn_data_trim_test_wide)[2],sep=" "))
  
  tmp_full_miss<-t_conn_data_trim_test_wide[,colSums(is.na(t_conn_data_trim_test_wide))==
                                              dim(t_conn_data_trim_test_wide)[1]]#Here I'm calculating how many NAs are in a column. If the number of NAs in a column == the dim[2] of the table, 
  # it means the entire column is filled with NAs. RECALL THAT I REMOVED SUBTHRESHOLD CONNS ALREADY. I'M HAVING TO DO THIS AGAIN BECAUSE SOME METRICS HAD CONNS THAT SURVIVED THRESHOLD WHEREAS SOME DIDN'T. 
  # WHEN I RESHAPE TO WIDE, IT'S EXPECTING EVERYTHING TO BE REPEATED (E.G.: IF I HAVE METRIC=FA, CONN=10_11, IT'S EXPECTING ME TO HAVE METRIC=NUFO, CONN=10_11. IF IT'S MISSING, IT WILL CREATE IT IN THE WIDE TABLE AND FILL IT WITH NANS.
  # I NEED TO FLAG THESE AND REMOVE THEM.
  
  message("Due to lack of some conns in specific metrics but not others, ",dim(tmp_full_miss)[2])
  message("Conns that were introduced by error and now removed: ",colnames(tmp_full_miss))
  t_conn_data_trim_test_wide<-t_conn_data_trim_test_wide[,!colSums(is.na(t_conn_data_trim_test_wide))==dim(t_conn_data_trim_test_wide)[1]] #SO HERE I'M KEEPING THE COLUMNS THAT WEREN'T INTRODUCED BY MISTAKE DURING THE RESHAPING.
  
  ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#@@@@@@@@@@@@@@@@@@@@@@@
  
  ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@IMPUTING MISSING CONN DATA: TESTING SET#@@@@@@@@@@@@@@@@@@@@@@@  
  if (thresh<1){
    message("Imputing missing data - testing set")
    
    t_conn_data_trim_test_imput<-t_conn_data_trim_test_wide[,-1]#Needs to be the full wide table. First remove the subjects.
    
    for (i in 1:dim(t_conn_data_trim_test_imput)[2]){
      idx_miss<-rownames(t_conn_data_trim_test_imput[is.na(t_conn_data_trim_test_imput[,i])==1,])
      idx_miss<-data.frame(idx_miss)
      idx_comp<-rownames(t_conn_data_trim_test_imput[is.na(t_conn_data_trim_test_imput[,i])==0,])
      idx_comp<-data.frame(idx_comp)
      if (dim(idx_miss)[1]!=0 && dim(idx_miss)[1]<dim(idx_comp)[1]){
        imp_size=dim(as.data.frame(idx_miss))[1]
        
        set.seed(i)
        idx_rep<-sample(as.vector(idx_comp$idx_comp),size=imp_size)
        t_conn_data_trim_test_imput[rownames(t_conn_data_trim_test_imput)==as.numeric(as.character(idx_miss$idx_miss)),i]<-t_conn_data_trim_test_imput[idx_rep,i]
      }
    }
    
  } else if (thresh==1){
    t_conn_data_trim_test_imput<-t_conn_data_trim_test_wide
  }
  
  
  t_conn_data_trim_test_imput$subjectkey<-t_conn_data_trim_test_wide$sub#Add back subjectkey to perform the merge of the nuisance vars
  ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#@@@@@@@@@@@@@@@@@@@@@@@
  
  
  ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@SPLITING SYMPTOM DATA: TESTING SET#@@@@@@@@@@@@@@@@@@@@@@@  
  message("Splitting symptoms into train/test sets")  
  symptoms_test<-symp_data[as.character(symp_data$subjectkey) %in% t_conn_data_trim_test_imput$subjectkey,] # filter symptom data by looking at subs that are in training set
  
  # symptoms_test<-symptoms_mTBI[as.character(symptoms_mTBI$subjectkey) %in% t_conn_data_trim_test_imput$subjectkey,] # filter symptom data by looking at subs that are in training set
  ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#@@@@@@@@@@@@@@@@@@@@@@@
  
  ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@IMPUTING MISSING SYMPTOM DATA: TESTING SET#@@@@@@@@@@@@@@@@@@@@@@@  
  for (i in 1:dim(symptoms_test)[2]){
    idx_miss<-rownames(symptoms_test[is.na(symptoms_test[,i])==1,])
    if (length(idx_miss)!=0){
      imp_size=dim(as.data.frame(idx_miss))[1]
      # set.seed(123) # Should we set a seed?
      idx_comp<-rownames(symptoms_test[is.na(symptoms_test[,i])==0,])
      set.seed(i)
      idx_rep<-sample(idx_comp,size=imp_size)
      symptoms_test[idx_miss,i]<-symptoms_test[idx_rep,i]
    }
  }
  ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#@@@@@@@@@@@@@@@@@@@@@@@
  
  ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@IMPUTING MISSING NUISANCE DATA: TESTING SET#@@@@@@@@@@@@@@@@@@@@@@@  
  nuisance_dat_test<-nuisance_dat[as.character(nuisance_dat$subjectkey) %in% t_conn_data_trim_test_imput$subjectkey,] # filter symptom data by looking at subs that are in training set
  
  # nuisance_dat_test<-nuisance_covs[as.character(nuisance_covs$subjectkey) %in% t_conn_data_trim_test_imput$subjectkey,] # filter symptom data by looking at subs that are in training set
  message("Imputing missing nuisance covariates")
  for (i in 1:dim(nuisance_dat_test)[2]){
    idx_miss<-rownames(nuisance_dat_test[is.na(nuisance_dat_test[,i])==1,])
    if (length(idx_miss)!=0){
      imp_size=dim(as.data.frame(idx_miss))[1]
      # set.seed(123) # Should we set a seed?
      idx_comp<-rownames(nuisance_dat_test[is.na(nuisance_dat_test[,i])==0,])
      set.seed(i)
      idx_rep<-sample(idx_comp,size=imp_size)
      nuisance_dat_test[idx_miss,i]<-nuisance_dat_test[idx_rep,i]
    }
  }
  
  ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#@@@@@@@@@@@@@@@@@@@@@@@
  
  
  ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@REGRESSING OUT NUISANCE DATA: TESTING SET#@@@@@@@@@@@@@@@@@@@@@@@  
  
  if (resume_regout==FALSE){
    message("Regressing out nuisance vars - step 1: merging nuisance vars")
    t_conn_data_trim_test_nuisance<-merge(t_conn_data_trim_test_imput,nuisance_dat_test)
    
    message("There are",dim(data.frame(unique(t_conn_data_trim_test_nuisance$subjectkey)))[1]," subs (",
            dim(data.frame(unique(t_conn_data_trim_test_imput$subjectkey)))[1]-dim(data.frame(unique(t_conn_data_trim_test_nuisance$subjectkey)))[1],"lost) after incorporating nuisance info")
    
    ##CHECKING ON MISSING DATA##
    miss_nuis_dat<-t_conn_data_trim_test_imput[!t_conn_data_trim_test_imput$subjectkey%in%t_conn_data_trim_test_nuisance$subjectkey,]
    message(dim(miss_nuis_dat)[1]," subs are missing nuisance dat: ",miss_nuis_dat$subjectkey)
    
    
    ##DO THE SAME FOR SYMPTOMS##  
    symptoms_test_nuisance<-merge(symptoms_test,nuisance_dat_test)#By merging like this I'm just keeping the nuisance vars that are in train set
    
    message("There are",dim(data.frame(unique(symptoms_test_nuisance$subjectkey)))[1]," subs (",
            dim(data.frame(unique(symptoms_test$subjectkey)))[1]-dim(data.frame(unique(symptoms_test_nuisance$subjectkey)))[1],"lost) after incorporating nuisance info")
    
    
    message("Regressing out site - step 2: regress out site")
    message("Brain")
    # nuisance_covs[,c(tmpt)]
    #NOTE: BEFORE I WAS REGRESSING OUT SITE ON THE LONG TABLE. BUT I COULD JUST REGRESS OUT SITE FOR EACH OF THESE COLUMNS IN THE WIDE DIRECTORY.
    t_conn_data_trim_test_regout<-matrix(data = NA,nrow = dim(t_conn_data_trim_test_nuisance)[1],ncol = dim(t_conn_data_trim_test_nuisance[,-which(colnames(t_conn_data_trim_test_nuisance) %in% 
                                                                                                                                                     c("subjectkey","interview_age","t_since_inj",nuisance_vars))])[2])
    t_conn_data_trim_test_regout<-as.data.frame(t_conn_data_trim_test_regout)
    # dim(t_conn_data_trim_test_regout_tmp)
    colnames(t_conn_data_trim_test_regout)<-colnames(t_conn_data_trim_test_nuisance[,-which(colnames(t_conn_data_trim_test_nuisance) %in% c("subjectkey","interview_age","t_since_inj",nuisance_vars))])
    for (i in 1:dim(t_conn_data_trim_test_regout)[2]){
      #Loop through the data frame but only the connections
      message("Regressing out ", nuisance_vars," for ",colnames(t_conn_data_trim_test_regout)[i]," [",i,"/",dim(t_conn_data_trim_test_regout)[2],"]")
      in_dat<-t_conn_data_trim_test_nuisance[,c(colnames(t_conn_data_trim_test_regout)[i],nuisance_vars)]
      colnames(in_dat)<-c("mean",nuisance_vars)
      # Iteratively change your input data for the model. Since the col names in the output table have been preset, you can use those names to select the input conn. I have to do this because you can't input variable names
      #in lm apparently. But I don't want to write out the nuisance vars instead we want to change them. So by doing it this way, we can use variable nuisance vars. The '.' in the formula means: use all other columns.
      lmod<-lm(mean ~ .,data = in_dat)#
      t_conn_data_trim_test_regout[,i]<-residuals(lmod)
    }
    
    t_conn_data_trim_test_regout$subjectkey<-t_conn_data_trim_test_nuisance$subjectkey#Add back subjectkey
    t_conn_data_trim_test_regout<-as.data.frame(t_conn_data_trim_test_regout)
    
    # write.table(t_conn_data_trim_test_regout,paste(out_dir,"conn_features_thresh_",thresh,"_regd_test.txt",sep=""),quote = F)
    
    
    message("Symptoms")
    # nuisance_vars_symps = c("sex","ehi_y_ss_scoreb","pubertal_stage")
    symptoms_test_regout<-matrix(data = NA,nrow = dim(symptoms_test_nuisance)[1],ncol = dim(symptoms_test_nuisance[,-which(colnames(symptoms_test_nuisance) %in% 
                                                                                                                             c("subjectkey","interview_age","t_since_inj","mri_info_deviceserialnumber",nuisance_vars_symps))])[2])
    symptoms_test_regout<-as.data.frame(symptoms_test_regout)
    dim(symptoms_test_regout)
    colnames(symptoms_test_regout)<-colnames(symptoms_test_nuisance[,-which(colnames(symptoms_test_nuisance) %in% c("subjectkey","interview_age","t_since_inj","mri_info_deviceserialnumber",nuisance_vars_symps))])
    for (i in 1:dim(symptoms_test_regout)[2]){
      #Loop through the data frame but only the connections
      message("Regressing out ", nuisance_vars_symps," for ",colnames(symptoms_test_regout)[i]," [",i,"/",dim(symptoms_test_regout)[2],"]")
      in_dat<-symptoms_test_nuisance[,c(colnames(symptoms_test_regout)[i],nuisance_vars_symps)]
      colnames(in_dat)<-c("mean",nuisance_vars_symps)
      # Iteratively change your input data for the model. Since the col names in the output table have been preset, you can use those names to select the input conn. I have to do this because you can't input variable names
      #in lm apparently. But I don't want to write out the nuisance vars instead we want to change them. So by doing it this way, we can use variable nuisance vars. The '.' in the formula means: use all other columns.
      lmod<-lm(mean ~ .,data = in_dat)
      symptoms_test_regout[,i]<-residuals(lmod)
    }
    
    symptoms_test_regout$subjectkey<-symptoms_test_nuisance$subjectkey#Add back subjectkey
    symptoms_test_regout<-as.data.frame(symptoms_test_regout)
    
    write.table(symptoms_test_regout,paste(out_dir,"symp_features_thresh_",thresh,"_regd_test.txt",sep=""),quote = F)
    
    
  } else if (resume_regout==TRUE){
    message("Resuming from regressing out site step")
    t_conn_data_trim_test_regout<-read.table(regout_file, sep=' ', header=T)####IMPORTANT NOTE TO SELF: IM LOADING THE SAME FILE FOR TRAIN AND TEST. I DOUBT IT'S TOO IMPORTANT BECAUSE THE SCRIPT RUNS MUCH FASTER AFTER THRESHOLDING SO I PROBABLY WON'T BE SAVING THE OUTPUT OF REGOUT
    t_conn_data_trim_test_regout<-data.frame(t_conn_data_trim_test_regout)
    
    #NOTE: these lines below are wrong since I'm importing the same file as above. It's just a placeholder, most likely I will be removing this option of resuming here.
    symptoms_test_regout<-read.table(regout_file, sep=' ', header=T)
    symptoms_test_regout<-data.frame(t_conn_data_trim_test_regout)
  }
  
  ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#@@@@@@@@@@@@@@@@@@@@@@@
  
  
  
  
  
  ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ADDITIONAL DATA TRANSFORMATIONS: TESTING SET#@@@@@@@@@@@@@@@@@@@@@@@
  #Check if there are any NANs
  #And for the moment, set the entire column to 0, as another way to trim, note that there might be mistakes so this trim needs to be double checked
  message("Checking that NaNs have been removed")
  tmp_wide_off<-t_conn_data_trim_test_regout
  for (i in 1:dim(tmp_wide_off)[2]){
    if (sum(is.na(tmp_wide_off[,i]))>0){
      sub<-tmp_wide_off[is.na(tmp_wide_off[,i])==T,c("subjectkey")]
      col<-colnames(t_conn_data_trim_test_regout)[i]
      message("Found missing value in col ",col," for sub ",sub," remember to threshold again!")
      # tmp_wide[is.na(tmp_wide[,i])==T,i]<-0
      # tmp_wide[,i]<-NA
    }
  }
  if (sum(is.na(tmp_wide_off))==0){
    message("You got rid of all NaNs, you may proceed ...")
  }
  
  
  
  message("Quick sanity check: Are symptoms and connectivity features the same dims?")
  message("Conn feats: ",dim(t_conn_data_trim_test_regout)[1],"Symps:",dim(symptoms_test_regout)[1])
  ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#@@@@@@@@@@@@@@@@@@@@@@@
  
  
  ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@SETTING UP DF OUTPUT: TRAINING SET#@@@@@@@@@@@@@@@@@@@@@@@
  message("Setting up DF output - Testing set")
  ##I"M CLOSE. THE THING IS, IF I DON'T SPLIT THE METRICS HERE I'll have to do it outside the function and that will add code
  conn_feat_names<-as.data.frame(colnames(t_conn_data_trim_test_regout))
  colnames(conn_feat_names)<-"names"
  
  split_conn_feats_test<-vector(mode = "list", length = dim(metrics_df)[1])
  for(i in 1:dim(metrics_df)[1]){
    mtr_find<-apply(conn_feat_names,1,function(x){regexpr(metrics_df[i,],x)})  
    split_conn_feats_test[[i]]<-t_conn_data_trim_test_regout[,mtr_find>0]
  }
  
  names(split_conn_feats_test)<-metrics
  
  mydata<-list(yy_train=symptoms_train_regout,xx_train=split_conn_feats_train,nuisance_train=nuisance_dat_train,yy_test=symptoms_test_regout,xx_test=split_conn_feats_test,nuisance_test=nuisance_dat_test)
  return(mydata)
  ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#@@@@@@@@@@@@@@@@@@@@@@@
}
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
pre_proc_noScanRO <- function(conn_data,symp_data,nuisance_dat,metrics,nuisance_vars,nuisance_vars_symps,thresh,tr_ratio,out_dir,resume_thresh,thresh_file,resume_regout,regout_file,resume_scale,scale_file){
  # function to do the pre-processing without the step for regressing out site, done for validation analyses
  # conn_data:the connectivity data after QC removal
  # symp_data:the symptoms data
  # site_dat: site data
  # thresh: threshold 
  # tr_ratio: training/testing group ratio
  
  metrics_df<-as.data.frame(metrics)#Turning metrics into a data frame (will be used further down)
  colnames(metrics_df)<-"metric"
  
  conn_data<-conn_data[conn_data$metric %in% metrics,]#Restricting the conn_data to the metrics the user inputed
  conn_data$metric<-factor(conn_data$metric, levels=metrics)
  
  message("We are starting with ",dim(data.frame(unique(conn_data$sub)))[1]," subs, processing metrics ",metrics,sep=" ")
  
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@THRESHOLDING@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  if (resume_thresh==FALSE){
    # t_str_conn_feats<-thresh_conn(str_conn_feats_QCpass,str_conn_feats_QCpass$conn,str_conn_feats_QCpass$sub,str_conn_feats_QCpass$metric,1)#For troubleshoot
    message("Thresholding conns")
    t_conn_data<-thresh_conn(conn_data,conn_data$conn,conn_data$sub,metrics,thresh)
    
    message("Flagging thresholded conns")
    t_conn_data$keep<-0
    t_conn_data[is.na(t_conn_data$mean)==0,c("keep")]<-1#HERE I FLAG CONNS TO REMOVE. I'M FLAGGING THEM BY LOOKING AT THE ONES THAT HAVE NANS. REMEMBER, THOSE THAT ARE MISSING BUT FROM CONNS THAT PASSED THRESHOLD, THEY ARE 0, NOT NANS.
    
    message(paste("After thresholding with t=",thresh,", ",sum(t_conn_data$keep),"/",dim(t_conn_data)[1]," features were kept"))
    
    write.table(t_conn_data,paste(out_dir,"/conn_features_thresh_",thresh,sep=""),quote = F) ##REMEMBER: CHANGE THE OUTPUT NAME HERE ONCE YOU ARE SURE THE SCRIPT WORKS##
    
    t_conn_data_trim<-t_conn_data[t_conn_data$keep==1,]#Trimming first, since we don't want to use missing data to regress out site. HERE IS WHERE I'M REMOVING ROWS BELONGING TO CONNS THAT DID NOT PASS THRESHOLD.
    message(paste("There are ",dim(t_conn_data)[1]-dim(t_conn_data_trim)[1]," features that were removed after thresholding"))
    
  } else if (resume_thresh==TRUE){
    message("Resuming from threshold step")
    t_conn_data<-read.table(thresh_file, sep=' ', header=T)
    t_conn_data<-data.frame(t_conn_data)
    
    t_conn_data_trim<-t_conn_data[t_conn_data$keep==1,]#Trimming first, since we don't want to use missing data to regress out site. 
    message(paste("There are ",dim(t_conn_data)[1]-dim(t_conn_data_trim)[1]," features that were removed after thresholding"))
    
    # 
    # message("Resuming from threshold step")
    # t_conn_data<-read.table('/Users/Guido/Desktop/Projects/PCHStudy/Data/conn_features_thresh_0.9_newesttest.txt', sep=' ', header=T)
    # t_conn_data<-data.frame(t_conn_data)
    # 
    # t_conn_data_trim<-t_conn_data[t_conn_data$keep==1,]#Trimming first, since we don't want to use missing data to regress out site. 
    # message(paste("There are ",dim(t_conn_data)[1]-dim(t_conn_data_trim)[1]," features that were removed after thresholding"))
    # 
    # 
  }
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@SPLITTING DATASETS@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  message("Splitting dataset into training and testing sets")
  # colnames(t_conn_data_trim_site_tmp)<-c("subjectkey","mean","conn","metric","keep","mri_info_deviceserialnumber","residuals")
  subs<-data.frame(unique(t_conn_data_trim$sub))
  # dim(subs)
  colnames(subs)<-"subjectkey"
  smp_size=floor(tr_ratio*nrow(subs))
  set.seed(123)
  train_idx<-sample(seq_len(nrow(subs)),size=smp_size)
  subs$grp<-2
  subs[train_idx,c("grp")]<-1
  
  t_conn_data_trim<-merge(t_conn_data_trim,subs,by.x = "sub",by.y="subjectkey")
  
  ##HERE IS WHERE THE DATA IS ACTUALLY SPLIT IN TWO SETS##
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@SPLIT CONN DATA: TRAINING SET@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  t_conn_data_trim_train<-t_conn_data_trim[t_conn_data_trim$grp==1,]
  message(paste("There are ",dim(as.data.frame(unique(t_conn_data_trim_train$sub)))[1],"subs in the train set, ",
                dim(as.data.frame(unique(t_conn_data_trim$sub)))[1]-dim(as.data.frame(unique(t_conn_data_trim_train$sub)))[1],"in the test set"))
  
  #Changing data structure#
  message("Changing data structure - long to wide")
  t_conn_data_trim_train_halfwide<-reshape(t_conn_data_trim_train[,c("sub","mean","conn","metric")],idvar = c("sub","metric"),timevar="conn",direction="wide")
  
  t_conn_data_trim_train_wide<-reshape(t_conn_data_trim_train_halfwide,idvar = "sub",timevar="metric",direction="wide")
  message(paste("Data in wide has the following dimensions:",dim(t_conn_data_trim_train_wide)[1],dim(t_conn_data_trim_train_wide)[2],sep=" "))
  
  tmp_full_miss<-t_conn_data_trim_train_wide[,colSums(is.na(t_conn_data_trim_train_wide))==
                                               dim(t_conn_data_trim_train_wide)[1]]#Here I'm calculating how many NAs are in a column. If the number of NAs in a column == the dim[2] of the table,
  # it means the entire column is filled with NAs. RECALL THAT I REMOVED SUBTHRESHOLD CONNS ALREADY. I'M HAVING TO DO THIS AGAIN BECAUSE SOME METRICS HAD CONNS THAT SURVIVED THRESHOLD WHEREAS SOME DIDN'T.
  # WHEN I RESHAPE TO WIDE, IT'S EXPECTING EVERYTHING TO BE REPEATED (E.G.: IF I HAVE METRIC=FA, CONN=10_11, IT'S EXPECTING ME TO HAVE METRIC=NUFO, CONN=10_11. IF IT'S MISSING, IT WILL CREATE IT IN THE WIDE TABLE AND FILL IT WITH NANS.
  # I NEED TO FLAG THESE AND REMOVE THEM.
  
  message("Due to lack of some conns in specific metrics but not others, ",dim(tmp_full_miss)[2])
  message("Conns that were introduced by error and now removed: ",colnames(tmp_full_miss))
  t_conn_data_trim_train_wide<-t_conn_data_trim_train_wide[,!colSums(is.na(t_conn_data_trim_train_wide))==dim(t_conn_data_trim_train_wide)[1]] #SO HERE I'M KEEPING THE COLUMNS THAT WEREN'T INTRODUCED BY MISTAKE DURING THE RESHAPING.
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  
  
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@IMPUTING CONN DATA: TRAINING SET@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  if (thresh<1){
    message("Imputing missing conn data - training set")
    
    t_conn_data_trim_train_imput<-t_conn_data_trim_train_wide[,-1]#Needs to be the full wide table. First remove the subjects.
    
    for (i in 1:dim(t_conn_data_trim_train_imput)[2]){
      idx_miss<-rownames(t_conn_data_trim_train_imput[is.na(t_conn_data_trim_train_imput[,i])==1,])
      idx_miss<-data.frame(idx_miss)
      idx_comp<-rownames(t_conn_data_trim_train_imput[is.na(t_conn_data_trim_train_imput[,i])==0,])
      idx_comp<-data.frame(idx_comp)
      if (dim(idx_miss)[1]!=0 && dim(idx_miss)[1]<dim(idx_comp)[1]){
        imp_size=dim(as.data.frame(idx_miss))[1]
        set.seed(i)
        idx_rep<-sample(as.vector(idx_comp$idx_comp),size=imp_size)
        t_conn_data_trim_train_imput[rownames(t_conn_data_trim_train_imput)==as.numeric(as.character(idx_miss$idx_miss)),i]<-t_conn_data_trim_train_imput[idx_rep,i]
      }
    }
  } else if (thresh==1){
    t_conn_data_trim_train_imput<-t_conn_data_trim_train_wide
  }
  
  
  t_conn_data_trim_train_imput$subjectkey<-t_conn_data_trim_train_wide$sub#Add back subjectkey to perform the merge of the nuisance vars
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#@@@@@@@@@@@@@
  
  
  ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@SPLIT SYMPTOM DATA: TRAINING SET#@@@@@@@@@@@@@@@@@@@@@@@
  message("Splitting symptoms into train/test sets")
  symptoms_train<-symp_data[as.character(symp_data$subjectkey) %in% t_conn_data_trim_train_imput$subjectkey,] # filter symptom data by looking at subs that are in training set
  
  ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@IMPUTING MISSING SYMPTOM DATA: TRAINING SET#@@@@@@@@@@@@@@@@@@@@@@@
  message("Imputing missing symp data - training set")
  for (i in 1:dim(symptoms_train)[2]){
    idx_miss<-rownames(symptoms_train[is.na(symptoms_train[,i])==1,])
    if (length(idx_miss)!=0){
      imp_size=dim(as.data.frame(idx_miss))[1]
      # set.seed(123) # Should we set a seed?
      idx_comp<-rownames(symptoms_train[is.na(symptoms_train[,i])==0,])
      set.seed(i)
      idx_rep<-sample(idx_comp,size=imp_size)
      symptoms_train[idx_miss,i]<-symptoms_train[idx_rep,i]
    }
  }
  
  ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#@@@@@@@@@@@@@@@@@@@@@@@
  
  ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@IMPUTING MISSING NUISANCE DATA: TRAINING SET#@@@@@@@@@@@@@@@@@@@@@@@
  nuisance_dat_train<-nuisance_dat[as.character(nuisance_dat$subjectkey) %in% t_conn_data_trim_train_imput$subjectkey,] # filter symptom data by looking at subs that are in training set
  message("Imputing missing nuisance covariates - training set")
  for (i in 1:dim(nuisance_dat_train)[2]){
    idx_miss<-rownames(nuisance_dat_train[is.na(nuisance_dat_train[,i])==1,])
    if (length(idx_miss)!=0){
      imp_size=dim(as.data.frame(idx_miss))[1]
      # set.seed(123) # Should we set a seed?
      idx_comp<-rownames(nuisance_dat_train[is.na(nuisance_dat_train[,i])==0,])
      set.seed(i)
      idx_rep<-sample(idx_comp,size=imp_size)
      nuisance_dat_train[idx_miss,i]<-nuisance_dat_train[idx_rep,i]
    }
  }
  
  ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
  ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@REGRESSING OUT NUISANCE VARS: TRAINING SET#@@@@@@@@@@@@@@@@@@@@@@@
  if (resume_regout==FALSE){
    message("Regressing out nuisance vars - step 1: merging nuisance vars")
    
    t_conn_data_trim_train_nuisance<-merge(t_conn_data_trim_train_imput,nuisance_dat_train)#By merging like this I'm just keeping the nuisance vars that are in train set
    
    message("There are",dim(data.frame(unique(t_conn_data_trim_train_nuisance$subjectkey)))[1]," subs (",
            dim(data.frame(unique(t_conn_data_trim_train_imput$subjectkey)))[1]-dim(data.frame(unique(t_conn_data_trim_train_nuisance$subjectkey)))[1],"lost) after incorporating nuisance info")
    
    ##DO THE SAME FOR SYMPTOMS##
    symptoms_train_nuisance<-merge(symptoms_train,nuisance_dat_train)#By merging like this I'm just keeping the nuisance vars that are in train set
    
    message("There are",dim(data.frame(unique(symptoms_train_nuisance$subjectkey)))[1]," subs (",
            dim(data.frame(unique(symptoms_train$subjectkey)))[1]-dim(data.frame(unique(symptoms_train_nuisance$subjectkey)))[1],"lost) after incorporating nuisance info")
    
    
    
    ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@REGRESSING OUT NUISANCE VARS: TRAINING SET#@@@@@@@@@@@@@@@@@@@@@@@
    message("Regressing out site - step 2: regress out site")
    message("Brain")
    t_conn_data_trim_train_regout<-matrix(data = NA,nrow = dim(t_conn_data_trim_train_nuisance)[1],ncol = dim(t_conn_data_trim_train_nuisance[,-which(colnames(t_conn_data_trim_train_nuisance) %in%
                                                                                                                                                        c("subjectkey","interview_age","t_since_inj",nuisance_vars))])[2])
    t_conn_data_trim_train_regout<-as.data.frame(t_conn_data_trim_train_regout)
    dim(t_conn_data_trim_train_regout)
    colnames(t_conn_data_trim_train_regout)<-colnames(t_conn_data_trim_train_nuisance[,-which(colnames(t_conn_data_trim_train_nuisance) %in% c("subjectkey","interview_age","t_since_inj",nuisance_vars))])
    for (i in 1:dim(t_conn_data_trim_train_regout)[2]){
      #Loop through the data frame but only the connections
      message("Regressing out ", nuisance_vars," for ",colnames(t_conn_data_trim_train_regout)[i]," [",i,"/",dim(t_conn_data_trim_train_regout)[2],"]")
      in_dat<-t_conn_data_trim_train_nuisance[,c(colnames(t_conn_data_trim_train_regout)[i],nuisance_vars)]
      colnames(in_dat)<-c("mean",nuisance_vars)
      # Iteratively change your input data for the model. Since the col names in the output table have been preset, you can use those names to select the input conn. I have to do this because you can't input variable names
      #in lm apparently. But I don't want to write out the nuisance vars instead we want to change them. So by doing it this way, we can use variable nuisance vars. The '.' in the formula means: use all other columns.
      lmod<-lm(mean ~ .,data = in_dat)#
      t_conn_data_trim_train_regout[,i]<-residuals(lmod)
    }
    
    t_conn_data_trim_train_regout$subjectkey<-t_conn_data_trim_train_nuisance$subjectkey#Add back subjectkey
    t_conn_data_trim_train_regout<-as.data.frame(t_conn_data_trim_train_regout)
    
    # write.table(t_conn_data_trim_train_regout,paste(out_dir,"conn_features_thresh_",thresh,"_regd_train.txt",sep=""),quote = F)
    
    message("Symptoms")
    symptoms_train_regout<-matrix(data = NA,nrow = dim(symptoms_train_nuisance)[1],ncol = dim(symptoms_train_nuisance[,-which(colnames(symptoms_train_nuisance) %in%
                                                                                                                                c("subjectkey","interview_age","t_since_inj","mri_info_deviceserialnumber",nuisance_vars_symps))])[2])
    symptoms_train_regout<-as.data.frame(symptoms_train_regout)
    dim(symptoms_train_regout)
    #NOTE: notice that I'm regressing out different things. For symptoms, I'm not regressing out mri site. So that's why I have nuisance_vars_symps
    colnames(symptoms_train_regout)<-colnames(symptoms_train_nuisance[,-which(colnames(symptoms_train_nuisance) %in% c("subjectkey","interview_age","t_since_inj","mri_info_deviceserialnumber",nuisance_vars_symps))])
    for (i in 1:dim(symptoms_train_regout)[2]){
      #Loop through the data frame but only the connections
      message("Regressing out ", nuisance_vars_symps," for ",colnames(symptoms_train_regout)[i]," [",i,"/",dim(symptoms_train_regout)[2],"]")
      in_dat<-symptoms_train_nuisance[,c(colnames(symptoms_train_regout)[i],nuisance_vars_symps)]
      colnames(in_dat)<-c("mean",nuisance_vars_symps)
      # Iteratively change your input data for the model. Since the col names in the output table have been preset, you can use those names to select the input conn. I have to do this because you can't input variable names
      #in lm apparently. But I don't want to write out the nuisance vars instead we want to change them. So by doing it this way, we can use variable nuisance vars. The '.' in the formula means: use all other columns.
      lmod<-lm(mean ~ .,data = in_dat)#
      symptoms_train_regout[,i]<-residuals(lmod)
    }
    
    symptoms_train_regout$subjectkey<-symptoms_train_nuisance$subjectkey#Add back subjectkey
    symptoms_train_regout<-as.data.frame(symptoms_train_regout)
    
    write.table(symptoms_train_regout,paste(out_dir,"symp_features_thresh_",thresh,"_regd_train.txt",sep=""),quote = F)
    
    
    
    
  } else if (resume_regout==TRUE){
    message("Resuming from regressing out site step")
    
    t_conn_data_trim_train_regout<-read.table(regout_file, sep=' ', header=T)
    t_conn_data_trim_train_regout<-data.frame(t_conn_data_trim_train_regout)
    
    #NOTE: these lines below are wrong since I'm importing the same file as above. It's just a placeholder, most likely I will be removing this option of resuming here.
    symptoms_train_regout<-read.table(regout_file, sep=' ', header=T)
    symptoms_train_regout<-data.frame(t_conn_data_trim_train_regout)
    
  }
  ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
  ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ADDITIONAL DATA TRANSFORMATIONS: TRAINING SET#@@@@@@@@@@@@@@@@@@@@@@@
  #Check if there are any NANs
  #And for the moment, set the entire column to 0, as another way to trim, note that there might be mistakes so this trim needs to be double checked
  message("Checking that NaNs have been removed")
  tmp_wide_off<-t_conn_data_trim_train_regout
  for (i in 1:dim(tmp_wide_off)[2]){
    if (sum(is.na(tmp_wide_off[,i]))>0){
      sub<-tmp_wide_off[is.na(tmp_wide_off[,i])==T,c("subjectkey")]
      col<-colnames(t_conn_data_trim_regout)[i]
      message("Found missing value in col ",col," for sub ",sub," remember to threshold again!")
      # tmp_wide[is.na(tmp_wide[,i])==T,i]<-0
      # tmp_wide[,i]<-NA
    }
  }
  if (sum(is.na(tmp_wide_off))==0){
    message("You got rid of all NaNs, you may proceed ...")
  }
  
  message("Quick sanity check: Are symptoms and connectivity features the same dims?")
  message("Conn feats: ",dim(t_conn_data_trim_train_regout)[1],"Symps:",dim(symptoms_train_regout)[1])
  ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#@@@@@@@@@@@@@@@@@@@@@@@
  
  
  ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@SETTING UP DF OUTPUT: TRAINING SET#@@@@@@@@@@@@@@@@@@@@@@@
  message("Setting up DF output - Training set")
  ##I"M CLOSE. THE THING IS, IF I DON'T SPLIT THE METRICS HERE I'll have to do it outside the function and that will add code
  conn_feat_names<-as.data.frame(colnames(t_conn_data_trim_train_regout))
  colnames(conn_feat_names)<-"names"
  split_conn_feats_train<-vector(mode = "list", length = dim(metrics_df)[1])
  for(i in 1:dim(metrics_df)[1]){
    mtr_find<-apply(conn_feat_names,1,function(x){regexpr(metrics_df[i,],x)})
    split_conn_feats_train[[i]]<-t_conn_data_trim_train_regout[,mtr_find>0]
  }
  
  names(split_conn_feats_train)<-metrics
  ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#@@@@@@@@@@@@@@@@@@@@@@@
  
  
  
  ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@SPLIT CONN DATA: TESTING SET#@@@@@@@@@@@@@@@@@@@@@@@
  t_conn_data_trim_test<-t_conn_data_trim[t_conn_data_trim$grp==2,]
  
  #Changing data structure#
  message("Changing data structure - long to wide") 
  #NOTE: had to proceed with this two-step reshape because I had multiple repeating units
  t_conn_data_trim_test_halfwide<-reshape(t_conn_data_trim_test[,c("sub","mean","conn","metric")],idvar = c("sub","metric"),timevar="conn",direction="wide")
  
  t_conn_data_trim_test_wide<-reshape(t_conn_data_trim_test_halfwide,idvar = "sub",timevar="metric",direction="wide")
  
  message(paste("Data in wide has the following dimensions:",dim(t_conn_data_trim_test_wide)[1],dim(t_conn_data_trim_test_wide)[2],sep=" "))
  
  tmp_full_miss<-t_conn_data_trim_test_wide[,colSums(is.na(t_conn_data_trim_test_wide))==
                                              dim(t_conn_data_trim_test_wide)[1]]#Here I'm calculating how many NAs are in a column. If the number of NAs in a column == the dim[2] of the table, 
  # it means the entire column is filled with NAs. RECALL THAT I REMOVED SUBTHRESHOLD CONNS ALREADY. I'M HAVING TO DO THIS AGAIN BECAUSE SOME METRICS HAD CONNS THAT SURVIVED THRESHOLD WHEREAS SOME DIDN'T. 
  # WHEN I RESHAPE TO WIDE, IT'S EXPECTING EVERYTHING TO BE REPEATED (E.G.: IF I HAVE METRIC=FA, CONN=10_11, IT'S EXPECTING ME TO HAVE METRIC=NUFO, CONN=10_11. IF IT'S MISSING, IT WILL CREATE IT IN THE WIDE TABLE AND FILL IT WITH NANS.
  # I NEED TO FLAG THESE AND REMOVE THEM.
  
  message("Due to lack of some conns in specific metrics but not others, ",dim(tmp_full_miss)[2])
  message("Conns that were introduced by error and now removed: ",colnames(tmp_full_miss))
  t_conn_data_trim_test_wide<-t_conn_data_trim_test_wide[,!colSums(is.na(t_conn_data_trim_test_wide))==dim(t_conn_data_trim_test_wide)[1]] #SO HERE I'M KEEPING THE COLUMNS THAT WEREN'T INTRODUCED BY MISTAKE DURING THE RESHAPING.
  
  ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#@@@@@@@@@@@@@@@@@@@@@@@
  
  ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@IMPUTING MISSING CONN DATA: TESTING SET#@@@@@@@@@@@@@@@@@@@@@@@  
  if (thresh<1){
    message("Imputing missing data - testing set")
    
    t_conn_data_trim_test_imput<-t_conn_data_trim_test_wide[,-1]#Needs to be the full wide table. First remove the subjects.
    
    for (i in 1:dim(t_conn_data_trim_test_imput)[2]){
      idx_miss<-rownames(t_conn_data_trim_test_imput[is.na(t_conn_data_trim_test_imput[,i])==1,])
      idx_miss<-data.frame(idx_miss)
      idx_comp<-rownames(t_conn_data_trim_test_imput[is.na(t_conn_data_trim_test_imput[,i])==0,])
      idx_comp<-data.frame(idx_comp)
      if (dim(idx_miss)[1]!=0 && dim(idx_miss)[1]<dim(idx_comp)[1]){
        imp_size=dim(as.data.frame(idx_miss))[1]
        
        set.seed(i)
        idx_rep<-sample(as.vector(idx_comp$idx_comp),size=imp_size)
        t_conn_data_trim_test_imput[rownames(t_conn_data_trim_test_imput)==as.numeric(as.character(idx_miss$idx_miss)),i]<-t_conn_data_trim_test_imput[idx_rep,i]
      }
    }
    
  } else if (thresh==1){
    t_conn_data_trim_test_imput<-t_conn_data_trim_test_wide
  }
  
  
  t_conn_data_trim_test_imput$subjectkey<-t_conn_data_trim_test_wide$sub#Add back subjectkey to perform the merge of the nuisance vars
  ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#@@@@@@@@@@@@@@@@@@@@@@@
  
  
  ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@SPLITING SYMPTOM DATA: TESTING SET#@@@@@@@@@@@@@@@@@@@@@@@  
  message("Splitting symptoms into train/test sets")  
  symptoms_test<-symp_data[as.character(symp_data$subjectkey) %in% t_conn_data_trim_test_imput$subjectkey,] # filter symptom data by looking at subs that are in training set
  
  # symptoms_test<-symptoms_mTBI[as.character(symptoms_mTBI$subjectkey) %in% t_conn_data_trim_test_imput$subjectkey,] # filter symptom data by looking at subs that are in training set
  ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#@@@@@@@@@@@@@@@@@@@@@@@
  
  ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@IMPUTING MISSING SYMPTOM DATA: TESTING SET#@@@@@@@@@@@@@@@@@@@@@@@  
  for (i in 1:dim(symptoms_test)[2]){
    idx_miss<-rownames(symptoms_test[is.na(symptoms_test[,i])==1,])
    if (length(idx_miss)!=0){
      imp_size=dim(as.data.frame(idx_miss))[1]
      # set.seed(123) # Should we set a seed?
      idx_comp<-rownames(symptoms_test[is.na(symptoms_test[,i])==0,])
      set.seed(i)
      idx_rep<-sample(idx_comp,size=imp_size)
      symptoms_test[idx_miss,i]<-symptoms_test[idx_rep,i]
    }
  }
  ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#@@@@@@@@@@@@@@@@@@@@@@@
  
  ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@IMPUTING MISSING NUISANCE DATA: TESTING SET#@@@@@@@@@@@@@@@@@@@@@@@  
  nuisance_dat_test<-nuisance_dat[as.character(nuisance_dat$subjectkey) %in% t_conn_data_trim_test_imput$subjectkey,] # filter symptom data by looking at subs that are in training set
  
  # nuisance_dat_test<-nuisance_covs[as.character(nuisance_covs$subjectkey) %in% t_conn_data_trim_test_imput$subjectkey,] # filter symptom data by looking at subs that are in training set
  message("Imputing missing nuisance covariates")
  for (i in 1:dim(nuisance_dat_test)[2]){
    idx_miss<-rownames(nuisance_dat_test[is.na(nuisance_dat_test[,i])==1,])
    if (length(idx_miss)!=0){
      imp_size=dim(as.data.frame(idx_miss))[1]
      # set.seed(123) # Should we set a seed?
      idx_comp<-rownames(nuisance_dat_test[is.na(nuisance_dat_test[,i])==0,])
      set.seed(i)
      idx_rep<-sample(idx_comp,size=imp_size)
      nuisance_dat_test[idx_miss,i]<-nuisance_dat_test[idx_rep,i]
    }
  }
  
  ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#@@@@@@@@@@@@@@@@@@@@@@@
  
  
  ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@REGRESSING OUT NUISANCE DATA: TESTING SET#@@@@@@@@@@@@@@@@@@@@@@@  
  
  if (resume_regout==FALSE){
    message("Regressing out nuisance vars - step 1: merging nuisance vars")
    t_conn_data_trim_test_nuisance<-merge(t_conn_data_trim_test_imput,nuisance_dat_test)
    
    message("There are",dim(data.frame(unique(t_conn_data_trim_test_nuisance$subjectkey)))[1]," subs (",
            dim(data.frame(unique(t_conn_data_trim_test_imput$subjectkey)))[1]-dim(data.frame(unique(t_conn_data_trim_test_nuisance$subjectkey)))[1],"lost) after incorporating nuisance info")
    
    ##CHECKING ON MISSING DATA##
    miss_nuis_dat<-t_conn_data_trim_test_imput[!t_conn_data_trim_test_imput$subjectkey%in%t_conn_data_trim_test_nuisance$subjectkey,]
    message(dim(miss_nuis_dat)[1]," subs are missing nuisance dat: ",miss_nuis_dat$subjectkey)
    
    
    ##DO THE SAME FOR SYMPTOMS##  
    symptoms_test_nuisance<-merge(symptoms_test,nuisance_dat_test)#By merging like this I'm just keeping the nuisance vars that are in train set
    
    message("There are",dim(data.frame(unique(symptoms_test_nuisance$subjectkey)))[1]," subs (",
            dim(data.frame(unique(symptoms_test$subjectkey)))[1]-dim(data.frame(unique(symptoms_test_nuisance$subjectkey)))[1],"lost) after incorporating nuisance info")
    
    
    message("Regressing out site - step 2: regress out site")
    message("Brain")
    # nuisance_covs[,c(tmpt)]
    #NOTE: BEFORE I WAS REGRESSING OUT SITE ON THE LONG TABLE. BUT I COULD JUST REGRESS OUT SITE FOR EACH OF THESE COLUMNS IN THE WIDE DIRECTORY.
    t_conn_data_trim_test_regout<-matrix(data = NA,nrow = dim(t_conn_data_trim_test_nuisance)[1],ncol = dim(t_conn_data_trim_test_nuisance[,-which(colnames(t_conn_data_trim_test_nuisance) %in% 
                                                                                                                                                     c("subjectkey","interview_age","t_since_inj",nuisance_vars))])[2])
    t_conn_data_trim_test_regout<-as.data.frame(t_conn_data_trim_test_regout)
    # dim(t_conn_data_trim_test_regout_tmp)
    colnames(t_conn_data_trim_test_regout)<-colnames(t_conn_data_trim_test_nuisance[,-which(colnames(t_conn_data_trim_test_nuisance) %in% c("subjectkey","interview_age","t_since_inj",nuisance_vars))])
    for (i in 1:dim(t_conn_data_trim_test_regout)[2]){
      #Loop through the data frame but only the connections
      message("Regressing out ", nuisance_vars," for ",colnames(t_conn_data_trim_test_regout)[i]," [",i,"/",dim(t_conn_data_trim_test_regout)[2],"]")
      in_dat<-t_conn_data_trim_test_nuisance[,c(colnames(t_conn_data_trim_test_regout)[i],nuisance_vars)]
      colnames(in_dat)<-c("mean",nuisance_vars)
      # Iteratively change your input data for the model. Since the col names in the output table have been preset, you can use those names to select the input conn. I have to do this because you can't input variable names
      #in lm apparently. But I don't want to write out the nuisance vars instead we want to change them. So by doing it this way, we can use variable nuisance vars. The '.' in the formula means: use all other columns.
      lmod<-lm(mean ~ .,data = in_dat)#
      t_conn_data_trim_test_regout[,i]<-residuals(lmod)
    }
    
    t_conn_data_trim_test_regout$subjectkey<-t_conn_data_trim_test_nuisance$subjectkey#Add back subjectkey
    t_conn_data_trim_test_regout<-as.data.frame(t_conn_data_trim_test_regout)
    
    # write.table(t_conn_data_trim_test_regout,paste(out_dir,"conn_features_thresh_",thresh,"_regd_test.txt",sep=""),quote = F)
    
    
    message("Symptoms")
    # nuisance_vars_symps = c("sex","ehi_y_ss_scoreb","pubertal_stage")
    symptoms_test_regout<-matrix(data = NA,nrow = dim(symptoms_test_nuisance)[1],ncol = dim(symptoms_test_nuisance[,-which(colnames(symptoms_test_nuisance) %in% 
                                                                                                                             c("subjectkey","interview_age","t_since_inj","mri_info_deviceserialnumber",nuisance_vars_symps))])[2])
    symptoms_test_regout<-as.data.frame(symptoms_test_regout)
    dim(symptoms_test_regout)
    colnames(symptoms_test_regout)<-colnames(symptoms_test_nuisance[,-which(colnames(symptoms_test_nuisance) %in% c("subjectkey","interview_age","t_since_inj","mri_info_deviceserialnumber",nuisance_vars_symps))])
    for (i in 1:dim(symptoms_test_regout)[2]){
      #Loop through the data frame but only the connections
      message("Regressing out ", nuisance_vars_symps," for ",colnames(symptoms_test_regout)[i]," [",i,"/",dim(symptoms_test_regout)[2],"]")
      in_dat<-symptoms_test_nuisance[,c(colnames(symptoms_test_regout)[i],nuisance_vars_symps)]
      colnames(in_dat)<-c("mean",nuisance_vars_symps)
      # Iteratively change your input data for the model. Since the col names in the output table have been preset, you can use those names to select the input conn. I have to do this because you can't input variable names
      #in lm apparently. But I don't want to write out the nuisance vars instead we want to change them. So by doing it this way, we can use variable nuisance vars. The '.' in the formula means: use all other columns.
      lmod<-lm(mean ~ .,data = in_dat)
      symptoms_test_regout[,i]<-residuals(lmod)
    }
    
    symptoms_test_regout$subjectkey<-symptoms_test_nuisance$subjectkey#Add back subjectkey
    symptoms_test_regout<-as.data.frame(symptoms_test_regout)
    
    write.table(symptoms_test_regout,paste(out_dir,"symp_features_thresh_",thresh,"_regd_test.txt",sep=""),quote = F)
    
    
  } else if (resume_regout==TRUE){
    message("Resuming from regressing out site step")
    t_conn_data_trim_test_regout<-read.table(regout_file, sep=' ', header=T)####IMPORTANT NOTE TO SELF: IM LOADING THE SAME FILE FOR TRAIN AND TEST. I DOUBT IT'S TOO IMPORTANT BECAUSE THE SCRIPT RUNS MUCH FASTER AFTER THRESHOLDING SO I PROBABLY WON'T BE SAVING THE OUTPUT OF REGOUT
    t_conn_data_trim_test_regout<-data.frame(t_conn_data_trim_test_regout)
    
    #NOTE: these lines below are wrong since I'm importing the same file as above. It's just a placeholder, most likely I will be removing this option of resuming here.
    symptoms_test_regout<-read.table(regout_file, sep=' ', header=T)
    symptoms_test_regout<-data.frame(t_conn_data_trim_test_regout)
  }
  
  ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#@@@@@@@@@@@@@@@@@@@@@@@
  
  
  
  
  
  ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ADDITIONAL DATA TRANSFORMATIONS: TESTING SET#@@@@@@@@@@@@@@@@@@@@@@@
  #Check if there are any NANs
  #And for the moment, set the entire column to 0, as another way to trim, note that there might be mistakes so this trim needs to be double checked
  message("Checking that NaNs have been removed")
  tmp_wide_off<-t_conn_data_trim_test_regout
  for (i in 1:dim(tmp_wide_off)[2]){
    if (sum(is.na(tmp_wide_off[,i]))>0){
      sub<-tmp_wide_off[is.na(tmp_wide_off[,i])==T,c("subjectkey")]
      col<-colnames(t_conn_data_trim_test_regout)[i]
      message("Found missing value in col ",col," for sub ",sub," remember to threshold again!")
      # tmp_wide[is.na(tmp_wide[,i])==T,i]<-0
      # tmp_wide[,i]<-NA
    }
  }
  if (sum(is.na(tmp_wide_off))==0){
    message("You got rid of all NaNs, you may proceed ...")
  }
  
  
  
  message("Quick sanity check: Are symptoms and connectivity features the same dims?")
  message("Conn feats: ",dim(t_conn_data_trim_test_regout)[1],"Symps:",dim(symptoms_test_regout)[1])
  ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#@@@@@@@@@@@@@@@@@@@@@@@
  
  
  ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@SETTING UP DF OUTPUT: TRAINING SET#@@@@@@@@@@@@@@@@@@@@@@@
  message("Setting up DF output - Testing set")
  ##I"M CLOSE. THE THING IS, IF I DON'T SPLIT THE METRICS HERE I'll have to do it outside the function and that will add code
  conn_feat_names<-as.data.frame(colnames(t_conn_data_trim_test_regout))
  colnames(conn_feat_names)<-"names"
  
  split_conn_feats_test<-vector(mode = "list", length = dim(metrics_df)[1])
  for(i in 1:dim(metrics_df)[1]){
    mtr_find<-apply(conn_feat_names,1,function(x){regexpr(metrics_df[i,],x)})  
    split_conn_feats_test[[i]]<-t_conn_data_trim_test_regout[,mtr_find>0]
  }
  
  names(split_conn_feats_test)<-metrics
  
  mydata<-list(yy_train=symptoms_train_regout,xx_train=split_conn_feats_train,nuisance_train=nuisance_dat_train,yy_test=symptoms_test_regout,xx_test=split_conn_feats_test,nuisance_test=nuisance_dat_test)
  return(mydata)
  
}
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
row_permute<-function(data,seed_vect){
  perm_data<-matrix(data=NA,nrow=dim(data)[1],ncol=dim(data)[2])
  for (i in 1:dim(perm_data)[1]){
    set.seed(seed_vect[i])
    shuf_inds<-shuffle(dim(data)[2])
    perm_data[i,]<-data[i,shuf_inds]
  }
  return(perm_data)
}
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
uvar_feat_sel<-function(Yvar,Xvar,Nfeat,dir){
  univar_corrs<-cor(Yvar,Xvar)
  univar_corrs<-t(univar_corrs)
  # colnames(univar_corrs_tmp)<-colnames(X_tmp)
  univar_corrs<-melt(univar_corrs)
  colnames(univar_corrs)<-c("conn","symptom","corr")
  univar_corrs$corr<-abs(univar_corrs$corr)
  
  #SOLVING REPEATS HERE, SO THAT I CAN JUST SELECT THE TOPN AND NOT USE A BUFFER
  #The idea is to, for every conn, replace all its corrs with symptoms by the max corr with a symptom, that way we can get rid of duplicates and then select the topn
  #if I want to pick the bottom feats, I should select the reverse
  tmp_conn_nums<-as.data.frame(unique(univar_corrs$conn)) 
  colnames(tmp_conn_nums)<-"conn"
  for (i in 1:dim(tmp_conn_nums)[1]){
    tmp<-subset(univar_corrs,conn==tmp_conn_nums$conn[i])
    if (dir=="top"){
      univar_corrs[univar_corrs$conn==tmp_conn_nums$conn[i],c("corr")]<-max(tmp$corr)#Replace all corrs by the max corr, it'll be easier to get rid of duplications that way  
    } else if (dir=="bottom") {
      univar_corrs[univar_corrs$conn==tmp_conn_nums$conn[i],c("corr")]<-min(tmp$corr)#Replace all corrs by the min corr, it'll be easier to get rid of duplications that way  
    }
  }
  
  univar_corrs<-univar_corrs %>% distinct(conn, .keep_all = TRUE)
  
  #Get the rownames of the top/bottom 200 corrs
  if (dir=="top"){
    top.corrs.idx<-univar_corrs %>% rownames_to_column() %>% top_n(Nfeat, corr) %>% pull(rowname)
    
    top.corrs.X<-univar_corrs$conn[as.numeric(top.corrs.idx)]
    
    #Get the raw feat corresponding to top corrs
    top.corr.raw<-Xvar[,colnames(Xvar) %in% univar_corrs$conn[as.numeric(top.corrs.idx)]]  
  } else if (dir=="bottom"){
    bot.corrs.idx<-univar_corrs %>% rownames_to_column() %>% top_n(-topn, corr) %>% pull(rowname)
    #NOTE: that's how you select the bottom, it's top_n but put a negative for the number
    bot.corrs.X<-univar_corrs$conn[as.numeric(bot.corrs.idx)]
    
    #Get the raw feat corresponding to bottom corrs
    bot.corr.raw<-Xvar[,colnames(Xvar) %in% univar_corrs$conn[as.numeric(bot.corrs.idx)]]
  }
  
}
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
get_pval <- function(real, null_dist, better="smaller"){
  if (better == "smaller"){
    rank <- sum(real < null_dist) + 1
  }
  if (better == "bigger"){
    rank <- sum(real > null_dist) + 1
  }
  pval <- rank / (length(null_dist) + 1)
  return(pval)
}
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
plsc_and_permute<-function(X,Y,num_comps,nperms,metrics,topn,out_flag){
  
  ##
  registerDoMC(cores=2) # to run it multicore
  
  print("Univariate feature selection")
  
  top.corr.raw<-uvar_feat_sel(Y,X,topn,"top")
  
  print("performing PLSc")
  real_model<-tepPLS(top.corr.raw,Y,k = min(dim(top.corr.raw)[2],dim(Y)[2]),center1 = TRUE,center2 = TRUE,scale1 = TRUE,scale2 = TRUE,graphs = FALSE)
  real_eigs<-real_model$TExPosition.Data$eigs
  
  shuffled_xindexes <- sapply(1:nperms, function(x){
    shuffle(1:nrow(X))
  }
  )
  null_results <- foreach(i=1:nperms) %dopar% {
    null_model <- apply_pls(top.corr.raw[shuffled_xindexes[,i],],Y,ncomps=num_comps,scaleflag = T,permflag = T,outflag = 1)
    list(null_model)
  }
  
  # transform null results lists to data frame
  null_outputs <- lapply(null_results, function(x){return(x[[1]])})
  
  null_outputs <- as.data.frame(do.call(rbind, null_outputs))
  
  print("computing p-values")
  pvals_cancor <- mapply(function(real, null_dist){
    get_pval(real, null_dist, better="smaller")},
    real_eigs,
    null_outputs)
  
  if (out_flag==1){
    return(list(pvals_cancor,real_eigs,null_outputs,real_model))
  } else if (out_flag==2){
    return(real_model)
  }
}
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
plsc_and_permute_ver2<-function(X,Y,nperms,num_comps,metrics,topn,out_flag,pls_flag){
  
  #Function to test what happens if we permute BEFORE the univariate feature selection step. Part of the sensitivity analyses.
  
  registerDoMC(cores=2) # to run it multicore
  
  print("Univariate feature selection")
  top.correlations.raw<-uvar_feat_sel(Y,X,topn,"top") #Note: I made sure that the name is different than what is defined in the function
  
  print("performing PLSc")
  real_model<-tepPLS(top.correlations.raw,Y,k = min(dim(top.correlations.raw)[2],dim(Y)[2]),center1 = TRUE,center2 = TRUE,scale1 = TRUE,scale2 = TRUE,graphs = FALSE)
  real_eigs<-real_model$TExPosition.Data$eigs
  
  print("permutation testing")
  
  
  shuffled_xindexes <- sapply(1:nperms, function(x){
    shuffle(1:nrow(X))
  }
  )
  
  
  null_results <- foreach(i=1:nperms) %dopar% {
    top.correlations.shuffled<-uvar_feat_sel(Y,X[shuffled_xindexes[,i],],topn,"top") #SHUFFLE PRIOR TO UVAR SELECTION
    #Shuffling indexes for conns. So now the univar feat sel selects rows that are already shuffled
    null_model <- apply_pls(top.correlations.shuffled,Y,ncomps=num_comps,scaleflag = T,permflag = T,outflag = 1)
    list(null_model)
  }
  
  # transform null results lists to data frame
  null_eigs <- lapply(null_results, function(x){return(x[[1]])})
  null_eigs <- as.data.frame(do.call(rbind, null_eigs))
  
  print("computing p-values")
  pvals_eigs <- mapply(function(real, null_dist){
    get_pval(real, null_dist, better="smaller")},
    real_eigs,
    null_eigs)
  
  if (out_flag==1){
    return(list(pvals_eigs,real_eigs,null_eigs,real_model))
  } else if (out_flag==2){
    return(real_model)
  }
}
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
get_cormat_pvals<-function(x,y){
  pval_mat<-matrix(data=NA,nrow=dim(y)[2],ncol = dim(x)[2])
  rownames(pval_mat)<-colnames(y)
  colnames(pval_mat)<-colnames(x)
  for (i in 1:dim(x)[2]){
    for (j in 1:dim(y)[2]){
      tmpcor<-cor.test(x[,i],y[,j])
      pval_mat[j,i]<-tmpcor$p.value
    }
  }
  return(pval_mat) 
}
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
create_polar_plots<-function(indat,var,color,colorpal,ylims,boot,lbldat){
  #Note: the variables need to be stored in a column called "variable"
  if (boot==TRUE){
    lbl_dat<-merge(lbldat,indat[indat$variable==var,c("gen_prob","value","boottest")])
    
    # calculate the alignment of labels: right or left
    # If I am on the left part of the plot, my labels have currently an angle < -90
    lbl_dat$hjust<-ifelse( angle < -90, 1, 0)
    
    # flip angle BY to make them readable
    lbl_dat$angle<-ifelse(angle < -90, angle+180, angle)    
  }else {
    lbl_dat<-merge(label_data,indat[indat$variable==var,c("gen_prob","value")])
  }
  
  plt<-ggplot(data=indat[indat$variable==var,],aes(x=gen_prob,y=value,fill=gen_prob)) +
    geom_bar(position = 'dodge', stat='identity')+scale_fill_manual(values=colorpal)+
    coord_polar(start=0)+ylim(ylims)+
    geom_bar(stat = "identity", aes(x = gen_prob, y = 0.01), fill = "white") +
    geom_bar(stat = "identity", aes(x = gen_prob, y = -0.04), fill = "white") +
    geom_bar(stat = "identity", aes(x = gen_prob, y = -0.03), fill = color) +
    theme(legend.position="none",
          axis.text = element_blank(),
          axis.title.y = element_blank(),
          panel.grid = element_blank(),
          axis.ticks=element_blank(),
          panel.background = element_rect(fill="white"))+
    xlab(var)
  if (boot==TRUE){
    plt<-plt+geom_text(data=lbl_dat, aes(x=id, y=ifelse(value<0,value-.05,value+.05), label=ifelse(boottest==TRUE,"*",""), hjust=hjust), color="black",size=2, angle= lbl_dat$angle, inherit.aes = FALSE )
  }
  return(plt)
}


####Import Data####
####################################################SYMPTOMS####################################################################
##    1. Headache - CBCL Headaches (cbcl_q56b_p)                                                                              ##
##    2. Nausea - CBCL Nausea, feels sick (cbcl_q56c_p)                                                                       ##
##    3. Vomiting - CBCL Vomiting, throwing up (cbcl_q56g_p)                                                                  ##
##    4. Dizziness - CBCL Feels dizzy or lightheaded (cbcl_q51_p)                                                             ##
##    5. Fatigue - CBCL Overtired without good reason (cbcl_q54_p)                                                            ##
##    6. Drowsiness - SDS Daytime sleepiness (sleepdisturb25_p)                                                               ##
##    7. Trouble falling asleep - SDS Difficulty getting to sleep (sleepdisturb4_p)                                           ##
##    8. Sleep more than usual - CBCL Sleeps more than most kids (cbcl_q77_p)                                                 ##
##    9. Sleep less than usual - CBCL Sleeps less than most kids (cbcl_q76_p)                                                 ##
##    10. Sadness - CBCL Depression (DSM) T score (cbcl_scr_dsm5_depress_t)                                                   ##
##    11. Nervousness - CBCL Anxiety Disorder (DSM) T score (cbcl_scr_dsm5_anxdisord_t)                                       ##
##    12. Diffculty concentrating - CBCL Attention problems T score (cbcl_scr_syn_attention_t)                                ##
##    13. Picture Sequence Memory Test Fully Corrected T score (nihtbx_picture_fc)                                            ##
##    14. List Sorting Working Memory Test – Fully-Corrected T-score (nihtbx_list_fc)                                         ##
##    15. RAVLT Short Delay Trial VI – Total Correct (pea_ravlt_sd_trial_vi_tc)                                               ##
##    16. RAVLT Long Delay Trial VII – Total Correct (pea_ravlt_ld_trial_vii_tc)                                              ##
##    17. Processing speed - Pattern Comparison Processing Speed Test – Fully-Corrected T-score (nihtbx_pattern_fc)           ##
##    18. Executive function - Dimensional Change Card Sort Test – Fully-Corrected T-score (nihtbx_cardsort_fc)               ##
##    19. Aggression - CBCL Aggression T score (cbcl_scr_syn_aggressive_t)                                                    ##
################################################################################################################################
####IMPORT SUBS EXCLUDED BY POST-PROC QC####
subs_excl<-read.table('/Users/Guido/Desktop/Projects/PCHStudy/Data/TBI_QC_fail_final.txt', sep='\t', header=F)
colnames(subs_excl)<-"subjectkey"
subs_excl$subjectkey<-switch_names(subs_excl$subjectkey,Rem = "sub-NDAR",Add = "NDAR_",separ = "",side="L")
dim(subs_excl)


####IMPORT SYMPTOMS####
#@@@@@@@@@@CBCL Summary@@@@@@@@@@@
##12, 13, 15, 22
cbcl_summary <- read.table('/Users/Guido/Desktop/GeneralData/ABCD/GeneralData/abcd_cbcls01.txt', sep='\t', header=T)#for troubleshooting
datadic_cbcl_summary <- cbcl_summary[1,]
datadic_cbcl_summary <-t(datadic_cbcl_summary)
cbcl_summary <-cbcl_summary[-1,]

cbcl_summary_clean<-cbcl_summary[cbcl_summary$eventname=="baseline_year_1_arm_1",c("subjectkey","sex","cbcl_scr_dsm5_depress_t",
                                                                                   "cbcl_scr_dsm5_anxdisord_t","cbcl_scr_syn_attention_t",
                                                                                   "cbcl_scr_syn_aggressive_t")]
dim(cbcl_summary_clean)

#@@@@@@@@@@CBCL Raw@@@@@@@@@@@@@@
##1, 2, 3, 4, 5, 8, 9, 10, 11, 14, 16
cbcl_raw <- read.table('/Users/Guido/Desktop/GeneralData/ABCD/GeneralData/abcd_cbcl01.txt', sep='\t', header=T)#for troubleshooting
datadic_cbcl_raw <- cbcl_raw[1,]
datadic_cbcl_raw <-t(datadic_cbcl_raw)
cbcl_raw <-cbcl_raw[-1,]
cbcl_raw_clean<-cbcl_raw[cbcl_raw$eventname=="baseline_year_1_arm_1",c("subjectkey","sex","cbcl_q56b_p","cbcl_q56c_p","cbcl_q56g_p","cbcl_q51_p",
                                                                       "cbcl_q54_p","cbcl_q77_p","cbcl_q76_p")]
dim(cbcl_raw_clean)


#@@@@@@@@@@@Sleep Disturbance Scale@@@@@@@@@@
##6, 7
sds <- read.table('/Users/Guido/Desktop/GeneralData/ABCD/GeneralData/abcd_sds01.txt', sep='\t', header=T)#for troubleshooting
datadic_sds <- sds[1,]
datadic_sds <-t(datadic_sds)
sds <-sds[-1,]
sds_clean<-sds[sds$eventname=="baseline_year_1_arm_1",c("subjectkey","sex","sleepdisturb25_p","sleepdisturb4_p")]
dim(sds_clean)



#@@@@@@@@@@@@@Pearson@@@@@@@@@@@@@
##17.3, 17.4
pearson <- read.table('/Users/Guido/Desktop/GeneralData/ABCD/GeneralData/abcd_ps01.txt', sep='\t', header=T)#for troubleshooting
datadic_pearson <- pearson[1,]
datadic_pearson <-t(datadic_pearson)
pearson <-pearson[-1,]
pearson_clean<-pearson[pearson$eventname=="baseline_year_1_arm_1",c("subjectkey","sex","pea_ravlt_sd_trial_vi_tc","pea_ravlt_ld_trial_vii_tc")]
dim(pearson_clean)


#@@@@@@@@@@NIH ToolBox Summary Scores@@@@@@@@@@@@@@
##17.1, 17.2, 20, 21
tbss <- read.table('/Users/Guido/Desktop/GeneralData/ABCD/GeneralData/abcd_tbss01.txt', sep='\t', header=T)#for troubleshooting
datadic_tbss <- tbss[1,]
datadic_tbss <-t(datadic_tbss)
tbss <-tbss[-1,]
tbss_clean<-tbss[tbss$eventname=="baseline_year_1_arm_1",c("subjectkey","sex","nihtbx_picture_fc","nihtbx_list_fc",
                                                           "nihtbx_pattern_fc","nihtbx_cardsort_fc")]
dim(tbss_clean)

#############MERGE SYMPTOMS##############
tmp1<-merge(cbcl_summary_clean,cbcl_raw_clean,by="subjectkey")
dim(tmp1)
tmp2<-merge(tmp1,sds_clean,by="subjectkey")
dim(tmp2)
tmp3<-merge(tmp2,pearson_clean,by="subjectkey")
dim(tmp3)
symptoms_all<-merge(tmp3,tbss_clean,by="subjectkey")
dim(symptoms_all)

#############IMPORT mTBI SUB LIST (processed ones only)###################
# mtbi_subs <- read.table('/Users/guigub/Documents/Research/PhD/Projects/BiotypingStudy/Data/mTBI_IDs.txt', sep='\t', header=T)#for troubleshooting
mtbi_subs <- read.table('/Users/Guido/Desktop/Projects/PCHStudy/Data/mTBI_IDs.txt', sep='\t', header=F)#for troubleshooting
# mtbi_subs<-as.data.frame(mtbi_subs[-1,])
colnames(mtbi_subs)<-"subjectkey"
dim(mtbi_subs)


mtbi_subs$subjectkey<-switch_names(mtbi_subs$subjectkey,"sub-NDAR","NDAR","_","L")
symptoms_mTBI<-symptoms_all[symptoms_all$subjectkey %in% mtbi_subs$subjectkey,]
dim(symptoms_mTBI)

symptoms_mTBI<-symptoms_mTBI[! symptoms_mTBI$subjectkey %in% subs_excl$subjectkey,]
dim(symptoms_mTBI)
dim(subs_excl)

# View(summary(symptoms_mTBI))
symptoms_mTBI<-symptoms_mTBI[,-which(names(symptoms_mTBI) %in% c("sex.x","sex.y"))]

#Convert symptoms to numeric
for (i in 1:dim(symptoms_mTBI)[2]){
  if (colnames(symptoms_mTBI)[i]!="subjectkey"&colnames(symptoms_mTBI)[i]!="sex"){
    symptoms_mTBI[,i]<-as.numeric(as.character(symptoms_mTBI[,i]))
  }
}


####REVERSE CODING SYMPTOMS (ON ENTIRE DATASET)####
symptoms_mTBI$pea_ravlt_sd_trial_vi_tc<-max(symptoms_mTBI$pea_ravlt_sd_trial_vi_tc,na.rm = T)-symptoms_mTBI$pea_ravlt_sd_trial_vi_tc
symptoms_mTBI$pea_ravlt_ld_trial_vii_tc<-max(symptoms_mTBI$pea_ravlt_ld_trial_vii_tc,na.rm = T)-symptoms_mTBI$pea_ravlt_ld_trial_vii_tc
symptoms_mTBI$nihtbx_picture_fc<-max(symptoms_mTBI$nihtbx_picture_fc,na.rm = T)-symptoms_mTBI$nihtbx_picture_fc
symptoms_mTBI$nihtbx_list_fc<-max(symptoms_mTBI$nihtbx_list_fc,na.rm = T)-symptoms_mTBI$nihtbx_list_fc
symptoms_mTBI$nihtbx_pattern_fc<-max(symptoms_mTBI$nihtbx_pattern_fc,na.rm = T)-symptoms_mTBI$nihtbx_pattern_fc
symptoms_mTBI$nihtbx_cardsort_fc<-max(symptoms_mTBI$nihtbx_cardsort_fc,na.rm = T)-symptoms_mTBI$nihtbx_cardsort_fc


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@STRUCTURAL CONN FEATURES@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
####LOAD LABEL LIST####
#Load complete list of lbls (will be required often)
lbls<-read.table('/Users/Guido/Desktop/Projects/PCHStudy/Data/DKT_lbls.txt',fill = T,header = T)
colnames(lbls)<-c("lbls","num","full_name")

####LOAD ALL STR CONN FEATURES####
str_conn_feats<-read.csv('/Users/Guido/Desktop/Projects/PCHStudy/Data/mTBI_results_for_Biotyping_study/post_proc_results/conn_features.csv')
str_conn_feats$sub<-switch_names(str_conn_feats$sub,Rem = "_ses",Add="",separ="",side="R")
str_conn_feats$sub<-switch_names(str_conn_feats$sub,Rem = "sub-NDAR",Add = "NDAR_",separ = "",side="L")
head(str_conn_feats)
tail(str_conn_feats)


####LOAD NUISANCE DATS: MRI SITES####
MRI <- read.table('/Users/Guido/Desktop/GeneralData/ABCD/GeneralData/abcd_mri01.txt', sep='\t', header=T)
MRI <- MRI[c('src_subject_id', 'mri_info_deviceserialnumber')]
MRI<-MRI[-1,]
colnames(MRI)<-c("subjectkey","mri_info_deviceserialnumber")


message(paste("There are currently",dim(data.frame(unique(str_conn_feats$sub)))[1],"subs",sep=" "))

####ONLY KEEP CONNS FROM SUBS THAT PASSED QC####
str_conn_feats_QCpass<-str_conn_feats[! str_conn_feats$sub %in% subs_excl$subjectkey,]
dim(str_conn_feats_QCpass)
####ONLY KEEP SYMPTOMS FORM SUBS THAT PASSED QC####
symptoms_mTBI<-symptoms_mTBI[symptoms_mTBI$subjectkey %in% str_conn_feats_QCpass$sub,]
dim(symptoms_mTBI)

####DEFINE METRICS####
mtrs<-as.data.frame(unique(str_conn_feats$metric))
colnames(mtrs)<-c("metric")


####INCORPORATE TBI CHARACTERISTICS####
#Incorporate t_since_injury and Sex in MRI to be able to regress them out
#Incorporate information to improve interpretability
TBI_info<-read.table('/Users/Guido/Desktop/GeneralData/ABCD/GeneralData/abcd_tbi01.txt', sep='\t', header=T)#for troubleshooting
datadic_tbi <- TBI_info[1,]
datadic_tbi <-t(datadic_tbi)
TBI_info <-TBI_info[-1,]

#tbi_ss_ntbiloc: n of TBIs with LOC
#tbi_ss_nmrpi n of periods with multiple TBIs
#tbi_ss_agefirst: age at first TBI
#tbi_ss_ntbiloc30: n TBIs with LOC >-30mins

TBI_info_clean<-TBI_info[TBI_info$eventname=="baseline_year_1_arm_1",c("subjectkey","sex","tbi_ss_ntbiloc","tbi_ss_nmrpi","tbi_ss_agefirst","tbi_ss_ntbiloc30")]
dim(TBI_info_clean)

TBI_info_clean_study<-TBI_info_clean[ TBI_info_clean$subjectkey %in% symptoms_mTBI$subjectkey,]
summary(TBI_info_clean_study)
dim(TBI_info_clean_study)

TBI_info_raw<-read.table('/Users/Guido/Desktop/GeneralData/ABCD/GeneralData/abcd_otbi01.txt', sep='\t', header=T)#for troubleshooting
datadic_tbi_raw <- TBI_info_raw[1,]
datadic_tbi_raw <-t(datadic_tbi_raw)
TBI_info_raw <-TBI_info_raw[-1,]

TBI_info_raw_study<-TBI_info_raw[TBI_info_raw$subjectkey %in% symptoms_mTBI$subjectkey,]
TBI_raw_ages<-TBI_info_raw_study[,c("subjectkey","tbi_1d","tbi_2d","tbi_3d","tbi_4d","tbi_5d")]
t_since_inj<-TBI_info_raw_study[,c("subjectkey","interview_age")]
t_since_inj$interview_age<-floor(as.numeric(as.character(t_since_inj$interview_age))/12)#way to get current age, given in months
for (j in 2:dim(TBI_raw_ages)[2]){
  TBI_raw_ages[,j]<-as.numeric(as.character(TBI_raw_ages[,j]))
}
TBI_raw_ages[is.na(TBI_raw_ages)]<- -1 #Used -1 to replace NAs because 0 is a plausible number here
t_since_inj$t_since_inj<-NA
for (i in 1:dim(TBI_raw_ages)[1]){
  t_since_inj[i,c("t_since_inj")]<-t_since_inj[i,c("interview_age")]-max(TBI_raw_ages[i,c("tbi_1d","tbi_2d","tbi_3d","tbi_4d","tbi_5d")])
  if (t_since_inj[i,c("t_since_inj")]<0){
    message(paste("subject",t_since_inj[i,c("subjectkey")],"has negative t since inj (",t_since_inj[i,c("interview_age")],
                  "interview age - ",max(TBI_raw_ages[i,c("tbi_1d","tbi_2d","tbi_3d","tbi_4d","tbi_5d")]),"age at most recent TBI), setting to 0",sep=" "))
    t_since_inj[i,c("t_since_inj")]<-0
  }
}


####INCORPORATE PUBERTAL STAGE####
pubertal_stage <- read.table('/Users/Guido/Desktop/GeneralData/ABCD/GeneralData/abcd_ssphp01.txt', sep='\t', header=T)
pubertal_stage<-pubertal_stage[-1,]
pubertal_stage<-subset(pubertal_stage,pubertal_stage$eventname=='baseline_year_1_arm_1')
pubertal_stage <- pubertal_stage[c('subjectkey','sex', 'pds_p_ss_female_category','pds_p_ss_male_category')]
# subset(pubertal_stage,pubertal_stage=="2 1")#somebody has this value
# subset(pubertal_stage,subjectkey=="NDAR_INV00HEV6HB")#it's this person

pubertal_stage$pubertal_stage<-paste(pubertal_stage$pds_p_ss_female_category,pubertal_stage$pds_p_ss_male_category)
pubertal_stage$pubertal_stage<-as.numeric(as.character(pubertal_stage$pubertal_stage))
pubertal_stage <- pubertal_stage[c('subjectkey','sex', 'pubertal_stage')]
pubertal_stage_study<-pubertal_stage[pubertal_stage$subjectkey %in% symptoms_mTBI$subjectkey,]


####INCORPORATE HANDEDNESS####
hand <- read.table('/Users/Guido/Desktop/GeneralData/ABCD/GeneralData/abcd_ehis01.txt', sep='\t', header=T)
hand <- hand[,c('subjectkey', 'ehi_y_ss_scoreb')]
hand_study<-hand[hand$subjectkey %in% symptoms_mTBI$subjectkey,]


####MERGE NUISANCE COVS####
nuisance_covs<-merge(MRI,t_since_inj,all.y=TRUE)

nuisance_covs<-merge(nuisance_covs,pubertal_stage_study)
nuisance_covs<-merge(nuisance_covs,hand_study)


####Preproc####

mydata_t09<-pre_proc(conn_data = str_conn_feats_QCpass,symp_data = symptoms_mTBI,nuisance_dat = nuisance_covs,metrics=c("fa","ad","md","rd","nufo","afdf"),
                     nuisance_vars = c("sex","ehi_y_ss_scoreb","mri_info_deviceserialnumber","pubertal_stage"),nuisance_vars_symps =c("sex","ehi_y_ss_scoreb","pubertal_stage"),
                     thresh = 0.85,tr_ratio = .70,out_dir = "/Users/Guido/Desktop/Projects/PCHStudy/Data/",resume_thresh = FALSE,resume_regout = FALSE,resume_scale = FALSE)

#Running after resuming thresholding
mydata_t09<-pre_proc(conn_data = str_conn_feats_QCpass,symp_data = symptoms_mTBI,nuisance_dat = nuisance_covs,metrics=c("fa","ad","md","rd","nufo","afdf"),
                     nuisance_vars = c("sex","ehi_y_ss_scoreb","mri_info_deviceserialnumber","pubertal_stage"),nuisance_vars_symps =c("sex","ehi_y_ss_scoreb","pubertal_stage"),
                     thresh = 0.9,tr_ratio = .70,out_dir = "/Users/Guido/Desktop/Projects/PCHStudy/Data/",resume_thresh = TRUE,thresh_file = '/Users/Guido/Desktop/Projects/PCHStudy/Data/conn_features_thresh_0.9',
                     resume_regout = FALSE,resume_scale = FALSE)

mydata_t085<-pre_proc(conn_data = str_conn_feats_QCpass,symp_data = symptoms_mTBI,nuisance_dat = nuisance_covs,metrics=c("fa","ad","md","rd","nufo","afdf"),
                      nuisance_vars = c("sex","ehi_y_ss_scoreb","mri_info_deviceserialnumber","pubertal_stage"),nuisance_vars_symps =c("sex","ehi_y_ss_scoreb","pubertal_stage"),
                      thresh = 0.85,tr_ratio = .70,out_dir = "/Users/Guido/Desktop/Projects/PCHStudy/Data/",resume_thresh = TRUE,thresh_file = '/Users/Guido/Desktop/Projects/PCHStudy/Data/conn_features_thresh_0.85',
                      resume_regout = FALSE,resume_scale = FALSE)

mydata_t1<-pre_proc(conn_data = str_conn_feats_QCpass,symp_data = symptoms_mTBI,nuisance_dat = nuisance_covs,metrics=c("fa","ad","md","rd","nufo","afdf"),
                    nuisance_vars = c("sex","ehi_y_ss_scoreb","mri_info_deviceserialnumber","pubertal_stage"),nuisance_vars_symps =c("sex","ehi_y_ss_scoreb","pubertal_stage"),
                    thresh = 1,tr_ratio = .70,out_dir = "/Users/Guido/Desktop/Projects/PCHStudy/Data/",resume_thresh = TRUE,thresh_file = '/Users/Guido/Desktop/Projects/PCHStudy/Data/conn_features_thresh_1',
                    resume_regout = FALSE,resume_scale = FALSE)

mydata_t095<-pre_proc(conn_data = str_conn_feats_QCpass,symp_data = symptoms_mTBI,nuisance_dat = nuisance_covs,metrics=c("fa","ad","md","rd","nufo","afdf"),
                      nuisance_vars = c("sex","ehi_y_ss_scoreb","mri_info_deviceserialnumber","pubertal_stage"),nuisance_vars_symps =c("sex","ehi_y_ss_scoreb","pubertal_stage"),
                      thresh = 0.95,tr_ratio = .70,out_dir = "/Users/Guido/Desktop/Projects/PCHStudy/Data/",resume_thresh = TRUE,thresh_file = '/Users/Guido/Desktop/Projects/PCHStudy/Data/conn_features_thresh_0.95',
                      resume_regout = FALSE,resume_scale = FALSE)

mydata_t09_backup<-mydata_t09###@@@@@@@@@@@@@@@@@RUN TO BACK UP###
mydata_t09<-mydata_t09_backup##@@@@@@@@@@@@@@@@@@RUN TO RECOVER T09###

mydata_t09<-mydata_t1##@@@@@@@@@@@@@@@@@@@@RUN TO AVOID HAVING TO REWRITE EVERYTHING WITH A DIFFERENT THRESHOLD (useful later for sensitivity analyses###



####Clean myData and stack all metrics together####
#HERE I NEED TO FIND OUT IF ANY METRIC HAD LESS CONNS, AND ONLY KEEP THE SMALLEST NUMBER OF CONNS
extract_colnames<-function(input_colnames){
  # tmpbigcols<-substr(input_colnames,regexpr("\\.",input_colnames)+1,20)#for some reason, this line doesn't extract the conn number alone anymore
  tmpbigcols<-substr(input_colnames,regexpr("\\.",input_colnames)+1,20)
  tmpbigcols<-substr(tmpbigcols,1,regexpr("\\.",tmpbigcols)-1)
  return(tmpbigcols)
}


tmp_dims<-c(dim(mydata_t09$xx_train$fa)[2],dim(mydata_t09$xx_train$md)[2],dim(mydata_t09$xx_train$rd)[2],dim(mydata_t09$xx_train$ad)[2],dim(mydata_t09$xx_train$nufo)[2],dim(mydata_t09$xx_train$afdf)[2])
tmp_mtrs<-c("fa","md","rd","ad","nufo","afdf")
tmp_dim_mtrs<-rbind(tmp_dims,tmp_mtrs)
min_mtr<-tmp_dim_mtrs[rownames(tmp_dim_mtrs)=="tmp_mtrs",which(tmp_dim_mtrs[1,]==min(tmp_dim_mtrs[1,]))][1]
max_mtr<-tmp_dim_mtrs[rownames(tmp_dim_mtrs)=="tmp_mtrs",which(tmp_dim_mtrs[1,]==max(tmp_dim_mtrs[1,]))][1]
# tmpbigcols<-extract_colnames(input_colnames = colnames(as.data.frame(mydata_t09$xx_train[max_mtr])))

tmpbigcols<-substr(colnames(as.data.frame(mydata_t09$xx_train[max_mtr])),regexpr("\\.",colnames(as.data.frame(mydata_t09$xx_train[max_mtr])))+1,20)#for some reason, this line doesn't extract the conn number alone anymore
tmpbigcols<-extract_colnames(tmpbigcols)

# tmpbigcols<-substr(tmpbigcols,regexpr("\\.",tmpbigcols)+1,20)
# tmpbigcols<-substr(tmpbigcols,1,regexpr("\\.",tmpbigcols)-1)


tmp_misscols<-cbind(as.vector(tmpbigcols[!tmpbigcols%in%extract_colnames(colnames(mydata_t09$xx_train$fa))]),as.vector(tmpbigcols[!tmpbigcols%in%extract_colnames(colnames(mydata_t09$xx_train$md))]),
                    as.vector(tmpbigcols[!tmpbigcols%in%extract_colnames(colnames(mydata_t09$xx_train$rd))]),as.vector(tmpbigcols[!tmpbigcols%in%extract_colnames(colnames(mydata_t09$xx_train$ad))]),
                    as.vector(tmpbigcols[!tmpbigcols%in%extract_colnames(colnames(mydata_t09$xx_train$nufo))]),as.vector(tmpbigcols[!tmpbigcols%in%extract_colnames(colnames(mydata_t09$xx_train$afdf))]))
misscols<-unique(melt(tmp_misscols)$value)
tmpsmallcols<-tmpbigcols[!tmpbigcols%in%misscols]

# tmpsmallcols<-substr(colnames(as.data.frame(mydata_t09$xx_train[min_mtr])),regexpr("\\.",colnames(as.data.frame(mydata_t09$xx_train[min_mtr])))+1,20)


tmp_submat<-matrix(data=NA,ncol = length(tmpsmallcols),nrow = length(mydata_t09$yy_train$subjectkey))
colnames(tmp_submat)<-tmpsmallcols
tmp_submat<-as.data.frame(tmp_submat)
tmp_submat$subjectkey<-mydata_t09$yy_train$subjectkey
tmp_submat_melt<-melt(tmp_submat,id.vars = "subjectkey")

if (length(tmpbigcols)>length(tmpsmallcols)){
  #Note that here I'm removing those "extra" conns that metrics other than NuFO had.
  
  length(mydata_t09$xx_train$fa[,which(extract_colnames(colnames(mydata_t09$xx_train$fa))%in%tmpsmallcols)])
  length(mydata_t09$xx_train$md[,which(extract_colnames(colnames(mydata_t09$xx_train$md))%in%tmpsmallcols)])
  length(mydata_t09$xx_train$rd[,which(extract_colnames(colnames(mydata_t09$xx_train$rd))%in%tmpsmallcols)])
  length(mydata_t09$xx_train$ad[,which(extract_colnames(colnames(mydata_t09$xx_train$ad))%in%tmpsmallcols)])
  length(mydata_t09$xx_train$nufo[,which(extract_colnames(colnames(mydata_t09$xx_train$nufo))%in%tmpsmallcols)])
  length(mydata_t09$xx_train$afdf[,which(extract_colnames(colnames(mydata_t09$xx_train$afd))%in%tmpsmallcols)])
  
  mydata_t09$xx_train$all_metrics<-melt(mydata_t09$xx_train$fa[,which(extract_colnames(colnames(mydata_t09$xx_train$fa))%in%tmpsmallcols)])
  mydata_t09$xx_train$all_metrics$md<-melt(mydata_t09$xx_train$md[,which(extract_colnames(colnames(mydata_t09$xx_train$md))%in%tmpsmallcols)])$value
  mydata_t09$xx_train$all_metrics$rd<-melt(mydata_t09$xx_train$rd[,which(extract_colnames(colnames(mydata_t09$xx_train$rd))%in%tmpsmallcols)])$value
  mydata_t09$xx_train$all_metrics$ad<-melt(mydata_t09$xx_train$ad[,which(extract_colnames(colnames(mydata_t09$xx_train$ad))%in%tmpsmallcols)])$value
  mydata_t09$xx_train$all_metrics$nufo<-melt(mydata_t09$xx_train$nufo[,which(extract_colnames(colnames(mydata_t09$xx_train$nufo))%in%tmpsmallcols)])$value
  mydata_t09$xx_train$all_metrics$afdf<-melt(mydata_t09$xx_train$afdf[,which(extract_colnames(colnames(mydata_t09$xx_train$afdf))%in%tmpsmallcols)])$value
  mydata_t09$xx_train$all_metrics$subjectkey<-tmp_submat_melt$subjectkey
  mydata_t09$xx_train$all_metrics$variable<-tmp_submat_melt$variable
  colnames(mydata_t09$xx_train$all_metrics)<-c("conn","fa","md","rd","ad","nufo","afdf","subjectkey")
  
} else if (length(tmpbigcols)==length(tmpsmallcols)){
  #Note that here I'm removing those "extra" conns that metrics other than NuFO had.
  mydata_t09$xx_train$all_metrics<-melt(mydata_t09$xx_train$fa)
  mydata_t09$xx_train$all_metrics$md<-melt(mydata_t09$xx_train$md)$value
  mydata_t09$xx_train$all_metrics$rd<-melt(mydata_t09$xx_train$rd)$value
  mydata_t09$xx_train$all_metrics$ad<-melt(mydata_t09$xx_train$ad)$value
  mydata_t09$xx_train$all_metrics$nufo<-melt(mydata_t09$xx_train$nufo)$value
  mydata_t09$xx_train$all_metrics$afdf<-melt(mydata_t09$xx_train$afdf)$value
  mydata_t09$xx_train$all_metrics$subjectkey<-tmp_submat_melt$subjectkey
  colnames(mydata_t09$xx_train$all_metrics)<-c("conn","fa","md","rd","ad","nufo","afdf","subjectkey")
}


####COMMIT####
#LOAD COMMIT2 WEIGHTS
commit1_weights<-read.csv('/Users/Guido/Desktop/Projects/PCHStudy/Data/mTBI_results_for_Biotyping_study/commit1_outputs/commit1_weights.csv')
head(commit1_weights)
tail(commit1_weights)

sum(is.na(commit1_weights$mean))


#Remove brainstem connections
commit1_weights$start<-substr(commit1_weights$conn,1,regexpr("_",commit1_weights$conn)-1)
commit1_weights$end<-substr(commit1_weights$conn,regexpr("_",commit1_weights$conn)+1,20)
commit1_weights$conn

commit1_weights_nobs<-commit1_weights[!(commit1_weights$start==16|commit1_weights$end==16),]
length(unique(commit1_weights_nobs$conn)) #2850
#This number is the number of possible permutations given 76 labels, divided by 2 given that reverse order is the same


##Now filter all connections by COMMIT1 weights
#1. Check how many have non-zero c1 weights
selected_conns<-unique(mydata_t09$xx_train$all_metrics$conn)
selected_conns<-as.data.frame(selected_conns)
dim(selected_conns)#1026

#Select only the connections that passed the thresholding
commit1_weights_nobs<-commit1_weights_nobs[commit1_weights_nobs$conn%in%selected_conns$selected_conns,]
length(unique(commit1_weights_nobs$conn))#1026, same number of conns selected after thresholding

#Set the right levels
commit1_weights_nobs$conn<-factor(commit1_weights_nobs$conn,levels=levels(selected_conns$selected_conns))

selected_conns$c1w_pos<-0
crit=90
for (i in 1:dim(selected_conns)[1]){
  tmp<-commit1_weights_nobs[commit1_weights_nobs$conn==selected_conns[i,c("selected_conns")],]
  if (sum(tmp$mean>0)/length(unique(commit1_weights_nobs$sub))*100>=crit){
    selected_conns[i,c("c1w_pos")]<-1
  }  
}
sum(selected_conns$c1w_pos)#629/1026
###IMPORTANT NOTE: THIS COMMIT WAS APPLIED ACROSS THE ENTIRE 306 SUBJECTS, NOT DONE SEPARATELY BY SET!!!###

####PCA####
pca.across.metrics<-prcomp(mydata_t09$xx_train$all_metrics[mydata_t09$xx_train$all_metrics$conn%in%selected_conns[selected_conns$c1w_pos==1,c("selected_conns")],
                                                           -which(colnames(mydata_t09$xx_train$all_metrics)%in%c("subjectkey","conn"))],center = TRUE,scale. = TRUE)

corrplot(pca.across.metrics$rotation, is.corr=FALSE)
summary(pca.across.metrics)
#0.8089 0.1626

####EXTRACT FIRST 2 PCs####
#Extract first 2 PCS, they account for 97% of variance
PC1.scores<-as.data.frame(pca.across.metrics$x[,1])
colnames(PC1.scores)<-"score"
PC2.scores<-as.data.frame(pca.across.metrics$x[,2])
colnames(PC2.scores)<-"score"
PC1.scores$conns<-mydata_t09$xx_train$all_metrics[mydata_t09$xx_train$all_metrics$conn%in%selected_conns[selected_conns$c1w_pos==1,c("selected_conns")],c("conn")]
PC1.scores$subjectkey<-mydata_t09$xx_train$all_metrics[mydata_t09$xx_train$all_metrics$conn%in%selected_conns[selected_conns$c1w_pos==1,c("selected_conns")],c("subjectkey")]
PC2.scores$conns<-mydata_t09$xx_train$all_metrics[mydata_t09$xx_train$all_metrics$conn%in%selected_conns[selected_conns$c1w_pos==1,c("selected_conns")],c("conn")]
PC2.scores$subjectkey<-mydata_t09$xx_train$all_metrics[mydata_t09$xx_train$all_metrics$conn%in%selected_conns[selected_conns$c1w_pos==1,c("selected_conns")],c("subjectkey")]


####RESTRUCTURE PCS BACK TO WIDE FORMAT####
PC1.scores.wide<-reshape(PC1.scores,idvar = "subjectkey",timevar = "conns",direction = "wide")
colnames(PC1.scores.wide)<-substr(colnames(PC1.scores.wide),regexpr("\\.",colnames(PC1.scores.wide))+1,25)

PC2.scores.wide<-reshape(PC2.scores,idvar = "subjectkey",timevar = "conns",direction = "wide")
colnames(PC2.scores.wide)<-substr(colnames(PC2.scores.wide),regexpr("\\.",colnames(PC2.scores.wide))+1,25)



####PLSC####
plsc.t09.pc1<-plsc_and_permute(PC1.scores.wide[,-which(colnames(PC1.scores.wide)%in%c("subjectkey"))],mydata_t09$yy_train[,-which(colnames(mydata_t09$yy_train)%in%c("subjectkey","sex"))],
                               num_comps = 20,nperms = 2000,topn = 200,out_flag = 1)

plsc.t09.pc2<-plsc_and_permute(PC2.scores.wide[,-which(colnames(PC2.scores.wide)%in%c("subjectkey"))],mydata_t09$yy_train[,-which(colnames(mydata_t09$yy_train)%in%c("subjectkey","sex"))],
                               num_comps = 20,nperms = 2000,topn = 200,out_flag = 1)


####BootstrapTests####
boot.pc1<-Boot4PLSC(PC1.scores.wide[,which(colnames(PC1.scores.wide)%in%c(rownames(plsc.t09.pc1[[4]]$TExPosition.Data$pdq$p)))],mydata_t09$yy_train[,-which(colnames(mydata_t09$yy_train)%in%c("subjectkey","sex"))],
                    center1 = TRUE,center2 = TRUE,scale1 = TRUE,scale2 = TRUE,nf2keep = 19,nIter = 2000,critical.value = 1.96,eig = TRUE)
boot.pc1$bootRatios.j[, plsc.t09.pc1[[1]]<.05]
boot.pc1$bootRatiosSignificant.j[, plsc.t09.pc1[[1]]<.05]

boot.pc2<-Boot4PLSC(PC2.scores.wide[,which(colnames(PC2.scores.wide)%in%c(rownames(plsc.t09.pc2[[4]]$TExPosition.Data$pdq$p)))],mydata_t09$yy_train[,-which(colnames(mydata_t09$yy_train)%in%c("subjectkey","sex"))],
                    center1 = TRUE,center2 = TRUE,scale1 = TRUE,scale2 = TRUE,nf2keep = 19,nIter = 2000,critical.value = 1.96,eig = TRUE)
boot.pc2$bootRatios.j[, plsc.t09.pc2[[1]]<.05]
boot.pc2$bootRatiosSignificant.j[, plsc.t09.pc2[[1]]<.05]


####IncorporateExternalVariables####
#@@@@@@@@@@@@@@@@@@@@@@START INCORPORATING EXTERNAL VARIABLES@@@@@@@@@@@@@@@@@@@@@@
####KSADS####
# ## ## ## ## ## ## ## ## KSADS DIAGNOSES USED IN "ADVERSE PSYCH OUTCOME" # ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
# ksads_1_843_p: Diagnosis - Persistent Depressive Disorder (Dysthymia) PRESENT F34.1
# ksads_1_840_p: Diagnosis - Major Depressive Disorder Present
# ksads_1_846_p: Diagnosis - Unspecified Depressive Disorder Current (F32.9)
# ksads_2_835_p: Diagnosis - Bipolar II Disorder currently hypomanic ξF31.81 
# ksads_2_836_p: Diagnosis - Bipolar II Disorder currently depressed F31.81 
# ksads_2_831_p: Diagnosis - Bipolar I Disorder current episode depressed F31.3x 
# ksads_2_832_p: Diagnosis - Bipolar I Disorder currently hypomanic  F31.0 
# ksads_2_830_p: Diagnosis - Bipolar I Disorder current episode manic (F31.1x) 
# ksads_2_838_p: Diagnosis - Unspecified Bipolar and Related Disorder current ξ(F31.9) 
# ksads_3_848_p: Diagnosis - Disruptive Mood Dysregulation Disorder (DMDD) Current (F34.8) 
# ksads_4_851_p: Diagnosis - Unspecified Schizophrenia Spectrum and Other Psychotic Disorder F29 (Current) 
# ksads_4_826_p: Diagnosis - Hallucinations (Present) 
# ksads_4_828_p: Diagnosis - Delusions (Present)
# ksads_4_849_p: Diagnosis - Associated Psychotic Symptoms \u008a\u0097\u0096 Current 
# ksads_5_906_p: Diagnosis - Other Specified Anxiety Disorder (Panic Disorder impairment does not meet full criteria) F41.8 
# ksads_5_857_p: Diagnosis - Panic Disorder (F41.0) PRESENT 
# ksads_6_908_p: Diagnosis - Other Specified Anxiety Disorder (Agoraphobia impairment does not meet full criteria) F41.8 
# ksads_6_859_p: Diagnosis - Agoraphobia (F40.00) PRESENT 
# ksads_7_861_p: Diagnosis - Separation Anxiety Disorder (F93.00) PRESENT 
# ksads_7_909_p: Diagnosis - Other Specified Anxiety Disorder (Separation Anxiety Disorder impairment does not meet full criteria) F41.8 
# ksads_8_863_p: Diagnosis - Social Anxiety Disorder (F40.10) PRESENT 
# ksads_8_911_p: Diagnosis - Other Specified Anxiety Disorder (Social Anxiety Disorder impairment does not meet minimum duration) F41.8
# ksads_9_867_p: Diagnosis - Specific Phobia ξPRESENT (F40.2XX)
# ksads_10_913_p: Diagnosis - Other Specified Anxiety Disorder (Generalized Anxiety Disorder impairment does not meet minimum duration) F41.8 
# ksads_10_869_p: Diagnosis - Generalized Anxiety Disorder \u008a\u0097\u0096Present (F41.1)
# ksads_11_917_p: Diagnosis - Obsessive-Compulsive Disorder \u008a\u0097\u0096Present (F42)
# ksads_11_919_p: Diagnosis - Other Specified Obsessive-Compulsive and Related Disorder present  does not meet full criteria (F42)
# ksads_12_927_p: Diagnosis - Encopresis Present (F98.1) 
# ksads_12_925_p: Diagnosis - Enuresis Present (F98.0)
# ksads_13_938_p: Diagnosis - Binge-Eating Disorder (F50.8) CURRENT 
# ksads_13_929_p: Diagnosis - Anorexia Nervosa (F50.02) ξBinge eating/purging subtype PRESENT 
# ksads_13_932_p: Diagnosis - Anorexia Nervosa (F50.01) ξRestricting subtype PRESENT
# ksads_13_935_p: Diagnosis - Bulimia Nervosa (F50.2) PRESENT
# ksads_13_942_p: Diagnosis - Other Specified Feeding or Eating Disorder Bulimia Nervosa current does not meet full criteria (F50.8) 
# ksads_13_944_p: Diagnosis - Other Specified Feeding or Eating Disorder Binge Eating Disorder present does not meet full criteria (F50.8) 
# ksads_13_941_p: Diagnosis - Other Specified Feeding or Eating Disorder Anorexia Nervosa current does not meet full criteria (F50.8) 
# ksads_14_856_p: Diagnosis - Unspecified Attention-Deficit/Hyperactivity Disorder (F90.9)
# ksads_14_853_p: Diagnosis - Attention-Deficit/Hyperactivity Disorder Present
# ksads_15_901_p: Diagnosis - Oppositional Defiant Disorder Present F91.3
# ksads_16_897_p: Diagnosis - Conduct Disorder present childhood onset (F91.1)
# ksads_16_898_p: Diagnosis - Conduct Disorder present adolescent onset (F91.2)
# ksads_17_904_p: Diagnosis - Unspecified Tic Disorder present ξ(F95.9) 
# ksads_18_903_p: Diagnosis - Other Specified Neurodevelopmental Disorder Autism Spectrum Disorder full criteria not assessed (F88.0) 
# ksads_19_895_p: Diagnosis - Unspecified Alcohol-Related Disorder present (F10.99)
# ksads_19_891_p: Diagnosis - Alcohol Use Disorder Present 
# ksads_20_893_p: Diagnosis - Unspecified Substance Related Disorder Present (F10.99)
# ksads_20_874_p: Diagnosis - Stimulant Use Disorder Present: Cocaine 
# ksads_20_872_p: Diagnosis - Stimulant Use Disorder Present: ξAmphetamine-type substance
# ksads_20_889_p: Diagnosis - Substance Use Disorder CURRENT
# ksads_20_878_p: Diagnosis - Inhalant Use Disorder Present 
# ksads_20_877_p: Diagnosis - Phencycllidine (PCP) Use Disorder Present 
# ksads_20_875_p: Diagnosis - Opiod Use Disorder Present 
# ksads_20_876_p: Diagnosis - Other Hallucinagen Use Disorder Present 
# ksads_20_879_p: Diagnosis - Other Drugs Use Disorder Present 
# ksads_20_873_p: Diagnosis - Sedative Hypnotic or Anxiolytic Use Disorder Present 
# ksads_20_871_p: Diagnosis - Cannabis Use Disorder Present 
# ksads_21_923_p: Diagnosis - Other Specified Trauma-and Stressor-Related Disorder present (PTSD impairment does not meet full criteria (F43.8) 
# ksads_21_921_p: Diagnosis - Post-Traumatic Stress Disorder PRESENT (F94.1) 
# ksads_22_969_p: Diagnosis - SLEEP PROBLEMS Present 
# ksads_23_946_p: Diagnosis - SuicidalideationPassivePresent 
# ksads_23_954_p: Diagnosis - SuicideAttemptPresent 
# ksads_23_945_p: Diagnosis - SelfInjuriousBehaviorwithoutsuicidalintentPresent 
# ksads_23_950_p: Diagnosis - SuicidalideationActiveplanPresent 
# ksads_23_947_p: Diagnosis - SuicidalideationActivenonspecificPresent 
# ksads_23_948_p: Diagnosis - SuicidalideationActivemethodPresent 
# ksads_23_949_p: Diagnosis - SuicidalideationActiveintentPresent 
# ksads_23_952_p: Diagnosis - InterruptedAttemptPresent 
# ksads_23_955_p: Diagnosis - NosuicidalideationorbehaviorPresent 
# ksads_23_951_p: Diagnosis - PreparatoryActionstowardimminentSuicidalbehaviorPresent 
# ksads_23_953_p: Diagnosis - AbortedAttemptPresent 
# ksads_24_967_p: Diagnosis - HOMICIDAL IDEATION AND BEHAVIOR Present 
# ksads_25_915_p: Diagnosis - Other Specified Anxiety Disorder (Selective Mutism  does not meet minimum duration) F41.8 
# ksads_25_865_p: Diagnosis - Selective Mutism (F94.0) PRESENT 



ksads <- read.table('/Users/Guido/Desktop/GeneralData/ABCD/GeneralData/abcd_ksad01.txt', sep='\t', header=T)
datadic_ksads_summary <- ksads[1,]
datadic_ksads_summary <-t(datadic_ksads_summary)
ksads <-ksads[-1,]

datadic_ksads_summary<-as.data.frame(datadic_ksads_summary)

# ksads_15_901_p = Diagnosis - Oppositional Defiant Disorder Present F91.3
# ksads_15_902_p = Diagnosis - Oppositional Defiant Disorder Past F91.3
# ksads_15_97_p = Duration (at least 6 months) ODD, Past
# ksads_16_900_p = Diagnosis - Conduct Disorder past adolescent onset (F91.2)
# ksads_16_897_p = Diagnosis - Conduct Disorder present childhood onset (F91.1)
# ksads_16_899_p = Diagnosis - Conduct Disorder past childhood onset (F91.1)
# ksads_16_898_p = Diagnosis - Conduct Disorder present adolescent onset
#extract baseline
ksads_clean<-ksads[ksads$eventname=="baseline_year_1_arm_1",] 
#extract diagnoses
ksads_dx<-ksads_clean[,regexpr("Diagnosis",datadic_ksads_summary[,1])>0]
ksads_dx$subjectkey<-ksads_clean$subjectkey
datadic_ksads_dx<-datadic_ksads_summary[regexpr("Diagnosis",datadic_ksads_summary[,1])>0,]

ksads_dx_present<-ksads_dx[,c("subjectkey","ksads_1_843_p","ksads_1_840_p","ksads_1_846_p","ksads_2_835_p","ksads_2_836_p","ksads_2_831_p","ksads_2_832_p","ksads_2_830_p","ksads_2_838_p","ksads_3_848_p",
                              "ksads_4_851_p","ksads_4_826_p","ksads_4_828_p","ksads_4_849_p","ksads_5_906_p","ksads_5_857_p","ksads_6_908_p","ksads_6_859_p","ksads_7_861_p","ksads_7_909_p","ksads_8_863_p","ksads_8_911_p",
                              "ksads_9_867_p","ksads_10_913_p","ksads_10_869_p","ksads_11_917_p","ksads_11_919_p","ksads_12_927_p","ksads_12_925_p","ksads_13_938_p","ksads_13_929_p","ksads_13_932_p","ksads_13_935_p",
                              "ksads_13_942_p","ksads_13_944_p","ksads_13_941_p","ksads_14_856_p","ksads_14_853_p","ksads_15_901_p","ksads_16_897_p","ksads_16_898_p","ksads_17_904_p","ksads_18_903_p","ksads_19_895_p",
                              "ksads_19_891_p","ksads_20_893_p","ksads_20_874_p","ksads_20_872_p","ksads_20_889_p","ksads_20_878_p","ksads_20_877_p","ksads_20_875_p","ksads_20_876_p","ksads_20_879_p","ksads_20_873_p",
                              "ksads_20_871_p","ksads_21_923_p","ksads_21_921_p","ksads_22_969_p","ksads_23_946_p","ksads_23_954_p","ksads_23_945_p","ksads_23_950_p","ksads_23_947_p","ksads_23_948_p","ksads_23_949_p",
                              "ksads_23_952_p","ksads_23_955_p","ksads_23_951_p","ksads_23_953_p","ksads_24_967_p","ksads_25_915_p","ksads_25_865_p")]


ksads_dx_present<-as.data.frame(ksads_dx_present)
current_dx<-apply(ksads_dx_present[,-which(colnames(ksads_dx_present)%in%c("subjectkey"))],1,function(x){sum(as.numeric(x),na.rm = TRUE)})

ksads_dx_present$current_dx<-current_dx
# ksads_dx_present$adverse_psych_outcome<-ifelse(ksads_dx_present$current_dx>0,1,0)
ksads_numeric<-ksads_dx_present[,!colnames(ksads_dx_present)%in%c("current_dx","adverse_psych_outcome")]
ksads_numeric<-ksads_numeric[,c("subjectkey","ksads_18_903_p","ksads_14_853_p")]



####INJURY CAUSES####
#"tbi_2": MVA
#"tbi_3": Fall or hit by object (accidentally)
#"tbi_4": Fight, hit intentionally, or shaken
#"tbi_5": Explosion or blast

#causes:
# 1=MVA
# 2=Fall or hit by object (accidentally)
# 3=Fight, hit intentionally, or shaken
# 4=Explosion or blast
# 5=Other
# 6=Multiple

TBI_causes<-TBI_info_raw_study[,c("subjectkey","tbi_2","tbi_3","tbi_4","tbi_5")]
TBI_causes$causes<-NA
for (i in 1:dim(TBI_causes)[1]){
  if (sum(as.numeric(TBI_causes[i,colnames(TBI_causes)!=c("subjectkey","causes")])-2)==1){
    #ugly hack because as.numeric keeps turning my factors into numbers that are +2 away from the real number
    TBI_causes[i,c("causes")]<-ifelse(TBI_causes[i,c("tbi_2")]==1,yes = 1,
                                      no = ifelse(TBI_causes[i,c("tbi_3")]==1,yes=2, no=ifelse(TBI_causes[i,c("tbi_4")]==1,yes=3,
                                                                                               no=ifelse(TBI_causes[i,c("tbi_5")]==1,yes=4,no=0))))
  } else if (sum(as.numeric(TBI_causes[i,colnames(TBI_causes)!=c("subjectkey","causes")])-2)==0){
    TBI_causes[i,c("causes")]<-5
  }else{
    TBI_causes[i,c("causes")]<-6
  }
  
}


TBI_causes$totalTBIs<-apply(TBI_causes[,colnames(TBI_causes)!=c("subjectkey","causes")],1,function(x){sum(as.numeric(x))})
TBI_causes$causes<-as.factor(TBI_causes$causes)
TBI_causes<-TBI_causes[,c("subjectkey","causes","totalTBIs")]


####SOCIODEMOGRAPHIC CHARACTERISTICS: PARENT DEMOGRAPHICS SURVEY####
sdm<-read.table('/Users/Guido/Desktop/GeneralData/ABCD/GeneralData/pdem02.txt',sep='\t',header=T)
# ## ## ## ## ## ## ## ## NOTE: THIS FILE CONTAINS# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
#           demo_comb_income_v2:	What is your TOTAL COMBINED FAMILY INCOME for the past 12 months?                   ##
#           demo_prtnr_ed_v2:	What is the highest grade or level of school your partner completed                     ##
#           demo_prnt_ed_v2:	What is the highest grade or level of school you have completed                         ##
#           demo_race_a_p___0:	What race do you consider the child to be? Please check all that apply                ##
#           FROM THE LATTER 3 WE BUILD:                                                                               ##
#           highest_perental_edu: highest parental education (highest edu between two parents)                        ##
#           race-ethnicity: race ethnicity                                                                            ## 
#             According to Eman:                                                                                      ## 
#                Asian Indian, Chinese, Fillipino, Hapanese, Korean, Vietnamese, Other Asian =Asian                   ## 
#                American Indian/Native American, Alaska Native = AIAN                                                ## 
#                Native Hawaiian, Guamanian, Samoan, Other Pacific Islander = NHPI                                    ## 
#                All others are left non-modified                                                                     ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##  ## 
datadic_sdm <- sdm[1,]
datadic_sdm <-t(datadic_sdm)
sdm <-sdm[-1,]
datadic_sdm<-as.data.frame(datadic_sdm)

sdm_clean<-sdm[sdm$eventname=="baseline_year_1_arm_1",]
#Construct a race/ethnicity var
#NOTE: Code for race ethnicity was provided by Eman Nishat
# create dummys first 
sdm_clean$white <- 0 
sdm_clean$white[sdm_clean$demo_race_a_p___10==1] <- 1
#check - values should be same
table(sdm_clean$white) 
table(sdm_clean$demo_race_a_p___10)

# black/ african american 
sdm_clean$blk <- 0 
sdm_clean$blk [sdm_clean$demo_race_a_p___11==1] <- 1
table(sdm_clean$blk) 

# asian 
#asian indian, chinese, fillipino, japanese, korean, vietnamese, other asian
sdm_clean$asian <- 0 
sdm_clean$asian [sdm_clean$demo_race_a_p___18==1 | sdm_clean$demo_race_a_p___19==1 | sdm_clean$demo_race_a_p___20==1 | sdm_clean$demo_race_a_p___21==1 |
                   sdm_clean$demo_race_a_p___22==1 | sdm_clean$demo_race_a_p___23==1 | sdm_clean$demo_race_a_p___24==1  ] <- 1
table(sdm_clean$asian) 

# aian 
#american indian/native american, alaska native 
sdm_clean$aian <- 0 
sdm_clean$aian [sdm_clean$demo_race_a_p___12==1 | sdm_clean$demo_race_a_p___13==1] <- 1
table(sdm_clean$aian) 

# nhpi 
#native hawaiian, guamanian, samoan, other pacific islander
sdm_clean$nhpi <- 0 
sdm_clean$nhpi [sdm_clean$demo_race_a_p___14==1 | sdm_clean$demo_race_a_p___15==1 | sdm_clean$demo_race_a_p___16==1 | sdm_clean$demo_race_a_p___17==1] <- 1
table(sdm_clean$nhpi) 

# other 
sdm_clean$oth <- 0 
sdm_clean$oth [sdm_clean$demo_race_a_p___25==1 ] <-1 
table(sdm_clean$oth) 



# race counts with certain groups set to 1 only (asian aian nhpi)
sdm_clean$racecount2 <- 0
sdm_clean$racecount2 <- (sdm_clean$white + sdm_clean$blk + sdm_clean$asian + sdm_clean$aian + sdm_clean$nhpi+ sdm_clean$oth)
table(sdm_clean$racecount2) 

## indicator of mixed race
#create column
sdm_clean$mixed1 <- NA
#not mixed
sdm_clean$mixed1 [sdm_clean$racecount2 <=1]  <- 0 
#anyone with value more than 1 is mixed race
sdm_clean$mixed1 [sdm_clean$racecount2 > 1] <- 1 
table (sdm_clean$mixed1) 

sdm_clean$race_eth <- NA
sdm_clean$race_eth [sdm_clean$demo_race_a_p___10==1 ] <- 1
sdm_clean$race_eth [sdm_clean$demo_race_a_p___11==1 ] <- 2
sdm_clean$race_eth [sdm_clean$demo_race_a_p___12==1 | sdm_clean$demo_race_a_p___13==1 ] <- 5
sdm_clean$race_eth [sdm_clean$demo_race_a_p___14==1 | sdm_clean$demo_race_a_p___15==1 | sdm_clean$demo_race_a_p___16==1 | sdm_clean$demo_race_a_p___17==1 ] <- 6
sdm_clean$race_eth [sdm_clean$demo_race_a_p___18==1 | sdm_clean$demo_race_a_p___19==1 | sdm_clean$demo_race_a_p___20==1 | sdm_clean$demo_race_a_p___21==1 | 
                      sdm_clean$demo_race_a_p___22==1 | sdm_clean$demo_race_a_p___23==1 | sdm_clean$demo_race_a_p___24==1 ] <- 4
sdm_clean$race_eth [sdm_clean$demo_race_a_p___25==1 ] <- 8
sdm_clean$race_eth [sdm_clean$demo_race_a_p___99==1 | sdm_clean$demo_race_a_p___77==1] <- NA
sdm_clean$race_eth [sdm_clean$mixed1==1] <- 9
sdm_clean$race_eth [sdm_clean$demo_ethn_v2==1] <- 3 

# table(sdm_clean$race_ethc,exclude=NULL)

sdm_clean$race_ethc <- factor(sdm_clean$race_eth, levels = c(1,2,3,4,5,6,8,9), labels = c("Non-Hispanic White", "Non-Hispanic Black", "Hispanic", "Asian" ,"AIAN", "NHPI", "OTHER" ,"MULTIPLE"))

#Construct highest_parental_edu#
#Make hpe numeric to facilitate comparison of size
sdm_clean$demo_prnt_ed_v2<-as.numeric(as.character(sdm_clean$demo_prnt_ed_v2))
sdm_clean$demo_prtnr_ed_v2<-as.numeric(as.character(sdm_clean$demo_prtnr_ed_v2))
#Replacing missing values with 0, this way it's the lowest value possible, which allows us to compare the two.
sdm_clean$demo_prnt_ed_v2[is.na(sdm_clean$demo_prnt_ed_v2)==TRUE]<-0
sdm_clean[sdm_clean$demo_prnt_ed_v2==777|sdm_clean$demo_prnt_ed_v2==999,c("demo_prnt_ed_v2")]<-0

sdm_clean$demo_prtnr_ed_v2[is.na(sdm_clean$demo_prtnr_ed_v2)==TRUE]<-0
sdm_clean[sdm_clean$demo_prtnr_ed_v2==777|sdm_clean$demo_prtnr_ed_v2==999,c("demo_prtnr_ed_v2")]<-0
sdm_clean$highest_parental_edu<- ifelse(as.numeric(as.character(sdm_clean$demo_prnt_ed_v2)) >= as.numeric(as.character(sdm_clean$demo_prtnr_ed_v2)), 
                                        as.numeric(as.character(sdm_clean$demo_prnt_ed_v2)), as.numeric(as.character(sdm_clean$demo_prtnr_ed_v2)))
sdm_clean<-sdm_clean[,c("subjectkey","demo_comb_income_v2","race_ethc","highest_parental_edu")]

#Change the 999 in income to NA
sdm_clean<-sdm_clean[,c("subjectkey","demo_comb_income_v2","race_ethc","highest_parental_edu")]
sdm_clean[sdm_clean$demo_comb_income_v2==999|sdm_clean$demo_comb_income_v2==777,c("demo_comb_income_v2")]<-NA


####INCORPORATE "EXTERNAL" VARIABLES####
external_vars<-merge(TBI_causes,ksads_numeric,all.x=TRUE)#Incorporate TBI causes and adverse psych outcome
external_vars<-merge(external_vars,sdm_clean)#Incorporate sdm
external_vars<-merge(external_vars,nuisance_covs[,-which(colnames(nuisance_covs)%in%c("interview_age"))],by="subjectkey")


#SPLIT EXTERNAL VARS INTO TRAIN AND TEST#
external_vars$set<-ifelse(external_vars$subjectkey%in%mydata_t09$yy_train$subjectkey,1,2)

#Fill empty cells with NAs
tmp_ext<-external_vars[rowSums(external_vars=="",na.rm = T)>0,colSums(external_vars=="",na.rm = T)>0] #Store the columns that contain a missing value
for (i in 1:length(colnames(tmp_ext))){ #loop through the number of cols with empty values
  external_vars[rownames(tmp_ext)[tmp_ext[,colnames(tmp_ext)[i]]==""],colnames(tmp_ext)[i]]<-NA #from external_vars, take the rownames corresponding to the rows in tmp_ext that have an empty value
  
}


####IMPUTE MISSING EXTERNAL VARS####
message("Imputing missing external variables")
#This script should be imputing within each set
for (j in 1:max(external_vars$set)){#loop through the set
  tmp_external<-external_vars[external_vars$set==j,]#temporarily store data of one set
  for (i in 1:dim(tmp_external)[2]){
    idx_miss<-rownames(tmp_external[is.na(tmp_external[,i])==1,])#works because row names are preserved
    if (length(idx_miss)!=0){
      imp_size=dim(as.data.frame(idx_miss))[1]
      # set.seed(123) # Should we set a seed?
      idx_comp<-rownames(tmp_external[is.na(tmp_external[,i])==0,])
      set.seed(i)
      idx_rep<-sample(idx_comp,size=imp_size)
      external_vars[idx_miss,i]<-external_vars[idx_rep,i]
    }
  }
}

####Analyses on Replication Set####
####1. To show replicability, we want to redo analyses on testing set and compare the results. This tests robustness of pipeline to different inputs####
#This question requires redoing commit on these subjects, and redoing PCA, PLSc, etc.
#Note: data from these analyses only serves for a small number of comparisons
#1.1. How do the PCA results compare?

###Apply PCA on testing set###

colnames(mydata_t09$xx_test$fa)<-substr(colnames(mydata_t09$xx_test$fa),regexpr("mean",colnames(mydata_t09$xx_test$fa))+5,regexpr("fa",colnames(mydata_t09$xx_test$fa))-2)
colnames(mydata_t09$xx_test$md)<-substr(colnames(mydata_t09$xx_test$md),regexpr("mean",colnames(mydata_t09$xx_test$md))+5,regexpr("md",colnames(mydata_t09$xx_test$md))-2)
colnames(mydata_t09$xx_test$rd)<-substr(colnames(mydata_t09$xx_test$rd),regexpr("mean",colnames(mydata_t09$xx_test$rd))+5,regexpr("rd",colnames(mydata_t09$xx_test$rd))-2)
colnames(mydata_t09$xx_test$ad)<-substr(colnames(mydata_t09$xx_test$ad),regexpr("mean",colnames(mydata_t09$xx_test$ad))+5,regexpr("ad",colnames(mydata_t09$xx_test$ad))-2)
colnames(mydata_t09$xx_test$nufo)<-substr(colnames(mydata_t09$xx_test$nufo),regexpr("mean",colnames(mydata_t09$xx_test$nufo))+5,regexpr("nufo",colnames(mydata_t09$xx_test$nufo))-2)
colnames(mydata_t09$xx_test$afdf)<-substr(colnames(mydata_t09$xx_test$afdf),regexpr("mean",colnames(mydata_t09$xx_test$afdf))+5,regexpr("afdf",colnames(mydata_t09$xx_test$afdf))-2)

#HERE I NEED TO FIND OUT IF ANY METRIC HAD LESS CONNS, AND ONLY KEEP THE SMALLEST NUMBER OF CONNS
tmp_dims<-c(dim(mydata_t09$xx_test$fa)[2],dim(mydata_t09$xx_test$md)[2],dim(mydata_t09$xx_test$rd)[2],dim(mydata_t09$xx_test$ad)[2],dim(mydata_t09$xx_test$nufo)[2],dim(mydata_t09$xx_test$afdf)[2])
tmp_mtrs<-c("fa","md","rd","ad","nufo","afdf")
tmp_dim_mtrs<-rbind(tmp_dims,tmp_mtrs)
min_mtr<-tmp_dim_mtrs[rownames(tmp_dim_mtrs)=="tmp_mtrs",which(tmp_dim_mtrs[1,]==min(tmp_dim_mtrs[1,]))][1]
max_mtr<-tmp_dim_mtrs[rownames(tmp_dim_mtrs)=="tmp_mtrs",which(tmp_dim_mtrs[1,]==max(tmp_dim_mtrs[1,]))][1]
tmpbigcols<-substr(colnames(as.data.frame(mydata_t09$xx_test[max_mtr])),regexpr("\\.",colnames(as.data.frame(mydata_t09$xx_test[max_mtr])))+1,20)

tmp_misscols<-cbind(as.vector(tmpbigcols[!tmpbigcols%in%colnames(mydata_t09$xx_test$fa)]),as.vector(tmpbigcols[!tmpbigcols%in%colnames(mydata_t09$xx_test$md)]),
                    as.vector(tmpbigcols[!tmpbigcols%in%colnames(mydata_t09$xx_test$rd)]),as.vector(tmpbigcols[!tmpbigcols%in%colnames(mydata_t09$xx_test$ad)]),
                    as.vector(tmpbigcols[!tmpbigcols%in%colnames(mydata_t09$xx_test$nufo)]),as.vector(tmpbigcols[!tmpbigcols%in%colnames(mydata_t09$xx_test$afdf)]))
misscols<-unique(melt(tmp_misscols)$value)
tmpsmallcols<-tmpbigcols[!tmpbigcols%in%misscols]

# tmpsmallcols<-substr(colnames(as.data.frame(mydata_t09$xx_test[min_mtr])),regexpr("\\.",colnames(as.data.frame(mydata_t09$xx_test[min_mtr])))+1,20)

tmp_submat<-matrix(data=NA,ncol = length(tmpsmallcols),nrow = length(mydata_t09$yy_test$subjectkey))
colnames(tmp_submat)<-tmpsmallcols
tmp_submat<-as.data.frame(tmp_submat)
tmp_submat$subjectkey<-mydata_t09$yy_test$subjectkey
tmp_submat_melt<-melt(tmp_submat,id.vars = "subjectkey")

#HERE I NEED TO ADAPT THE SCRIPT SO THAT IT CAN HANDLE THE FACT THAT DEPENDING ON SOME THRESHOLDS, THERE WON"T BE A MISMATCH BETWEEN NUFO AND OTHER METRICS
if (length(tmpbigcols)>length(tmpsmallcols)){
  #Note that here I'm removing those "extra" conns that metrics other than NuFO had.
  # 
  mydata_t09$xx_test$all_metrics<-melt(mydata_t09$xx_test$fa[,which(colnames(mydata_t09$xx_test$fa)%in%tmpsmallcols)])
  mydata_t09$xx_test$all_metrics$md<-melt(mydata_t09$xx_test$md[,which(colnames(mydata_t09$xx_test$md)%in%tmpsmallcols)])$value
  mydata_t09$xx_test$all_metrics$rd<-melt(mydata_t09$xx_test$rd[,which(colnames(mydata_t09$xx_test$rd)%in%tmpsmallcols)])$value
  mydata_t09$xx_test$all_metrics$ad<-melt(mydata_t09$xx_test$ad[,which(colnames(mydata_t09$xx_test$ad)%in%tmpsmallcols)])$value
  mydata_t09$xx_test$all_metrics$nufo<-melt(mydata_t09$xx_test$nufo[,which(colnames(mydata_t09$xx_test$nufo)%in%tmpsmallcols)])$value
  mydata_t09$xx_test$all_metrics$afdf<-melt(mydata_t09$xx_test$afdf[,which(colnames(mydata_t09$xx_test$afdf)%in%tmpsmallcols)])$value
  mydata_t09$xx_test$all_metrics$subjectkey<-tmp_submat_melt$subjectkey
  colnames(mydata_t09$xx_test$all_metrics)<-c("conn","fa","md","rd","ad","nufo","afdf","subjectkey")
  
} else if (length(tmpbigcols)==length(tmpsmallcols)){
  #Note that here I'm removing those "extra" conns that metrics other than NuFO had.
  mydata_t09$xx_test$all_metrics<-melt(mydata_t09$xx_test$fa)
  mydata_t09$xx_test$all_metrics$md<-melt(mydata_t09$xx_test$md)$value
  mydata_t09$xx_test$all_metrics$rd<-melt(mydata_t09$xx_test$rd)$value
  mydata_t09$xx_test$all_metrics$ad<-melt(mydata_t09$xx_test$ad)$value
  mydata_t09$xx_test$all_metrics$nufo<-melt(mydata_t09$xx_test$nufo)$value
  mydata_t09$xx_test$all_metrics$afdf<-melt(mydata_t09$xx_test$afdf)$value
  mydata_t09$xx_test$all_metrics$subjectkey<-tmp_submat_melt$subjectkey
  colnames(mydata_t09$xx_test$all_metrics)<-c("conn","fa","md","rd","ad","nufo","afdf","subjectkey")
}


pca.across.metrics.test<-prcomp(mydata_t09$xx_test$all_metrics[mydata_t09$xx_test$all_metrics$conn%in%selected_conns[selected_conns$c1w_pos==1,c("selected_conns")],
                                                               -which(colnames(mydata_t09$xx_test$all_metrics)%in%c("subjectkey","conn"))],center = TRUE,scale. = TRUE)

corrplot(pca.across.metrics.test$rotation, is.corr=FALSE)
summary(pca.across.metrics.test)
# 0.7847 0.1712

#Comparison of loadings
cor.test(pca.across.metrics$rotation[,1],pca.across.metrics.test$rotation[,1])
# r= -0.9997, p<0.001

cor.test(pca.across.metrics$rotation[,2],pca.across.metrics.test$rotation[,2])
# r= -0.996, p<0.001

#1.2. How many of the same connections were selected by the uvarfeatsel?
#Note: to do this I need to already run the PLSc, since the uvarfeatsel is implemented inside

###EXTRACT FIRST 2 PCs###
#Extract first 2 PCS, they account for 95% of variance
PC1.scores.test<-as.data.frame(pca.across.metrics.test$x[,1])
colnames(PC1.scores.test)<-"score"
PC2.scores.test<-as.data.frame(pca.across.metrics.test$x[,2])
colnames(PC2.scores.test)<-"score"
PC1.scores.test$conns<-mydata_t09$xx_test$all_metrics[mydata_t09$xx_test$all_metrics$conn%in%selected_conns[selected_conns$c1w_pos==1,c("selected_conns")],c("conn")]
PC1.scores.test$subjectkey<-mydata_t09$xx_test$all_metrics[mydata_t09$xx_test$all_metrics$conn%in%selected_conns[selected_conns$c1w_pos==1,c("selected_conns")],c("subjectkey")]
PC2.scores.test$conns<-mydata_t09$xx_test$all_metrics[mydata_t09$xx_test$all_metrics$conn%in%selected_conns[selected_conns$c1w_pos==1,c("selected_conns")],c("conn")]
PC2.scores.test$subjectkey<-mydata_t09$xx_test$all_metrics[mydata_t09$xx_test$all_metrics$conn%in%selected_conns[selected_conns$c1w_pos==1,c("selected_conns")],c("subjectkey")]


###RESTRUCTURE PCS BACK TO WIDE FORMAT###
PC1.scores.test.wide<-reshape(PC1.scores.test,idvar = "subjectkey",timevar = "conns",direction = "wide")
colnames(PC1.scores.test.wide)<-substr(colnames(PC1.scores.test.wide),regexpr("\\.",colnames(PC1.scores.test.wide))+1,25)

PC2.scores.test.wide<-reshape(PC2.scores.test,idvar = "subjectkey",timevar = "conns",direction = "wide")
colnames(PC2.scores.test.wide)<-substr(colnames(PC2.scores.test.wide),regexpr("\\.",colnames(PC2.scores.test.wide))+1,25)

##PLSc##
plsc.t09.pc1.test<-plsc_and_permute(PC1.scores.test.wide[,-which(colnames(PC1.scores.test.wide)%in%c("subjectkey"))],mydata_t09$yy_test[,-which(colnames(mydata_t09$yy_test)%in%c("subjectkey","sex"))],
                                    num_comps = 20,nperms = 2000,topn = 200,out_flag = 1)

plsc.t09.pc2.test<-plsc_and_permute(PC2.scores.test.wide[,-which(colnames(PC2.scores.test.wide)%in%c("subjectkey"))],mydata_t09$yy_test[,-which(colnames(mydata_t09$yy_test)%in%c("subjectkey","sex"))],
                                    num_comps = 20,nperms = 2000,topn = 200,out_flag = 1)

#Now check how many connections were also selected in this dataset
sum(rownames(plsc.t09.pc1[[4]]$TExPosition.Data$pdq$p)%in%rownames(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$p))#*73* of the connections from discovery set were selected by uvarfeatsel in replication set. 

#1.3. Comparisons of loadings (correlations + polar plots)
common_cons<-rownames(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$p)[rownames(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$p)%in%rownames(plsc.t09.pc1[[4]]$TExPosition.Data$pdq$p)]

#cor of the loadings OF THE COMMON CONNECTIONS
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$p[rownames(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$p)%in%common_cons,1],plsc.t09.pc1[[4]]$TExPosition.Data$pdq$p[rownames(plsc.t09.pc1[[4]]$TExPosition.Data$pdq$p)%in%common_cons,1])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$p[rownames(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$p)%in%common_cons,2],plsc.t09.pc1[[4]]$TExPosition.Data$pdq$p[rownames(plsc.t09.pc1[[4]]$TExPosition.Data$pdq$p)%in%common_cons,2])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$p[rownames(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$p)%in%common_cons,3],plsc.t09.pc1[[4]]$TExPosition.Data$pdq$p[rownames(plsc.t09.pc1[[4]]$TExPosition.Data$pdq$p)%in%common_cons,3])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$p[rownames(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$p)%in%common_cons,4],plsc.t09.pc1[[4]]$TExPosition.Data$pdq$p[rownames(plsc.t09.pc1[[4]]$TExPosition.Data$pdq$p)%in%common_cons,4])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$p[rownames(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$p)%in%common_cons,5],plsc.t09.pc1[[4]]$TExPosition.Data$pdq$p[rownames(plsc.t09.pc1[[4]]$TExPosition.Data$pdq$p)%in%common_cons,5])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$p[rownames(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$p)%in%common_cons,6],plsc.t09.pc1[[4]]$TExPosition.Data$pdq$p[rownames(plsc.t09.pc1[[4]]$TExPosition.Data$pdq$p)%in%common_cons,6])#r=0.006, p=0.957 (min)
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$p[rownames(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$p)%in%common_cons,7],plsc.t09.pc1[[4]]$TExPosition.Data$pdq$p[rownames(plsc.t09.pc1[[4]]$TExPosition.Data$pdq$p)%in%common_cons,7])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$p[rownames(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$p)%in%common_cons,8],plsc.t09.pc1[[4]]$TExPosition.Data$pdq$p[rownames(plsc.t09.pc1[[4]]$TExPosition.Data$pdq$p)%in%common_cons,8])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$p[rownames(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$p)%in%common_cons,9],plsc.t09.pc1[[4]]$TExPosition.Data$pdq$p[rownames(plsc.t09.pc1[[4]]$TExPosition.Data$pdq$p)%in%common_cons,9])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$p[rownames(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$p)%in%common_cons,10],plsc.t09.pc1[[4]]$TExPosition.Data$pdq$p[rownames(plsc.t09.pc1[[4]]$TExPosition.Data$pdq$p)%in%common_cons,10])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$p[rownames(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$p)%in%common_cons,11],plsc.t09.pc1[[4]]$TExPosition.Data$pdq$p[rownames(plsc.t09.pc1[[4]]$TExPosition.Data$pdq$p)%in%common_cons,11])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$p[rownames(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$p)%in%common_cons,12],plsc.t09.pc1[[4]]$TExPosition.Data$pdq$p[rownames(plsc.t09.pc1[[4]]$TExPosition.Data$pdq$p)%in%common_cons,12])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$p[rownames(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$p)%in%common_cons,13],plsc.t09.pc1[[4]]$TExPosition.Data$pdq$p[rownames(plsc.t09.pc1[[4]]$TExPosition.Data$pdq$p)%in%common_cons,13])#r=0.227, p=0.053 (max)
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$p[rownames(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$p)%in%common_cons,14],plsc.t09.pc1[[4]]$TExPosition.Data$pdq$p[rownames(plsc.t09.pc1[[4]]$TExPosition.Data$pdq$p)%in%common_cons,14])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$p[rownames(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$p)%in%common_cons,15],plsc.t09.pc1[[4]]$TExPosition.Data$pdq$p[rownames(plsc.t09.pc1[[4]]$TExPosition.Data$pdq$p)%in%common_cons,15])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$p[rownames(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$p)%in%common_cons,16],plsc.t09.pc1[[4]]$TExPosition.Data$pdq$p[rownames(plsc.t09.pc1[[4]]$TExPosition.Data$pdq$p)%in%common_cons,16])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$p[rownames(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$p)%in%common_cons,17],plsc.t09.pc1[[4]]$TExPosition.Data$pdq$p[rownames(plsc.t09.pc1[[4]]$TExPosition.Data$pdq$p)%in%common_cons,17])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$p[rownames(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$p)%in%common_cons,18],plsc.t09.pc1[[4]]$TExPosition.Data$pdq$p[rownames(plsc.t09.pc1[[4]]$TExPosition.Data$pdq$p)%in%common_cons,18])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$p[rownames(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$p)%in%common_cons,19],plsc.t09.pc1[[4]]$TExPosition.Data$pdq$p[rownames(plsc.t09.pc1[[4]]$TExPosition.Data$pdq$p)%in%common_cons,19])


#cor of the loadings OF THE symptoms
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$q[,1],plsc.t09.pc1[[4]]$TExPosition.Data$pdq$q[,1])#r=0.810, p<0.001
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$q[,2],plsc.t09.pc1[[4]]$TExPosition.Data$pdq$q[,2])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$q[,3],plsc.t09.pc1[[4]]$TExPosition.Data$pdq$q[,3])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$q[,4],plsc.t09.pc1[[4]]$TExPosition.Data$pdq$q[,4])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$q[,5],plsc.t09.pc1[[4]]$TExPosition.Data$pdq$q[,5])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$q[,6],plsc.t09.pc1[[4]]$TExPosition.Data$pdq$q[,6])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$q[,7],plsc.t09.pc1[[4]]$TExPosition.Data$pdq$q[,7])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$q[,8],plsc.t09.pc1[[4]]$TExPosition.Data$pdq$q[,8])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$q[,9],plsc.t09.pc1[[4]]$TExPosition.Data$pdq$q[,9])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$q[,10],plsc.t09.pc1[[4]]$TExPosition.Data$pdq$q[,10])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$q[,11],plsc.t09.pc1[[4]]$TExPosition.Data$pdq$q[,11])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$q[,12],plsc.t09.pc1[[4]]$TExPosition.Data$pdq$q[,12])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$q[,13],plsc.t09.pc1[[4]]$TExPosition.Data$pdq$q[,13])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$q[,14],plsc.t09.pc1[[4]]$TExPosition.Data$pdq$q[,14])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$q[,15],plsc.t09.pc1[[4]]$TExPosition.Data$pdq$q[,15])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$q[,16],plsc.t09.pc1[[4]]$TExPosition.Data$pdq$q[,16])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$q[,17],plsc.t09.pc1[[4]]$TExPosition.Data$pdq$q[,17])#r=-0.087, p=0.722 (min)
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$q[,18],plsc.t09.pc1[[4]]$TExPosition.Data$pdq$q[,18])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$q[,19],plsc.t09.pc1[[4]]$TExPosition.Data$pdq$q[,19])

####BootstrapTests####
boot.pc1.t09.test<-Boot4PLSC(PC1.scores.test.wide[,-which(colnames(PC1.scores.test.wide)%in%c("subjectkey"))],mydata_t09$yy_test[,-which(colnames(mydata_t09$yy_test)%in%c("subjectkey","sex"))],
                         center1 = TRUE,center2 = TRUE,scale1 = TRUE,scale2 = TRUE,nf2keep = 19,nIter = 2000,critical.value = 1.96,eig = TRUE)
boot.pc1.t09.test$bootRatios.j[, plsc.t09.pc1.test[[1]]<.05]
boot.pc1.t09.test$bootRatiosSignificant.j[, plsc.t09.pc1.test[[1]]<.05]


####Illustrate all sig LFs####
plot_data<-plsc.t09.pc1.test
boot_data<-boot.pc1.t09.test

###APPLIES TO ALL POLAR PLOTS###
#Store loadings
df<-as.data.frame(plot_data[[4]]$TExPosition.Data$pdq$q)
#Change column names
colnames(df)<-paste("Y",substr(colnames(df),2,5),sep="_")
#Only keep the lfs that were significant
df<-df[,plot_data[[1]]<.05]
#Store symptoms as a new column
df$symptom<-rownames(df)
#Create a column to specify symptom domains (not really used atm)
df$domain<-c("mood","mood","cognitive","mood","somatic","somatic","somatic","somatic","somatic","somatic","somatic","somatic","somatic","cognitive","cognitive","cognitive","cognitive","cognitive","cognitive")
#Create a column for the common name of the symptom
df$gen_prob<-c("depression","anxiety","attention","aggression","headache","nausea","vomiting","dizziness","fatigue","sleep+","sleep-","drowsiness","troublesleeping","s-recall","l-recall","seq-memory","work-memory","proc-speed","cardsorting")
#Flag significant symptoms based on bootstrap tests
df$boottest<-boot_data$bootRatiosSignificant.j[,plot_data[[1]]<.05]
#Melt df
df_melt<-melt(df)

#Create color palette
c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "sienna", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)


#Create df for labels
label_data <- as.data.frame(sort(df$gen_prob))
colnames(label_data)<-"gen_prob"
label_data$id<-seq(1,19)
# calculate the ANGLE of the labels
number_of_bar <- nrow(label_data)
angle <-  90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)

###@@@###

#Note: will reverse a few vars for visualization purposes
df_melt[df_melt$variable=="Y_1",c("value")]<-df_melt[df_melt$variable=="Y_1",c("value")]*-1
df_melt[df_melt$variable=="Y_2",c("value")]<-df_melt[df_melt$variable=="Y_2",c("value")]*-1
df_melt[df_melt$variable=="Y_3",c("value")]<-df_melt[df_melt$variable=="Y_3",c("value")]*-1
df_melt[df_melt$variable=="Y_7",c("value")]<-df_melt[df_melt$variable=="Y_7",c("value")]*-1

###SPECIFIC POLAR PLOTS###
#Note: the commands that are commented out correspond to the latent factors that were not found to be significant

# p1<-create_polar_plots(df_melt,var = "Y_1",color = "black",colorpal = c25,ylims = c(-0.5,0.8),boot=TRUE,label_data)
# p2<-create_polar_plots(df_melt,var = "Y_2",color = "black",colorpal = c25,ylims = c(-0.9,0.9),boot=TRUE,label_data)
# p3<-create_polar_plots(df_melt,var = "Y_3",color = "black",colorpal = c25,ylims = c(-0.6,0.8),boot=TRUE,label_data)
# p4<-create_polar_plots(df_melt,var = "Y_4",color = "black",colorpal = c25,ylims = c(-0.7,0.8),boot=TRUE,label_data)
p5<-create_polar_plots(df_melt,var = "Y_5",color = "black",colorpal = c25,ylims = c(-0.6,0.8),boot=TRUE,label_data)
p6<-create_polar_plots(df_melt,var = "Y_6",color = "black",colorpal = c25,ylims = c(-0.7,0.8),boot=TRUE,label_data)
p7<-create_polar_plots(df_melt,var = "Y_7",color = "black",colorpal = c25,ylims = c(-0.8,0.8),boot=TRUE,label_data)
p8<-create_polar_plots(df_melt,var = "Y_8",color = "black",colorpal = c25,ylims = c(-0.6,0.8),boot=TRUE,label_data)
p9<-create_polar_plots(df_melt,var = "Y_9",color = "black",colorpal = c25,ylims = c(-0.8,0.8),boot=TRUE,label_data)
p10<-create_polar_plots(df_melt,var = "Y_10",color = "black",colorpal = c25,ylims = c(-0.7,0.8),boot=TRUE,label_data)
# p11<-create_polar_plots(df_melt,var = "Y_11",color = "black",colorpal = c25,ylims = c(-0.6,0.8),boot=TRUE,label_data)
p12<-create_polar_plots(df_melt,var = "Y_12",color = "black",colorpal = c25,ylims = c(-0.8,0.8),boot=TRUE,label_data)
# p13<-create_polar_plots(df_melt,var = "Y_13",color = "black",colorpal = c25,ylims = c(-0.8,0.8),boot=TRUE,label_data)
# p14<-create_polar_plots(df_melt,var = "Y_14",color = "black",colorpal = c25,ylims = c(-0.8,0.8),boot=TRUE,label_data)
# p15<-create_polar_plots(df_melt,var = "Y_15",color = "black",colorpal = c25,ylims = c(-0.8,0.8),boot=TRUE,label_data)
# p16<-create_polar_plots(df_melt,var = "Y_16",color = "black",colorpal = c25,ylims = c(-0.8,0.9),boot=TRUE,label_data)
# p17<-create_polar_plots(df_melt,var = "Y_17",color = "black",colorpal = c25,ylims = c(-0.9,0.9),boot=TRUE,label_data)
# p18<-create_polar_plots(df_melt,var = "Y_18",color = "black",colorpal = c25,ylims = c(-0.8,0.8),boot=TRUE,label_data)
# p19<-create_polar_plots(df_melt,var = "Y_19",color = "black",colorpal = c25,ylims = c(-0.8,0.8),boot=TRUE,label_data)

plot_grid(p5,p6,p7,p8,p9,p10,p12)






####2. To test robustness of results, we want to project the replication data onto the latent spaces derived on the discovery set####

test.scaled <- scale(mydata_t09$xx_test$all_metrics[mydata_t09$xx_test$all_metrics$conn%in%selected_conns[selected_conns$c1w_pos==1,c("selected_conns")],-which(colnames(mydata_t09$xx_test$all_metrics)%in%c("subjectkey","conn"))], 
                     center= pca.across.metrics$center,scale=pca.across.metrics$scale)
test.projected <- test.scaled %*% pca.across.metrics$rotation
test.projected<-as.data.frame(test.projected)

#Extract first 2 PCS, they account for 97% of variance
##NOTE: "Ptest" because these are not the scores from the PCA on replication ("test") set, but the scores from the *P*rojected "test" data
PC1.scores.Ptest<-as.data.frame(test.projected$PC1)
colnames(PC1.scores.Ptest)<-"score"
PC2.scores.Ptest<-as.data.frame(test.projected$PC2)
colnames(PC2.scores.Ptest)<-"score"

PC1.scores.Ptest$conns<-mydata_t09$xx_test$all_metrics[mydata_t09$xx_test$all_metrics$conn%in%selected_conns[selected_conns$c1w_pos==1,c("selected_conns")],c("conn")]
PC1.scores.Ptest$subjectkey<-mydata_t09$xx_test$all_metrics[mydata_t09$xx_test$all_metrics$conn%in%selected_conns[selected_conns$c1w_pos==1,c("selected_conns")],c("subjectkey")]
PC2.scores.Ptest$conns<-mydata_t09$xx_test$all_metrics[mydata_t09$xx_test$all_metrics$conn%in%selected_conns[selected_conns$c1w_pos==1,c("selected_conns")],c("conn")]
PC2.scores.Ptest$subjectkey<-mydata_t09$xx_test$all_metrics[mydata_t09$xx_test$all_metrics$conn%in%selected_conns[selected_conns$c1w_pos==1,c("selected_conns")],c("subjectkey")]


#Restructure PCs back to wide format
PC1.scores.Ptest.wide<-reshape(PC1.scores.Ptest,idvar = "subjectkey",timevar = "conns",direction = "wide")
colnames(PC1.scores.Ptest.wide)<-substr(colnames(PC1.scores.Ptest.wide),regexpr("\\.",colnames(PC1.scores.Ptest.wide))+1,25)

PC2.scores.Ptest.wide<-reshape(PC2.scores.Ptest,idvar = "subjectkey",timevar = "conns",direction = "wide")
colnames(PC2.scores.Ptest.wide)<-substr(colnames(PC2.scores.Ptest.wide),regexpr("\\.",colnames(PC2.scores.Ptest.wide))+1,25)


#Scaling and centering to the latent spaces
#PC1
pc1.xscores.Ptest<-as.matrix(PC1.scores.Ptest.wide[,colnames(PC1.scores.Ptest.wide)%in%rownames(plsc.t09.pc1[[4]]$TExPosition.Data$pdq$p)])#We do this because the plsc does uvarfeatsel, so we don't really use all the conns.
scaled.pc1.xscores.Ptest<-scale(pc1.xscores.Ptest, center= plsc.t09.pc1[[4]]$TExPosition.Data$data1.norm$center,scale = plsc.t09.pc1[[4]]$TExPosition.Data$data1.norm$scale)

pc1.yscores.Ptest<-as.matrix(mydata_t09$yy_test[,-which(colnames(mydata_t09$yy_test)%in%c("subjectkey"))])
scaled.pc1.yscores.Ptest<-scale(pc1.yscores.Ptest, center= plsc.t09.pc1[[4]]$TExPosition.Data$data2.norm$center,scale= plsc.t09.pc1[[4]]$TExPosition.Data$data2.norm$scale)

#PC2
pc2.xscores.Ptest<-as.matrix(PC2.scores.Ptest.wide[,colnames(PC2.scores.Ptest.wide)%in%rownames(plsc.t09.pc2[[4]]$TExPosition.Data$pdq$p)])
scaled.pc2.xscores.Ptest<-scale(pc2.xscores.Ptest, center= plsc.t09.pc2[[4]]$TExPosition.Data$data1.norm$center,scale = plsc.t09.pc2[[4]]$TExPosition.Data$data1.norm$scale)

pc2.yscores.Ptest<-as.matrix(mydata_t09$yy_test[,-which(colnames(mydata_t09$yy_test)%in%c("subjectkey"))])
scaled.pc2.yscores.Ptest<-scale(pc2.yscores.Ptest, center= plsc.t09.pc2[[4]]$TExPosition.Data$data2.norm$center,scale= plsc.t09.pc2[[4]]$TExPosition.Data$data2.norm$scale)

#Projecting to latent spaces
#PC1
projected.xscores.Ptest.pc1<-as.data.frame(apply(plsc.t09.pc1[[4]]$TExPosition.Data$pdq$p,2,function(x){scaled.pc1.xscores.Ptest %*% as.vector(x)}))
projected.yscores.Ptest.pc1<-as.data.frame(apply(plsc.t09.pc1[[4]]$TExPosition.Data$pdq$q,2,function(x){scaled.pc1.yscores.Ptest %*% as.vector(x)}))

#PC2
projected.xscores.Ptest.pc2<-as.data.frame(apply(plsc.t09.pc2[[4]]$TExPosition.Data$pdq$p,2,function(x){scaled.pc2.xscores.Ptest %*% as.vector(x)}))
projected.yscores.Ptest.pc2<-as.data.frame(apply(plsc.t09.pc2[[4]]$TExPosition.Data$pdq$q,2,function(x){scaled.pc2.yscores.Ptest %*% as.vector(x)}))


#Select significant LFs
projected.xscores.Ptest.pc1<-as.data.frame(projected.xscores.Ptest.pc1[,plsc.t09.pc1[[1]]<0.05])
projected.yscores.Ptest.pc1<-as.data.frame(projected.yscores.Ptest.pc1[,plsc.t09.pc1[[1]]<0.05])

projected.xscores.Ptest.pc2<-as.data.frame(projected.xscores.Ptest.pc2[,plsc.t09.pc2[[1]]<0.05])
projected.yscores.Ptest.pc2<-as.data.frame(projected.yscores.Ptest.pc2[,plsc.t09.pc2[[1]]<0.05])

#Set column names
colnames(projected.xscores.Ptest.pc1)<-paste("X",which(plsc.t09.pc1[[1]]<0.05),sep="_")
colnames(projected.yscores.Ptest.pc1)<-paste("Y",which(plsc.t09.pc1[[1]]<0.05),sep="_")

colnames(projected.xscores.Ptest.pc2)<-paste("X",which(plsc.t09.pc2[[1]]<0.05),sep="_")
colnames(projected.yscores.Ptest.pc2)<-paste("Y",which(plsc.t09.pc2[[1]]<0.05),sep="_")

#Bind X and Y projections
projected.scores.Ptest.pc1<-cbind(projected.xscores.Ptest.pc1,projected.yscores.Ptest.pc1)

projected.scores.Ptest.pc2<-cbind(projected.xscores.Ptest.pc2,projected.yscores.Ptest.pc2)


#Add subjectkey
projected.scores.Ptest.pc1$subjectkey<-mydata_t09$yy_test$subjectkey

projected.scores.Ptest.pc2$subjectkey<-mydata_t09$yy_test$subjectkey

#Merge projections with external vars
all.data.Ptest.pc1<-merge(projected.scores.Ptest.pc1,external_vars)

all.data.Ptest.pc2<-merge(projected.scores.Ptest.pc2,external_vars)
#Now that you stored all the Replication set data + external vars, you need to build the big dataset with all the data

#Get the replication set outputs
projected.xscores.train.pc1<-as.data.frame(plsc.t09.pc1[[4]]$TExPosition.Data$lx)
projected.yscores.train.pc1<-as.data.frame(plsc.t09.pc1[[4]]$TExPosition.Data$ly)

projected.xscores.train.pc2<-as.data.frame(plsc.t09.pc2[[4]]$TExPosition.Data$lx)
projected.yscores.train.pc2<-as.data.frame(plsc.t09.pc2[[4]]$TExPosition.Data$ly)

#Select significant LFs
projected.xscores.train.pc1<-as.data.frame(projected.xscores.train.pc1[,plsc.t09.pc1[[1]]<0.05])
projected.yscores.train.pc1<-as.data.frame(projected.yscores.train.pc1[,plsc.t09.pc1[[1]]<0.05])

projected.xscores.train.pc2<-as.data.frame(projected.xscores.train.pc2[,plsc.t09.pc2[[1]]<0.05])
projected.yscores.train.pc2<-as.data.frame(projected.yscores.train.pc2[,plsc.t09.pc2[[1]]<0.05])

#Set column names
colnames(projected.xscores.train.pc1)<-paste("X",which(plsc.t09.pc1[[1]]<0.05),sep="_")
colnames(projected.yscores.train.pc1)<-paste("Y",which(plsc.t09.pc1[[1]]<0.05),sep="_")

colnames(projected.xscores.train.pc2)<-paste("X",which(plsc.t09.pc2[[1]]<0.05),sep="_")
colnames(projected.yscores.train.pc2)<-paste("Y",which(plsc.t09.pc2[[1]]<0.05),sep="_")


#Bind X and Y projections
projected.scores.train.pc1<-cbind(projected.xscores.train.pc1,projected.yscores.train.pc1)

projected.scores.train.pc2<-cbind(projected.xscores.train.pc2,projected.yscores.train.pc2)

#Add subjectkey
projected.scores.train.pc1$subjectkey<-mydata_t09$yy_train$subjectkey
projected.scores.train.pc2$subjectkey<-mydata_t09$yy_train$subjectkey

all.data.train.pc1<-merge(projected.scores.train.pc1,external_vars)

all.data.train.pc2<-merge(projected.scores.train.pc2,external_vars)

####MERGE ALL DATA TOGETHER####
all.data.pc1<-rbind(all.data.train.pc1,all.data.Ptest.pc1)#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@YOUR BIG DATA FRAME WITH LFS AND EXTERNAL VARS

all.data.pc2<-rbind(all.data.train.pc2,all.data.Ptest.pc2)#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


#2.1. How much conn/symp covariance is explained by MTMS1 from discovery set (make sure to compare with permutation testing)?
#PC1
#CALCULATE % COVARIANCE EXPLAINED IN REPLICATION SET
Ptest.set.proj.covs.pc1<-diag(cov(projected.xscores.Ptest.pc1,projected.yscores.Ptest.pc1))
pc1.cove<-abs(Ptest.set.proj.covs.pc1)/sum(abs(Ptest.set.proj.covs.pc1))
barplot(abs(Ptest.set.proj.covs.pc1/sum(abs(Ptest.set.proj.covs.pc1))))

#PC2
Ptest.set.proj.covs.pc2<-diag(cov(projected.xscores.Ptest.pc2,projected.yscores.Ptest.pc2))
pc2.cove<-abs(Ptest.set.proj.covs.pc2)/sum(abs(Ptest.set.proj.covs.pc2))


###RUN THIS THROUGH PERMUTATION TESTING###

compute_percCOVe<-function(x,y){
  proj_covs<-diag(cov(x,y))
  abs(proj_covs)/sum(abs(proj_covs))
}

set.seed(123)
# shuffle within scan location
nperms=2000


shuffled_indexes <- sapply(1:nperms, function(x){
  shuffle(1:nrow(projected.xscores.Ptest.pc1))})

null_results <- foreach(i=1:nperms) %dopar% {
  null_model <- compute_percCOVe(projected.xscores.Ptest.pc1[shuffled_indexes[,i],],
                                 projected.yscores.Ptest.pc1)
  
  list(null_model)
  
}

# transform null results lists to data frame
null_dist_cancor <- lapply(null_results, function(x){return(x[[1]])})

null_dist_cancor <- as.data.frame(do.call(rbind, null_dist_cancor))

print("computing p-values")
pvals_cancor <- mapply(function(real, null_dist){
  get_pval(real, null_dist, better="smaller")},
  pc1.cove,
  null_dist_cancor)
#All NS  

length(pc1.cove)
length(null_dist_cancor)

shuffled_indexes <- sapply(1:nperms, function(x){
  shuffle(1:nrow(projected.xscores.Ptest.pc1))})

null_results <- foreach(i=1:nperms) %dopar% {
  null_model <- compute_percCOVe(projected.xscores.Ptest.pc2[shuffled_indexes[,i],],
                                 projected.yscores.Ptest.pc2)
  
  list(null_model)
  
}

# transform null results lists to data frame
null_dist_cancor <- lapply(null_results, function(x){return(x[[1]])})

null_dist_cancor <- as.data.frame(do.call(rbind, null_dist_cancor))

print("computing p-values")
pvals_cancor <- mapply(function(real, null_dist){
  get_pval(real, null_dist, better="smaller")},
  pc2.cove,
  null_dist_cancor)
#All NS



#2.2. Does MTMS expression differ by symptom set? Do this for all MTMS
#TAKE THE PROJECTED SCORES AND COMPUTE CORRELATIONS BETWEEN LF EXPRESSION AND A BINARY VARIABLE FOR SET
#Y
cor.test(all.data.pc1[,c("Y_1")],as.numeric(as.character(all.data.pc1[,c("set")])))
cor.test(all.data.pc1[,c("Y_2")],as.numeric(as.character(all.data.pc1[,c("set")])))
cor.test(all.data.pc1[,c("Y_3")],as.numeric(as.character(all.data.pc1[,c("set")])))
# cor.test(all.data.pc1[,c("Y_4")],as.numeric(as.character(all.data.pc1[,c("set")])))
cor.test(all.data.pc1[,c("Y_5")],as.numeric(as.character(all.data.pc1[,c("set")])))
# cor.test(all.data.pc1[,c("Y_6")],as.numeric(as.character(all.data.pc1[,c("set")])))
cor.test(all.data.pc1[,c("Y_7")],as.numeric(as.character(all.data.pc1[,c("set")])))
cor.test(all.data.pc1[,c("Y_8")],as.numeric(as.character(all.data.pc1[,c("set")])))
cor.test(all.data.pc1[,c("Y_9")],as.numeric(as.character(all.data.pc1[,c("set")])))
cor.test(all.data.pc1[,c("Y_10")],as.numeric(as.character(all.data.pc1[,c("set")])))
cor.test(all.data.pc1[,c("Y_11")],as.numeric(as.character(all.data.pc1[,c("set")])))
cor.test(all.data.pc1[,c("Y_12")],as.numeric(as.character(all.data.pc1[,c("set")])))
cor.test(all.data.pc1[,c("Y_13")],as.numeric(as.character(all.data.pc1[,c("set")])))
cor.test(all.data.pc1[,c("Y_14")],as.numeric(as.character(all.data.pc1[,c("set")])))
cor.test(all.data.pc1[,c("Y_15")],as.numeric(as.character(all.data.pc1[,c("set")])))
cor.test(all.data.pc1[,c("Y_16")],as.numeric(as.character(all.data.pc1[,c("set")])))
cor.test(all.data.pc1[,c("Y_17")],as.numeric(as.character(all.data.pc1[,c("set")])))
cor.test(all.data.pc1[,c("Y_18")],as.numeric(as.character(all.data.pc1[,c("set")])))
# cor.test(all.data.pc1[,c("Y_19")],as.numeric(as.character(all.data.pc1[,c("set")])))
#All NS

#X
cor.test(all.data.pc1[,c("X_1")],as.numeric(as.character(all.data.pc1[,c("set")])))
cor.test(all.data.pc1[,c("X_2")],as.numeric(as.character(all.data.pc1[,c("set")])))
cor.test(all.data.pc1[,c("X_3")],as.numeric(as.character(all.data.pc1[,c("set")])))
# cor.test(all.data.pc1[,c("X_4")],as.numeric(as.character(all.data.pc1[,c("set")])))
cor.test(all.data.pc1[,c("X_5")],as.numeric(as.character(all.data.pc1[,c("set")])))
# cor.test(all.data.pc1[,c("X_6")],as.numeric(as.character(all.data.pc1[,c("set")])))
cor.test(all.data.pc1[,c("X_7")],as.numeric(as.character(all.data.pc1[,c("set")])))
cor.test(all.data.pc1[,c("X_8")],as.numeric(as.character(all.data.pc1[,c("set")])))
cor.test(all.data.pc1[,c("X_9")],as.numeric(as.character(all.data.pc1[,c("set")])))
cor.test(all.data.pc1[,c("X_10")],as.numeric(as.character(all.data.pc1[,c("set")])))
cor.test(all.data.pc1[,c("X_11")],as.numeric(as.character(all.data.pc1[,c("set")])))
cor.test(all.data.pc1[,c("X_12")],as.numeric(as.character(all.data.pc1[,c("set")])))
cor.test(all.data.pc1[,c("X_13")],as.numeric(as.character(all.data.pc1[,c("set")])))
cor.test(all.data.pc1[,c("X_14")],as.numeric(as.character(all.data.pc1[,c("set")])))
cor.test(all.data.pc1[,c("X_15")],as.numeric(as.character(all.data.pc1[,c("set")])))
cor.test(all.data.pc1[,c("X_16")],as.numeric(as.character(all.data.pc1[,c("set")])))
cor.test(all.data.pc1[,c("X_17")],as.numeric(as.character(all.data.pc1[,c("set")])))
cor.test(all.data.pc1[,c("X_18")],as.numeric(as.character(all.data.pc1[,c("set")])))
# cor.test(all.data.pc1[,c("X_19")],as.numeric(as.character(all.data.pc1[,c("set")])))
#All NS

#2.3.: Comparisons of MT scores
#3.3.: Comparisons of MS scores

cors<-cor(plsc.t09.pc1[[4]]$TExPosition.Data$pdq$p[rownames(plsc.t09.pc1[[4]]$TExPosition.Data$pdq$p)%in%rownames(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$p),],
          plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$p[rownames(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$p)%in%rownames(plsc.t09.pc1[[4]]$TExPosition.Data$pdq$p),])
diag(cors)
#0.048085438 -0.074629158 -0.084037238 -0.046562193  0.050440294  0.006446025  0.050195472  0.021330575 -0.097690606 -0.071046031  0.141633263  0.175294018  0.226778135  0.033274578  0.094791551  0.004942270  0.116038584 -0.159708019 -0.170950404
#Low correlations, suggesting that the weights for conns are not all that similar. Take it with a grain of salt though, the connection set is different.

#Symptoms
cors<-cor(plsc.t09.pc1[[4]]$TExPosition.Data$pdq$q,plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$q)
diag(cors)
#Corrs are much better. This follows the polar plots, showing that MS feats are similar
# 0.80997214 -0.29013264 -0.64511601 -0.11465685  0.25482008 -0.20460342 -0.39935171 -0.12547420 -0.18267021 -0.20430239  0.39677681 -0.04353353 -0.28945218 -0.15687225  0.17486670  0.39171682 -0.08737018 -0.17811484  0.60249408

#Compute correlations between projected scores and current scores
#Conns
cors<-cor(projected.xscores.Ptest.pc1,plsc.t09.pc1.test[[4]]$TExPosition.Data$lx)
diag(cors)
# 0.669688512  0.322694089 -0.035847160  0.100370148 -0.439709467 -0.230274688 -0.025631168  0.307768580  0.090118604  0.010773336 -0.331161557  0.004282037  0.031078666 -0.073539445  0.253913499 -0.221513033
#These correlations are much better. The first few LFs have fairly similar projections, and this despite having pretty different connection sets. Note also this pattern whereby the highest corr is for LF1. This 
#is one important support for the idea that LF1 is generally the most consistent LF.

#Symptoms
cors<-cor(projected.yscores.Ptest.pc1,plsc.t09.pc1.test[[4]]$TExPosition.Data$ly)
diag(cors)
# 0.97863031 -0.41394424 -0.69016705 -0.59243769 -0.11797236 -0.40869419 -0.36763883  0.27210343 -0.20640915  0.18358989  0.12886723 -0.23127818  0.26245724  0.42061233 -0.17112541  0.02641009
#The first few correlations are very high, suggesting again that the multi symptom feats are very similar to each other. 



#1.4. Comparisons of expressions
common_cons<-rownames(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$p)[rownames(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$p)%in%rownames(plsc.t09.pc1[[4]]$TExPosition.Data$pdq$p)]

#cor of the expression of the LFs that were significant in the PLSC1
plsc.t09.pc1[[1]]<0.05
#1,2,3,5,7,8,10,11,12,13,14,15,16,17,18

cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$lx[,1],projected.xscores.Ptest.pc1[,1])#.670,p<0.001
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$lx[,2],projected.xscores.Ptest.pc1[,2])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$lx[,3],projected.xscores.Ptest.pc1[,3])#-0.036,0.734
# cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$lx[,4],projected.xscores.Ptest.pc1[,4])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$lx[,5],projected.xscores.Ptest.pc1[,4])
# cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$lx[,6],projected.xscores.Ptest.pc1[,6])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$lx[,7],projected.xscores.Ptest.pc1[,5])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$lx[,8],projected.xscores.Ptest.pc1[,6])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$lx[,9],projected.xscores.Ptest.pc1[,7])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$lx[,10],projected.xscores.Ptest.pc1[,8])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$lx[,11],projected.xscores.Ptest.pc1[,9])#0.036,p=0.737
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$lx[,12],projected.xscores.Ptest.pc1[,10])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$lx[,13],projected.xscores.Ptest.pc1[,11])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$lx[,14],projected.xscores.Ptest.pc1[,12])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$lx[,15],projected.xscores.Ptest.pc1[,13])#0.015,p=0.888
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$lx[,16],projected.xscores.Ptest.pc1[,14])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$lx[,17],projected.xscores.Ptest.pc1[,15])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$lx[,18],projected.xscores.Ptest.pc1[,16])
# cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$lx[,19],projected.xscores.Ptest.pc1[,17])


#cor of the expression OF THE symptoms
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$ly[,1],projected.yscores.Ptest.pc1[,1])#.979,p<0.001
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$ly[,2],projected.yscores.Ptest.pc1[,2])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$ly[,3],projected.yscores.Ptest.pc1[,3])
# cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$lx[,4],projected.xscores.Ptest.pc1[,4])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$ly[,5],projected.yscores.Ptest.pc1[,4])
# cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$lx[,6],projected.xscores.Ptest.pc1[,6])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$ly[,7],projected.yscores.Ptest.pc1[,5])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$ly[,8],projected.yscores.Ptest.pc1[,6])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$ly[,9],projected.yscores.Ptest.pc1[,7])#0.044,0.676
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$ly[,10],projected.yscores.Ptest.pc1[,8])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$ly[,11],projected.yscores.Ptest.pc1[,9])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$ly[,12],projected.yscores.Ptest.pc1[,10])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$ly[,13],projected.yscores.Ptest.pc1[,11])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$ly[,14],projected.yscores.Ptest.pc1[,12])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$ly[,15],projected.yscores.Ptest.pc1[,13])#0.010,p=0.921
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$ly[,16],projected.yscores.Ptest.pc1[,14])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$ly[,17],projected.yscores.Ptest.pc1[,15])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$ly[,18],projected.yscores.Ptest.pc1[,16])
# cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$lx[,19],projected.xscores.Ptest.pc1[,17])





#1.4. Comparisons of expressions
common_cons<-rownames(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$p)[rownames(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$p)%in%rownames(plsc.t09.pc1[[4]]$TExPosition.Data$pdq$p)]

#cor of the loadings OF THE COMMON CONNECTIONS
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$lx[,1],plsc.t09.pc1[[4]]$TExPosition.Data$lx[,1])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$lx[,2],plsc.t09.pc1[[4]]$TExPosition.Data$lx[,2])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$lx[,3],plsc.t09.pc1[[4]]$TExPosition.Data$lx[,3])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$lx[,4],plsc.t09.pc1[[4]]$TExPosition.Data$lx[,4])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$lx[,5],plsc.t09.pc1[[4]]$TExPosition.Data$lx[,5])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$lx[,6],plsc.t09.pc1[[4]]$TExPosition.Data$lx[,6])#r=0.006, p=0.957 (min)
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$lx[,7],plsc.t09.pc1[[4]]$TExPosition.Data$lx[,7])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$lx[,8],plsc.t09.pc1[[4]]$TExPosition.Data$lx[,8])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$lx[,9],plsc.t09.pc1[[4]]$TExPosition.Data$lx[,9])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$lx[,10],plsc.t09.pc1[[4]]$TExPosition.Data$lx[,10])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$lx[,11],plsc.t09.pc1[[4]]$TExPosition.Data$lx[,11])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$lx[,12],plsc.t09.pc1[[4]]$TExPosition.Data$lx[,12])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$lx[,13],plsc.t09.pc1[[4]]$TExPosition.Data$lx[,13])#r=0.227, p=0.053 (max)
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$lx[,14],plsc.t09.pc1[[4]]$TExPosition.Data$lx[,14])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$lx[,15],plsc.t09.pc1[[4]]$TExPosition.Data$lx[,15])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$lx[,16],plsc.t09.pc1[[4]]$TExPosition.Data$lx[,16])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$lx[,17],plsc.t09.pc1[[4]]$TExPosition.Data$lx[,17])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$lx[,18],plsc.t09.pc1[[4]]$TExPosition.Data$lx[,18])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$lx[,19],plsc.t09.pc1[[4]]$TExPosition.Data$lx[,19])


#cor of the loadings OF THE symptoms
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$q[,1],plsc.t09.pc1[[4]]$TExPosition.Data$pdq$q[,1])#r=0.810, p<0.001
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$q[,2],plsc.t09.pc1[[4]]$TExPosition.Data$pdq$q[,2])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$q[,3],plsc.t09.pc1[[4]]$TExPosition.Data$pdq$q[,3])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$q[,4],plsc.t09.pc1[[4]]$TExPosition.Data$pdq$q[,4])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$q[,5],plsc.t09.pc1[[4]]$TExPosition.Data$pdq$q[,5])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$q[,6],plsc.t09.pc1[[4]]$TExPosition.Data$pdq$q[,6])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$q[,7],plsc.t09.pc1[[4]]$TExPosition.Data$pdq$q[,7])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$q[,8],plsc.t09.pc1[[4]]$TExPosition.Data$pdq$q[,8])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$q[,9],plsc.t09.pc1[[4]]$TExPosition.Data$pdq$q[,9])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$q[,10],plsc.t09.pc1[[4]]$TExPosition.Data$pdq$q[,10])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$q[,11],plsc.t09.pc1[[4]]$TExPosition.Data$pdq$q[,11])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$q[,12],plsc.t09.pc1[[4]]$TExPosition.Data$pdq$q[,12])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$q[,13],plsc.t09.pc1[[4]]$TExPosition.Data$pdq$q[,13])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$q[,14],plsc.t09.pc1[[4]]$TExPosition.Data$pdq$q[,14])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$q[,15],plsc.t09.pc1[[4]]$TExPosition.Data$pdq$q[,15])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$q[,16],plsc.t09.pc1[[4]]$TExPosition.Data$pdq$q[,16])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$q[,17],plsc.t09.pc1[[4]]$TExPosition.Data$pdq$q[,17])#r=-0.087, p=0.722 (min)
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$q[,18],plsc.t09.pc1[[4]]$TExPosition.Data$pdq$q[,18])
cor.test(plsc.t09.pc1.test[[4]]$TExPosition.Data$pdq$q[,19],plsc.t09.pc1[[4]]$TExPosition.Data$pdq$q[,19])


####Testing influence of external variables####
####Analyses of sociodemographic variables (Supplementary Figure 2)####
p1.1<-ggplot(data = all.data.pc1[all.data.pc1$set==1,])+
  geom_point(aes(x=X_1*-1,y=Y_1*-1,color=as.numeric(demo_comb_income_v2)),size=2)+
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill="white",color="white"),
        axis.line =element_line(color = "black", size = 1, linetype = "solid"),
        axis.text = element_text(face = "bold", color = "black", 
                                 size = 15),
        axis.title = element_text(face = "bold", color = "black", 
                                  size = 20),
        axis.ticks = element_blank(),
        legend.position="none")+labs(y="PLSc1 - MSF1", x="PLSc1 - MCF1")+scale_color_continuous(type="viridis")
# p1.1

p1.2<-ggplot(data = all.data.pc1[all.data.pc1$set==1,])+
  geom_point(aes(x=X_1*-1,y=Y_1*-1,color=race_ethc),size=2)+
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill="white",color="white"),
        axis.line =element_line(color = "black", size = 1, linetype = "solid"),
        axis.text = element_text(face = "bold", color = "black", 
                                 size = 15),
        axis.title = element_text(face = "bold", color = "black", 
                                  size = 20),
        axis.ticks = element_blank(),
        legend.position="none")+labs(y="PLSc1 - MSF1",x="PLSc1 - MCF1")+scale_color_brewer(palette="Dark2")
# p1.2

p1.3<-ggplot(data = all.data.pc1[all.data.pc1$set==1,])+
  geom_point(aes(x=X_1*-1,y=Y_1*-1,color=sex),size=2)+
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill="white",color="white"),
        axis.line =element_line(color = "black", size = 1, linetype = "solid"),
        axis.text = element_text(face = "bold", color = "black", 
                                 size = 15),
        axis.title = element_text(face = "bold", color = "black", 
                                  size = 20),
        axis.ticks = element_blank(),legend.position="none")+labs(y="PLSc1 - MSF1", x="PLSc1 - MCF1")+scale_color_brewer(palette="Dark2")
# p1.3


#PC2 LF3
p2.1<-ggplot(data = all.data.pc2[all.data.pc1$set==1,])+
  geom_point(aes(x=X_3,y=Y_3,color=as.numeric(demo_comb_income_v2)),size=2)+
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill="white",color="white"),
        axis.line =element_line(color = "black", size = 1, linetype = "solid"),
        axis.text = element_text(face = "bold", color = "black", 
                                 size = 15),
        axis.title = element_text(face = "bold", color = "black", 
                                  size = 20),
        axis.ticks = element_blank(),legend.position="none")+labs(y="PLSc2 - MSF3",x="PLSc2 - MCF3")+scale_color_continuous(type="viridis")
# p2.1

p2.2<-ggplot(data = all.data.pc2[all.data.pc1$set==1,])+
  geom_point(aes(x=X_3,y=Y_3,color=race_ethc),size=2)+
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill="white",color="white"),
        axis.line =element_line(color = "black", size = 1, linetype = "solid"),
        axis.text = element_text(face = "bold", color = "black", 
                                 size = 15),
        axis.title = element_text(face = "bold", color = "black", 
                                  size = 20),
        axis.ticks = element_blank(),legend.position="none")+labs(y="PLSc2 - MSF3",x="PLSc2 - MCF3")+scale_color_brewer(palette="Dark2")
# p2.2

p2.3<-ggplot(data = all.data.pc2[all.data.pc1$set==1,])+
  geom_point(aes(x=X_3,y=Y_3,color=sex),size=2)+
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill="white",color="white"),
        axis.line =element_line(color = "black", size = 1, linetype = "solid"),
        axis.text = element_text(face = "bold", color = "black", 
                                 size = 15),
        axis.title = element_text(face = "bold", color = "black", 
                                  size = 20),
        axis.ticks = element_blank(),legend.position="none")+labs(y="PLSc2 - MSF3",x="PLSc2 - MCF3")+scale_color_brewer(palette="Dark2")
# p2.3



p3.1<-ggplot(data = all.data.pc1[all.data.pc1$set==1,])+
  geom_point(aes(x=X_7*-1,y=Y_7*-1,color=as.numeric(demo_comb_income_v2)),size=2)+
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill="white",color="white"),
        axis.line =element_line(color = "black", size = 1, linetype = "solid"),
        axis.text = element_text(face = "bold", color = "black", 
                                 size = 15),
        axis.title = element_text(face = "bold", color = "black", 
                                  size = 20),
        axis.ticks = element_blank(),
        legend.title = element_text(face="bold",color="black",size=10))+labs(y="PLSc1 - MSF7", 
                                                                             x="PLSc1 - MCF7",color="Family Income")+scale_color_continuous(type="viridis")
# p3.1

p3.2<-ggplot(data = all.data.pc1[all.data.pc1$set==1,])+
  geom_point(aes(x=X_7*-1,y=Y_7*-1,color=race_ethc),size=2)+
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill="white",color="white"),
        axis.line =element_line(color = "black", size = 1, linetype = "solid"),
        axis.text = element_text(face = "bold", color = "black", 
                                 size = 15),
        axis.title = element_text(face = "bold", color = "black", 
                                  size = 20),
        axis.ticks = element_blank(),
        legend.title = element_text(face="bold",color="black",size=10))+labs(y="PLSc1 - MSF7", 
                                                                             x="PLSc1 - MCF7",color="Race/Ethnicity")+scale_color_brewer(palette="Dark2")
# p3.2

p3.3<-ggplot(data = all.data.pc1[all.data.pc1$set==1,])+
  geom_point(aes(x=X_7*-1,y=Y_7*-1,color=sex),size=2)+
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill="white",color="white"),
        axis.line =element_line(color = "black", size = 1, linetype = "solid"),
        axis.text = element_text(face = "bold", color = "black", 
                                 size = 15),
        axis.title = element_text(face = "bold", color = "black", 
                                  size = 20),
        axis.ticks = element_blank(),
        legend.title = element_text(face="bold",color="black",size=10))+labs(y="PLSc1 - MSF7", 
                                                                             x="PLSc1 - MCF7",color="Sex")+scale_color_brewer(palette="Dark2")
# p3.3

# 
# # library(cowplot)
# plot_grid(
#   p1.1,p1.2,p1.3,p2.1,p2.2,p2.3,p3.1,p3.2,p3.3, ncol = 3
# )


plot_grid(
  p1.1,p2.1,p3.1,p1.2,p2.2,p3.2,p1.3,p2.3,p3.3, ncol = 3
)


###Analyses to add to this supplementary figure###
#1. 
#X FOR FIGURE
cor.test(all.data.pc1[all.data.pc1$set==1,c("X_1")]*-1,as.numeric(as.character(all.data.pc1[all.data.pc1$set==1,c("demo_comb_income_v2")])))
cor.test(all.data.pc2[all.data.pc2$set==1,c("X_3")],as.numeric(as.character(all.data.pc2[all.data.pc2$set==1,c("demo_comb_income_v2")])))
cor.test(all.data.pc1[all.data.pc1$set==1,c("X_7")]*-1,as.numeric(as.character(all.data.pc1[all.data.pc1$set==1,c("demo_comb_income_v2")])))

#Y FOR FIGURE
cor.test(all.data.pc1[all.data.pc1$set==1,c("Y_1")]*-1,as.numeric(as.character(all.data.pc1[all.data.pc1$set==1,c("demo_comb_income_v2")]))) #sig
cor.test(all.data.pc2[all.data.pc2$set==1,c("Y_3")],as.numeric(as.character(all.data.pc2[all.data.pc2$set==1,c("demo_comb_income_v2")])))
cor.test(all.data.pc1[all.data.pc1$set==1,c("Y_7")]*-1,as.numeric(as.character(all.data.pc1[all.data.pc1$set==1,c("demo_comb_income_v2")])))



#2.
#numerize sex
sex_train<-ifelse(all.data.pc1[all.data.pc1$set==1,c("sex")]=="F",1,0)

#X FOR FIGURE
cor.test(all.data.pc1[all.data.pc1$set==1,c("X_1")]*-1,sex_train)
cor.test(all.data.pc2[all.data.pc2$set==1,c("X_3")],sex_train)
cor.test(all.data.pc1[all.data.pc1$set==1,c("X_7")]*-1,sex_train)

#Y FOR FIGURE
cor.test(all.data.pc1[all.data.pc1$set==1,c("Y_1")]*-1,sex_train)
cor.test(all.data.pc2[all.data.pc2$set==1,c("Y_3")],sex_train)
cor.test(all.data.pc1[all.data.pc1$set==1,c("Y_7")]*-1,sex_train)


#3
#Dummy-code ethnicity
ethnicity_dc<-as.data.frame(ifelse(all.data.pc1[all.data.pc1$set==1,c("race_ethc")]=="Non-Hispanic White",1,0))
colnames(ethnicity_dc)<-"dceth_w"

ethnicity_dc$dceth_h<-ifelse(all.data.pc1[all.data.pc1$set==1,c("race_ethc")]=="Hispanic",1,0)
ethnicity_dc$dceth_b<-ifelse(all.data.pc1[all.data.pc1$set==1,c("race_ethc")]=="Non-Hispanic Black",1,0)
ethnicity_dc$dceth_as<-ifelse(all.data.pc1[all.data.pc1$set==1,c("race_ethc")]=="Asian",1,0)
ethnicity_dc$dceth_ai<-ifelse(all.data.pc1[all.data.pc1$set==1,c("race_ethc")]=="AIAN",1,0)
ethnicity_dc$dceth_n<-ifelse(all.data.pc1[all.data.pc1$set==1,c("race_ethc")]=="NHPI",1,0)
ethnicity_dc$dceth_o<-ifelse(all.data.pc1[all.data.pc1$set==1,c("race_ethc")]=="OTHER",1,0)
ethnicity_dc$dceth_m<-ifelse(all.data.pc1[all.data.pc1$set==1,c("race_ethc")]=="MULTIPLE",1,0)

#X FOR FIGURE:
#Note: commented out subgroups with n=0
c1<-cor(all.data.pc1[all.data.pc1$set==1,c("X_1")]*-1,ethnicity_dc$dceth_w)
c2<-cor(all.data.pc1[all.data.pc1$set==1,c("X_1")]*-1,ethnicity_dc$dceth_h)
c3<-cor(all.data.pc1[all.data.pc1$set==1,c("X_1")]*-1,ethnicity_dc$dceth_b)
c4<-cor(all.data.pc1[all.data.pc1$set==1,c("X_1")]*-1,ethnicity_dc$dceth_as)
# c5<-cor(all.data.pc1[all.data.pc1$set==1,c("X_1")],ethnicity_dc$dceth_ai) #because there are 0 cases
# c6<-cor(all.data.pc1[all.data.pc1$set==1,c("X_1")],ethnicity_dc$dceth_n) #because there are 0 cases
c7<-cor(all.data.pc1[all.data.pc1$set==1,c("X_1")]*-1,ethnicity_dc$dceth_o)
c8<-cor(all.data.pc1[all.data.pc1$set==1,c("X_1")]*-1,ethnicity_dc$dceth_m)
cors<-rbind(c1,c2,c3,c4,c7,c8)
# cors<-ifelse(is.na(cors),0,cors)
barplot(t(cors))



c1<-cor(all.data.pc1[all.data.pc2$set==1,c("X_3")],ethnicity_dc$dceth_w)
c2<-cor(all.data.pc1[all.data.pc2$set==1,c("X_3")],ethnicity_dc$dceth_h)
c3<-cor(all.data.pc1[all.data.pc2$set==1,c("X_3")],ethnicity_dc$dceth_b)
c4<-cor(all.data.pc1[all.data.pc2$set==1,c("X_3")],ethnicity_dc$dceth_as)
# c5<-cor(all.data.pc1[all.data.pc2$set==1,c("X_3")],ethnicity_dc$dceth_ai)
# c6<-cor(all.data.pc1[all.data.pc2$set==1,c("X_3")],ethnicity_dc$dceth_n)
c7<-cor(all.data.pc1[all.data.pc2$set==1,c("X_3")],ethnicity_dc$dceth_o)
c8<-cor(all.data.pc1[all.data.pc2$set==1,c("X_3")],ethnicity_dc$dceth_m)
cors<-rbind(c1,c2,c3,c4,c7,c8)
# cors<-ifelse(is.na(cors),0,cors)
barplot(t(cors))


c1<-cor(all.data.pc1[all.data.pc1$set==1,c("X_7")]*-1,ethnicity_dc$dceth_w)
c2<-cor(all.data.pc1[all.data.pc1$set==1,c("X_7")]*-1,ethnicity_dc$dceth_h)
c3<-cor(all.data.pc1[all.data.pc1$set==1,c("X_7")]*-1,ethnicity_dc$dceth_b)
c4<-cor(all.data.pc1[all.data.pc1$set==1,c("X_7")]*-1,ethnicity_dc$dceth_as)
# c5<-cor(all.data.pc1[all.data.pc1$set==1,c("X_7")],ethnicity_dc$dceth_ai)
# c6<-cor(all.data.pc1[all.data.pc1$set==1,c("X_7")],ethnicity_dc$dceth_n)
c7<-cor(all.data.pc1[all.data.pc1$set==1,c("X_7")]*-1,ethnicity_dc$dceth_o)
c8<-cor(all.data.pc1[all.data.pc1$set==1,c("X_7")]*-1,ethnicity_dc$dceth_m)
cors<-rbind(c1,c2,c3,c4,c7,c8)
# cors<-ifelse(is.na(cors),0,cors)
barplot(t(cors))

#Y FOR FIGURE:
c1<-cor(all.data.pc1[all.data.pc1$set==1,c("Y_1")]*-1,ethnicity_dc$dceth_w)
c2<-cor(all.data.pc1[all.data.pc1$set==1,c("Y_1")]*-1,ethnicity_dc$dceth_h)
c3<-cor(all.data.pc1[all.data.pc1$set==1,c("Y_1")]*-1,ethnicity_dc$dceth_b)
c4<-cor(all.data.pc1[all.data.pc1$set==1,c("Y_1")]*-1,ethnicity_dc$dceth_as)
# c5<-cor(all.data.pc1[all.data.pc1$set==1,c("Y_1")],ethnicity_dc$dceth_ai)
# c6<-cor(all.data.pc1[all.data.pc1$set==1,c("Y_1")],ethnicity_dc$dceth_n)
c7<-cor(all.data.pc1[all.data.pc1$set==1,c("Y_1")]*-1,ethnicity_dc$dceth_o)
c8<-cor(all.data.pc1[all.data.pc1$set==1,c("Y_1")]*-1,ethnicity_dc$dceth_m)
cors<-rbind(c1,c2,c3,c4,c7,c8)
# cors<-ifelse(is.na(cors),0,cors)
barplot(t(cors))



c1<-cor(all.data.pc1[all.data.pc2$set==1,c("Y_3")],ethnicity_dc$dceth_w)
c2<-cor(all.data.pc1[all.data.pc2$set==1,c("Y_3")],ethnicity_dc$dceth_h)
c3<-cor(all.data.pc1[all.data.pc2$set==1,c("Y_3")],ethnicity_dc$dceth_b)
c4<-cor(all.data.pc1[all.data.pc2$set==1,c("Y_3")],ethnicity_dc$dceth_as)
# c5<-cor(all.data.pc1[all.data.pc2$set==1,c("Y_3")],ethnicity_dc$dceth_ai)
# c6<-cor(all.data.pc1[all.data.pc2$set==1,c("Y_3")],ethnicity_dc$dceth_n)
c7<-cor(all.data.pc1[all.data.pc2$set==1,c("Y_3")],ethnicity_dc$dceth_o)
c8<-cor(all.data.pc1[all.data.pc2$set==1,c("Y_3")],ethnicity_dc$dceth_m)
cors<-rbind(c1,c2,c3,c4,c7,c8)
# cors<-ifelse(is.na(cors),0,cors)
barplot(t(cors))


c1<-cor(all.data.pc1[all.data.pc1$set==1,c("Y_7")]*-1,ethnicity_dc$dceth_w)
c2<-cor(all.data.pc1[all.data.pc1$set==1,c("Y_7")]*-1,ethnicity_dc$dceth_h)
c3<-cor(all.data.pc1[all.data.pc1$set==1,c("Y_7")]*-1,ethnicity_dc$dceth_b)
c4<-cor(all.data.pc1[all.data.pc1$set==1,c("Y_7")]*-1,ethnicity_dc$dceth_as)
# c5<-cor(all.data.pc1[all.data.pc1$set==1,c("Y_7")],ethnicity_dc$dceth_ai)
# c6<-cor(all.data.pc1[all.data.pc1$set==1,c("Y_7")],ethnicity_dc$dceth_n)
c7<-cor(all.data.pc1[all.data.pc1$set==1,c("Y_7")]*-1,ethnicity_dc$dceth_o)
c8<-cor(all.data.pc1[all.data.pc1$set==1,c("Y_7")]*-1,ethnicity_dc$dceth_m)
cors<-rbind(c1,c2,c3,c4,c7,c8)
# cors<-ifelse(is.na(cors),0,cors)
barplot(t(cors))

####Analyses of other external variables####
##NOTE: all of these analyses (and the scatter plots) show entire dataset (306)
#create df that only contains the expression of conn LFs and t since inj

#T Since Injury
#PC1
cor.test(all.data.pc1$X_1*-1,all.data.pc1$t_since_inj)
cor.test(all.data.pc1$X_2,all.data.pc1$t_since_inj)
cor.test(all.data.pc1$X_3,all.data.pc1$t_since_inj)
# cor.test(all.data.pc1$X_4,all.data.pc1$t_since_inj)
cor.test(all.data.pc1$X_5,all.data.pc1$t_since_inj)
# cor.test(all.data.pc1$X_6,all.data.pc1$t_since_inj)
cor.test(all.data.pc1$X_7,all.data.pc1$t_since_inj)
cor.test(all.data.pc1$X_8,all.data.pc1$t_since_inj)
cor.test(all.data.pc1$X_9,all.data.pc1$t_since_inj)#r=0.11, p=0.048
cor.test(all.data.pc1$X_10,all.data.pc1$t_since_inj)
cor.test(all.data.pc1$X_11,all.data.pc1$t_since_inj)
cor.test(all.data.pc1$X_12,all.data.pc1$t_since_inj)
cor.test(all.data.pc1$X_13,all.data.pc1$t_since_inj)
cor.test(all.data.pc1$X_14,all.data.pc1$t_since_inj)
cor.test(all.data.pc1$X_15,all.data.pc1$t_since_inj)
cor.test(all.data.pc1$X_16,all.data.pc1$t_since_inj)
cor.test(all.data.pc1$X_17,all.data.pc1$t_since_inj)
cor.test(all.data.pc1$X_18,all.data.pc1$t_since_inj)
# cor.test(all.data.pc1$X_19,all.data.pc1$t_since_inj)

#PC2
# cor.test(all.data.pc2$X_1,all.data.pc1$t_since_inj)
cor.test(all.data.pc2$X_2,all.data.pc1$t_since_inj)
cor.test(all.data.pc2$X_3,all.data.pc1$t_since_inj)
# cor.test(all.data.pc2$X_4,all.data.pc1$t_since_inj)
cor.test(all.data.pc2$X_5,all.data.pc1$t_since_inj)
# cor.test(all.data.pc1$X_6,all.data.pc1$t_since_inj)
# cor.test(all.data.pc2$X_7,all.data.pc1$t_since_inj)
cor.test(all.data.pc2$X_8,all.data.pc1$t_since_inj)
cor.test(all.data.pc2$X_9,all.data.pc1$t_since_inj)
cor.test(all.data.pc2$X_10,all.data.pc1$t_since_inj)
cor.test(all.data.pc2$X_11,all.data.pc1$t_since_inj)
# cor.test(all.data.pc2$X_12,all.data.pc1$t_since_inj)
# cor.test(all.data.pc2$X_13,all.data.pc1$t_since_inj)
# cor.test(all.data.pc2$X_14,all.data.pc1$t_since_inj)
# cor.test(all.data.pc2$X_15,all.data.pc1$t_since_inj)
# cor.test(all.data.pc2$X_16,all.data.pc1$t_since_inj)
# cor.test(all.data.pc2$X_17,all.data.pc1$t_since_inj)
# cor.test(all.data.pc2$X_18,all.data.pc1$t_since_inj)
# cor.test(all.data.pc1$X_19,all.data.pc1$t_since_inj)



#Injury Mechanisms
tmp_Xcols<-colnames(all.data.pc1)[regexpr("X",colnames(all.data.pc1))>0]
tmp.pc1.xLFs.nuisancevars<-all.data.pc1[all.data.pc1$set==1,c(tmp_Xcols,"causes")]
apply(tmp.pc1.xLFs.nuisancevars,2,function(x){kruskal.test(x~causes,data=tmp.pc1.xLFs.nuisancevars)})

#X_2
#X_15

ggplot(tmp.pc1.xLFs.nuisancevars, aes(x = factor(causes), y = X_2)) + 
  geom_bar(stat = "summary", fun = "mean",fill="black")+
  xlab("Injury Mechanism")+
  ylab("Mean")+
  theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size=15,face="bold"),
        axis.title.x=element_text(size=20,face="bold"),
        axis.text.y = element_text(size=15,face="bold"),
        axis.title.y=element_text(size=20,face="bold"),
        legend.position = "top"
  ) 



ggplot(tmp.pc1.xLFs.nuisancevars, aes(x = factor(causes), y = X_15)) + 
  geom_bar(stat = "summary", fun = "mean",fill="black")+
  xlab("Injury Mechanism")+
  ylab("Mean")+
  theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size=15,face="bold"),
        axis.title.x=element_text(size=20,face="bold"),
        axis.text.y = element_text(size=15,face="bold"),
        axis.title.y=element_text(size=20,face="bold"),
        legend.position = "top"
  ) 


tmp.pc1.xLFs.nuisancevars<-all.data.pc1[all.data.pc1$set==1,c(tmp_Xcols,"totalTBIs")]

apply(tmp.pc1.xLFs.nuisancevars,2,function(x){kruskal.test(x~totalTBIs,data=tmp.pc1.xLFs.nuisancevars)})


#X_15

#Switch to PC1 before running
ggplot(tmp.pc1.xLFs.nuisancevars, aes(x = factor(totalTBIs), y = X_15)) + 
  geom_bar(stat = "summary", fun = "mean",fill="black")+
  xlab("Total TBIs")+
  ylab("Mean")+
  theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size=15,face="bold"),
        axis.title.x=element_text(size=20,face="bold"),
        axis.text.y = element_text(size=15,face="bold"),
        axis.title.y=element_text(size=20,face="bold"),
        legend.position = "top"
  ) 



##NOW DO AGAIN FOR PC2
#Switch to PC2 before running
tmp_Xcols<-colnames(all.data.pc2)[regexpr("X",colnames(all.data.pc2))>0]
tmp.pc2.xLFs.nuisancevars<-all.data.pc2[all.data.pc2$set==1,c(tmp_Xcols,"causes")]
apply(tmp.pc2.xLFs.nuisancevars,2,function(x){kruskal.test(x~causes,data=tmp.pc2.xLFs.nuisancevars)})
#Non-sig


tmp.pc2.xLFs.nuisancevars<-all.data.pc2[all.data.pc2$set==1,c(tmp_Xcols,"totalTBIs")]

apply(tmp.pc2.xLFs.nuisancevars,2,function(x){kruskal.test(x~totalTBIs,data=tmp.pc2.xLFs.nuisancevars)})
#Non-sig



####ADJACENCY MATRIX ANALYSIS (Figure 5)####
#30_36 (lpOPER-lPreC) : 1018_1024
#61_65 (rpOPER-rpostC): 2018_2022
#67_71 (rpreC-rSFG): 2024_2028
#61_74 (rpOPER-rSMAR): 2018_2031

xMFConns<-boot.pc1$bootRatiosSignificant.i
Yboots<-boot.pc1$bootRatiosSignificant.j


adjmat_CbyS<-expand_grid(rownames(xMFConns),rownames(Yboots))
colnames(adjmat_CbyS)<-c("xConns","ySymps")

adjmat_CbyS$gen_prob<-ifelse(adjmat_CbyS$ySymps=="cbcl_q56b_p","headache",
                             ifelse(adjmat_CbyS$ySymps=="cbcl_q56c_p","nausea",
                                    ifelse(adjmat_CbyS$ySymps=="cbcl_q56g_p","vomiting",
                                           ifelse(adjmat_CbyS$ySymps=="cbcl_q51_p","dizziness",
                                                  ifelse(adjmat_CbyS$ySymps=="cbcl_q54_p","fatigue",
                                                         ifelse(adjmat_CbyS$ySymps=="sleepdisturb25_p","drowsiness",
                                                                ifelse(adjmat_CbyS$ySymps=="sleepdisturb4_p","troublesleeping",
                                                                       ifelse(adjmat_CbyS$ySymps=="cbcl_q77_p","sleep+",
                                                                              ifelse(adjmat_CbyS$ySymps=="cbcl_q76_p","sleep-",
                                                                                     ifelse(adjmat_CbyS$ySymps=="cbcl_scr_dsm5_anxdisord_t","anxiety",
                                                                                            ifelse(adjmat_CbyS$ySymps=="cbcl_scr_dsm5_depress_t","depression",
                                                                                                   ifelse(adjmat_CbyS$ySymps=="cbcl_scr_syn_aggressive_t","aggression",
                                                                                                          ifelse(adjmat_CbyS$ySymps=="cbcl_scr_syn_attention_t","attention",
                                                                                                                 ifelse(adjmat_CbyS$ySymps=="nihtbx_cardsort_fc","cardsorting",
                                                                                                                        ifelse(adjmat_CbyS$ySymps=="nihtbx_list_fc","work-memory",
                                                                                                                               ifelse(adjmat_CbyS$ySymps=="nihtbx_pattern_fc","proc-speed",
                                                                                                                                      ifelse(adjmat_CbyS$ySymps=="nihtbx_picture_fc","seq-memory",
                                                                                                                                             ifelse(adjmat_CbyS$ySymps=="pea_ravlt_sd_trial_vi_tc","s-recall",
                                                                                                                                                    ifelse(adjmat_CbyS$ySymps=="pea_ravlt_ld_trial_vii_tc","l-recall",
                                                                                                                                                           "NANA")))))))))))))))))))
#LEGEND:
#blue=somatic
#red=sleep
#green=mood
#purple=cognitive
adjmat_CbyS$prob_domain<-ifelse(adjmat_CbyS$ySymps=="cbcl_q56b_p",1,
                             ifelse(adjmat_CbyS$ySymps=="cbcl_q56c_p",1,
                                    ifelse(adjmat_CbyS$ySymps=="cbcl_q56g_p",1,
                                           ifelse(adjmat_CbyS$ySymps=="cbcl_q51_p",1,
                                                  ifelse(adjmat_CbyS$ySymps=="cbcl_q54_p",1,
                                                         ifelse(adjmat_CbyS$ySymps=="sleepdisturb25_p",2,
                                                                ifelse(adjmat_CbyS$ySymps=="sleepdisturb4_p",2,
                                                                       ifelse(adjmat_CbyS$ySymps=="cbcl_q77_p",2,
                                                                              ifelse(adjmat_CbyS$ySymps=="cbcl_q76_p",2,
                                                                                     ifelse(adjmat_CbyS$ySymps=="cbcl_scr_dsm5_anxdisord_t",3,
                                                                                            ifelse(adjmat_CbyS$ySymps=="cbcl_scr_dsm5_depress_t",3,
                                                                                                   ifelse(adjmat_CbyS$ySymps=="cbcl_scr_syn_aggressive_t",3,
                                                                                                          ifelse(adjmat_CbyS$ySymps=="cbcl_scr_syn_attention_t",4,
                                                                                                                 ifelse(adjmat_CbyS$ySymps=="nihtbx_cardsort_fc",4,
                                                                                                                        ifelse(adjmat_CbyS$ySymps=="nihtbx_list_fc",4,
                                                                                                                               ifelse(adjmat_CbyS$ySymps=="nihtbx_pattern_fc",4,
                                                                                                                                      ifelse(adjmat_CbyS$ySymps=="nihtbx_picture_fc",4,
                                                                                                                                             ifelse(adjmat_CbyS$ySymps=="pea_ravlt_sd_trial_vi_tc",4,
                                                                                                                                                    ifelse(adjmat_CbyS$ySymps=="pea_ravlt_ld_trial_vii_tc",4,
                                                                                                                                                           5)))))))))))))))))))

adjmat_CbyS$adj<-0

for (i in 1:dim(xMFConns)[2]){
  #Loop through LFs
  tmpX<-rownames(xMFConns)[xMFConns[,i]]#for each lf, store the names of the sig conns
  tmpY<-rownames(Yboots)[Yboots[,i]]#for each lf, store the names of the symptoms
  adjmat_CbyS[adjmat_CbyS$xConns%in%tmpX&adjmat_CbyS$ySymps%in%tmpY,c("adj")]<-adjmat_CbyS[adjmat_CbyS$xConns%in%tmpX&adjmat_CbyS$ySymps%in%tmpY,c("adj")]+1 #subselect the conn-symp combo, and add 1 onto itself
}

adjmat_CbyS$start<-substr(adjmat_CbyS$xConns,1,regexpr("_",adjmat_CbyS$xConns)-1)

adjmat_CbyS<-adjmat_CbyS[order(adjmat_CbyS$prob_domain,as.numeric(as.character(adjmat_CbyS$start))),]

adjmat_CbyS$xConns <- factor(adjmat_CbyS$xConns, levels = unique(adjmat_CbyS$xConns))
adjmat_CbyS$ySymps <- factor(adjmat_CbyS$ySymps, levels = unique(adjmat_CbyS$ySymps))
adjmat_CbyS$gen_prob <- factor(adjmat_CbyS$gen_prob, levels = unique(adjmat_CbyS$gen_prob))


ggplot(adjmat_CbyS,aes(gen_prob,xConns,fill=as.factor(adj)))+geom_tile()+scale_fill_manual(values=c("white", "pink", "red","darkred"))+
# +scale_fill_distiller(palette = "Greys")+
theme(axis.text.x = element_text(angle = 315,vjust = 0.5, hjust=0,size = 20),
      axis.title = element_text(face = "bold", color = "black", 
                                size = 20),
      axis.ticks = element_blank())+
  labs(y="Connections",x="Symptoms")




#1008_1015 (left inferior parietal - left middle temporal) is a connection that is implicated in 3 MTMS pairs, along with trouble sleeping. 
sum(boot.pc1$bootRatiosSignificant.i[rownames(boot.pc1$bootRatiosSignificant.i)%in%"1008_1015",])#3 MTMS pairs
sum(boot.pc1$bootRatiosSignificant.j[rownames(boot.pc1$bootRatiosSignificant.j)%in%"sleepdisturb4_p",])#7 MTMS pairs

#Highlighting certain connection/symptom pairs
#1. Broad, non-specific symptom category connection correspondences
boot.pc1$bootRatiosSignificant.i[rownames(boot.pc1$bootRatiosSignificant.i)%in%c("2011_2025"),]#MTMS3; #Conn 2011_2025 (right lateral occipital - right precuneus) associated with most cognitive problems
boot.pc1$bootRatiosSignificant.i[rownames(boot.pc1$bootRatiosSignificant.i)%in%c("12_1026"),]#MTMS3; lPUT_lrACC (left putamen - left rostral anterior cingulate)
boot.pc1$bootRatiosSignificant.i[rownames(boot.pc1$bootRatiosSignificant.i)%in%c("12_1013"),]#MTMS3; lPUT_lLING (left putamen - left lingual)

#Visualize the polar plots for Figure 5
plot_data<-plsc.t09.pc1
boot_data<-boot.pc1

###APPLIES TO ALL POLAR PLOTS###
#Store loadings
df<-as.data.frame(plot_data[[4]]$TExPosition.Data$pdq$q)
#Change column names
colnames(df)<-paste("Y",substr(colnames(df),2,5),sep="_")
#Only keep the lfs that were significant
df<-df[,plot_data[[1]]<.05]
#Store symptoms as a new column
df$symptom<-rownames(df)
#Create a column to specify symptom domains (not really used atm)
df$domain<-c("mood","mood","cognitive","mood","somatic","somatic","somatic","somatic","somatic","somatic","somatic","somatic","somatic","cognitive","cognitive","cognitive","cognitive","cognitive","cognitive")
#Create a column for the common name of the symptom
df$gen_prob<-c("depression","anxiety","attention","aggression","headache","nausea","vomiting","dizziness","fatigue","sleep+","sleep-","drowsiness","troublesleeping","s-recall","l-recall","seq-memory","work-memory","proc-speed","cardsorting")
#Flag significant symptoms based on bootstrap tests
df$boottest<-boot_data$bootRatiosSignificant.j[,plot_data[[1]]<.05]
#Melt df
df_melt<-melt(df)

#Color palette
c25 <- c(
  "black", "black", # red
  "black",
  "#6A3D9A", # purple
  "black", # orange
  "black", "black",
  "black", "black", # lt pink
  "palegreen2",
  "black", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "black", "black", "black", "black", "steelblue4",
  "black", "black", "black", "black",
  "black", "black"
)

#Create df for labels
label_data <- as.data.frame(sort(df$gen_prob))
colnames(label_data)<-"gen_prob"
label_data$id<-seq(1,19)
# calculate the ANGLE of the labels
number_of_bar <- nrow(label_data)
angle <-  90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)

###@@@###

#Note: will reverse a few vars for visualization purposes
df_melt[df_melt$variable=="Y_3",c("value")]<-df_melt[df_melt$variable=="Y_3",c("value")]*-1

###SPECIFIC POLAR PLOTS###
p3<-create_polar_plots(df_melt,var = "Y_3",color = "black",colorpal = c25,ylims = c(-0.6,0.8),boot=TRUE,label_data)



#Color palette
c25 <- c(
  "black", "black", # red
  "green4",
  "black", # purple
  "black", # orange
  "black", "black",
  "black", "black", # lt pink
  "black",
  "black", # lt purple
  "black", # lt orange
  "black", "black",
  "black", "black", "black", "black", "black",
  "black", "black", "black", "black",
  "black", "black"
)

#Note: will reverse a few vars for visualization purposes
df_melt[df_melt$variable=="Y_1",c("value")]<-df_melt[df_melt$variable=="Y_1",c("value")]*-1
df_melt[df_melt$variable=="Y_7",c("value")]<-df_melt[df_melt$variable=="Y_7",c("value")]*-1

###SPECIFIC POLAR PLOTS###
p1<-create_polar_plots(df_melt,var = "Y_1",color = "black",colorpal = c25,ylims = c(-0.5,0.8),boot=TRUE,label_data)
p7<-create_polar_plots(df_melt,var = "Y_7",color = "black",colorpal = c25,ylims = c(-0.8,0.8),boot=TRUE,label_data)
p8<-create_polar_plots(df_melt,var = "Y_8",color = "black",colorpal = c25,ylims = c(-0.6,0.8),boot=TRUE,label_data)
p16<-create_polar_plots(df_melt,var = "Y_16",color = "black",colorpal = c25,ylims = c(-0.8,0.9),boot=TRUE,label_data)
plot_grid(p1,p7,p8,p16)


####Create visualization of all LFs####
plot_data<-plsc.t09.pc1
boot_data<-boot.pc1

###APPLIES TO ALL POLAR PLOTS###
#Store loadings
df<-as.data.frame(plot_data[[4]]$TExPosition.Data$pdq$q)
#Change column names
colnames(df)<-paste("Y",substr(colnames(df),2,5),sep="_")
#Only keep the lfs that were significant
df<-df[,plot_data[[1]]<.05]
#Store symptoms as a new column
df$symptom<-rownames(df)
#Create a column to specify symptom domains (not really used atm)
df$domain<-c("mood","mood","cognitive","mood","somatic","somatic","somatic","somatic","somatic","somatic","somatic","somatic","somatic","cognitive","cognitive","cognitive","cognitive","cognitive","cognitive")
#Create a column for the common name of the symptom
df$gen_prob<-c("depression","anxiety","attention","aggression","headache","nausea","vomiting","dizziness","fatigue","sleep+","sleep-","drowsiness","troublesleeping","s-recall","l-recall","seq-memory","work-memory","proc-speed","cardsorting")
#Flag significant symptoms based on bootstrap tests
df$boottest<-boot_data$bootRatiosSignificant.j[,plot_data[[1]]<.05]
#Melt df
df_melt<-melt(df)

#Create color palette
c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "sienna", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)


#Create df for labels
label_data <- as.data.frame(sort(df$gen_prob))
colnames(label_data)<-"gen_prob"
label_data$id<-seq(1,19)
# calculate the ANGLE of the labels
number_of_bar <- nrow(label_data)
angle <-  90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)

###@@@###

#Note: will reverse a few vars for visualization purposes
df_melt[df_melt$variable=="Y_1",c("value")]<-df_melt[df_melt$variable=="Y_1",c("value")]*-1
df_melt[df_melt$variable=="Y_2",c("value")]<-df_melt[df_melt$variable=="Y_2",c("value")]*-1
df_melt[df_melt$variable=="Y_3",c("value")]<-df_melt[df_melt$variable=="Y_3",c("value")]*-1
df_melt[df_melt$variable=="Y_7",c("value")]<-df_melt[df_melt$variable=="Y_7",c("value")]*-1

###SPECIFIC POLAR PLOTS###
p1<-create_polar_plots(df_melt,var = "Y_1",color = "black",colorpal = c25,ylims = c(-0.5,0.8),boot=TRUE,label_data)
p2<-create_polar_plots(df_melt,var = "Y_2",color = "black",colorpal = c25,ylims = c(-0.9,0.9),boot=TRUE,label_data)
p3<-create_polar_plots(df_melt,var = "Y_3",color = "black",colorpal = c25,ylims = c(-0.6,0.8),boot=TRUE,label_data)
# p4<-create_polar_plots(df_melt,var = "Y_4",color = "black",colorpal = c25,ylims = c(-0.7,0.8),boot=TRUE,label_data)
p5<-create_polar_plots(df_melt,var = "Y_5",color = "black",colorpal = c25,ylims = c(-0.6,0.8),boot=TRUE,label_data)
# p6<-create_polar_plots(df_melt,var = "Y_6",color = "black",colorpal = c25,ylims = c(-0.7,0.8),boot=TRUE,label_data)
p7<-create_polar_plots(df_melt,var = "Y_7",color = "black",colorpal = c25,ylims = c(-0.8,0.8),boot=TRUE,label_data)
p8<-create_polar_plots(df_melt,var = "Y_8",color = "black",colorpal = c25,ylims = c(-0.6,0.8),boot=TRUE,label_data)
p9<-create_polar_plots(df_melt,var = "Y_9",color = "black",colorpal = c25,ylims = c(-0.8,0.8),boot=TRUE,label_data)
p10<-create_polar_plots(df_melt,var = "Y_10",color = "black",colorpal = c25,ylims = c(-0.7,0.8),boot=TRUE,label_data)
p11<-create_polar_plots(df_melt,var = "Y_11",color = "black",colorpal = c25,ylims = c(-0.6,0.8),boot=TRUE,label_data)
p12<-create_polar_plots(df_melt,var = "Y_12",color = "black",colorpal = c25,ylims = c(-0.6,0.8),boot=TRUE,label_data)
p13<-create_polar_plots(df_melt,var = "Y_13",color = "black",colorpal = c25,ylims = c(-0.8,0.8),boot=TRUE,label_data)
p14<-create_polar_plots(df_melt,var = "Y_14",color = "black",colorpal = c25,ylims = c(-0.8,0.8),boot=TRUE,label_data)
p15<-create_polar_plots(df_melt,var = "Y_15",color = "black",colorpal = c25,ylims = c(-0.8,0.8),boot=TRUE,label_data)
p16<-create_polar_plots(df_melt,var = "Y_16",color = "black",colorpal = c25,ylims = c(-0.8,0.9),boot=TRUE,label_data)
p17<-create_polar_plots(df_melt,var = "Y_17",color = "black",colorpal = c25,ylims = c(-0.9,0.9),boot=TRUE,label_data)
p18<-create_polar_plots(df_melt,var = "Y_18",color = "black",colorpal = c25,ylims = c(-0.8,0.8),boot=TRUE,label_data)
# p19<-create_polar_plots(df_melt,var = "Y_19",color = "black",colorpal = c25,ylims = c(-0.8,0.8),boot=TRUE,label_data)

plot_grid(p1,p2,p3,p5,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18)



#PC2
plot_data<-plsc.t09.pc2
boot_data<-boot.pc2

###APPLIES TO ALL POLAR PLOTS###
#Store loadings
df<-as.data.frame(plot_data[[4]]$TExPosition.Data$pdq$q)
#Change column names
colnames(df)<-paste("Y",substr(colnames(df),2,5),sep="_")
#Only keep the lfs that were significant
df<-df[,plot_data[[1]]<.05]
#Store symptoms as a new column
df$symptom<-rownames(df)
#Create a column to specify symptom domains (not really used atm)
df$domain<-c("mood","mood","cognitive","mood","somatic","somatic","somatic","somatic","somatic","somatic","somatic","somatic","somatic","cognitive","cognitive","cognitive","cognitive","cognitive","cognitive")
#Create a column for the common name of the symptom
df$gen_prob<-c("depression","anxiety","attention","aggression","headache","nausea","vomiting","dizziness","fatigue","sleep+","sleep-","drowsiness","troublesleeping","s-recall","l-recall","seq-memory","work-memory","proc-speed","cardsorting")
#Flag significant symptoms based on bootstrap tests
df$boottest<-boot_data$bootRatiosSignificant.j[,plot_data[[1]]<.05]
#Melt df
df_melt<-melt(df)

#Create color palette
c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "sienna", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)


#Create df for labels
label_data <- as.data.frame(sort(df$gen_prob))
colnames(label_data)<-"gen_prob"
label_data$id<-seq(1,19)
# calculate the ANGLE of the labels
number_of_bar <- nrow(label_data)
angle <-  90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)

###@@@###

#Note: will reverse a few vars for visualization purposes
df_melt[df_melt$variable=="Y_1",c("value")]<-df_melt[df_melt$variable=="Y_1",c("value")]*-1
df_melt[df_melt$variable=="Y_2",c("value")]<-df_melt[df_melt$variable=="Y_2",c("value")]*-1
df_melt[df_melt$variable=="Y_3",c("value")]<-df_melt[df_melt$variable=="Y_3",c("value")]*-1
df_melt[df_melt$variable=="Y_7",c("value")]<-df_melt[df_melt$variable=="Y_7",c("value")]*-1

###SPECIFIC POLAR PLOTS###
# p1<-create_polar_plots(df_melt,var = "Y_1",color = "black",colorpal = c25,ylims = c(-0.5,0.8),boot=TRUE,label_data)
p2<-create_polar_plots(df_melt,var = "Y_2",color = "black",colorpal = c25,ylims = c(-0.9,0.9),boot=TRUE,label_data)
p3<-create_polar_plots(df_melt,var = "Y_3",color = "black",colorpal = c25,ylims = c(-0.6,0.8),boot=TRUE,label_data)
# p4<-create_polar_plots(df_melt,var = "Y_4",color = "black",colorpal = c25,ylims = c(-0.7,0.8),boot=TRUE,label_data)
p5<-create_polar_plots(df_melt,var = "Y_5",color = "black",colorpal = c25,ylims = c(-0.6,0.8),boot=TRUE,label_data)
# p6<-create_polar_plots(df_melt,var = "Y_6",color = "black",colorpal = c25,ylims = c(-0.7,0.8),boot=TRUE,label_data)
# p7<-create_polar_plots(df_melt,var = "Y_7",color = "black",colorpal = c25,ylims = c(-0.8,0.8),boot=TRUE,label_data)
p8<-create_polar_plots(df_melt,var = "Y_8",color = "black",colorpal = c25,ylims = c(-0.9,0.9),boot=TRUE,label_data)
p9<-create_polar_plots(df_melt,var = "Y_9",color = "black",colorpal = c25,ylims = c(-0.8,0.8),boot=TRUE,label_data)
p10<-create_polar_plots(df_melt,var = "Y_10",color = "black",colorpal = c25,ylims = c(-0.7,0.8),boot=TRUE,label_data)
p11<-create_polar_plots(df_melt,var = "Y_11",color = "black",colorpal = c25,ylims = c(-0.6,0.8),boot=TRUE,label_data)
# p12<-create_polar_plots(df_melt,var = "Y_12",color = "black",colorpal = c25,ylims = c(-0.6,0.8),boot=TRUE,label_data)
# p13<-create_polar_plots(df_melt,var = "Y_13",color = "black",colorpal = c25,ylims = c(-0.8,0.8),boot=TRUE,label_data)
# p14<-create_polar_plots(df_melt,var = "Y_14",color = "black",colorpal = c25,ylims = c(-0.8,0.8),boot=TRUE,label_data)
# p15<-create_polar_plots(df_melt,var = "Y_15",color = "black",colorpal = c25,ylims = c(-0.8,0.8),boot=TRUE,label_data)
# p16<-create_polar_plots(df_melt,var = "Y_16",color = "black",colorpal = c25,ylims = c(-0.8,0.9),boot=TRUE,label_data)
# p17<-create_polar_plots(df_melt,var = "Y_17",color = "black",colorpal = c25,ylims = c(-0.9,0.9),boot=TRUE,label_data)
# p18<-create_polar_plots(df_melt,var = "Y_18",color = "black",colorpal = c25,ylims = c(-0.8,0.8),boot=TRUE,label_data)
# p19<-create_polar_plots(df_melt,var = "Y_19",color = "black",colorpal = c25,ylims = c(-0.8,0.8),boot=TRUE,label_data)

plot_grid(p2,p3,p5,p8,p9,p10,p11)


####Figure 3 Analyses####
###Histogram of %COV explained###
perc.cove<-as.data.frame(plsc.t09.pc2[[4]]$TExPosition.Data$t)
colnames(perc.cove)<-"perc_cov"
perc.cove$lf<-seq(1:19)

p<-ggplot(data=perc.cove, aes(x=lf, y=perc_cov)) +
  geom_bar(color="gray",fill="gray",stat="identity")+
  labs(x="Multi-Tract Multi-Symptom Pair",y="%COV Explained")+
  theme_classic(base_size =20)+
  theme(axis.title=element_text(face="bold"),
        axis.text=element_text(face="bold"))
p

#PART BELOW WILL BE NECESSARY WHEN CREATING THE FIGURES
lfint<-3
col<-"turquoise"
p+annotate("bar", x = lfint, y = perc.cove$perc_cov[lfint], colour = col,fill=col)




#Scatter 1
#Create color column
all.data.pc1$psych_col<-ifelse(all.data.pc1[,c("ksads_14_853_p")]==1,"black","blueviolet")

#Create plot
cords<-c("X_1","Y_1")
p1<-ggplot(all.data.pc1[all.data.pc1$set==1,], aes(x=X_1*-1, y=Y_1*-1))+
  # geom_point(colour="blueviolet",aes(size=4.5,alpha=0.5,shape=as.factor(set)))+
  geom_point(colour=all.data.pc1[all.data.pc1$set==1,c("psych_col")],size=7.5,alpha=0.7)+
  # p1<-ggplot(all.data.pc2, aes(x=X_1, y=Y_1))+
  # geom_point(colour="blueviolet",aes(size=4.5,alpha=0.5))+
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill="white",color="white"),
        axis.line =element_line(color = "black", size = 1, linetype = "solid"),
        axis.text = element_text(face = "bold", color = "black", 
                                 size = 15),
        axis.title = element_text(face = "bold", color = "black", 
                                  size = 20),
        axis.ticks = element_blank())+
  labs(y="Multi-Symptom Feature 1", 
       x="Multi-Tract Connectivity Feature 1")+
  
  geom_text(data = all.data.pc1[all.data.pc1$subjectkey=="NDAR_INV4TB7WE3E",cords], label="1",col = 'white',size=6)+
  geom_text(data = all.data.pc1[all.data.pc1$subjectkey=="NDAR_INVPGE70ZD4",cords], label="2",col = 'white',size=6)+
  geom_text(data = all.data.pc1[all.data.pc1$subjectkey=="NDAR_INVAH6KK71Y",cords], label="3",col = 'white',size=6)+
  geom_text(data = all.data.pc1[all.data.pc1$subjectkey=="NDAR_INV0MPBK7TU",cords], label="4",col = 'white',size=6)+
  geom_text(data = all.data.pc1[all.data.pc1$subjectkey=="NDAR_INVTY9UFVZW",cords], label="5",col = 'white',size=6)+
  geom_text(data = all.data.pc1[all.data.pc1$subjectkey=="NDAR_INV0UMM15GY",cords], label="6",col = 'white',size=6)+
  theme(legend.position="none")

p1

#THE ONE ABOVE WORKS WELL


#Scatter 2
#Since there's a graphical outlier in x_1, I'll remove it here (see line below). Remember that it's for illustrative purposes only, we're not actually removing it for any analyses
# all.data.pc2.psych<-all.data.pc2.psych[all.data.pc2.psych$X_1<20,]

#Create color column
all.data.pc2$psych_col<-ifelse(all.data.pc2[,c("ksads_14_853_p")]==1,"black","turquoise")

cords<-c("X_3","Y_3")
p1<-ggplot(all.data.pc2[all.data.pc2$set==1,], aes(x=X_3, y=Y_3))+
  # geom_point(colour="blueviolet",aes(size=4.5,alpha=0.5,shape=as.factor(set)))+
  geom_point(colour=all.data.pc2[all.data.pc2$set==1,c("psych_col")],size=7.5,alpha=0.7)+
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill="white",color="white"),
        axis.line =element_line(color = "black", size = 1, linetype = "solid"),
        axis.text = element_text(face = "bold", color = "black", 
                                 size = 15),
        axis.title = element_text(face = "bold", color = "black", 
                                  size = 20),
        axis.ticks = element_blank())+
  labs(y="Multi-Symptom Feature 3", 
       x="Multi-Tract Connectivity Feature 3")+
  
  
  geom_text(data = all.data.pc2[all.data.pc2$subjectkey=="NDAR_INVPGE70ZD4",cords], label="1",col = 'white',size=6)+
  geom_text(data = all.data.pc2[all.data.pc2$subjectkey=="NDAR_INV4TB7WE3E",cords], label="2",col = 'white',size=6)+
  geom_text(data = all.data.pc2[all.data.pc2$subjectkey=="NDAR_INVAH6KK71Y",cords], label="3",col = 'white',size=6)+
  geom_text(data = all.data.pc2[all.data.pc2$subjectkey=="NDAR_INV0MPBK7TU",cords], label="4",col = 'white',size=6)+
  geom_text(data = all.data.pc2[all.data.pc2$subjectkey=="NDAR_INVTY9UFVZW",cords], label="5",col = 'white',size=6)+
  geom_text(data = all.data.pc2[all.data.pc2$subjectkey=="NDAR_INV0UMM15GY",cords], label="6",col = 'white',size=6)+
  theme(legend.position="none")
p1          



#PLSC1_LF7
all.data.pc1$psych_col<-ifelse(all.data.pc1[,c("ksads_14_853_p")]==1,"black","greenyellow")
#NOTE: removing one outlier for visualization purposes only. Not removed for the actual analyses.

cords<-c("X_7","Y_7")
p1<-ggplot(all.data.pc1[all.data.pc1$set==1&all.data.pc1$X_7<7,], aes(x=X_7*-1, y=Y_7*-1))+
  # geom_point(colour="blueviolet",aes(size=4.5,alpha=0.5,shape=as.factor(set)))+
  geom_point(colour=all.data.pc1[all.data.pc1$set==1&all.data.pc1$X_7<7,c("psych_col")],size=7.5,alpha=0.7)+
  # p1<-ggplot(all.data.pc2, aes(x=X_1, y=Y_1))+
  # geom_point(colour="blueviolet",aes(size=4.5,alpha=0.5))+
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill="white",color="white"),
        axis.line =element_line(color = "black", size = 1, linetype = "solid"),
        axis.text = element_text(face = "bold", color = "black", 
                                 size = 15),
        axis.title = element_text(face = "bold", color = "black", 
                                  size = 20),
        axis.ticks = element_blank())+
  labs(y="Multi-Symptom Feature 7", 
       x="Multi-Tract Connectivity Feature 7")+
  
  
  geom_text(data = all.data.pc1[all.data.pc1$subjectkey=="NDAR_INVPGE70ZD4",cords], label="1",col = 'white',size=6)+
  geom_text(data = all.data.pc1[all.data.pc1$subjectkey=="NDAR_INV4TB7WE3E",cords], label="2",col = 'white',size=6)+
  geom_text(data = all.data.pc1[all.data.pc1$subjectkey=="NDAR_INVAH6KK71Y",cords], label="3",col = 'white',size=6)+
  geom_text(data = all.data.pc1[all.data.pc1$subjectkey=="NDAR_INV0MPBK7TU",cords], label="4",col = 'white',size=6)+
  geom_text(data = all.data.pc1[all.data.pc1$subjectkey=="NDAR_INVTY9UFVZW",cords], label="5",col = 'white',size=6)+
  geom_text(data = all.data.pc1[all.data.pc1$subjectkey=="NDAR_INV0UMM15GY",cords], label="6",col = 'white',size=6)+
  theme(legend.position="none")
p1    


##ADD INDIVIDUAL SYMPTOM PROFILES##

#Scaling them here because they are scaled within the PLSC algorithm, but it would be helpful to show them scaled in this figure
raw_symps.train.scaled<-apply(mydata_t09$yy_train[,-which(colnames(mydata_t09$yy_train)%in%c("subjectkey"))],2,scale)
raw_symps.train.scaled<-as.data.frame(raw_symps.train.scaled)
raw_symps.train.scaled$subjectkey<-mydata_t09$yy_train$subjectkey

#List of set=1 subs who have ADHD

#CHANGE SUB HERE TO GET INDIVID SUB SYMPTOMS#

subint="NDAR_INVAH6KK71Y"
indiv_symps<-raw_symps.train.scaled[raw_symps.train.scaled$subjectkey==subint,-which(colnames(raw_symps.train.scaled[raw_symps.train.scaled$subjectkey==subint,])%in%c("subjectkey"))]
indiv_symps_melt<-melt(indiv_symps)
colnames(indiv_symps_melt)<-c("symptom","value")
indiv_symps_melt<-merge(indiv_symps_melt,df[,c("symptom","gen_prob")])
p2<-ggplot(indiv_symps_melt,aes(x=gen_prob,y=value,fill=gen_prob))+
  geom_bar(position = 'dodge', stat='identity')+scale_fill_manual(values=c25)+
  geom_hline(yintercept=1, linetype="dashed", color = "black")+
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill="white",color="white"),
        legend.position = "none",
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks=element_blank()
  )
p2


####Figure 3. Illustration of LFs (POLAR PLOTS)####
df<-as.data.frame(plsc.t09.pc1[[4]]$TExPosition.Data$pdq$q)
colnames(df)<-paste("Y",substr(colnames(df),2,5),sep="_")
df<-df[,plsc.t09.pc1[[1]]<.05]
# df<-as.data.frame(plsc.t09.pc1[[4]]$TExPosition.Data$pdq$q[,plsc.t09.pc1[[1]]<.05])
df$symptom<-rownames(df)
df$domain<-c("mood","mood","cognitive","mood","somatic","somatic","somatic","somatic","somatic","somatic","somatic","somatic","somatic","cognitive","cognitive","cognitive","cognitive","cognitive","cognitive")
df$gen_prob<-c("depression","anxiety","attention","aggression","headache","nausea","vomiting","dizziness","fatigue","sleep+","sleep-","drowsiness","troublesleeping","s-recall","l-recall","seq-memory","work-memory","proc-speed","cardsorting")
df$boottest<-boot.pc1$bootRatiosSignificant.j[,plsc.t09.pc1[[1]]<.05]
df_melt<-melt(df)



c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "sienna", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)


#Figure 3: polar plots showing loadings for LFs
# Get the name and the y position of each label

label_data <- as.data.frame(sort(df$gen_prob))
colnames(label_data)<-"gen_prob"
label_data$id<-seq(1,19)
# calculate the ANGLE of the labels
number_of_bar <- nrow(label_data)
angle <-  90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data<-merge(label_data,df_melt[df_melt$variable=="Y_1",c("gen_prob","value","boottest")])
# label_data$value<-label_data$value*-1
# calculate the alignment of labels: right or left
# If I am on the left part of the plot, my labels have currently an angle < -90
label_data$hjust<-ifelse( angle < -90, 1, 0)

# flip angle BY to make them readable
label_data$angle<-ifelse(angle < -90, angle+180, angle)


# p<-ggplot(data=df_melt[df_melt$variable=="V1",],aes(x=gen_prob,y=value*-1,fill=gen_prob)) + #IF MULTIPLYING BY -1
p<-ggplot(data=df_melt[df_melt$variable=="Y_1",],aes(x=gen_prob,y=value*-1,fill=gen_prob)) +
  geom_bar(position = 'dodge', stat='identity')+scale_fill_manual(values=c25)+
  coord_polar(start=0)+ylim(-0.4,0.75)+
  geom_bar(stat = "identity", aes(x = gen_prob, y = 0.01), fill = "white") +
  geom_bar(stat = "identity", aes(x = gen_prob, y = -0.04), fill = "white") +
  geom_bar(stat = "identity", aes(x = gen_prob, y = -0.03), fill = "blueviolet") +
  theme(legend.position="none",
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        axis.ticks=element_blank(),
        panel.background = element_rect(fill="white"))+
  geom_text(data=label_data, aes(x=id, y=.415, label=gen_prob, hjust=hjust,color=gen_prob), fontface="bold", size=10.5, angle= label_data$angle, inherit.aes = FALSE )+ scale_colour_manual(values=c25)
# +geom_text(data=label_data, aes(x=id, y=ifelse(value<0,value+.05,value-.05), label=ifelse(boottest==TRUE,"*",""), hjust=hjust), color="white",size=7.8, angle= label_data$angle, inherit.aes = FALSE )
p




label_data <- as.data.frame(sort(df$gen_prob))
colnames(label_data)<-"gen_prob"
label_data$id<-seq(1,19)
# calculate the ANGLE of the labels
number_of_bar <- nrow(label_data)
angle <-  90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data<-merge(label_data,df_melt[df_melt$variable=="Y_7",c("gen_prob","value","boottest")])
label_data$value<-label_data$value
# calculate the alignment of labels: right or left
# If I am on the left part of the plot, my labels have currently an angle < -90
label_data$hjust<-ifelse( angle < -90, 1, 0)

# flip angle BY to make them readable
label_data$angle<-ifelse(angle < -90, angle+180, angle)

p<-ggplot(data=df_melt[df_melt$variable=="Y_7",],aes(x=gen_prob,y=value*-1,fill=gen_prob)) +
  # p<-ggplot(data=df_melt[df_melt$variable=="Y_3",],aes(x=gen_prob,y=value,fill=gen_prob)) +
  geom_bar(position = 'dodge', stat='identity')+scale_fill_manual(values=c25)+
  coord_polar(start=0)+ylim(-0.8,1)+
  geom_bar(stat = "identity", aes(x = gen_prob, y = 0.01), fill = "white") +
  geom_bar(stat = "identity", aes(x = gen_prob, y = -0.04), fill = "white") +
  geom_bar(stat = "identity", aes(x = gen_prob, y = -0.03), fill = "greenyellow") +
  theme(legend.position="none",
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        axis.ticks=element_blank(),
        panel.background = element_rect(fill="white",color="white"))+
  geom_text(data=label_data, aes(x=id, y=.6, label=gen_prob, hjust=hjust,color=gen_prob), fontface="bold", size=10.5, angle= label_data$angle, inherit.aes = FALSE )+ scale_colour_manual(values=c25)
# +geom_text(data=label_data, aes(x=id, y=ifelse(value<0,value+.05,value-.06), label=ifelse(boottest==TRUE,"*",""), hjust=hjust), color="white",size=7.8, angle= label_data$angle, inherit.aes = FALSE )
p




df<-as.data.frame(plsc.t09.pc2[[4]]$TExPosition.Data$pdq$q)
colnames(df)<-paste("Y",substr(colnames(df),2,5),sep="_")
df<-df[,plsc.t09.pc2[[1]]<.05]
df$symptom<-rownames(df)
df$domain<-c("mood","mood","cognitive","mood","somatic","somatic","somatic","somatic","somatic","somatic","somatic","somatic","somatic","cognitive","cognitive","cognitive","cognitive","cognitive","cognitive")
df$gen_prob<-c("depression","anxiety","attention","aggression","headache","nausea","vomiting","dizziness","fatigue","sleep+","sleep-","drowsiness","troublesleeping","s-recall","l-recall","seq-memory","work-memory","proc-speed","cardsorting")
df$boottest<-boot.pc2$bootRatiosSignificant.j[,plsc.t09.pc2[[1]]<.05]
df_melt<-melt(df)

label_data <- as.data.frame(sort(df$gen_prob))
colnames(label_data)<-"gen_prob"
label_data$id<-seq(1,19)
# calculate the ANGLE of the labels
number_of_bar <- nrow(label_data)
angle <-  90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data<-merge(label_data,df_melt[df_melt$variable=="Y_3",c("gen_prob","value","boottest")])
label_data$value<-label_data$value*-1
# calculate the alignment of labels: right or left
# If I am on the left part of the plot, my labels have currently an angle < -90
label_data$hjust<-ifelse( angle < -90, 1, 0)

# flip angle BY to make them readable
label_data$angle<-ifelse(angle < -90, angle+180, angle)

p<-ggplot(data=df_melt[df_melt$variable=="Y_3",],aes(x=gen_prob,y=value,fill=gen_prob)) +
  # p<-ggplot(data=df_melt[df_melt$variable=="Y_4",],aes(x=gen_prob,y=value,fill=gen_prob)) +
  geom_bar(position = 'dodge', stat='identity')+scale_fill_manual(values=c25)+
  coord_polar(start=0)+ylim(-0.7,1)+
  geom_bar(stat = "identity", aes(x = gen_prob, y = 0.01), fill = "white") +
  geom_bar(stat = "identity", aes(x = gen_prob, y = -0.04), fill = "white") +
  geom_bar(stat = "identity", aes(x = gen_prob, y = -0.03), fill = "turquoise") +
  theme(legend.position="none",
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        axis.ticks=element_blank(),
        panel.background = element_rect(fill="white",color="white"))+
  geom_text(data=label_data, aes(x=id, y=.5, label=gen_prob, hjust=hjust,color=gen_prob), fontface="bold", size=10.5, angle= label_data$angle, inherit.aes = FALSE )+ scale_colour_manual(values=c25)+
geom_text(data=label_data, aes(x=id, y=ifelse(value<0,value+.05,value-.08), label=ifelse(boottest==TRUE,"*",""), hjust=hjust), color="black",size=7.8, angle= label_data$angle, inherit.aes = FALSE )
p




####Correlations with ADHD####
#PC1 xLF
tmp1<-cor.test(all.data.pc1[all.data.pc1$set==1,c("X_1")]*-1,as.numeric(as.character(all.data.pc1[all.data.pc1$set==1,c("ksads_14_853_p")])))#
tmp11<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="X_1")
tmp1<-cor.test(all.data.pc1[all.data.pc1$set==1,c("X_2")],as.numeric(as.character(all.data.pc1[all.data.pc1$set==1,c("ksads_14_853_p")])))#
tmp21<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="X_2")
tmp1<-cor.test(all.data.pc1[all.data.pc1$set==1,c("X_3")],as.numeric(as.character(all.data.pc1[all.data.pc1$set==1,c("ksads_14_853_p")])))#
tmp31<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="X_3")
# tmp1<-cor.test(all.data.pc1[all.data.pc1$set==1,c("X_4")],as.numeric(as.character(all.data.pc1[all.data.pc1$set==1,c("ksads_14_853_p")])))#
tmp1<-cor.test(all.data.pc1[all.data.pc1$set==1,c("X_5")],as.numeric(as.character(all.data.pc1[all.data.pc1$set==1,c("ksads_14_853_p")])))#
tmp41<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="X_5")
# tmp1<-cor.test(all.data.pc1[all.data.pc1$set==1,c("X_6")],as.numeric(as.character(all.data.pc1[all.data.pc1$set==1,c("ksads_14_853_p")])))#
# tmp51<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="X_6")
tmp1<-cor.test(all.data.pc1[all.data.pc1$set==1,c("X_7")]*-1,as.numeric(as.character(all.data.pc1[all.data.pc1$set==1,c("ksads_14_853_p")])))#
tmp61<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="X_7")
tmp1<-cor.test(all.data.pc1[all.data.pc1$set==1,c("X_8")],as.numeric(as.character(all.data.pc1[all.data.pc1$set==1,c("ksads_14_853_p")])))#
tmp71<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="X_8")
tmp1<-cor.test(all.data.pc1[all.data.pc1$set==1,c("X_9")],as.numeric(as.character(all.data.pc1[all.data.pc1$set==1,c("ksads_14_853_p")])))#
tmp81<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="X_9")
tmp1<-cor.test(all.data.pc1[all.data.pc1$set==1,c("X_10")],as.numeric(as.character(all.data.pc1[all.data.pc1$set==1,c("ksads_14_853_p")])))#
tmp91<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="X_10")
tmp1<-cor.test(all.data.pc1[all.data.pc1$set==1,c("X_11")],as.numeric(as.character(all.data.pc1[all.data.pc1$set==1,c("ksads_14_853_p")])))#
tmp101<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="X_11")
tmp1<-cor.test(all.data.pc1[all.data.pc1$set==1,c("X_12")],as.numeric(as.character(all.data.pc1[all.data.pc1$set==1,c("ksads_14_853_p")])))#
tmp111<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="X_12")
tmp1<-cor.test(all.data.pc1[all.data.pc1$set==1,c("X_13")],as.numeric(as.character(all.data.pc1[all.data.pc1$set==1,c("ksads_14_853_p")])))#
tmp121<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="X_13")
tmp1<-cor.test(all.data.pc1[all.data.pc1$set==1,c("X_14")],as.numeric(as.character(all.data.pc1[all.data.pc1$set==1,c("ksads_14_853_p")])))#
tmp131<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="X_14")
tmp1<-cor.test(all.data.pc1[all.data.pc1$set==1,c("X_15")],as.numeric(as.character(all.data.pc1[all.data.pc1$set==1,c("ksads_14_853_p")])))#
tmp141<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="X_15")
tmp1<-cor.test(all.data.pc1[all.data.pc1$set==1,c("X_16")],as.numeric(as.character(all.data.pc1[all.data.pc1$set==1,c("ksads_14_853_p")])))#
tmp151<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="X_16")
tmp1<-cor.test(all.data.pc1[all.data.pc1$set==1,c("X_17")],as.numeric(as.character(all.data.pc1[all.data.pc1$set==1,c("ksads_14_853_p")])))#
tmp161<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="X_17")
tmp1<-cor.test(all.data.pc1[all.data.pc1$set==1,c("X_18")],as.numeric(as.character(all.data.pc1[all.data.pc1$set==1,c("ksads_14_853_p")])))#
tmp171<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="X_18")
# tmp1<-cor.test(all.data.pc1[all.data.pc1$set==1,c("X_19")],as.numeric(as.character(all.data.pc1[all.data.pc1$set==1,c("ksads_14_853_p")])))#
adhd_cors_pc1_xLF<-rbind(tmp11,tmp21,tmp31,tmp41,tmp61,tmp71,tmp81,tmp91,tmp101,tmp111,tmp121,tmp131,tmp141,tmp151,tmp161,tmp171)


#PC1 yLF
tmp1<-cor.test(all.data.pc1[all.data.pc1$set==1,c("Y_1")]*-1,as.numeric(as.character(all.data.pc1[all.data.pc1$set==1,c("ksads_14_853_p")])))#
tmp11<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="Y_1")
tmp1<-cor.test(all.data.pc1[all.data.pc1$set==1,c("Y_2")],as.numeric(as.character(all.data.pc1[all.data.pc1$set==1,c("ksads_14_853_p")])))#
tmp21<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="Y_2")
tmp1<-cor.test(all.data.pc1[all.data.pc1$set==1,c("Y_3")],as.numeric(as.character(all.data.pc1[all.data.pc1$set==1,c("ksads_14_853_p")])))#
tmp31<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="Y_3")
# tmp1<-cor.test(all.data.pc1[all.data.pc1$set==1,c("X_4")],as.numeric(as.character(all.data.pc1[all.data.pc1$set==1,c("ksads_14_853_p")])))#
tmp1<-cor.test(all.data.pc1[all.data.pc1$set==1,c("Y_5")],as.numeric(as.character(all.data.pc1[all.data.pc1$set==1,c("ksads_14_853_p")])))#
tmp41<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="Y_5")
# tmp1<-cor.test(all.data.pc1[all.data.pc1$set==1,c("Y_6")],as.numeric(as.character(all.data.pc1[all.data.pc1$set==1,c("ksads_14_853_p")])))#
# tmp51<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="Y_6")
tmp1<-cor.test(all.data.pc1[all.data.pc1$set==1,c("Y_7")]*-1,as.numeric(as.character(all.data.pc1[all.data.pc1$set==1,c("ksads_14_853_p")])))#
tmp61<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="Y_7")
tmp1<-cor.test(all.data.pc1[all.data.pc1$set==1,c("Y_8")],as.numeric(as.character(all.data.pc1[all.data.pc1$set==1,c("ksads_14_853_p")])))#
tmp71<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="Y_8")
tmp1<-cor.test(all.data.pc1[all.data.pc1$set==1,c("Y_9")],as.numeric(as.character(all.data.pc1[all.data.pc1$set==1,c("ksads_14_853_p")])))#
tmp81<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="Y_9")
tmp1<-cor.test(all.data.pc1[all.data.pc1$set==1,c("Y_10")],as.numeric(as.character(all.data.pc1[all.data.pc1$set==1,c("ksads_14_853_p")])))#
tmp91<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="Y_10")
tmp1<-cor.test(all.data.pc1[all.data.pc1$set==1,c("Y_11")],as.numeric(as.character(all.data.pc1[all.data.pc1$set==1,c("ksads_14_853_p")])))#
tmp101<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="Y_11")
tmp1<-cor.test(all.data.pc1[all.data.pc1$set==1,c("Y_12")],as.numeric(as.character(all.data.pc1[all.data.pc1$set==1,c("ksads_14_853_p")])))#
tmp111<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="Y_12")
tmp1<-cor.test(all.data.pc1[all.data.pc1$set==1,c("Y_13")],as.numeric(as.character(all.data.pc1[all.data.pc1$set==1,c("ksads_14_853_p")])))#
tmp121<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="Y_13")
tmp1<-cor.test(all.data.pc1[all.data.pc1$set==1,c("Y_14")],as.numeric(as.character(all.data.pc1[all.data.pc1$set==1,c("ksads_14_853_p")])))#
tmp131<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="Y_14")
tmp1<-cor.test(all.data.pc1[all.data.pc1$set==1,c("Y_15")],as.numeric(as.character(all.data.pc1[all.data.pc1$set==1,c("ksads_14_853_p")])))#
tmp141<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="Y_15")
tmp1<-cor.test(all.data.pc1[all.data.pc1$set==1,c("Y_16")],as.numeric(as.character(all.data.pc1[all.data.pc1$set==1,c("ksads_14_853_p")])))#
tmp151<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="Y_16")
tmp1<-cor.test(all.data.pc1[all.data.pc1$set==1,c("Y_17")],as.numeric(as.character(all.data.pc1[all.data.pc1$set==1,c("ksads_14_853_p")])))#
tmp161<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="Y_17")
tmp1<-cor.test(all.data.pc1[all.data.pc1$set==1,c("Y_18")],as.numeric(as.character(all.data.pc1[all.data.pc1$set==1,c("ksads_14_853_p")])))#
tmp171<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="Y_18")
# tmp1<-cor.test(all.data.pc1[all.data.pc1$set==1,c("X_19")],as.numeric(as.character(all.data.pc1[all.data.pc1$set==1,c("ksads_14_853_p")])))#
adhd_cors_pc1_yLF<-rbind(tmp11,tmp21,tmp31,tmp41,tmp61,tmp71,tmp81,tmp91,tmp101,tmp111,tmp121,tmp131,tmp141,tmp151,tmp161,tmp171)






#PC2 xLF
# tmp1<-cor.test(all.data.pc2[all.data.pc2$set==1,c("X_1")],as.numeric(as.character(all.data.pc2[all.data.pc2$set==1,c("ksads_14_853_p")])))#
# tmp11<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="X_1")
tmp1<-cor.test(all.data.pc2[all.data.pc2$set==1,c("X_2")],as.numeric(as.character(all.data.pc2[all.data.pc2$set==1,c("ksads_14_853_p")])))#
tmp21<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="X_2")
tmp1<-cor.test(all.data.pc2[all.data.pc2$set==1,c("X_3")],as.numeric(as.character(all.data.pc2[all.data.pc2$set==1,c("ksads_14_853_p")])))#
tmp31<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="X_3")
# tmp1<-cor.test(all.data.pc1[all.data.pc2$set==1,c("X_4")],as.numeric(as.character(all.data.pc2[all.data.pc1$set==1,c("ksads_14_853_p")])))#
tmp1<-cor.test(all.data.pc2[all.data.pc2$set==1,c("X_5")],as.numeric(as.character(all.data.pc2[all.data.pc2$set==1,c("ksads_14_853_p")])))#
tmp41<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="X_5")
# tmp1<-cor.test(all.data.pc1[all.data.pc2$set==1,c("X_6")],as.numeric(as.character(all.data.pc2[all.data.pc1$set==1,c("ksads_14_853_p")])))#
# tmp51<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="X_6")
# tmp1<-cor.test(all.data.pc2[all.data.pc2$set==1,c("X_7")],as.numeric(as.character(all.data.pc2[all.data.pc2$set==1,c("ksads_14_853_p")])))#
# tmp61<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="X_7")
tmp1<-cor.test(all.data.pc2[all.data.pc2$set==1,c("X_8")],as.numeric(as.character(all.data.pc2[all.data.pc2$set==1,c("ksads_14_853_p")])))#
tmp71<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="X_8")
tmp1<-cor.test(all.data.pc2[all.data.pc2$set==1,c("X_9")],as.numeric(as.character(all.data.pc2[all.data.pc2$set==1,c("ksads_14_853_p")])))#
tmp81<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="X_9")
tmp1<-cor.test(all.data.pc2[all.data.pc2$set==1,c("X_10")],as.numeric(as.character(all.data.pc2[all.data.pc2$set==1,c("ksads_14_853_p")])))#
tmp91<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="X_10")
tmp1<-cor.test(all.data.pc2[all.data.pc2$set==1,c("X_11")],as.numeric(as.character(all.data.pc2[all.data.pc2$set==1,c("ksads_14_853_p")])))#
tmp101<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="X_11")
# tmp1<-cor.test(all.data.pc2[all.data.pc2$set==1,c("X_12")],as.numeric(as.character(all.data.pc2[all.data.pc2$set==1,c("ksads_14_853_p")])))#
# tmp111<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="X_12")
# tmp1<-cor.test(all.data.pc2[all.data.pc2$set==1,c("X_13")],as.numeric(as.character(all.data.pc2[all.data.pc2$set==1,c("ksads_14_853_p")])))#
# tmp121<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="X_13")
# tmp1<-cor.test(all.data.pc2[all.data.pc2$set==1,c("X_14")],as.numeric(as.character(all.data.pc2[all.data.pc2$set==1,c("ksads_14_853_p")])))#
# tmp131<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="X_14")
# tmp1<-cor.test(all.data.pc2[all.data.pc2$set==1,c("X_15")],as.numeric(as.character(all.data.pc2[all.data.pc2$set==1,c("ksads_14_853_p")])))#
# tmp141<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="X_15")
# tmp1<-cor.test(all.data.pc2[all.data.pc2$set==1,c("X_16")],as.numeric(as.character(all.data.pc2[all.data.pc2$set==1,c("ksads_14_853_p")])))#
# tmp151<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="X_16")
# tmp1<-cor.test(all.data.pc2[all.data.pc2$set==1,c("X_17")],as.numeric(as.character(all.data.pc2[all.data.pc2$set==1,c("ksads_14_853_p")])))#
# tmp161<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="X_17")
# tmp1<-cor.test(all.data.pc2[all.data.pc2$set==1,c("X_18")],as.numeric(as.character(all.data.pc2[all.data.pc2$set==1,c("ksads_14_853_p")])))#
# tmp171<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="X_18")
# tmp1<-cor.test(all.data.pc2[all.data.pc2$set==1,c("X_19")],as.numeric(as.character(all.data.pc2[all.data.pc2$set==1,c("ksads_14_853_p")])))#
adhd_cors_pc2_xLF<-rbind(tmp21,tmp31,tmp41,tmp71,tmp81,tmp91,tmp101)



#PC2 yLF
# tmp1<-cor.test(all.data.pc2[all.data.pc2$set==1,c("Y_1")],as.numeric(as.character(all.data.pc2[all.data.pc2$set==1,c("ksads_14_853_p")])))#
# tmp11<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="Y_1")
tmp1<-cor.test(all.data.pc2[all.data.pc2$set==1,c("Y_2")],as.numeric(as.character(all.data.pc2[all.data.pc2$set==1,c("ksads_14_853_p")])))#
tmp21<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="Y_2")
tmp1<-cor.test(all.data.pc2[all.data.pc2$set==1,c("Y_3")],as.numeric(as.character(all.data.pc2[all.data.pc2$set==1,c("ksads_14_853_p")])))#
tmp31<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="Y_3")
# tmp1<-cor.test(all.data.pc1[all.data.pc1$set==1,c("X_4")],as.numeric(as.character(all.data.pc1[all.data.pc1$set==1,c("ksads_14_853_p")])))#
tmp1<-cor.test(all.data.pc2[all.data.pc2$set==1,c("Y_5")],as.numeric(as.character(all.data.pc2[all.data.pc2$set==1,c("ksads_14_853_p")])))#
tmp41<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="Y_5")
# tmp1<-cor.test(all.data.pc1[all.data.pc1$set==1,c("X_6")],as.numeric(as.character(all.data.pc1[all.data.pc1$set==1,c("ksads_14_853_p")])))#
# tmp51<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="X_6")
# tmp1<-cor.test(all.data.pc2[all.data.pc2$set==1,c("Y_7")],as.numeric(as.character(all.data.pc2[all.data.pc2$set==1,c("ksads_14_853_p")])))#
# tmp61<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="Y_7")
tmp1<-cor.test(all.data.pc2[all.data.pc2$set==1,c("Y_8")],as.numeric(as.character(all.data.pc2[all.data.pc2$set==1,c("ksads_14_853_p")])))#
tmp71<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="Y_8")
tmp1<-cor.test(all.data.pc2[all.data.pc2$set==1,c("Y_9")],as.numeric(as.character(all.data.pc2[all.data.pc2$set==1,c("ksads_14_853_p")])))#
tmp81<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="Y_9")
tmp1<-cor.test(all.data.pc2[all.data.pc2$set==1,c("Y_10")],as.numeric(as.character(all.data.pc2[all.data.pc2$set==1,c("ksads_14_853_p")])))#
tmp91<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="Y_10")
tmp1<-cor.test(all.data.pc2[all.data.pc2$set==1,c("Y_11")],as.numeric(as.character(all.data.pc2[all.data.pc2$set==1,c("ksads_14_853_p")])))#
tmp101<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="Y_11")
# tmp1<-cor.test(all.data.pc2[all.data.pc2$set==1,c("Y_12")],as.numeric(as.character(all.data.pc2[all.data.pc2$set==1,c("ksads_14_853_p")])))#
# tmp111<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="Y_12")
# tmp1<-cor.test(all.data.pc2[all.data.pc2$set==1,c("Y_13")],as.numeric(as.character(all.data.pc2[all.data.pc2$set==1,c("ksads_14_853_p")])))#
# tmp121<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="Y_13")
# tmp1<-cor.test(all.data.pc2[all.data.pc2$set==1,c("Y_14")],as.numeric(as.character(all.data.pc2[all.data.pc2$set==1,c("ksads_14_853_p")])))#
# tmp131<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="Y_14")
# tmp1<-cor.test(all.data.pc2[all.data.pc2$set==1,c("Y_15")],as.numeric(as.character(all.data.pc2[all.data.pc2$set==1,c("ksads_14_853_p")])))#
# tmp141<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="Y_15")
# tmp1<-cor.test(all.data.pc2[all.data.pc2$set==1,c("Y_16")],as.numeric(as.character(all.data.pc2[all.data.pc2$set==1,c("ksads_14_853_p")])))#
# tmp151<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="Y_16")
# tmp1<-cor.test(all.data.pc2[all.data.pc2$set==1,c("Y_17")],as.numeric(as.character(all.data.pc2[all.data.pc2$set==1,c("ksads_14_853_p")])))#
# tmp161<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="Y_17")
# tmp1<-cor.test(all.data.pc2[all.data.pc2$set==1,c("Y_18")],as.numeric(as.character(all.data.pc2[all.data.pc2$set==1,c("ksads_14_853_p")])))#
# tmp171<-data.frame(cor=tmp1$estimate, pval=tmp1$p.value,lf="Y_18")
# tmp1<-cor.test(all.data.pc1[all.data.pc1$set==1,c("Y_19")],as.numeric(as.character(all.data.pc1[all.data.pc1$set==1,c("ksads_14_853_p")])))#
adhd_cors_pc2_yLF<-rbind(tmp21,tmp31,tmp41,tmp71,tmp81,tmp91,tmp101)


adhd_cors_pc1_xLF
adhd_cors_pc1_yLF
adhd_cors_pc2_xLF
adhd_cors_pc2_yLF


####Figure 4 % Overlap Analyses####
####UNIVARIATE BRAIN COMPS VS MULTIVARIATE LFS (line graphs)####
#Store the subs with clinical-level scores
clinsubs<-all.data.pc1[all.data.pc1$ksads_14_853_p==1,c("subjectkey")]

PC1.scores.clincomps<-PC1.scores
PC1.scores.clincomps$adhd<-ifelse(PC1.scores.clincomps$subjectkey%in%clinsubs,1,0)


stat.test <- PC1.scores.clincomps %>%
  group_by(conns) %>%
  t_test(score ~ adhd) %>%
  add_significance()
u_conns<-as.data.frame(stat.test[,c("conns","p")])
u_conns$sig<-ifelse(u_conns$p<0.05&u_conns$p>0,1,0) 

u_conns_pc1<-u_conns

u_conns_mat<-reshape_conn_mat(u_conns,lbls,"p")

u_conns_mat_pc1<-u_conns_mat

u_conns_mat_pc1.sig<-ifelse(u_conns_mat_pc1<0.05&u_conns_mat_pc1>0,1,0)


loadings.pc1<-as.data.frame(plsc.t09.pc1[[4]]$TExPosition.Data$pdq$p)
colnames(loadings.pc1)<-paste("X",substr(colnames(as.data.frame(plsc.t09.pc1[[4]]$TExPosition.Data$pdq$p)),2,5),sep="_")
loadings.pc1$conn<-rownames(loadings.pc1)


xLF1.pc1<-reshape_conn_mat(loadings.pc1[,c("X_1","conn")],lbls,"X_1")
xLF2.pc1<-reshape_conn_mat(loadings.pc1[,c("X_2","conn")],lbls,"X_2")
xLF3.pc1<-reshape_conn_mat(loadings.pc1[,c("X_3","conn")],lbls,"X_3")
xLF4.pc1<-reshape_conn_mat(loadings.pc1[,c("X_4","conn")],lbls,"X_4")
xLF5.pc1<-reshape_conn_mat(loadings.pc1[,c("X_5","conn")],lbls,"X_5")
xLF6.pc1<-reshape_conn_mat(loadings.pc1[,c("X_6","conn")],lbls,"X_6")
xLF7.pc1<-reshape_conn_mat(loadings.pc1[,c("X_7","conn")],lbls,"X_7")
xLF8.pc1<-reshape_conn_mat(loadings.pc1[,c("X_8","conn")],lbls,"X_8")
xLF9.pc1<-reshape_conn_mat(loadings.pc1[,c("X_9","conn")],lbls,"X_9")
xLF10.pc1<-reshape_conn_mat(loadings.pc1[,c("X_10","conn")],lbls,"X_10")
xLF11.pc1<-reshape_conn_mat(loadings.pc1[,c("X_11","conn")],lbls,"X_11")
xLF12.pc1<-reshape_conn_mat(loadings.pc1[,c("X_12","conn")],lbls,"X_12")
xLF13.pc1<-reshape_conn_mat(loadings.pc1[,c("X_13","conn")],lbls,"X_13")
xLF14.pc1<-reshape_conn_mat(loadings.pc1[,c("X_14","conn")],lbls,"X_14")
xLF15.pc1<-reshape_conn_mat(loadings.pc1[,c("X_15","conn")],lbls,"X_15")
xLF16.pc1<-reshape_conn_mat(loadings.pc1[,c("X_16","conn")],lbls,"X_16")
xLF17.pc1<-reshape_conn_mat(loadings.pc1[,c("X_17","conn")],lbls,"X_17")
xLF18.pc1<-reshape_conn_mat(loadings.pc1[,c("X_18","conn")],lbls,"X_18")
xLF19.pc1<-reshape_conn_mat(loadings.pc1[,c("X_19","conn")],lbls,"X_19")

# #FLIP THE FIRST THREE SO THEY'RE CONSISTENT WITH THE SCATTER PLOTS
# xLF1<-xLF1*-1
# xLF3<-xLF3*-1
# xLF4<-xLF4*-1

#NOW ONLY STORE THE SIGNIFICANT ONES


xLF1.pc1.bootsig<-as.data.frame(boot.pc1$bootRatiosSignificant.i[,1])
colnames(xLF1.pc1.bootsig)<-"sig"
xLF1.pc1.bootsig$conn<-rownames(xLF1.pc1.bootsig)
xLF1.pc1.bootsig$signum<-ifelse(xLF1.pc1.bootsig$sig==TRUE,1,0)
xLF1.pc1.bootsig.mat<-reshape_conn_mat(xLF1.pc1.bootsig[,c("signum","conn")],lbls,"signum")#This works
xLF1.pc1.loadings.sig<-xLF1.pc1*xLF1.pc1.bootsig.mat
View(xLF1.pc1.loadings.sig)


xLF2.pc1.bootsig<-as.data.frame(boot.pc1$bootRatiosSignificant.i[,2])
colnames(xLF2.pc1.bootsig)<-"sig"
xLF2.pc1.bootsig$conn<-rownames(xLF2.pc1.bootsig)
xLF2.pc1.bootsig$signum<-ifelse(xLF2.pc1.bootsig$sig==TRUE,1,0)
xLF2.pc1.bootsig.mat<-reshape_conn_mat(xLF2.pc1.bootsig[,c("signum","conn")],lbls,"signum")#This works
xLF2.pc1.loadings.sig<-xLF2.pc1*xLF2.pc1.bootsig.mat

xLF3.pc1.bootsig<-as.data.frame(boot.pc1$bootRatiosSignificant.i[,3])
colnames(xLF3.pc1.bootsig)<-"sig"
xLF3.pc1.bootsig$conn<-rownames(xLF3.pc1.bootsig)
xLF3.pc1.bootsig$signum<-ifelse(xLF3.pc1.bootsig$sig==TRUE,1,0)
xLF3.pc1.bootsig.mat<-reshape_conn_mat(xLF3.pc1.bootsig[,c("signum","conn")],lbls,"signum")#This works
xLF3.pc1.loadings.sig<-xLF3.pc1*xLF3.pc1.bootsig.mat

xLF4.pc1.bootsig<-as.data.frame(boot.pc1$bootRatiosSignificant.i[,4])
colnames(xLF4.pc1.bootsig)<-"sig"
xLF4.pc1.bootsig$conn<-rownames(xLF4.pc1.bootsig)
xLF4.pc1.bootsig$signum<-ifelse(xLF4.pc1.bootsig$sig==TRUE,1,0)
xLF4.pc1.bootsig.mat<-reshape_conn_mat(xLF4.pc1.bootsig[,c("signum","conn")],lbls,"signum")#This works
xLF4.pc1.loadings.sig<-xLF4.pc1*xLF4.pc1.bootsig.mat

xLF5.pc1.bootsig<-as.data.frame(boot.pc1$bootRatiosSignificant.i[,5])
colnames(xLF5.pc1.bootsig)<-"sig"
xLF5.pc1.bootsig$conn<-rownames(xLF5.pc1.bootsig)
xLF5.pc1.bootsig$signum<-ifelse(xLF5.pc1.bootsig$sig==TRUE,1,0)
xLF5.pc1.bootsig.mat<-reshape_conn_mat(xLF5.pc1.bootsig[,c("signum","conn")],lbls,"signum")#This works
xLF5.pc1.loadings.sig<-xLF5.pc1*xLF5.pc1.bootsig.mat

xLF6.pc1.bootsig<-as.data.frame(boot.pc1$bootRatiosSignificant.i[,6])
colnames(xLF6.pc1.bootsig)<-"sig"
xLF6.pc1.bootsig$conn<-rownames(xLF6.pc1.bootsig)
xLF6.pc1.bootsig$signum<-ifelse(xLF6.pc1.bootsig$sig==TRUE,1,0)
xLF6.pc1.bootsig.mat<-reshape_conn_mat(xLF6.pc1.bootsig[,c("signum","conn")],lbls,"signum")#This works
xLF6.pc1.loadings.sig<-xLF6.pc1*xLF6.pc1.bootsig.mat


xLF7.pc1.bootsig<-as.data.frame(boot.pc1$bootRatiosSignificant.i[,7])
colnames(xLF7.pc1.bootsig)<-"sig"
xLF7.pc1.bootsig$conn<-rownames(xLF7.pc1.bootsig)
xLF7.pc1.bootsig$signum<-ifelse(xLF7.pc1.bootsig$sig==TRUE,1,0)
xLF7.pc1.bootsig.mat<-reshape_conn_mat(xLF7.pc1.bootsig[,c("signum","conn")],lbls,"signum")#This works
xLF7.pc1.loadings.sig<-xLF7.pc1*xLF7.pc1.bootsig.mat

xLF8.pc1.bootsig<-as.data.frame(boot.pc1$bootRatiosSignificant.i[,8])
colnames(xLF8.pc1.bootsig)<-"sig"
xLF8.pc1.bootsig$conn<-rownames(xLF8.pc1.bootsig)
xLF8.pc1.bootsig$signum<-ifelse(xLF8.pc1.bootsig$sig==TRUE,1,0)
xLF8.pc1.bootsig.mat<-reshape_conn_mat(xLF8.pc1.bootsig[,c("signum","conn")],lbls,"signum")#This works
xLF8.pc1.loadings.sig<-xLF8.pc1*xLF8.pc1.bootsig.mat

xLF9.pc1.bootsig<-as.data.frame(boot.pc1$bootRatiosSignificant.i[,9])
colnames(xLF9.pc1.bootsig)<-"sig"
xLF9.pc1.bootsig$conn<-rownames(xLF9.pc1.bootsig)
xLF9.pc1.bootsig$signum<-ifelse(xLF9.pc1.bootsig$sig==TRUE,1,0)
xLF9.pc1.bootsig.mat<-reshape_conn_mat(xLF9.pc1.bootsig[,c("signum","conn")],lbls,"signum")#This works
xLF9.pc1.loadings.sig<-xLF9.pc1*xLF9.pc1.bootsig.mat

xLF10.pc1.bootsig<-as.data.frame(boot.pc1$bootRatiosSignificant.i[,10])
colnames(xLF10.pc1.bootsig)<-"sig"
xLF10.pc1.bootsig$conn<-rownames(xLF10.pc1.bootsig)
xLF10.pc1.bootsig$signum<-ifelse(xLF10.pc1.bootsig$sig==TRUE,1,0)
xLF10.pc1.bootsig.mat<-reshape_conn_mat(xLF10.pc1.bootsig[,c("signum","conn")],lbls,"signum")#This works
xLF10.pc1.loadings.sig<-xLF10.pc1*xLF10.pc1.bootsig.mat

xLF11.pc1.bootsig<-as.data.frame(boot.pc1$bootRatiosSignificant.i[,11])
colnames(xLF11.pc1.bootsig)<-"sig"
xLF11.pc1.bootsig$conn<-rownames(xLF11.pc1.bootsig)
xLF11.pc1.bootsig$signum<-ifelse(xLF11.pc1.bootsig$sig==TRUE,1,0)
xLF11.pc1.bootsig.mat<-reshape_conn_mat(xLF11.pc1.bootsig[,c("signum","conn")],lbls,"signum")#This works
xLF11.pc1.loadings.sig<-xLF11.pc1*xLF11.pc1.bootsig.mat

xLF12.pc1.bootsig<-as.data.frame(boot.pc1$bootRatiosSignificant.i[,12])
colnames(xLF12.pc1.bootsig)<-"sig"
xLF12.pc1.bootsig$conn<-rownames(xLF12.pc1.bootsig)
xLF12.pc1.bootsig$signum<-ifelse(xLF12.pc1.bootsig$sig==TRUE,1,0)
xLF12.pc1.bootsig.mat<-reshape_conn_mat(xLF12.pc1.bootsig[,c("signum","conn")],lbls,"signum")#This works
xLF12.pc1.loadings.sig<-xLF12.pc1*xLF12.pc1.bootsig.mat

xLF13.pc1.bootsig<-as.data.frame(boot.pc1$bootRatiosSignificant.i[,13])
colnames(xLF13.pc1.bootsig)<-"sig"
xLF13.pc1.bootsig$conn<-rownames(xLF13.pc1.bootsig)
xLF13.pc1.bootsig$signum<-ifelse(xLF13.pc1.bootsig$sig==TRUE,1,0)
xLF13.pc1.bootsig.mat<-reshape_conn_mat(xLF13.pc1.bootsig[,c("signum","conn")],lbls,"signum")#This works
xLF13.pc1.loadings.sig<-xLF13.pc1*xLF13.pc1.bootsig.mat

xLF14.pc1.bootsig<-as.data.frame(boot.pc1$bootRatiosSignificant.i[,14])
colnames(xLF14.pc1.bootsig)<-"sig"
xLF14.pc1.bootsig$conn<-rownames(xLF14.pc1.bootsig)
xLF14.pc1.bootsig$signum<-ifelse(xLF14.pc1.bootsig$sig==TRUE,1,0)
xLF14.pc1.bootsig.mat<-reshape_conn_mat(xLF14.pc1.bootsig[,c("signum","conn")],lbls,"signum")#This works
xLF14.pc1.loadings.sig<-xLF14.pc1*xLF14.pc1.bootsig.mat

xLF15.pc1.bootsig<-as.data.frame(boot.pc1$bootRatiosSignificant.i[,15])
colnames(xLF15.pc1.bootsig)<-"sig"
xLF15.pc1.bootsig$conn<-rownames(xLF15.pc1.bootsig)
xLF15.pc1.bootsig$signum<-ifelse(xLF15.pc1.bootsig$sig==TRUE,1,0)
xLF15.pc1.bootsig.mat<-reshape_conn_mat(xLF15.pc1.bootsig[,c("signum","conn")],lbls,"signum")#This works
xLF15.pc1.loadings.sig<-xLF15.pc1*xLF15.pc1.bootsig.mat


xLF16.pc1.bootsig<-as.data.frame(boot.pc1$bootRatiosSignificant.i[,16])
colnames(xLF16.pc1.bootsig)<-"sig"
xLF16.pc1.bootsig$conn<-rownames(xLF16.pc1.bootsig)
xLF16.pc1.bootsig$signum<-ifelse(xLF16.pc1.bootsig$sig==TRUE,1,0)
xLF16.pc1.bootsig.mat<-reshape_conn_mat(xLF16.pc1.bootsig[,c("signum","conn")],lbls,"signum")#This works
xLF16.pc1.loadings.sig<-xLF16.pc1*xLF16.pc1.bootsig.mat

xLF17.pc1.bootsig<-as.data.frame(boot.pc1$bootRatiosSignificant.i[,17])
colnames(xLF17.pc1.bootsig)<-"sig"
xLF17.pc1.bootsig$conn<-rownames(xLF17.pc1.bootsig)
xLF17.pc1.bootsig$signum<-ifelse(xLF17.pc1.bootsig$sig==TRUE,1,0)
xLF17.pc1.bootsig.mat<-reshape_conn_mat(xLF17.pc1.bootsig[,c("signum","conn")],lbls,"signum")#This works
xLF17.pc1.loadings.sig<-xLF17.pc1*xLF17.pc1.bootsig.mat

xLF18.pc1.bootsig<-as.data.frame(boot.pc1$bootRatiosSignificant.i[,18])
colnames(xLF18.pc1.bootsig)<-"sig"
xLF18.pc1.bootsig$conn<-rownames(xLF18.pc1.bootsig)
xLF18.pc1.bootsig$signum<-ifelse(xLF18.pc1.bootsig$sig==TRUE,1,0)
xLF18.pc1.bootsig.mat<-reshape_conn_mat(xLF18.pc1.bootsig[,c("signum","conn")],lbls,"signum")#This works
xLF18.pc1.loadings.sig<-xLF18.pc1*xLF18.pc1.bootsig.mat

xLF19.pc1.bootsig<-as.data.frame(boot.pc1$bootRatiosSignificant.i[,19])
colnames(xLF19.pc1.bootsig)<-"sig"
xLF19.pc1.bootsig$conn<-rownames(xLF19.pc1.bootsig)
xLF19.pc1.bootsig$signum<-ifelse(xLF19.pc1.bootsig$sig==TRUE,1,0)
xLF19.pc1.bootsig.mat<-reshape_conn_mat(xLF19.pc1.bootsig[,c("signum","conn")],lbls,"signum")#This works
xLF19.pc1.loadings.sig<-xLF19.pc1*xLF19.pc1.bootsig.mat



sim_mat1<-matrix(data = 0,nrow = dim(lbls)[1],ncol = dim(lbls)[1])
sim_mat2<-matrix(data = 0,nrow = dim(lbls)[1],ncol = dim(lbls)[1])
sim_mat3<-matrix(data = 0,nrow = dim(lbls)[1],ncol = dim(lbls)[1])
sim_mat4<-matrix(data = 0,nrow = dim(lbls)[1],ncol = dim(lbls)[1])
sim_mat5<-matrix(data = 0,nrow = dim(lbls)[1],ncol = dim(lbls)[1])
sim_mat6<-matrix(data = 0,nrow = dim(lbls)[1],ncol = dim(lbls)[1])
sim_mat7<-matrix(data = 0,nrow = dim(lbls)[1],ncol = dim(lbls)[1])
sim_mat8<-matrix(data = 0,nrow = dim(lbls)[1],ncol = dim(lbls)[1])
sim_mat9<-matrix(data = 0,nrow = dim(lbls)[1],ncol = dim(lbls)[1])
sim_mat10<-matrix(data = 0,nrow = dim(lbls)[1],ncol = dim(lbls)[1])
sim_mat11<-matrix(data = 0,nrow = dim(lbls)[1],ncol = dim(lbls)[1])
sim_mat12<-matrix(data = 0,nrow = dim(lbls)[1],ncol = dim(lbls)[1])
sim_mat13<-matrix(data = 0,nrow = dim(lbls)[1],ncol = dim(lbls)[1])
sim_mat14<-matrix(data = 0,nrow = dim(lbls)[1],ncol = dim(lbls)[1])
sim_mat15<-matrix(data = 0,nrow = dim(lbls)[1],ncol = dim(lbls)[1])
sim_mat16<-matrix(data = 0,nrow = dim(lbls)[1],ncol = dim(lbls)[1])
sim_mat17<-matrix(data = 0,nrow = dim(lbls)[1],ncol = dim(lbls)[1])
sim_mat18<-matrix(data = 0,nrow = dim(lbls)[1],ncol = dim(lbls)[1])
sim_mat19<-matrix(data = 0,nrow = dim(lbls)[1],ncol = dim(lbls)[1])

# 

tmp_conns<-xLF1.pc1.bootsig[xLF1.pc1.bootsig$signum==1,c("conn")]
commonsig<-sum(u_conns[as.character(u_conns$conns)%in%tmp_conns,c("sig")])
perc_sim1<-(commonsig/(sum(u_conns$sig)+length(tmp_conns)-commonsig))*100

tmp_conns<-xLF2.pc1.bootsig[xLF2.pc1.bootsig$signum==1,c("conn")]
commonsig<-sum(u_conns[as.character(u_conns$conns)%in%tmp_conns,c("sig")])
perc_sim2<-(commonsig/(sum(u_conns$sig)+length(tmp_conns)-commonsig))*100

tmp_conns<-xLF3.pc1.bootsig[xLF3.pc1.bootsig$signum==1,c("conn")]
commonsig<-sum(u_conns[as.character(u_conns$conns)%in%tmp_conns,c("sig")])
perc_sim3<-(commonsig/(sum(u_conns$sig)+length(tmp_conns)-commonsig))*100

tmp_conns<-xLF4.pc1.bootsig[xLF4.pc1.bootsig$signum==1,c("conn")]
commonsig<-sum(u_conns[as.character(u_conns$conns)%in%tmp_conns,c("sig")])
perc_sim4<-(commonsig/(sum(u_conns$sig)+length(tmp_conns)-commonsig))*100

tmp_conns<-xLF5.pc1.bootsig[xLF5.pc1.bootsig$signum==1,c("conn")]
commonsig<-sum(u_conns[as.character(u_conns$conns)%in%tmp_conns,c("sig")])
perc_sim5<-(commonsig/(sum(u_conns$sig)+length(tmp_conns)-commonsig))*100

tmp_conns<-xLF6.pc1.bootsig[xLF6.pc1.bootsig$signum==1,c("conn")]
commonsig<-sum(u_conns[as.character(u_conns$conns)%in%tmp_conns,c("sig")])
perc_sim6<-(commonsig/(sum(u_conns$sig)+length(tmp_conns)-commonsig))*100

tmp_conns<-xLF7.pc1.bootsig[xLF7.pc1.bootsig$signum==1,c("conn")]
commonsig<-sum(u_conns[as.character(u_conns$conns)%in%tmp_conns,c("sig")])
perc_sim7<-(commonsig/(sum(u_conns$sig)+length(tmp_conns)-commonsig))*100

tmp_conns<-xLF8.pc1.bootsig[xLF8.pc1.bootsig$signum==1,c("conn")]
commonsig<-sum(u_conns[as.character(u_conns$conns)%in%tmp_conns,c("sig")])
perc_sim8<-(commonsig/(sum(u_conns$sig)+length(tmp_conns)-commonsig))*100

tmp_conns<-xLF9.pc1.bootsig[xLF9.pc1.bootsig$signum==1,c("conn")]
commonsig<-sum(u_conns[as.character(u_conns$conns)%in%tmp_conns,c("sig")])
perc_sim9<-(commonsig/(sum(u_conns$sig)+length(tmp_conns)-commonsig))*100

tmp_conns<-xLF10.pc1.bootsig[xLF10.pc1.bootsig$signum==1,c("conn")]
commonsig<-sum(u_conns[as.character(u_conns$conns)%in%tmp_conns,c("sig")])
perc_sim10<-(commonsig/(sum(u_conns$sig)+length(tmp_conns)-commonsig))*100

tmp_conns<-xLF11.pc1.bootsig[xLF11.pc1.bootsig$signum==1,c("conn")]
commonsig<-sum(u_conns[as.character(u_conns$conns)%in%tmp_conns,c("sig")])
perc_sim11<-(commonsig/(sum(u_conns$sig)+length(tmp_conns)-commonsig))*100

tmp_conns<-xLF12.pc1.bootsig[xLF12.pc1.bootsig$signum==1,c("conn")]
commonsig<-sum(u_conns[as.character(u_conns$conns)%in%tmp_conns,c("sig")])
perc_sim12<-(commonsig/(sum(u_conns$sig)+length(tmp_conns)-commonsig))*100

tmp_conns<-xLF13.pc1.bootsig[xLF13.pc1.bootsig$signum==1,c("conn")]
commonsig<-sum(u_conns[as.character(u_conns$conns)%in%tmp_conns,c("sig")])
perc_sim13<-(commonsig/(sum(u_conns$sig)+length(tmp_conns)-commonsig))*100

tmp_conns<-xLF14.pc1.bootsig[xLF14.pc1.bootsig$signum==1,c("conn")]
commonsig<-sum(u_conns[as.character(u_conns$conns)%in%tmp_conns,c("sig")])
perc_sim14<-(commonsig/(sum(u_conns$sig)+length(tmp_conns)-commonsig))*100

tmp_conns<-xLF15.pc1.bootsig[xLF15.pc1.bootsig$signum==1,c("conn")]
commonsig<-sum(u_conns[as.character(u_conns$conns)%in%tmp_conns,c("sig")])
perc_sim15<-(commonsig/(sum(u_conns$sig)+length(tmp_conns)-commonsig))*100

tmp_conns<-xLF16.pc1.bootsig[xLF16.pc1.bootsig$signum==1,c("conn")]
commonsig<-sum(u_conns[as.character(u_conns$conns)%in%tmp_conns,c("sig")])
perc_sim16<-(commonsig/(sum(u_conns$sig)+length(tmp_conns)-commonsig))*100

tmp_conns<-xLF17.pc1.bootsig[xLF17.pc1.bootsig$signum==1,c("conn")]
commonsig<-sum(u_conns[as.character(u_conns$conns)%in%tmp_conns,c("sig")])
perc_sim17<-(commonsig/(sum(u_conns$sig)+length(tmp_conns)-commonsig))*100

tmp_conns<-xLF18.pc1.bootsig[xLF18.pc1.bootsig$signum==1,c("conn")]
commonsig<-sum(u_conns[as.character(u_conns$conns)%in%tmp_conns,c("sig")])
perc_sim18<-(commonsig/(sum(u_conns$sig)+length(tmp_conns)-commonsig))*100

tmp_conns<-xLF19.pc1.bootsig[xLF19.pc1.bootsig$signum==1,c("conn")]
commonsig<-sum(u_conns[as.character(u_conns$conns)%in%tmp_conns,c("sig")])
perc_sim19<-(commonsig/(sum(u_conns$sig)+length(tmp_conns)-commonsig))*100


perc_sims<-as.data.frame(rbind(perc_sim1,perc_sim2,perc_sim3,perc_sim4,perc_sim5,perc_sim6,perc_sim7,perc_sim8,perc_sim9,perc_sim10,perc_sim11,perc_sim12,perc_sim13,perc_sim14,perc_sim15,perc_sim16,perc_sim17,perc_sim18,perc_sim19))
colnames(perc_sims)<-"perc_overlap"
perc_sims$lf<-seq(1:19)
perc_sims$pc<-1
perc_sims$sig<-ifelse(plsc.t09.pc1[[1]]<0.05,1,0)

#PC2

clinsubs<-all.data.pc2[all.data.pc2$ksads_14_853_p==1,c("subjectkey")]

PC2.scores.clincomps<-PC2.scores
PC2.scores.clincomps$adhd<-ifelse(PC2.scores.clincomps$subjectkey%in%clinsubs,1,0)


stat.test <- PC2.scores.clincomps %>%
  group_by(conns) %>%
  t_test(score ~ adhd) %>%
  add_significance()
u_conns<-as.data.frame(stat.test[,c("conns","p")])
u_conns$sig<-ifelse(u_conns$p<0.05&u_conns$p>0,1,0)


u_conns_pc2<-u_conns

u_conns_mat<-reshape_conn_mat(u_conns,lbls,"p")

u_conns_mat_pc2<-u_conns_mat

u_conns_mat_pc2.sig<-ifelse(u_conns_mat_pc2<0.05&u_conns_mat_pc2>0,1,0)


loadings.pc2<-as.data.frame(plsc.t09.pc2[[4]]$TExPosition.Data$pdq$p)
colnames(loadings.pc2)<-paste("X",substr(colnames(as.data.frame(plsc.t09.pc2[[4]]$TExPosition.Data$pdq$p)),2,5),sep="_")
loadings.pc2$conn<-rownames(loadings.pc2)

xLF1.pc2<-reshape_conn_mat(loadings.pc2[,c("X_1","conn")],lbls,"X_1")#This works
xLF2.pc2<-reshape_conn_mat(loadings.pc2[,c("X_2","conn")],lbls,"X_2")#This works
xLF3.pc2<-reshape_conn_mat(loadings.pc2[,c("X_3","conn")],lbls,"X_3")#This works
xLF4.pc2<-reshape_conn_mat(loadings.pc2[,c("X_4","conn")],lbls,"X_4")#This works
xLF5.pc2<-reshape_conn_mat(loadings.pc2[,c("X_5","conn")],lbls,"X_5")#This works
xLF6.pc2<-reshape_conn_mat(loadings.pc2[,c("X_6","conn")],lbls,"X_6")#This works
xLF7.pc2<-reshape_conn_mat(loadings.pc2[,c("X_7","conn")],lbls,"X_7")#This works
xLF8.pc2<-reshape_conn_mat(loadings.pc2[,c("X_8","conn")],lbls,"X_8")#This works
xLF9.pc2<-reshape_conn_mat(loadings.pc2[,c("X_9","conn")],lbls,"X_9")#This works
xLF10.pc2<-reshape_conn_mat(loadings.pc2[,c("X_10","conn")],lbls,"X_10")#This works
xLF11.pc2<-reshape_conn_mat(loadings.pc2[,c("X_11","conn")],lbls,"X_11")#This works
xLF12.pc2<-reshape_conn_mat(loadings.pc2[,c("X_12","conn")],lbls,"X_12")#This works
xLF13.pc2<-reshape_conn_mat(loadings.pc2[,c("X_13","conn")],lbls,"X_13")#This works
xLF14.pc2<-reshape_conn_mat(loadings.pc2[,c("X_14","conn")],lbls,"X_14")#This works
xLF15.pc2<-reshape_conn_mat(loadings.pc2[,c("X_15","conn")],lbls,"X_15")#This works
xLF16.pc2<-reshape_conn_mat(loadings.pc2[,c("X_16","conn")],lbls,"X_16")#This works
xLF17.pc2<-reshape_conn_mat(loadings.pc2[,c("X_17","conn")],lbls,"X_17")#This works
xLF18.pc2<-reshape_conn_mat(loadings.pc2[,c("X_18","conn")],lbls,"X_18")#This works
xLF19.pc2<-reshape_conn_mat(loadings.pc2[,c("X_19","conn")],lbls,"X_19")#This works

# #FLIP THE FIRST THREE SO THEY'RE CONSISTENT WITH THE SCATTER PLOTS
# xLF1<-xLF1*-1
# xLF3<-xLF3*-1
# xLF4<-xLF4*-1

#NOW ONLY STORE THE SIGNIFICANT ONES



xLF1.pc2.bootsig<-as.data.frame(boot.pc2$bootRatiosSignificant.i[,1])
colnames(xLF1.pc2.bootsig)<-"sig"
xLF1.pc2.bootsig$conn<-rownames(xLF1.pc2.bootsig)
xLF1.pc2.bootsig$signum<-ifelse(xLF1.pc2.bootsig$sig==TRUE,1,0)
xLF1.pc2.bootsig.mat<-reshape_conn_mat(xLF1.pc2.bootsig[,c("signum","conn")],lbls,"signum")#This works
xLF1.pc2.loadings.sig<-xLF1.pc2*xLF1.pc2.bootsig.mat

xLF2.pc2.bootsig<-as.data.frame(boot.pc2$bootRatiosSignificant.i[,2])
colnames(xLF2.pc2.bootsig)<-"sig"
xLF2.pc2.bootsig$conn<-rownames(xLF2.pc2.bootsig)
xLF2.pc2.bootsig$signum<-ifelse(xLF2.pc2.bootsig$sig==TRUE,1,0)
xLF2.pc2.bootsig.mat<-reshape_conn_mat(xLF2.pc2.bootsig[,c("signum","conn")],lbls,"signum")#This works
xLF2.pc2.loadings.sig<-xLF2.pc2*xLF2.pc2.bootsig.mat

xLF3.pc2.bootsig<-as.data.frame(boot.pc2$bootRatiosSignificant.i[,3])
colnames(xLF3.pc2.bootsig)<-"sig"
xLF3.pc2.bootsig$conn<-rownames(xLF3.pc2.bootsig)
xLF3.pc2.bootsig$signum<-ifelse(xLF3.pc2.bootsig$sig==TRUE,1,0)
xLF3.pc2.bootsig.mat<-reshape_conn_mat(xLF3.pc2.bootsig[,c("signum","conn")],lbls,"signum")#This works
xLF3.pc2.loadings.sig<-xLF3.pc2*xLF3.pc2.bootsig.mat

xLF4.pc2.bootsig<-as.data.frame(boot.pc2$bootRatiosSignificant.i[,4])
colnames(xLF4.pc2.bootsig)<-"sig"
xLF4.pc2.bootsig$conn<-rownames(xLF4.pc2.bootsig)
xLF4.pc2.bootsig$signum<-ifelse(xLF4.pc2.bootsig$sig==TRUE,1,0)
xLF4.pc2.bootsig.mat<-reshape_conn_mat(xLF4.pc2.bootsig[,c("signum","conn")],lbls,"signum")#This works
xLF4.pc2.loadings.sig<-xLF4.pc2*xLF4.pc2.bootsig.mat

xLF5.pc2.bootsig<-as.data.frame(boot.pc2$bootRatiosSignificant.i[,5])
colnames(xLF5.pc2.bootsig)<-"sig"
xLF5.pc2.bootsig$conn<-rownames(xLF5.pc2.bootsig)
xLF5.pc2.bootsig$signum<-ifelse(xLF5.pc2.bootsig$sig==TRUE,1,0)
xLF5.pc2.bootsig.mat<-reshape_conn_mat(xLF5.pc2.bootsig[,c("signum","conn")],lbls,"signum")#This works
xLF5.pc2.loadings.sig<-xLF5.pc2*xLF5.pc2.bootsig.mat

xLF6.pc2.bootsig<-as.data.frame(boot.pc2$bootRatiosSignificant.i[,6])
colnames(xLF6.pc2.bootsig)<-"sig"
xLF6.pc2.bootsig$conn<-rownames(xLF6.pc2.bootsig)
xLF6.pc2.bootsig$signum<-ifelse(xLF6.pc2.bootsig$sig==TRUE,1,0)
xLF6.pc2.bootsig.mat<-reshape_conn_mat(xLF6.pc2.bootsig[,c("signum","conn")],lbls,"signum")#This works
xLF6.pc2.loadings.sig<-xLF6.pc2*xLF6.pc2.bootsig.mat


xLF7.pc2.bootsig<-as.data.frame(boot.pc2$bootRatiosSignificant.i[,7])
colnames(xLF7.pc2.bootsig)<-"sig"
xLF7.pc2.bootsig$conn<-rownames(xLF7.pc2.bootsig)
xLF7.pc2.bootsig$signum<-ifelse(xLF7.pc2.bootsig$sig==TRUE,1,0)
xLF7.pc2.bootsig.mat<-reshape_conn_mat(xLF7.pc2.bootsig[,c("signum","conn")],lbls,"signum")#This works
xLF7.pc2.loadings.sig<-xLF7.pc2*xLF7.pc2.bootsig.mat

xLF8.pc2.bootsig<-as.data.frame(boot.pc2$bootRatiosSignificant.i[,8])
colnames(xLF8.pc2.bootsig)<-"sig"
xLF8.pc2.bootsig$conn<-rownames(xLF8.pc2.bootsig)
xLF8.pc2.bootsig$signum<-ifelse(xLF8.pc2.bootsig$sig==TRUE,1,0)
xLF8.pc2.bootsig.mat<-reshape_conn_mat(xLF8.pc2.bootsig[,c("signum","conn")],lbls,"signum")#This works
xLF8.pc2.loadings.sig<-xLF8.pc2*xLF8.pc2.bootsig.mat

xLF9.pc2.bootsig<-as.data.frame(boot.pc2$bootRatiosSignificant.i[,9])
colnames(xLF9.pc2.bootsig)<-"sig"
xLF9.pc2.bootsig$conn<-rownames(xLF9.pc2.bootsig)
xLF9.pc2.bootsig$signum<-ifelse(xLF9.pc2.bootsig$sig==TRUE,1,0)
xLF9.pc2.bootsig.mat<-reshape_conn_mat(xLF9.pc2.bootsig[,c("signum","conn")],lbls,"signum")#This works
xLF9.pc2.loadings.sig<-xLF9.pc2*xLF9.pc2.bootsig.mat

xLF10.pc2.bootsig<-as.data.frame(boot.pc2$bootRatiosSignificant.i[,10])
colnames(xLF10.pc2.bootsig)<-"sig"
xLF10.pc2.bootsig$conn<-rownames(xLF10.pc2.bootsig)
xLF10.pc2.bootsig$signum<-ifelse(xLF10.pc2.bootsig$sig==TRUE,1,0)
xLF10.pc2.bootsig.mat<-reshape_conn_mat(xLF10.pc2.bootsig[,c("signum","conn")],lbls,"signum")#This works
xLF10.pc2.loadings.sig<-xLF10.pc2*xLF10.pc2.bootsig.mat

xLF11.pc2.bootsig<-as.data.frame(boot.pc2$bootRatiosSignificant.i[,11])
colnames(xLF11.pc2.bootsig)<-"sig"
xLF11.pc2.bootsig$conn<-rownames(xLF11.pc2.bootsig)
xLF11.pc2.bootsig$signum<-ifelse(xLF11.pc2.bootsig$sig==TRUE,1,0)
xLF11.pc2.bootsig.mat<-reshape_conn_mat(xLF11.pc2.bootsig[,c("signum","conn")],lbls,"signum")#This works
xLF11.pc2.loadings.sig<-xLF11.pc2*xLF11.pc2.bootsig.mat

xLF12.pc2.bootsig<-as.data.frame(boot.pc2$bootRatiosSignificant.i[,12])
colnames(xLF12.pc2.bootsig)<-"sig"
xLF12.pc2.bootsig$conn<-rownames(xLF12.pc2.bootsig)
xLF12.pc2.bootsig$signum<-ifelse(xLF12.pc2.bootsig$sig==TRUE,1,0)
xLF12.pc2.bootsig.mat<-reshape_conn_mat(xLF12.pc2.bootsig[,c("signum","conn")],lbls,"signum")#This works
xLF12.pc2.loadings.sig<-xLF12.pc2*xLF12.pc2.bootsig.mat

xLF13.pc2.bootsig<-as.data.frame(boot.pc2$bootRatiosSignificant.i[,13])
colnames(xLF13.pc2.bootsig)<-"sig"
xLF13.pc2.bootsig$conn<-rownames(xLF13.pc2.bootsig)
xLF13.pc2.bootsig$signum<-ifelse(xLF13.pc2.bootsig$sig==TRUE,1,0)
xLF13.pc2.bootsig.mat<-reshape_conn_mat(xLF13.pc2.bootsig[,c("signum","conn")],lbls,"signum")#This works
xLF13.pc2.loadings.sig<-xLF13.pc2*xLF13.pc2.bootsig.mat

xLF14.pc2.bootsig<-as.data.frame(boot.pc2$bootRatiosSignificant.i[,14])
colnames(xLF14.pc2.bootsig)<-"sig"
xLF14.pc2.bootsig$conn<-rownames(xLF14.pc2.bootsig)
xLF14.pc2.bootsig$signum<-ifelse(xLF14.pc2.bootsig$sig==TRUE,1,0)
xLF14.pc2.bootsig.mat<-reshape_conn_mat(xLF14.pc2.bootsig[,c("signum","conn")],lbls,"signum")#This works
xLF14.pc2.loadings.sig<-xLF14.pc2*xLF14.pc2.bootsig.mat

xLF15.pc2.bootsig<-as.data.frame(boot.pc2$bootRatiosSignificant.i[,15])
colnames(xLF15.pc2.bootsig)<-"sig"
xLF15.pc2.bootsig$conn<-rownames(xLF15.pc2.bootsig)
xLF15.pc2.bootsig$signum<-ifelse(xLF15.pc2.bootsig$sig==TRUE,1,0)
xLF15.pc2.bootsig.mat<-reshape_conn_mat(xLF15.pc2.bootsig[,c("signum","conn")],lbls,"signum")#This works
xLF15.pc2.loadings.sig<-xLF15.pc2*xLF15.pc2.bootsig.mat

xLF16.pc2.bootsig<-as.data.frame(boot.pc2$bootRatiosSignificant.i[,16])
colnames(xLF16.pc2.bootsig)<-"sig"
xLF16.pc2.bootsig$conn<-rownames(xLF16.pc2.bootsig)
xLF16.pc2.bootsig$signum<-ifelse(xLF16.pc2.bootsig$sig==TRUE,1,0)
xLF16.pc2.bootsig.mat<-reshape_conn_mat(xLF16.pc2.bootsig[,c("signum","conn")],lbls,"signum")#This works
xLF16.pc2.loadings.sig<-xLF16.pc2*xLF16.pc2.bootsig.mat

xLF17.pc2.bootsig<-as.data.frame(boot.pc2$bootRatiosSignificant.i[,17])
colnames(xLF17.pc2.bootsig)<-"sig"
xLF17.pc2.bootsig$conn<-rownames(xLF17.pc2.bootsig)
xLF17.pc2.bootsig$signum<-ifelse(xLF17.pc2.bootsig$sig==TRUE,1,0)
xLF17.pc2.bootsig.mat<-reshape_conn_mat(xLF17.pc2.bootsig[,c("signum","conn")],lbls,"signum")#This works
xLF17.pc2.loadings.sig<-xLF17.pc2*xLF17.pc2.bootsig.mat

xLF18.pc2.bootsig<-as.data.frame(boot.pc2$bootRatiosSignificant.i[,18])
colnames(xLF18.pc2.bootsig)<-"sig"
xLF18.pc2.bootsig$conn<-rownames(xLF18.pc2.bootsig)
xLF18.pc2.bootsig$signum<-ifelse(xLF18.pc2.bootsig$sig==TRUE,1,0)
xLF18.pc2.bootsig.mat<-reshape_conn_mat(xLF18.pc2.bootsig[,c("signum","conn")],lbls,"signum")#This works
xLF18.pc2.loadings.sig<-xLF18.pc2*xLF18.pc2.bootsig.mat

xLF19.pc2.bootsig<-as.data.frame(boot.pc2$bootRatiosSignificant.i[,19])
colnames(xLF19.pc2.bootsig)<-"sig"
xLF19.pc2.bootsig$conn<-rownames(xLF19.pc2.bootsig)
xLF19.pc2.bootsig$signum<-ifelse(xLF19.pc2.bootsig$sig==TRUE,1,0)
xLF19.pc2.bootsig.mat<-reshape_conn_mat(xLF19.pc2.bootsig[,c("signum","conn")],lbls,"signum")#This works
xLF19.pc2.loadings.sig<-xLF19.pc2*xLF19.pc2.bootsig.mat



sim_mat1<-matrix(data = 0,nrow = dim(lbls)[1],ncol = dim(lbls)[1])
sim_mat2<-matrix(data = 0,nrow = dim(lbls)[1],ncol = dim(lbls)[1])
sim_mat3<-matrix(data = 0,nrow = dim(lbls)[1],ncol = dim(lbls)[1])
sim_mat4<-matrix(data = 0,nrow = dim(lbls)[1],ncol = dim(lbls)[1])
sim_mat5<-matrix(data = 0,nrow = dim(lbls)[1],ncol = dim(lbls)[1])
sim_mat6<-matrix(data = 0,nrow = dim(lbls)[1],ncol = dim(lbls)[1])
sim_mat7<-matrix(data = 0,nrow = dim(lbls)[1],ncol = dim(lbls)[1])
sim_mat8<-matrix(data = 0,nrow = dim(lbls)[1],ncol = dim(lbls)[1])
sim_mat9<-matrix(data = 0,nrow = dim(lbls)[1],ncol = dim(lbls)[1])
sim_mat10<-matrix(data = 0,nrow = dim(lbls)[1],ncol = dim(lbls)[1])
sim_mat11<-matrix(data = 0,nrow = dim(lbls)[1],ncol = dim(lbls)[1])
sim_mat12<-matrix(data = 0,nrow = dim(lbls)[1],ncol = dim(lbls)[1])
sim_mat13<-matrix(data = 0,nrow = dim(lbls)[1],ncol = dim(lbls)[1])
sim_mat14<-matrix(data = 0,nrow = dim(lbls)[1],ncol = dim(lbls)[1])
sim_mat15<-matrix(data = 0,nrow = dim(lbls)[1],ncol = dim(lbls)[1])
sim_mat16<-matrix(data = 0,nrow = dim(lbls)[1],ncol = dim(lbls)[1])
sim_mat17<-matrix(data = 0,nrow = dim(lbls)[1],ncol = dim(lbls)[1])
sim_mat18<-matrix(data = 0,nrow = dim(lbls)[1],ncol = dim(lbls)[1])
sim_mat19<-matrix(data = 0,nrow = dim(lbls)[1],ncol = dim(lbls)[1])



tmp_conns<-xLF1.pc2.bootsig[xLF1.pc2.bootsig$signum==1,c("conn")]
commonsig<-sum(u_conns[as.character(u_conns$conns)%in%tmp_conns,c("sig")])
perc_sim1<-(commonsig/(sum(u_conns$sig)+length(tmp_conns)-commonsig))*100

tmp_conns<-xLF2.pc2.bootsig[xLF2.pc2.bootsig$signum==1,c("conn")]
commonsig<-sum(u_conns[as.character(u_conns$conns)%in%tmp_conns,c("sig")])
perc_sim2<-(commonsig/(sum(u_conns$sig)+length(tmp_conns)-commonsig))*100

tmp_conns<-xLF3.pc2.bootsig[xLF3.pc2.bootsig$signum==1,c("conn")]
commonsig<-sum(u_conns[as.character(u_conns$conns)%in%tmp_conns,c("sig")])
perc_sim3<-(commonsig/(sum(u_conns$sig)+length(tmp_conns)-commonsig))*100

tmp_conns<-xLF4.pc2.bootsig[xLF4.pc2.bootsig$signum==1,c("conn")]
commonsig<-sum(u_conns[as.character(u_conns$conns)%in%tmp_conns,c("sig")])
perc_sim4<-(commonsig/(sum(u_conns$sig)+length(tmp_conns)-commonsig))*100

tmp_conns<-xLF5.pc2.bootsig[xLF5.pc2.bootsig$signum==1,c("conn")]
commonsig<-sum(u_conns[as.character(u_conns$conns)%in%tmp_conns,c("sig")])
perc_sim5<-(commonsig/(sum(u_conns$sig)+length(tmp_conns)-commonsig))*100

tmp_conns<-xLF6.pc2.bootsig[xLF6.pc2.bootsig$signum==1,c("conn")]
commonsig<-sum(u_conns[as.character(u_conns$conns)%in%tmp_conns,c("sig")])
perc_sim6<-(commonsig/(sum(u_conns$sig)+length(tmp_conns)-commonsig))*100

tmp_conns<-xLF7.pc2.bootsig[xLF7.pc2.bootsig$signum==1,c("conn")]
commonsig<-sum(u_conns[as.character(u_conns$conns)%in%tmp_conns,c("sig")])
perc_sim7<-(commonsig/(sum(u_conns$sig)+length(tmp_conns)-commonsig))*100

tmp_conns<-xLF8.pc2.bootsig[xLF8.pc2.bootsig$signum==1,c("conn")]
commonsig<-sum(u_conns[as.character(u_conns$conns)%in%tmp_conns,c("sig")])
perc_sim8<-(commonsig/(sum(u_conns$sig)+length(tmp_conns)-commonsig))*100

tmp_conns<-xLF9.pc2bootsig[xLF9.pc2.bootsig$signum==1,c("conn")]
commonsig<-sum(u_conns[as.character(u_conns$conns)%in%tmp_conns,c("sig")])
perc_sim9<-(commonsig/(sum(u_conns$sig)+length(tmp_conns)-commonsig))*100

tmp_conns<-xLF10.pc2.bootsig[xLF10.pc2.bootsig$signum==1,c("conn")]
commonsig<-sum(u_conns[as.character(u_conns$conns)%in%tmp_conns,c("sig")])
perc_sim10<-(commonsig/(sum(u_conns$sig)+length(tmp_conns)-commonsig))*100

tmp_conns<-xLF11.pc2.bootsig[xLF11.pc2.bootsig$signum==1,c("conn")]
commonsig<-sum(u_conns[as.character(u_conns$conns)%in%tmp_conns,c("sig")])
perc_sim11<-(commonsig/(sum(u_conns$sig)+length(tmp_conns)-commonsig))*100

tmp_conns<-xLF12.pc2.bootsig[xLF12.pc2.bootsig$signum==1,c("conn")]
commonsig<-sum(u_conns[as.character(u_conns$conns)%in%tmp_conns,c("sig")])
perc_sim12<-(commonsig/(sum(u_conns$sig)+length(tmp_conns)-commonsig))*100

tmp_conns<-xLF13.pc2.bootsig[xLF13.pc2.bootsig$signum==1,c("conn")]
commonsig<-sum(u_conns[as.character(u_conns$conns)%in%tmp_conns,c("sig")])
perc_sim13<-(commonsig/(sum(u_conns$sig)+length(tmp_conns)-commonsig))*100

tmp_conns<-xLF14.pc2.bootsig[xLF14.pc2.bootsig$signum==1,c("conn")]
commonsig<-sum(u_conns[as.character(u_conns$conns)%in%tmp_conns,c("sig")])
perc_sim14<-(commonsig/(sum(u_conns$sig)+length(tmp_conns)-commonsig))*100

tmp_conns<-xLF15.pc2.bootsig[xLF15.pc2.bootsig$signum==1,c("conn")]
commonsig<-sum(u_conns[as.character(u_conns$conns)%in%tmp_conns,c("sig")])
perc_sim15<-(commonsig/(sum(u_conns$sig)+length(tmp_conns)-commonsig))*100

tmp_conns<-xLF16.pc2.bootsig[xLF16.pc2.bootsig$signum==1,c("conn")]
commonsig<-sum(u_conns[as.character(u_conns$conns)%in%tmp_conns,c("sig")])
perc_sim16<-(commonsig/(sum(u_conns$sig)+length(tmp_conns)-commonsig))*100

tmp_conns<-xLF17.pc2.bootsig[xLF17.pc2.bootsig$signum==1,c("conn")]
commonsig<-sum(u_conns[as.character(u_conns$conns)%in%tmp_conns,c("sig")])
perc_sim17<-(commonsig/(sum(u_conns$sig)+length(tmp_conns)-commonsig))*100

tmp_conns<-xLF18.pc2.bootsig[xLF18.pc2.bootsig$signum==1,c("conn")]
commonsig<-sum(u_conns[as.character(u_conns$conns)%in%tmp_conns,c("sig")])
perc_sim18<-(commonsig/(sum(u_conns$sig)+length(tmp_conns)-commonsig))*100

tmp_conns<-xLF19.pc2.bootsig[xLF19.pc2.bootsig$signum==1,c("conn")]
commonsig<-sum(u_conns[as.character(u_conns$conns)%in%tmp_conns,c("sig")])
perc_sim19<-(commonsig/(sum(u_conns$sig)+length(tmp_conns)-commonsig))*100


tmp_perc_sims<-as.data.frame(rbind(perc_sim1,perc_sim2,perc_sim3,perc_sim4,perc_sim5,perc_sim6,perc_sim7,perc_sim8,perc_sim9,perc_sim10,perc_sim11,perc_sim12,perc_sim13,perc_sim14,perc_sim15,perc_sim16,perc_sim17,perc_sim18,perc_sim19))
colnames(tmp_perc_sims)<-"perc_overlap"
tmp_perc_sims$lf<-seq(1:19)
tmp_perc_sims$pc<-2
tmp_perc_sims$sig<-ifelse(plsc.t09.pc2[[1]]<0.05,1,0)

perc_sims<-rbind(perc_sims,tmp_perc_sims)
perc_sims$perc_overlap<-ifelse(is.na(perc_sims$perc_overlap),0,perc_sims$perc_overlap)


p<-ggplot(data=perc_sims, aes(x=lf, y=perc_overlap,group=as.factor(pc))) +
  geom_line(size=1,aes(linetype=as.factor(pc)))+
  geom_point(size=2)+ylim(0,15)+
  labs(x="Multi-Tract Connectivity Feature",y="% Overlap")+scale_linetype_discrete(name="PLSc", labels=c("PLSc1","PLSc2"))+
  scale_x_continuous(breaks = seq(1,19))+
  theme_classic(base_size = 20)+
  theme(axis.title=element_text(face="bold"),
        axis.text=element_text(face="bold"))
p



####Figure 4 brain graphs + brain graphs for supp material####

u_conns_mat_pc1.bin<-ifelse(u_conns_mat_pc1.sig==0,0,1)
u_conns_mat_pc2.bin<-ifelse(u_conns_mat_pc2.sig==0,0,1)

#First binarize your multivar mats, you only care about which connections are identified
xLF1.pc1.loadings.sig.bin<-ifelse(xLF1.pc1.loadings.sig==0,0,1)
xLF2.pc1.loadings.sig.bin<-ifelse(xLF2.pc1.loadings.sig==0,0,1)
xLF3.pc1.loadings.sig.bin<-ifelse(xLF3.pc1.loadings.sig==0,0,1)
xLF4.pc1.loadings.sig.bin<-ifelse(xLF4.pc1.loadings.sig==0,0,1)
xLF5.pc1.loadings.sig.bin<-ifelse(xLF5.pc1.loadings.sig==0,0,1)
xLF6.pc1.loadings.sig.bin<-ifelse(xLF6.pc1.loadings.sig==0,0,1)
xLF7.pc1.loadings.sig.bin<-ifelse(xLF7.pc1.loadings.sig==0,0,1)
xLF8.pc1.loadings.sig.bin<-ifelse(xLF8.pc1.loadings.sig==0,0,1)
xLF9.pc1.loadings.sig.bin<-ifelse(xLF9.pc1.loadings.sig==0,0,1)
xLF10.pc1.loadings.sig.bin<-ifelse(xLF10.pc1.loadings.sig==0,0,1)
xLF11.pc1.loadings.sig.bin<-ifelse(xLF11.pc1.loadings.sig==0,0,1)
xLF12.pc1.loadings.sig.bin<-ifelse(xLF12.pc1.loadings.sig==0,0,1)
xLF13.pc1.loadings.sig.bin<-ifelse(xLF13.pc1.loadings.sig==0,0,1)
xLF14.pc1.loadings.sig.bin<-ifelse(xLF14.pc1.loadings.sig==0,0,1)
xLF15.pc1.loadings.sig.bin<-ifelse(xLF15.pc1.loadings.sig==0,0,1)
xLF16.pc1.loadings.sig.bin<-ifelse(xLF16.pc1.loadings.sig==0,0,1)
xLF17.pc1.loadings.sig.bin<-ifelse(xLF17.pc1.loadings.sig==0,0,1)
xLF18.pc1.loadings.sig.bin<-ifelse(xLF18.pc1.loadings.sig==0,0,1)
xLF19.pc1.loadings.sig.bin<-ifelse(xLF19.pc1.loadings.sig==0,0,1)


xLF3.pc2.loadings.sig.bin<-ifelse(xLF3.pc2.loadings.sig==0,0,1)

#Store conn mats
write.table(u_conns_mat_pc1.bin,file = '/Users/Guido/Desktop/Projects/PCHStudy/Figures/DraftOfFinalFigures/FinalDraft/Round3-Resubmission/COMMITBeforePCA/uconnspc1.txt',quote = F,col.names = F,row.names = F)
write.table(u_conns_mat_pc2.bin,file = '/Users/Guido/Desktop/Projects/PCHStudy/Figures/DraftOfFinalFigures/FinalDraft/Round3-Resubmission/COMMITBeforePCA/uconnspc2.txt',quote = F,col.names = F,row.names = F)

write.table(xLF1.pc1.loadings.sig.bin,file = '/Users/Guido/Desktop/Projects/PCHStudy/Figures/DraftOfFinalFigures/FinalDraft/Round3-Resubmission/COMMITBeforePCA/xLF1pc1loadingssig.txt',quote = F,col.names = F,row.names = F)
write.table(xLF2.pc1.loadings.sig.bin,file = '/Users/Guido/Desktop/Projects/PCHStudy/Figures/DraftOfFinalFigures/FinalDraft/Round3-Resubmission/COMMITBeforePCA/xLF2pc1loadingssig.txt',quote = F,col.names = F,row.names = F)
write.table(xLF3.pc1.loadings.sig.bin,file = '/Users/Guido/Desktop/Projects/PCHStudy/Figures/DraftOfFinalFigures/FinalDraft/Round3-Resubmission/COMMITBeforePCA/xLF3pc1loadingssig.txt',quote = F,col.names = F,row.names = F)
write.table(xLF4.pc1.loadings.sig.bin,file = '/Users/Guido/Desktop/Projects/PCHStudy/Figures/DraftOfFinalFigures/FinalDraft/Round3-Resubmission/COMMITBeforePCA/xLF4pc1loadingssig.txt',quote = F,col.names = F,row.names = F)
write.table(xLF5.pc1.loadings.sig.bin,file = '/Users/Guido/Desktop/Projects/PCHStudy/Figures/DraftOfFinalFigures/FinalDraft/Round3-Resubmission/COMMITBeforePCA/xLF5pc1loadingssig.txt',quote = F,col.names = F,row.names = F)
write.table(xLF6.pc1.loadings.sig.bin,file = '/Users/Guido/Desktop/Projects/PCHStudy/Figures/DraftOfFinalFigures/FinalDraft/Round3-Resubmission/COMMITBeforePCA/xLF6pc1loadingssig.txt',quote = F,col.names = F,row.names = F)
write.table(xLF7.pc1.loadings.sig.bin,file = '/Users/Guido/Desktop/Projects/PCHStudy/Figures/DraftOfFinalFigures/FinalDraft/Round3-Resubmission/COMMITBeforePCA/xLF7pc1loadingssig.txt',quote = F,col.names = F,row.names = F)
write.table(xLF8.pc1.loadings.sig.bin,file = '/Users/Guido/Desktop/Projects/PCHStudy/Figures/DraftOfFinalFigures/FinalDraft/Round3-Resubmission/COMMITBeforePCA/xLF8pc1loadingssig.txt',quote = F,col.names = F,row.names = F)
write.table(xLF9.pc1.loadings.sig.bin,file = '/Users/Guido/Desktop/Projects/PCHStudy/Figures/DraftOfFinalFigures/FinalDraft/Round3-Resubmission/COMMITBeforePCA/xLF9pc1loadingssig.txt',quote = F,col.names = F,row.names = F)
write.table(xLF10.pc1.loadings.sig.bin,file = '/Users/Guido/Desktop/Projects/PCHStudy/Figures/DraftOfFinalFigures/FinalDraft/Round3-Resubmission/COMMITBeforePCA/xLF10pc1loadingssig.txt',quote = F,col.names = F,row.names = F)
write.table(xLF11.pc1.loadings.sig.bin,file = '/Users/Guido/Desktop/Projects/PCHStudy/Figures/DraftOfFinalFigures/FinalDraft/Round3-Resubmission/COMMITBeforePCA/xLF11pc1loadingssig.txt',quote = F,col.names = F,row.names = F)
write.table(xLF12.pc1.loadings.sig.bin,file = '/Users/Guido/Desktop/Projects/PCHStudy/Figures/DraftOfFinalFigures/FinalDraft/Round3-Resubmission/COMMITBeforePCA/xLF12pc1loadingssig.txt',quote = F,col.names = F,row.names = F)
write.table(xLF13.pc1.loadings.sig.bin,file = '/Users/Guido/Desktop/Projects/PCHStudy/Figures/DraftOfFinalFigures/FinalDraft/Round3-Resubmission/COMMITBeforePCA/xLF13pc1loadingssig.txt',quote = F,col.names = F,row.names = F)
write.table(xLF14.pc1.loadings.sig.bin,file = '/Users/Guido/Desktop/Projects/PCHStudy/Figures/DraftOfFinalFigures/FinalDraft/Round3-Resubmission/COMMITBeforePCA/xLF14pc1loadingssig.txt',quote = F,col.names = F,row.names = F)
write.table(xLF15.pc1.loadings.sig.bin,file = '/Users/Guido/Desktop/Projects/PCHStudy/Figures/DraftOfFinalFigures/FinalDraft/Round3-Resubmission/COMMITBeforePCA/xLF15pc1loadingssig.txt',quote = F,col.names = F,row.names = F)
write.table(xLF16.pc1.loadings.sig.bin,file = '/Users/Guido/Desktop/Projects/PCHStudy/Figures/DraftOfFinalFigures/FinalDraft/Round3-Resubmission/COMMITBeforePCA/xLF16pc1loadingssig.txt',quote = F,col.names = F,row.names = F)
write.table(xLF17.pc1.loadings.sig.bin,file = '/Users/Guido/Desktop/Projects/PCHStudy/Figures/DraftOfFinalFigures/FinalDraft/Round3-Resubmission/COMMITBeforePCA/xLF17pc1loadingssig.txt',quote = F,col.names = F,row.names = F)
write.table(xLF18.pc1.loadings.sig.bin,file = '/Users/Guido/Desktop/Projects/PCHStudy/Figures/DraftOfFinalFigures/FinalDraft/Round3-Resubmission/COMMITBeforePCA/xLF18pc1loadingssig.txt',quote = F,col.names = F,row.names = F)
write.table(xLF19.pc1.loadings.sig.bin,file = '/Users/Guido/Desktop/Projects/PCHStudy/Figures/DraftOfFinalFigures/FinalDraft/Round3-Resubmission/COMMITBeforePCA/xLF19pc1loadingssig.txt',quote = F,col.names = F,row.names = F)

write.table(xLF3.pc2.loadings.sig.bin,file = '/Users/Guido/Desktop/Projects/PCHStudy/Figures/DraftOfFinalFigures/FinalDraft/Round3-Resubmission/COMMITBeforePCA/xLF3pc2loadingssig.txt',quote = F,col.names = F,row.names = F)

#Compute a matrix of connection frequency 
xConnSus.pc1.loadings.sig.bin<-xLF1.pc1.loadings.sig.bin+xLF2.pc1.loadings.sig.bin+xLF3.pc1.loadings.sig.bin+xLF4.pc1.loadings.sig.bin+xLF5.pc1.loadings.sig.bin+
  xLF6.pc1.loadings.sig.bin+xLF7.pc1.loadings.sig.bin+xLF8.pc1.loadings.sig.bin+xLF9.pc1.loadings.sig.bin+xLF10.pc1.loadings.sig.bin+xLF11.pc1.loadings.sig.bin+
  xLF12.pc1.loadings.sig.bin+xLF13.pc1.loadings.sig.bin+xLF14.pc1.loadings.sig.bin+xLF15.pc1.loadings.sig.bin+xLF16.pc1.loadings.sig.bin+
  xLF17.pc1.loadings.sig.bin+xLF18.pc1.loadings.sig.bin+xLF19.pc1.loadings.sig.bin



write.table(xConnSus.pc1.loadings.sig.bin,file = '/Users/Guido/Desktop/Projects/PCHStudy/Figures/DraftOfFinalFigures/FinalDraft/Round3-Resubmission/COMMITBeforePCA/xConnSusPC1.txt',quote = F,col.names = F,row.names = F)  

max(xConnSus.pc1.loadings.sig.bin)#4
sum(xConnSus.pc1.loadings.sig.bin==4)#4 but remember it's a symmetric matrix (so 2)
#So the max number of MTF that a single conn was repeated was 4 and 2 conns were among the most found

xConnSus.pc1.loadings.sig.bin.df<-as.data.frame(xConnSus.pc1.loadings.sig.bin)
colnames(xConnSus.pc1.loadings.sig.bin.df)[colSums(xConnSus.pc1.loadings.sig.bin.df==4)>0]#"V30" "V36" "V61" "V65" "V67" "V71" "V74"
rownames(xConnSus.pc1.loadings.sig.bin.df)[rowSums(xConnSus.pc1.loadings.sig.bin.df==4)>0]#"30" "36" "61" "65" "67" "71" "74"

colnames(xConnSus.pc1.loadings.sig.bin.df[61,])[xConnSus.pc1.loadings.sig.bin.df[61,]==4]#2018_2022 & 2018_2031
colnames(xConnSus.pc1.loadings.sig.bin.df[65,])[xConnSus.pc1.loadings.sig.bin.df[65,]==4]#61_65 (2018_2022) (rpOPER-rpostC)
colnames(xConnSus.pc1.loadings.sig.bin.df[74,])[xConnSus.pc1.loadings.sig.bin.df[74,]==4]#61_74 (2018_2031) (rpOPER-rSMAR)

###Alternative way of doing the same analysis as above, just to verify###
##Find out which pairs implicated them##

colnames(xLF1.pc1.bootsig)<-c("sig.X1","conn","signum.X1")
colnames(xLF2.pc1.bootsig)<-c("sig.X2","conn","signum.X2")
colnames(xLF3.pc1.bootsig)<-c("sig.X3","conn","signum.X3")
colnames(xLF4.pc1.bootsig)<-c("sig.X4","conn","signum.X4")
colnames(xLF5.pc1.bootsig)<-c("sig.X5","conn","signum.X5")
colnames(xLF6.pc1.bootsig)<-c("sig.X6","conn","signum.X6")
colnames(xLF7.pc1.bootsig)<-c("sig.X7","conn","signum.X7")
colnames(xLF8.pc1.bootsig)<-c("sig.X8","conn","signum.X8")
colnames(xLF9.pc1.bootsig)<-c("sig.X9","conn","signum.X9")
colnames(xLF10.pc1.bootsig)<-c("sig.X10","conn","signum.X10")
colnames(xLF11.pc1.bootsig)<-c("sig.X11","conn","signum.X11")
colnames(xLF12.pc1.bootsig)<-c("sig.X12","conn","signum.X12")
colnames(xLF13.pc1.bootsig)<-c("sig.X13","conn","signum.X13")
colnames(xLF14.pc1.bootsig)<-c("sig.X14","conn","signum.X14")
colnames(xLF15.pc1.bootsig)<-c("sig.X15","conn","signum.X15")
colnames(xLF16.pc1.bootsig)<-c("sig.X16","conn","signum.X16")
colnames(xLF17.pc1.bootsig)<-c("sig.X17","conn","signum.X17")
colnames(xLF18.pc1.bootsig)<-c("sig.X18","conn","signum.X18")
colnames(xLF19.pc1.bootsig)<-c("sig.X19","conn","signum.X19")

xConnSus.pc1.bootsig<-cbind(xLF1.pc1.bootsig,xLF2.pc1.bootsig,xLF3.pc1.bootsig,xLF4.pc1.bootsig,xLF5.pc1.bootsig,xLF6.pc1.bootsig,
                            xLF7.pc1.bootsig,xLF8.pc1.bootsig,xLF9.pc1.bootsig,xLF10.pc1.bootsig,xLF11.pc1.bootsig,xLF12.pc1.bootsig,
                            xLF13.pc1.bootsig,xLF14.pc1.bootsig,xLF15.pc1.bootsig,xLF16.pc1.bootsig,xLF17.pc1.bootsig,xLF18.pc1.bootsig,xLF19.pc1.bootsig)

xConnSus.pc1.bootsig$sums<-xConnSus.pc1.bootsig$signum.X1+xConnSus.pc1.bootsig$signum.X2+xConnSus.pc1.bootsig$signum.X3+xConnSus.pc1.bootsig$signum.X4+xConnSus.pc1.bootsig$signum.X5+xConnSus.pc1.bootsig$signum.X6+xConnSus.pc1.bootsig$signum.X7+
  xConnSus.pc1.bootsig$signum.X8+xConnSus.pc1.bootsig$signum.X9+xConnSus.pc1.bootsig$signum.X10+xConnSus.pc1.bootsig$signum.X11+xConnSus.pc1.bootsig$signum.X12+xConnSus.pc1.bootsig$signum.X13+xConnSus.pc1.bootsig$signum.X14+
  xConnSus.pc1.bootsig$signum.X15+xConnSus.pc1.bootsig$signum.X16+xConnSus.pc1.bootsig$signum.X17+xConnSus.pc1.bootsig$signum.X18+xConnSus.pc1.bootsig$signum.X19

xConnSus.pc1.bootsig[xConnSus.pc1.bootsig$sums==max(xConnSus.pc1.bootsig$sums),]

xConnSus.pc1.bootsig[xConnSus.pc1.bootsig$conn%in%c("2018_2022","2018_2031"),]
# 2018_2022 (rpOPER_rpostC): LF1, LF2, LF7, LF8
# 2018_2031 (rpOPER_rSMAR): LF1, LF8, LF16, LF18
tmp<-boot.pc1$bootRatiosSignificant.j[rownames(boot.pc1$bootRatiosSignificant.j)=="cbcl_scr_syn_attention_t",]
#ATTENTION: LF1, LF7, LF8, LF16

xConnSus.pc1.bootsig[xConnSus.pc1.bootsig$sums==3,]


sum(xConnSus.pc1.loadings.sig.bin==3)#26 but remember it's a symmetric matrix (so 13)
colnames(xConnSus.pc1.loadings.sig.bin.df)[colSums(xConnSus.pc1.loadings.sig.bin.df==3)>0]#"V11" "V12" "V15" "V16" "V20" "V21" "V24" "V27" "V34" "V36" "V37" "V40" "V42" "V45" "V46" "V47" "V50" "V54" "V61" "V70" "V71" "V72"
rownames(xConnSus.pc1.loadings.sig.bin.df)[rowSums(xConnSus.pc1.loadings.sig.bin.df==3)>0]#"11" "12" "15" "16" "20" "21" "24" "27" "34" "36" "37" "40" "42" "45" "46" "47" "50" "54" "61" "70" "71" "72"

colnames(xConnSus.pc1.loadings.sig.bin.df[11,])[xConnSus.pc1.loadings.sig.bin.df[11,]==3]#11_12
colnames(xConnSus.pc1.loadings.sig.bin.df[12,])[xConnSus.pc1.loadings.sig.bin.df[12,]==3]
colnames(xConnSus.pc1.loadings.sig.bin.df[15,])[xConnSus.pc1.loadings.sig.bin.df[15,]==3]#15_36
colnames(xConnSus.pc1.loadings.sig.bin.df[16,])[xConnSus.pc1.loadings.sig.bin.df[16,]==3]#16_71 (callosal)
colnames(xConnSus.pc1.loadings.sig.bin.df[20,])[xConnSus.pc1.loadings.sig.bin.df[20,]==3]#20_27
colnames(xConnSus.pc1.loadings.sig.bin.df[21,])[xConnSus.pc1.loadings.sig.bin.df[21,]==3]#21_40
colnames(xConnSus.pc1.loadings.sig.bin.df[24,])[xConnSus.pc1.loadings.sig.bin.df[24,]==3]#24_45
colnames(xConnSus.pc1.loadings.sig.bin.df[27,])[xConnSus.pc1.loadings.sig.bin.df[27,]==3]#27_42
colnames(xConnSus.pc1.loadings.sig.bin.df[34,])[xConnSus.pc1.loadings.sig.bin.df[34,]==3]#34_72 (callosal)
colnames(xConnSus.pc1.loadings.sig.bin.df[36,])[xConnSus.pc1.loadings.sig.bin.df[36,]==3]#
colnames(xConnSus.pc1.loadings.sig.bin.df[37,])[xConnSus.pc1.loadings.sig.bin.df[37,]==3]#37_54 (callosal)
colnames(xConnSus.pc1.loadings.sig.bin.df[40,])[xConnSus.pc1.loadings.sig.bin.df[40,]==3]#40_47 (callosal)
colnames(xConnSus.pc1.loadings.sig.bin.df[42,])[xConnSus.pc1.loadings.sig.bin.df[42,]==3]#
colnames(xConnSus.pc1.loadings.sig.bin.df[45,])[xConnSus.pc1.loadings.sig.bin.df[45,]==3]#
colnames(xConnSus.pc1.loadings.sig.bin.df[46,])[xConnSus.pc1.loadings.sig.bin.df[46,]==3]#46_61
colnames(xConnSus.pc1.loadings.sig.bin.df[47,])[xConnSus.pc1.loadings.sig.bin.df[47,]==3]#
colnames(xConnSus.pc1.loadings.sig.bin.df[50,])[xConnSus.pc1.loadings.sig.bin.df[50,]==3]#50_54
colnames(xConnSus.pc1.loadings.sig.bin.df[54,])[xConnSus.pc1.loadings.sig.bin.df[54,]==3]#
colnames(xConnSus.pc1.loadings.sig.bin.df[61,])[xConnSus.pc1.loadings.sig.bin.df[61,]==3]#
colnames(xConnSus.pc1.loadings.sig.bin.df[70,])[xConnSus.pc1.loadings.sig.bin.df[70,]==3]#70_71
colnames(xConnSus.pc1.loadings.sig.bin.df[71,])[xConnSus.pc1.loadings.sig.bin.df[71,]==3]#
colnames(xConnSus.pc1.loadings.sig.bin.df[72,])[xConnSus.pc1.loadings.sig.bin.df[72,]==3]#

#So, the most common conns are not callosal. There are some callosal conns that are less common across LFs.

rownames(boot.pc1$bootRatiosSignificant.i)[boot.pc1$bootRatiosSignificant.i[,1]==TRUE]
####How many conns are callosal in UPC1 UPC2 LF1 AND LF7?####

#Calculate how many callosal connections were found
#Note: the analyses below are not reported in the manuscript, they were just out of interest
#uPC1
u_conns_pc1$end1<-substr(u_conns_pc1$conn,1,regexpr("_",u_conns_pc1$conn)-1)
u_conns_pc1$end2<-substr(u_conns_pc1$conn,regexpr("_",u_conns_pc1$conn)+1,20)
u_conns_pc1$call<-ifelse(((u_conns_pc1$end1>1000&u_conns_pc1$end1<2000)&u_conns_pc1$end2>2000)|
                           ((u_conns_pc1$end1>9&u_conns_pc1$end1<27)&u_conns_pc1$end2>2000)|
                           ((u_conns_pc1$end1>9&u_conns_pc1$end1<27)&(u_conns_pc1$end2>48&u_conns_pc1$end2<59))|
                           ((u_conns_pc1$end1>48&u_conns_pc1$end1<59)&(u_conns_pc1$end2>1000&u_conns_pc1$end2<2000)),1,0)
#uPC2
u_conns_pc2$end1<-substr(u_conns_pc2$conn,1,regexpr("_",u_conns_pc2$conn)-1)
u_conns_pc2$end2<-substr(u_conns_pc2$conn,regexpr("_",u_conns_pc2$conn)+1,20)
u_conns_pc2$call<-ifelse(((u_conns_pc2$end1>1000&u_conns_pc2$end1<2000)&u_conns_pc2$end2>2000)|
                           ((u_conns_pc2$end1>9&u_conns_pc2$end1<27)&u_conns_pc2$end2>2000)|
                           ((u_conns_pc2$end1>9&u_conns_pc2$end1<27)&(u_conns_pc2$end2>48&u_conns_pc2$end2<59))|
                           ((u_conns_pc2$end1>48&u_conns_pc2$end1<59)&(u_conns_pc2$end2>1000&u_conns_pc2$end2<2000)),1,0)


colnames(xLF1.pc1.bootsig)<-c("sig","conn","signum")

#LF1
xLF1.pc1.bootsig$end1<-substr(xLF1.pc1.bootsig$conn,1,regexpr("_",xLF1.pc1.bootsig$conn)-1)
xLF1.pc1.bootsig$end2<-substr(xLF1.pc1.bootsig$conn,regexpr("_",xLF1.pc1.bootsig$conn)+1,20)
xLF1.pc1.bootsig$call<-ifelse(((xLF1.pc1.bootsig$end1>1000&xLF1.pc1.bootsig$end1<2000)&xLF1.pc1.bootsig$end2>2000)|
                                ((xLF1.pc1.bootsig$end1>9&xLF1.pc1.bootsig$end1<27)&xLF1.pc1.bootsig$end2>2000)|
                                ((xLF1.pc1.bootsig$end1>9&xLF1.pc1.bootsig$end1<27)&(xLF1.pc1.bootsig$end2>48&xLF1.pc1.bootsig$end2<59))|
                                ((xLF1.pc1.bootsig$end1>48&xLF1.pc1.bootsig$end1<59)&(xLF1.pc1.bootsig$end2>1000&xLF1.pc1.bootsig$end2<2000)),1,0)
#If cortical left -> cortical right OR subcort left to cort right OR subcort left to subcort right OR subcort right to cort left


colnames(xLF7.pc1.bootsig)<-c("sig","conn","signum")
#LF7
xLF7.pc1.bootsig$end1<-substr(xLF7.pc1.bootsig$conn,1,regexpr("_",xLF7.pc1.bootsig$conn)-1)
xLF7.pc1.bootsig$end2<-substr(xLF7.pc1.bootsig$conn,regexpr("_",xLF7.pc1.bootsig$conn)+1,20)
xLF7.pc1.bootsig$call<-ifelse(((xLF7.pc1.bootsig$end1>1000&xLF7.pc1.bootsig$end1<2000)&xLF7.pc1.bootsig$end2>2000)|
                                ((xLF7.pc1.bootsig$end1>9&xLF7.pc1.bootsig$end1<27)&xLF7.pc1.bootsig$end2>2000)|
                                ((xLF7.pc1.bootsig$end1>9&xLF7.pc1.bootsig$end1<27)&(xLF7.pc1.bootsig$end2>48&xLF7.pc1.bootsig$end2<59))|
                                ((xLF7.pc1.bootsig$end1>48&xLF7.pc1.bootsig$end1<59)&(xLF7.pc1.bootsig$end2>1000&xLF7.pc1.bootsig$end2<2000)),1,0)


#Calculate proportion of callosal tracts over sig tracts
sum(u_conns_pc1[u_conns_pc1$sig==1,c("call")])/sum(u_conns_pc1$sig)*100#7/39
sum(u_conns_pc2[u_conns_pc2$sig==1,c("call")])/sum(u_conns_pc2$sig)*100#18/50
sum(xLF1.pc1.bootsig[xLF1.pc1.bootsig$sig==1,c("call")])/sum(xLF1.pc1.bootsig$sig)*100#6/24
sum(xLF7.pc1.bootsig[xLF7.pc1.bootsig$sig==1,c("call")])/sum(xLF7.pc1.bootsig$sig)*100#2/11



tmp1<-perc_sims[1:19,c("perc_overlap")]
tmp2<-colSums(boot.pc1$bootRatiosSignificant.i)
cbind(tmp1,tmp2)



####Sensitivity analsyes: Impact of Reshuffling strategy####
#plsc_and_permute_ver2 is the modified function. Note that I'm computing it on the basis of covariance explained, shuffling the labels for the symptoms
plsc.t09.pc1.ver2<-plsc_and_permute_ver2(PC1.scores.wide[,-which(colnames(PC1.scores.wide)%in%c("subjectkey"))],mydata_t09$yy_train[,-which(colnames(mydata_t09$yy_train)%in%c("subjectkey","sex"))],
                                         nperms = 2000,topn = 200,out_flag = 1,num_comps = 20)



plsc.t09.pc2.ver2<-plsc_and_permute_ver2(PC2.scores.wide[,-which(colnames(PC2.scores.wide)%in%c("subjectkey"))],mydata_t09$yy_train[,-which(colnames(mydata_t09$yy_train)%in%c("subjectkey","sex"))],
                                         nperms = 2000,topn = 200,out_flag = 1, num_comps=20)




####Sensitivity analysis: Why 200 features?####
Xinput<-uvar_feat_sel(mydata_t09$yy_train[,-which(colnames(mydata_t09$yy_train)%in%c("subjectkey","sex"))],PC1.scores.wide[,-which(colnames(PC1.scores.wide)%in%c("subjectkey"))],200,"top")
Yinput<-mydata_t09$yy_train[,-which(colnames(mydata_t09$yy_train)%in%c("subjectkey","sex"))]
dimX<-dim(Xinput)[2]
dimY<-dim(Yinput)[2]
dimsubs<-length(mydata_t09$yy_train$subjectkey)

#PLS1
cov_exps <- foreach(i=min(dimX,dimY):200) %dopar% {
  #Note: I'm starting with 19 which is the same number of symptoms. This avoids having to deal with a changing number of LFs.
  covs <- apply_pls(Xinput[,1:i],Yinput,ncomps = 20,scaleflag = TRUE,permflag = TRUE,outflag = 2)
  list(covs)
}

# transform null results lists to data frame
all_covs <- lapply(cov_exps, function(x){return(x[[1]])})
all_covs <- do.call(rbind, all_covs)

cove_data<-data.frame(cove=all_covs[,1],feats=min(dimX,dimY):max(dimX,dimY))#LF1

#CREATE FIGURE
ggplot(cove_data, aes(x=feats)) +
  geom_line( aes(y=cove), size=2, color="black") + 
  scale_y_continuous(
    # Features of the first axis
    name = "% Covariance Explained"
  ) +
  scale_x_continuous(
    # Features of the first axis
    name = "Number of Features Selected"
  ) +
  theme_classic(base_size = 20)+
  theme(axis.title=element_text(face="bold",colour = "black"),
        axis.text=element_text(face="bold",colour = "black"))
##

####Sensitivity analysis: Impact of Regressing out site####
mydata_t09_noScanRO<-pre_proc(conn_data = str_conn_feats_QCpass,symp_data = symptoms_mTBI,nuisance_dat = nuisance_covs,metrics=c("fa","ad","md","rd","nufo","afdf"),
                              nuisance_vars = c("sex","ehi_y_ss_scoreb","pubertal_stage"),nuisance_vars_symps =c("sex","ehi_y_ss_scoreb","pubertal_stage"),
                              thresh = 0.9,tr_ratio = .70,out_dir = "/Users/Guido/Desktop/Projects/PCHStudy/Data/",resume_thresh = TRUE,thresh_file = '/Users/Guido/Desktop/Projects/PCHStudy/Data/conn_features_thresh_0.9',
                              resume_regout = FALSE,resume_scale = FALSE)

mydata_t09_noScanRO$xx_train$fa
mydata_t09_noScanRO$nuisance_train


tmp_dims<-c(dim(mydata_t09_noScanRO$xx_train$fa)[2],dim(mydata_t09_noScanRO$xx_train$md)[2],dim(mydata_t09_noScanRO$xx_train$rd)[2],dim(mydata_t09_noScanRO$xx_train$ad)[2],
            dim(mydata_t09_noScanRO$xx_train$nufo)[2],dim(mydata_t09_noScanRO$xx_train$afdf)[2])
tmp_mtrs<-c("fa","md","rd","ad","nufo","afdf")
tmp_dim_mtrs<-rbind(tmp_dims,tmp_mtrs)
min_mtr<-tmp_dim_mtrs[rownames(tmp_dim_mtrs)=="tmp_mtrs",which(tmp_dim_mtrs[1,]==min(tmp_dim_mtrs[1,]))][1]
max_mtr<-tmp_dim_mtrs[rownames(tmp_dim_mtrs)=="tmp_mtrs",which(tmp_dim_mtrs[1,]==max(tmp_dim_mtrs[1,]))][1]
# tmpbigcols<-extract_colnames(input_colnames = colnames(as.data.frame(mydata_t09$xx_train[max_mtr])))

tmpbigcols<-substr(colnames(as.data.frame(mydata_t09_noScanRO$xx_train[max_mtr])),regexpr("\\.",colnames(as.data.frame(mydata_t09_noScanRO$xx_train[max_mtr])))+1,20)#for some reason, this line doesn't extract the conn number alone anymore
tmpbigcols<-extract_colnames(tmpbigcols)

# tmpbigcols<-substr(tmpbigcols,regexpr("\\.",tmpbigcols)+1,20)
# tmpbigcols<-substr(tmpbigcols,1,regexpr("\\.",tmpbigcols)-1)


tmp_misscols<-cbind(as.vector(tmpbigcols[!tmpbigcols%in%extract_colnames(colnames(mydata_t09_noScanRO$xx_train$fa))]),as.vector(tmpbigcols[!tmpbigcols%in%extract_colnames(colnames(mydata_t09_noScanRO$xx_train$md))]),
                    as.vector(tmpbigcols[!tmpbigcols%in%extract_colnames(colnames(mydata_t09_noScanRO$xx_train$rd))]),as.vector(tmpbigcols[!tmpbigcols%in%extract_colnames(colnames(mydata_t09_noScanRO$xx_train$ad))]),
                    as.vector(tmpbigcols[!tmpbigcols%in%extract_colnames(colnames(mydata_t09_noScanRO$xx_train$nufo))]),as.vector(tmpbigcols[!tmpbigcols%in%extract_colnames(colnames(mydata_t09_noScanRO$xx_train$afdf))]))
misscols<-unique(melt(tmp_misscols)$value)
tmpsmallcols<-tmpbigcols[!tmpbigcols%in%misscols]

# tmpsmallcols<-substr(colnames(as.data.frame(mydata_t09$xx_train[min_mtr])),regexpr("\\.",colnames(as.data.frame(mydata_t09$xx_train[min_mtr])))+1,20)


tmp_submat<-matrix(data=NA,ncol = length(tmpsmallcols),nrow = length(mydata_t09_noScanRO$yy_train$subjectkey))
colnames(tmp_submat)<-tmpsmallcols
tmp_submat<-as.data.frame(tmp_submat)
tmp_submat$subjectkey<-mydata_t09_noScanRO$yy_train$subjectkey
tmp_submat_melt<-melt(tmp_submat,id.vars = "subjectkey")

if (length(tmpbigcols)>length(tmpsmallcols)){
  #Note that here I'm removing those "extra" conns that metrics other than NuFO had.
  
  length(mydata_t09_noScanRO$xx_train$fa[,which(extract_colnames(colnames(mydata_t09_noScanRO$xx_train$fa))%in%tmpsmallcols)])
  length(mydata_t09_noScanRO$xx_train$md[,which(extract_colnames(colnames(mydata_t09_noScanRO$xx_train$md))%in%tmpsmallcols)])
  length(mydata_t09_noScanRO$xx_train$rd[,which(extract_colnames(colnames(mydata_t09_noScanRO$xx_train$rd))%in%tmpsmallcols)])
  length(mydata_t09_noScanRO$xx_train$ad[,which(extract_colnames(colnames(mydata_t09_noScanRO$xx_train$ad))%in%tmpsmallcols)])
  length(mydata_t09_noScanRO$xx_train$nufo[,which(extract_colnames(colnames(mydata_t09_noScanRO$xx_train$nufo))%in%tmpsmallcols)])
  length(mydata_t09_noScanRO$xx_train$afdf[,which(extract_colnames(colnames(mydata_t09_noScanRO$xx_train$afd))%in%tmpsmallcols)])
  
  mydata_t09_noScanRO$xx_train$all_metrics<-melt(mydata_t09_noScanRO$xx_train$fa[,which(extract_colnames(colnames(mydata_t09_noScanRO$xx_train$fa))%in%tmpsmallcols)])
  mydata_t09_noScanRO$xx_train$all_metrics$md<-melt(mydata_t09_noScanRO$xx_train$md[,which(extract_colnames(colnames(mydata_t09_noScanRO$xx_train$md))%in%tmpsmallcols)])$value
  mydata_t09_noScanRO$xx_train$all_metrics$rd<-melt(mydata_t09_noScanRO$xx_train$rd[,which(extract_colnames(colnames(mydata_t09_noScanRO$xx_train$rd))%in%tmpsmallcols)])$value
  mydata_t09_noScanRO$xx_train$all_metrics$ad<-melt(mydata_t09_noScanRO$xx_train$ad[,which(extract_colnames(colnames(mydata_t09_noScanRO$xx_train$ad))%in%tmpsmallcols)])$value
  mydata_t09_noScanRO$xx_train$all_metrics$nufo<-melt(mydata_t09_noScanRO$xx_train$nufo[,which(extract_colnames(colnames(mydata_t09_noScanRO$xx_train$nufo))%in%tmpsmallcols)])$value
  mydata_t09_noScanRO$xx_train$all_metrics$afdf<-melt(mydata_t09_noScanRO$xx_train$afdf[,which(extract_colnames(colnames(mydata_t09_noScanRO$xx_train$afdf))%in%tmpsmallcols)])$value
  mydata_t09_noScanRO$xx_train$all_metrics$subjectkey<-tmp_submat_melt$subjectkey
  mydata_t09_noScanRO$xx_train$all_metrics$variable<-tmp_submat_melt$variable
  colnames(mydata_t09_noScanRO$xx_train$all_metrics)<-c("conn","fa","md","rd","ad","nufo","afdf","subjectkey")
  
} else if (length(tmpbigcols)==length(tmpsmallcols)){
  #Note that here I'm removing those "extra" conns that metrics other than NuFO had.
  mydata_t09$xx_train$all_metrics<-melt(mydata_t09$xx_train$fa)
  mydata_t09$xx_train$all_metrics$md<-melt(mydata_t09$xx_train$md)$value
  mydata_t09$xx_train$all_metrics$rd<-melt(mydata_t09$xx_train$rd)$value
  mydata_t09$xx_train$all_metrics$ad<-melt(mydata_t09$xx_train$ad)$value
  mydata_t09$xx_train$all_metrics$nufo<-melt(mydata_t09$xx_train$nufo)$value
  mydata_t09$xx_train$all_metrics$afdf<-melt(mydata_t09$xx_train$afdf)$value
  mydata_t09$xx_train$all_metrics$subjectkey<-tmp_submat_melt$subjectkey
  colnames(mydata_t09$xx_train$all_metrics)<-c("conn","fa","md","rd","ad","nufo","afdf","subjectkey")
}

###PCA ACROSS MODALITIES###

pca.across.metrics<-prcomp(mydata_t09$xx_train$all_metrics[mydata_t09$xx_train$all_metrics$conn%in%selected_conns[selected_conns$c1w_pos==1,c("selected_conns")],
                                                           -which(colnames(mydata_t09$xx_train$all_metrics)%in%c("subjectkey","conn"))],center = TRUE,scale. = TRUE)

#Notice that I'm subsetting by the other mydata_t09, because they should have the same structure.
pca.across.metrics.noScanRO<-prcomp(mydata_t09_noScanRO$xx_train$all_metrics[mydata_t09$xx_train$all_metrics$conn%in%selected_conns[selected_conns$c1w_pos==1,c("selected_conns")],
                                                                    -which(colnames(mydata_t09_noScanRO$xx_train$all_metrics)%in%c("subjectkey","conn"))],center = TRUE,scale. = TRUE)

corrplot(pca.across.metrics.noScanRO$rotation, is.corr=FALSE)

###EXTRACT FIRST 2 PCs###
#Extract first 2 PCS, they account for 97% of variance (EVEN WITHOUT SCAN RO)
PC1.scores.noScanRO<-as.data.frame(pca.across.metrics.noScanRO$x[,1])
colnames(PC1.scores.noScanRO)<-"score"
PC2.scores.noScanRO<-as.data.frame(pca.across.metrics.noScanRO$x[,2])
colnames(PC2.scores.noScanRO)<-"score"
PC1.scores.noScanRO$conns<-mydata_t09_noScanRO$xx_train$all_metrics[mydata_t09$xx_train$all_metrics$conn%in%selected_conns[selected_conns$c1w_pos==1,c("selected_conns")],c("conn")]
PC1.scores.noScanRO$subjectkey<-mydata_t09_noScanRO$xx_train$all_metrics[mydata_t09$xx_train$all_metrics$conn%in%selected_conns[selected_conns$c1w_pos==1,c("selected_conns")],c("subjectkey")]
PC2.scores.noScanRO$conns<-mydata_t09_noScanRO$xx_train$all_metrics[mydata_t09$xx_train$all_metrics$conn%in%selected_conns[selected_conns$c1w_pos==1,c("selected_conns")],c("conn")]
PC2.scores.noScanRO$subjectkey<-mydata_t09_noScanRO$xx_train$all_metrics[mydata_t09$xx_train$all_metrics$conn%in%selected_conns[selected_conns$c1w_pos==1,c("selected_conns")],c("subjectkey")]



###RESTRUCTURE PCS BACK TO WIDE FORMAT###
PC1.scores.wide.noScanRO<-reshape(PC1.scores.noScanRO,idvar = "subjectkey",timevar = "conns",direction = "wide")
colnames(PC1.scores.wide.noScanRO)<-substr(colnames(PC1.scores.wide.noScanRO),regexpr("\\.",colnames(PC1.scores.wide.noScanRO))+1,25)

PC2.scores.wide.noScanRO<-reshape(PC2.scores.noScanRO,idvar = "subjectkey",timevar = "conns",direction = "wide")
colnames(PC2.scores.wide.noScanRO)<-substr(colnames(PC2.scores.wide.noScanRO),regexpr("\\.",colnames(PC2.scores.wide.noScanRO))+1,25)


###PLSc ON PCs### 
plsc.t09.pc1.noScanRO<-plsc_and_permute(PC1.scores.wide.noScanRO[,-which(colnames(PC1.scores.wide.noScanRO)%in%c("subjectkey"))],
                                        mydata_t09_noScanRO$yy_train[,-which(colnames(mydata_t09_noScanRO$yy_train)%in%c("subjectkey","sex"))],
                                        num_comps = 20,nperms = 2000,topn = 200,out_flag = 1)

plsc.t09.pc2.noScanRO<-plsc_and_permute(PC2.scores.wide.noScanRO[,-which(colnames(PC2.scores.wide.noScanRO)%in%c("subjectkey"))],
                                        mydata_t09_noScanRO$yy_train[,-which(colnames(mydata_t09_noScanRO$yy_train)%in%c("subjectkey","sex"))],
                                        num_comps = 20,nperms = 2000,topn = 200,out_flag = 1)


mydata_t09$nuisance_train[,c("subjectkey","mri_info_deviceserialnumber")]
mydata_t09_noScanRO$nuisance_train[,c("subjectkey","mri_info_deviceserialnumber")]

cor(plsc.t09.pc1.noScanRO[[4]]$TExPosition.Data$lx,plsc.t09.pc1[[4]]$TExPosition.Data$lx)

SiteRO_x_projections<-plsc.t09.pc1[[4]]$TExPosition.Data$lx
colnames(SiteRO_x_projections)<-paste("X_",1:19,sep = "")
SiteRO_y_projections<-plsc.t09.pc1[[4]]$TExPosition.Data$ly
colnames(SiteRO_y_projections)<-paste("Y_",1:19,sep = "")

NoSiteRO_x_projections<-plsc.t09.pc1.noScanRO[[4]]$TExPosition.Data$lx
colnames(NoSiteRO_x_projections)<-paste("X_",1:19,sep = "")
NoSiteRO_y_projections<-plsc.t09.pc1.noScanRO[[4]]$TExPosition.Data$ly
colnames(NoSiteRO_y_projections)<-paste("Y_",1:19,sep = "")

SiteRO_projections<-cbind(SiteRO_x_projections,SiteRO_y_projections,mydata_t09$nuisance_train[,c("subjectkey","mri_info_deviceserialnumber")])
NoSiteRO_projections<-cbind(NoSiteRO_x_projections,NoSiteRO_y_projections,mydata_t09_noScanRO$nuisance_train[,c("subjectkey","mri_info_deviceserialnumber")])

summSiteRO_projections<-SiteRO_projections %>% group_by(mri_info_deviceserialnumber) %>% summarize(mean=mean(X_1))
summNoSiteRO_projections<-NoSiteRO_projections %>% group_by(mri_info_deviceserialnumber) %>% summarize(mean=mean(X_1))

summNoSiteRO_projections[summNoSiteRO_projections$mean==max(summNoSiteRO_projections$mean),c("mri_info_deviceserialnumber")]
summNoSiteRO_projections[summNoSiteRO_projections$mean==min(summNoSiteRO_projections$mean),c("mri_info_deviceserialnumber")]

# HASH5b2fcf80	-1.8783003
# HASH4b0b8b05	-1.7329407
# HASHc3bf3d9c	-1.5038814
# HASHfeb7e81a	-1.4208568
# HASH03db707f	-1.2354671
# HASHa3e45734	-0.8362058
# HASHd7cb4c6d	-0.8111347
# HASH69f406fa	-0.7420654
# HASH5b0cf1bb	-0.1840217
# HASHe3ce02d3	-0.1234141
# HASH7911780b	-0.1185814
# HASH3935c89e	0.1847400
# HASHe4f6957a	0.1970060
# HASH1314a204	0.2398370
# HASH4d1ed7b1	0.3652730
# HASH311170b9	0.4596933
# HASHdb2589d4	0.5776911
# HASHd422be27	0.6920907
# HASHc9398971	0.9123946
# HASHb640a1b8	0.9485583
# HASH96a0c182	0.9755831
# HASH65b39280	1.0009283
# HASH7f91147d	1.0442244
# HASH6b4422a7	1.7313422
# HASH11ad4ed5	1.8638408

#Note: I selected two where the effect was more obvious

SiteRO_projections$col<-ifelse(SiteRO_projections$mri_info_deviceserialnumber=="HASH11ad4ed5","blue",ifelse(SiteRO_projections$mri_info_deviceserialnumber=="HASH5b2fcf80","red","gray"))
NoSiteRO_projections$col<-ifelse(NoSiteRO_projections$mri_info_deviceserialnumber=="HASH11ad4ed5","blue",ifelse(NoSiteRO_projections$mri_info_deviceserialnumber=="HASH5b2fcf80","red","gray"))

p1<-ggplot(SiteRO_projections,aes(x=X_1*-1, y=Y_1*-1,color="gray"))+
  geom_point(colour=SiteRO_projections$col,size=4.5,alpha=0.7)+
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill="white",color="white"),
        axis.line =element_line(color = "black", size = 1, linetype = "solid"),
        axis.text = element_text(face = "bold", color = "black", 
                                 size = 15),
        axis.title = element_text(face = "bold", color = "black", 
                                  size = 20),
        axis.ticks = element_blank())+
  labs(y="Multi-Symptom Feature 1", 
       x="Multi-Tract Connectivity Feature 1")+
  theme(legend.position="none")
p1    


p1<-ggplot(NoSiteRO_projections,aes(x=X_1*-1, y=Y_1*-1,color="gray"))+
  geom_point(colour=NoSiteRO_projections$col,size=4.5,alpha=0.7)+
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill="white",color="white"),
        axis.line =element_line(color = "black", size = 1, linetype = "solid"),
        axis.text = element_text(face = "bold", color = "black", 
                                 size = 15),
        axis.title = element_text(face = "bold", color = "black", 
                                  size = 20),
        axis.ticks = element_blank())+
  labs(y="Multi-Symptom Feature 1", 
       x="Multi-Tract Connectivity Feature 1")+
  theme(legend.position="none")
p1    


p<-ggplot(data=summNoSiteRO_projections,aes(x=mri_info_deviceserialnumber, y=mean*-1)) +
  geom_bar(stat="identity",fill=factor(ifelse(summNoSiteRO_projections$mri_info_deviceserialnumber=="HASH11ad4ed5","blue",ifelse(summNoSiteRO_projections$mri_info_deviceserialnumber=="HASH5b2fcf80","red","gray"))))+ylim(-3,3)+
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill="white",color="white"),
        axis.line =element_line(color = "black", size = 1, linetype = "solid"),
        axis.text = element_text(face = "bold", color = "black", 
                                 size = 15),
        axis.text.x = element_text(angle = 315,vjust = 0.5, hjust=0),
        axis.title = element_text(face = "bold", color = "black", 
                                  size = 20),
        axis.ticks = element_blank(),
        plot.margin = margin(1,2.5,1,1,"cm")
  )+
  labs(y="Average Expression", 
       x="Scanner")+
  theme(legend.position="none")
p


####Sensitivity Analysis: Different Thresholds####

mydata_t085<-pre_proc(conn_data = str_conn_feats_QCpass,symp_data = symptoms_mTBI,nuisance_dat = nuisance_covs,metrics=c("fa","ad","md","rd","nufo","afdf"),
                      nuisance_vars = c("sex","ehi_y_ss_scoreb","mri_info_deviceserialnumber","pubertal_stage"),nuisance_vars_symps =c("sex","ehi_y_ss_scoreb","pubertal_stage"),
                      thresh = 0.85,tr_ratio = .70,out_dir = "/Users/Guido/Desktop/Projects/PCHStudy/Data/",resume_thresh = TRUE,thresh_file = '/Users/Guido/Desktop/Projects/PCHStudy/Data/conn_features_thresh_0.85',
                      resume_regout = FALSE,resume_scale = FALSE)

mydata_t1<-pre_proc(conn_data = str_conn_feats_QCpass,symp_data = symptoms_mTBI,nuisance_dat = nuisance_covs,metrics=c("fa","ad","md","rd","nufo","afdf"),
                    nuisance_vars = c("sex","ehi_y_ss_scoreb","mri_info_deviceserialnumber","pubertal_stage"),nuisance_vars_symps =c("sex","ehi_y_ss_scoreb","pubertal_stage"),
                    thresh = 1,tr_ratio = .70,out_dir = "/Users/Guido/Desktop/Projects/PCHStudy/Data/",resume_thresh = TRUE,thresh_file = '/Users/Guido/Desktop/Projects/PCHStudy/Data/conn_features_thresh_1',
                    resume_regout = FALSE,resume_scale = FALSE)

mydata_t095<-pre_proc(conn_data = str_conn_feats_QCpass,symp_data = symptoms_mTBI,nuisance_dat = nuisance_covs,metrics=c("fa","ad","md","rd","nufo","afdf"),
                      nuisance_vars = c("sex","ehi_y_ss_scoreb","mri_info_deviceserialnumber","pubertal_stage"),nuisance_vars_symps =c("sex","ehi_y_ss_scoreb","pubertal_stage"),
                      thresh = 0.95,tr_ratio = .70,out_dir = "/Users/Guido/Desktop/Projects/PCHStudy/Data/",resume_thresh = TRUE,thresh_file = '/Users/Guido/Desktop/Projects/PCHStudy/Data/conn_features_thresh_0.95',
                      resume_regout = FALSE,resume_scale = FALSE)

mydata_t09_backup<-mydata_t09###@@@@@@@@@@@@@@@@@RUN TO BACK UP###
mydata_t09<-mydata_t09_backup##@@@@@@@@@@@@@@@@@@RUN TO RECOVER T09###

mydata_t09<-mydata_t1##@@@@@@@@@@@@@@@@@@@@RUN TO AVOID HAVING TO REWRITE EVERYTHING WITH A DIFFERENT THRESHOLD###


####@85%@####
####@@SET UP ALLDATA@@####
tmp_dims<-c(dim(mydata_t085$xx_train$fa)[2],dim(mydata_t085$xx_train$md)[2],dim(mydata_t085$xx_train$rd)[2],dim(mydata_t085$xx_train$ad)[2],dim(mydata_t085$xx_train$nufo)[2],dim(mydata_t085$xx_train$afdf)[2])
tmp_mtrs<-c("fa","md","rd","ad","nufo","afdf")
tmp_dim_mtrs<-rbind(tmp_dims,tmp_mtrs)
min_mtr<-tmp_dim_mtrs[rownames(tmp_dim_mtrs)=="tmp_mtrs",which(tmp_dim_mtrs[1,]==min(tmp_dim_mtrs[1,]))][1]
max_mtr<-tmp_dim_mtrs[rownames(tmp_dim_mtrs)=="tmp_mtrs",which(tmp_dim_mtrs[1,]==max(tmp_dim_mtrs[1,]))][1]
# tmpbigcols<-extract_colnames(input_colnames = colnames(as.data.frame(mydata_t09$xx_train[max_mtr])))

tmpbigcols<-substr(colnames(as.data.frame(mydata_t085$xx_train[max_mtr])),regexpr("\\.",colnames(as.data.frame(mydata_t085$xx_train[max_mtr])))+1,20)#for some reason, this line doesn't extract the conn number alone anymore
tmpbigcols<-extract_colnames(tmpbigcols)

# tmpbigcols<-substr(tmpbigcols,regexpr("\\.",tmpbigcols)+1,20)
# tmpbigcols<-substr(tmpbigcols,1,regexpr("\\.",tmpbigcols)-1)


tmp_misscols<-cbind(as.vector(tmpbigcols[!tmpbigcols%in%extract_colnames(colnames(mydata_t085$xx_train$fa))]),as.vector(tmpbigcols[!tmpbigcols%in%extract_colnames(colnames(mydata_t085$xx_train$md))]),
                    as.vector(tmpbigcols[!tmpbigcols%in%extract_colnames(colnames(mydata_t085$xx_train$rd))]),as.vector(tmpbigcols[!tmpbigcols%in%extract_colnames(colnames(mydata_t085$xx_train$ad))]),
                    as.vector(tmpbigcols[!tmpbigcols%in%extract_colnames(colnames(mydata_t085$xx_train$nufo))]),as.vector(tmpbigcols[!tmpbigcols%in%extract_colnames(colnames(mydata_t085$xx_train$afdf))]))
misscols<-unique(melt(tmp_misscols)$value)
tmpsmallcols<-tmpbigcols[!tmpbigcols%in%misscols]

# tmpsmallcols<-substr(colnames(as.data.frame(mydata_t09$xx_train[min_mtr])),regexpr("\\.",colnames(as.data.frame(mydata_t09$xx_train[min_mtr])))+1,20)


tmp_submat<-matrix(data=NA,ncol = length(tmpsmallcols),nrow = length(mydata_t085$yy_train$subjectkey))
colnames(tmp_submat)<-tmpsmallcols
tmp_submat<-as.data.frame(tmp_submat)
tmp_submat$subjectkey<-mydata_t085$yy_train$subjectkey
tmp_submat_melt<-melt(tmp_submat,id.vars = "subjectkey")

if (length(tmpbigcols)>length(tmpsmallcols)){
  #Note that here I'm removing those "extra" conns that metrics other than NuFO had.
  # mydata_t09$xx_train$all_metrics<-melt(mydata_t09$xx_train$fa[,which(colnames(mydata_t09$xx_train$fa)%in%tmpsmallcols)|
  #                                                                which(colnames(mydata_t09$xx_train$fa)%in%"subjectkey")],id.vars = "subjectkey")
  
  length(mydata_t085$xx_train$fa[,which(extract_colnames(colnames(mydata_t085$xx_train$fa))%in%tmpsmallcols)])
  length(mydata_t085$xx_train$md[,which(extract_colnames(colnames(mydata_t085$xx_train$md))%in%tmpsmallcols)])
  length(mydata_t085$xx_train$rd[,which(extract_colnames(colnames(mydata_t085$xx_train$rd))%in%tmpsmallcols)])
  length(mydata_t085$xx_train$ad[,which(extract_colnames(colnames(mydata_t085$xx_train$ad))%in%tmpsmallcols)])
  length(mydata_t085$xx_train$nufo[,which(extract_colnames(colnames(mydata_t085$xx_train$nufo))%in%tmpsmallcols)])
  length(mydata_t085$xx_train$afdf[,which(extract_colnames(colnames(mydata_t085$xx_train$afd))%in%tmpsmallcols)])
  
  mydata_t085$xx_train$all_metrics<-melt(mydata_t085$xx_train$fa[,which(extract_colnames(colnames(mydata_t085$xx_train$fa))%in%tmpsmallcols)])
  mydata_t085$xx_train$all_metrics$md<-melt(mydata_t085$xx_train$md[,which(extract_colnames(colnames(mydata_t085$xx_train$md))%in%tmpsmallcols)])$value
  mydata_t085$xx_train$all_metrics$rd<-melt(mydata_t085$xx_train$rd[,which(extract_colnames(colnames(mydata_t085$xx_train$rd))%in%tmpsmallcols)])$value
  mydata_t085$xx_train$all_metrics$ad<-melt(mydata_t085$xx_train$ad[,which(extract_colnames(colnames(mydata_t085$xx_train$ad))%in%tmpsmallcols)])$value
  mydata_t085$xx_train$all_metrics$nufo<-melt(mydata_t085$xx_train$nufo[,which(extract_colnames(colnames(mydata_t085$xx_train$nufo))%in%tmpsmallcols)])$value
  mydata_t085$xx_train$all_metrics$afdf<-melt(mydata_t085$xx_train$afdf[,which(extract_colnames(colnames(mydata_t085$xx_train$afdf))%in%tmpsmallcols)])$value
  mydata_t085$xx_train$all_metrics$subjectkey<-tmp_submat_melt$subjectkey
  mydata_t085$xx_train$all_metrics$variable<-tmp_submat_melt$variable
  colnames(mydata_t085$xx_train$all_metrics)<-c("conn","fa","md","rd","ad","nufo","afdf","subjectkey")
  
} else if (length(tmpbigcols)==length(tmpsmallcols)){
  #Note that here I'm removing those "extra" conns that metrics other than NuFO had.
  mydata_t085$xx_train$all_metrics<-melt(mydata_t085$xx_train$fa)
  mydata_t085$xx_train$all_metrics$md<-melt(mydata_t085$xx_train$md)$value
  mydata_t085$xx_train$all_metrics$rd<-melt(mydata_t085$xx_train$rd)$value
  mydata_t085$xx_train$all_metrics$ad<-melt(mydata_t085$xx_train$ad)$value
  mydata_t085$xx_train$all_metrics$nufo<-melt(mydata_t085$xx_train$nufo)$value
  mydata_t085$xx_train$all_metrics$afdf<-melt(mydata_t085$xx_train$afdf)$value
  mydata_t085$xx_train$all_metrics$subjectkey<-tmp_submat_melt$subjectkey
  colnames(mydata_t085$xx_train$all_metrics)<-c("conn","fa","md","rd","ad","nufo","afdf","subjectkey")
}

#Now modify the connection names because they have the metric in them. In all_metrics only .fa appears, so I need to replace that one only
if(regexpr("fa",mydata_t085$xx_train$all_metrics$conn)>1){
  mydata_t085$xx_train$all_metrics$conn<-substr(mydata_t085$xx_train$all_metrics$conn,regexpr("\\.",mydata_t085$xx_train$all_metrics$conn)+1,40)
  mydata_t085$xx_train$all_metrics$conn<-substr(mydata_t085$xx_train$all_metrics$conn,1,regexpr("\\.",mydata_t085$xx_train$all_metrics$conn)-1)
}

###PCA###

pca.across.metrics.t085<-prcomp(mydata_t085$xx_train$all_metrics[mydata_t085$xx_train$all_metrics$conn%in%selected_conns[selected_conns$c1w_pos==1,c("selected_conns")],
                                                           -which(colnames(mydata_t085$xx_train$all_metrics)%in%c("subjectkey","conn"))],center = TRUE,scale. = TRUE)

corrplot(pca.across.metrics.t085$rotation, is.corr=FALSE)
summary(pca.across.metrics.t085)
#0.8089 0.1626


###EXTRACT FIRST 2 PCs###
#Extract first 2 PCS, they account for 97% of variance
PC1.scores.t085<-as.data.frame(pca.across.metrics.t085$x[,1])
colnames(PC1.scores.t085)<-"score"
PC2.scores.t085<-as.data.frame(pca.across.metrics.t085$x[,2])
colnames(PC2.scores.t085)<-"score"
PC1.scores.t085$conns<-mydata_t085$xx_train$all_metrics[mydata_t085$xx_train$all_metrics$conn%in%selected_conns[selected_conns$c1w_pos==1,c("selected_conns")],c("conn")]
PC1.scores.t085$subjectkey<-mydata_t085$xx_train$all_metrics[mydata_t085$xx_train$all_metrics$conn%in%selected_conns[selected_conns$c1w_pos==1,c("selected_conns")],c("subjectkey")]
PC2.scores.t085$conns<-mydata_t085$xx_train$all_metrics[mydata_t085$xx_train$all_metrics$conn%in%selected_conns[selected_conns$c1w_pos==1,c("selected_conns")],c("conn")]
PC2.scores.t085$subjectkey<-mydata_t085$xx_train$all_metrics[mydata_t085$xx_train$all_metrics$conn%in%selected_conns[selected_conns$c1w_pos==1,c("selected_conns")],c("subjectkey")]


###RESTRUCTURE PCS BACK TO WIDE FORMAT###
PC1.scores.wide.t085<-reshape(PC1.scores.t085,idvar = "subjectkey",timevar = "conns",direction = "wide")
colnames(PC1.scores.wide.t085)<-substr(colnames(PC1.scores.wide.t085),regexpr("\\.",colnames(PC1.scores.wide.t085))+1,25)

PC2.scores.wide.t085<-reshape(PC2.scores.t085,idvar = "subjectkey",timevar = "conns",direction = "wide")
colnames(PC2.scores.wide.t085)<-substr(colnames(PC2.scores.wide.t085),regexpr("\\.",colnames(PC2.scores.wide.t085))+1,25)




####95%####
###SET UP ALLDATA###
tmp_dims<-c(dim(mydata_t095$xx_train$fa)[2],dim(mydata_t095$xx_train$md)[2],dim(mydata_t095$xx_train$rd)[2],dim(mydata_t095$xx_train$ad)[2],dim(mydata_t095$xx_train$nufo)[2],dim(mydata_t095$xx_train$afdf)[2])
tmp_mtrs<-c("fa","md","rd","ad","nufo","afdf")
tmp_dim_mtrs<-rbind(tmp_dims,tmp_mtrs)
min_mtr<-tmp_dim_mtrs[rownames(tmp_dim_mtrs)=="tmp_mtrs",which(tmp_dim_mtrs[1,]==min(tmp_dim_mtrs[1,]))][1]
max_mtr<-tmp_dim_mtrs[rownames(tmp_dim_mtrs)=="tmp_mtrs",which(tmp_dim_mtrs[1,]==max(tmp_dim_mtrs[1,]))][1]
# tmpbigcols<-extract_colnames(input_colnames = colnames(as.data.frame(mydata_t09$xx_train[max_mtr])))

tmpbigcols<-substr(colnames(as.data.frame(mydata_t095$xx_train[max_mtr])),regexpr("\\.",colnames(as.data.frame(mydata_t095$xx_train[max_mtr])))+1,20)#for some reason, this line doesn't extract the conn number alone anymore
tmpbigcols<-extract_colnames(tmpbigcols)

# tmpbigcols<-substr(tmpbigcols,regexpr("\\.",tmpbigcols)+1,20)
# tmpbigcols<-substr(tmpbigcols,1,regexpr("\\.",tmpbigcols)-1)


tmp_misscols<-cbind(as.vector(tmpbigcols[!tmpbigcols%in%extract_colnames(colnames(mydata_t095$xx_train$fa))]),as.vector(tmpbigcols[!tmpbigcols%in%extract_colnames(colnames(mydata_t095$xx_train$md))]),
                    as.vector(tmpbigcols[!tmpbigcols%in%extract_colnames(colnames(mydata_t095$xx_train$rd))]),as.vector(tmpbigcols[!tmpbigcols%in%extract_colnames(colnames(mydata_t095$xx_train$ad))]),
                    as.vector(tmpbigcols[!tmpbigcols%in%extract_colnames(colnames(mydata_t095$xx_train$nufo))]),as.vector(tmpbigcols[!tmpbigcols%in%extract_colnames(colnames(mydata_t095$xx_train$afdf))]))
misscols<-unique(melt(tmp_misscols)$value)
tmpsmallcols<-tmpbigcols[!tmpbigcols%in%misscols]

# tmpsmallcols<-substr(colnames(as.data.frame(mydata_t09$xx_train[min_mtr])),regexpr("\\.",colnames(as.data.frame(mydata_t09$xx_train[min_mtr])))+1,20)



#Note: I think I'm encountering an issue. With T=1, there's an even bigger mismatch. afdf has 254, nufo has 256. But there's a chance that they're not missing the same conns. so I need 
#to make this tmpsmallcols include missing columns from all metrics. Maybe compare the tmpbigcols against all others. And then just create tmpsmallcols by "subtracting" missing cols from tmpbigcols


tmp_submat<-matrix(data=NA,ncol = length(tmpsmallcols),nrow = length(mydata_t095$yy_train$subjectkey))
colnames(tmp_submat)<-tmpsmallcols
tmp_submat<-as.data.frame(tmp_submat)
tmp_submat$subjectkey<-mydata_t095$yy_train$subjectkey
tmp_submat_melt<-melt(tmp_submat,id.vars = "subjectkey")

#HERE I NEED TO ADAPT THE SCRIPT SO THAT IT CAN HANDLE THE FACT THAT DEPENDING ON SOME THRESHOLDS, THERE WON"T BE A MISMATCH BETWEEN NUFO AND OTHER METRICS
if (length(tmpbigcols)>length(tmpsmallcols)){
  #Note that here I'm removing those "extra" conns that metrics other than NuFO had.
  # mydata_t09$xx_train$all_metrics<-melt(mydata_t09$xx_train$fa[,which(colnames(mydata_t09$xx_train$fa)%in%tmpsmallcols)|
  #                                                                which(colnames(mydata_t09$xx_train$fa)%in%"subjectkey")],id.vars = "subjectkey")
  
  length(mydata_t095$xx_train$fa[,which(extract_colnames(colnames(mydata_t095$xx_train$fa))%in%tmpsmallcols)])
  length(mydata_t095$xx_train$md[,which(extract_colnames(colnames(mydata_t095$xx_train$md))%in%tmpsmallcols)])
  length(mydata_t095$xx_train$rd[,which(extract_colnames(colnames(mydata_t095$xx_train$rd))%in%tmpsmallcols)])
  length(mydata_t095$xx_train$ad[,which(extract_colnames(colnames(mydata_t095$xx_train$ad))%in%tmpsmallcols)])
  length(mydata_t095$xx_train$nufo[,which(extract_colnames(colnames(mydata_t095$xx_train$nufo))%in%tmpsmallcols)])
  length(mydata_t095$xx_train$afdf[,which(extract_colnames(colnames(mydata_t095$xx_train$afd))%in%tmpsmallcols)])
  
  mydata_t095$xx_train$all_metrics<-melt(mydata_t095$xx_train$fa[,which(extract_colnames(colnames(mydata_t095$xx_train$fa))%in%tmpsmallcols)])
  mydata_t095$xx_train$all_metrics$md<-melt(mydata_t095$xx_train$md[,which(extract_colnames(colnames(mydata_t095$xx_train$md))%in%tmpsmallcols)])$value
  mydata_t095$xx_train$all_metrics$rd<-melt(mydata_t095$xx_train$rd[,which(extract_colnames(colnames(mydata_t095$xx_train$rd))%in%tmpsmallcols)])$value
  mydata_t095$xx_train$all_metrics$ad<-melt(mydata_t095$xx_train$ad[,which(extract_colnames(colnames(mydata_t095$xx_train$ad))%in%tmpsmallcols)])$value
  mydata_t095$xx_train$all_metrics$nufo<-melt(mydata_t095$xx_train$nufo[,which(extract_colnames(colnames(mydata_t095$xx_train$nufo))%in%tmpsmallcols)])$value
  mydata_t095$xx_train$all_metrics$afdf<-melt(mydata_t095$xx_train$afdf[,which(extract_colnames(colnames(mydata_t095$xx_train$afdf))%in%tmpsmallcols)])$value
  mydata_t095$xx_train$all_metrics$subjectkey<-tmp_submat_melt$subjectkey
  mydata_t095$xx_train$all_metrics$variable<-tmp_submat_melt$variable
  colnames(mydata_t095$xx_train$all_metrics)<-c("conn","fa","md","rd","ad","nufo","afdf","subjectkey")
  
} else if (length(tmpbigcols)==length(tmpsmallcols)){
  #Note that here I'm removing those "extra" conns that metrics other than NuFO had.
  mydata_t095$xx_train$all_metrics<-melt(mydata_t095$xx_train$fa)
  mydata_t095$xx_train$all_metrics$md<-melt(mydata_t095$xx_train$md)$value
  mydata_t095$xx_train$all_metrics$rd<-melt(mydata_t095$xx_train$rd)$value
  mydata_t095$xx_train$all_metrics$ad<-melt(mydata_t095$xx_train$ad)$value
  mydata_t095$xx_train$all_metrics$nufo<-melt(mydata_t095$xx_train$nufo)$value
  mydata_t095$xx_train$all_metrics$afdf<-melt(mydata_t095$xx_train$afdf)$value
  mydata_t095$xx_train$all_metrics$subjectkey<-tmp_submat_melt$subjectkey
  colnames(mydata_t095$xx_train$all_metrics)<-c("conn","fa","md","rd","ad","nufo","afdf","subjectkey")
}

#Now modify the connection names because they have the metric in them. In all_metrics only .fa appears, so I need to replace that one only
if(regexpr("fa",mydata_t095$xx_train$all_metrics$conn)>1){
  mydata_t095$xx_train$all_metrics$conn<-substr(mydata_t095$xx_train$all_metrics$conn,regexpr("\\.",mydata_t095$xx_train$all_metrics$conn)+1,40)
  mydata_t095$xx_train$all_metrics$conn<-substr(mydata_t095$xx_train$all_metrics$conn,1,regexpr("\\.",mydata_t095$xx_train$all_metrics$conn)-1)
}

####@@PCA@@####

pca.across.metrics.t095<-prcomp(mydata_t095$xx_train$all_metrics[mydata_t095$xx_train$all_metrics$conn%in%selected_conns[selected_conns$c1w_pos==1,c("selected_conns")],
                                                                 -which(colnames(mydata_t095$xx_train$all_metrics)%in%c("subjectkey","conn"))],center = TRUE,scale. = TRUE)

corrplot(pca.across.metrics.t095$rotation, is.corr=FALSE)
summary(pca.across.metrics.t095)
#0.8089 0.1626


####@@EXTRACT FIRST 2 PCs@@####
#Extract first 2 PCS, they account for 97% of variance
PC1.scores.t095<-as.data.frame(pca.across.metrics.t095$x[,1])
colnames(PC1.scores.t095)<-"score"
PC2.scores.t095<-as.data.frame(pca.across.metrics.t095$x[,2])
colnames(PC2.scores.t095)<-"score"
PC1.scores.t095$conns<-mydata_t095$xx_train$all_metrics[mydata_t095$xx_train$all_metrics$conn%in%selected_conns[selected_conns$c1w_pos==1,c("selected_conns")],c("conn")]
PC1.scores.t095$subjectkey<-mydata_t095$xx_train$all_metrics[mydata_t095$xx_train$all_metrics$conn%in%selected_conns[selected_conns$c1w_pos==1,c("selected_conns")],c("subjectkey")]
PC2.scores.t095$conns<-mydata_t095$xx_train$all_metrics[mydata_t095$xx_train$all_metrics$conn%in%selected_conns[selected_conns$c1w_pos==1,c("selected_conns")],c("conn")]
PC2.scores.t095$subjectkey<-mydata_t095$xx_train$all_metrics[mydata_t095$xx_train$all_metrics$conn%in%selected_conns[selected_conns$c1w_pos==1,c("selected_conns")],c("subjectkey")]


####RESTRUCTURE PCS BACK TO WIDE FORMAT####
PC1.scores.wide.t095<-reshape(PC1.scores.t095,idvar = "subjectkey",timevar = "conns",direction = "wide")
colnames(PC1.scores.wide.t095)<-substr(colnames(PC1.scores.wide.t095),regexpr("\\.",colnames(PC1.scores.wide.t095))+1,25)

PC2.scores.wide.t095<-reshape(PC2.scores.t095,idvar = "subjectkey",timevar = "conns",direction = "wide")
colnames(PC2.scores.wide.t095)<-substr(colnames(PC2.scores.wide.t095),regexpr("\\.",colnames(PC2.scores.wide.t095))+1,25)



####PLSC####
plsc.t095.pc1<-plsc_and_permute(PC1.scores.wide.t095[,-which(colnames(PC1.scores.wide.t095)%in%c("subjectkey"))],mydata_t095$yy_train[,-which(colnames(mydata_t095$yy_train)%in%c("subjectkey","sex"))],
                                num_comps = 20,nperms = 2000,topn = 200,out_flag = 1)

plsc.t095.pc2<-plsc_and_permute(PC2.scores.wide.t095[,-which(colnames(PC2.scores.wide.t095)%in%c("subjectkey"))],mydata_t095$yy_train[,-which(colnames(mydata_t095$yy_train)%in%c("subjectkey","sex"))],
                                num_comps = 20,nperms = 2000,topn = 200,out_flag = 1)




####@100%@####
####@@SET UP ALLDATA@@####
tmp_dims<-c(dim(mydata_t1$xx_train$fa)[2],dim(mydata_t1$xx_train$md)[2],dim(mydata_t1$xx_train$rd)[2],dim(mydata_t1$xx_train$ad)[2],dim(mydata_t1$xx_train$nufo)[2],dim(mydata_t1$xx_train$afdf)[2])
tmp_mtrs<-c("fa","md","rd","ad","nufo","afdf")
tmp_dim_mtrs<-rbind(tmp_dims,tmp_mtrs)
min_mtr<-tmp_dim_mtrs[rownames(tmp_dim_mtrs)=="tmp_mtrs",which(tmp_dim_mtrs[1,]==min(tmp_dim_mtrs[1,]))][1]
max_mtr<-tmp_dim_mtrs[rownames(tmp_dim_mtrs)=="tmp_mtrs",which(tmp_dim_mtrs[1,]==max(tmp_dim_mtrs[1,]))][1]
# tmpbigcols<-extract_colnames(input_colnames = colnames(as.data.frame(mydata_t09$xx_train[max_mtr])))

tmpbigcols<-substr(colnames(as.data.frame(mydata_t1$xx_train[max_mtr])),regexpr("\\.",colnames(as.data.frame(mydata_t1$xx_train[max_mtr])))+1,20)#for some reason, this line doesn't extract the conn number alone anymore
tmpbigcols<-extract_colnames(tmpbigcols)

# tmpbigcols<-substr(tmpbigcols,regexpr("\\.",tmpbigcols)+1,20)
# tmpbigcols<-substr(tmpbigcols,1,regexpr("\\.",tmpbigcols)-1)


tmp_misscols<-cbind(as.vector(tmpbigcols[!tmpbigcols%in%extract_colnames(colnames(mydata_t1$xx_train$fa))]),as.vector(tmpbigcols[!tmpbigcols%in%extract_colnames(colnames(mydata_t1$xx_train$md))]),
                    as.vector(tmpbigcols[!tmpbigcols%in%extract_colnames(colnames(mydata_t1$xx_train$rd))]),as.vector(tmpbigcols[!tmpbigcols%in%extract_colnames(colnames(mydata_t1$xx_train$ad))]),
                    as.vector(tmpbigcols[!tmpbigcols%in%extract_colnames(colnames(mydata_t1$xx_train$nufo))]),as.vector(tmpbigcols[!tmpbigcols%in%extract_colnames(colnames(mydata_t1$xx_train$afdf))]))
misscols<-unique(melt(tmp_misscols)$value)
tmpsmallcols<-tmpbigcols[!tmpbigcols%in%misscols]

# tmpsmallcols<-substr(colnames(as.data.frame(mydata_t09$xx_train[min_mtr])),regexpr("\\.",colnames(as.data.frame(mydata_t09$xx_train[min_mtr])))+1,20)



#Note: I think I'm encountering an issue. With T=1, there's an even bigger mismatch. afdf has 254, nufo has 256. But there's a chance that they're not missing the same conns. so I need 
#to make this tmpsmallcols include missing columns from all metrics. Maybe compare the tmpbigcols against all others. And then just create tmpsmallcols by "subtracting" missing cols from tmpbigcols


tmp_submat<-matrix(data=NA,ncol = length(tmpsmallcols),nrow = length(mydata_t1$yy_train$subjectkey))
colnames(tmp_submat)<-tmpsmallcols
tmp_submat<-as.data.frame(tmp_submat)
tmp_submat$subjectkey<-mydata_t1$yy_train$subjectkey
tmp_submat_melt<-melt(tmp_submat,id.vars = "subjectkey")

#HERE I NEED TO ADAPT THE SCRIPT SO THAT IT CAN HANDLE THE FACT THAT DEPENDING ON SOME THRESHOLDS, THERE WON"T BE A MISMATCH BETWEEN NUFO AND OTHER METRICS
if (length(tmpbigcols)>length(tmpsmallcols)){
  #Note that here I'm removing those "extra" conns that metrics other than NuFO had.
  # mydata_t09$xx_train$all_metrics<-melt(mydata_t09$xx_train$fa[,which(colnames(mydata_t09$xx_train$fa)%in%tmpsmallcols)|
  #                                                                which(colnames(mydata_t09$xx_train$fa)%in%"subjectkey")],id.vars = "subjectkey")
  
  length(mydata_t1$xx_train$fa[,which(extract_colnames(colnames(mydata_t1$xx_train$fa))%in%tmpsmallcols)])
  length(mydata_t1$xx_train$md[,which(extract_colnames(colnames(mydata_t1$xx_train$md))%in%tmpsmallcols)])
  length(mydata_t1$xx_train$rd[,which(extract_colnames(colnames(mydata_t1$xx_train$rd))%in%tmpsmallcols)])
  length(mydata_t1$xx_train$ad[,which(extract_colnames(colnames(mydata_t1$xx_train$ad))%in%tmpsmallcols)])
  length(mydata_t1$xx_train$nufo[,which(extract_colnames(colnames(mydata_t1$xx_train$nufo))%in%tmpsmallcols)])
  length(mydata_t1$xx_train$afdf[,which(extract_colnames(colnames(mydata_t1$xx_train$afd))%in%tmpsmallcols)])
  
  mydata_t1$xx_train$all_metrics<-melt(mydata_t1$xx_train$fa[,which(extract_colnames(colnames(mydata_t1$xx_train$fa))%in%tmpsmallcols)])
  mydata_t1$xx_train$all_metrics$md<-melt(mydata_t1$xx_train$md[,which(extract_colnames(colnames(mydata_t1$xx_train$md))%in%tmpsmallcols)])$value
  mydata_t1$xx_train$all_metrics$rd<-melt(mydata_t1$xx_train$rd[,which(extract_colnames(colnames(mydata_t1$xx_train$rd))%in%tmpsmallcols)])$value
  mydata_t1$xx_train$all_metrics$ad<-melt(mydata_t1$xx_train$ad[,which(extract_colnames(colnames(mydata_t1$xx_train$ad))%in%tmpsmallcols)])$value
  mydata_t1$xx_train$all_metrics$nufo<-melt(mydata_t1$xx_train$nufo[,which(extract_colnames(colnames(mydata_t1$xx_train$nufo))%in%tmpsmallcols)])$value
  mydata_t1$xx_train$all_metrics$afdf<-melt(mydata_t1$xx_train$afdf[,which(extract_colnames(colnames(mydata_t1$xx_train$afdf))%in%tmpsmallcols)])$value
  mydata_t1$xx_train$all_metrics$subjectkey<-tmp_submat_melt$subjectkey
  mydata_t1$xx_train$all_metrics$variable<-tmp_submat_melt$variable
  colnames(mydata_t1$xx_train$all_metrics)<-c("conn","fa","md","rd","ad","nufo","afdf","subjectkey")
  
} else if (length(tmpbigcols)==length(tmpsmallcols)){
  #Note that here I'm removing those "extra" conns that metrics other than NuFO had.
  mydata_t1$xx_train$all_metrics<-melt(mydata_t1$xx_train$fa)
  mydata_t1$xx_train$all_metrics$md<-melt(mydata_t1$xx_train$md)$value
  mydata_t1$xx_train$all_metrics$rd<-melt(mydata_t1$xx_train$rd)$value
  mydata_t1$xx_train$all_metrics$ad<-melt(mydata_t1$xx_train$ad)$value
  mydata_t1$xx_train$all_metrics$nufo<-melt(mydata_t1$xx_train$nufo)$value
  mydata_t1$xx_train$all_metrics$afdf<-melt(mydata_t1$xx_train$afdf)$value
  mydata_t1$xx_train$all_metrics$subjectkey<-tmp_submat_melt$subjectkey
  colnames(mydata_t1$xx_train$all_metrics)<-c("conn","fa","md","rd","ad","nufo","afdf","subjectkey")
}

#Now modify the connection names because they have the metric in them. In all_metrics only .fa appears, so I need to replace that one only
if(regexpr("fa",mydata_t1$xx_train$all_metrics$conn)>1){
  mydata_t1$xx_train$all_metrics$conn<-substr(mydata_t1$xx_train$all_metrics$conn,regexpr("\\.",mydata_t1$xx_train$all_metrics$conn)+1,40)
  mydata_t1$xx_train$all_metrics$conn<-substr(mydata_t1$xx_train$all_metrics$conn,1,regexpr("\\.",mydata_t1$xx_train$all_metrics$conn)-1)
}

####@@PCA@@####

pca.across.metrics.t1<-prcomp(mydata_t1$xx_train$all_metrics[mydata_t1$xx_train$all_metrics$conn%in%selected_conns[selected_conns$c1w_pos==1,c("selected_conns")],
                                                                 -which(colnames(mydata_t1$xx_train$all_metrics)%in%c("subjectkey","conn"))],center = TRUE,scale. = TRUE)

corrplot(pca.across.metrics.t1$rotation, is.corr=FALSE)
summary(pca.across.metrics.t1)
#0.784 0.175


####@@EXTRACT FIRST 2 PCs@@####
#Extract first 2 PCS, they account for 97% of variance
PC1.scores.t1<-as.data.frame(pca.across.metrics.t1$x[,1])
colnames(PC1.scores.t1)<-"score"
PC2.scores.t1<-as.data.frame(pca.across.metrics.t1$x[,2])
colnames(PC2.scores.t1)<-"score"
PC1.scores.t1$conns<-mydata_t1$xx_train$all_metrics[mydata_t1$xx_train$all_metrics$conn%in%selected_conns[selected_conns$c1w_pos==1,c("selected_conns")],c("conn")]
PC1.scores.t1$subjectkey<-mydata_t1$xx_train$all_metrics[mydata_t1$xx_train$all_metrics$conn%in%selected_conns[selected_conns$c1w_pos==1,c("selected_conns")],c("subjectkey")]
PC2.scores.t1$conns<-mydata_t1$xx_train$all_metrics[mydata_t1$xx_train$all_metrics$conn%in%selected_conns[selected_conns$c1w_pos==1,c("selected_conns")],c("conn")]
PC2.scores.t1$subjectkey<-mydata_t1$xx_train$all_metrics[mydata_t1$xx_train$all_metrics$conn%in%selected_conns[selected_conns$c1w_pos==1,c("selected_conns")],c("subjectkey")]


####RESTRUCTURE PCS BACK TO WIDE FORMAT####
PC1.scores.wide.t1<-reshape(PC1.scores.t1,idvar = "subjectkey",timevar = "conns",direction = "wide")
colnames(PC1.scores.wide.t1)<-substr(colnames(PC1.scores.wide.t1),regexpr("\\.",colnames(PC1.scores.wide.t1))+1,25)

PC2.scores.wide.t1<-reshape(PC2.scores.t1,idvar = "subjectkey",timevar = "conns",direction = "wide")
colnames(PC2.scores.wide.t1)<-substr(colnames(PC2.scores.wide.t1),regexpr("\\.",colnames(PC2.scores.wide.t1))+1,25)

####PLSC####
plsc.t1.pc1<-plsc_and_permute(PC1.scores.wide.t1[,-which(colnames(PC1.scores.wide.t1)%in%c("subjectkey"))],mydata_t1$yy_train[,-which(colnames(mydata_t1$yy_train)%in%c("subjectkey","sex"))],
                                num_comps = 20,nperms = 2000,topn = 200,out_flag = 1)

plsc.t1.pc2<-plsc_and_permute(PC2.scores.wide.t1[,-which(colnames(PC2.scores.wide.t1)%in%c("subjectkey"))],mydata_t1$yy_train[,-which(colnames(mydata_t1$yy_train)%in%c("subjectkey","sex"))],
                                num_comps = 20,nperms = 2000,topn = 200,out_flag = 1)






##WAIT: uvar feat sel is performed after PCA. BUT, COMMIT is performed before. That doesn't guarantee they'll have the same conns, if thresholding removed some of them.
#JUST MAKE SURE THE CORRELATION BETWEEN PROJECTIONS IS NOT 100%. THERE SHOULD BE SOME SMALL TINY DIFFERENCES AT LEAST
dim(mydata_t085$xx_train$fa)
dim(mydata_t09$xx_train$fa)
dim(mydata_t095$xx_train$fa)
dim(mydata_t1$xx_train$fa)#258

dim(pca.across.metrics.t085$x)
dim(pca.across.metrics$x)
dim(pca.across.metrics.t095$x)
dim(pca.across.metrics.t1$x)
#t1 is the only one that lost some conns. Which is why it's really the ONLY one that is not identical.
#Comparisons between thresholds are unnecesary therefore. The only thing worth reporting is the fact that except for t1, all other connections selected were the same, which means that 
#there seems to be some agreement between connections that exist across most people and COMMIT weights. It's reassuring to see that COMMIT filtering is mostly agreeing with thresholding.
#Remember that COMMIT comes after thresholding. So you run COMMIT on whatever is left. But what is reassuring is that uvarfeatsel selected most of the connections that would go on to survive COMMIT. 

#ARE THERE DIFFERENCES IN THE CONNECTIONS SELECTED IN BOTH ANALYSES?

conns.t085<-unique(PC1.scores.t085$conns)
conns.t09<-unique(PC1.scores$conns)
conns.t095<-unique(PC1.scores.t095$conns)
conns.t1<-unique(PC1.scores.t1$conns)

length(conns.t1)#252
#So there were only 6 connections that were selected by thresholding that were not selected by COMMIT. 

sum(conns.t085%in%conns.t09)
sum(conns.t09%in%conns.t095)
sum(conns.t09%in%conns.t1)



###Since T1 is the only threshold that yields different connections, I'll only proceed with the T1 connections
###BootstrapTests###
boot.pc1.t1<-Boot4PLSC(PC1.scores.wide.t1[,which(colnames(PC1.scores.wide.t1)%in%c(rownames(plsc.t1.pc1[[4]]$TExPosition.Data$pdq$p)))],mydata_t1$yy_train[,-which(colnames(mydata_t1$yy_train)%in%c("subjectkey","sex"))],
                         center1 = TRUE,center2 = TRUE,scale1 = TRUE,scale2 = TRUE,nf2keep = 19,nIter = 2000,critical.value = 1.96,eig = TRUE)
# boot.pc1.t1$bootRatios.j[, plsc.t1.pc1[[1]]<.05]
# boot.pc1.t1$bootRatiosSignificant.j[, plsc.t1.pc1[[1]]<.05]


###Illustrate all sig LFs###
plot_data<-plsc.t1.pc1
boot_data<-boot.pc1.t1

###APPLIES TO ALL POLAR PLOTS###
#Store loadings
df<-as.data.frame(plot_data[[4]]$TExPosition.Data$pdq$q)
#Change column names
colnames(df)<-paste("Y",substr(colnames(df),2,5),sep="_")
#Only keep the lfs that were significant
# df<-df[,plot_data[[1]]<.05]###NOTE: IM PLOTTING EVEN THE ONES THAT ARE NONSIG, TO COMPARE AGAINST T90. BUT I HAVE TO SPECIFY THAT THEY'RE NOT ALL SIG.
#Store symptoms as a new column
df$symptom<-rownames(df)
#Create a column to specify symptom domains (not really used atm)
df$domain<-c("mood","mood","cognitive","mood","somatic","somatic","somatic","somatic","somatic","somatic","somatic","somatic","somatic","cognitive","cognitive","cognitive","cognitive","cognitive","cognitive")
#Create a column for the common name of the symptom
df$gen_prob<-c("depression","anxiety","attention","aggression","headache","nausea","vomiting","dizziness","fatigue","sleep+","sleep-","drowsiness","troublesleeping","s-recall","l-recall","seq-memory","work-memory","proc-speed","cardsorting")
#Flag significant symptoms based on bootstrap tests
# df$boottest<-boot_data$bootRatiosSignificant.j[,plot_data[[1]]<.05]
df$boottest<-boot_data$bootRatiosSignificant.j
#Melt df
df_melt<-melt(df)

#Create color palette
c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "sienna", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)


#Create df for labels
label_data <- as.data.frame(sort(df$gen_prob))
colnames(label_data)<-"gen_prob"
label_data$id<-seq(1,19)
# calculate the ANGLE of the labels
number_of_bar <- nrow(label_data)
angle <-  90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)

###@@@###

#Note: will reverse a few vars for visualization purposes
df_melt[df_melt$variable=="Y_1",c("value")]<-df_melt[df_melt$variable=="Y_1",c("value")]*-1
df_melt[df_melt$variable=="Y_2",c("value")]<-df_melt[df_melt$variable=="Y_2",c("value")]*-1
df_melt[df_melt$variable=="Y_3",c("value")]<-df_melt[df_melt$variable=="Y_3",c("value")]*-1
df_melt[df_melt$variable=="Y_7",c("value")]<-df_melt[df_melt$variable=="Y_7",c("value")]*-1


###SPECIFIC POLAR PLOTS###
p1<-create_polar_plots(df_melt,var = "Y_1",color = "black",colorpal = c25,ylims = c(-0.5,0.8),boot=TRUE,label_data)
p2<-create_polar_plots(df_melt,var = "Y_2",color = "black",colorpal = c25,ylims = c(-0.9,0.9),boot=TRUE,label_data)
p3<-create_polar_plots(df_melt,var = "Y_3",color = "black",colorpal = c25,ylims = c(-0.6,0.8),boot=TRUE,label_data)
# p4<-create_polar_plots(df_melt,var = "Y_4",color = "black",colorpal = c25,ylims = c(-0.7,0.8),boot=TRUE,label_data)
p5<-create_polar_plots(df_melt,var = "Y_5",color = "black",colorpal = c25,ylims = c(-0.6,0.8),boot=TRUE,label_data)
# p6<-create_polar_plots(df_melt,var = "Y_6",color = "black",colorpal = c25,ylims = c(-0.7,0.8),boot=TRUE,label_data)
p7<-create_polar_plots(df_melt,var = "Y_7",color = "black",colorpal = c25,ylims = c(-0.8,0.8),boot=TRUE,label_data)
p8<-create_polar_plots(df_melt,var = "Y_8",color = "black",colorpal = c25,ylims = c(-0.6,0.8),boot=TRUE,label_data)
p9<-create_polar_plots(df_melt,var = "Y_9",color = "black",colorpal = c25,ylims = c(-0.8,0.8),boot=TRUE,label_data)
p10<-create_polar_plots(df_melt,var = "Y_10",color = "black",colorpal = c25,ylims = c(-0.7,0.8),boot=TRUE,label_data)
p11<-create_polar_plots(df_melt,var = "Y_11",color = "black",colorpal = c25,ylims = c(-0.6,0.8),boot=TRUE,label_data)
p12<-create_polar_plots(df_melt,var = "Y_12",color = "black",colorpal = c25,ylims = c(-0.6,0.8),boot=TRUE,label_data)
p13<-create_polar_plots(df_melt,var = "Y_13",color = "black",colorpal = c25,ylims = c(-0.8,0.8),boot=TRUE,label_data)
p14<-create_polar_plots(df_melt,var = "Y_14",color = "black",colorpal = c25,ylims = c(-0.8,0.8),boot=TRUE,label_data)
p15<-create_polar_plots(df_melt,var = "Y_15",color = "black",colorpal = c25,ylims = c(-0.8,0.8),boot=TRUE,label_data)
p16<-create_polar_plots(df_melt,var = "Y_16",color = "black",colorpal = c25,ylims = c(-0.8,0.9),boot=TRUE,label_data)
p17<-create_polar_plots(df_melt,var = "Y_17",color = "black",colorpal = c25,ylims = c(-0.9,0.9),boot=TRUE,label_data)
p18<-create_polar_plots(df_melt,var = "Y_18",color = "black",colorpal = c25,ylims = c(-0.8,0.8),boot=TRUE,label_data)
# p19<-create_polar_plots(df_melt,var = "Y_19",color = "black",colorpal = c25,ylims = c(-0.8,0.8),boot=TRUE,label_data)

plot_grid(p1,p2,p3,p5,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18)


###Comparisons between thresholds###
#Comparisons are unnecessary, except against t1. It's the only one I'll do. I'll compare all others against t1. 

t085_lx<-as.data.frame(plsc.t085.pc1[[4]]$TExPosition.Data$lx)
colnames(t085_lx)<-paste("t085-",substr(colnames(t085_lx),2,5),sep="")  

t09_lx<-as.data.frame(plsc.t09.pc1[[4]]$TExPosition.Data$lx)
colnames(t09_lx)<-paste("t090-",substr(colnames(t09_lx),2,5),sep="")  

t095_lx<-as.data.frame(plsc.t095.pc1[[4]]$TExPosition.Data$lx)
colnames(t095_lx)<-paste("t095-",substr(colnames(t095_lx),2,5),sep="")  

t1_lx<-as.data.frame(plsc.t1.pc1[[4]]$TExPosition.Data$lx)
colnames(t1_lx)<-paste("t1-",substr(colnames(t1_lx),2,5),sep="")  


lxs<-cbind(t09_lx,t1_lx)
corrplot(cor(lxs),method="color",type="lower",tl.col = "black",tl.pos = "l")

