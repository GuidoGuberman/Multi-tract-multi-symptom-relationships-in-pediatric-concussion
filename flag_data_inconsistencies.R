###########################################################
##        SCRIPT TO HANDLE DATA INCONSISTENCIES          ##
###########################################################

##############
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
#############


# library(gmp)
args = commandArgs(trailingOnly=TRUE)
library(MatchIt)
library(dplyr)


##########Handle Philips rev_b0 problem by computing mutual information######
image_info_BIDS_noMI <- read.table(paste(args[1],'/Data/BIDS/image_info_BIDS.txt',sep=''), sep='', header=T)
# image_info_BIDS_noMI <- read.table('/Users/Guido/Desktop/GeneralData/ABCD/Harmonization_Controls/Data/BIDS/image_info_BIDS.txt', sep='', header=T)#For troubleshooting
# image_info_BIDS_noMI <- read.table('/Users/Guido/Desktop/GeneralData/ABCD/mTBI_Group/Data/BIDS/image_info_BIDS.txt', sep='', header=T)#For troubleshooting
# image_info_BIDS_noMI <- read.table('/Users/Guido/Desktop/Projects/TBI-ConductDisordersStudy/Data/mTBI_Group/Data/BIDS/image_info_BIDS.txt', sep='', header=T)#For troubleshooting
grp=args[2]####GROUP (Defined in Master Script)#######
mutual_info <- read.table(paste(args[1],'/Data/BIDS/extra/mutual_info_Philips.txt',sep=''), sep='', header=F)
# mutual_info <- read.table('/Users/Guido/Desktop/GeneralData/ABCD/Harmonization_Controls/Data/BIDS/extra/mutual_info_Philips.txt', sep='', header=F)#For troubleshooting
# mutual_info <- read.table('/Users/Guido/Desktop/GeneralData/ABCD/mTBI_Group/Data/BIDS/extra/mutual_info_Philips.txt', sep='', header=F)#For troubleshooting
# mutual_info <- read.table('/Users/Guido/Desktop/Projects/TBI-ConductDisordersStudy/Data/mTBI_Group/Data/BIDS/extra/mutual_info_Philips.txt', sep='', header=F)#For troubleshooting

#Name the columns in mutual info file
names(mutual_info)[names(mutual_info) == "V1"] <- "subjectkey"
names(mutual_info)[names(mutual_info) == "V2"] <- "FileName"
names(mutual_info)[names(mutual_info) == "V3"] <- "ProtocolName"
names(mutual_info)[names(mutual_info) == "V4"] <- "MutualInformation"

#NOTE: There may have been a change in the json files provided by ABCD, because I can't find PatientName anymore. I have to change my script because the lack of PatientName is the reason image_info_BIDS_noMI is coming out empty

#Change the name of the files to refer to json files
mutual_info$FileName<-switch_names(mutual_info$FileName,".nii.gz",".json",'',"R")


#Change the way we refer to subs
# mutual_info$subjectkey<-switch_names(mutual_info$subjectkey,"sub-NDAR","NDAR",'_',"L")
image_info_BIDS<-merge(image_info_BIDS_noMI,mutual_info,all=T)

# dim(image_info_BIDS_noMI)
# dim(mutual_info)
# dim(image_info_BIDS)
# Correct PE direction for Philips (remember: Philips scans don't have PE Direction, just PE Axis)
for (i in 1:dim(image_info_BIDS)[1]){
  # Find all scans for a given person
  tmp1<-subset(image_info_BIDS,subjectkey==image_info_BIDS$subjectkey[i])
  
  # If it's a Philips scan
  if (tmp1$Manufacturer[1]=="Philips"){
    
    #Find b0s
    b0s<-subset(tmp1,tmp1$ProtocolName=='ABCD_dMRI_DistortionMap_AP'|tmp1$ProtocolName=='ABCD_dMRI_DistortionMap_PA'|
                  tmp1$ProtocolName=='ABCD-Diffusion-FM_GE_original_(baseline_year_1_arm_1)'|tmp1$ProtocolName=='DTI_Fieldmap_A'|
                  tmp1$ProtocolName=='DTI_Fieldmap_P'|tmp1$ProtocolName=='WIP_DTI_Fieldmap_A'|tmp1$ProtocolName=='WIP_DTI_Fieldmap_P')
    ###Identify rev_b0 as those that have the lowest mutual information###
    rev_b0<-b0s[b0s$MutualInformation==min(b0s$MutualInformation),]
    
    #In the column for the PED, add a j- to signify a rev_b0 image
    image_info_BIDS[row.names(rev_b0),c('PhaseEncodingDirection')] <- 'j-'
  }
}


# Flag b0s that are not rev_b0
for (i in 1:dim(image_info_BIDS)[1]){
    tmp1<-subset(image_info_BIDS,subjectkey==image_info_BIDS$subjectkey[i])
    dwi<-subset(tmp1,tmp1$ProtocolName=='ABCD_dMRI'|tmp1$ProtocolName=='ABCD-DTI_GE_original_(baseline_year_1_arm_1)'|tmp1$ProtocolName=='DTI_1'|
                     tmp1$ProtocolName=='DTI_2'|tmp1$ProtocolName=='WIP_DTI_1'|tmp1$ProtocolName=='WIP_DTI_2')
    #Identify b0s
    b0s<-subset(tmp1,tmp1$ProtocolName=='ABCD_dMRI_DistortionMap_AP'|tmp1$ProtocolName=='ABCD_dMRI_DistortionMap_PA'|
                  tmp1$ProtocolName=='ABCD-Diffusion-FM_GE_original_(baseline_year_1_arm_1)'|tmp1$ProtocolName=='DTI_Fieldmap_A'|
                  tmp1$ProtocolName=='DTI_Fieldmap_P'|tmp1$ProtocolName=='WIP_DTI_Fieldmap_A'|tmp1$ProtocolName=='WIP_DTI_Fieldmap_P')
    #Find b0s whose PED is equal to the PED of DWI images (i.e.: they're not rev_b0)
    non_rev_b0<-subset(b0s,b0s$PhaseEncodingDirection==unique(dwi$PhaseEncodingDirection))
    
    image_info_BIDS[row.names(non_rev_b0),c('non_rev_b0')]<-1
}

# Flag subjects without any rev_b0
for (i in 1:dim(image_info_BIDS)[1]){
  tmp1<-subset(image_info_BIDS,subjectkey==image_info_BIDS$subjectkey[i])
  dwi<-subset(tmp1,tmp1$ProtocolName=='ABCD_dMRI'|tmp1$ProtocolName=='ABCD-DTI_GE_original_(baseline_year_1_arm_1)'|tmp1$ProtocolName=='DTI_1'|
                tmp1$ProtocolName=='DTI_2'|tmp1$ProtocolName=='WIP_DTI_1'|tmp1$ProtocolName=='WIP_DTI_2')
  
  b0s<-subset(tmp1,tmp1$ProtocolName=='ABCD_dMRI_DistortionMap_AP'|tmp1$ProtocolName=='ABCD_dMRI_DistortionMap_PA'|
                tmp1$ProtocolName=='ABCD-Diffusion-FM_GE_original_(baseline_year_1_arm_1)'|tmp1$ProtocolName=='DTI_Fieldmap_A'|
                tmp1$ProtocolName=='DTI_Fieldmap_P'|tmp1$ProtocolName=='WIP_DTI_Fieldmap_A'|tmp1$ProtocolName=='WIP_DTI_Fieldmap_P')
  
  #Identify the rev_b0 by finding images that have PED that is different from the DWI's PED
  rev_b0<-subset(b0s,b0s$PhaseEncodingDirection!=unique(dwi$PhaseEncodingDirection))
  if (dim(rev_b0)[1]==0){
    #If there are no rev_b0s
    image_info_BIDS[row.names(tmp1),c('missing_rev_b0')]<-1
  }
}

tmp<-subset(image_info_BIDS,image_info_BIDS$missing_rev_b0==1) 
subs_missing_rev_b0<-data.frame(unique(tmp$subjectkey))
names(subs_missing_rev_b0) <- "subjectkey"

#dim(subs_missing_rev_b0)[1]#mTBI: 0 subs with a missing rev_b0

#Remove subjects without any rev_b0
image_info_BIDS_norev0miss_pre <- image_info_BIDS[! image_info_BIDS$subjectkey %in% subs_missing_rev_b0$subjectkey, ]
# dim(data.frame(unique(image_info_BIDS_norev0miss_pre$subjectkey)))#mTBI=351 

##############INTEGRATE INFO TO FIND TENSOR PEAK ERRORS################
file_dims <- read.table(paste(args[1],'/Data/BIDS/extra/verifying_file_dimensions.txt',sep=''), sep='', header=T)
# file_dims <- read.table('/Users/Guido/Desktop/GeneralData/ABCD/mTBI_Group/Data/BIDS/extra/verifying_file_dimensions.txt', sep='', header=T)#For troubleshooting
# file_dims <- read.table('/Users/Guido/Desktop/GeneralData/ABCD/Harmonization_Controls/Data/BIDS/extra/verifying_file_dimensions.txt', sep='', header=T)#For troubleshooting
# file_dims <- read.table('/Users/Guido/Desktop/Projects/TBI-ConductDisordersStudy/Data/mTBI_Group/Data/BIDS/extra/verifying_file_dimensions.txt', sep='', header=T)#For troubleshooting
fmaps<-subset(file_dims,file_dims$ImageType=="fmap")

#b0s are a single volume so should only have 3dims
errors<-subset(fmaps,fmaps$Obtained_Dim>3)

group_n<-dim(data.frame(unique(errors$subjectkey)))[1]
cat('There are',group_n,'subjects with Tensor Peak errors ','\n')
#NOTE: Having 4Ds does not necessarily guarantee that it's a tensor-peak error, but it definitely means something is wrong with the fmap

#add a nifti line to the larger list.
image_info_BIDS_norev0miss_pre$FileName_nifti<-switch_names(image_info_BIDS_norev0miss_pre$FileName,".json",".nii.gz",'',"R")

# dim(data.frame(unique(image_info_BIDS_norev0miss_pre$subjectkey)))

image_info_BIDS_norev0miss<-image_info_BIDS_norev0miss_pre[! image_info_BIDS_norev0miss_pre$FileName_nifti %in% errors$FileName_nifti, ]


image_info_BIDS_norev0miss$subjectkey<-switch_names(image_info_BIDS_norev0miss$subjectkey,"NDAR_","sub-NDAR",'',"L")

image_info_BIDS_norev0miss[ image_info_BIDS_norev0miss$subjectkey %in% errors$subjectkey,]#after removing images with dimension inconsistencies, which subjects become "subjects with missing data"?

# dim(data.frame(unique(image_info_BIDS_norev0miss$subjectkey)))
#Here I need to once again check the number of people who have missing data.

# unique(image_info_BIDS_norev0miss$ProtocolName)




#Now that we've removed rows corresponding to images that had missing QC scores or had failed QC, we will count again which subjects have complete data.
tmp <- image_info_BIDS_norev0miss %>% group_by(subjectkey) %>% summarise(DTI_FM_verif = sum(ProtocolName == 'ABCD-Diffusion-FM_GE_original_(baseline_year_1_arm_1)'), 
                                                                         DTI_FM_PA_verif = sum(ProtocolName == 'DTI_Fieldmap_P'|ProtocolName == 'WIP_DTI_Fieldmap_P'|ProtocolName == 'ABCD_dMRI_DistortionMap_PA'),
                                                                         DTI_FM_AP_verif = sum(ProtocolName == 'DTI_Fieldmap_A'|ProtocolName == 'WIP_DTI_Fieldmap_A'|ProtocolName == 'ABCD_dMRI_DistortionMap_AP'),
                                                                         ABCD_DTI_verif = sum(ProtocolName == 'DTI_1'|ProtocolName == 'ABCD_dMRI'|ProtocolName == 'ABCD-DTI_GE_original_(baseline_year_1_arm_1)'|ProtocolName == 'WIP_DTI_1'), 
                                                                         ABCD_T1_verif = sum(ProtocolName =='WIP_3D_T1'|ProtocolName =='ABCD_T1w_MPR_vNav'|ProtocolName =='ABCD-T1_GE_original_(baseline_year_1_arm_1)'|ProtocolName =='ABCD_T1w_MPR_vNav_repeat'))
Verif <- merge(tmp, image_info_BIDS_norev0miss, by='subjectkey', all=T)
# dim(tmp)

PHILIPS_verif <- subset(Verif, Manufacturer=='Philips')#I've removed the irrelevant b0s so now Philips shouldn't have both. I've also concatenated DWIs
PHILIPS_verif$missingdata_verif <- ifelse(PHILIPS_verif$ABCD_DTI_verif<1|PHILIPS_verif$DTI_FM_PA_verif<1|PHILIPS_verif$DTI_FM_AP_verif<1|PHILIPS_verif$ABCD_T1_verif<1, 1,0)
PHILIPS_verif$T1REPEAT_verif <- ifelse(PHILIPS_verif$ABCD_T1_verif>1, 1, 0)
PHILIPS_verif$DWIREPEAT_verif <- ifelse(PHILIPS_verif$ABCD_DTI_verif>2 | PHILIPS_verif$DTI_FM_PA_verif>1| PHILIPS_verif$DTI_FM_AP_verif>1, 1,0)

SIEMENS_verif <- subset(Verif, Manufacturer=='Siemens')
SIEMENS_verif$missingdata_verif <- ifelse(SIEMENS_verif$ABCD_DTI_verif<1|SIEMENS_verif$DTI_FM_PA_verif<1|SIEMENS_verif$DTI_FM_AP_verif<1|SIEMENS_verif$ABCD_T1_verif<1, 1,0)
SIEMENS_verif$T1REPEAT_verif <- ifelse(SIEMENS_verif$ABCD_T1_verif>1, 1, 0)
SIEMENS_verif$DWIREPEAT_verif <- ifelse(SIEMENS_verif$ABCD_DTI_verif>1 | SIEMENS_verif$DTI_FM_PA_verif| SIEMENS_verif$DTI_FM_AP_verif>1, 1,0)

GE_verif <- subset(Verif, Manufacturer=='GE')
GE_verif$missingdata_verif <- ifelse(GE_verif$ABCD_DTI_verif<1|GE_verif$DTI_FM_verif<1|GE_verif$ABCD_T1_verif<1, 1,0)
GE_verif$T1REPEAT_verif <- ifelse(GE_verif$ABCD_T1_verif>1, 1, 0)
GE_verif$DWIREPEAT_verif <- ifelse(GE_verif$ABCD_DTI_verif>1 | GE_verif$DTI_FM_verif>1, 1,0)



FINAL_inconsistencies_verif <- rbind(GE_verif,SIEMENS_verif,PHILIPS_verif) #Here's the important table for the next steps.

subs_missingdata_verif<-data.frame(unique(FINAL_inconsistencies_verif$subjectkey[FINAL_inconsistencies_verif$missingdata_verif==1]))
colnames(subs_missingdata_verif)<-'subjectkey'

FINAL_data_clean<-image_info_BIDS_norev0miss[! image_info_BIDS_norev0miss$subjectkey %in% subs_missingdata_verif$subjectkey, ]

group_nn<-dim(data.frame(subs_missingdata_verif))[1]#should be 350 (mtbi) or 460 (Match Controls)
group_n<-dim(data.frame(unique(FINAL_data_clean$subjectkey)))[1]#mTBI: 339 (consistent with what I had before) MC: 433 subs (460-27)
cat('There are',group_n,'subjects who have complete data after removing potential tensor peak errors (',group_nn,'now became subjects with missing data)','\n')



################MAIN LOOP TO GET THE FINAL LIST OF SUBJECTS AND SCANS####################
#An iterative loop that recovers all these cases + solves the repeats. It is based on a hierarchical procedure.
#The loop tries finding as much complete data as it can within a session, and whatever it can't find, it looks for in the next available
#session. It prioritizes getting data that is as temporally close as possible. But repeats within a session are handled by selecting the
#last acquired (as per ABCD's recommendation).
#1. Look within a session, select the last-acquired scans if data is repeated, and look for complete data. 2. If data is not complete,
#look for the next available session. Apply the same within-session criterion (select last-acquired scan if data is repeated) and see if you
#can complete data. Prioritize having the most data come from the same session (i.e.: if both scans have a T1, prioritize the T1 scan that
#comes from the same session as most of the other data (maybe minimize the differences in sessionUID))

image_info_BIDS_norev0miss<-FINAL_data_clean
image_info_BIDS_norev0miss[,c('complete_data_same_session_ASR')]<-NA
#ASR=After Solving for Repeats
subs<-data.frame(unique(image_info_BIDS_norev0miss$subjectkey))
colnames(subs)<-'subjectkey'
for (i in 1:dim(subs)[1]){
  tmp1<-subset(image_info_BIDS_norev0miss,subjectkey==subs$subjectkey[i])#Get data for a given person
  tmp2<-data.frame(unique(tmp1$StudyInstanceUID)) #Then take all study UIDs (unique for each session)
  for (l in 1:dim(tmp2)[1]){
    tmp1<-subset(image_info_BIDS_norev0miss,subjectkey==subs$subjectkey[i])
    #redefine tmp1 within this loop, because tmp1 should change as this loop runs
    
    comp_win_sess_dat_check<-subset(tmp1, tmp1$complete_data_same_session_ASR==1) #a variable to verify if the person already has been found to have complete data in a session

    if (dim(comp_win_sess_dat_check)[1]==0){
      tmp3<-subset(tmp1,tmp1$StudyInstanceUID==tmp2[l,]) #take scans from the same session

      dwi<-subset(tmp3,tmp3$ProtocolName=='ABCD_dMRI'|tmp3$ProtocolName=='ABCD-DTI_GE_original_(baseline_year_1_arm_1)'|tmp3$ProtocolName=='DTI_1'|
                      tmp3$ProtocolName=='DTI_2'|tmp3$ProtocolName=='WIP_DTI_1'|tmp3$ProtocolName=='WIP_DTI_2')

      #apply your criterion to select between repeat scans
      if (dim(dwi)[1]>1){
        seriestime<-data.frame(matrix(NA_character_, nrow = dim(dwi)[1], ncol = 1))
        colnames(seriestime)<-'seriestime'
        seriestime[]<-as.character(1000000)#Filling this empty matrix with a large character
        
        tmp_str<-data.frame(strsplit(as.character(dwi$SeriesInstanceUID), ".",fixed=T)) #SeriesInstanceUID is unique for each scan. 
        #It contains '.' which is why I had to split them this way and treat as character because it couldn't handle numbers this large
        
        for (h in 1:dim(tmp_str)[2]){
          seriestime[h,1]<-paste(tmp_str[,h],collapse = "")#Collapse to obtain a single character for SeriesInstanceUID
        }
        dwi$seriestime<-seriestime$seriestime
          dwi<-dwi[dwi$seriestime==max(dwi$seriestime),] #(I hope SeriesInstanceUID really does reflect the time when it was acquired)
          ######This line above is where the 'last acquired' criterion is applied. This part can be changed if a different criterion is desired.
      }

      #Repeat the same procedure with other image types
      t1<-subset(tmp3,tmp3$ProtocolName=='ABCD_T1w_MPR_vNav'|tmp3$ProtocolName=='ABCD_T1w_MPR_vNav_repeat'|tmp3$ProtocolName=='ABCD-T1_GE_original_(baseline_year_1_arm_1)'|
                   tmp3$ProtocolName=='WIP_3D_T1')

      if (dim(t1)[1]>1){
        seriestime<-data.frame(matrix(NA_character_, nrow = dim(t1)[1], ncol = 1))
        colnames(seriestime)<-'seriestime'
        seriestime[]<-as.character(1000000)
        tmp_str<-data.frame(strsplit(as.character(t1$SeriesInstanceUID), ".",fixed=T))
        for (h in 1:dim(tmp_str)[2]){
          seriestime[h,1]<-paste(tmp_str[,h],collapse = "")
        }
        t1$seriestime<-seriestime$seriestime
        t1<-t1[t1$seriestime==max(t1$seriestime),]
        ######This line above is where the 'last acquired' criterion is applied. This part can be changed if a different criterion is desired.
      }

      b0s<-subset(tmp3,tmp3$ProtocolName=='ABCD_dMRI_DistortionMap_AP'|tmp3$ProtocolName=='ABCD_dMRI_DistortionMap_PA'|
                    tmp3$ProtocolName=='ABCD-Diffusion-FM_GE_original_(baseline_year_1_arm_1)'|tmp3$ProtocolName=='DTI_Fieldmap_A'|
                    tmp3$ProtocolName=='DTI_Fieldmap_P'|tmp3$ProtocolName=='WIP_DTI_Fieldmap_A'|tmp3$ProtocolName=='WIP_DTI_Fieldmap_P')


      rev_b0<-subset(b0s,b0s$PhaseEncodingDirection!=unique(dwi$PhaseEncodingDirection))



      if (dim(rev_b0)[1]>1){
        seriestime<-data.frame(matrix(NA_character_, nrow = dim(rev_b0)[1], ncol = 1))
        colnames(seriestime)<-'seriestime'
        seriestime[]<-as.character(1000000)

        tmp_str<-data.frame(strsplit(as.character(rev_b0$SeriesInstanceUID), ".",fixed=T))
        for (h in 1:dim(tmp_str)[2]){
          
          seriestime[h,1]<-paste(tmp_str[,h],collapse = "")
        }
        rev_b0$seriestime<-seriestime$seriestime

        rev_b0<-rev_b0[rev_b0$seriestime==max(rev_b0$seriestime),]
        ######This line above is where the 'last acquired' criterion is applied. This part can be changed if a different criterion is desired.
      }

      if (dim(dwi)[1]==1 & dim(rev_b0)[1]==1 & dim(t1)[1]==1){
        #If you find they are complete within a session and no previous session had complete data, then keep them. If a previous session had
        #complete data, they will not be kept.

        image_info_BIDS_norev0miss[c(row.names(dwi),row.names(rev_b0),row.names(t1)),c('complete_data_same_session_ASR')]<-1
      }
    }
  }
}


#second loop for data that is only complete between sessions
for (i in 1:dim(subs)[1]){
  tmp1<-subset(image_info_BIDS_norev0miss,subjectkey==subs$subjectkey[i])
  comp_win_sess_dat_check<-subset(tmp1, tmp1$complete_data_same_session_ASR==1) #a variable to verify if the person already has been found to have complete data in a session
  if (dim(comp_win_sess_dat_check)[1]==0){

  tmp2<-data.frame(unique(tmp1$StudyInstanceUID)) #Then take all study UIDs

  #order them such that you will first look at the session with the most data.
  tmp3_sizes<-data.frame(matrix(NA_integer_, nrow = dim(tmp2)[1], ncol = 1))
  colnames(tmp3_sizes)<-'sizes'
  for (q in 1:dim(tmp2)[1]){
    tmp3<-subset(tmp1,tmp1$StudyInstanceUID==tmp2[q,]) #take scans from the same session
    tmp3_sizes[q,1]<-dim(tmp3)[1]
  }
  sorted_idxs<-sort(tmp3_sizes$sizes,decreasing=T,index.return=T)$ix
  sorted_tmp2<-data.frame(tmp2[sorted_idxs,1])#in case of a tie, which session will be taken first is arbitrary

  for (l in 1:dim(sorted_tmp2)[1]){ #this loop tries to prioritize data from the same session

    tmp1<-subset(image_info_BIDS_norev0miss,subjectkey==subs$subjectkey[i])
    #redefine tmp1 within this loop, because tmp1 should change as this loop runs

      #define these scans outside to check for data that was flagged as 2 (data that was not complete in a single session)
      dwi_sub<-subset(tmp1,tmp1$ProtocolName=='ABCD_dMRI'|tmp1$ProtocolName=='ABCD-DTI_GE_original_(baseline_year_1_arm_1)'|tmp1$ProtocolName=='DTI_1'|
                        tmp1$ProtocolName=='DTI_2'|tmp1$ProtocolName=='WIP_DTI_1'|tmp1$ProtocolName=='WIP_DTI_2')
      t1_sub<-subset(tmp1,tmp1$ProtocolName=='ABCD_T1w_MPR_vNav'|tmp1$ProtocolName=='ABCD_T1w_MPR_vNav_repeat'|tmp1$ProtocolName=='ABCD-T1_GE_original_(baseline_year_1_arm_1)'|
                       tmp1$ProtocolName=='WIP_3D_T1')
      b0s_sub<-subset(tmp1,tmp1$ProtocolName=='ABCD_dMRI_DistortionMap_AP'|tmp1$ProtocolName=='ABCD_dMRI_DistortionMap_PA'|
                        tmp1$ProtocolName=='ABCD-Diffusion-FM_GE_original_(baseline_year_1_arm_1)'|tmp1$ProtocolName=='DTI_Fieldmap_A'|
                        tmp1$ProtocolName=='DTI_Fieldmap_P'|tmp1$ProtocolName=='WIP_DTI_Fieldmap_A'|tmp1$ProtocolName=='WIP_DTI_Fieldmap_P')
      rev_b0_sub<-subset(b0s_sub,b0s_sub$PhaseEncodingDirection!=unique(dwi_sub$PhaseEncodingDirection))

      dwi_comp<-subset(dwi_sub, dwi_sub$complete_data_same_session_ASR==2)
      t1_comp<-subset(t1_sub, t1_sub$complete_data_same_session_ASR==2)
      rev_b0_comp<-subset(rev_b0_sub, rev_b0_sub$complete_data_same_session_ASR==2)


      tmp3<-subset(tmp1,tmp1$StudyInstanceUID==sorted_tmp2[l,]) #take scans from the same session

      if (dim(dwi_comp)[1]==0){
        dwi<-subset(tmp3,tmp3$ProtocolName=='ABCD_dMRI'|tmp3$ProtocolName=='ABCD-DTI_GE_original_(baseline_year_1_arm_1)'|tmp3$ProtocolName=='DTI_1'|
                      tmp3$ProtocolName=='DTI_2'|tmp3$ProtocolName=='WIP_DTI_1'|tmp3$ProtocolName=='WIP_DTI_2')

        #apply your criterion to select between repeat scans
        if (dim(dwi)[1]>1){
          seriestime<-data.frame(matrix(NA_character_, nrow = dim(dwi)[1], ncol = 1))
          colnames(seriestime)<-'seriestime'
          seriestime[]<-as.character(1000000)
          tmp_str<-data.frame(strsplit(as.character(dwi$SeriesInstanceUID), ".",fixed=T))
          for (h in 1:dim(tmp_str)[2]){
            seriestime[h,1]<-paste(tmp_str[,h],collapse = "")
          }
          dwi$seriestime<-seriestime$seriestime
          dwi<-dwi[dwi$seriestime==max(dwi$seriestime),] #(I hope SeriesInstanceUID really does reflect the time when it was acquired)
          ######This line above is where the 'last acquired' criterion is applied. This part can be changed if a different criterion is desired.
        }

      }

      if (dim(t1_comp)[1]==0){
        t1<-subset(tmp3,tmp3$ProtocolName=='ABCD_T1w_MPR_vNav'|tmp3$ProtocolName=='ABCD_T1w_MPR_vNav_repeat'|tmp3$ProtocolName=='ABCD-T1_GE_original_(baseline_year_1_arm_1)'|
                     tmp3$ProtocolName=='WIP_3D_T1')

        if (dim(t1)[1]>1){
          seriestime<-data.frame(matrix(NA_character_, nrow = 2, ncol = 1))
          colnames(seriestime)<-'seriestime'
          seriestime[]<-as.character(1000000)
          tmp_str<-data.frame(strsplit(as.character(t1$SeriesInstanceUID), ".",fixed=T))
          for (h in 1:dim(tmp_str)[2]){
            seriestime[h,1]<-paste(tmp_str[,h],collapse = "")
          }
          t1$seriestime<-seriestime$seriestime
          t1<-t1[t1$seriestime==max(t1$seriestime),]
          ######This line above is where the 'last acquired' criterion is applied. This part can be changed if a different criterion is desired.
        }
      }


      if (dim(rev_b0_comp)[1]==0){
        b0s<-subset(tmp3,tmp3$ProtocolName=='ABCD_dMRI_DistortionMap_AP'|tmp3$ProtocolName=='ABCD_dMRI_DistortionMap_PA'|
                      tmp3$ProtocolName=='ABCD-Diffusion-FM_GE_original_(baseline_year_1_arm_1)'|tmp3$ProtocolName=='DTI_Fieldmap_A'|
                      tmp3$ProtocolName=='DTI_Fieldmap_P'|tmp3$ProtocolName=='WIP_DTI_Fieldmap_A'|tmp3$ProtocolName=='WIP_DTI_Fieldmap_P')


        rev_b0<-subset(b0s,b0s$PhaseEncodingDirection!=unique(dwi$PhaseEncodingDirection))



        if (dim(rev_b0)[1]>1){
          seriestime<-data.frame(matrix(NA_character_, nrow = 2, ncol = 1))
          colnames(seriestime)<-'seriestime'
          seriestime[]<-as.character(1000000)

          tmp_str<-data.frame(strsplit(as.character(rev_b0$SeriesInstanceUID), ".",fixed=T))
          for (h in 1:dim(tmp_str)[2]){
            seriestime[h,1]<-paste(tmp_str[,h],collapse = "")
          }
          rev_b0$seriestime<-seriestime$seriestime

          rev_b0<-rev_b0[rev_b0$seriestime==max(rev_b0$seriestime),]
          ######This line above is where the 'last acquired' criterion is applied. This part can be changed if a different criterion is desired.
        }
      }
        #When you try storing a new variable using an empty row.name, it will just ignore that, which is good.
        #These cases are people that have complete data but in separate sessions. So far, this line below can therefore
        #store the data that IS available within a session. The next step is to go look into the next session for the missing data.
        image_info_BIDS_norev0miss[c(row.names(dwi),row.names(rev_b0),row.names(t1)),c('complete_data_same_session_ASR')]<-2
        #stores the data that it does have from the session with the most data. Anything else will be taken from the next largest session.
    }
  }
}

#Isolate all the scans to be kept
final_data_list<-subset(image_info_BIDS_norev0miss,image_info_BIDS_norev0miss$complete_data_same_session_ASR==1|image_info_BIDS_norev0miss$complete_data_same_session_ASR==2)

#Storing it in tmp because I will use it separately
tmp<-subset(image_info_BIDS_norev0miss,image_info_BIDS_norev0miss$complete_data_same_session_ASR==1|image_info_BIDS_norev0miss$complete_data_same_session_ASR==2)

# tmp<-substr(final_data_list$FileName, regexpr('sub', final_data_list$FileName), regexpr('.json', final_data_list$FileName)-1)
final_data_list$FileNameFinal<-switch_names(final_data_list$FileName,'.json','.nii.gz','','R')

final_sub_list<-data.frame(unique(final_data_list$subjectkey))


group_n<-dim(final_sub_list)[1]
cat('There are',group_n,'subjects with complete data','\n')


#Loop to add the image type back to final list
for (i in 1:dim(final_data_list)[1]){


dwi<-subset(final_data_list[i,],final_data_list[i,]$ProtocolName=='ABCD_dMRI'|final_data_list[i,]$ProtocolName=='ABCD-DTI_GE_original_(baseline_year_1_arm_1)'|final_data_list[i,]$ProtocolName=='DTI_1'|
                  final_data_list[i,]$ProtocolName=='DTI_2'|final_data_list[i,]$ProtocolName=='WIP_DTI_1'|final_data_list[i,]$ProtocolName=='WIP_DTI_2')
t1<-subset(final_data_list[i,],final_data_list[i,]$ProtocolName=='ABCD_T1w_MPR_vNav'|final_data_list[i,]$ProtocolName=='ABCD_T1w_MPR_vNav_repeat'|final_data_list[i,]$ProtocolName=='ABCD-T1_GE_original_(baseline_year_1_arm_1)'|
                 final_data_list[i,]$ProtocolName=='WIP_3D_T1')
b0s<-subset(final_data_list[i,],final_data_list[i,]$ProtocolName=='ABCD_dMRI_DistortionMap_AP'|final_data_list[i,]$ProtocolName=='ABCD_dMRI_DistortionMap_PA'|
                  final_data_list[i,]$ProtocolName=='ABCD-Diffusion-FM_GE_original_(baseline_year_1_arm_1)'|final_data_list[i,]$ProtocolName=='DTI_Fieldmap_A'|
                  final_data_list[i,]$ProtocolName=='DTI_Fieldmap_P'|final_data_list[i,]$ProtocolName=='WIP_DTI_Fieldmap_A'|final_data_list[i,]$ProtocolName=='WIP_DTI_Fieldmap_P')

  if (dim(dwi)[1]!=0){
    final_data_list[i,c("ProtocolFinal")]<-"dwi"
  } else if (dim(t1)[1]!=0){
    final_data_list[i,c("ProtocolFinal")]<-"anat"
  } else if (dim(b0s)[1]!=0){
    final_data_list[i,c("ProtocolFinal")]<-"fmap"
  }
}



final_data_list_clean <- final_data_list[,c('subjectkey','FileName','FileNameFinal','ProtocolFinal')]
# tmp<-substr(final_data_list_clean$subjectkey, regexpr('_', final_data_list_clean$subjectkey)+1, regexpr('sub_', final_data_list_clean$subjectkey)+30)
# final_data_list_clean$subjectkey<-paste("sub-NDAR", tmp, sep="")
###########################################################
###########################################################

##########This part handles group-specific issues#####################
if (grp==1){
  #mTBI Group
  group_dir<-'/mTBI_Group'
  final_data_list_clean_clean<-final_data_list_clean
} else if (grp==2){
  #Harmonization Group
  #This loop handles redundancies left after handling data inconsistecies (reference group lost subs + we had matched using a larger ratio to have a buffer)
  #It is almost exactly like the Select_HC_Participants script except that the match pool is now the subjects left after checking for data inconsistencies, and the ratio is 1:1
  group_dir<-'/Harmonization_Controls'
  harmonization_list<-read.table(paste(args[1],'/Harmonization_Data.txt',sep=''), sep='', header=T)
  #NOTE (17/01/2020): I don't remember why the Harmonization_Data.txt ends up with a col of sub ids at the end
  # harmonization_list<-read.table('/Users/Guido/Desktop/GeneralData/ABCD/Harmonization_Controls/Harmonization_Data.txt', sep='', header=F)#For troubleshooting
  
  
  colnames(harmonization_list)<-c('subjectkey','mri_info_deviceserialnumber','pubertal_stage','sex','ehi_y_ss_scoreb','IQ_crystalized','group','distance','weights')
    
  harmonization_list$subjectkey<-switch_names(harmonization_list$subjectkey,Rem = 'NDAR_',Add = 'sub-NDAR',separ = '',side = 'L')
  
  complete_data<-merge(harmonization_list,final_data_list_clean,by = 'subjectkey')

  sites<-unique(data.frame(complete_data$mri_info_deviceserialnumber))
  colnames(sites)<-'mri_info_deviceserialnumber'
  
    ref_site<-"HASHe76e6d72"
    # ref_site<-'HASHc9398971'#this target site was selected because it has the lowest average z score (for FA, MD, AD, RD, AFD_Max, AFD_total, NuFO) WILL HAVE TO IMPLEMENT THIS PART LATER.
    ref_site_data<-subset(complete_data,complete_data$mri_info_deviceserialnumber==ref_site)
    # 'HASHe76e6d72' Reference site
    desired_dim<-dim(unique(data.frame(ref_site_data$subjectkey)))[1]#Defined based on how many people are left in the reference group
    cat('There are',desired_dim,'subjects with complete data in reference group',ref_site,'\n')
    
    
    sites$complete<-1 #Prefill this column
    for (k in 1:dim(sites)[1]){
      tmp<-subset(complete_data,complete_data$mri_info_deviceserialnumber==sites[k,]$mri_info_deviceserialnumber)
      if (dim(unique(data.frame(tmp$subjectkey)))[1]<desired_dim){
        # If there are less subjects than needed
        sites[k,]$complete<-0
      } else if (dim(unique(data.frame(tmp$subjectkey)))[1]>desired_dim){
        # If there are more subjects than needed
        sites[k,]$complete<-2
        # Everybody else will be left with a 1=exactly the number of subjects needed
      }
    }

    # TD <- complete_data[c('subjectkey', 'mri_info_deviceserialnumber', 'interview_age', 'sex', 'ehi_y_ss_scoreb', 'IQ_crystalized')]
    # TD <- complete_data[c('subjectkey', 'mri_info_deviceserialnumber', 'pubertal_stage', 'sex', 'ehi_y_ss_scoreb', 'IQ_crystalized')]
    TD <- complete_data[c('subjectkey', 'mri_info_deviceserialnumber', 'pubertal_stage', 'sex', 'ehi_y_ss_scoreb', 'IQ_crystalized','group')]
    # TD <- complete_data[,c('subjectkey','FileName','FileNameFinal','ProtocolFinal')]#12/02/2020: TRYING NOW TO PROCESS THE LARGE POOL OF SUBS, QC THEM, AND THEN EITHER MATCH AGAIN ON THE QCED SUBS, OR PROCEED WITH HARMONIZATION WITH THE ONES THAT PASS QC (DEPENDING ON HOW MANY LEFT)
    TD<-unique(TD)# TD had 3 rows per subject, I only need one
    #TD <- TD_full[c('subjectkey', 'mri_info_deviceserialnumber', 'interview_age.x', 'sex', 'ehi_y_ss_scoreb', 'IQ_crystalized')]
    
    #Get site sizes
    sites$sizes<-1 #Prefill this column
    for(q in 1:dim(sites)[1]){ #Loop through the numbers of ID's instead of the ID's
      site<-sites[q,]$mri_info_deviceserialnumber
      site_dat<-subset(TD,TD$mri_info_deviceserialnumber==site)
      sites[q,]$sizes<-dim(site_dat)[1]
    }
    
    ############################################################
    ##                                                        ##
    ##               T Y M    2    M A T C H                  ##
    ##                                                        ##
    ##                                                        ##
    ##                        MATCH ON:                       ##
    ##              age, gender, iq, handedness               ##
    ##                                                        ##
    ##           HASHe76e6d72 - is the reference site         ##
    ##                ^^ bc is the smollest                   ##
    ##                                                        ##
    ############################################################
    # table(TD$mri_info_deviceserialnumber=='HASHe76e6d72')
    
    ## to verify that the hashed serial number is the number of scanners i'm going to...
    ## ## use the publication 'Detecting and harmonizing scanner differences in the ABCD study - annual release 1.0' - Nielson, 2018 that said they used hashed serial numbers
    
    ## make a dumby variable for group - 0: reference, 1: target
    target<-data.frame(subset(TD,TD$group==1))
    
    ## force age to be a number!
    # TD$interview_age <- as.numeric(as.character(TD$interview_age))
    
    all_sites<-unique(data.frame(TD$mri_info_deviceserialnumber))
    # dim(all_sites)
    #Note: there will be 30, because some people don't have a serial number
    
    colnames(all_sites)<-'sites'
    
    #remove target from the list of sites
    match_sites<-subset(all_sites,all_sites$sites!='HASHe76e6d72')
    # match_sites<-subset(all_sites,all_sites$sites!='HASHc9398971')
    
    
    
    # dim(match_sites)
    #Note: 28 without counting target site. Matching occurs in pairs with the reference, there are 29 unique serial numbers therefore 28 groups to match
    
    ## Pre-define an empty matrix, whose dimensions will be defined by the amount of sites, data in target, and match ratio
    #Amount of people per match
    step_size=(as.numeric(dim(target)[1])+as.numeric(as.numeric(args[3])*as.numeric(dim(target)[1])))
    # step_size=(as.numeric(dim(target)[1])+as.numeric(as.numeric(2)*as.numeric(dim(target)[1])))#For troubleshoot
    #NOTE: Number of people in target site + number of people in match site (target * match ratio)
    
    # step_size=(as.numeric(dim(target)[1])+as.numeric(2*as.numeric(dim(target)[1])))#For troubleshooting
    
    #Create an empty data frame the size of #pairs (i.e.: nsites-1) * step_size
    full_size = as.numeric(dim(match_sites)[1])*as.numeric(step_size)
    
    # To define number of columns
    variables = as.numeric(dim(TD)[2])+3
    #NOTE: the 7 original variables in TD + the two that matchit always adds (distance, weights)+match_iteration
    
    output <- data.frame(matrix(NA,ncol=variables, nrow=full_size))
    # colnames(output)<-c('subjectkey','mri_info_deviceserialnumber','interview_age','sex','ehi_y_ss_scoreb','IQ_crystalized','group','distance','weights')
    colnames(output)<-c('subjectkey','mri_info_deviceserialnumber','pubertal_stage','sex','ehi_y_ss_scoreb','IQ_crystalized','group','distance','weights','match_idx')
    
    TD<-na.omit(TD)#REMEMBER: TD IS NOW THE CLEANED POOL OF SUBS WITH COMPLETE DATA
    cnt=0#start a counter to index which matching iteration you're doing (because the target site will be duplicated a bunch of times)
    for(i in 1:dim(match_sites)[1]){
      cnt=cnt+1
      #Pool of data from which to match (contains target site data + one other site's data)
      match_pool <- subset(TD, mri_info_deviceserialnumber=='HASHe76e6d72' | mri_info_deviceserialnumber==match_sites[i,])
      # match_pool <- subset(TD, mri_info_deviceserialnumber=='HASHc9398971' | mri_info_deviceserialnumber==match_sites[i,])
      # match.it <- matchit(group~ interview_age + sex + ehi_y_ss_scoreb + IQ_crystalized , data=match_pool, method='nearest', ratio=as.numeric(args[3]))
      # # match.it <- matchit(group~ interview_age + sex + ehi_y_ss_scoreb + IQ_crystalized , data=match_pool, method='nearest', ratio=2) #For troubleshooting
      
      match.it <- matchit(group~ pubertal_stage + sex + ehi_y_ss_scoreb + IQ_crystalized , data=match_pool, method='nearest', ratio=as.numeric(args[3]))
      # match.it <- matchit(group~ pubertal_stage + sex + ehi_y_ss_scoreb + IQ_crystalized , data=match_pool, method='nearest', ratio=1) #For troubleshooting
      
      data<-match.data(match.it)
      data$match_idx<-cnt
      ## Variables to create a moving window in which to store the data
      new_ss<-dim(data)[1]
      start_idx<-(1+(as.numeric(step_size)*i)-as.numeric(step_size))#Starts at beginning of each window. If the previous one was shorter, it will skip a few lines
      end_idx<-as.numeric(start_idx)+(as.numeric(new_ss))-1#Ends at start+size of matched data. Ensures that the index always works even if there's not enough data to respect the ratio
      
      # dim(data)#For troubleshooting
      # dim(data.frame(start_idx:end_idx))#For troubleshooting
      
      #Had to make sure these were treated as characters, they were being converted into numbers
      data$subjectkey<-as.character(data$subjectkey)
      data$mri_info_deviceserialnumber<-as.character(data$mri_info_deviceserialnumber)
      
      #Store data
      output[start_idx:end_idx,] <- data
      
    }
    
    
    
    output<-subset(output,is.na(output$subjectkey)==F)

    #Remove duplicated data (the target site will be duplicated a bunch of times)
    DATAFORHARMONIZATION <- output[!duplicated(output$subjectkey), ]
    DATAFORHARMONIZATION<-subset(DATAFORHARMONIZATION,is.na(DATAFORHARMONIZATION$subjectkey)==F)
    
    
  
    
    #If site has >20, randomly select 20 subjects from them. Otherwise, keep all the subjects from a site
    listofdfs <- list() #Create a list in which you intend to save your df's.
    for(l in 1:dim(sites)[1]){
      site<-sites[l,]$mri_info_deviceserialnumber
      site_dat<-subset(TD,TD$mri_info_deviceserialnumber==site)
      if (dim(site_dat)[1]>20){
        df<-site_dat[sample(nrow(site_dat),20), ]#RANDOM SAMPLING DONE HERE
      }else{
        df<-site_dat
      }
      listofdfs[[l]] <- df # save your dataframes into the list
    }
    
    TD_new<-bind_rows(listofdfs,.id="column_label")
    dim(TD_new)
    
    #Do a first check to see if sites are matched
    # Compute the analysis of variance
    res.aov <- aov(sex ~ mri_info_deviceserialnumber, data = TD_new)
    # summary(res.aov)
    sex<-summary(res.aov)[[1]][["Pr(>F)"]][1]
    
    res.aov <- aov(ehi_y_ss_scoreb ~ mri_info_deviceserialnumber, data = TD_new)
    # summary(res.aov)
    hand<-summary(res.aov)[[1]][["Pr(>F)"]][1]
    
    res.aov <- aov(IQ_crystalized ~ mri_info_deviceserialnumber, data = TD_new)
    # summary(res.aov)
    IQ<-summary(res.aov)[[1]][["Pr(>F)"]][1]
    
    # res.aov <- aov(interview_age ~ mri_info_deviceserialnumber, data = TD_new)
    res.aov <- aov(pubertal_stage ~ mri_info_deviceserialnumber, data = TD_new)
    # summary(res.aov)
    age<-summary(res.aov)[[1]][["Pr(>F)"]][1]
    
    
    

    if (sex>0.05&hand>0.05&IQ>0.05&age>0.05){
      #If sites are matched, tell the person and just proceed
      cat('Sites remain matched after removing data inconsistencies and random selection','\n')
    } else {
      #If they are not, tell the person and then try randomly selecting a pre-defined number of tries. If that still doesn't work, tell the person to go back and change the initial matching ratio to get a bigger buffer.
      cat('Sites are not matched, redoing random selection','\n')
      cnt=0
      chk=0
      # while(cnt<=args[3]&chk==0){ #Try X times, but stop if you've found a matched result
      while(cnt<=3&chk==0){ #for troubleshoot
        
        #If site has >20, randomly select 20 subjects from them. Otherwise, keep all the subjects from a site
        listofdfs <- list() #Create a list in which you intend to save your df's.
        for(l in 1:dim(sites)[1]){
          site<-sites[l,]$mri_info_deviceserialnumber
          site_dat<-subset(TD,TD$mri_info_deviceserialnumber==site)
          if (dim(site_dat)[1]>20){
            df<-site_dat[sample(nrow(site_dat),20), ]
          }else{
            df<-site_dat
          }
          listofdfs[[l]] <- df # save your dataframes into the list
        }
        
        TD_new<-bind_rows(listofdfs,.id="column_label")
        
        res.aov <- aov(sex ~ mri_info_deviceserialnumber, data = TD_new)
        sex<-summary(res.aov)[[1]][["Pr(>F)"]][1]
        
        res.aov <- aov(ehi_y_ss_scoreb ~ mri_info_deviceserialnumber, data = TD_new)
        hand<-summary(res.aov)[[1]][["Pr(>F)"]][1]
        
        res.aov <- aov(IQ_crystalized ~ mri_info_deviceserialnumber, data = TD_new)
        IQ<-summary(res.aov)[[1]][["Pr(>F)"]][1]
        
        res.aov <- aov(interview_age ~ mri_info_deviceserialnumber, data = TD_new)
        age<-summary(res.aov)[[1]][["Pr(>F)"]][1]
        
        if (sex>0.05&hand>0.05&IQ>0.05&age>0.05){
          chk=1
        }
        cnt=cnt+1#Counter goes up on every try
      }
    }
    
   
    
    DATAFORHARMONIZATION<-TD_new
    group_n<-dim(unique(data.frame(DATAFORHARMONIZATION$subjectkey)))[1]
    cat('There are',group_n,'subjects for harmonization after second round of matching','\n')

    
    complete_data_matched<-complete_data[complete_data$subjectkey %in% DATAFORHARMONIZATION$subjectkey, ]
    # dim(unique(data.frame(complete_data_matched$subjectkey))) # 348, 12 subjects per site (29 sites)
    # dim(complete_data_matched)[1]#1044, 3 files per person
    #
    complete_data_matched_clean <- complete_data_matched[,c('subjectkey','FileName','FileNameFinal','ProtocolFinal')]
    final_data_list_clean_clean<-complete_data_matched_clean
} else if (grp==3){
  #Matched Control Group
} else if (grp==4){
  #Twin group
}











# dim(complete_data_matched_clean)
write.table(final_data_list_clean_clean, file = paste(args[1],'/Data/BIDS/final_data_list.txt',sep=''), sep = ' ',eol = '\n',
            row.names = FALSE, col.names = FALSE, quote=FALSE)
# write.table(TD, file = '/Users/Guido/Desktop/GeneralData/ABCD/Harmonization_Controls/Data/BIDS/final_data_list.txt', sep = ' ',eol = '\n',
#             row.names = FALSE, col.names = FALSE, quote=FALSE)



############EXTRA#############
# write.table(final_data_list_clean, file = "/Users/guigub/Documents/Research/PhD/GeneralData/ABCD/MC_Group/final_data_list.txt", sep = ' ',eol = '\n',
#             row.names = FALSE, col.names = FALSE, quote=FALSE)

# #writing a table with subject ID and manufacturer for comparing FRFs between vendors (in another script)
# tmp<-subset(final_data_list,is.na(final_data_list$PhaseEncodingDirection)==T)
# IDs_manufacturers<-tmp[,c('subjectkey','Manufacturer')]

# write.table(IDs_manufacturers, file = "/Users/guigub/Documents/Research/PhD/GeneralData/ABCD/mTBIGroup/ID_manufacturers.txt", sep = ' ',eol = '\n',
            # row.names = FALSE, col.names = FALSE, quote=FALSE)