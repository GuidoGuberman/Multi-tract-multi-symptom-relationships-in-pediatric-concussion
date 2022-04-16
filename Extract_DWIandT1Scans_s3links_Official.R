##########################################################
##        SCRIPT TO DO FIRST CHECK FOR MISSING DATA     ##
##########################################################

args = commandArgs(trailingOnly=TRUE)
library(dplyr)
image_info <- read.table(paste(args[1],'/GeneralData/image03.txt',sep =''), sep='\t', header=T)
# image_info <- read.table('/Users/Guido/Desktop/GeneralData/ABCD/GeneralData/image03.txt', sep='\t', header=T)#for troubleshooting
datadic_image_info <- image_info[1,]
datadic_image_info <-t(datadic_image_info)
image_info <-image_info[-1,]


###########GROUP (Defined in Master Script)################
grp=args[2]

if (grp==1){
  #mTBI Group
  group_dir<-'/mTBI_Group'
} else if (grp==2){
  #Harmonization Group
  group_dir<-'/Harmonization_Controls'
} else if (grp==3){
  #Matched Control Group
} else if (grp==4){
  #Twin group
}

##########################################################


# Load IDs from the selection step
subs<-read.table(paste(args[1],group_dir,'/IDs.txt',sep =''), sep='\t', header=T)
# subs<-read.table('/Users/Guido/Desktop/GeneralData/ABCD/mTBI_Group/IDs.txt', sep='\t', header=T)#for troubleshooting
#Line with specific path, for troubleshooting 
# subs<-read.table('/Users/Guido/Desktop/GeneralData/ABCD/Harmonization_Controls/IDs.txt', sep='\t', header=T)
  
colnames(subs)<-'subjectkey' #Name the column
  
image_info_subs<-image_info[image_info$subjectkey %in% subs$subjectkey, ]#Get, from image_info, subs from the group of interest

group_n<-dim(unique(data.frame(image_info_subs$subjectkey)))[1]
cat('There are',group_n,'subjects with image info','\n')
  

########EXTRACT DWI AND T1 SCANS #########
#First select baseline scans only
image_info_subs_baseline <- subset(image_info_subs, visit == 'baseline_year_1_arm_1')
group_n<- dim(unique(data.frame(image_info_subs_baseline$subjectkey)))[1]
cat('There are',group_n,'subjects with baseline data','\n')

#Then look for necessary images
DWI_and_T1_scans_baseline <- subset(image_info_subs_baseline, image_description== 'ABCD-Diffusion-FM-PA' | image_description== 'ABCD-Diffusion-FM-AP' | 
                                      image_description== 'ABCD-DTI' | image_description== 'ABCD-Diffusion-FM' | image_description == 'ABCD-T1')
  
group_n<- dim(unique(data.frame(DWI_and_T1_scans_baseline$subjectkey)))[1]
cat('There are',group_n,'subjects with minimal image requirements','\n')
  
# #Clean up this table
DWI_and_T1_scans_baseline_clean<-DWI_and_T1_scans_baseline[c('subjectkey','interview_date','interview_age','sex','image03_id','image_file',
                                                             'image_description','scanner_manufacturer_pd','visit')]

#Extract series time from large image_info table (will be used later)
DWI_and_T1_scans_baseline_clean$iqc_seriestime=substr(DWI_and_T1_scans_baseline_clean$image_file,
                                                      regexpr('.tgz', DWI_and_T1_scans_baseline_clean$image_file)-6,
                                                      regexpr('.tgz', DWI_and_T1_scans_baseline_clean$image_file)-1)
  
  
###EXTRACT SCANNER MANUFACTURER###
scanner <- DWI_and_T1_scans_baseline_clean[,c('subjectkey','scanner_manufacturer_pd','image_description')]

#####HANDLE DATA DUPLICATIONS####
#this part flags them
  for (i in 1:dim(data.frame(DWI_and_T1_scans_baseline_clean$subjectkey))[1]){
    #Find all scans for a subject
    tmp1<-subset(DWI_and_T1_scans_baseline_clean,subjectkey==DWI_and_T1_scans_baseline_clean$subjectkey[i])
    #Find all scans with the same image_description
    scans<-subset(tmp1,image_description==DWI_and_T1_scans_baseline_clean[i,]$image_description)
    #Find how many images of a particular type this subject has
    #Note: all images have a unique image_id, even if duplicated
    tmp2<-dim(data.frame(unique(scans$image03_id)))[1]
    #Find how many images with unique series time there are (duplicated images will have the same series time)
    tmp3<-dim(data.frame(unique(scans$iqc_seriestime)))[1]
    #If there are less unique series time than image_ids, it means some images were duplicated
    DWI_and_T1_scans_baseline_clean$duplicated[i] <- ifelse(tmp2!=tmp3, 1,0)
  }
  
#this part ensures that only one of the duplications is kept
  for (i in 1:dim(data.frame(DWI_and_T1_scans_baseline_clean$subjectkey))[1]){
    #Find all scans for a subject
    tmp1<-subset(DWI_and_T1_scans_baseline_clean,subjectkey==DWI_and_T1_scans_baseline_clean$subjectkey[i])
    
    #Find all scans with the same image_description
    scans_old<-subset(tmp1,image_description==DWI_and_T1_scans_baseline_clean[i,]$image_description)
    
    #Find all scans with the same series time
    scans<-subset(scans_old,iqc_seriestime==DWI_and_T1_scans_baseline_clean[i,]$iqc_seriestime)
    
    #The number of duplicated images = number of images with the same series time
    if (sum(scans$duplicated)==dim(scans)[1]){
      #The first time it will encounter duplicated images, this condition will be met. Next time, it will not. So only one of the 
      #duplicated images will be labelled with 0 to be kept
      DWI_and_T1_scans_baseline_clean$duplicated[i] = 0
    }
  }
  
  DWI_and_T1_scans_baseline_clean_noDups <-subset(DWI_and_T1_scans_baseline_clean,duplicated==0)
  
  group_n<- dim(data.frame(unique(DWI_and_T1_scans_baseline_clean_noDups$subjectkey)))[1]
  cat('There are',group_n,'subjects after removing duplications','\n')
  
#####FLAG DATA INCONSISTENCIES####
#Repeats: same image type but not the same exact image; missing: missing either T1, DWI, or PA/FM)

############################################################
##                                                        ##
##               R E Q U I R E M E N T S                  ##
##        Philips: 2 DWI, 1 T1, 1 FM-AP and FM-PA         ##
##        GE: 1 DWI, 1 T1, 1 FM                           ##
##        Siemens: 1 DWI, 1 T1, 1 FM-AP and FM-PA         ##
##                                                        ##
##                                                        ##
##                                                        ##
##                                                        ##
##  Note: AP and PA are both required, since labels are   ##
##       not to be trusted                                ##
############################################################

  temp <- DWI_and_T1_scans_baseline_clean_noDups %>% group_by(subjectkey) %>% summarise(DTI_FM = sum(image_description == 'ABCD-Diffusion-FM'), DTI_FM_PA = sum(image_description == 'ABCD-Diffusion-FM-PA'),
                                                                                          DTI_FM_AP = sum(image_description == 'ABCD-Diffusion-FM-AP'), ABCD_DTI = sum(image_description == 'ABCD-DTI'), ABCD_T1 = sum(image_description =='ABCD-T1' ))
  FINAL <- merge(temp, DWI_and_T1_scans_baseline_clean_noDups, by='subjectkey', all=T)
    
  PHILIPS <- subset(FINAL, scanner_manufacturer_pd=='Philips Medical Systems')
  PHILIPS$missingdata <- ifelse(PHILIPS$ABCD_DTI<2|PHILIPS$DTI_FM_PA<1|PHILIPS$DTI_FM_AP<1|PHILIPS$ABCD_T1<1, 1,0)
  PHILIPS$T1REPEAT <- ifelse(PHILIPS$ABCD_T1>1, 1, 0)
  PHILIPS$DWIREPEAT <- ifelse(PHILIPS$ABCD_DTI>2 | PHILIPS$DTI_FM_PA>1 | PHILIPS$DTI_FM_AP>1, 1,0)
    
  SIEMENS <- subset(FINAL, scanner_manufacturer_pd=='SIEMENS')
  SIEMENS$missingdata <- ifelse(SIEMENS$ABCD_DTI<1|SIEMENS$DTI_FM_PA<1|SIEMENS$DTI_FM_AP<1|SIEMENS$ABCD_T1<1, 1,0)
  SIEMENS$T1REPEAT <- ifelse(SIEMENS$ABCD_T1>1, 1, 0)
  SIEMENS$DWIREPEAT <- ifelse(SIEMENS$ABCD_DTI>1 | SIEMENS$DTI_FM_PA>1 | SIEMENS$DTI_FM_AP>1, 1,0)
    
  GE <- subset(FINAL, scanner_manufacturer_pd=='GE MEDICAL SYSTEMS')
  GE$missingdata <- ifelse(GE$ABCD_DTI<1|GE$DTI_FM<1|GE$ABCD_T1<1, 1,0)
  GE$T1REPEAT <- ifelse(GE$ABCD_T1>1, 1, 0)
  GE$DWIREPEAT <- ifelse(GE$ABCD_DTI>1 | GE$DTI_FM>1, 1,0)
  FINAL_inconsistencies <- rbind(GE,SIEMENS,PHILIPS) #Here's the important table for the next steps.
    
  
#Checking number of subjects lost so far
  subs_missingdata<-data.frame(unique(FINAL_inconsistencies$subjectkey[FINAL_inconsistencies$missingdata==1]))
  colnames(subs_missingdata)<-'subjectkey'
  group_nn<-dim(subs_missingdata)[1]
  
  FINAL_inconsistencies<-subset(FINAL_inconsistencies,missingdata!=1)
  subs_FINAL_inconsistencies<-data.frame(unique(FINAL_inconsistencies$subjectkey))
  colnames(subs_FINAL_inconsistencies)<-'subjectkey' #Name the column
  
  group_n<-dim(subs_FINAL_inconsistencies)[1]#should be 350 (mtbi) or 460 (Match Controls)
  cat('There are',group_n,'subjects with complete data (',group_nn,'missing) after removing subjects without minimum image requirements','\n')
  

  
######QC CHECK#########
  QC_file <- read.table(paste(args[1],'/GeneralData/mriqcrp102.txt',sep =''), sep='\t', header=T)
  
  #File with path, for troubleshooting
  # QC_file <- read.table('/Users/Guido/Desktop/GeneralData/ABCD/GeneralData/mriqcrp102.txt', sep='\t', header=T)
  
  datadic_QC_file <- QC_file[1,]
  datadic_QC_file <-t(datadic_QC_file)
  QC_file <-QC_file[-1,]
  
  QC_file_subs<-QC_file[QC_file$subjectkey %in% subs_FINAL_inconsistencies$subjectkey, ]#Get from the selected subs so far, the subs from the group of interest
  
  group_n<-dim(data.frame(unique(QC_file_subs$subjectkey)))[1]
  cat('There are',group_n,'subjects in QC file','\n')
  
  
  T1_QC <- data.frame(QC_file_subs[,c('subjectkey','iqc_t1_1_qc_score','iqc_t1_1_seriestime','iqc_t1_2_qc_score','iqc_t1_2_seriestime',
                                             'iqc_t1_3_qc_score','iqc_t1_3_seriestime')])
  
  
  DWI_QC <- data.frame(QC_file_subs[,c('subjectkey','iqc_dmri_1_qc_score','iqc_dmri_1_seriestime','iqc_dmri_1_fm_missing',
                                              'iqc_dmri_2_qc_score','iqc_dmri_2_seriestime','iqc_dmri_2_fm_missing',
                                              'iqc_dmri_3_qc_score','iqc_dmri_3_seriestime','iqc_dmri_3_fm_missing',
                                              'iqc_dmri_4_qc_score','iqc_dmri_4_seriestime','iqc_dmri_4_fm_missing',
                                              'iqc_dmri_5_qc_score','iqc_dmri_5_seriestime','iqc_dmri_5_fm_missing','iqc_dmri_6_qc_score','iqc_dmri_6_seriestime','iqc_dmri_6_fm_missing')])
  
  
  #Integrate QC info to image dataset
  #Renaming columns in preparation for switching from large to long format
  names(DWI_QC)[names(DWI_QC) == "iqc_dmri_1_qc_score"] <- "iqc_dmri_qc_scores.1"
  names(DWI_QC)[names(DWI_QC) == "iqc_dmri_2_qc_score"] <- "iqc_dmri_qc_scores.2"
  names(DWI_QC)[names(DWI_QC) == "iqc_dmri_3_qc_score"] <- "iqc_dmri_qc_scores.3"
  names(DWI_QC)[names(DWI_QC) == "iqc_dmri_4_qc_score"] <- "iqc_dmri_qc_scores.4"
  names(DWI_QC)[names(DWI_QC) == "iqc_dmri_5_qc_score"] <- "iqc_dmri_qc_scores.5"
  names(DWI_QC)[names(DWI_QC) == "iqc_dmri_6_qc_score"] <- "iqc_dmri_qc_scores.6"
  names(DWI_QC)[names(DWI_QC) == "iqc_dmri_1_seriestime"] <- "iqc_dmri_seriestime.1"
  names(DWI_QC)[names(DWI_QC) == "iqc_dmri_2_seriestime"] <- "iqc_dmri_seriestime.2"
  names(DWI_QC)[names(DWI_QC) == "iqc_dmri_3_seriestime"] <- "iqc_dmri_seriestime.3"
  names(DWI_QC)[names(DWI_QC) == "iqc_dmri_4_seriestime"] <- "iqc_dmri_seriestime.4"
  names(DWI_QC)[names(DWI_QC) == "iqc_dmri_5_seriestime"] <- "iqc_dmri_seriestime.5"
  names(DWI_QC)[names(DWI_QC) == "iqc_dmri_6_seriestime"] <- "iqc_dmri_seriestime.6"
  names(DWI_QC)[names(DWI_QC) == "iqc_dmri_1_fm_missing"] <- "iqc_dmri_fm_missing.1"
  names(DWI_QC)[names(DWI_QC) == "iqc_dmri_2_fm_missing"] <- "iqc_dmri_fm_missing.2"
  names(DWI_QC)[names(DWI_QC) == "iqc_dmri_3_fm_missing"] <- "iqc_dmri_fm_missing.3"
  names(DWI_QC)[names(DWI_QC) == "iqc_dmri_4_fm_missing"] <- "iqc_dmri_fm_missing.4"
  names(DWI_QC)[names(DWI_QC) == "iqc_dmri_5_fm_missing"] <- "iqc_dmri_fm_missing.5"
  names(DWI_QC)[names(DWI_QC) == "iqc_dmri_6_fm_missing"] <- "iqc_dmri_fm_missing.6"
  
  names(T1_QC)[names(T1_QC) == "iqc_t1_1_qc_score"] <- "iqc_t1_qc_scores.1"
  names(T1_QC)[names(T1_QC) == "iqc_t1_2_qc_score"] <- "iqc_t1_qc_scores.2"
  names(T1_QC)[names(T1_QC) == "iqc_t1_3_qc_score"] <- "iqc_t1_qc_scores.3"
  names(T1_QC)[names(T1_QC) == "iqc_t1_1_seriestime"] <- "iqc_t1_seriestime.1"
  names(T1_QC)[names(T1_QC) == "iqc_t1_2_seriestime"] <- "iqc_t1_seriestime.2"
  names(T1_QC)[names(T1_QC) == "iqc_t1_3_seriestime"] <- "iqc_t1_seriestime.3"
  
  #Ensuring all numerical values will be treated as such
  T1_QC$iqc_t1_qc_scores.1<-as.numeric(as.vector(T1_QC$iqc_t1_qc_scores.1))
  T1_QC$iqc_t1_qc_scores.2<-as.numeric(as.vector(T1_QC$iqc_t1_qc_scores.2))
  T1_QC$iqc_t1_qc_scores.3<-as.numeric(as.vector(T1_QC$iqc_t1_qc_scores.3))
  
  T1_QC$iqc_t1_seriestime.1<-as.numeric(as.vector(T1_QC$iqc_t1_seriestime.1))
  T1_QC$iqc_t1_seriestime.2<-as.numeric(as.vector(T1_QC$iqc_t1_seriestime.2))
  T1_QC$iqc_t1_seriestime.3<-as.numeric(as.vector(T1_QC$iqc_t1_seriestime.3))
  
  
  DWI_QC$iqc_dmri_qc_scores.1<-as.numeric(as.vector(DWI_QC$iqc_dmri_qc_scores.1))
  DWI_QC$iqc_dmri_qc_scores.2<-as.numeric(as.vector(DWI_QC$iqc_dmri_qc_scores.2))
  DWI_QC$iqc_dmri_qc_scores.3<-as.numeric(as.vector(DWI_QC$iqc_dmri_qc_scores.3))
  DWI_QC$iqc_dmri_qc_scores.4<-as.numeric(as.vector(DWI_QC$iqc_dmri_qc_scores.4))
  DWI_QC$iqc_dmri_qc_scores.5<-as.numeric(as.vector(DWI_QC$iqc_dmri_qc_scores.5))
  DWI_QC$iqc_dmri_qc_scores.6<-as.numeric(as.vector(DWI_QC$iqc_dmri_qc_scores.6))
  
  DWI_QC$iqc_dmri_seriestime.1<-as.numeric(as.vector(DWI_QC$iqc_dmri_seriestime.1))
  DWI_QC$iqc_dmri_seriestime.2<-as.numeric(as.vector(DWI_QC$iqc_dmri_seriestime.2))
  DWI_QC$iqc_dmri_seriestime.3<-as.numeric(as.vector(DWI_QC$iqc_dmri_seriestime.3))
  DWI_QC$iqc_dmri_seriestime.4<-as.numeric(as.vector(DWI_QC$iqc_dmri_seriestime.4))
  DWI_QC$iqc_dmri_seriestime.5<-as.numeric(as.vector(DWI_QC$iqc_dmri_seriestime.5))
  DWI_QC$iqc_dmri_seriestime.6<-as.numeric(as.vector(DWI_QC$iqc_dmri_seriestime.6))
  
  
  #NOTE: the way this script works starting here is that images (in image_info) with missing QC data, or QC fails will be removed.
  # Since some subjects might have redundancies, the number of scans per person will be tallied up again at the end and that is 
  # when subjects will once again be checked for missing data, and removed if necessary
  
  #Reshape T1 QC file
  T1_QC_long<-reshape(T1_QC, varying = list(c("iqc_t1_qc_scores.1", "iqc_t1_qc_scores.2","iqc_t1_qc_scores.3"),
                                                        c("iqc_t1_seriestime.1", "iqc_t1_seriestime.2","iqc_t1_seriestime.3")),
                            v.names = c("iqc_t1_qc_scores","iqc_seriestime"),
                            idvar = "subjectkey", ids = 1:NROW(T1_QC),timevar="acquisition",
                            times=1:3,direction="long")
  row.names(T1_QC_long) <- NULL
  
  T1_QC_long$image_description="ABCD-T1"
  
  ####################################################################################
  #NOTE: some images appeated in QC file but not in image_info, and some vice-versa.
  #The sections below try to first remove data without QC scores (from either source)
  #and then check the images that have a passing QC score.
  ####################################################################################
  group_nn<-dim(data.frame(unique(T1_QC_long$subjectkey)))[1]#num subs before removal of missing T1 QC
  T1_QC_long<-subset(T1_QC_long,T1_QC_long$iqc_t1_qc_score!="")#Remove images with missing QC score (from QC file)
  group_n<-dim(data.frame(unique(T1_QC_long$subjectkey)))[1]#num subs before after of missing T1 QC
  cat('There are',group_nn-group_n,'subjects with missing T1 QC scores','\n')
  
  
  #Floor series times for T1 scans from QC file (because some have decimal places for some reason ¯\_(ツ)_/¯ )
  T1_QC_long$iqc_seriestime<-floor(as.numeric(as.vector(T1_QC_long$iqc_seriestime)))
  
  #Reshape DWI QC file
  DWI_QC_long<-reshape(DWI_QC, varying = list(c("iqc_dmri_qc_scores.1", "iqc_dmri_qc_scores.2","iqc_dmri_qc_scores.3","iqc_dmri_qc_scores.4","iqc_dmri_qc_scores.5","iqc_dmri_qc_scores.6"),
                                                                c("iqc_dmri_seriestime.1", "iqc_dmri_seriestime.2","iqc_dmri_seriestime.3","iqc_dmri_seriestime.4","iqc_dmri_seriestime.5","iqc_dmri_seriestime.6"),
                                                                c("iqc_dmri_fm_missing.1", "iqc_dmri_fm_missing.2","iqc_dmri_fm_missing.3","iqc_dmri_fm_missing.4","iqc_dmri_fm_missing.5","iqc_dmri_fm_missing.6")),
                                v.names = c("iqc_dmri_qc_scores","iqc_seriestime","iqc_dmri_fm_missing"),
                                idvar = "subjectkey", ids = 1:NROW(DWI_QC),timevar="acquisition",
                                times=1:6,direction="long")
  row.names(DWI_QC_long) <- NULL
  DWI_QC_long<-DWI_QC_long[,c('subjectkey','acquisition','iqc_dmri_qc_scores','iqc_seriestime','iqc_dmri_fm_missing')]
  DWI_QC_long$image_description="ABCD-DTI"

  
  group_nn<-dim(data.frame(unique(DWI_QC_long$subjectkey)))[1]#num subs before removal of missing DWI QC
  DWI_QC_long<-subset(DWI_QC_long,DWI_QC_long$iqc_dmri_qc_score!="")#Remove images missing QC score (from QC file)
  group_n<-dim(data.frame(unique(DWI_QC_long$subjectkey)))[1]#num subs before after of missing DWI QC
  cat('There are',group_nn-group_n,'subjects with missing DWI QC scores','\n')
  
  FINAL_inconsistencies$iqc_seriestime<-as.numeric(as.vector(FINAL_inconsistencies$iqc_seriestime))
  
  #Merge in main table the T1 QC scores
  FINAL_inconsistencies_tmp<-merge(FINAL_inconsistencies,T1_QC_long,all=T)
  
  #Floor series times for T1 scans from QC file (because some have decimal places for some reason ¯\_(ツ)_/¯ )
  DWI_QC_long$iqc_seriestime<-floor(as.numeric(as.vector(DWI_QC_long$iqc_seriestime)))
  
  FINAL_inconsistencies$iqc_seriestime<-as.numeric(as.vector(FINAL_inconsistencies$iqc_seriestime))
  
  #Merge in main table the DWI QC scores
  FINAL_inconsistencies_tmp<-merge(FINAL_inconsistencies_tmp,DWI_QC_long,all=T)
  
  #Because the ABCD folks misnamed a few scans (gave them filenames with slightly incorrect 'series time'), this merge has
  #a few mismatches (rows mostly populated by NA). I will take those mismatches to "correct" the real 'series time' (i.e.: make
  #it match the incorrect filename) and then I will merge again.
  
  for (i in 1:dim(FINAL_inconsistencies_tmp)[1]){
    #If image id is missing (all rows from image_info have it), it means it was a mismatched row
    if (is.na(FINAL_inconsistencies_tmp[i,]$image03_id)==T){
      #Find all data for this subject
      tmp1<-subset(FINAL_inconsistencies_tmp,subjectkey==FINAL_inconsistencies_tmp$subjectkey[i])
      #Find all scans of a given image_description
      scans<-subset(tmp1,image_description==FINAL_inconsistencies_tmp[i,]$image_description)
      
      #Loop through all unique series times
      for (l in 1:dim(data.frame(unique(scans$iqc_seriestime)))[1]){
        #Find instances of rows that refer to the same scan
        #These are identified by taking the difference between the series time from image_info and the one from QC file
        #The assumption is that if there's a typo, this difference will be small (arbitrarily decided it would be about 50)
        same_scans<-scans[(as.numeric(as.vector(scans$iqc_seriestime))-as.numeric(as.vector(FINAL_inconsistencies_tmp[i,]$iqc_seriestime))<50&
                                as.numeric(as.vector(scans$iqc_seriestime))-as.numeric(as.vector(FINAL_inconsistencies_tmp[i,]$iqc_seriestime))>-50),]
        
        #If there are duplicated rows
        if (dim(same_scans)[1]>1){
          #Take one row that contains that duplicated data. This one will be the row from QC file (remember first if statement)
          input_idx<-row.names(FINAL_inconsistencies_tmp[i,])
          
          #Find the row, from same_scans, that also contains this duplicated data (but that it's not the same as input_idx, so should be the one from image_info)
          target_idx<-row.names(same_scans)[row.names(same_scans)!=input_idx]
          
          #Now replace image_info row's series time and acquisition with the one from the QC file (which has the correct one)
          FINAL_inconsistencies_tmp[target_idx,]$iqc_seriestime<-FINAL_inconsistencies_tmp[input_idx,]$iqc_seriestime
          FINAL_inconsistencies_tmp[target_idx,]$acquisition<-FINAL_inconsistencies_tmp[input_idx,]$acquisition
          FINAL_inconsistencies_tmp[target_idx,]$iqc_t1_qc_scores<-FINAL_inconsistencies_tmp[input_idx,]$iqc_t1_qc_scores
          FINAL_inconsistencies_tmp[target_idx,]$iqc_dmri_qc_scores<-FINAL_inconsistencies_tmp[input_idx,]$iqc_dmri_qc_scores
          FINAL_inconsistencies_tmp[target_idx,]$iqc_dmri_fm_missing<-FINAL_inconsistencies_tmp[input_idx,]$iqc_dmri_fm_missing
        }
      }
    }
  }
  
  #Now that image_info series times have been corrected, remove the rows with NA
  FINAL_inconsistencies_tmp<-subset(FINAL_inconsistencies_tmp,is.na(image03_id)==F)
  #NOTE: I'm using image03_id because this variable should be present in all rows that come from image_info. 
  #If it's not, it signals a row from the QC table that was not properly matched.
  
  FINAL_inconsistencies_tmp$QC_missing<-0
  for (i in 1:dim(FINAL_inconsistencies_tmp)[1]){
    if ((FINAL_inconsistencies_tmp[i,]$image_description=='ABCD-T1'&is.na(FINAL_inconsistencies_tmp[i,]$iqc_t1_qc_scores)==T) |
        (FINAL_inconsistencies_tmp[i,]$image_description=='ABCD-DTI'& is.na(FINAL_inconsistencies_tmp[i,]$iqc_dmri_qc_scores)==T)){
      #flag if there is a T1 that doesn't have T1-QC score and same thing with DWI.
      FINAL_inconsistencies_tmp[i,]$QC_missing<-1
      #NOTE: here I'm flagging rows from image_info with missing data (because some were not even listed in QC file)
      }
  }
  
  missing_QC<-subset(FINAL_inconsistencies_tmp,QC_missing==1)
  
  subs_missing_QC<-data.frame(unique(missing_QC$subjectkey))
  dim(subs_missing_QC)[1]
  # For your own use if you want to see how many subs had some QC data missing
  
  #Remove images that didn't have QC data.
  FINAL_inconsistencies_QC<-subset(FINAL_inconsistencies_tmp,FINAL_inconsistencies_tmp$QC_missing==0) 
  
  #How many images failed QC
  failed_T1_QC<-subset(FINAL_inconsistencies_QC,iqc_t1_qc_scores==0)
  # dim(failed_T1_QC)#For your own use
  #How many subs had a T1 QC fail
  subs_T1_QC_fail<-unique(failed_T1_QC$subjectkey)
  # dim(data.frame(subs_T1_QC_fail))#For your own use
  
  #How many images failed QC
  failed_DWI_QC<-subset(FINAL_inconsistencies_QC,iqc_dmri_qc_scores==0)
  # dim(failed_DWI_QC)#For your own use
  
  #How many subs had a DWI QC fail
  subs_DWI_QC_fail<-unique(failed_DWI_QC$subjectkey)
  # dim(data.frame(subs_DWI_QC_fail))#For your own use

  #verify if you got rid of mismatched by searching for ABCD_T1 = NA (just me being paranoid)
  tmp<-FINAL_inconsistencies_QC[is.na(FINAL_inconsistencies_QC$ABCD_T1)==T,]
  # table(is.na(FINAL_inconsistencies_QC$ABCD_T1))
  
  
  #Remove DWI and T1 scans that didn't pass QC
  FINAL_data_QCPass_tmp<-subset(FINAL_inconsistencies_QC,iqc_t1_qc_scores==1|is.na(iqc_t1_qc_scores)==T)
  #NOTE: Take people who passed T1 QC or who had NA (now, NA will only be inrows where T1 QC is irrelevant (e.g.: for a DWI image))

  FINAL_data_QCPass<-subset(FINAL_data_QCPass_tmp,iqc_dmri_qc_scores==1|is.na(iqc_dmri_qc_scores)==T)
  #Same thing for DWI
  
  FINAL_data <- FINAL_data_QCPass
  
  #Now that we've removed rows corresponding to images that had missing QC scores or had failed QC, we will count again which subjects have complete data.
  tmp <- FINAL_data %>% group_by(subjectkey) %>% summarise(DTI_FM_verif = sum(image_description == 'ABCD-Diffusion-FM'), DTI_FM_PA_verif = sum(image_description == 'ABCD-Diffusion-FM-PA'), DTI_FM_AP_verif = sum(image_description == 'ABCD-Diffusion-FM-AP'),
                                                                                        ABCD_DTI_verif = sum(image_description == 'ABCD-DTI'), ABCD_T1_verif = sum(image_description =='ABCD-T1' ))
  Verif <- merge(tmp, FINAL_data, by='subjectkey', all=T)
  
  PHILIPS_verif <- subset(Verif, scanner_manufacturer_pd=='Philips Medical Systems')
  PHILIPS_verif$missingdata_verif <- ifelse(PHILIPS_verif$ABCD_DTI_verif<2|PHILIPS_verif$DTI_FM_PA_verif<1|PHILIPS_verif$DTI_FM_AP_verif<1|PHILIPS_verif$ABCD_T1_verif<1, 1,0)
  PHILIPS_verif$T1REPEAT_verif <- ifelse(PHILIPS_verif$ABCD_T1_verif>1, 1, 0)
  PHILIPS_verif$DWIREPEAT_verif <- ifelse(PHILIPS_verif$ABCD_DTI_verif>2 | PHILIPS_verif$DTI_FM_PA_verif>1| PHILIPS_verif$DTI_FM_AP_verif>1, 1,0)
  
  SIEMENS_verif <- subset(Verif, scanner_manufacturer_pd=='SIEMENS')
  SIEMENS_verif$missingdata_verif <- ifelse(SIEMENS_verif$ABCD_DTI_verif<1|SIEMENS_verif$DTI_FM_PA_verif<1|SIEMENS_verif$DTI_FM_AP_verif<1|SIEMENS_verif$ABCD_T1_verif<1, 1,0)
  SIEMENS_verif$T1REPEAT_verif <- ifelse(SIEMENS_verif$ABCD_T1_verif>1, 1, 0)
  SIEMENS_verif$DWIREPEAT_verif <- ifelse(SIEMENS_verif$ABCD_DTI_verif>1 | SIEMENS_verif$DTI_FM_PA_verif| SIEMENS_verif$DTI_FM_AP_verif>1, 1,0)
  
  GE_verif <- subset(Verif, scanner_manufacturer_pd=='GE MEDICAL SYSTEMS')
  GE_verif$missingdata_verif <- ifelse(GE_verif$ABCD_DTI_verif<1|GE_verif$DTI_FM_verif<1|GE_verif$ABCD_T1_verif<1, 1,0)
  GE_verif$T1REPEAT_verif <- ifelse(GE_verif$ABCD_T1_verif>1, 1, 0)
  GE_verif$DWIREPEAT_verif <- ifelse(GE_verif$ABCD_DTI_verif>1 | GE_verif$DTI_FM_verif>1, 1,0)
  
  
  
  FINAL_inconsistencies_verif <- rbind(GE_verif,SIEMENS_verif,PHILIPS_verif) #Here's the important table for the next steps.
  
  subs_missingdata_verif<-data.frame(unique(FINAL_inconsistencies_verif$subjectkey[FINAL_inconsistencies_verif$missingdata_verif==1]))
  colnames(subs_missingdata_verif)<-'subjectkey'
  
  FINAL_data_clean<-FINAL_data[! FINAL_data$subjectkey %in% subs_missingdata_verif$subjectkey, ]
  
  group_nn<-dim(data.frame(subs_missingdata_verif))[1]#should be 350 (mtbi) or 460 (Match Controls)
  group_n<-dim(data.frame(unique(FINAL_data_clean$subjectkey)))[1]#mTBI: 339 (consistent with what I had before) MC: 433 subs (460-27)
  cat('There are',group_n,'subjects with available and passing QC data (',group_nn,'missing or failed)','\n')
  
  #save the S3 links in a text file
  write.table(FINAL_data_clean$image_file,paste(args[1],group_dir,'/S3links.txt',sep=''),sep = '\n',row.names = FALSE, col.names = FALSE, quote=FALSE)
  

#########SECTION THAT WILL BE RELEVANT TO A QC SCRIPT I'M WRITING###########
  #Text file containing IDs and scanner type
  # FINAL_subs_scanners<-FINAL_data_clean[,c('subjectkey','scanner_manufacturer_pd')]
  # write.table(FINAL_subs_scanners,paste(args[1],'/mTBI_Group/Final_subs_scanners.txt',sep = ' ',eol='\n',row.names = FALSE, col.names = FALSE, quote=FALSE)
  # write.table(FINAL_subs_scanners,"/Users/Guido/Desktop/GeneralData/ABCD/mTBI_Group/Final_subs_scanners.txt",sep = ' ',eol='\n',row.names = FALSE, col.names = FALSE, quote=FALSE)
# Tensor_Peak_error <- read.table('/Users/Guido/Desktop/GeneralData/ABCD/Control_Group_Harmonization/Data/BIDS/final_data_for_upload/MC_subs_weird_errors.txt', sep='\t', header=T)
# colnames(Tensor_Peak_error)<-'subjectkey'
