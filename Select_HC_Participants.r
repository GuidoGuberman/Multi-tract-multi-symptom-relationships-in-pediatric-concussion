##########################################################
##        SCRIPT TO SELECT HC PARTICIPANTS              ##
##########################################################

args = commandArgs(trailingOnly=TRUE)
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



############################################################
##                    READ IN ALL DATA                    ##
############################################################

allimage <- read.table(paste(args[1],'/GeneralData/image03.txt',sep=''), sep='\t', header=T)
# allimage <- read.table('/Users/Guido/Desktop/GeneralData/ABCD/GeneralData/image03.txt', sep='\t', header=T)
allimage <-allimage[-1,]

## TBI Info
TBI <- read.table(paste(args[1],'/GeneralData/abcd_tbi01.txt',sep=''), sep='\t', header=T)
# TBI <- read.table('/Users/Guido/Desktop/GeneralData/ABCD/GeneralData/abcd_tbi01.txt', sep='\t', header=T)
datadic <- TBI[1,]
datadic <-t(datadic)
TBI <-TBI[-1,]

###############################################################################
##                          TBI CRITERIA                                     ##
##    1. Improbable TBI (no TBI or TBI w/o LOC or memory loss)               ##
##    2. Possible mild TBI (TBI w/o LOC but memory loss)                     ##
##    3. Mild TBI (TBI w/LOC ≤ 30 min)                                       ##
##    4. Moderate TBI (TBI w/LOC  30 min - 24 hrs)                           ##
##    5. Severe TBI (TBI w/ LOC ≥ 24 hrs)                                    ##
###############################################################################
TBI <- TBI[c('subjectkey', 'tbi_ss_worst_overall')]

## ABCD Parent Medical History Questionnaire
parent_MEDHIST <- read.table(paste(args[1],'/GeneralData/abcd_mx01.txt',sep=''), sep='\t', header=T)
# parent_MEDHIST <- read.table('/Users/Guido/Desktop/GeneralData/ABCD/GeneralData/abcd_mx01.txt', sep='\t', header=T)
# parent_MEDHIST <- read.table('/Users/atepavac/Documents/SONJA/ABCD/metadata/raw/MEDICALHISTORYPARENT/abcd_mx01.txt', sep='\t', header=T)

## Site information
MRI <- read.table(paste(args[1],'/GeneralData/abcd_mri01.txt',sep=''), sep='\t', header=T)
# MRI <- read.table('/Users/Guido/Desktop/GeneralData/ABCD/GeneralData/abcd_mri01.txt', sep='\t', header=T)
# MRI <- read.table('/Users/atepavac/Documents/SONJA/ABCD/metadata/raw/MRIINFO/abcd_mri01.txt', sep='\t', header=T)
MRI <- MRI[c('src_subject_id', 'mri_info_deviceserialnumber')]
#NOTE: there are people with no device serial number ...?

## HANDEDNESS
hand <- read.table(paste(args[1],'/GeneralData/abcd_ehis01.txt',sep=''), sep='\t', header=T)
# hand <- read.table('/Users/Guido/Desktop/GeneralData/ABCD/GeneralData/abcd_ehis01.txt', sep='\t', header=T)
# hand <- read.table('/Users/atepavac/Documents/SONJA/ABCD/metadata/raw/HANDEDNESS/abcd_ehis01.txt', sep='\t', header=T)
hand <- hand[c('subjectkey', 'ehi_y_ss_scoreb')]

## IQ
# Given recent controversies regarding the association between substance misuse and IQ decline over time (Jackson et al., 2016; Meier et al., 2012), 
#     the Consortium recognized the importance of including measures of baseline fluid (Matrix Reasoning) 
#     and crystallized (NIH Toolbox Picture Vocabulary™) reasoning with the goal of re-assessing these skills over time.
# https://www.sciencedirect.com/science/article/pii/S1878929317302384

## nihtbx_picvocab_agecorrected - pearson wisc matrix reasoning total scaled score - PROXY FOR FLUID IQ  (https://www.sciencedirect.com/science/article/pii/S1878929317302384)
# cryst_iq <- read.table('/Users/atepavac/Documents/SONJA/ABCD/metadata/raw/IQ_CRYSTALIZED/abcd_tbss01.txt', sep='\t', header=T)
cryst_iq <- read.table(paste(args[1],'/GeneralData/abcd_tbss01.txt',sep=''), sep='\t', header=T)
# cryst_iq <- read.table('/Users/Guido/Desktop/GeneralData/ABCD/GeneralData/abcd_tbss01.txt', sep='\t', header=T)
cryst_iq <- cryst_iq[c('subjectkey', 'nihtbx_picvocab_agecorrected')]
cryst_iq$IQ_crystalized <- as.numeric(as.character(cryst_iq$nihtbx_picvocab_agecorrected))
#NOTE: Removed people with missing data


## pea_wiscv_tss - pearson wisc matrix reasoning total scaled score - PROXY FOR FLUID IQ  (https://www.sciencedirect.com/science/article/pii/S1878929317302384)
# fluid_iq <- read.table('/Users/atepavac/Documents/SONJA/ABCD/metadata/raw/IQ_FLUID/abcd_ps01.txt', sep='\t', header=T)
fluid_iq <- read.table(paste(args[1],'/GeneralData/abcd_ps01.txt',sep=''), sep='\t', header=T)
# fluid_iq <- read.table('/Users/Guido/Desktop/GeneralData/ABCD/GeneralData/abcd_ps01.txt', sep='\t', header=T)
fluid_iq <- fluid_iq[c('subjectkey', 'pea_wiscv_tss')]
fluid_iq$IQ_fluid <- as.numeric(as.character(fluid_iq$pea_wiscv_tss))
#NOTE: Removed people with missing data

pubertal_stage <- read.table(paste(args[1],'/GeneralData/abcd_ssphp01.txt',sep=''), sep='\t', header=T)
# pubertal_stage <- read.table('/Users/Guido/Desktop/GeneralData/ABCD/GeneralData/abcd_ssphp01.txt', sep='\t', header=T)
pubertal_stage<-pubertal_stage[-1,]
pubertal_stage<-subset(pubertal_stage,pubertal_stage$eventname=='baseline_year_1_arm_1')
pubertal_stage <- pubertal_stage[c('subjectkey', 'pds_p_ss_female_category','pds_p_ss_male_category')]
# subset(pubertal_stage,pubertal_stage=="2 1")#somebody has this value
# subset(pubertal_stage,subjectkey=="NDAR_INV00HEV6HB")#it's this person

pubertal_stage$pubertal_stage<-paste(pubertal_stage$pds_p_ss_female_category,pubertal_stage$pds_p_ss_male_category)
pubertal_stage$pubertal_stage<-as.numeric(as.character(pubertal_stage$pubertal_stage))
pubertal_stage <- pubertal_stage[c('subjectkey', 'pubertal_stage')]

table(pubertal_stage$pubertal_stage)
#NOTE: Removed people with missing data

#25/02/2020: Comparing age against AFD (REMOVE ONCE DONE)


afds <- read.table('/Users/Guido/Desktop/GeneralData/ABCD/Harmonization_Controls/Data/BIDS/Tests_of_AFD/b3000/metrics/afds.txt', sep=' ', header=F)
colnames(afds)<-c('subjectkey','NA','afd_max_b3000','afd_max_b3000_val','afd_max_ms','afd_max_ms_val')
tmp<-switch_names(afds$subjectkey,Rem = "_ses",Add = "",separ = '',side = 'R')
tmp<-data.frame(tmp)
colnames(tmp)<-c('subjectkey')
afds$subjectkey<-switch_names(tmp$subjectkey,Rem='sub-NDAR',Add='NDAR_',separ='',side='L')

afds<-merge(afds,pubertal_stage,by = 'subjectkey')
afds<-afds[!afds$subjectkey=='NDAR_INV06DX0E7R',]


p<-cor(afds[,c('afd_max_b3000_val','afd_max_ms_val','pubertal_stage')])
plot(afds[,c('afd_max_b3000_val','pubertal_stage')])
abline(lm(afds$pubertal_stage~afds$afd_max_b3000_val))
plot(afds[,c('afd_max_ms_val','pubertal_stage')])
abline(lm(afds$pubertal_stage~afds$afd_max_ms_val))





#############################################################
##                    CLEAN IMAGE INFO                     ##
##               RISH HARM REQUIRES DWI ONLY               ##
##        BUT TO PROCESS DWI, TRACTOFLOW STILL NEEDS       ##
##            ALL OTHER IMAGES (T1, REV_B0)                ##
#############################################################
scans <- subset(allimage, image_description== 'ABCD-Diffusion-FM-PA' | image_description== 'ABCD-Diffusion-FM-AP' | 
                             image_description== 'ABCD-DTI' | image_description== 'ABCD-Diffusion-FM'| image_description == 'ABCD-T1')
scans <- subset(scans, visit=='baseline_year_1_arm_1' )

## seperate out scanner type to use for the function later
scanner <- scans[c('subjectkey', 'scanner_manufacturer_pd')]
scanner <- unique(scanner)

# for this to work I'm going to need to rotate the axis ... from more than one value for image_description along the y axis to many columns along the x axis to then apply a function
## so yea take DWI_scans and hit it with tidyr ... spread? idk. who knows? not me.)
library(tidyr)
library(dplyr)
# install.packages("MatchIt")
library(MatchIt)


temp <- scans %>% group_by(subjectkey) %>% summarise(DTI_FM = sum(image_description == 'ABCD-Diffusion-FM'), DTI_FM_PA = sum(image_description == 'ABCD-Diffusion-FM-PA'), ABCD_DTI = sum(image_description == 'ABCD-DTI'), DTI_FM_AP = sum(image_description == 'ABCD-Diffusion-FM-AP'))
FINAL <- merge(temp, scanner, by='subjectkey', all=T)

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

PHILIPS <- subset(FINAL, scanner_manufacturer_pd=='Philips Medical Systems')
PHILIPS$missingdata <- ifelse(PHILIPS$ABCD_DTI<2|PHILIPS$DTI_FM_PA<1|PHILIPS$DTI_FM_AP<1, 1,0)


SIEMENS <- subset(FINAL, scanner_manufacturer_pd=='SIEMENS')
SIEMENS$missingdata <- ifelse(SIEMENS$ABCD_DTI<1|SIEMENS$DTI_FM_PA<1|SIEMENS$DTI_FM_AP<1, 1,0)


GE <- subset(FINAL, scanner_manufacturer_pd=='GE MEDICAL SYSTEMS')
GE$missingdata <- ifelse(GE$ABCD_DTI<1|GE$DTI_FM<1, 1,0)
allimage_withdata <- rbind(GE,SIEMENS,PHILIPS)
allimage_withdata <- subset(allimage_withdata, missingdata==0 )


#NOTE: it might seem redundant to check for missing data twice, but here it was done so as to not match people who will be later removed.
#The script that comes after also checks for missing data, but additionally removes data with missing or failing QC.

############################################################
##                    MERGE ALL DATA                      ##
############################################################
temp <- merge (allimage_withdata, TBI, by = 'subjectkey', all.x=T )
temp <- merge(temp, MRI, by.x = 'subjectkey', by.y='src_subject_id', all.x = T)
temp <- merge(temp, hand, by.x = 'subjectkey', by.y='subjectkey', all.x = T)
temp <- merge(temp, cryst_iq, by.x = 'subjectkey', by.y='subjectkey', all.x = T)
temp <- merge(temp, fluid_iq, by.x = 'subjectkey', by.y='subjectkey', all.x = T)
temp <- merge(temp, pubertal_stage, by.x = 'subjectkey', by.y='subjectkey', all.x = T)
TD <- merge (temp, parent_MEDHIST, by = 'subjectkey', all.x=T )



############################################################
##              E X C L U S I O N   C R I T               ##
##                Epilsepsy                               ##
##                Lead Poisoning                          ##
##                Multiple Sclerosis                      ##
##                Cerebral Palsy                          ##
##                TBI                                     ##
############################################################

## EPILEPSY
TD$epilepsy <- (as.numeric(as.character(TD$medhx_2h)))

## LEAD POISONING
TD$leadpoisoning <- (as.numeric(as.character(TD$medhx_2k)))

## MULTIPLE SCLEROSIS
TD$multscler <- (as.numeric(as.character(TD$medhx_2m)))

## CEREBRAL PALSY
TD$cerebralpals <- (as.numeric(as.character(TD$medhx_2f)))

## remove participants with these diagnoses
TD_sub1 <- subset(TD, epilepsy==0)
TD_sub2 <- subset(TD_sub1, leadpoisoning==0)
TD_sub3 <- subset(TD_sub2, multscler==0)
TD_sub4 <- subset(TD_sub3, cerebralpals==0)
## only keep particiapnts without TBI
TD_sub5 <- subset(TD_sub4, tbi_ss_worst_overall==1)
## remove participants who don't have a serial number
TD_sub6 <- subset(TD_sub5, is.na(mri_info_deviceserialnumber)==FALSE)
## remove particiapnts who don't have information on handedness
TD_sub7 <- subset(TD_sub6, is.na(ehi_y_ss_scoreb)=='FALSE')
## remove particiapnts who don't have information on IQ
## more people are missing information on fluid IQ than crystallized, so this is what we will use.
TD_sub8 <- subset(TD_sub7, is.na(pubertal_stage)=='FALSE')
TD_full <- subset(TD_sub8, is.na(IQ_crystalized)=='FALSE')

# TD <- TD_full[c('subjectkey', 'mri_info_deviceserialnumber', 'interview_age', 'sex', 'ehi_y_ss_scoreb', 'IQ_crystalized')]
TD <- TD_full[c('subjectkey', 'mri_info_deviceserialnumber', 'pubertal_stage', 'sex', 'ehi_y_ss_scoreb', 'IQ_crystalized')] #Trying to match by pubertal stage this time

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
TD$group <- ifelse(TD$mri_info_deviceserialnumber=='HASHe76e6d72', 1, 0)#######################17/01/2020: here's the place where I will likely have to modify
target<-data.frame(subset(TD,TD$group==1))

## force age to be a number!
# TD$interview_age <- as.numeric(as.character(TD$interview_age))

all_sites<-unique(data.frame(TD$mri_info_deviceserialnumber))
# dim(all_sites)
#Note: there will be 30, because some people don't have a serial number

colnames(all_sites)<-'sites'

#remove target from the list of sites
match_sites<-subset(all_sites,all_sites$sites!='HASHe76e6d72')

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
#NOTE: the 6 original variables in TD + the three that matchit always adds (group, distance, weights)

output <- data.frame(matrix(NA,ncol=variables, nrow=full_size))
# colnames(output)<-c('subjectkey','mri_info_deviceserialnumber','interview_age','sex','ehi_y_ss_scoreb','IQ_crystalized','group','distance','weights')
colnames(output)<-c('subjectkey','mri_info_deviceserialnumber','pubertal_stage','sex','ehi_y_ss_scoreb','IQ_crystalized','group','distance','weights')

TD<-na.omit(TD)

for(i in 1:dim(match_sites)[1]){
  #Pool of data from which to match (contains target site data + one other site's data)
  match_pool <- subset(TD, mri_info_deviceserialnumber=='HASHe76e6d72' | mri_info_deviceserialnumber==match_sites[i,])
  # match.it <- matchit(group~ interview_age + sex + ehi_y_ss_scoreb + IQ_crystalized , data=match_pool, method='nearest', ratio=as.numeric(args[3]))
  # # match.it <- matchit(group~ interview_age + sex + ehi_y_ss_scoreb + IQ_crystalized , data=match_pool, method='nearest', ratio=2) #For troubleshooting
  
  match.it <- matchit(group~ pubertal_stage + sex + ehi_y_ss_scoreb + IQ_crystalized , data=match_pool, method='nearest', ratio=as.numeric(args[3]))
  # match.it <- matchit(group~ pubertal_stage + sex + ehi_y_ss_scoreb + IQ_crystalized , data=match_pool, method='nearest', ratio=2) #For troubleshooting
  
  data<-match.data(match.it)
  
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




#Remove duplicated data (the target site will be duplicated a bunch of times)
DATAFORHARMONIZATION <- output[!duplicated(output$subjectkey), ]
DATAFORHARMONIZATION<-subset(DATAFORHARMONIZATION,is.na(DATAFORHARMONIZATION$subjectkey)==F)

group_n<-dim(data.frame(unique(DATAFORHARMONIZATION$subjectkey)))[1]
cat('There are',group_n,'Harmonization Control subjects','\n')

#Verify that data is matched
res.aov <- aov(sex ~ mri_info_deviceserialnumber, data = DATAFORHARMONIZATION)
sex<-summary(res.aov)[[1]][["Pr(>F)"]][1]

res.aov <- aov(ehi_y_ss_scoreb ~ mri_info_deviceserialnumber, data = DATAFORHARMONIZATION)
hand<-summary(res.aov)[[1]][["Pr(>F)"]][1]

res.aov <- aov(IQ_crystalized ~ mri_info_deviceserialnumber, data = DATAFORHARMONIZATION)
IQ<-summary(res.aov)[[1]][["Pr(>F)"]][1]

res.aov <- aov(pubertal_stage ~ mri_info_deviceserialnumber, data = DATAFORHARMONIZATION)
#Quick look at the pubertal ages
# d<-aggregate(DATAFORHARMONIZATION$pubertal_stage, list(DATAFORHARMONIZATION$mri_info_deviceserialnumber), mean)
# xrange<-dim(d)[1]
# yrange<-max(d$x)
# plot(xrange, yrange, type="n", xlab="Site",
#      ylab="Pubertal stage" )
# barplot(d$x)

pubertal_stage<-summary(res.aov)[[1]][["Pr(>F)"]][1]

if (sex>0.05&hand>0.05&IQ>0.05&pubertal_stage>0.05){
  #If sites are matched, tell the person and just proceed
  cat('Sites are matched!','\n')
} else {
  #If they are not, tell the person and then try randomly selecting a pre-defined number of tries. If that still doesn't work, tell the person to go back and change the initial matching ratio to get a bigger buffer.
  cat('Sites are not matched, you should redo the matching','\n')
}

write.table(DATAFORHARMONIZATION$subjectkey, paste(args[1],group_dir,'/IDs.txt',sep=''), quote = F, col.names = F, row.names = F )

#Table used to check if sites are matched after checking for data inconsistencies (later script)
write.table(DATAFORHARMONIZATION, paste(args[1],group_dir,'/Harmonization_Data.txt',sep=''), quote = F, col.names = F, row.names = F )



