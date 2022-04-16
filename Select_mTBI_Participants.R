##########################################################
##        SCRIPT TO SELECT mTBI PARTICIPANTS              ##
##########################################################
args = commandArgs(trailingOnly=TRUE) # Line to get user input

#ABCD TBI Summary Scores
tbi_summary <- read.table(paste(args[1],'/GeneralData/abcd_tbi01.txt',sep = ''), sep='\t', header=T)
#tbi_summary <- read.table('/Users/Guido/Desktop/GeneralData/ABCD/GeneralData/abcd_tbi01.txt', sep='\t', header=T)#for troubleshooting
datadic_tbi_summary <- tbi_summary[1,]
datadic_tbi_summary <-t(datadic_tbi_summary)
tbi_summary <-tbi_summary[-1,]

# Select those who have mild TBI as their worst injury 
###############################################################################
##                          TBI CRITERIA                                     ##
##    1. Improbable TBI (no TBI or TBI w/o LOC or memory loss)               ##
##    2. Possible mild TBI (TBI w/o LOC but memory loss)                     ##
##    3. Mild TBI (TBI w/LOC ≤ 30 min)                                       ##
##    4. Moderate TBI (TBI w/LOC  30 min - 24 hrs)                           ##
##    5. Severe TBI (TBI w/ LOC ≥ 24 hrs)                                    ##
###############################################################################

mTBI <- subset(tbi_summary, tbi_ss_worst_overall==2 | tbi_ss_worst_overall==3 )


#ABCD Parent Medical History Questionnaire
parent_MEDHIST <- read.table(paste(args[1],'/GeneralData/abcd_mx01.txt',sep = ''), sep='\t', header=T)
# parent_MEDHIST <- read.table('/Users/Guido/Desktop/GeneralData/ABCD/GeneralData/abcd_mx01.txt', sep='\t', header=T)#for troubleshooting

############################################################
##              E X C L U S I O N   C R I T               ##
##                Epilsepsy                               ##
##                Lead Poisoning                          ##
##                Multiple Sclerosis                      ##
##                Cerebral Palsy                          ##
############################################################
mTBI_medhx <- merge (mTBI, parent_MEDHIST, by = 'subjectkey', all.x=T )

## EPILEPSY
mTBI_medhx$epilepsy <- (as.numeric(as.character(mTBI_medhx$medhx_2h)))

## LEAD POISONING
mTBI_medhx$leadpoisoning <- (as.numeric(as.character(mTBI_medhx$medhx_2k)))

## MULTIPLE SCLEROSIS
mTBI_medhx$multscler <- (as.numeric(as.character(mTBI_medhx$medhx_2m)))

## CEREBRAL PALSY
mTBI_medhx$cerebralpals <- (as.numeric(as.character(mTBI_medhx$medhx_2f)))

## remove participants with these diagnoses
mTBI_medhx_sub1 <- subset(mTBI_medhx, epilepsy==0)
mTBI_medhx_sub2 <- subset(mTBI_medhx_sub1, leadpoisoning==0)
mTBI_medhx_sub3 <- subset(mTBI_medhx_sub2, multscler==0)
mTBI_medhx_excluded <- subset(mTBI_medhx_sub3, cerebralpals==0)


group_n<-dim(data.frame(mTBI_medhx_excluded$subjectkey))[1]
cat('There are',group_n,'mTBI subjects','\n')

write.table(mTBI_medhx_excluded$subjectkey,paste(args[1],'/mTBI_Group/IDs.txt',sep=''),row.names = F, col.names = F, quote = F)
# write.table(mTBI_medhx_excluded$subjectkey,'/Users/Guido/Desktop/GeneralData/ABCD/mTBI_Group/IDs.txt',row.names = F, col.names = F, quote = F)#for troubleshooting