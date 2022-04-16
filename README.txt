########################################################################################
            Multi-Tract Multi-Symptom Relationships in Pediatric Concussion
      Guberman GI, Stojanovski S, Nishat E, Ptito A, Bdzok D, Wheeler A, Descoteaux M

                                  _---~~(~~-_.__
                                _{        )     )
                              ,   ) -~~- ( ,-' ' )_
                             (  ''-,_..''., )-- '_,)
                            ( '' _)  (  -~( -_ ',  }
                            (_-  _  ~_-~~~~',  ,' )
                              '~ -^(    __;-,((()))
                                    ~~~~ {_ -_(())
                                           '\  }
                                             { }      SJW


#########################################################################################

This README describes each of the scripts that were used in the Multi-tract multi-symptom relationships in pediatric concussions study.

When using this code, please cite:
Guberman, G. I., Stojanovski, S., Nishat, E., Ptito, A., Bzdok, D., Wheeler, A., & Descoteaux, M. (2022). Multi-tract multi-symptom relationships in pediatric concussions. eLife.


Overall structure: this project runs mainly using two large scripts, and a few supporting scripts
1. MasterScript.sh: this bash script handles most of the data organizing and cleaning. The requirements are as follows:

Data from ABCD:
-gen_dir/GeneralData
   |- abcd_tbi01.txt
   |- abcd_mx01.txt
   |- image03.txt
   |- abcd_mri01.txt
   |- abcd_ehis01.txt
   |- abcd_tbss01.txt
   |- abcd_ps01.txt
   |- mriqcrp102.txt
   |- abcd_ssphp01.txt
   (note: see MasterScript.doc for details)

R PACKAGES:
 - MatchIt
 - dplyr

SCRIPTS (added to path):
 - organize_ABCD_data.sh
 - Select_mTBI_Participants.R
 - Select_HC_Participants.R
 - flag_data_inconsistencies.R
 - mutual_information.py
 - Extract_DWIandT1Scans_s3links_Official.R
 - downloadcmd (see MasterScript.doc for details)

OTHER TOOLS:
 - SCILPY (https://github.com/scilus/scilpy)
 - MRtrix3 (https://www.mrtrix.org)
 - FSL (https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/)
 - Rstudio (https://www.rstudio.com)
 - NDA Download Manager (see MasterScript.doc for details)

2. MTMSAnalysesOfficial.R: this R script handles the post-processing and statistical analysis. All packages and libraries are installed and loaded within the script.

As input, it requires:

- Data processed using Tractoflow (Theaud, G., Houde, J. C., Boré, A., Rheault, F., Morency, F., & Descoteaux, M. (2020). TractoFlow:
A robust, efficient and reproducible diffusion MRI pipeline leveraging Nextflow & Singularity. NeuroImage, 218, 116889) implemented in a ComputeCanada cluster. The script to call Tractoflow is called process_ABCD_mTBI.sh.

- DKT parcellations (Klein A, Tourville J. 101 labeled brain images and a consistent human cortical labeling protocol. Front Neurosci. 2012;6:171) fitted with Freesurfer implemented in McGill's C-Brain platform.
The file is called DKT_lbls.txt (provided in Github)


- Connectivity matrices built using Connectoflow (Rheault, F., Houde, J. C., Sidhu, J., Obaid, S., Guberman, G., Daducci, A., & Descoteaux, M. Connectoflow: A cutting-edge Nextflow pipeline for structural connectomics),
implemented in a ComputeCanada cluster. The script to call Connectoflow is called call_connectoflow.sh. This script also uses an MNI152_T1_1mm_brain.nii template available online.

- TBI_QC_fail_final.txt : a txt file with the subject IDs of individuals that were deemed to have failed visual quality control after processing. Given that this file contains ABCD data, it cannot be shared publicly without permission from the ABCD Study team

(from ABCD)
- abcd_cbcls01.txt
- abcd_cbcl01.txt
- abcd_sds01.txt
- abcd_ps01.txt
- abcd_tbss01.txt
- abcd_mri01.txt
- abcd_tbi01.txt
- abcd_otbi01.txt
- abcd_ssphp01.txt
- abcd_ehis01.txt
- abcd_ksad01.txt
- pdem02.txt
