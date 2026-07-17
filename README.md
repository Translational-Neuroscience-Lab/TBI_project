# TBI Project APP NL-F

This repository contains all of the codes used in the processing and analysis of the data corresponding to each published paper so that it can be replicated by anyone who so desires. 

## Description 
The main codes are used for:
* Seizure scoring and analysis
* Power spectral density (PSD) pre-processing and computation from EEG data
* Multiscale Entropy (MSE) curves pre-processing and computation from EEG data
* Mean-centered partial least squares statistical analysis for PSD and MSE

The project has been split into many papers and not all papers use the same codes. Therefore, please refer to the appropriate folder for each paper. Folders are divided as follows:

| Folder Name | Description |
| --- | --- | 
|**APPNL-F_ Subacute**| Contains the codes for the paper titled "Repeated mild traumatic brain injury does not affect sleep or epileptiform activity one-month post-injury in a knock-in mouse model of Alzheimer’s disease" |
|**APPNL-F_ Chronic**| Contains the cordes for the paper including the 6 and 12 months old mice analyzed at 15 months. Title TBD. |
|**APPNL-F_ Acute**| Contains the codes for the paper including the 6 and 12 months old mice analyzed immidiately after the TBI. Title TBD. |
|**APPNL-F_MSE**| Contains all the MSE codes for the MSE paper. Title TBD. |

## Getting Started

You will need Matlab to run all the codes. We used version R2023a but should be able to run in any version.

### PSD and Seizure codes

For the PSD and seizure codes you will need two files for each mouse/subject: the EDF file with all the raw EEG data and the TSV files with the sleep scoring of each epoch. 
*NOTE: we used the Sirenia Sleep software to acquire the EEG data and to score sleep, therefore, our files came from that software. If using another acquisition/scoring software you will still need to convert your data to EDF/TSV files to be compatible with this codes.

### MSE

In construction

### PLS
For the PLS analysis please refer to this repository from our collaborators: https://github.com/McIntosh-Lab/PLS/
You can donwload the PLS package to run our PLS code.

For all other codes, no download or installation is necessary.

The data for the PLS analysis and plotting was obtained from our previous PSD pre-processing code that converted all the aggregate data into excel files. The excel files were organized as shown below with frequencies going all the way to 30Hz:

<img width="738" height="41" alt="image" src="https://github.com/user-attachments/assets/304b61c5-a918-4949-af87-0e44051f231c" />

This requires some manual manipulation of the data before it can be used in the PLS code. The data was separated by vigilance state (NREM, REM, Wake) and by light period (on/off).
*NOTE: the "Sex PLS" that considers both treatment and sex, is constructed so that your data should be organized as such: sham females, sham males, TBI females, TBI males. If altering the ordering of the data you will also need to alter the code.

More details on how to use each code will be available soon under the "README files" folder.
