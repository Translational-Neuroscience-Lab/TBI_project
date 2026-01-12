# Seizure detection analysis for TBI study

## Overview: 
This code base supports the Seizure analysis of mouse EEG data. The main MATLAB script (Seizure_detection.m) processes EEG and sleep scoring data to detect epileptiform activity (EA) using a manually selected threshold and plot epochs with suspected EA for manual inspection.

It is tailored for data collected using Sirenia Acquisition (Pinnacle Technology) and assumes:

- EEG is exported in .edf format (24hs of interest)
- Sleep scoring is in .tsv format
- Scoring uses 10-second epochs
- Each mouse is processed separately
- Instructions on how we named our files will be included here

## Input files:

.edf file: Contains raw EEG time series
.tsv file: Contains epoch-level sleep scoring (Wake = 1, NREM = 2, REM = 3, Unscored = 255)

## What the script does:

- Loads EEG data (3 channels - EEGEEG1A_B, EEGEEG2A_B, EMGEMG- at 400Hz sampling rate) and #sleep scores
- Breaks the .edf time series into 10s epochs
- #Aligns sleep scorings with corresponding 10s EEG epochs
- #Applies a bandpass filter (0.5–100 Hz)
- #Removes artifacts (epochs with signal amplitude > 399 µV) - maybe I can try this plus clipping?
- It uses the Hilbert transform function to calculate the mean envelope amplitude for each epoch
- Detect epochs with "peaks" greater than x standard deviations (SD) from the mean (only in channels EEG1 and EEG2)
- Plots the first 6 epochs with EA detected and provides information about seizure epochs in each channel
- Allows the user to decide whether to keep or discard

## How to use it:

- Change the filepath according to the file that needs processing
- Set the threshold accordingly (we recommend a minimum of 5 SD) in the "detectionResults, epochStats" after "fs" and under "visualizeDetectionResults" after "seizures"



