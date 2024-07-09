# Novel _ex vivo_ Culturing Techniques to Study Spatial Patterns in the Human Intratumoral Immune Response to Immune Checkpoint Blockade


This project focuses on deciphering immune responses to immune checkpoint blockade (ICB) in a spatial manner. To investigate this, we have developed two platforms to achieve this: the patient-derived tumour fragment (PDTF) platform, established in Voabil et al. (2021), and the patient-derived tumour slice (PDTS) platform established in this project. This repository contains the R scripts used to analyse flow cytometry and LEGENDPLEX data gathered in this project.

The scripts in this repository are defined as follows:

_FlowJo_SLICE.Rmd_: This file contains code to convert FlowJo gating data from all experiment in the PDTS/SLICE project (exported as an Excel file) to a data table containing numbers and percentages of each gated cell population per sample. Furthermore, this file contains scripts for the visualisation and export of this flow cytometry data.

_FlowJo_TLS.Rmd_: This file contains code to convert FlowJo gating data in the PDTF/TLS project (exported as an Excel file) to a data table containing numbers and percentages of each gated cell population per sample. Furthermore, this file contains scripts for the visualisation and export of this flow cytometry data. This script also performs filtering of data based on cell numbers, as well as assigning fragments B cell high/intermediate/low scores. Any analyses involving cellular populations in the TLS project are in this file.

_Legendplex.Rmd_: This is a general-use Rmd file which can be used to calculate several metrics from Legendplex data: it takes a file containing raw cytokine and chemokine values and calculates corrected values, deltas and PDTF scores and outputs this per tumour individually.

_Legendplex_from_template.Rmd_: This is an Rmd file to handle Legendplex data from the TLS project specifically. It can perform all functions of the Legendplex.Rmd file, but in addition can also perform normalisation and quality control based on concentrations and MFI.

_legendplex_all_tumours.R_: This file is used to handle all output form the Legendplex_from_template.Rmd script. It is used to calculate z-scores and cumulative scores for many different comparisons, fold-changes, outlier correction and heatmap plotting. Any analysis involving Legendplex data in the TLS project is handled in this file.

_Functions.R_: A master file containing several, more complex, functions that were designed to aid in the analysis of all data. Other scripts, such as _legendplex_all_tumours.R_, rely on some of the functions specified in this file and will call these in using the source() function from R.


If there are any questions regarding the scripts included in this project, place them in the Discussion page for this GitHub repository.
Happy coding!

Sidney
