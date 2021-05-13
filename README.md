# POAR-range-limits

This repo contains all the data and code necessary to reproduce the analyses in: Miller and Compagnoni, "Two-sex demography, sexual niche differentiation, and the geographic range limits of Texas bluegrass".

## In the code folder:
* The script https://github.com/texmiller/POAR-range-limits/blob/master/code/POAR_range_limits_analysis.R includes all statistical and demographic modeling and produces all of the data figures in the manuscript. This script downloads data from files stored in the cloud. Some sections of this script takes hours-to-days to run; these are commented out and we link straight to the saved objects that they produce. You could always recreate those objects by uncommenting and running that code. 

* The script https://github.com/texmiller/POAR-range-limits/blob/master/code/twosexMPM.R includes all of the source functions of the two-sex MPM. This script is sourced by POAR_range_limits_analysis.R.

* The script https://github.com/texmiller/POAR-range-limits/blob/master/code/POAR_size_distribution_appendix.R includes the analysis presented in Appendix C of the manuscript

* There are several sub-folders under "code". The "format" folder includes all of our data manipulation steps (these need not be re-run to reproduce our analysis -- POAR_range_limits_analysis.R works with derived data products produced in the format scripts).

* The "stan" folder includes stan model scripts. POAR_range_limits_analysis.R calls one of these scripts.

## In the Manuscript folder:

* https://github.com/texmiller/POAR-range-limits/blob/master/Manuscript/POAR_range_limits_AmNat_submission1.Rnw is the manuscript file. The manuscript was written in knitr and it sources figures and data objects generated by POAR_range_limits_analysis.R. Compiling this file will generate a .tex file of the same name (in this folder), and compiling that one will generate a .pdf file of the same name (also in this folder). 

## Any questions, feel free to contact me at tom.miller@rice.edu
