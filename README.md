# FMD-brachialtools-analyzer
R script for the analysis of FMD including viscosity and shear stress calculations

The FMD-R Scripts are designed for the automated analysis of FMD and blood viscosity using data from Brachial Tools Software and Labchart

Requires:
 R ([link](https://mirror.rcg.sfu.ca/mirror/CRAN/)) and RStudio ([link](https://posit.co/downloads/)) 
## Creating your project folder


Into this folder, you can now drop your working data

## File Preparation

All files should be in the same folder, and a popup dialog will appear under your R program to select the directory where these are located
 1. File 1: Baseline diameter data report exported from brachial tools (.txt file)
- Settings for export:
	Must include Reader ID, Participant ID, and StudyType as condition, export report as .csv
- rename BL_DIA.txt
2. File 2: Post-occlusion data report exported from brachial tools (.txt file)
- Settings for export:
	Must include Reader ID, Participant ID, and StudyType as condition, export report as .csv
- rename FMD_DIA.txt

3. File 3: Baseline labchart data exported from labchart - 
- Settings for export:
	- channels: Flow velocity, Finger pressure <- Channel order must match
	- export selection(highlight baseline data)
	- downsample by 333 for 10ksample rate (33 for 1k sample rate)
	- blockheaderON,
	- Time selected as true (always seconds)
	- everything else unchecked
- rename BL_LC.txt
4. File 4: Post-occlusion labchart data exported from labchart - 
- Settings for export:
	- See above
- rename FMD_LC.txt
5. File 5: Viscometer report (.csv)
- rename VISC.csv
- if no visco report omit this file and run the no_visc version of the script.
## Analyzing Files
Run FMD launcher and select version required using command
 ```R
 source("https://github.com/nkcheung95/FMD-brachialtools-analyzer/blob/main/FMDBTOOLSAnalysis-loader.R?raw=TRUE")
 ```

 Select Output folder using button, and run desired program version




