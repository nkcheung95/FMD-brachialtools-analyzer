# FMD-brachialtools-analyzer
R script for the analysis of FMD including viscosity and shear stress calculations

The FMD-R Scripts are designed for the automated analysis of FMD and blood viscosity using data from Brachial Tools Software and Labchart

Requires:
 R ([link](https://mirror.rcg.sfu.ca/mirror/CRAN/)) and RStudio ([link](https://posit.co/downloads/)) 
## Creating your project folder

Begin setting up for analysis by creating a project in RStudio (File > New Directory > New Project).
In the new project console, copy and paste the following code:
```R
source("https://github.com/nkcheung95/FMD-brachialtools-analyzer/blob/main/FMD_R_mod_visc.r?raw=TRUE")
```

You should now have a folder labelled "data" nested in the working folder.
Into this folder, you can now drop your working data

## File Preparation

Place the following files into the "data" folder

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
Run Script version required using given command

NEW* Launcher here
```R
source("https://github.com/nkcheung95/FMD-brachialtools-analyzer/blob/main/FMDBTOOLSAnalysis-loader.R?raw=TRUE")
```

 - FMD-R mod_visc
	 - For application of a 2-phase exponential decay for shear-dependent viscosity modelling
```R
source("https://github.com/nkcheung95/FMD-brachialtools-analyzer/blob/main/FMD_R_mod_visc.r?raw=TRUE")
```
 - FMD-R sp_visc
	 - Single-point viscosity at a shear rate of 225
```R
source("https://github.com/nkcheung95/FMD-brachialtools-analyzer/blob/main/FMD_R_sp_visc.r?raw=TRUE")
```
 - FMD-R no_visc
	 - For running FMD's without viscosity
```R
source("https://github.com/nkcheung95/FMD-brachialtools-analyzer/blob/main/FMD_R_no_visc.r?raw=TRUE")
```
- FMD-R mod_visc_FLIPPED
	- FMD-R mod_visc with the labchart channels in opposite order (finger pressure, flow velocity)
```R
source("https://github.com/nkcheung95/FMD-brachialtools-analyzer/blob/main/FMD_R_mod_visc_FLIPPED.r?raw=TRUE")
```
- FMD-R no_visc_FLIPPED
	- FMD-R no_visc with the labchart channels in opposite order (finger pressure, flow velocity)
```R
source("https://github.com/nkcheung95/FMD-brachialtools-analyzer/blob/main/FMD_R_no_visc_FLIPPED.r?raw=TRUE")
```
Results immediately available for QC in working directory and saved under participant IDs from the brachial tools output.


