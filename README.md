# FMD-brachialtools-analyzer
R script for the analysis of FMD including viscosity and shear stress calculations
FMD_R ANALYSIS README

It is imperative that files are set up properly before use with FMD_R 

If running for the first time, COPY .r file into directory you want the files to populate 
- This folder will hold all data for the entire project/analysis
- Run the file by opening it in RStudio. select all the text from the r file (should load in top left workspace) and click Run

You should now have a folder labelled "data" nested in the folder alongside your .r file
Into this folder, you can now drop your working data

Place the following files into the "data" folder

File 1: Baseline diameter data report exported from brachial tools (.txt file)
-Settings for export:
	Must include Reader ID, Participant ID, and StudyType as condition, export report as .csv
- rename BL_DIA.txt
File 2: Post-occlusion data report exported from brachial tools (.txt file)
-Settings for export:
	Must include Reader ID, Participant ID, and StudyType as condition, export report as .csv
- rename FMD_DIA.txt
File 3: Baseline labchart data exported from labchart - 
- Settings for export:
	channels:QDAT flow, Finger Pressure <- Channel order must match
	export selection(highlight baseline data)
	downsample by 333 for 10ksample rate (33 for 1k sample rate)
	blockheaderON,
	Time selected as true (always seconds)
	everything else unchecked
- rename BL_LC.txt
File 4: Post-occlusion labchart data exported from labchart - 
- Settings for export:
	channels:QDAT flow, Finger Pressure <- Channel order must match
	export selection(highlight post-occlusion data)
	downsample by 333 for 10ksample rate (33 for 1k sample rate)
	blockheaderON,
	Time selected as true
	everything else unchecked
- rename BL_LC.txt
File 5: Viscometer report (.csv)
- rename VISC.csv
- if no visco report omit this file and run the NO_VISC version of the script.

Run R script
Results immediately available for QC in working directory and saved under participant IDs
