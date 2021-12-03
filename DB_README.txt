README FILE: Database_System
All coding: CGrazal, 2021
Format: RKumar

LOCATION:
Prod3:
/data/cgrazal/Database_System

##################################################################################################################
A. LAYOUT/DESCRIPTION: 
	- Within the Database_System, there are 4 main folders:
		-DATABASES
		-REPORTS
		-SCRIPTS
		-DATA
B. SCRIPT HOW TO RUN

####################################################################################################################
A. LAYOUT/DESCRIPTION
- DATABASES: 
	- Where the databases are stored. At present, there are three (GEIS_DB_forReport.db, 
	GEIS_FinalDB_storage, and OtherProjects_DB).
	
	- For COVID/GEIS project, the "forReport" db is used to generate an automated preliminary report. Consider it a "single use".
	- The "FinalDB_storage" db is used as a storage component for all samples/runs. Consider it appended of all Geis CGS samples/running db.
	- These database pull from Joshua _Assigned and _Metadata csv sheets with sample information,
		and also includes all VSalign information, and full path to fasta files on /mnt/CGS_2021/Covid_Data_Analysis.

	- OtherProjects_DB is used for any other project. It includes only VSALIGN output.
		-Metadata table is the "override.txt" file from VSALIGN.

- REPORTS:
	- This folder is where generated reports are saved to.

- SCRIPTS:
	- Three versions are stored currently, corresponding to GEIS_forReport, GEIS_storage,
	and OtherProjects.
	- For GEIS, use the "DBParser_forReport" to create a single use db used to generate the automated report.
	- For GEIS Storage to append to ongoing storage db, use "DBParser_finalDB_storage" script.
	- DBParser_OtherProjects.py is used for all other projects.
	- There is a CREATE_REPORT folder with script to create automated preliminary report.

- DATA:
	- The data used for this is copied from P_drive/VGB to the DATA folder.
	- For Otherprojects, the _Assigned and _Metadata folders are irrelevant; it only pulls from the CGS folder.
	- For Covid/GEIS, the _Assigned and _Metadata folders are relevant.
		- "GEIS_Assigned_samples.xlsx" and "GEIS_Metadata.xlsx" are copied to here from VGB.

#####################################################################################################################
B. HOW TO RUN

1. Copy CGS_000### file from VGB to /data/cgrazal/Database_System/DATA, making sure to follow correct file format.
	- Be mindful of the organization of the CGS folder. 
		- It should be CGS_### > analyzed_data > ###k_yyyymmdd > Sample folder names > VSalign output for each
	- Check previous CGS numbers to make sure a new run has the same layout, as the parser depends on that layout.

2. If Covid/GEIS, copy "GEIS_Assigned_Samples.xlsx" and "GEIS_Metadata.xlsx" from VGB to /data/cgrazal/Database_System/DATA/_Assigned or _Metadata.
	- Edit files: Open and copy only working data to a new file, save as .csv
	- For GEIS_Assigned_samples, remove the top header noting, "Trizol", "Extraction", etc.
	- Working .csv sheet should have just the names "GEIS_Assigned.Samples.csv" or "GEIS_Metadata.csv"
	- Previous versions are kept here with the date run in the name.
	- Delete the .xlsx version, not needed.

3. Navigate to SCRIPTS folder.
	--> RUN:  python DBParser_...version..._.py /data/cgrazal/Database_System/DATA/CGS_000###
	- Make sure no "/" at end of CGS_000### 

4. Once steps 1-3 are done and database is updated, create report.
	- Version is "GEIS_create_report.sh" or "OtherProjects_create_report.sh"
	- Edit the /SCRIPTS/CREATE_REPORT/file.sh
		1. Under "CREATE VIEW REPORT1 AS...", WHERE CGS_NUMBER = "#####" ; UPDATE this with CGS_000### of current run
		2. At the end of file, change the ".output ../REPORTS/######" to the name of the Preliminary report you want to create.
			- Good template: PRELIMINARY_REPORT_CGS_000###_YYYYMMDD.csv"	
	--> TO RUN: Navigate to DATABASES folder; Copy and paste entire "create_report.sh" file into terminal.

5. Final report is automatically saved under REPORTS.

6. NOTES: If the report is blank- Know that for Covid/GEIS samples, it depends on the samples being listed on GEIS_Metadata.csv or GEIS_Assigned_Samples.csv,
with the correct CGS Number. 9/10 times check those data sheets and it is not listed. Make sure the samples are added on .csv sheets first, then try re-running.


DONE! 
THANKS!
