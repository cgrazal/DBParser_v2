#Parser to create report from original files.
#!/usr/bin/python
import DBParser_functions_appended as DB
import sqlite3
import sys
import os
from os import walk
from os.path import dirname
import numpy as np

###################### GET NECESSARY INFO ######################
#Get database
database_input = "/mnt/CGS_2021/users/cgrazal/Database_v2/WGS_Database_appended.db"

#Data source is CGS# from _GEIS_2021 folder
data_source=str(sys.argv[-1])
print(data_source)
cgs_number = data_source.rsplit('/',1)[-1]
random_int = np.random.randint(1,9999,1)

#Get GEIS sample sheets from _GEIS_2021
directory = dirname(data_source)
metadata = directory + "/" + "_Metadata/GEIS_Metadata.csv" 
assigned = directory + "/" + "_Assigned/GEIS_Assigned_samples.csv"

#Get snp and metrics for each sample grouping.
analyzed_data = data_source+"/analyzed_data"
vsalign_input = walk(analyzed_data).next()[1]

vsalign_input_list_path = []
for each in vsalign_input:
	foldernames = analyzed_data + "/" + each
	vsalign_input_list_path.append(foldernames)

vsalign_standard_folders = []
for each in vsalign_input_list_path:
	samplefolders = walk(each).next()[1]
	vsalign_standard_folders.append(samplefolders)

outputs = [[] for i in range(len(vsalign_input_list_path))]
for i in range(len(vsalign_input_list_path)):	
	for each in vsalign_standard_folders[i]:
		info = vsalign_input_list_path[i]+"/"+(each)
		outputs[i].append(info)

snps = [[] for i in range(len(outputs))]
metrics = [[] for i in range(len(outputs))]
for i in range(len(outputs)):
	for each in outputs[i]:
		snpFolder = each + "/" + "_snps"
		metricsFolder = each + "/" + "_metrics"
		snps[i].append(snpFolder)
		metrics[i].append(metricsFolder)

#Get the pango and nextclade outputs for each grouping.
pango_list = [[] for i in range(len(outputs))]
nextclade_list = [[] for i in range(len(outputs))]
for i in range(len(outputs)):
	for each in outputs[i]:
		p_name = [file for file in os.listdir (each) if file.startswith("results")]
		n_name = [file for file in os.listdir(each) if file.startswith("nextclade")]
		p_name = str(p_name[0])
		n_name = str(n_name[0])
		pango_list[i].append(p_name)
		nextclade_list[i].append(n_name)
pangolin=[]
nextclade=[]

for (a, b) in zip (outputs, pango_list):
	pango_paths = [e + "/" + f for e, f in zip(a, b)]
	pangolin.append(pango_paths)
for (c, d) in zip (outputs, nextclade_list):
	nc_paths = [g + "/" + h for g, h in zip(c, d)]
	nextclade.append(nc_paths)


######################INPUT DATA#############################

############Row IDs############
DB.GetRowID_Analysis(database_input, random_int)
DB.GetRowID_snp(database_input, random_int)

############Analyzed_data############
for (a, b, c) in zip (metrics, pangolin, nextclade):
	DB.GetVSALIGN(a, database_input, cgs_number)
	DB.GetPANGO(b)
	DB.GetNEXTCLADE(c)
	DB.AnalysisResults(a, database_input, cgs_number, b, c, random_int)

print("\nANALYZED_DATA input complete. [1/3]\n")

############GEIS Metadata############
DB.MetadataTable(metadata, assigned, database_input, random_int)

print("\nMETADATA input complete. [2/3]\n")

############Sample_variants############
for item in snps:
	DB.GetVariants(item, database_input)
	DB.VariantsResults(item, database_input, random_int)

print("\nSAMPLE_VARIANTS_TABLE input complete. [3/3]\n")
print("\nDATABASE UPDATED.")
