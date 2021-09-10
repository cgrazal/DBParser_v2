#Parser to create report from original files.
#!/usr/bin/python
import DBParser_functions as DB
import sqlite3
import sys
import os
from os import walk
from os.path import dirname

#Data source is CGS# from _GEIS_2021 folder
data_source=str(sys.argv[-1])
print(data_source)
cgs_number = data_source.rsplit('/',1)[-1]

#/Users/Mine./Desktop/USAMRIID/DummyData

#Get GEIS sample sheets from _GEIS_2021
directory = dirname(data_source)
metadata = directory + "/" + "_Metadata/GEIS_Metadata.csv" 
assigned = directory + "/" + "_Assigned/GEIS_Assigned_samples.csv"

#Get snp and metrics
analyzed_data = data_source+"/analyzed_data"
vsalign_input = walk(analyzed_data).next()[1]

vsalign_input_list_path = []
for each in vsalign_input:
	foldernames = analyzed_data+"/"+each
	vsalign_input_list_path.append(foldernames)

vsalign_standard_folders = []
for each in vsalign_input_list_path:
	samplefolders = walk(each).next()[1]
	vsalign_standard_folders.append(samplefolders)

outputs = []
for i in range(len(vsalign_input_list_path)):
	for each in vsalign_standard_folders[i]:
		info = vsalign_input_list_path[i]+"/"+(each)
		outputs.append(info)

snps = []
metrics = []

for each in outputs:
	snpFolder = each + "/" + "_snps"
	metricsFolder = each + "/" + "_metrics"
	snps.append(snpFolder)
	metrics.append(metricsFolder)

pango_list = []
nextclade_list = []
for each in outputs:
	p_name = [file for file in os.listdir (each) if file.startswith("results")]
	n_name = [file for file in os.listdir(each) if file.startswith("nextclade")]
	pango_list.append(p_name)
	nextclade_list.append(n_name)

pangolin=[]
for i in range(len(outputs)):
	for each in pango_list[i]:
		pango_paths = outputs[i] + "/" + each
		pangolin.append(pango_paths)

nextclade=[]
for i in range(len(outputs)):
	for each in nextclade_list[i]:
		nextclade_paths = outputs[i] + "/" + each
		nextclade.append(nextclade_paths)
#Get database
database_input = "WGS_Database_v2.db"
print(nextclade)
#INPUT DATA
#Analyzed_data
DB.GetVSALIGN(metrics, database_input, cgs_number)
DB.GetPANGO(pangolin)
DB.GetNEXTCLADE(nextclade)
DB.AnalysisResults(metrics, database_input, cgs_number, pangolin, nextclade)

#GEIS Metadata
DB.MetadataTable(metadata, assigned, database_input)

#Sample_variants
DB.GetVariants(snps, database_input)
DB.VariantsResults(snps, database_input)

