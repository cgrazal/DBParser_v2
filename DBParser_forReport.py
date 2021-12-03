#Parser to create report from original files.
#!/usr/bin/python
import functions_forReport as DB
import sqlite3
import sys
import os
from os import walk
from os.path import dirname
import numpy as np
from datetime import datetime

###################### GET NECESSARY INFO ######################
#Get database
database_input = "/data/scripts/cgrazal/Database_System_backup/DATABASES/GEIS_forReport.db"

#Data source is CGS# from _GEIS_2021 folder
data_source=str(sys.argv[-1])
print(data_source)

cgs_number = data_source.rsplit('/',1)[-1]
random_int = np.random.randint(1,9999,1)

file_report_name = "PRELIMINARY_REPORT_" + datetime.now().strftime("%Y%m%d") + "_" + cgs_number + ".csv"
#print(file_report_name)

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
		print(p_name)
		print(each)
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

########Get .fasta full path from /mnt/Covid_Data_Analysis############
fasta_path = "/mnt/CGS_2021/Covid19_Data_Analysis/" + cgs_number + "/Analyzed_data"		
inputAmt = walk(fasta_path).next()[1]

stndInput = []
for each in inputAmt:
	folders = fasta_path + "/" + each
	stndInput.append(folders)

standard_folders = []
for each in stndInput:
	sampfolders = walk(each).next()[1]
	standard_folders.append(sampfolders)

out_puts = [[] for i in range(len(standard_folders))]
for i in range(len(standard_folders)):	
	for each in standard_folders[i]:
		info = stndInput[i]+"/"+(each)
		out_puts[i].append(info)

for i in out_puts:
	if i == []:	
		out_puts.remove(i)

genomes = [[] for i in range(len(out_puts))]
for i in range(len(out_puts)):
	for each in out_puts[i]:
		genome_folder = each + "/" + "_genomes"
		genomes[i].append(genome_folder)

genome_list = [[] for i in range(len(genomes))]
for i in range(len(genomes)):
	for each in genomes[i]:
		gen_fasta = [file for file in os.listdir (each) if file.startswith("S") and file.endswith(".fasta")]
		genome_list[i].append(gen_fasta)

gen_full_paths=[]
genomes_working = genomes[0]
genome_list_working = genome_list[0]
for i in range(len(genomes_working)):
	for each in genome_list_working[i]:
		g_path = genomes_working[i] + "/" + each
		gen_full_paths.append(g_path)

######################INPUT DATA#############################

############Row IDs############
DB.GetRowID_Analysis(database_input, random_int)
DB.GetRowID_snp(database_input, random_int)

############Analyzed_data############
for (a, b, c) in zip (metrics, pangolin, nextclade):
	DB.GetVSALIGN(a, database_input, cgs_number)
	DB.GetPANGO(b)
	DB.GetNEXTCLADE(c)
	DB.AnalysisResults(a, database_input, cgs_number, b, c, random_int, gen_full_paths)
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
