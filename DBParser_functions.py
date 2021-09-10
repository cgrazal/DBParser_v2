import shutil
from shutil import copyfile
import sys
from sys import exit 
import os
import sqlite3
from sqlite3 import Error
import pandas as pd 
import csv
import numpy as np
#from pathlib import Path
from datetime import date
random_int = np.random.randint(1,9999,1)

#Create mysqlite3 DB
def create_connection(sqliteDB):
    conn=None
    try:
        conn=sqlite3.connect(sqliteDB)
        return conn
    except Error as e:
        print(e)
    
    return conn

def create_table(conn, create_table_sql):
	try:
		c=conn.cursor()
		c.execute(create_table_sql)
	except Error as e:
		print (e)

def main(DBname):
	database = DBname

	sql_create_Analyzed_data="""CREATE TABLE IF NOT EXISTS ANALYZED_DATA_TABLE(
	FILE_ID INTEGER,
	ROW_ID INTEGER,
	ATTR_ID TEXT,
	ATTR_VAL INTEGER
	);"""

	sql_create_Sample_variants="""CREATE TABLE IF NOT EXISTS SAMPLE_VARIANTS_TABLE(
	FILE_ID INTEGER,
	ROW_ID INTEGER,
	ATTR_ID TEXT,
	ATTR_VAL INTEGER
	);"""

	sql_create_Metadata="""CREATE TABLE IF NOT EXISTS SAMPLE_METADATA_TABLE(
	FILE_ID INTEGER,
	ROW_ID INTEGER,
	ATTR_ID TEXT,
	ATTR_VAL INTEGER
	);"""

	conn=create_connection(database)

	if conn is not None:
		create_table(conn, sql_create_Analyzed_data)
		create_table(conn, sql_create_Sample_variants)
		create_table(conn, sql_create_Metadata)
		print ("\nConnected to database!\n")
        
	else:
		print ("Error- Cannot create the database connection.")

def GetVSALIGN(mets_paths, dbname, CGSNumber):
	vsalign_results = []

	for each in mets_paths:
		Cleaned = each + "/01_Cleaning.met"
		Aligned = each + "/02_Aligned.met"
		SNPTable = each + "/03_SNP_Table.met"

		with open(Cleaned) as csv_file:
			cl_data = pd.read_csv(csv_file, sep='\t')

		with open(Aligned) as csv_file:
			al_data = pd.read_csv(csv_file,  sep='\t')

		with open(SNPTable) as csv_file:
			snp_data = pd.read_csv(csv_file, sep='\t')


		cl_data.sort_values(by=['SAMPLE'], inplace=True)
		al_data.sort_values(by=['SAMPLE'], inplace=True)
		snp_data.sort_values(by=['SAMPLE'], inplace=True)

    
		#take average of paired reads in Aligned file
		al_data2 = al_data.groupby('SAMPLE', as_index=False).mean()

		cl_data.rename(columns={'SAMPLE':'SAMPLE_ID'}, inplace=True)
		al_data2.rename(columns={'SAMPLE':'SAMPLE_ID'}, inplace=True)
		snp_data.rename(columns={'SAMPLE':'SAMPLE_ID'}, inplace=True)

		#combine dfs
		merged=pd.merge(left=cl_data, right=al_data2, how='outer', left_on="SAMPLE_ID", right_on="SAMPLE_ID")
		Results_int= pd.merge(left=merged, right=snp_data, how='outer', left_on="SAMPLE_ID", right_on="SAMPLE_ID")

		Results_int.sort_values(by=['SAMPLE_ID'], inplace=True)

		Results_int['CGS_Number']=CGSNumber

		vsalign_results.append(Results_int)

		vsalign_results_all = pd.concat(vsalign_results)
	
	return(vsalign_results_all)

def GetPANGO(pango):
	pango_results = []
	for each in pango:
		with open (each) as csv_file:
			pangolin_lineages = pd.read_csv(csv_file)
			print(pangolin_lineages.columns)
		#pangolin_lineages.rename(columns={'?Sequence name':'SAMPLE_ID'}, inplace=True)
		pango_results.append(pangolin_lineages)
		pangolin_all = pd.concat(pango_results)
	return(pangolin_all)

def GetNEXTCLADE(nc):
	nc_results = []
	for each in nc:
		with open (each) as csv_file:
			nc_lineages = pd.read_csv(csv_file, sep=';')
		#nc_lineages.rename(columns={'seqName': 'SAMPLE_ID'}, inplace=True)
		nc_results.append(nc_lineages)
		nc_all = pd.concat(nc_results)
	return(nc_all)

def AnalysisResults (mets_paths, dbname, CGSNumber, pango, nc):
	vsalign = GetVSALIGN(mets_paths, dbname, CGSNumber)
	pangolin = GetPANGO(pango)
	nextclade = GetNEXTCLADE(nc)

	lineages = pd.merge(left = pangolin, right = nextclade, how = 'outer', left_on = 'Sequence name', right_on = 'seqName')
	lineage_cols_to_include = ['SAMPLE_ID', 'clade', 'qc.overallScore', 'totalSubstitutions', 'totalAminoacidSubstitutions', 
	'substitutions', 'deletions', 'insertions', 'aaSubstitutions', 'aaDeletions', 'qc.privateMutations.total',
	'Lineage', 'Ambiguity', 'Scorpio call', 'Scorpio support', 'pangolin version', 'pangoLEARN version']

	lineages = lineages.dropna(axis=0, subset=['SAMPLE_ID'])
	lineages.drop(lineages.index[lineages['SAMPLE_ID'] == 'NC_045512.2'], inplace=True)
	lineages.drop(lineages.index[lineages['SAMPLE_ID'] == 'NC_045512.2_v4'], inplace=True)
	lineages_revised = lineages.drop(columns = [col for col in lineages if col not in lineage_cols_to_include])

	analyzed_data = pd.merge(left = vsalign, right = lineages_revised, how = 'outer', left_on = 'SAMPLE_ID', right_on = 'SAMPLE_ID')
	Samples = analyzed_data['SAMPLE_ID'].apply(str)

	spl1 = []
	for each in Samples:
		value = each.split("_", 1)[1]
		spl1.append(value)
	
	samples = []
	for each in spl1:
		value = each.rsplit("_",1)[0]
		samples.append(value)

	analyzed_data['SAMPLE_ID'] = samples

	analyzed_data['ROW_ID'] = analyzed_data.index + 1

	cols = analyzed_data.columns.tolist()
	cols.pop(-1)
	
	Analyzed_Results = pd.melt(analyzed_data, id_vars=["ROW_ID"],value_vars=cols, var_name='ATTR_ID',value_name="ATTR_VAL")

	Analyzed_Results.loc[:, 'FILE_ID'] = random_int
	Analyzed_Results = Analyzed_Results[['FILE_ID', 'ROW_ID', 'ATTR_ID', 'ATTR_VAL']]

	#Insert into database
	databasename = dbname
	conn=sqlite3.connect(databasename)
	cur=conn.cursor()
	Analyzed_Results.to_sql("ANALYZED_DATA_TABLE", conn, if_exists='append', index=False)
	conn.commit()

	print("\nANALYZED_DATA input complete. [1/3]\n")

def MetadataTable(meta, assign, dbname):
	with open(meta) as csv_file:
		metadata = pd.read_csv(csv_file)

	with open(assign) as csv_file:
		assigned_data = pd.read_csv(csv_file)

	final_meta = pd.merge(left = metadata, right = assigned_data, how = "outer", left_on = "SAMPLE_ID", right_on = "SAMPLE_ID")

	final_meta = final_meta.dropna(axis=0, subset=['SAMPLE_ID'])

	today = date.today()
	final_meta["TIMESTAMP"] = today.strftime("%Y%m%d")

	final_meta['ROW_ID'] = final_meta.index + 1

	cols = final_meta.columns.tolist()
	cols.pop(-1)
	
	Metadata_results = pd.melt(final_meta, id_vars=["ROW_ID"],value_vars=cols, var_name='ATTR_ID',value_name="ATTR_VAL")

	Metadata_results.loc[:, 'FILE_ID'] = random_int
	Metadata_results = Metadata_results[['FILE_ID', 'ROW_ID', 'ATTR_ID', 'ATTR_VAL']]
		
	#Insert into database
	databasename = dbname
	conn=sqlite3.connect(databasename)
	conn.text_factory = str
	cur=conn.cursor()
	Metadata_results.to_sql("SAMPLE_METADATA_TABLE", conn, if_exists='append', index=False)
	conn.commit()

	print("\nMETADATA input complete. [2/3]\n")

def GetVariants(SV_paths, dbname):
	fullList = []
	for each in SV_paths:
		files = []
		files += [name for name in os.listdir(each) if name.endswith('.snp')]
		fullList.append(files)

	fullPathsnp = []
	for i in range(len(SV_paths)):
		for each in fullList[i]:
			paths = SV_paths[i] + "/" + each
			fullPathsnp.append(paths)

	for item in fullPathsnp:
		if item.endswith("snpSummary.snp"):
			fullPathsnp.remove(item)

	snpData = []
	for each in fullPathsnp:
		with open(each) as csv_file:
			data = pd.read_csv(csv_file, sep='\t', low_memory=False)
		
		name_path = os.path.basename(each)
		name = os.path.splitext(name_path)[0]
		data["SAMPLE_ID"] = name
		snpData.append(data)
		variants = pd.concat(snpData)
	
	return(variants)

def VariantsResults(SV_paths, dbname):
	allvariants = GetVariants(SV_paths, dbname)

	colList = ["SAMPLE_ID","TYPE", "BASE", "SNP","FREQ", "CODON", "FEATURE", "DEPTH","DEL", "INS", "SCD", "NOTE"]
	final_variants = allvariants.drop (columns = [col for col in allvariants if col not in colList])
	final_variants.sort_values(by = ['SAMPLE_ID'], inplace = True)
	final_variants['ROW_ID'] = (final_variants['SAMPLE_ID']!= final_variants['SAMPLE_ID'].shift()).cumsum()

	final_variants = final_variants[["SAMPLE_ID","TYPE", "BASE", "SNP","FREQ", "CODON", "FEATURE", "DEPTH","DEL", "INS", "SCD", "NOTE", "ROW_ID"]]
	cols = final_variants.columns.tolist()
	cols.pop(-1)

	variants_results = pd.melt(final_variants, id_vars=["ROW_ID"],value_vars=cols, var_name='ATTR_ID',value_name="ATTR_VAL")
			
	variants_results.loc[:, 'FILE_ID'] = random_int
	variants_results = variants_results[['FILE_ID', 'ROW_ID', 'ATTR_ID', 'ATTR_VAL']]

	databasename = dbname
	conn = sqlite3.connect(databasename)
	cur = conn.cursor()
	variants_results.to_sql ("SAMPLE_VARIANTS_TABLE", conn, if_exists='append', index=False)
	conn.commit()

	print("\nSAMPLE_VARIANTS_TABLE input complete. [3/3]\n")

	print("\nDATABASE UPDATED.")

#View updated


