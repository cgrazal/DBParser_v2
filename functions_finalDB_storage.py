import shutil
from shutil import copyfile
import sys
from sys import exit 
import os
import sqlite3
from sqlite3 import Error
import pandas as pd 
import csv
#from pathlib import Path
from datetime import datetime

######### Create mysqlite3 DB #########
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
	conn.text_factory = str

	if conn is not None:
		create_table(conn, sql_create_Analyzed_data)
		create_table(conn, sql_create_Sample_variants)
		create_table(conn, sql_create_Metadata)
		print ("\nConnected to database!\n")
        
	else:
		print ("Error- Cannot create the database connection.")

######### Input incrememnting row Id's; New file start with last row ID and add +1. #########
def GetRowID_Analysis(dbname, r_int):
	databasename = dbname
	conn=sqlite3.connect(databasename)
	cur=conn.cursor()

	tablequery = cur.execute("""SELECT name from sqlite_master where type = 'table' and name = 'ANALYZED_DATA_TABLE';""").fetchall()
	rInt = ' '.join([str(elem) for elem in r_int])
	if tablequery == []:
		maximum_analysis = 1
	else:
		rowquery = cur.execute("""SELECT MAX(ROW_ID) as maximum from ANALYZED_DATA_TABLE where FILE_ID =?;""", (rInt,)).fetchall()
		
		if rowquery == [(None,)]:
			maximum_analysis = 1
		else:
			cur.execute("SELECT MAX(ROW_ID) as maximum from ANALYZED_DATA_TABLE where FILE_ID = ?;", (rInt,))
			result = cur.fetchall()
			for i in result:
				maximum_analysis = float(i[0]) + 1
	return (maximum_analysis)

def GetRowID_snp(dbname, r_int):
	databasename = dbname
	conn=sqlite3.connect(databasename)
	cur=conn.cursor()

	tablequery = cur.execute("""SELECT name from sqlite_master where type = 'table' and name = 'SAMPLE_VARIANTS_TABLE';""").fetchall()
	rInt = ' '.join([str(elem) for elem in r_int])
	if tablequery == []:
		maximum_snp = 1
	else:
		rowquery = cur.execute("""SELECT MAX(ROW_ID) as maximum from SAMPLE_VARIANTS_TABLE where FILE_ID =?;""", (rInt,)).fetchall()
		
		if rowquery == [(None,)]:
			maximum_snp = 1
		else:
			cur.execute("SELECT MAX(ROW_ID) as maximum from SAMPLE_VARIANTS_TABLE where FILE_ID = ?;", (rInt,))
			result = cur.fetchall()
			for i in result:
				maximum_snp = float(i[0]) + 1
	return (maximum_snp)

##########################################################################################
############################## Analyzed Data Table ##############################

######### Get VSAlign output #########

def GetVSALIGN(mets_paths, dbname, CGSNumber):
	vsalign_results = []

	for each in mets_paths:
		Cleaned = each + "/01_Cleaning.met"
		Aligned = each + "/02_Aligned.met"
		SNPTable = each + "/03_SNP_Table.met"

		with open(Cleaned) as csv_file:
			cl_data = pd.read_csv(csv_file, sep='\t', low_memory=False)

		with open(Aligned) as csv_file:
			al_data = pd.read_csv(csv_file,  sep='\t', low_memory=False)

		with open(SNPTable) as csv_file:
			snp_data = pd.read_csv(csv_file, sep='\t', low_memory=False)


		cl_data.sort_values(by=['SAMPLE'], inplace=True)
		cl_data.drop(cl_data.index[cl_data['SAMPLE'] == 'Sample'], inplace=True)
		cl_data.drop(cl_data.index[cl_data['SAMPLE'] == 'Sequence'], inplace=True)
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

		gbPath = os.path.split(each)[0] + "/"
		items = os.listdir(gbPath)
		gbFile = []
		for item in items:
			if item.endswith(".gb"):
				gbFile.append(item)
		refFile = ''.join(gbFile)
		Results_int['GENBANK_REFERENCE_FILE'] = refFile


		append = vsalign_results.append
		append(Results_int)
		vsalign_results_all = pd.concat(vsalign_results)
	
	return(vsalign_results_all)

######### Get Pangolin and Nextclade #########

def GetPANGO(pango):
	pango_results = []
	for each in pango:
		with open (each) as csv_file:
			pangolin_lineages = pd.read_csv(csv_file, encoding = 'utf-8-sig', low_memory=False)
		pangolin_lineages.rename(columns={'Sequence name':'SAMPLE_ID'}, inplace=True)
		pango_results.append(pangolin_lineages)
		pangolin_all = pd.concat(pango_results)

	return(pangolin_all)

def GetNEXTCLADE(nc):
	nc_results = []
	for each in nc:
		with open (each) as csv_file2:
			nc_lineages = pd.read_csv(csv_file2, sep='\t', encoding = 'utf-8-sig', low_memory=False)
		nc_lineages.rename(columns={'seqName': 'SAMPLE_ID'}, inplace=True)
		nc_results.append(nc_lineages)
		nc_all = pd.concat(nc_results)

	return(nc_all)

######### INPUT Analysis output (VSalign and Pangolin/Nextclade) #########

def AnalysisResults (mets_paths, dbname, CGSNumber, pango, nc, r_int):
	vsalign = GetVSALIGN(mets_paths, dbname, CGSNumber)
	pangolin = GetPANGO(pango)
	nextclade = GetNEXTCLADE(nc)

	lineages = pd.merge(left = pangolin, right = nextclade, how = 'outer', left_on = 'SAMPLE_ID', right_on = 'SAMPLE_ID')
	lineages = lineages.dropna(axis=0, subset=['SAMPLE_ID'])
	lineages.drop(lineages.index[lineages['SAMPLE_ID'] == 'NC_045512.2'], inplace=True)
	lineages.drop(lineages.index[lineages['SAMPLE_ID'] == 'NC_045512.2_v4'], inplace=True)

	included_cols = ['SAMPLE_ID', 'clade', 'qc.overallScore', 'totalSubstitutions', 'totalAminoacidSubstitutions', 
	'substitutions', 'deletions', 'insertions', 'aaSubstitutions', 'aaDeletions', 'qc.privateMutations.total',
	'Lineage', 'Ambiguity score', 'Scorpio call', 'Scorpio support', 'pangolin version', 'pangoLEARN version']

	lineages_revised = lineages[included_cols]
	analyzed_data = pd.merge(left = vsalign, right = lineages_revised, how = 'outer', left_on = 'SAMPLE_ID', right_on = 'SAMPLE_ID')

	Samples = analyzed_data['SAMPLE_ID'].apply(str)
	spl1 = []
	for each in Samples:
		value = each.split("_", 1)[1]
		append1 = spl1.append
		append1(value)
	
	samples = []
	for each in spl1:
		value = each.rsplit("_",1)[0]
		append2 = samples.append
		append2(value)
	
	analyzed_data = analyzed_data.drop('SAMPLE_ID', 1)
	analyzed_data['SAMPLE_ID'] = samples

	maximumVal = GetRowID_Analysis(dbname, r_int)
	analyzed_data['ROW_ID'] = analyzed_data.index + maximumVal

	cols = analyzed_data.columns.tolist()
	cols.pop(-1)
	Analyzed_Results = pd.melt(analyzed_data, id_vars=["ROW_ID"],value_vars=cols, var_name='ATTR_ID',value_name="ATTR_VAL")
	Analyzed_Results.loc[:, 'FILE_ID'] = r_int
	Analyzed_Results = Analyzed_Results[['FILE_ID', 'ROW_ID', 'ATTR_ID', 'ATTR_VAL']]

	#Insert into database
	databasename = dbname
	conn=sqlite3.connect(databasename)
	conn.text_factory = str
	cur=conn.cursor()
	Analyzed_Results.to_sql("ANALYZED_DATA_TABLE", conn, if_exists='append', index=False)
	conn.commit()

##########################################################################################
############################## Metadata Table ##############################

def MetadataTable(meta, assign, dbname, r_int):
	with open(meta) as csv_file:
		metadata = pd.read_csv(csv_file, encoding = "unicode_escape", low_memory=False)

	with open(assign) as csv_file:
		assigned_data = pd.read_csv(csv_file, low_memory=False)

	now = datetime.now()

	metadata['SAMPLE_ID'] = metadata['SAMPLE_ID'].str.replace(r'[^A-Za-z0-9_]+', '')
	metadata["FileName"] = "GEIS_Metadata.csv"
	metadata["TIMESTAMP"] = now.strftime("%Y%m%d, %H:%M:%S")
	metadata = metadata.reset_index(drop=True)
	metadata['ROW_ID'] = metadata.index
	cols_metadata = metadata.columns.tolist()
	cols_metadata.pop(-1)
	Metadata_results = pd.melt(metadata, id_vars=["ROW_ID"],value_vars=cols_metadata, var_name='ATTR_ID',value_name="ATTR_VAL")
	Metadata_results.loc[:, 'FILE_ID'] = r_int
	Metadata_results = Metadata_results[['FILE_ID', 'ROW_ID', 'ATTR_ID', 'ATTR_VAL']]
		
	assigned_data['SAMPLE_ID'] = assigned_data['SAMPLE_ID'].str.replace(r'[^A-Za-z0-9_]+', '')
	assigned_data["FileName"] = "GEIS_Assigned_samples.csv"
	assigned_data["TIMESTAMP"] = now.strftime("%Y%m%d, %H:%M:%S")
	assigned_data = assigned_data.reset_index(drop=True)
	assigned_data['ROW_ID'] = assigned_data.index + len(metadata.index)
	cols_assigned_data = assigned_data.columns.tolist()
	cols_assigned_data.pop(-1)
	Assigned_results = pd.melt(assigned_data, id_vars=["ROW_ID"],value_vars=cols_assigned_data, var_name='ATTR_ID',value_name="ATTR_VAL")
	Assigned_results.loc[:, 'FILE_ID'] = r_int
	Assigned_results = Assigned_results[['FILE_ID', 'ROW_ID', 'ATTR_ID', 'ATTR_VAL']]

	#Insert into database
	databasename = dbname
	conn=sqlite3.connect(databasename)
	conn.text_factory = str
	cur=conn.cursor()
	Metadata_results.to_sql("SAMPLE_METADATA_TABLE", conn, if_exists='append', index=False)
	Assigned_results.to_sql("SAMPLE_METADATA_TABLE", conn, if_exists='append', index=False)
	conn.commit()

##########################################################################################
############################## Sample Variants Table ##############################

######### Get sample variants from output #########

def GetVariants(SV_paths, dbname):
	fullList = []
	for each in SV_paths:
		files = []
		files += [name for name in os.listdir(each) if name.endswith('.snp')]
		append3 = fullList.append
		append3(files)

	fullPathsnp = []
	for i in range(len(SV_paths)):
		for each in fullList[i]:
			paths = SV_paths[i] + "/" + each
			append4 = fullPathsnp.append
			append4(paths)

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
		append5 = snpData.append
		append5(data)
		variants = pd.concat(snpData)
	
	return(variants)

######### Input variants info. #########

def VariantsResults(SV_paths, dbname, r_int):
	allvariants = GetVariants(SV_paths, dbname)

	colList = ["SAMPLE_ID","TYPE", "BASE", "SNP","FREQ", "CODON", "FEATURE", "DEPTH","DEL", "INS", "SCD", "NOTE"]
	allvariants = allvariants[colList]
	allvariants.sort_values(by = ['SAMPLE_ID'], inplace = True)
	
	maximumVal = GetRowID_snp(dbname, r_int)
	allvariants['ROW_ID'] = ((allvariants['SAMPLE_ID']!= allvariants['SAMPLE_ID'].shift()).cumsum()) + maximumVal

	final_variants = allvariants[["SAMPLE_ID","TYPE", "BASE", "SNP","FREQ", "CODON", "FEATURE", "DEPTH","DEL", "INS", "SCD", "NOTE", "ROW_ID"]]
	cols = final_variants.columns.tolist()
	cols.pop(-1)

	variants_results = pd.melt(final_variants, id_vars=["ROW_ID"],value_vars=cols, var_name='ATTR_ID',value_name="ATTR_VAL")
			
	variants_results.loc[:, 'FILE_ID'] = r_int
	variants_results = variants_results[['FILE_ID', 'ROW_ID', 'ATTR_ID', 'ATTR_VAL']]

	databasename = dbname
	conn = sqlite3.connect(databasename)
	cur = conn.cursor()
	variants_results.to_sql ("SAMPLE_VARIANTS_TABLE", conn, if_exists='append', index=False)
	conn.commit()

#View updated

