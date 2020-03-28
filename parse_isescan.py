'''
This program runs the program isescan with a genome as input 
and parses de most important fields in the resulting gff file.
Input: genome in fasta format
Output: file parsed as csv

'''
import os
import pandas as pd
import numpy as np

def run_isescan(genome_file):
	
	# A command to configure isescan program
	ssw_command = 'export LD_LIBRARY_PATH=./ISEScan-master:$LD_LIBRARY_PATH'
	os.system(ssw_command)
	
	# Run the isescan program
	isescan_command = 'python ./ISEScan-master/isescan.py ' + str(genome_file) + ' ./ISEScan-master/proteome ./ISEScan-master/hmm'
	os.system(isescan_command)
	
def remove_files_isescan():
	#DO NOT USE
	# This function deletes all output files from isescan that are not usefull,
	# only the gff file remains
	delete_command = 'find . -type f ! -name "*.gff" -delete'
	os.system(delete_command)
	
def parse_isescan(gff_file):
	'''
	This function takes as input the gff file from the isescan program and creates a new parsed file with the
	folowing tab delimited fields:
	id:IS_X
	Family
	Cluster
	Element type
	Start coordinate
	Finish coordinate
	Orientation
	'''
	gff_arr = pd.read_table(gff_file,
							names = ['accesion', 'program', 'type', 
									'start', 'end','score', 'strand', 'frame', 'attribute'],
							skiprows = [0]
							)
							
	#Filter the datafram to just is_sequence position
	is_seq_arr = gff_arr[gff_arr['type'] == 'insertion_sequence']
	
	# Create a new column with the corresponding id in format IS_1, IS_2 ....
	gff_arr['id'] = gff_arr['attribute'].str.extract('(IS_\d+)')
	
	#Create new column with corresponding family in format ISX, ISY ....
	gff_arr['family'] = gff_arr['attribute'].str.split(';', expand = True)[1]
	gff_arr['family'] = gff_arr['family'].str.extract('(IS\d+)')
	#gff_arr.loc[:,'family'][gff_arr['family'] == 'IS'] = np.nan
	
	#Create new column with corresponding cluster in format ISX, ISY ....
	gff_arr['cluster'] = gff_arr['attribute'].str.split(';', expand = True)[2]
	gff_arr['cluster'] = gff_arr['cluster'].str.extract('(IS\d+_\d+)')
	#gff_arr.loc[:,'cluster'][gff_arr['cluster'] == 'IS'] = np.nan
	
	#Add an interval column to later merge dataframes
	gff_arr['interval'] = [pd.Interval(start, end, closed='both') for (start, end) in zip(gff_arr['start'], gff_arr['end'])]
	
	#Create the dataframe with the wanted order
	gff = gff_arr[['id','family','cluster','type','start','end','strand', 'interval']]
	
	
	return gff
	
def write_as_csv(gff_df, genome_name =''):
	'''This function takes a dataframe and outputs a csv_file'''
	gff_df.to_csv('./results/' + genome_name + '_is.csv')
	

if __name__ == '__main__':
	#run_isescan('NCGM2.S1.fna')
	#remove_files_isescan()
	NCGM2 = parse_isescan('./prediction/NCGM2.S1.fna.gff')
	write_as_csv(NCGM2, 'NCGM2.S1')
	
