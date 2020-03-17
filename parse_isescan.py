'''
This program runs the program isescan with a genome as input 
and parses de most important fields in the resulting gff file.
Input: genome in fasta format
Output: file parsed as tab delimited

'''
import os
import pandas as pd
from io import StringIO

def run_isescan(genome_file):
	
	# A command to configure isescan program
	ssw_command = 'export LD_LIBRARY_PATH=./ISEScan-master:$LD_LIBRARY_PATH'
	os.system(ssw_command)
	
	# Run the isescan program
	isescan_command = 'python ./ISEScan-master/isescan.py ' + str(genome_file) + ' proteome hmm'
	os.system(isescan_command)
	
def remove_files_isescan():
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
	is_seq_arr = gff_arr[gff_arr['type'] == 'insertion_sequence']
	gff_arr['attribute'].str.replace('.*','IS_\d*')
	print(gff_arr['attribute'])

if __name__ == '__main__':
	#run_isescan('PAO1.fna')
	#remove_files_isescan()
	parse_isescan('./prediction/PAO1.fna.gff')
	
