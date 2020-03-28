'''
This program recieves a gff file of a genome 
and parses de most important fields in the resulting gff file.
Input: genome in fasta format
Output: file parsed as csv
'''
import os
import pandas as pd
import numpy as np

def parse_genome(gff_file, only_genes = True):

	'''
	This function takes as input a gff_file from a genome and creates a pandas df
	with most relevant features: ID, start, finish, gene_name, orientation, function, locus tag
	Also, we only use the gene locations, the rest genomic elements are not used
	'''
	# Read gff file and convert it to pd dataframe with all its columns
	gff_arr = pd.read_table(gff_file,
							names = ['accesion', 'program', 'type', 
									'start', 'end','score', 'strand', 'frame', 'attribute'],
							skiprows = [0],
							comment = '#'
							)
							
	#Keep only the rows with gene products
	if only_genes:
		gff_arr = gff_arr[gff_arr['type'] == 'gene']

	# Create a separate series with all the attributes to later parse	(?<= = ).* 				
	attribute_arr = gff_arr['attribute']
	
	#Parse ID (match a pattern that starts with ID= until;
	gff_arr['id'] = attribute_arr.str.extract('ID=(.*?);')
	
	#Parse gene name
	#gff_arr['product'] = attribute_arr[4].str.extract('gene=(.*)')
	gff_arr['gene_name'] = attribute_arr.str.extract('Name=(.*?);')
	
	#TODO:Parse protein function
	
	#TODO:Extract the corresponding protein attribute from the gene it belongs
	
	#Parse locus tag
	gff_arr['locus_tag'] = attribute_arr.str.extract('locus_tag=(.*)')
	
	#Add an interval column to later merge dataframes
	gff_arr['interval'] = [pd.Interval(start, end, closed='both') for (start, end) in zip(gff_arr['start'], gff_arr['end'])]
	#Create the dataframe with the wanted order
	gff = gff_arr[['id', 'start','end','gene_name','strand', 'locus_tag', 'interval']]
	
	return gff

def write_as_csv(gff_df, genome_name =''):
	'''This function takes a dataframe and outputs a csv_file'''
	gff_df.to_csv('./results/' + genome_name + '_genome.csv')

if __name__ == '__main__':
	NCGM2 = parse_genome('./genomes/NCGM2.S1.gff')
	write_as_csv(NCGM2, 'NCGM2.S1')
