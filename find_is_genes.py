'''
This program imports a csv file with is sequence locations and gene locations
and produces a csv file with the genes that are interrupted by is_sequences
also the ones that between 2 genes
'''
from parse_isescan import parse_isescan
from parse_genome import parse_genome
import numpy as np
import pandas as pd
import pyranges as pr

def create_csv_files(genome_name):
	'''
	TODO:This function calls the IS and genome parsers and produces the respective csv files
	'''
	parse_isescan(genome_name)
	parse_genome(genome_name)

def import_csv_files(genome_name):
	is_df = pd.read_csv('./results/' + genome_name + '_is.csv')
	genome_df = pd.read_csv('./results/' + genome_name + '_genome.csv')
	
	return is_df, genome_df

def find_is_genes(genome_name):
	
	is_df, genome_df = import_csv_files(genome_name)
	
	#Select only insertion sequence type
	is_df = is_df[is_df['type'] == 'insertion_sequence']
	
	
	is_range = pr.PyRanges(chromosomes=genome_name,starts=is_df['start'], ends=is_df['end'],strands=is_df['strand'])
	
	genome_range = pr.PyRanges(chromosomes=genome_name,starts=genome_df['start'], ends=genome_df['end'],strands=genome_df['strand'])
	genome_range.id = genome_df['id']
	
	intersected_df = genome_range.overlap(is_range).as_df()
	
	merged_df = pd.merge(genome_df, intersected_df, on='id')
	print('There are ' + str(len(merged_df)) + ' IS sequences interrupting genes in ' + genome_name)
	print('The genes are:')
	print(merged_df['gene_name'])
	
	
	
if __name__ == '__main__':
	find_is_genes('PAO1')
	
