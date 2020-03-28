'''
This class finds the proteins and genes in a genome that are interrupted by an IS
'''
import pandas as pd
import numpy as np
import pyranges as pr

class is_parser():
	
	def __init__(self, genome_name, is_gff_file, genome_gff_file):
		
		self.is_gff_file = is_gff_file
		self.genome_gff_file = genome_gff_file
		self.genome_name = genome_name
		self.is_df = self.parse_is_file()
		self.genome_df = self.parse_genome_file()
		
	
	def read_gff_file(self, gff_file):
		'''
		This function reads a gff file and converts it to pandas df
		'''
		gff_df = pd.read_table(gff_file,
								names = ['accesion', 'program', 'type', 
								'start', 'end','score', 'strand', 'frame', 'attribute'],
								comment = '#'
								)
		return gff_df
		
	def parse_is_file(self):
		
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
		#Read the insertion sequence gff file
		gff_arr = self.read_gff_file(self.is_gff_file)
								
		#Filter the dataframe to just is_sequence position
		gff_arr = gff_arr[gff_arr['type'] == 'insertion_sequence']
		
		# Create a separate series with all the attributes to later parse	(?<= = ).* 				
		attribute_arr = gff_arr['attribute']
		
		#Parse ID (match a pattern that starts with ID= until;
		gff_arr['id'] = attribute_arr.str.extract('ID=.*?(IS_.*?);')
		
		#Parse family (match a pattern that starts with family= until;
		gff_arr['family'] = attribute_arr.str.extract('family=(.*?);')
		
		#Parse cluster (match a pattern that starts with ID= until the end
		gff_arr['cluster'] = attribute_arr.str.extract('cluster=(.*)')
		
		#Add an interval column to later merge dataframes
		gff_arr['interval'] = [pd.Interval(start, end, closed='both') for (start, end) in zip(gff_arr['start'], gff_arr['end'])]
		
		#Create the dataframe with the wanted order
		gff = gff_arr[['id','family','cluster','type','start','end','strand', 'interval']]
		
		return gff
		
		
	def parse_genome_file(self):
		'''
		This function takes as input a gff_file from a genome and creates a pandas df
		with most relevant features: ID, start, finish, gene_name, orientation, function, locus tag
		Also, we only use the gene locations, the rest genomic elements are not used
		'''
		# Read genome gff file and convert it to pd dataframe with all its columns
		gff_arr = self.read_gff_file(self.genome_gff_file)
								
		#Select just the CDS of the file
		gff_arr = gff_arr[gff_arr['type'] == 'CDS']

		# Create a separate series with all the attributes to later parse	(?<= = ).* 				
		attribute_arr = gff_arr['attribute']
		
		#Parse ID (match a pattern that starts with ID= until;
		gff_arr['id'] = attribute_arr.str.extract('protein_id=(.*?);')
		
		#Parse gene name
		gff_arr['product'] = attribute_arr.str.extract('product=(.*?);')
		
		gff_arr['gene_name'] = attribute_arr.str.extract('gene=(.*?);')
		
		gff_arr['parent_gene'] = attribute_arr.str.extract('Parent=(.*?);')
		
		#Parse locus tag
		gff_arr['locus_tag'] = attribute_arr.str.extract('locus_tag=(.*?);')
		
		#Add an interval column to later merge dataframes
		gff_arr['interval'] = [pd.Interval(start, end, closed='both') for (start, end) in zip(gff_arr['start'], gff_arr['end'])]
		#Create the dataframe with the wanted order
		gff = gff_arr[['id', 'start','end','strand','product','gene_name','parent_gene', 'locus_tag', 'interval']]
		
		return gff
		
		
		
	def find_is_genes(self):
		
			is_range = pr.PyRanges(chromosomes=self.genome_name,starts=self.is_df['start'], ends=self.is_df['end'],strands=self.is_df['strand'])
			
			genome_range = pr.PyRanges(chromosomes=self.genome_name,starts=self.genome_df['start'], ends=self.genome_df['end'],strands=self.genome_df['strand'])
			
			genome_range.id = self.genome_df['id']
			
			intersected_df = genome_range.overlap(is_range).as_df()
			
			merged_df = pd.merge(self.genome_df, intersected_df, on='id')
			print('There are ' + str(len(merged_df)) + ' IS sequences interrupting genes in ' + self.genome_name)
			print('The genes are:')
			print(merged_df[['product','gene_name']])
			
	def is_families(self):
		'''
		This method returns all the is families present in the genome
		'''
		return self.is_df['family']
	
	def group_by_is_family(self):
		'''
		this method returns a df with all the families and gene products
		'''
		for name, group in self.is_df.groupby(['family'])['id']:
			print(name)
			print(group)
			print('\n')
			
	def test(self):
		print(self.group_by_is_family())
	
		
			
			
if __name__ == '__main__':
	PAO1_parser = is_parser('PAO1','./prediction/PAO1.fna.gff','./genomes/PAO1_genome.gff')
	NCGM2_parser = is_parser('NCGM2.S1','./prediction/NCGM2.S1.fna.gff','./genomes/NCGM2.S1.gff')
	PAO1_parser.test()
	NCGM2_parser.test()
		
