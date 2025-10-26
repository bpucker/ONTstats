### Boas Pucker ###
### pucker@uni-bonn.de ###

__version__ = "v0.1"

__usage__ = """
			Evaluation of all ONT sequencing runs (""" + __version__ + """)
			python3 evaluate_ONT_P2_runs.py
			--in <INPUT_FOLDER>
			--out <OUTPUT_FOLDER>
			"""


import os, sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from operator import itemgetter
# --- end of imports --- #


def load_data( input_file, masking ):
	"""! @brief load data from given input file """
	
	data = []
	data_per_project = {}
	data_per_flowcell = {}
	counter = 1	#only relevant for masking FC IDs
	fc_mapping_table = {}	#only relevant for masking FC IDs
	
	with open( input_file, "r" ) as f:
		headers = f.readline().strip().split('\t')
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			if len( parts ) >= 10:	#check if complete line
				try:
					# --- load all data from a row --- #
					entry = {}
					for head in headers:
						if head in [ "ReadN50[kb]", "Output[Gbp]" ]:
							try:
								entry.update( { head: int( parts[ headers.index( head ) ] ) } )
							except:
								try:
									entry.update( { head: float( parts[ headers.index( head ) ] ) } )
								except:
									entry.update( { head: 0 } )
						else:
							entry.update( { head: parts[ headers.index( head ) ] } )
					data.append( entry )
					
					# --- get values --- #
					try:
						output = float( parts[5] )
					except ValueError:
						output = 0
					try:
						N50 = float( parts[6] )
					except ValueError:
						N50 = 0.0
					try:
						reads = int( parts[7] )
					except ValueError:
						reads = 0
					
					# --- load data per project --- #
					try:	#project ID is given in column1
						data_per_project[ parts[1] ]['output'].append( output )
						data_per_project[ parts[1] ]['N50'].append( N50 )
						data_per_project[ parts[1] ]['reads'].append( reads )
					except KeyError:
						data_per_project.update( { parts[1]: { 	'output': [ output ],
																'N50': [ N50 ],
																'reads': [ reads ]
																} } )
					# --- load data per flowcell --- #
					if masking:
						try:
							fcid = fc_mapping_table[ parts[3] ]
						except KeyError:
							fcid = "FC" + str( counter )
							fc_mapping_table.update( { parts[3]: fcid } )
							counter += 1
					else:
						fcid = parts[3] + ""
					
					try:	#flowcell ID is given in column4
						data_per_flowcell[ fcid ]['output'].append( output )
						data_per_flowcell[ fcid ]['N50'].append( N50 )
						data_per_flowcell[ fcid ]['reads'].append( reads )
					except KeyError:
						data_per_flowcell.update( { fcid: { 	'output': [ output ],
																	'N50': [ N50 ],
																	'reads': [ reads ]
																	} } )
				except:
					print( "ERROR: " + line )
			line = f.readline()
	return data, data_per_project, data_per_flowcell


def generate_output_figure( dataset, figfilename, parameter ):
	"""! @brief generate stacked barplots """
	
	sample_names = sorted( list( dataset.keys() ) )
	
	data = []
	max_len = 0
	for sample in sample_names:
		data.append( [ sample ] + dataset[ sample ][ parameter ] )
		if len( [ sample ] + dataset[ sample ][ parameter ] ) > max_len:
			max_len = len( [ sample ] + dataset[ sample ][ parameter ] )
	
	df = pd.DataFrame( data, columns=[ 'Sample' ] + list( map( str, range( max_len ) ) )[1:] )
	if parameter == "output":
		ax = df.plot( x="Sample", kind="bar", stacked=True, title=parameter, legend=False )
		ax.set_ylabel( "Output [Gbp]" )
	elif parameter == "reads":
		ax = df.plot( x="Sample", kind="bar", stacked=True, title=parameter, legend=False )
		ax.set_ylabel( "thousands reads" )
	elif parameter == "N50":
		ax = df.plot( x="Sample", kind="bar", stacked=True, title=parameter, legend=False )
		ax.set_ylabel( "N50[kbp]" )
	else:
		ax = df.plot( x="Sample", kind="bar", stacked=True, title=parameter, legend=False )
	plt.savefig( figfilename, bbox_inches='tight', dpi=300 )


def highscore( data, highscore_file ):
	"""! @brief identify best performing runs, FCs, and more """
	
	with open( highscore_file, "w" ) as out:
		highest_output_runs = sorted( data, key=itemgetter( "Output[Gbp]" ) )[::-1]
		for each in highest_output_runs[:5]:
			out.write( each["#Run"] + "\t" + str( each["Output[Gbp]"] ) + "\n" )
		out.write( "\n\n" )
		highest_N50_runs = sorted( data, key=itemgetter( "ReadN50[kb]" ) )[::-1]
		for each in highest_N50_runs[:5]:
			out.write( each["#Run"] + "\t" + str( each["ReadN50[kb]"] ) + "\n" )


def main( arguments ):
	"""! @brief run everything """
	
	input_file = arguments[ arguments.index('--in')+1 ]
	output_folder = arguments[ arguments.index('--out')+1 ]
	
	if '--nomasking' in arguments:
		masking = False
	else:
		masking = True
	
	if output_folder[-1] != "/":
		output_folder += "/"
	if not os.path.exists( output_folder ):
		os.makedirs( output_folder )
	
	data, data_per_project, data_per_flowcell = load_data( input_file, masking )
	
	project_output_figure = output_folder + "output_per_project.output.png"
	generate_output_figure( data_per_project, project_output_figure, "output" )
	#project_output_figure = output_folder + "output_per_project.reads.png"
	#generate_output_figure( data_per_project, project_output_figure, "reads" )
	#project_output_figure = output_folder + "output_per_project.N50.png"
	#generate_output_figure( data_per_project, project_output_figure, "N50" )
	
	flowcell_output_figure = output_folder + "output_per_flowcell.output.png"
	generate_output_figure( data_per_flowcell, flowcell_output_figure, "output" )
	#flowcell_output_figure = output_folder + "output_per_flowcell.reads.png"
	#generate_output_figure( data_per_flowcell, flowcell_output_figure, "reads" )
	#flowcell_output_figure = output_folder + "output_per_flowcell.N50.png"
	#generate_output_figure( data_per_flowcell, flowcell_output_figure, "N50" )
	
	highscore_file = output_folder + "highscore.txt"
	highscore( data, highscore_file )


if '--in' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
