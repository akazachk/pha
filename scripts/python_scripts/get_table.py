## Arguments
# First argument: boolean (0: do not recompute, 1: recompute, 2: recompute and save best runs), second argument is "test", "best", or "full"

## Set up proper path variables
from os import environ
PHA_DIR = environ['PHA_DIR']
#PHA_DIR = '../..'
ADD_TO_PATH = PHA_DIR + '/scripts/python_scripts'

from sys import path
path.append(ADD_TO_PATH)

## Now import all the relevant code
from utility import *
from lgreader import LGReader
from matrix2latex import matrix2latex
import csv

## Get arguments
from sys import argv
is_test = False
is_best = False
is_full = True
recompute = False
should_save_best_run_data = False

if (len(argv) > 1):
	try:
		tmp_recompute = int(argv[1])
	except:
		raise TypeError( "First argument should be an integer." )

	if (tmp_recompute not in [0,1,2]):
		raise ValueError( "recompute needs to be 0, 1, or 2, but given is: %d." % tmp_recompute )
	else:
		should_save_best_run_data = should_save_best_run_data or (tmp_recompute == 2)
		recompute = (tmp_recompute >= 1)
	
	if (len(argv) > 2):
		is_test = True if argv[2].lower() in ['test'] else False
		is_best = True if argv[2].lower() in ['best'] else False
		is_full = True if argv[2].lower() in ['full'] else False

## In and out directories
in_dir = PHA_DIR + '/results'
out_dir = PHA_DIR + '/results/tables'
if is_test:
	f_name = in_dir + '/pha-test.csv'
if is_best:
	f_name = in_dir + '/pha-best.csv'
if is_full:
	f_name = in_dir + '/pha-full.csv'

## Instances
zero = [
	'glass4',
	'misc01',
	'misc02',
	'misc03',
	'misc07',
	'p0201',
	'rgn'
]

instance = [
	'bell3a',
	'bell3b',
	'bell4',
	'bell5',
	'blend2',
	'bm23',
	'egout',
	'flugpl',
#	'glass4',
	'gt2',
	'k16x240',
	'lseu',
	'mas74',
	'mas76',
	'mas284',
#	'misc01',
#	'misc02',
#	'misc03',
	'misc05',
#	'misc07',
	'mod008',
	'mod013',
	'modglob',
	'p0033',
	'p0040',
#	'p0201',
	'p0282',
	'p0291',
	'pipex',
	'pp08a',
	'probportfolio',
#	'rgn',
	'sample2',
	'sentoy',
	'stein15_nosym',
	'stein27_nosym',
	'stein45_nosym',
	'timtab1',
	'vpm1',
	'vpm2'
]
instance.extend(zero)

instance_shortlist = [
	'bell3a',
	'bell3b',
	'bell4',
	'bell5',
	'blend2',
	'bm23',
	'egout',
	'flugpl', # 0 gap closed
	'gt2',
	'k16x240',
	'lseu',
	'mas74',
	'mas76',
	'mas284',
	'misc05',
	'mod008',
	'mod013',
	'modglob',
	'p0033',
	'p0040', # 0 gap closed
	'p0282',
	'p0291',
	'pipex',
	'pp08a',
	'probportfolio',
	'sample2',
	'sentoy',
	'stein15_nosym',
	'stein27_nosym',
	'stein45_nosym',
	'timtab1', # 0 gap closed
	'vpm1',
	'vpm2'
]
shortlist_indices = [
				i for i in range(len(instance)) 
				if instance[i] in instance_shortlist
		]
len_instancelist = len(instance)
len_shortlist = len(shortlist_indices)

## Table names
table_name_list = [
				"inst-info",	# 0
				"gap-closed",	# 1
				"obj-fails",	# 2 
				"point-analysis",	# 3
				"hplane-analysis",	# 4
				"cut-heur-analysis",	# 5
				"active-cuts",	# 6
				"param-analysis"	# 7
		]
num_tables = len(table_name_list)

## Tables to ignore; those that are commented out are NOT ignored
ignore_table_list = [
				0,	# 0, inst-info
				#1,	# 1, gap-closed
				#2,	# 2, obj-fails
				#3,	# 3, point-analysis
				#4,	# 4, hplane-analysis
				#5,	# 5, cut-heur-analysis
				#6,	# 6, active-cuts
				## DO NOT USE ## 
				7	 # 7, param-analysis
		]

if is_best:
	ignore_table_list.extend([2,3,4,5,6,7])

## Lists for making the tables
# Whether to use the short or long instance list
is_short_list = [
				False,	# 0, inst-info
				False,	# 1, gap-closed
				False,	# 2, obj-fails
				True,	 # 3, point-analysis
				True,	 # 4, hplane-analysis
				True,	 # 5, cut-heur-analysis
				True,	 # 6, active-cuts
				False,	# 7, param-analysis
		]

# Caption
caption_list = [
				"Set of instances used in experiments, along with the number of rows and columns for each instance, and the average amount of time.",	# 0, inst-info
				"Best percent gap closed and number of cuts.",	# 1, gap-closed
				"Number of cuts generated.",	# 2, obj-fails
				"Number of points and final points generated for instances with any improvement over SICs.",	# 3, point-analysis
				"Best percent gap closed by hyperplane activation choice for instances with any improvement over SICs. Some values differ in the thousandths digit.",	# 4, hplane-analysis
				#"Best percent gap closed by objective function used for instances with any improvement over SICs across all settings. Some values differ in the thousandths digit.",	# 5, cut-heur-analysis
				"Best percent gap closed by objective function used for the 30 instances with any improvement over SICs across all settings. \
				Highlighted cells indicate objectives (or combinations of objectives) achieving best result. \
		For the combinations, the cell is highlighted only if using the combination is actually necessary to achieve the result. \
				Some values differ in the thousandths digit.",	# 5, cut-heur-analysis
				"Active cuts for instances with improvements over SICs by heuristic used to generate them.",	# 6, active-cuts
				"Parameters."	# 7, param-analysis
		]

# Label
label_list = [
				"tab:"+table_name_list[i]
				for i in range(num_tables)
		]

# Alignment of columns
align_list = [
				r"lccc",	# 0, inst-info
				r"l*{2}{H}*{4}{H}cc>{\bfseries}c*{4}{c}H",	# 1, gap-closed
				r"lHcccccHHccccHccHcH",	# 2, obj-fails
				r"lccc",	# 3, point-analysis
				#r"lccH|*{3}{c}|*{5}{c}|*{4}{c}*{24}{H}",	# 4, hplane-analysis
				r"lccH|*{3}{c}|*{5}{c}|*{4}{c}",	# 4, hplane-analysis
				r"lc|*{4}{c}|*{6}{c}|*{4}{c}",	# 5, cut-heur-analysis
				r"lcccccccccccccc",	# 6, active-cuts
				r"lccccccccc"	# 7, param-analysis
		]

# Format for each of the columns
float_style0 = "%.0f"
float_style1 = "%.1f"
float_style2 = "%.2f"
int_style = "%d"
gen_style = "%g"
format_list = [
				[gen_style, int_style, int_style, float_style1],	# 0, inst-info
				[gen_style] + 2*[int_style] + 6*[float_style2] + [float_style2] + 4*[int_style] + [float_style0],	# 1, gap-closed
				int_style,	# 2, obj-fails 
				[gen_style] + 8*[int_style] + [float_style0],	# 3, point-analysis
				float_style2,	# 4, hplane-analyis
				float_style2,	# 5, cut-heur-analysis
				int_style,	# 6, active-cuts
				float_style2	# 7, param-analysis
		]

# Number of headers for the table (to decide where to put \midrule)
header_length_list = [
				1,	# 0, inst-info
				2,	# 1, gap-closed
				1,	# 2, obj-fails
				1,	# 3, point-analysis
				2,	# 4, hplane-analysis
				1,	# 5, cut-heur-analysis
				1,	# 6, active-cuts
				1	 # 7, param-analysis
		]

# Number of footers for the table (to decide where to put the last midrule, before the average rows)
summaryrows_list = [
				0,	# 0, inst-info
				1,	# 1, gap-closed *
				0,	# 2, obj-fails
				0,	# 3, point-analysis
				1,	# 4, hplane-analyis
				1,	# 5, cut-heur-analysis
				0,	# 6, active-cuts
				1	 # 7, param-analysis
		]

# Where to put other midrules
midruleIndex_list = [
				[],	# 0, inst-info
				[len(instance_shortlist), len(instance_shortlist)+1],	# 1, gap-closed
				[],	# 2, obj-fails
				[],	# 3, point-analyis
				[],	# 4, hplane-analysis
				[],	# 5, cut-heur-analysis
				[],	# 6, active-cuts
				[]	 # 7, param-analysis
		]

columnIndices_list = [
				[], # 0, inst-info
				[], # 1, gap-closed
				#[1,8,9,10,11,12,13,14]
				[], # 2, cut-analysis
				[], # 3, point-analysis
				range(0,16), # 4, hplane-analysis
				[], # 5, cut-heur-analysis
				[], # 6, active-cuts
				[], # 7, param-analysis
		]

avg_col_indices = [
				[],						 # 0, inst-info
				[7,8,9,14],		 # 1, gap-closed
				[],						 # 2, cut-analysis
				[],	#[8],			# 3, point-analysis
				range(1,16),		# 4, hplane-analysis
				range(1,16),		# 5, cut-heur-analysis
				[],						 # 6, active-cuts
				range(1,10)		 # 7, param-analysis
		]

def chain_elements_or_slices(*elements_or_slices):
		new_list = []
		for i in elements_or_slices:
				if isinstance(i, list):
						new_list.extend(i)
				else:
						new_list.append(i)
		return new_list

def save_tab(curr_tab, i):
		'''
		Save csv and tex documents
		'''
		start_table = header_length_list[i]
		curr_stub = out_dir + '/' + table_name_list[i]
		if is_short_list[i]:
				curr_stub += "_short"
				tmp_tab = chain_elements_or_slices(
								curr_tab[0:start_table], 
								#[curr_tab[j+start_table] for j in shortlist_indices]
								[curr_tab[j]
												for j in range(len(curr_tab)) if curr_tab[j][0] in instance_shortlist
										]
								)
				curr_tab = tmp_tab
		if (len(avg_col_indices[i]) > 0):
				# Add new average row
				len_list = len_shortlist if is_short_list[i] else len_instancelist
				avg = [
								sum([float(row[j]) for row in curr_tab[start_table:]]) / len_list
								if j in avg_col_indices[i]
								else ''
								for j in range(len(curr_tab[0]))
						]
				avg[0] = "Average"
				curr_tab.append(avg)
		out_name = curr_stub + '.csv'
		if (recompute):
				print_table(curr_tab, sep=',', outfile=out_name, overwrite=True)
		out_name = curr_stub + '.tex'
		with open(out_name, 'wb') as out_f:
				out_f.write(
						matrix2latex(curr_tab[start_table:], None,
						"singlespace", "adjustbox", "tabular",
						headerRow=curr_tab[0:start_table], 
						alignment=align_list[i], 
						label=label_list[i], 
						formatColumn=format_list[i],
						summaryrows = summaryrows_list[i],
						midruleIndex = midruleIndex_list[i],
						columnIndices = columnIndices_list[i],
						caption=caption_list[i])
				)

if __name__ == "__main__":
		tab = []
		if (recompute):
				reader = LGReader(f_name, compute_gap_closed = False)
				reader.set_param("CUT_LIMIT",1000)
				reader.set_param("CUT_PRESOLVE",0)
				reader.set_inst(instance, compute_gap_closed = False)
				reader.get_best_row()
				if (should_save_best_run_data):
					# Save the best run data
					print( "## Saving best runs to file ##" )
					reader.write_best_params(PHA_DIR + "/data/params")
					#with open(out_dir + '/' + 'best_runs.csv','wb') as myfile:
					#	wr = csv.writer(myfile)
					#	wr.writerow(reader._header[0])
					#	wr.writerow(reader._header[1])
					#	for inst in range(len(instance)):
					#		print( "\tBest run: instance %s." % instance[inst] )
					#		tmp = reader.get_row(reader._bestrow[inst])
					#		wr.writerow(tmp)
				if (0 in ignore_table_list):
						tab.append([])
				else:
						curr_tab = reader.inst_info()
						tab.append(curr_tab) 
						save_tab(curr_tab, 0)
				if (1 in ignore_table_list):
						tab.append([])
				else:
						curr_tab = reader.gap_closed_table()
						tab.append(curr_tab) 
						save_tab(curr_tab, 1)
				if (2 in ignore_table_list):
						tab.append([])
				else:
						curr_tab = reader.obj_fails_table()
						tab.append(curr_tab) 
						save_tab(curr_tab, 2)
				if (3 in ignore_table_list):
						tab.append([])
				else:
						curr_tab = reader.point_analysis()
						tab.append(curr_tab) 
						save_tab(curr_tab, 3)
				if (4 in ignore_table_list):
						tab.append([])
				else:
						curr_tab = reader.hplane_analysis()
						tab.append(curr_tab) 
						save_tab(curr_tab, 4)
				if (5 in ignore_table_list):
						tab.append([])
				else:
						curr_tab = reader.cut_heur_analysis()
						tab.append(curr_tab) 
						save_tab(curr_tab, 5)
				if (6 in ignore_table_list):
						tab.append([])
				else:
						curr_tab = reader.active_cuts_table()
						tab.append(curr_tab) 
						save_tab(curr_tab, 6)
				if (7 in ignore_table_list):
						tab.append([])
				else:
						curr_tab = reader.param_analysis()
						tab.append(curr_tab) 
						save_tab(curr_tab, 7)
		else:
				import re
				for i in range(num_tables):
						if (i in ignore_table_list):
								continue
						curr_stub = out_dir + '/' + table_name_list[i]
						if is_short_list[i]:
								curr_stub += "_short"
						curr_tab = csv2table(curr_stub + '.csv')
						save_tab(curr_tab, i)
