# File for reading output from csv files
# Specific to GIC output files by LiftGICs

from __future__ import print_function
import csv
import utility as util
from statenum import StatEnum
import copy
import param

class ProcessCSV(object):
		"""
		Class ProcessCSV:
				Tailored to GIC class but can be adapted based on your needs
				Would be better to use kwargs in the stat_col function...
		"""
		## Class variables
		default_num_header_lines = 2 # number of lines in the file that correspond to the header
		default_col_name_row = 1 # where to look for column names (in header lines), index starts at 0
		inst_col = 0 # where to look for instance names (row names)
		default_num_dec_places = 10 #2

#		format_float = lambda x: "{:." + str(x) + "f}"
		def __init__(self, f_name, inst_name = None, col_info = None, num_header_lines = None, 
									compute_gap_closed = True):
				""" Constructor """
				if __debug__:
						print( "\n## In ProcessCSV(): Opening file %s for reading ##" % f_name )
				self._f_name = f_name
				self._f = open(f_name, 'r')
				self._reader = csv.reader(self._f)
				
				# Number decimal places to use for floats
				self._num_dec_places = self.default_num_dec_places

				if (num_header_lines is not None):
						self._num_header_lines = num_header_lines
						self._col_name_row = self._num_header_lines - 1 # assuming the last line of the headers is the relevant one
				else:
						self._num_header_lines = self.default_num_header_lines
						self._col_name_row = self.default_col_name_row

				self.read_header()
				if __debug__:
						print( "## In ProcessCSV(): Read in %s, with %d header lines ##" % (f_name, self._num_header_lines) )

				# Set inst_name if it is specified 
				# Do not get all instances, in case user wishes to specify later
				self._param_container = param.Param()
				if (inst_name is not None):
						self.set_inst(inst_name, compute_gap_closed)

				# If column info specified, input it
				if col_info is not None:
						self.set_col(col_info)

		def read_header(self):
				""" Reads in header if at line 0, but does not reset to line 0 if not there currently """
				if (self._reader.line_num == 0):
						self._header = []
						for i in range(self._num_header_lines):
								self._header.append( self._reader.next() )
		
		def set_inst(self, inst_name = None, compute_gap_closed = True):
				""" Set instances to be read """
				if __debug__:
						print( "## Setting instance names ##" )
				if ((inst_name is None) or (inst_name == [])):
						# Read unique set of instances from column inst_col of the file
						tmp = []
						self.restart_reader(False)
						for row in self._reader:
								tmp.append(row[self.inst_col])
						inst_name = {}.fromkeys(tmp).keys()
						inst_name.sort()

				# Change instances in Param class
				self._param_container.set_param("INSTANCE", inst_name, True, str)
				self._nInst = len(self._param_container.get_param("INSTANCE"))
				
				if __debug__:
						print( "## Set instance list to ", end='' )
						print( self._param_container.get_param("INSTANCE"), end='' )
						print( " ##" )

		def get_inst(self, inst):
				""" Return name or index of inst """
				return self._param_container.get_inst(inst)

		def set_param(self, param_name, param_val):
				""" Set parameter param_name to param_val """
				if (param_name == "INSTANCE"):
						self.set_inst(param_val)
				else:
						self._param_container.set_param(param_name, param_val, True)
						if __debug__:
								print( "## Set parameter %s to " % param_name, end='' )
								print( self._param_container.get_param(param_name), end='' )
								print( " ##" )

		def get_param(self, which_param):
				if (type(which_param) is not str) and (type(which_param) is not int):
						raise TypeError(
								"Param name needs to be of type str or int, but is of type %s."
								% type(which_param)
						)

				return self._param_container.get_param(which_param)
				
		def set_col(self, col_info, set_val = None, is_secondary = False):
				""" 
				Find the index of the column for which statistics are being collected 
				Returns [nTargetCols, col_index, col_name] list
				"""
				if (set_val is None):
						set_val = True

				if (set_val):
						if __debug__:
								print( "## Setting column info ##" )

				# If nothing is specified, return the previously set col_info
				# Assumes this been set already
				if (col_info is None):
						if is_secondary:
								return [self._nTargetCols, self._secondary_col_index, self._secondary_col_name]
						else:
								return [self._nTargetCols, self._col_index, self._col_name]

				# In case header has not been read
				self.read_header()

				# If col_info is not a list, then make it one
				if (type(col_info) is not list):
						col_info = [col_info]

				nTargetCols = len(col_info)
				col_name = ["" for i in range(nTargetCols)]
				col_index = [-1 for i in range(nTargetCols)]
				for j in range(nTargetCols):
						# If col_info is a str, it is the name of the column and we must find the index
						if (type(col_info[j]) is str):
								col_name[j] = col_info[j]
								col_index[j] = self.get_col_index(col_name[j])
						elif (type(col_info[j]) is int):
								col_index[j] = col_info[j]
								col_name[j] = self.get_col_name(col_index[j])
						elif (type(col_info[j]) is list):	# may also be a list, in which case we repeat above for each col in list
								col_name[j] = [""] * len(col_info[j])
								col_index[j] = [-1] * len(col_info[j])
								for c in range(len(col_info[j])):
										if (type(col_info[j][c]) is str):
												col_name[j][c] = col_info[j][c]
												col_index[j][c] = self.get_col_index(col_name[j][c])
										elif (type(col_info[j][c]) is int):
												col_index[j][c] = col_info[j][c]
												col_name[j][c] = self.get_col_name(col_index[j][c])
										else:
												raise TypeError("Column specified is neither a string, for the column name, nor an integer, for the column index.")
						elif (is_secondary and col_info[j] is None):
								pass
						else:
								raise TypeError("Column specified is neither a string, for the column name, nor an integer, for the column index, nor a list of the above.")
						if (type(col_index[j]) is int and col_index[j] < 0):
								raise KeyError("Did not find column %s in header" % col_name[j])
						if __debug__:
								if (type(col_info[j]) is not list):
										print( '## Found column %s at index %d ##' % (col_name[j], col_index[j]) )
								else:
										print( '## Found columns (%s,%s) at indices (%d,%d) ##' % (col_name[j][0], col_name[j][1], col_index[j][0], col_index[j][1]) )

				if (set_val):
						if is_secondary:
								self._secondary_col_name = copy.deepcopy(col_name)
								self._secondary_col_index = copy.deepcopy(col_index)
						else:
								self._col_name = copy.deepcopy(col_name)
								self._col_index = copy.deepcopy(col_index)
								self._nTargetCols = len(col_info)
				return [nTargetCols, col_index, col_name]
						
		def init_stat(self, stat, nTargetCols, nInst, typecol, is_secondary = False):
				""" 
				Initialize stat we are keeping:
				avg: stat = 0
				min: stat = 1
				max: stat = 2
				first: stat = 3
				maxratio: stat = 4
				"""
				stat_default = ['' for j in range(nTargetCols)]
				statrow_default = ['' for j in range(nTargetCols)]
				for j in range(nTargetCols):
						if (stat[j] == StatEnum.AVG):	# avg
								stat_default[j] = typecol[j](0.0)
						elif (stat[j] == StatEnum.MIN): # min
								stat_default[j] = float('inf')
								statrow_default[j] = -1
						elif (stat[j] == StatEnum.MAX):	# max
								stat_default[j] = float('-inf')
								statrow_default[j] = -1
						elif (stat[j] == StatEnum.FIRST):	# first
								stat_default[j] = ''
								statrow_default[j] = -1
						elif (stat[j] == StatEnum.MAXRATIO):	# maxratio
								stat_default[j] = 0.0
								statrow_default[j] = -1
						elif (is_secondary and stat[j] is None):	# no secondary stat for this col
								pass
						
				if is_secondary:
						self._secondary_statcol = [copy.deepcopy(stat_default) for i in range(nInst)]
				else:
						self._statcol = [copy.deepcopy(stat_default) for i in range(nInst)]
						self._statrow = [copy.deepcopy(statrow_default) for i in range(nInst)]

		def update_stat_from_row(self, 
																col_index, nTargetCols, inst_name, 
																params, row, rowindex, 
																stat, typecol,
																secondary_stat, secondary_col_index,
																secondary_typecol):
				""" 
				Update stat we are keeping:
				avg: stat = 0
				min: stat = 1
				max: stat = 2
				first: stat = 3
				maxratio: stat = 4
				"""
				if (len(params[0]) > 0):
						inst_index = util.index_of(row[self.inst_col], inst_name)
				else:
						inst_index = 0
				self._nRows[inst_index] += 1 

				for j in range(nTargetCols):
						try:
								if (stat[j] is not StatEnum.MAXRATIO):
										if (str(row[col_index[j]]) is 'infty'):
												curr_val = float('inf')
										else:
												curr_val = typecol[j](row[col_index[j]])
								else:
										if ('infty' in [str(row[col_index[j][0]]), str(row[col_index[j][1]])]):
												curr_val = float('inf')
										elif (abs(float(row[col_index[j][1]])) < 1e-7):
												if (abs(float(row[col_index[j][0]])) < 1e-7):
														curr_val = float(0.0)
												else:
														curr_val = float('inf')
										else:
												curr_val = (float(row[col_index[j][0]]) / float(row[col_index[j][1]]))
						except ValueError:
								if (stat[j] is not StatEnum.MAXRATIO):
										raise ValueError("Value error accessing row_{%d}[%d] with value %s." % (rowindex, col_index[j], row[col_index[j]]))
								else:
										raise ValueError("Value error accessing row_{%d}[%d] with value %s, row_{%d}[%d] with value %s." 
														% (rowindex, col_index[j][0], row[col_index[j][0]],
															 rowindex, col_index[j][1], row[col_index[j][1]]))

						check_secondary = False
						update_secondary_anyway = False
						if (stat[j] == StatEnum.AVG):	# avg
								self._statcol[inst_index][j] += curr_val
						elif (stat[j] == StatEnum.MIN):	# min
								if (curr_val < self._statcol[inst_index][j]):
										self._statcol[inst_index][j] = curr_val
										self._statrow[inst_index][j] = rowindex
										update_secondary_anyway = type(secondary_stat) is list and secondary_stat[j] is not None
								elif (curr_val == self._statcol[inst_index][j]):
										check_secondary = type(secondary_stat) is list and secondary_stat[j] is not None
						elif (stat[j] == StatEnum.MAX):	# max
								if (curr_val > self._statcol[inst_index][j]):
										self._statcol[inst_index][j] = curr_val
										self._statrow[inst_index][j] = rowindex
										update_secondary_anyway = type(secondary_stat) is list and secondary_stat[j] is not None
								elif (curr_val == self._statcol[inst_index][j]):
										check_secondary = type(secondary_stat) is list and secondary_stat[j] is not None
						elif (stat[j] == StatEnum.FIRST):	# first
								if (self._statrow[inst_index][j] < 0):
										self._statcol[inst_index][j] = curr_val
										self._statrow[inst_index][j] = rowindex
						elif (stat[j] == StatEnum.MAXRATIO):	# maxratio
								if (curr_val > self._statcol[inst_index][j]):
										self._statcol[inst_index][j] = curr_val
										self._statrow[inst_index][j] = rowindex
										update_secondary_anyway = type(secondary_stat) is list and secondary_stat[j] is not None
								elif (curr_val == self._statcol[inst_index][j]):
										check_secondary = type(secondary_stat) is list and secondary_stat[j] is not None

						if update_secondary_anyway or check_secondary:
								try:
										if (secondary_stat[j] is not StatEnum.MAXRATIO):
												if (str(row[secondary_col_index[j]]) is 'infty'):
														curr_val = float('inf')
												else:
														curr_val = secondary_typecol[j](row[secondary_col_index[j]])
										else:
												if ('infty' in [str(row[secondary_col_index[j][0]]), str(row[secondary_col_index[j][1]])]):
														curr_val = float('inf')
												elif (abs(float(row[secondary_col_index[j][1]])) < 1e-7):
														if (abs(float(row[secondary_col_index[j][0]])) < 1e-7):
																curr_val = float(0.0)
														else:
																curr_val = float('inf')
												else:
														curr_val = (float(row[secondary_col_index[j][0]]) / float(row[secondary_col_index[j][1]]))
								except ValueError:
										if (secondary_stat[j] is not StatEnum.MAXRATIO):
												raise ValueError("Value error accessing row_{%d}[%d] with value %s." % (rowindex, secondary_col_index[j], row[secondary_col_index[j]]))
										else:
												raise ValueError("Value error accessing row_{%d}[%d] with value %s, row_{%d}[%d] with value %s." 
																% (rowindex, secondary_col_index[j][0], row[secondary_col_index[j][0]],
																	 rowindex, secondary_col_index[j][1], row[secondary_col_index[j][1]]))
												

								if update_secondary_anyway:
									self._secondary_statcol[inst_index][j] = curr_val
								elif (secondary_stat[j] == StatEnum.MIN):	# min
										if (curr_val < self._secondary_statcol[inst_index][j]):
												self._secondary_statcol[inst_index][j] = curr_val
												self._statrow[inst_index][j] = rowindex
								elif (secondary_stat[j] == StatEnum.MAX):	# max
										if (curr_val > self._secondary_statcol[inst_index][j]):
												self._secondary_statcol[inst_index][j] = curr_val
												self._statrow[inst_index][j] = rowindex
								elif (secondary_stat[j] == StatEnum.MAXRATIO):	# maxratio
										if (curr_val > self._secondary_statcol[inst_index][j]):
												self._secondary_statcol[inst_index][j] = curr_val
												self._statrow[inst_index][j] = rowindex


		def row_to_stat(self, row, params):
				""" Check if row should be used for calculating stats """
				for i in range(len(params)):
						if (len(params[i]) == 0):
								continue

						val = row[i]
						if (i > 0):
								val = float(val)
						param_index = util.index_of(val, params[i])
						if (param_index < 0):
								return False

				#if (len(row) > 2):
					#ind = self._param_container.index("CUT_LIMIT")
					#assert(int(row[ind]) == 1000)
					#ind = self._param_container.index("CUT_PRESOLVE")
					#assert(int(row[ind]) == 0)
				return True

		def stat_col(self, stat = None, col_info = None, typecol = None,
										inst_name = None, pha_act_option = None, 
										num_alg2_rounds = None, num_rays_cut = None, 
										max_frac_var = None, cut_limit = None,
										use_split_share = None, 
										num_cuts_iter_bilinear = None, 
										use_unit_vectors_heur = None,
										use_cut_vert_heur = None, 
										use_tight_points_heur = None,
										cut_presolve = None, 
										print_row = None,
										secondary_stat = None, secondary_col_info = None,
										secondary_typecol = None):
				""" 
				Find avg / min / max / first / maxratio of col_info
				Stat can be a list of stats to take for each column,
				or just the index of the stat you want 
				(if it is the same for all columns)
				avg: stat = 0
				min: stat = 1
				max: stat = 2
				first: stat = 3
				maxratio: stat = 4

				Returns out_table, table with found values
				"""
				## Defaults
				if (print_row is None):
						print_row = True
				# Set default stat
				if (stat is None):
						stat = StatEnum.FIRST 
				if (type(stat) is int):
						stat = [stat]
				elif (type(stat) is not list):
						raise TypeError("Type of stat must be int or list; given %s." % type(stat))

				if __debug__:
						print( "\n## In stat_col, for stat ", end='' )
						print( StatEnum.get_name_list(stat), end='' )
						print( " ##" )
				
				# Restart reader in case we have read something else in
				self.restart_reader()

				# Set col_info information from the header
				[nTargetCols, col_index, col_name] = self.set_col(col_info, set_val = False)

				# Make sure stat has the proper number of columns
				# If one column, assume meant for all columns
				if (len(stat) is not nTargetCols):
						if (len(stat) == 1):
								stat = [stat[0] for i in range(nTargetCols)]
						else:
								raise IndexError("Need to specify %d statistics (one for each column); given %d." % (nTargetCols,len(stat)))
				
				# Double check whichever stats are "ratio" to make sure num cols is 2
				ratio_stat_problems = [j for j in range(len(stat)) 
								if stat[j] is StatEnum.MAXRATIO and (type(col_index[j]) is not list or len(col_index[j]) is not 2)]
				if (len(ratio_stat_problems) > 0):
						raise IndexError("Need to specify 2 columns for every stat that is ratio.")

				if (len([j for j in range(nTargetCols) if stat[j] is not StatEnum.FIRST]) == 0):
						stat_is_first = True
				else:
						stat_is_first = False

				# Make sure typecol is set correctly
				if (typecol is None):
						typecol = [float for j in range(nTargetCols)]
				elif (type(typecol) is type):
						typecol = [typecol] * nTargetCols
				elif (type(typecol) is not list):
						print(typecol)
						raise TypeError("Type of typecol must be type or list; given %s." % type(typecol))

				if (len(typecol) is not nTargetCols):
						if (len(typecol) == 1):
								typecol = [typecol[0] for i in range(nTargetCols)]
						else:
								raise IndexError("Need to specify %d types for columns (one for each column); given %d." % (nTargetCols,len(typecol)))
				
				# Decide if there is secondary stat and set it up
				secondary_col_index = None
				secondary_col_name = None
				if (secondary_stat is not None):
						if (type(secondary_stat) is int):
								secondary_stat = [secondary_stat]
						elif (type(secondary_stat) is not list):
								raise TypeError("Type of secondary_stat must be int or list; given %s." % type(secondary_stat))

						if (secondary_col_info is None or len(secondary_col_info) is not nTargetCols):
								raise IndexError("Number of columns specified must be the same for the secondary stats as the first.")

						# None of the secondary stats can be average or first
						# All stats that are ratio need to have corresponding two columns specified
						for s in range(len(secondary_stat)):
								if (secondary_stat[s] is StatEnum.AVG or secondary_stat[s] is StatEnum.FIRST):
										raise TypeError("Type of secondary_stat cannot be average or first")
								if (secondary_stat[s] is StatEnum.MAXRATIO):
										if (type(secondary_col_info[j]) is not list or len(secondary_col_info[j]) is not 2):
												raise IndexError("Need to specify 2 columns for every secondary stat that is ratio.")
												
						[secondary_nTargetCols, secondary_col_index, secondary_col_name] = self.set_col(secondary_col_info, set_val = False, is_secondary = True)
						if (len(secondary_stat) is not nTargetCols):
								if (len(secondary_stat) == 1):
										secondary_stat = [secondary_stat[0] for i in range(nTargetCols)]
								else:
										raise IndexError("Secondary stat needed for %d columns; given %d." % (nTargetCols,len(secondary_stat)))
						
						# Set secondary_typecol
						if (secondary_typecol is None):
								secondary_typecol = [float for j in range(nTargetCols)]
						elif (type(secondary_typecol) is type):
								secondary_typecol = [secondary_typecol] * nTargetCols
						elif (type(secondary_typecol) is not list):
								print(secondary_typecol)
								raise TypeError("Type of secondary_typecol must be type or list; given %s." % type(secondary_typecol))
 
				if __debug__:
						print( "## Finding stat of cols: " )
						for j in range(nTargetCols):
								print( "\t(stat: %s, col: %s, index: %s" 
												% (StatEnum.get_name(stat[j]), str(col_name[j]), str(col_index[j])), end='' )
								if (type(secondary_stat) is list and secondary_stat[j] is not None):
										print( "; secondary_stat: %s, secondary_col: %s, secondary_index: %s"
														% (StatEnum.get_name(secondary_stat[j]), 
																str(secondary_col_name[j]), str(secondary_col_index[j])), end='' )
								print( ")" )
						print( "\t##" )

				# Default for inst_name and other parameters
				# instance, sics, pha, 
				# pha_act_option, num_alg2_rounds, 
				# num_rays_cut, max_frac_var, 
				# cut_limit,
				# use_split_share, num_cuts_iter_bil, use_unit_vectors_heur, 
				# use_cut_vert_heur, use_tight_points_heur,
				# cut_presolve, rounds, 
				# min_orthogonality, rayeps, eps, timelimit
				params = [
						inst_name, [], [],
						pha_act_option, num_alg2_rounds, 
						num_rays_cut, max_frac_var,
						cut_limit, 
						use_split_share, num_cuts_iter_bilinear, use_unit_vectors_heur, 
						use_cut_vert_heur, use_tight_points_heur,
						cut_presolve, [], 
						[], [], [], []
				]
				assert self._param_container is not None
				assert len(params) == self._param_container.num_params
				for i in range(len(params)):
						curr_param_name = self._param_container.param_names[i]
						params[i] = self._param_container.set_param(self._param_container.param_names[i], params[i], set_val = False)

				inst_name = copy.deepcopy(params[0])
				nInst = len(params[0])
				if (nInst == 0):	# In case it is the "all instances" object
						nInst = 1

				self._nRows = [0 for i in range(nInst)]
				self._nTotalRows = 0
				
				# Initialize statistic array
				self.init_stat(stat, nTargetCols, nInst, typecol)
				if (secondary_stat is not None):
						self.init_stat(secondary_stat, nTargetCols, nInst, typecol, is_secondary = True)

				# Go through each row and pick out values that are relevant
				for row in self._reader:
						# Skip this row if it is not information from a desired instance
						if (self.row_to_stat(row, params)):
								self.update_stat_from_row(col_index, nTargetCols, inst_name, params, row, 
												self._nTotalRows, stat, typecol, secondary_stat, secondary_col_index, secondary_typecol)

								if (stat_is_first):
										params[0].remove(row[self.inst_col])
						self._nTotalRows += 1		# Do this after update, since it's +1 from row index

						# Finish if stat == FIRST and no more instances left to read
						if (stat_is_first):
								if (len(params[0]) == 0):
										params[0] = inst_name
										break
				
				# Set out_table
				out_table_header = ["instance"]
				for j in range(nTargetCols):
						out_table_header.append(StatEnum.get_name(stat[j]) + " " + str(col_name[j]))
						if (print_row and (stat[j] in [StatEnum.MIN, StatEnum.MAX, StatEnum.FIRST, StatEnum.MAXRATIO])):
								out_table_header.append("row")
				if (print_row):
						out_table_header.append("num_entries")
				out_table = [out_table_header]
				for inst in range(nInst):
						if (len(params[0]) > 0):
								entry = [params[0][inst]]
						else:
								entry = ["all_instances"]
						if __debug__:
								print( "## Instance %s: Read %d rows ##" % (entry[0], self._nRows[inst]) )
						if (self._nRows[inst] == 0):
								if __debug__:
										print( "\tNo rows read for this instance, so appending dummy entries" )
						for j in range(nTargetCols):
								if (self._nRows[inst] == 0):
										if (typecol[j] is float):
											entry.append(util.format_float(self._num_dec_places).format(0.0))
										else:
											entry.append('0')
										if (stat[j] > 0):
											if (print_row):
												entry.append(str(-1))
								elif (stat[j] == StatEnum.AVG):	# avg
										if (self._nRows[inst] > 0):
												avg = float(self._statcol[inst][j]) / self._nRows[inst]
												if __debug__:
														print( "\tAverage of column %s (index %d) is %f" % (col_name[j], col_index[j], avg) )
										else:
												avg = float('-inf')
												if __debug__:
														print( "\tNo rows read for this instance, so average cannot be computed" )
										if (typecol[j] is float):
											entry.append(util.format_float(self._num_dec_places).format(avg))
										else:
											entry.append(str(avg))
								elif (stat[j] == StatEnum.MIN):	# min
										if __debug__:
												print( 
																"\tMinimum of column %s (index %d) is %f in row %d" 
																% (col_name[j], col_index[j], self._statcol[inst][j], self._statrow[inst][j]) 
														 )
										if (typecol[j] is float):
											entry.append(util.format_float(self._num_dec_places).format(self._statcol[inst][j]))
										else:
											entry.append(str(self._statcol[inst][j]))
										if (print_row):
												entry.append(str(self._statrow[inst][j]))
								elif (stat[j] == StatEnum.MAX):	# max
										if __debug__:
												print( 
																"\tMaximum of column %s (index %d) is %f in row %d" 
																% (col_name[j], col_index[j], self._statcol[inst][j], self._statrow[inst][j]) 
														 )
										if (typecol[j] is float):
											entry.append(util.format_float(self._num_dec_places).format(self._statcol[inst][j]))
										else:
											entry.append(str(self._statcol[inst][j]))
										if (print_row):
												entry.append(str(self._statrow[inst][j]))
								elif (stat[j] == StatEnum.FIRST):	# first
										if __debug__:
												print( 
																"\tFound instance %s in row %d and filled in relevant column information" 
																% (entry[0], self._statrow[inst][j]) 
														 )
										if (typecol[j] is float):
											entry.append(util.format_float(self._num_dec_places).format(self._statcol[inst][j]))
										else:
											entry.append(str(self._statcol[inst][j]))
										if (print_row):
												entry.append(str(self._statrow[inst][j]))
								elif (stat[j] == StatEnum.MAXRATIO):	# maxratio
										if __debug__:
												print( 
																"\tMaximum ratio of column %s over column %s (index %d / index %d) is %f in row %d" 
																% (col_name[j][0], col_name[j][1], col_index[j][0], col_index[j][1], self._statcol[inst][j], self._statrow[inst][j]) 
														 )
										entry.append(util.format_float(self._num_dec_places).format(self._statcol[inst][j]))
										if (print_row):
												entry.append(str(self._statrow[inst][j]))
								
						if (print_row):
								entry.append(self._nRows[inst])
						out_table.append(entry)

				return out_table

		def avg_col(self, col_info = None, typecol = None,
						inst_name = None, pha_act_option = None, 
						num_alg2_rounds = None, num_rays_cut = None, 
						max_frac_var = None, cut_limit = None,
						use_split_share = None, 
						num_cuts_iter_bilinear = None, 
						use_unit_vectors_heur = None,
						use_cut_vert_heur = None, 
						use_tight_points_heur = None,
						cut_presolve = None, 
						print_row = None):
				""" Find avg of col_info columns """
				return self.stat_col(StatEnum.AVG, 
								col_info, typecol, inst_name, pha_act_option,
								num_alg2_rounds, num_rays_cut, 
								max_frac_var, cut_limit,
								use_split_share, num_cuts_iter_bilinear,
								use_unit_vectors_heur,
								use_cut_vert_heur, cut_presolve, print_row)

		def min_col(self, col_info = None, typecol = None,
						inst_name = None, pha_act_option = None, 
						num_alg2_rounds = None, num_rays_cut = None, 
						max_frac_var = None, cut_limit = None,
						use_split_share = None, 
						num_cuts_iter_bilinear = None, 
						use_unit_vectors_heur = None,
						use_cut_vert_heur = None, 
						use_tight_points_heur = None,
						cut_presolve = None, 
						print_row = None):
				""" Find min of col_info columns """
				return self.stat_col(StatEnum.MIN, 
								col_info, typecol, inst_name, pha_act_option,
								num_alg2_rounds, num_rays_cut, 
								max_frac_var, cut_limit,
								use_split_share, num_cuts_iter_bilinear,
								use_unit_vectors_heur,
								use_cut_vert_heur, cut_presolve, print_row)
		
		def max_col(self, col_info = None, typecol = None,
						inst_name = None, pha_act_option = None, 
						num_alg2_rounds = None, num_rays_cut = None, 
						max_frac_var = None, cut_limit = None,
						use_split_share = None, 
						num_cuts_iter_bilinear = None, 
						use_unit_vectors_heur = None,
						use_cut_vert_heur = None, 
						use_tight_points_heur = None,
						cut_presolve = None, 
						print_row = None):
				""" Find max of col_info columns """
				return self.stat_col(StatEnum.MAX, 
								col_info, typecol, inst_name, pha_act_option,
								num_alg2_rounds, num_rays_cut, 
								max_frac_var, cut_limit,
								use_split_share, num_cuts_iter_bilinear,
								use_unit_vectors_heur,
								use_cut_vert_heur, cut_presolve, print_row)
		
		def find_first_val(self, col_info = None, typecol = None,
						inst_name = None, pha_act_option = None, 
						num_alg2_rounds = None, num_rays_cut = None, 
						max_frac_var = None, cut_limit = None,
						use_split_share = None, 
						num_cuts_iter_bilinear = None, 
						use_unit_vectors_heur = None,
						use_cut_vert_heur = None, 
						use_tight_points_heur = None,
						cut_presolve = None, 
						print_row = None):
				"""
				Finds the first value of the cols in col_info for each instance in inst_name
				Returns table with these values, as well as the row in which they were found
				"""
				return self.stat_col(StatEnum.FIRST, 
								col_info, typecol, inst_name, pha_act_option,
								num_alg2_rounds, num_rays_cut, 
								max_frac_var, cut_limit,
								use_split_share, num_cuts_iter_bilinear,
								use_unit_vectors_heur,
								use_cut_vert_heur, cut_presolve, print_row)

		def maxratio_col(self, col_info = None, typecol = None,
						inst_name = None, pha_act_option = None, 
						num_alg2_rounds = None, num_rays_cut = None, 
						max_frac_var = None, cut_limit = None,
						use_split_share = None, 
						num_cuts_iter_bilinear = None, 
						use_unit_vectors_heur = None,
						use_cut_vert_heur = None, 
						use_tight_points_heur = None,
						cut_presolve = None, 
						print_row = None):
				"""
				Find maxratio of the cols in col_info
				"""
				return self.stat_col(StatEnum.MAXRATIO, 
								col_info, typecol, inst_name, pha_act_option,
								num_alg2_rounds, num_rays_cut, 
								max_frac_var, cut_limit,
								use_split_share, num_cuts_iter_bilinear,
								use_unit_vectors_heur,
								use_cut_vert_heur, cut_presolve, print_row)

		def get_entry(self, row_ind, col_ind):
				""" Find entry (currently does not catch out-of-range errors) """
				self.restart_reader()

				for i in range(row_ind):
						self._reader.next()
				
				curr_row = self._reader.next()
				if (type(col_ind) is list):
					return [curr_row[i] for i in col_ind]
				else:
					return curr_row[col_ind]

		def get_row(self, row_ind):
				""" Prints specified row (and returns the row) """
				self.restart_reader()
				for i in range(row_ind):
						self._reader.next()
				curr_row = self._reader.next()
				return curr_row

		def get_col_name(self, col_index):
				""" Returns name of column in header at index col_index """
				assert type(col_index) is int
				return self._header[self._col_name_row][col_index]

		def get_col_index(self, col_name):
				""" Returns index for which col_name occurs in header (-1 if not found) """
				assert type(col_name) is str
				return util.index_of(col_name, self._header[self._col_name_row])

		def create_table(self, row_ind, col_ind):
				""" Creates table with values from the specified rows and columns """
				self.restart_reader()
				
				# col_ind could be an int, list of ints, string, or list of strings
				# Make it list of ints
				if (type(col_ind) is list):
						for i in range(len(col_ind)):
								if (type(col_ind[i]) is str):
										col_ind[i] = self.get_col_index(col_ind[i])
				elif (type(col_ind) is str):
						col_ind = [self.get_col_index(col_ind)]
				else:
						col_ind = [col_ind]

				# row_ind might be int or str or list; make it list
				if (type(row_ind) is not list):
						row_ind = [row_ind]

				tab = []
				# Now we loop through the rows and select the ones we want
				i = 0
				for row in self._reader:
						use_row = False
						for i_tmp in row_ind:
								if (type(i_tmp) is int):
										if (i == i_tmp):
												use_row = True
												break
								elif (type(i_tmp) is str):
										if (row[0] == i_tmp):
												use_row = True
												break
						if (use_row):
								tab.append([row[j] for j in col_ind])
						i = i + 1

				return tab

		def restart_reader(self, reset_header = None):
				""" Restart reader from line 0 + num header lines, and re-read the header if requested """
				if (reset_header is None):
						reset_header = False

				if (self._reader.line_num != self._num_header_lines):
						self._f.seek(0)
						self._reader = csv.reader(self._f)
						if (reset_header):
								self.read_header()
						else:
								# Assuming the header has already been read
								# Read the appropriate number of header lines
								for i in range(self._num_header_lines):
										self._reader.next()

		def reopen_reader(self):
				""" Reopen file (assumes is closed """
				self._f = open(f_name, 'r')
				self._reader = csv.reader(self._f)

		def close_f(self):
				""" Close file """
				self._f.close()

		def __del__(self):
				""" Destructor """
				self._f.close()


### MAIN METHOD ###
if __name__ == '__main__':
		# Variables for testing purposes
		from os import environ
		HOME_DIR = environ['HOME']
		results_dir = HOME_DIR + "/repos/pha2/results/saved_results"
		f_name = results_dir + "/pha.csv"
		col_name = 'Sub3'
		inst_name = ["bm23", "vpm1"]

		print( "Hello world!" )
		inst = "bm23"
		try:
				print( "Index of %s is %d." % (inst, util.index_of(inst, inst_name)) )
		except ValueError:
				print( "Instance not found in list." )
				
		inst = "Inst1"
		try:
				print( "Index of %s is %d." % (inst, util.index_of(inst, inst_name)) )
		except ValueError:
				print( "Instance not found in list." )

		csv_reader = ProcessCSV(f_name, col_info="Sub2")
		csv_reader.set_inst(inst_name)
#		with ProcessCSV(f_name, "Sub2") as csv_reader:
#		print( csv_reader._header )
		tab = csv_reader.avg_col()
		util.print_pretty_table(tab)
		#quit()
		csv_reader.set_col(col_name)
		csv_reader.avg_col()
		csv_reader.min_col()
		csv_reader.max_col()

		csv_reader = ProcessCSV(f_name, col_info="ALL CUTS")
		csv_reader.set_inst("bm23")
		util.print_pretty_table(csv_reader.max_col(pha_act_option = 1))
		tab = csv_reader.find_first_val("OSICS")
		util.print_pretty_table(tab)

		#print("INSTANCE DEFAULT")
		#print( csv_reader._param_container.get_param("INSTANCE") )

		tab = csv_reader.stat_col(inst_name = ["bm23", "vpm1"], col_info = ["OSICS", "ALL CUTS"], stat = [StatEnum.FIRST, StatEnum.MAX])
		util.print_pretty_table(tab)
