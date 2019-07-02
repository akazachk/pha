from os import environ
PHA_DIR = environ['PHA_DIR']
ADD_TO_PATH = PHA_DIR + '/scripts/python_scripts'

from sys import path
path.append(ADD_TO_PATH)

from processcsv import ProcessCSV
from processcsv import StatEnum
import sys
import utility as util
import copy
import csv
import numpy as np
from itertools import izip
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
from scipy import stats
import shutil
import os

class LGReader(ProcessCSV):
	"""
	class LGReader
	Performs various functions related to the output file from PHA code
	"""
	default_data_dir = os.environ['PHA_DIR'] + '/data'
	default_results_dir = os.environ['PHA_DIR'] + '/results'
	default_ip_opt_fname = default_data_dir + '/' + "ip_opt.csv"
	default_in_fname = default_results_dir + '/' + "pha.csv"
	max_num_dec_places = 10
	
	# Column names
	colname_inst						= "INSTANCE"
	colname_lpobj					 = "LP OBJ"
	colname_lpopt					 = "LP OBJ"
	colname_lpbound				 = "LP OBJ"
	colname_ipopt					 = "IP OBJ"
	colname_ipobj					 = "IP OBJ"
	colname_ipbound				 = "IP OBJ"
	colname_numrows				 = "NUM ROWS"
	colname_numcols				 = "NUM COLS"
	colname_fraccore				= "FRAC CORE"
	colname_numpoints			 = "NUM POINTS (TOTAL)"
	colname_numfinalpoints	= "NUM FINAL POINTS (TOTAL)"
	colname_numsics				 = "NUM SIC"
	colname_numpha					= "NUM PHA"
	colname_sicbound				= "SIC BOUND"
	colname_phabound				= "PHA BOUND"
	colname_allcutsbound		= "ALL CUTS (WITH SIC)"
	colname_activesics			= "ACTIVE SIC"
	colname_activegics			= "ACTIVE PHA"
	colname_objtried				= "OBJ"
	colname_dupsicfails		 = "DUPLICATE_SIC_FAILS"
	colname_dupgicfails		 = "DUPLICATE_GIC_FAILS"
	colname_totaltime			 = "TOTAL_TIME"
	colname_cutfails = [
			"DUAL_CUT_SOLVER_FAILS",
			"DUPLICATE_SIC_FAILS", 
			"DUPLICATE_GIC_FAILS", 
			"ORTHOGONALITY_FAILS",
			"TIMELIMIT_FAILS", 
			"ITERATION_FAILS", 
			"ABANDONED_FAILS",
			"NON_CUTTING_FAILS", 
			"SCALING_FAILS", 
			"UNSPEC_FAILS", 
			"CUT_LIMIT_FAILS",
			"PRIMAL_CUT_SOLVER_FAILS", 
			"PRIMAL_CUT_SOLVER_NO_OBJ_FAILS",
			"NUMERICAL_ISSUES_WARNING_NO_OBJ", 
			"NUMERICAL_ISSUES_NO_OBJ"
		]
	colname_numcuts = [
			"CUT_VERTICES_CUT_HEUR",
#			"DUMMY_OBJ_CUT_HEUR",
			"ITER_BILINEAR_CUT_HEUR",
			"UNIT_VECTORS_CUT_HEUR",
			"TIGHT_POINTS_CUT_HEUR",
			"SPLIT_SHARE_CUT_HEUR"
		]
	colname_activecuts = [
			"ACTIVE CUT_VERTICES_CUT_HEUR",
			"ACTIVE DUMMY_OBJ_CUT_HEUR",
			"ACTIVE ITER_BILINEAR_CUT_HEUR",
			"ACTIVE UNIT_VECTORS_CUT_HEUR",
			"ACTIVE TIGHT_POINTS_CUT_HEUR",
			"ACTIVE SPLIT_SHARE_CUT_HEUR"
		]
	

	def __init__(self, in_fname = None, hasipval = None, fill_ip_vals = None, 
							 ip_opt_fname = None, out_fname = None, 
							 inst_name = None, col_info = None, num_header_lines = None, 
							 compute_gap_closed = True):
		""" Constructor, sets reader object and instances to all instances """
		if (in_fname is None):
			self._in_fname = copy.deepcopy(self.default_in_fname)
		else:
			self._in_fname = copy.deepcopy(in_fname)
		if __debug__:
			print( "\n## In LGReader(): Opening file %s for reading ##" % self._in_fname )

		# If requested, set IP values if they are not available 
		if (fill_ip_vals is None):
			fill_ip_vals = False
		if (hasipval is None):
			self._hasipval = True	# This is so we do not needlessly recreate the ip file
		else:
			self._hasipval = hasipval
		if (ip_opt_fname is not None):
			self._ip_opt_fname = copy.deepcopy(ip_opt_fname)
		else:
			self._ip_opt_fname = copy.deepcopy(self.default_ip_opt_fname)

		if ((not self._hasipval) and fill_ip_vals):
			self.fill_ip_opt(out_fname = self._in_fname)	# Do overwrite the input file
			self._hasipval = True
		
		super(LGReader, self).__init__(self._in_fname, inst_name, col_info, num_header_lines, compute_gap_closed)
		
		# If the IP values exist, calculate gap closed
		# Note that we check the hasipval boolean, since otherwise gap_closed()
		# will grab the IP values whether or not fill_ip_vals is True
		if ((inst_name is not None) and (self._hasipval) and (not hasattr(self, '_ip_opt'))):
			self.get_ip_opt()
			if (compute_gap_closed):
				self.gap_closed()

	def set_inst(self, inst_name = None, compute_gap_closed = True):
		""" Override set_inst from parent class so class values will be reset """
		super(LGReader, self).set_inst(inst_name, compute_gap_closed)
		if ((inst_name is not None) and (self._hasipval)):
			self._ip_opt = None
			self.get_ip_opt()
			if (compute_gap_closed):
				self.gap_closed()

	def get_ip_opt(self, ip_opt_fname = None):
		"""
		Grabs IP optimum values from ip_opt file, only for relevant instances

		Saves the values internally as self._ip_opt, a numpy array
		"""
		# TODO fix bug here that ip_opt_fname might have changed... really not a bug
		if ((hasattr(self, '_ip_opt')) and (self._ip_opt is not None)):
			return self._ip_opt
		
		if (ip_opt_fname is None):
			ip_opt_fname = copy.deepcopy(self._ip_opt_fname)

		if __debug__:
			print( "\n## Reading IP file: %s ##" % ip_opt_fname )
		
		# Read IP opt file in
		ip_opt_reader = ProcessCSV(ip_opt_fname, num_header_lines = 1)
		ip_opt_reader._num_dec_places = self.max_num_dec_places

		inst_names = super(LGReader, self).get_param(self.colname_inst)
		self._ip_opt = ['' for i in range(len(inst_names))]
		for inst in range(len(inst_names)):
			curr_inst = inst_names[inst]

			# find_first_val returns a table, with a header row
			# The first row contains all the column information
			val_str = ip_opt_reader.find_first_val(col_info = self.colname_ipobj, inst_name = curr_inst)[1][1]
			if (len(val_str) > 0):
				curr_inst_ip_obj = float(val_str)
				self._ip_opt[inst] = curr_inst_ip_obj
				if __debug__:
					print( "Instance: %s\tIP obj: %f" % (curr_inst, curr_inst_ip_obj) )
			elif __debug__:
				print( "Instance %s not found in IP file" % curr_inst )

		del ip_opt_reader

		self._ip_opt = np.asarray(self._ip_opt)
		return self._ip_opt

	def inst_info(self, 
								inst_name = None, pha_act_option = None, 
								num_alg2_rounds = None, num_rays_cut = None,
								max_frac_var = None, cut_limit = None, 
								use_split_share = None, 
								num_cuts_iter_bilinear = None, 
								use_unit_vectors_heur = None,
								use_cut_vert_heur = None, 
								cut_presolve = None):
		""" Get instance basic information (rows, columns, avg time) """
		stat = [StatEnum.FIRST, StatEnum.FIRST, StatEnum.AVG]
		col_info = [self.colname_numrows, self.colname_numcols, self.colname_totaltime]
		typecol = [int, int, float]
		if (inst_name is None):
			inst_name = super(LGReader, self).get_param(self.colname_inst)

		tab = super(LGReader, self).stat_col(stat, col_info, typecol, inst_name, pha_act_option,
																		num_alg2_rounds, num_rays_cut, max_frac_var, cut_limit,
																		use_split_share, num_cuts_iter_bilinear, use_unit_vectors_heur,
																		use_cut_vert_heur, cut_presolve)

		# Remove the columns that are not relevant
		tab = [[tab[i][j] for j in [0,1,3,5]] for i in range(len(tab))]
		tab[0] = ["Instance", "Rows", "Cols", "Time"]

		# num_rows
		col_index = 1
		self._num_rows = np.asarray([int(round(float(tab[i][col_index]))) for i in range(1,len(tab))])

		# num_cols
		col_index += 1
		self._num_cols = np.asarray([int(round(float(tab[i][col_index]))) for i in range(1,len(tab))])
		
		util.print_pretty_table(tab)
		return tab

	def get_best_row(self,
								inst_name = None, pha_act_option = None, 
								num_alg2_rounds = None, num_rays_cut = None,
								max_frac_var = None, cut_limit = None, 
								use_split_share = None, 
								num_cuts_iter_bilinear = None, 
								use_unit_vectors_heur = None,
								use_cut_vert_heur = None, 
								use_tight_points_heur = None,
								cut_presolve = None):
		"""
		Get best-performing row per instance
		"""
		stat = [StatEnum.MAX]
		col_info = [self.colname_allcutsbound] #[self.colname_activesics, self.colname_numpha, self.colname_activegics] 
		typecol = [float] #[int, int, int]
		secondary_stat = [StatEnum.MAXRATIO]
		secondary_col_info = [[self.colname_activegics, self.colname_numpha]]
		tab = super(LGReader, self).stat_col(stat, col_info, typecol, inst_name, pha_act_option,
																		num_alg2_rounds, num_rays_cut, max_frac_var, cut_limit,
																		use_split_share, num_cuts_iter_bilinear,
																		use_unit_vectors_heur,
																		use_cut_vert_heur, use_tight_points_heur, cut_presolve,
																		secondary_stat = secondary_stat,
																		secondary_col_info = secondary_col_info)
		# Get best row for each instance
		self._best_row = [int(tab[i][2]) for i in range(1,len(tab))]

	def obj_fails_table(self,
								inst_name = None, pha_act_option = None, 
								num_alg2_rounds = None, num_rays_cut = None, 
								max_frac_var = None, cut_limit = None, 
								use_split_share = None, 
								num_cuts_iter_bilinear = None, 
								use_unit_vectors_heur = None,
								use_cut_vert_heur = None, 
								cut_presolve = None):
		"""
		Analyze number of cuts generated in relation to maximum possible,
		as well as some potential common problems (e.g., duplicate SICs in the cut LP)

		Returns 2D list with the information
			[instance, num splits, num sics, num active sics, num gics, num active gics, obj tried, fails...]
		"""
		if (inst_name is None):
			inst_name = super(LGReader, self).get_param(self.colname_inst)
		num_inst = len(inst_name)
		if (not hasattr(self, '_best_row')):
				self.get_best_row(inst_name, pha_act_option,
										num_alg2_rounds, num_rays_cut, 
										max_frac_var, cut_limit,
										use_split_share, num_cuts_iter_bilinear,
										use_unit_vectors_heur,
										use_cut_vert_heur, use_tight_points_heur, cut_presolve)
		#stat = [StatEnum.MAX]
		#col_info = [self.colname_activesics, self.colname_numpha, self.colname_activegics] 
		#typecol = [int, int, int]
		#tab = super(LGReader, self).stat_col(stat, col_info, typecol, inst_name, pha_act_option,
		#																num_alg2_rounds, num_rays_cut, cut_limit,
		#																use_split_share, num_cuts_iter_bilinear,
		#																use_cut_vert_heur, cut_presolve)
		
		out_tab = [[
			"Instance", "Splits", "SICs", "Active SICs", "GICs", "Active GICs", "Obj", 
			#"Duplicate tilt", "Unbounded tilt",
			"Unbounded", "Dup SIC", "Dup GIC", "Orth", "Time limit", "Iter limit", "Abandoned",
			"Non-cutting", "Scaling", "Unspec", "Cut limit", "Primal infeas", "Primal infeas (setup)",
			"Numerics warnings", "Numerics errors"]]
		outcol_obj_tried = out_tab[0].index("Obj");
		numcols_out_tab = len(out_tab[0]);
		
		np_out_tab = np.zeros(shape = (num_inst, numcols_out_tab), dtype=int)

		# Get column indices for each of the relevant stats
		index_numsplits = super(LGReader, self).get_col_index(self.colname_fraccore)
		index_numsics = super(LGReader, self).get_col_index(self.colname_numsics)
		index_activesics = super(LGReader, self).get_col_index(self.colname_activesics)
		index_numpha = super(LGReader, self).get_col_index(self.colname_numpha)
		index_activegics = super(LGReader, self).get_col_index(self.colname_activegics)
		index_objtried = super(LGReader, self).get_col_index(self.colname_objtried)
		index_cutfail = [-1 for i in range(len(self.colname_cutfails))]
		for i in range(len(self.colname_cutfails)):
			index_cutfail[i] = super(LGReader, self).get_col_index(self.colname_cutfails[i])

		for i in range(len(np_out_tab)):
			if __debug__:
				print( "## Obj_fails_table: Filling in information for instance %d with name %s ##" % (i,inst_name[i]) )
			if self._best_row[i] >= 0:
				row = super(LGReader, self).get_row(self._best_row[i])
				np_out_tab[i,1] = int(row[index_numsplits])
				np_out_tab[i,2] = int(row[index_numsics])
				np_out_tab[i,3] = int(row[index_activesics])
				np_out_tab[i,4] = int(row[index_numpha])
				np_out_tab[i,5] = int(row[index_activegics])
				obj_tried = int(row[index_objtried])
				np_out_tab[i,outcol_obj_tried] = obj_tried
				for j in range(len(index_cutfail)):
					curr_fail = int(row[index_cutfail[j]])
					np_out_tab[i,outcol_obj_tried+1+j] = 100. * curr_fail / obj_tried
			else:
				for j in range(len(index_cutfail)):
					np_out_tab[i,outcol_obj_tried+1+j] = 0
				
			
		out_tab.extend(np_out_tab.tolist())
		for i in range(1,num_inst+1):
			out_tab[i][0] = inst_name[i-1]

		util.print_pretty_table(out_tab)
		return out_tab

	def active_cuts_table(self,
								inst_name = None, pha_act_option = None, 
								num_alg2_rounds = None, num_rays_cut = None, 
								max_frac_var = None, cut_limit = None, 
								use_split_share = None, 
								num_cuts_iter_bilinear = None, 
								use_unit_vectors_heur = None,
								use_cut_vert_heur = None, 
								use_tight_points_heur = None,
								cut_presolve = None):
		"""
		Analyze which heuristics led to most active cuts 

		Returns 2D list with the information
			[instance, num sics, num active sics, num gics, num active gics, active from which heur...]
		"""
		if (inst_name is None):
			inst_name = super(LGReader, self).get_param(self.colname_inst)
		num_inst = len(inst_name)
		if (not hasattr(self, '_best_row')):
				self.get_best_row(inst_name, pha_act_option,
										num_alg2_rounds, num_rays_cut, max_frac_var, cut_limit,
										use_split_share, num_cuts_iter_bilinear, use_unit_vectors_heur,
										use_cut_vert_heur, use_tight_points_heur, cut_presolve)
		
		out_tab = [[
					"Instance", "SICs", "Active SICs", "GICs", "Active GICs",
					"V", "Active V", #
					"B", "Active B",
					"R", "Active R",
					"T", "Active T",
					"S", "Active S"
			]]
		numcols_out_tab = len(out_tab[0]);
		
		np_out_tab = np.zeros(shape = (num_inst, numcols_out_tab), dtype=int)
		
		# Get column indices for number of each type of cut
		# Number active of that type is the subsequent column (be aware could change)
		colindex_first = [super(LGReader, self).get_col_index(self.colname_numsics),
								super(LGReader, self).get_col_index(self.colname_activesics),
								super(LGReader, self).get_col_index(self.colname_numpha),
								super(LGReader, self).get_col_index(self.colname_activegics)]
		colindex_numcuts = [super(LGReader, self).get_col_index(self.colname_numcuts[i]) for i in range(len(self.colname_numcuts))]
		for i in range(len(np_out_tab)):
				if __debug__:
						print( "## Active_cuts_table: Filling in information for instance %d with name %s ##" % (i,inst_name[i]) )
				if self._best_row[i] >= 0:
					curr_row = super(LGReader, self).get_row(self._best_row[i])
					for j in range(len(colindex_first)):
							np_out_tab[i,j+1] = curr_row[colindex_first[j]]
					for j in range(len(colindex_numcuts)):
							np_out_tab[i,2*j+1+len(colindex_first)] = curr_row[colindex_numcuts[j]]
							np_out_tab[i,2*j+1+len(colindex_first)+1] = curr_row[colindex_numcuts[j]+1]
				else:
					for j in range(len(colindex_first)):
							np_out_tab[i,j+1] =	0
					for j in range(len(colindex_numcuts)):
							np_out_tab[i,2*j+1+len(colindex_first)] = 0
							np_out_tab[i,2*j+1+len(colindex_first)+1] = 0
					
		
		out_tab.extend(np_out_tab.tolist())
		for i in range(1,num_inst+1):
			out_tab[i][0] = inst_name[i-1]

		util.print_pretty_table(out_tab)
		return out_tab

	def write_best_params(self, out_dir = None):
		""" Writes best parameters for each instance to file """
		if (out_dir is None):
			out_dir = self.default_data_dir + "/params"
			try:
				os.makedirs(out_dir)
			except OSError:	# for the race condition, however unlikely
				if not os.path.isdir(out_dir):
					raise
		
		#num_params = super(LGReader, self)._param_container.num_params
		num_params = self._param_container.num_params
		inst_names = super(LGReader, self).get_param(self.colname_inst)
		if (not hasattr(self, '_best_row')):
				self.get_best_row(inst_names, cut_limit=[1000], cut_presolve=[0])  # these are current defaults we want to report for, but this can be confusing if we forget we set it...

		out_tab = [copy.deepcopy(self._header[0])]
		out_tab.append(copy.deepcopy(self._header[1]))
		for i in range(len(inst_names)):
			inst = inst_names[i]
			row = super(LGReader, self).get_row(self._best_row[i]) if self._best_row[i] >= 0 else [0 for i in range(num_params)]
			out_tab.append(row)

			curr_fname = out_dir + '/' + inst + "_params.txt"
			with open(curr_fname, 'wb') as out_f:
				if __debug__:
					print( "## Writing parameters for %s ##" % inst )
				for p in range(1,num_params):
					out_f.write(str(self._param_container.param_names[p]).lower() + ' ')
					curr_val = self._param_container.type_param[p](row[p])
					out_f.write(str(curr_val) + '\n')

		# Save parameter information
		with open(out_dir + '/' + "best_runs.csv", 'wb') as out_f:
			out_writer = csv.writer(out_f)
			out_writer.writerows(out_tab)
			del out_writer

	def gap_closed(self, inst_name = None, pha_act_option = None, 
							 num_alg2_rounds = None, num_rays_cut = None, 
							 max_frac_var = None, cut_limit = None, 
							 use_split_share = None, 
							 num_cuts_iter_bilinear = None, 
							 use_unit_vectors_heur = None,
							 use_cut_vert_heur = None, 
							 use_tight_points_heur = None,
							 cut_presolve = None, recompute_best_row = None):
		"""
		Adds gap closed information to the instance of LGReader
		In addition, keeps lp opt, ip opt, osic best, and all cuts best

		Defaults to all instances, unless inst_name is specified
		"""
		# Make sure that ip values are available 
		# (i.e., reread the file if necessary with ip values filled)
		if (not self._hasipval):
			self.fill_ip_opt(out_fname = self._in_fname)
			self._hasipval = True

		# Get best_row
		if (not hasattr(self, '_best_row')) or recompute_best_row:
				self.get_best_row(inst_name, pha_act_option,
										num_alg2_rounds, num_rays_cut, max_frac_var, cut_limit,
										use_split_share, num_cuts_iter_bilinear, use_unit_vectors_heur,
										use_cut_vert_heur, use_tight_points_heur, cut_presolve)

		# Set defaults
		col_info = [self.colname_lpobj, self.colname_sicbound, self.colname_allcutsbound]
		stat = [StatEnum.FIRST, StatEnum.MAX, StatEnum.MAX]
		typecol = [float, float, float]
		if (inst_name is None):
			inst_name = super(LGReader, self).get_param(self.colname_inst)

		# Save the current number of output decimal places, and set current outputted decimal places
		saved_num_dec_places = self._num_dec_places
		self._num_dec_places = self.max_num_dec_places

		if __debug__:
			print( "\n## Calculating gap_closed ##" )
		
		# tab will have columns 
		# [inst, lp opt, row, osic opt, row, all cut opt, row, num rows]
		tab = super(LGReader, self).stat_col(stat, col_info, typecol, inst_name, pha_act_option,
																		num_alg2_rounds, num_rays_cut, max_frac_var, cut_limit, 
																		use_split_share, num_cuts_iter_bilinear, use_unit_vectors_heur, 
																		use_cut_vert_heur, use_tight_points_heur, cut_presolve)
		self._num_dec_places = saved_num_dec_places
		#self._best_row = [int(tab[i][6]) for i in range(1,len(tab))]

		# ip opt
		if (not hasattr(self, '_ip_opt')):
			self.get_ip_opt()

		# Add information from instances
		lp_col_index = super(LGReader, self).get_col_index(self.colname_lpobj)
		sics_col_index = super(LGReader, self).get_col_index(self.colname_sicbound)
		allcuts_col_index = super(LGReader, self).get_col_index(self.colname_allcutsbound)
		self._lp_opt = []
		self._sic_opt= []
		self._sics_gap = []
		self._allcuts_opt = []
		self._allcuts_gap = []
		self._gap_closed = []
		for i in range(len(inst_name)):
				if self._best_row[i] >= 0:
					curr_row = super(LGReader, self).get_row(self._best_row[i])
					self._lp_opt.append(float(curr_row[lp_col_index]))
					self._sic_opt.append(float(curr_row[sics_col_index]))
					self._sics_gap.append(100 * (self._sic_opt[i] - self._lp_opt[i]) / (self._ip_opt[i] - self._lp_opt[i]))
					self._allcuts_opt.append(float(curr_row[allcuts_col_index]))
					self._allcuts_gap.append(100 * (self._allcuts_opt[i] - self._lp_opt[i]) / (self._ip_opt[i] - self._lp_opt[i]))
					self._gap_closed.append(100 * (self._allcuts_opt[i] - self._sic_opt[i]) / (self._ip_opt[i] - self._lp_opt[i]))
				else:
					self._lp_opt.append(0.0)
					self._sic_opt.append(0.0)
					self._sics_gap.append(0.0)
					self._allcuts_opt.append(0.0)
					self._allcuts_gap.append(0.0)
					self._gap_closed.append(0.0)

		return self._gap_closed

	def gap_closed_table(self, 
							 inst_name = None, pha_act_option = None,
							 num_alg2_rounds = None, num_rays_cut = None, 
							 max_frac_var = None, cut_limit = None, 
							 use_split_share = None, 
							 num_cuts_iter_bilinear = None, 
							 use_unit_vectors_heur = None,
							 use_cut_vert_heur = None, 
							 use_tight_points_heur = None,
							 cut_presolve = None,
							 recompute = False):
		""" Create table with gap closed information """
		if (not hasattr(self, '_gap_closed')) or recompute:
			self.gap_closed(inst_name, pha_act_option, num_alg2_rounds,
				num_rays_cut, max_frac_var, cut_limit, 
				use_split_share, num_cuts_iter_bilinear, use_unit_vectors_heur, use_cut_vert_heur,
				use_tight_points_heur, cut_presolve, recompute)
		if (not hasattr(self, '_num_rows') or not hasattr(self, '_num_cols')):
			self.inst_info()
		if (not hasattr(self, '_best_row')):
				self.get_best_row(inst_name, pha_act_option,
										num_alg2_rounds, num_rays_cut, max_frac_var, cut_limit,
										use_split_share, num_cuts_iter_bilinear, use_unit_vectors_heur,
										use_cut_vert_heur, use_tight_points_heur, cut_presolve)
		append_average = True
		eps = 1e-5

		nonzero_gap_indices = [i for i in range(len(self._gap_closed)) if (self._allcuts_gap[i] > eps)]
		zero_gap_indices = [i for i in range(len(self._gap_closed)) if (self._allcuts_gap[i] <= eps)]

		# Set up transpose of out_table
		out_tab_tr = [[super(LGReader, self).get_param(self.colname_inst)[i] for i in nonzero_gap_indices]]

		out_tab_tr.append(np.around([self._num_rows[i] for i in nonzero_gap_indices], decimals=self._num_dec_places).tolist())
		out_tab_tr.append(np.around([self._num_cols[i] for i in nonzero_gap_indices], decimals=self._num_dec_places).tolist())
		out_tab_tr.append(np.around([self._lp_opt[i] for i in nonzero_gap_indices], decimals=self._num_dec_places).tolist())
		out_tab_tr.append(np.around([self._ip_opt[i] for i in nonzero_gap_indices], decimals=self._num_dec_places).tolist())
		out_tab_tr.append(np.around([self._sic_opt[i] for i in nonzero_gap_indices], decimals=self._num_dec_places).tolist())
		out_tab_tr.append(np.around([self._allcuts_opt[i] for i in nonzero_gap_indices], decimals=self._num_dec_places).tolist())
		out_tab_tr.append(np.around([self._sics_gap[i] for i in nonzero_gap_indices], decimals=self._num_dec_places).tolist())
		out_tab_tr.append(np.around([self._allcuts_gap[i] for i in nonzero_gap_indices], decimals=self._num_dec_places).tolist())
		out_tab_tr.append(np.around([self._gap_closed[i] for i in nonzero_gap_indices], decimals=self._num_dec_places).tolist())

		# Also add num SICs, num active SICs, num GICs, num active GICs
		numsics_col_index = super(LGReader, self).get_col_index(self.colname_numsics)
		activesics_col_index = super(LGReader, self).get_col_index(self.colname_activesics)
		numpha_col_index = super(LGReader, self).get_col_index(self.colname_numpha)
		activegics_col_index = super(LGReader, self).get_col_index(self.colname_activegics)
		num_sics_tab = []
		num_active_sics_tab = []
		num_gics_tab = []
		num_active_gics_tab = []
		percent_active_tab = []
		for i in nonzero_gap_indices:
				if self._best_row[i] >= 0:
					curr_row = super(LGReader, self).get_row(self._best_row[i])
					num_sics_tab.append(int(curr_row[numsics_col_index]))
					num_active_sics_tab.append(int(curr_row[activesics_col_index]))
					num_gics_tab.append(int(curr_row[numpha_col_index]))
					num_active_gics_tab.append(int(curr_row[activegics_col_index]))
					num_pha = float(curr_row[numpha_col_index])
					percent_active_tab.append(
						(100. * float(curr_row[activegics_col_index]) / num_pha) if num_pha > 0 else 0)
				else:
					num_sics_tab.append(0)
					num_active_sics_tab.append(0)
					num_gics_tab.append(0)
					num_active_gics_tab.append(0)
					percent_active_tab.append(0.0)
					

		out_tab_tr.append(num_sics_tab) # num SICs
		out_tab_tr.append(num_active_sics_tab) # active SICs
		out_tab_tr.append(num_gics_tab) # num GICs
		out_tab_tr.append(num_active_gics_tab) # active GICs
		out_tab_tr.append(percent_active_tab) # % active

		# Header
		out_tab = [
						[
								'', '', '',
								"Opt", "Opt", "Opt", "Opt",
								"Best % gap closed", "Best % gap closed", "Best % gap closed",
								"# cuts", "# cuts", "# cuts", "# cuts",
								''
						],
						[
								"Instance", "Rows", "Cols",
								"LP", "IP", "SIC", "GIC+SIC", 
								"SIC", "GIC", "Diff",
								"SICs", "Active SICs", "GICs", "Active GICs", 
								"% active"
						]
				]
		out_tab.extend([list(t) for t in izip(*out_tab_tr)])

		if (append_average):
			out_tab.append(
				[
				"Average",'','','','','','',
				np.around(np.mean([self._sics_gap[i] for i in nonzero_gap_indices], dtype=np.float64), decimals=self._num_dec_places),
				np.around(np.mean([self._allcuts_gap[i] for i in nonzero_gap_indices], dtype=np.float64), decimals=self._num_dec_places),
				np.around(np.mean([self._gap_closed[i] for i in nonzero_gap_indices], dtype=np.float64), decimals=self._num_dec_places),
				'','','','',
				np.around(np.mean(percent_active_tab, dtype=np.float64), decimals=self._num_dec_places)
				]
			)
		
		out_tab_tr_zero = [[super(LGReader, self).get_param(self.colname_inst)[i] for i in zero_gap_indices]]
		out_tab_tr_zero.append(np.around([self._num_rows[i] for i in zero_gap_indices], decimals=self._num_dec_places).tolist())
		out_tab_tr_zero.append(np.around([self._num_cols[i] for i in zero_gap_indices], decimals=self._num_dec_places).tolist())
		out_tab_tr_zero.append(np.around([self._lp_opt[i] for i in zero_gap_indices], decimals=self._num_dec_places).tolist())
		out_tab_tr_zero.append(np.around([self._ip_opt[i] for i in zero_gap_indices], decimals=self._num_dec_places).tolist())
		out_tab_tr_zero.append(np.around([self._sic_opt[i] for i in zero_gap_indices], decimals=self._num_dec_places).tolist())
		out_tab_tr_zero.append(np.around([self._allcuts_opt[i] for i in zero_gap_indices], decimals=self._num_dec_places).tolist())
		out_tab_tr_zero.append(np.around([self._sics_gap[i] for i in zero_gap_indices], decimals=self._num_dec_places).tolist())
		out_tab_tr_zero.append(np.around([self._allcuts_gap[i] for i in zero_gap_indices], decimals=self._num_dec_places).tolist())
		out_tab_tr_zero.append(np.around([self._gap_closed[i] for i in zero_gap_indices], decimals=self._num_dec_places).tolist())
		num_sics_tab = []
		num_active_sics_tab = []
		num_gics_tab = []
		num_active_gics_tab = []
		percent_active_tab = []
		for i in zero_gap_indices:
				if self._best_row[i] >= 0:
					curr_row = super(LGReader, self).get_row(self._best_row[i])
					num_sics_tab.append(int(curr_row[numsics_col_index]))
					num_active_sics_tab.append(int(curr_row[activesics_col_index]))
					num_gics_tab.append(int(curr_row[numpha_col_index]))
					num_active_gics_tab.append(int(curr_row[activegics_col_index]))
					num_pha = float(curr_row[numpha_col_index])
					percent_active_tab.append(
						(100. * float(curr_row[activegics_col_index]) / num_pha) if num_pha > 0 else 0)
				else:
					num_sics_tab.append(0)
					num_active_sics_tab.append(0)
					num_gics_tab.append(0)
					num_active_gics_tab.append(0)
					percent_active_tab.append(0.0)
					

		out_tab_tr_zero.append(num_sics_tab) # num SICs
		out_tab_tr_zero.append(num_active_sics_tab) # active SICs
		out_tab_tr_zero.append(num_gics_tab) # num GICs
		out_tab_tr_zero.append(num_active_gics_tab) # active GICs
		out_tab_tr_zero.append(percent_active_tab) # % active

		out_tab.extend([list(t) for t in izip(*out_tab_tr_zero)])
		
		util.print_pretty_table(out_tab)
		return out_tab

	def hplane_analysis(self, hh_start = 0, hh_end = 2, num_act_start = 0, num_act_end = 4):
		"""
		Sees effect of the various hplane selection heuristics,
		as well as the number of activated hyperplanes

		Outputs 2D list with the data
		"""
		if (not hasattr(self, '_gap_closed')):
			self.gap_closed()

		col_name = [self.colname_allcutsbound]
		stat = StatEnum.MAX
		make_int = False
		
		inst_names = super(LGReader, self).get_param(self.colname_inst)

		numpoints_col_index = super(LGReader, self).get_col_index(self.colname_numpoints)
		numfinalpoints_col_index = super(LGReader, self).get_col_index(self.colname_numfinalpoints)
		points_tab_tr = []
		finalpoints_tab_tr = []

		out_tab_tr = [super(LGReader, self).get_param(self.colname_inst)]
		out_tab_tr.append(np.around(self._sics_gap, self._num_dec_places).tolist())
		out_tab_tr.append(np.around(self._allcuts_gap, self._num_dec_places).tolist())
		out_tab_tr.append(np.around(self._gap_closed, self._num_dec_places).tolist())

		append_average = True
		if (append_average):
			out_tab_tr[0].append("Average")
			out_tab_tr[1].append(np.around(np.mean(self._sics_gap, dtype=np.float64), self._num_dec_places))
			out_tab_tr[2].append(np.around(np.mean(self._allcuts_gap, dtype=np.float64), self._num_dec_places))
			out_tab_tr[3].append(np.around(np.mean(self._gap_closed, dtype=np.float64), self._num_dec_places))

		saved_num_dec_places = self._num_dec_places

		# First, best results for Alg 3 on its own
		for hh in range(hh_start,hh_end+1):
			self._num_dec_places = self.max_num_dec_places
			tab = super(LGReader, self).stat_col(col_info = col_name, pha_act_option = hh, stat = stat, num_alg2_rounds = 0)
			self._num_dec_places = saved_num_dec_places
			for c in range(len(col_name)):
				if (stat in [StatEnum.MIN, StatEnum.MAX, StatEnum.FIRST]):
					tab_col_index = 1 + 2 * c	# tab cols are [inst, stat, row, stat, row, ..., num_inst]
				elif (stat in [StatEnum.AVG]):
					tab_col_index = 1 + c	# tab cols are [inst, stat, stat, ..., num_inst]
				
				curr_gap = self.helper_col_analysis(col_name, stat, 
																						tab = tab, tab_col_index = tab_col_index,
																						append_average = append_average) 
				out_tab_tr.append(np.around(curr_gap, self._num_dec_places).tolist())

				curr_col = []
				for inst in range(len(inst_names)):
					curr_row = int(tab[inst+1][2])
					curr_col.append(super(LGReader, self).get_entry(curr_row, [numpoints_col_index, numfinalpoints_col_index]))
					
				points_tab_tr.append([int(curr_col[i][0]) for i in range(len(inst_names))])
				finalpoints_tab_tr.append([int(curr_col[i][1]) for i in range(len(inst_names))])
				
		# Now, best results for Alg 3 w/Alg 2, best among all hyperplane activation choices
		for num_act in range(num_act_start,num_act_end+1):
			self._num_dec_places = self.max_num_dec_places
			tab = super(LGReader, self).stat_col(col_info = col_name, num_alg2_rounds = num_act, stat = stat)
			self._num_dec_places = saved_num_dec_places
			for c in range(len(col_name)):
				if (stat in [StatEnum.MIN, StatEnum.MAX, StatEnum.FIRST]):
					tab_col_index = 1 + 2 * c	# tab cols are [inst, stat, row, stat, row, ..., num_inst]
				elif (stat in [StatEnum.AVG]):
					tab_col_index = 1 + c	# tab cols are [inst, stat, stat, ..., num_inst]
				
				curr_gap = self.helper_col_analysis(col_name, stat, 
																						tab = tab, tab_col_index = tab_col_index,
																						append_average = append_average) 
				out_tab_tr.append(np.around(curr_gap, self._num_dec_places).tolist())
				
				curr_col = []
				for inst in range(len(inst_names)):
					curr_row = int(tab[inst+1][2])
					curr_col.append(super(LGReader, self).get_entry(curr_row, [numpoints_col_index, numfinalpoints_col_index]))
					
				points_tab_tr.append([int(curr_col[i][0]) for i in range(len(inst_names))])
				finalpoints_tab_tr.append([int(curr_col[i][1]) for i in range(len(inst_names))])

		# Lastly, best results for Alg 2 on its own
		for num_act in range(num_act_start+1,num_act_end+1): # +1 because we do not want 0 alg 2 activations
			self._num_dec_places = self.max_num_dec_places
			tab = super(LGReader, self).stat_col(col_info = col_name, num_alg2_rounds = -1*num_act, stat = stat)
			self._num_dec_places = saved_num_dec_places
			for c in range(len(col_name)):
				if (stat in [StatEnum.MIN, StatEnum.MAX, StatEnum.FIRST]):
					tab_col_index = 1 + 2 * c	# tab cols are [inst, stat, row, stat, row, ..., num_inst]
				elif (stat in [StatEnum.AVG]):
					tab_col_index = 1 + c	# tab cols are [inst, stat, stat, ..., num_inst]
				
				curr_gap = self.helper_col_analysis(col_name, stat, 
																						tab = tab, tab_col_index = tab_col_index,
																						append_average = append_average) 
				out_tab_tr.append(np.around(curr_gap, self._num_dec_places).tolist())
				
				curr_col = []
				for inst in range(len(inst_names)):
					curr_row = int(tab[inst+1][2])
					curr_col.append(super(LGReader, self).get_entry(curr_row, [numpoints_col_index, numfinalpoints_col_index]))
					
				points_tab_tr.append([int(curr_col[i][0]) for i in range(len(inst_names))])
				finalpoints_tab_tr.append([int(curr_col[i][1]) for i in range(len(inst_names))])
		
		#out_tab_tr.extend(points_tab_tr)
		#out_tab_tr.extend(finalpoints_tab_tr)
		
		out_tab = [
				[ '' ]
				+ 3 * [ '' ]
				+ 3 * [ 'Only Alg. 3' ]
				+ 5 * [ 'Alg. 3 with Alg. 2' ]
				+ 4 * [ 'Only Alg. 2' ],
#				+ (hh_end-hh_start+1 + num_act_end-num_act_start+1 + 4) * ['Points']
#				+ (hh_end-hh_start+1 + num_act_end-num_act_start+1 + 4) * ['Final Points'],
				["Instance", "SIC", "Best", "Diff"] 
				+ ["H"+str(hh) for hh in range(hh_start+1,hh_end+2)] 
				+ [ "T+"+str(act) for act in range(num_act_start,num_act_end+1)]
				+ [ "+"+str(act) for act in range(1,4+1)]
#				+ ["H"+str(hh) for hh in range(hh_start+1,hh_end+2)] 
#				+ [ "T+"+str(act) for act in range(num_act_start,num_act_end+1)]
#				+ [ "+"+str(act) for act in range(1,4+1)]
#				+ ["H"+str(hh) for hh in range(hh_start+1,hh_end+2)] 
#				+ [ "T+"+str(act) for act in range(num_act_start,num_act_end+1)]
#				+ [ "+"+str(act) for act in range(1,4+1)]
		]

		out_tab.extend([list(x) for x in izip(*out_tab_tr)])

		util.print_pretty_table(out_tab)
		return out_tab
	
	def hplane_point_analysis(self, hh_start = 0, hh_end = 2, num_act_start = 0, num_act_end = 3):
		"""
		Sees effect of the various hplane selection heuristics,
		as well as the number of activated hyperplanes,
		on the number of points and final points

		Outputs 2D list with the data
		"""
		if (not hasattr(self, '_gap_closed')):
			self.gap_closed()

		col_name = [self.colname_allcutsbound]
		stat = StatEnum.MAX
		make_int = False

		numpoints_col_index = super(LGReader, self).get_col_index(self.colname_numpoints)
		numfinalpoints_col_index = super(LGReader, self).get_col_index(self.colname_numfinalpoints)

		inst_names = super(LGReader, self).get_param(self.colname_inst)
		out_tab_tr = [inst_names]
		points_tab_tr = []
		finalpoints_tab_tr = []

		saved_num_dec_places = self._num_dec_places
		for hh in range(hh_start,hh_end+1):
			self._num_dec_places = self.max_num_dec_places
			tab = super(LGReader, self).stat_col(col_info = col_name, pha_act_option = hh, stat = stat, num_alg2_rounds = 0)
			self._num_dec_places = saved_num_dec_places
			
			curr_col = []
			for inst in range(len(inst_names)):
				curr_row = int(tab[inst+1][2])
				curr_col.append(super(LGReader, self).get_entry(curr_row, [numpoints_col_index, numfinalpoints_col_index]))
			
			points_tab_tr.append([int(curr_col[i][0]) for i in range(len(inst_names))])
			finalpoints_tab_tr.append([int(curr_col[i][1]) for i in range(len(inst_names))])
				
		for num_act in range(num_act_start,num_act_end+1):
			self._num_dec_places = self.max_num_dec_places
			tab = super(LGReader, self).stat_col(col_info = col_name, num_alg2_rounds = num_act, stat = stat)
			self._num_dec_places = saved_num_dec_places
			
			curr_col = []
			for inst in range(len(inst_names)):
				curr_row = int(tab[inst+1][2])
				curr_col.append(super(LGReader, self).get_entry(curr_row, [numpoints_col_index, numfinalpoints_col_index]))
			
			points_tab_tr.append([int(curr_col[i][0]) for i in range(len(inst_names))])
			finalpoints_tab_tr.append([int(curr_col[i][1]) for i in range(len(inst_names))])

		for num_act in range(1,4+1):
			self._num_dec_places = self.max_num_dec_places
			tab = super(LGReader, self).stat_col(col_info = col_name, num_alg2_rounds = -1*num_act, stat = stat)
			self._num_dec_places = saved_num_dec_places
			
			curr_col = []
			for inst in range(len(inst_names)):
				curr_row = int(tab[inst+1][2])
				curr_col.append(super(LGReader, self).get_entry(curr_row, [numpoints_col_index, numfinalpoints_col_index]))
			
			points_tab_tr.append([int(curr_col[i][0]) for i in range(len(inst_names))])
			finalpoints_tab_tr.append([int(curr_col[i][1]) for i in range(len(inst_names))])

		out_tab_tr.extend(points_tab_tr)
		out_tab_tr.extend(finalpoints_tab_tr)
		
		out_tab = [
				[ '' ]
				+ (hh_end-hh_start+1 + num_act_end-num_act_start+1 + 4) * ['Points']
				+ (hh_end-hh_start+1 + num_act_end-num_act_start+1 + 4) * ['Final Points'],
				["Instance"] 
				+ ["HH"+str(hh) for hh in range(hh_start+1,hh_end+2)] 
				+ [ "+"+str(act)+"H" for act in range(num_act_start,num_act_end+1)]
				+ [ "-"+str(act)+"H" for act in range(1,4+1)]
				+ ["HH"+str(hh) for hh in range(hh_start+1,hh_end+2)] 
				+ [ "+"+str(act)+"H" for act in range(num_act_start,num_act_end+1)]
				+ [ "-"+str(act)+"H" for act in range(1,4+1)]]
		out_tab.extend([list(x) for x in izip(*out_tab_tr)])

		util.print_pretty_table(out_tab)
		return out_tab

	def cut_heur_analysis(self):
		"""
		Analyzes effect of different cut selection heuristics that were used

		Outputs 2D list with the data
		"""
		col_name = [self.colname_allcutsbound]
		stat = StatEnum.MAX
		out_tab_tr = [super(LGReader, self).get_param(self.colname_inst)]
		
		append_average = True
		if (append_average):
			out_tab_tr[0].append("Average")
		
		saved_num_dec_places = self._num_dec_places

		# Best
		out_tab_tr.append(np.around(self._allcuts_gap, self._num_dec_places).tolist())
		if (append_average):
			out_tab_tr[1].append(np.around(np.mean(self._allcuts_gap, dtype=np.float64), self._num_dec_places))

		# R
		self._num_dec_places = self.max_num_dec_places
		curr_gap = self.helper_col_analysis(col_name, stat, append_average = append_average,
						use_unit_vectors_heur = [], use_split_share = 0, num_cuts_iter_bilinear = 0, use_cut_vert_heur = 0, use_tight_points_heur = 0)
		self._num_dec_places = saved_num_dec_places
		out_tab_tr.append(np.around(curr_gap, self._num_dec_places).tolist())

		# S
		self._num_dec_places = self.max_num_dec_places
		curr_gap = self.helper_col_analysis(col_name, stat, append_average = append_average,
						use_unit_vectors_heur = 0, use_split_share = [], num_cuts_iter_bilinear = 0, use_cut_vert_heur = 0, use_tight_points_heur = 0)
		self._num_dec_places = saved_num_dec_places
		out_tab_tr.append(np.around(curr_gap, self._num_dec_places).tolist())

		# T
		self._num_dec_places = self.max_num_dec_places
		curr_gap = self.helper_col_analysis(col_name, stat, append_average = append_average,
						use_unit_vectors_heur = 0, use_split_share = 0, num_cuts_iter_bilinear = 0, use_cut_vert_heur = 0, use_tight_points_heur = [])
		self._num_dec_places = saved_num_dec_places
		out_tab_tr.append(np.around(curr_gap, self._num_dec_places).tolist())

		# V
		self._num_dec_places = self.max_num_dec_places
		curr_gap = self.helper_col_analysis(col_name, stat, append_average = append_average,
						use_unit_vectors_heur = 0, use_split_share = 0, num_cuts_iter_bilinear = 0, use_cut_vert_heur = [], use_tight_points_heur = 0)
		self._num_dec_places = saved_num_dec_places
		out_tab_tr.append(np.around(curr_gap, self._num_dec_places).tolist())

		if False:
				# B
				self._num_dec_places = self.max_num_dec_places
				curr_gap = self.helper_col_analysis(col_name, stat, append_average = append_average,
								use_unit_vectors_heur = 0, use_split_share = 0, num_cuts_iter_bilinear = [], use_cut_vert_heur = 0, use_tight_points_heur = 0)
				self._num_dec_places = saved_num_dec_places
				out_tab_tr.append(np.around(curr_gap, self._num_dec_places).tolist())

				# B+R
				self._num_dec_places = self.max_num_dec_places
				curr_gap = self.helper_col_analysis(col_name, stat, append_average = append_average,
								use_unit_vectors_heur = [], use_split_share = 0, num_cuts_iter_bilinear = [], use_cut_vert_heur = 0, use_tight_points_heur = 0)
				self._num_dec_places = saved_num_dec_places
				out_tab_tr.append(np.around(curr_gap, self._num_dec_places).tolist())

				# B+S
				self._num_dec_places = self.max_num_dec_places
				curr_gap = self.helper_col_analysis(col_name, stat, append_average = append_average,
								use_unit_vectors_heur = 0, use_split_share = [], num_cuts_iter_bilinear = [], use_cut_vert_heur = 0, use_tight_points_heur = 0)
				self._num_dec_places = saved_num_dec_places
				out_tab_tr.append(np.around(curr_gap, self._num_dec_places).tolist())

				# B+T
				self._num_dec_places = self.max_num_dec_places
				curr_gap = self.helper_col_analysis(col_name, stat, append_average = append_average,
								use_unit_vectors_heur = 0, use_split_share = 0, num_cuts_iter_bilinear = [], use_cut_vert_heur = 0, use_tight_points_heur = [])
				self._num_dec_places = saved_num_dec_places
				out_tab_tr.append(np.around(curr_gap, self._num_dec_places).tolist())

				# B+V
				self._num_dec_places = self.max_num_dec_places
				curr_gap = self.helper_col_analysis(col_name, stat, append_average = append_average,
								use_unit_vectors_heur = 0, use_split_share = 0, num_cuts_iter_bilinear = [], use_cut_vert_heur = [], use_tight_points_heur = 0)
				self._num_dec_places = saved_num_dec_places
				out_tab_tr.append(np.around(curr_gap, self._num_dec_places).tolist())

				# B+R+S
				self._num_dec_places = self.max_num_dec_places
				curr_gap = self.helper_col_analysis(col_name, stat, append_average = append_average,
								use_unit_vectors_heur = [], use_split_share = [], num_cuts_iter_bilinear = [], use_cut_vert_heur = 0, use_tight_points_heur = 0)
				self._num_dec_places = saved_num_dec_places
				out_tab_tr.append(np.around(curr_gap, self._num_dec_places).tolist())

				# B+R+T
				self._num_dec_places = self.max_num_dec_places
				curr_gap = self.helper_col_analysis(col_name, stat, append_average = append_average,
								use_unit_vectors_heur = [], use_split_share = 0, num_cuts_iter_bilinear = [], use_cut_vert_heur = 0, use_tight_points_heur = [])
				self._num_dec_places = saved_num_dec_places
				out_tab_tr.append(np.around(curr_gap, self._num_dec_places).tolist())
				
				# B+R+V
				self._num_dec_places = self.max_num_dec_places
				curr_gap = self.helper_col_analysis(col_name, stat, append_average = append_average,
								use_unit_vectors_heur = [], use_split_share = 0, num_cuts_iter_bilinear = [], use_cut_vert_heur = [], use_tight_points_heur = 0)
				self._num_dec_places = saved_num_dec_places
				out_tab_tr.append(np.around(curr_gap, self._num_dec_places).tolist())

		# R+S
		self._num_dec_places = self.max_num_dec_places
		curr_gap = self.helper_col_analysis(col_name, stat, append_average = append_average,
						use_unit_vectors_heur = [], use_split_share = [], num_cuts_iter_bilinear = 0, use_cut_vert_heur = 0, use_tight_points_heur = 0)
		self._num_dec_places = saved_num_dec_places
		out_tab_tr.append(np.around(curr_gap, self._num_dec_places).tolist())

		# R+T
		self._num_dec_places = self.max_num_dec_places
		curr_gap = self.helper_col_analysis(col_name, stat, append_average = append_average,
						use_unit_vectors_heur = [], use_split_share = 0, num_cuts_iter_bilinear = 0, use_cut_vert_heur = 0, use_tight_points_heur = [])
		self._num_dec_places = saved_num_dec_places
		out_tab_tr.append(np.around(curr_gap, self._num_dec_places).tolist())

		# R+V
		self._num_dec_places = self.max_num_dec_places
		curr_gap = self.helper_col_analysis(col_name, stat, append_average = append_average,
						use_unit_vectors_heur = [], use_split_share = 0, num_cuts_iter_bilinear = 0, use_cut_vert_heur = [], use_tight_points_heur = 0)
		self._num_dec_places = saved_num_dec_places
		out_tab_tr.append(np.around(curr_gap, self._num_dec_places).tolist())

		# S+T
		self._num_dec_places = self.max_num_dec_places
		curr_gap = self.helper_col_analysis(col_name, stat, append_average = append_average,
						use_unit_vectors_heur = 0, use_split_share = [], num_cuts_iter_bilinear = 0, use_cut_vert_heur = 0, use_tight_points_heur = [])
		self._num_dec_places = saved_num_dec_places
		out_tab_tr.append(np.around(curr_gap, self._num_dec_places).tolist())

		# S+V
		self._num_dec_places = self.max_num_dec_places
		curr_gap = self.helper_col_analysis(col_name, stat, append_average = append_average,
						use_unit_vectors_heur = 0, use_split_share = [], num_cuts_iter_bilinear = 0, use_cut_vert_heur = [], use_tight_points_heur = 0)
		self._num_dec_places = saved_num_dec_places
		out_tab_tr.append(np.around(curr_gap, self._num_dec_places).tolist())

		# T+V
		self._num_dec_places = self.max_num_dec_places
		curr_gap = self.helper_col_analysis(col_name, stat, append_average = append_average,
						use_unit_vectors_heur = 0, use_split_share = 0, num_cuts_iter_bilinear = 0, use_cut_vert_heur = [], use_tight_points_heur = [])
		self._num_dec_places = saved_num_dec_places
		out_tab_tr.append(np.around(curr_gap, self._num_dec_places).tolist())

		# R+S+T
		self._num_dec_places = self.max_num_dec_places
		curr_gap = self.helper_col_analysis(col_name, stat, append_average = append_average,
						use_unit_vectors_heur = [], use_split_share = [], num_cuts_iter_bilinear = 0, use_cut_vert_heur = 0, use_tight_points_heur = [])
		self._num_dec_places = saved_num_dec_places
		out_tab_tr.append(np.around(curr_gap, self._num_dec_places).tolist())

		# R+S+V
		self._num_dec_places = self.max_num_dec_places
		curr_gap = self.helper_col_analysis(col_name, stat, append_average = append_average,
						use_unit_vectors_heur = [], use_split_share = [], num_cuts_iter_bilinear = 0, use_cut_vert_heur = [], use_tight_points_heur = 0)
		self._num_dec_places = saved_num_dec_places
		out_tab_tr.append(np.around(curr_gap, self._num_dec_places).tolist())

		# R+T+V
		self._num_dec_places = self.max_num_dec_places
		curr_gap = self.helper_col_analysis(col_name, stat, append_average = append_average,
						use_unit_vectors_heur = [], use_split_share = 0, num_cuts_iter_bilinear = 0, use_cut_vert_heur = [], use_tight_points_heur = [])
		self._num_dec_places = saved_num_dec_places
		out_tab_tr.append(np.around(curr_gap, self._num_dec_places).tolist())
		
		# S+T+V
		self._num_dec_places = self.max_num_dec_places
		curr_gap = self.helper_col_analysis(col_name, stat, append_average = append_average,
						use_unit_vectors_heur = 0, use_split_share = [], num_cuts_iter_bilinear = 0, use_cut_vert_heur = [], use_tight_points_heur = [])
		self._num_dec_places = saved_num_dec_places
		out_tab_tr.append(np.around(curr_gap, self._num_dec_places).tolist())

		# Set header and transpose
		out_tab = [["Instance", "Best", 
								"R", "S", "T", "V", 
								#"B", 
								#"B+R", "B+S", "B+T", "B+V",
								#"B+R+S", "B+R+T", "B+R+V", "B+V+S", "B+V+T", "B+S+T",
								"R+S", "R+T", "R+V", 
								"S+T", "S+V", "T+V", 
								"R+S+T", "R+S+V", "R+T+V",
								"S+T+V"]]
		out_tab.extend([list(x) for x in izip(*out_tab_tr)])

		#util.print_pretty_table(out_tab_tr)
		util.print_pretty_table(out_tab)
		return out_tab

	def point_analysis(self):
		"""
		Analyzes number of points and final points

		Outputs 2D list with the data
		"""
		col_name = [self.colname_numpoints, self.colname_numfinalpoints]
		stat = StatEnum.MAX
		make_int = True
		out_tab_tr = [super(LGReader, self).get_param(self.colname_inst)]
		total_num_act = 4

		for num_act in range(0,total_num_act):
			tab = super(LGReader, self).stat_col(col_info = col_name, num_alg2_rounds = num_act, stat = stat)
			for c in range(len(col_name)):
				if (stat in [StatEnum.MIN, StatEnum.MAX, StatEnum.FIRST]):
					curr_index = 1 + 2 * c	# tab cols are [inst, stat, row, stat, row, ..., num_inst]
				elif (stat in [StatEnum.AVG]):
					curr_index = 1 + c	# tab cols are [inst, stat, stat, ..., num_inst]
				if (make_int):
					curr_col = [int(float(tab[i][curr_index])) for i in range(1,len(tab))]	# one header row
				else:
					curr_col = [float(tab[i][curr_index]) for i in range(1,len(tab))]	# one header row
				out_tab_tr.append(curr_col)
		
		# Rearrange (same columns should be next to each other)
		out_tab_tr = [
			out_tab_tr[0],
			out_tab_tr[1],
			out_tab_tr[3],
			out_tab_tr[5],
			out_tab_tr[7],
			out_tab_tr[2],
			out_tab_tr[4],
			out_tab_tr[6],
			out_tab_tr[8]
		]
		
		nInst = len(out_tab_tr[0])
		percent_final_points = [
						[
								round(100. * float(out_tab_tr[j][i]) / out_tab_tr[j-total_num_act][i], self._num_dec_places)
								if out_tab_tr[j-total_num_act][i] != 0 
								else 0 
								for j in range(5,9)
						]
						for i in range(nInst)
				]
		max_index = [percent_final_points[i].index(max(percent_final_points[i])) for i in range(nInst)]
		#percent_final_points =\
		#	[
		#		round(max(
		#						[100. * float(out_tab_tr[j][i]) / out_tab_tr[j-4][i] 
		#						if out_tab_tr[j-4][i] != 0 
		#						else 0 
		#						for j in range(5,9)]
		#				), self._num_dec_places) 
		#		for i in range(nInst)
		#	]
		out_tab = [["Instance", "Points", "Final", "% final"]]
		out_tab.extend([
						[out_tab_tr[0][i],
								out_tab_tr[1+max_index[i]][i],
								out_tab_tr[1+max_index[i]+total_num_act][i],
								percent_final_points[i][max_index[i]]]
						for i in range(nInst)
				])

		#out_tab = [
		#	#			['']+4*["Points"]+4*["Final points"]+[''],
		#				["Instance", "+0H", "+1H", "+2H", "+3H", 
		#				"+0H", "+1H", "+2H", "+3H", "% final"]
		#		]
		#out_tab.extend([list(x) for x in izip(*tmp_out_tab_tr)])

		util.print_pretty_table(out_tab)
		return out_tab

	def helper_col_analysis(self, col_name, stat, tab = None,
			tab_col_index = None, make_int = None, return_gap_closed = None, append_average = None,
			use_unit_vectors_heur = None,
			use_split_share = None, num_cuts_iter_bilinear = None, 
			use_cut_vert_heur = None, use_tight_points_heur = None, cut_limit = None):
		""" 
		Helper function (saves repeating of same lines)
		Calculates gap closed from best result for each instance using the given options
		If tab is given, then do not recalculate it

		Returns the numpy array with the gap closed for each instance or the column as a row
		"""
		# Defaults
		if (tab_col_index is None):
				tab_col_index = 1
		if (make_int is None):
			make_int = False
		if (return_gap_closed is None):
			return_gap_closed = True
		if (append_average is None):
			append_average = True
		if (tab is None):
			tab = super(LGReader, self).stat_col(col_info = col_name, stat = stat,
																			use_unit_vectors_heur = use_unit_vectors_heur,
																			use_split_share = use_split_share, 
																			num_cuts_iter_bilinear = num_cuts_iter_bilinear, 
																			use_cut_vert_heur = use_cut_vert_heur,
																			use_tight_points_heur = use_tight_points_heur,
                                                                                                                                                        cut_limit = cut_limit)

		nonzero_indices = [i-1 for i in range(1,len(tab)) if tab[i][3] > 0]

		if (make_int):
			curr_col = [
						int(round(float(tab[i][tab_col_index]))) 
						if tab[i][3] > 0
						else 0
						for i in range(1,len(tab))
				]	# one header row
		else:
			curr_col = [
						float(tab[i][tab_col_index]) 
						if tab[i][3] > 0
						else 0
						for i in range(1,len(tab))
				]	# one header row

		if (return_gap_closed):
			#out_col = 100 * np.true_divide(np.asarray(curr_col) - self._lp_opt, self._ip_opt - self._lp_opt) - self._sics_gap
			out_col = [
						100 * (curr_col[i] - self._lp_opt[i]) / (self._ip_opt[i] - self._lp_opt[i]) #- self._sics_gap[i]
						if tab[i+1][3] > 0
						else 0
						for i in range(0,len(tab)-1)
				]
		else:
			out_col = curr_col

		if (append_average):
			if (isinstance(out_col, np.ndarray)):
				#avg = np.around(np.mean(out_col, dtype=np.float64), self._num_dec_places)
				avg = np.around(np.mean([out_col[i] for i in nonzero_indices], dtype=np.float64), self._num_dec_places)
				out_col = np.append(out_col, avg)
			elif (type(out_col) is list):
				if (len(nonzero_indices) > 0):
					avg = float(sum([out_col[i] for i in nonzero_indices])) / len(nonzero_indices)
				else:
					avg = 0
				out_col.append(avg)
			else:
				raise TypeError("Somehow out_col is not list or numpy.ndarray. Has type %s." % type(out_col))

		if (make_int is not True):
			out_col = np.around(out_col, self._num_dec_places)
		return out_col
	
	def param_analysis(self, cut_limit = None):
		"""
		CURRENTLY NON FUNCTIONING
		Checks effect of # rays cut, limit cuts per split, and cut presolve

		Outputs 2D list with the data
		"""
		if (not hasattr(self, '_gap_closed')):
			self.gap_closed()

		col_name = [self.colname_allcutsbound]
		stat = StatEnum.MAX
		make_int = False
		out_tab_tr = [super(LGReader, self).get_param(self.colname_inst)]
		out_tab_tr.append(np.around(self._sics_gap, self._num_dec_places).tolist())
		out_tab_tr.append(np.around(self._allcuts_gap, self._num_dec_places).tolist())
		out_tab_tr.append(np.around(self._gap_closed, self._num_dec_places).tolist())

		append_average = True
		if (append_average):
			out_tab_tr[0].append("Average")
			out_tab_tr[1].append(np.around(np.mean(self._sics_gap, dtype=np.float64), self._num_dec_places))
			out_tab_tr[2].append(np.around(np.mean(self._allcuts_gap, dtype=np.float64), self._num_dec_places))
			out_tab_tr[3].append(np.around(np.mean(self._gap_closed, dtype=np.float64), self._num_dec_places))

		saved_num_dec_places = self._num_dec_places
		self._num_dec_places = self.max_num_dec_places
		for val in [0,1]:
				tab = super(LGReader, self).stat_col(col_info = col_name, stat = stat, limit_cuts_per_cgs = val, cut_limit = cut_limit)
				self._num_dec_places = saved_num_dec_places
				for c in range(len(col_name)):
					if (stat in [StatEnum.MIN, StatEnum.MAX, StatEnum.FIRST]):
						tab_col_index = 1 + 2 * c	# tab cols are [inst, stat, row, stat, row, ..., num_inst]
					elif (stat in [StatEnum.AVG]):
						tab_col_index = 1 + c	# tab cols are [inst, stat, stat, ..., num_inst]
						
					curr_gap = self.helper_col_analysis(col_name, stat, 
						tab = tab, tab_col_index = tab_col_index,
						append_average = append_average) 
					out_tab_tr.append(np.around(curr_gap, self._num_dec_places).tolist())

		for val in [0,1]:
				tab = super(LGReader, self).stat_col(col_info = col_name, stat = stat, cut_presolve = val, cut_limit = cut_limit)
				self._num_dec_places = saved_num_dec_places
				for c in range(len(col_name)):
					if (stat in [StatEnum.MIN, StatEnum.MAX, StatEnum.FIRST]):
						tab_col_index = 1 + 2 * c	# tab cols are [inst, stat, row, stat, row, ..., num_inst]
					elif (stat in [StatEnum.AVG]):
						tab_col_index = 1 + c	# tab cols are [inst, stat, stat, ..., num_inst]
						
					curr_gap = self.helper_col_analysis(col_name, stat, 
						tab = tab, tab_col_index = tab_col_index,
						append_average = append_average) 
					out_tab_tr.append(np.around(curr_gap, self._num_dec_places).tolist())

		out_tab = [["Instance", "SIC", "GIC", "Best", "Limit 0", "Limit 1", "Presolve 0", "Presolve 1"]]
		out_tab.extend([list(x) for x in izip(*out_tab_tr)])

		util.print_pretty_table(out_tab)
		return out_tab
	
	#def plotScatter(self, data_x, data_y, xTitle, yTitle, plotTitle, saveAsFile, folder):
	#	"""
	#	Makes scatter plot
	#	"""
	#	plt.plot(data_x, data_y, 'o')
	#	plt.xlabel(xTitle)
	#	plt.ylabel(yTitle)
	#	#plt.axis([0,100,0,100])
	#	plt.grid(False)
	#	plt.title(plotTitle + ' scatter')
	#	plt.savefig('/home/akazachk/repos/pointcuts4/data/figs/' + folder + '/' + saveAsFile + '_scatter.png')
	#	#plt.savefig('./figs/'+saveAsFile + '_scatter.png')
	#	plt.clf()

	def run_param_regression(self,
							 inst_list = None, pha_act_option = None,
							 num_alg2_rounds = None, num_rays_cut = None, 
							 cut_limit = None, 
							 limit_cuts_per_cgs = None, use_split_share = None, 
							 num_cuts_iter_bilinear = None, use_cut_vert_heur = None, 
							 use_tight_points_heur = None,
							 cut_presolve = None):
		"""
		Reports correlation and p-val of gap closed vs. various statistics
		(1) # FP, (2) # points, (3) % FP, (4) depth, (5) # cuts
		"""
		if (not hasattr(self, '_gap_closed')):
			self.gap_closed(inst_name, pha_act_option, num_alg2_rounds,
				num_rays_cut, cut_limit, limit_cuts_per_cgs,
				use_split_share, num_cuts_iter_bilinear, use_cut_vert_heur,
				use_tight_points_heur, cut_presolve)
		if (inst_list is None):
				inst_list = super(LGReader, self).get_param("INSTANCE")

		for inst_index in range(len(inst_list)):
				inst_name = inst_list[inst_index]

				if __debug__:
					print( "## Regressions for %s ##" % inst_name )

				# Points regression
				tab = super(LGReader, self).create_table([inst_name],['NUM FINAL POINTS (TOTAL)', 'NUM POINTS (TOTAL)', 'ALL CUTS (WITH SIC)', 'NUM PHA', 'SIC DEPTH POINTS (AVG)'])
				x_final_points = [(int)(tab[i][0]) for i in range(len(tab))]
				x_total_points = [(int)(tab[i][1]) for i in range(len(tab))]
				x_percent_fp = [100 * (float) (x_final_points[i])/x_total_points[i] for i in range (len(x_final_points))]
				x_num_pha = [(int)(tab[i][3]) for i in range(len(tab))]
				x_depth_points = [(float)(tab[i][4]) for i in range(len(tab))]

				curr_lp_opt = self._lp_opt[inst_index]
				curr_sic_opt = self._sic_opt[inst_index]
				curr_ip_opt = self._ip_opt[inst_index]
				y_gap = [100. * ((float)(tab[i][2]) - curr_lp_opt)/(curr_ip_opt - curr_lp_opt) for i in range(len(tab))]
				y_diff_gap = [100. * ((float)(tab[i][2]) - curr_sic_opt)/(curr_ip_opt - curr_lp_opt) for i in range(len(tab))]

				x = x_final_points
				y = y_gap
				self.plotScatter(x,y,'Num Final Points','% Gap Closed','Num FP vs Gap',inst_name,'NumFP')
				y = y_diff_gap
				self.plotScatter(x,y,'Num Final Points','Diff % Gap Closed','Num FP vs Diff Gap',inst_name + '_Diff','Diff/NumFP')
				
				x = x_total_points
				y = y_gap
				self.plotScatter(x,y,'Num Total Points','% Gap Closed','Num Points vs Gap',inst_name,'NumPoints')
				y = y_diff_gap
				self.plotScatter(x,y,'Num Total Points','Diff % Gap Closed','Num Points vs Diff Gap',inst_name + '_Diff','Diff/NumPoints')

				x = x_percent_fp
				y = y_gap
				self.plotScatter(x,y,'% Final Points','% Gap Closed','% FP vs Gap',inst_name,'PercentFP')
				y = y_diff_gap
				self.plotScatter(x,y,'% Final Points','Diff % Gap Closed','% FP vs Diff Gap',inst_name + '_Diff','Diff/PercentFP')

				x = x_num_pha
				y = y_gap
				self.plotScatter(x,y,'Num PHA','% Gap Closed','Num PHA vs Gap',inst_name,'NumCuts')
				y = y_diff_gap
				self.plotScatter(x,y,'Num PHA','Diff % Gap Closed','Num PHA vs Diff Gap',inst_name + '_Diff','Diff/NumCuts')

				x = x_depth_points
				y = y_gap
				self.plotScatter(x,y,'SIC Depth Points (Avg)','% Gap Closed','SIC Depth vs Gap',inst_name,'SICDepth')
				y = y_diff_gap
				self.plotScatter(x,y,'SIC Depth Points (Avg)','Diff % Gap Closed','SIC Depth vs Diff Gap',inst_name + '_Diff','Diff/SICDepth')

				#(slope,intercept,rval,pval,stderr)=stats.linregress(x,y)
				#print('inst=%s regression: slope=%f intercept=%f, rval=%f, pval=%f, std error= %f' % (inst_name,slope,intercept,rval,pval,stderr))


	def fill_ip_opt(self, ip_opt_fname = None, out_fname = None, overwrite = None):
		"""
		Fills in IP opt for each instance in the relevant column
		Creates a processcsv.ProcessCSV instance for lg_info and ip_opt
		Finds IP opt for each row in lg_info, and creates a new file (out_f) with all info
		"""
		if __debug__:
			print( "## Filling in IP values ##" )
		if (overwrite is None):
			overwrite = False
		
		# Process parameters
		assert self._in_fname is not None
		find_dot = util.index_of('.', self._in_fname)
		if (find_dot >= 0):
			in_fname_stub = self._in_fname[:find_dot]
		else:
			in_fname_stub = copy.deepcopy(self._in_fname)

		if (out_fname is None):
			out_fname = in_fname_stub + "_ip.csv"
		elif (out_fname == self._in_fname):
			overwrite = True
			out_fname = in_fname_stub + "_ip.csv"
		else:	# In case user mistakenly set overwrite = True
			overwrite = False
		
		if __debug__:
			print( "Infile: %s, Outfile: %s" % (self._in_fname, out_fname) )

		# Read IP opt file in
		self.get_ip_opt(ip_opt_fname)
		
		# Open out file
		out_f = open(out_fname, 'wb')
		output = csv.writer(out_f)

		# Write file line by line
		lg_ip_obj_col = super(LGReader, self).get_col_index(self.colname_ipobj)
		assert lg_ip_obj_col >= 0

		# Write header
		#for i in range(len(super(LGReader, self)._header)):
		#	output.writerow(super(LGReader, self)._header[i])
		for i in range(len(self._header)):
			output.writerow(self._header[i])

		# Write each row with filled-in value
		super(LGReader, self).restart_reader()
		#for row in super(LGReader, self)._reader:
		for row in self._reader:
			curr_inst = row[super(LGReader, self).inst_col]
			curr_inst_index = super(LGReader, self).get_inst(curr_inst)

			# find_first_val returns a table, with a header row
			# The first row contains all the column information
			#val_str = ip_opt_reader.find_first_val(col_info = self.colname_ipobj, inst_name = curr_inst)[1][1]
			print( "Curr inst: %s, curr_inst_index: %d" % (curr_inst, curr_inst_index) )
			val_str = str(self._ip_opt[curr_inst_index])	# Might be a float already, but might be the '' string
			if (len(val_str) > 0):
				curr_inst_ip_obj = float(val_str)
				if __debug__:
					print( "Instance: %s\tIP obj: %f" % (curr_inst, curr_inst_ip_obj) )
				row[lg_ip_obj_col] = curr_inst_ip_obj
			elif __debug__:
				print( "Instance %s not found in IP file" % curr_inst )

			output.writerow(row)

		# Close
		out_f.close()
		del output

		# Overwrite if out_fname and in_fname coincide
		if (overwrite):
			# Necessarily will need to recreate it later
			super(LGReader, self).close_f()
			# Overwite in_fname
			shutil.move(out_fname, self._in_fname)
			# Restart the file
			super(LGReader, self).reopen_reader()
	
	def get_command(self, row_info, inst_dir = None, out_dir = None):
		"""
		Get PointCuts command from row
		"""
		if (inst_dir is None):
			inst_dir = "../"
		elif (inst_dir[-1] is not '/'):
			inst_dir = inst_dir + '/'

		if (out_dir is None):
			out_dir = "."

		if (type(row_info) is int):
			row = super(LGReader, self).get_row(row_info)[0:self._param_container.num_params-2]
		elif (type(row_info) is list):
			assert(len(row_info) >= self._param_container.num_params-2)
			row = row_info[0:self._param_container.num_params-2]
		else:
			raise TypeError("Type of row is not int (for row number) or list (the row itself), but is %s." % type(row_info))

		return "$POINTCUTS" + " " + inst_dir + row[0] + ".mps" + " " + out_dir + " " + " ".join(row[1:])

	def __del__(self):
		""" Destructor """
		super(LGReader, self).__del__()


def main(argv):
	from os import environ
	HOME_DIR = os.environ['HOME']
	test_dir = HOME_DIR + "/repos/pha2/results/saved_results"
	test_name_stub = "pha.csv"
	test_name = test_dir + '/' + test_name_stub
	reader = LGReader(test_name);
	#reader = LGReader("/Users/akazachk/repos/pha2/results/pha.csv"); 
	reader.set_param("CUT_LIMIT",1000); reader.set_param("CUT_PRESOLVE",0); 
	reader.set_inst("mas74");
	#tab = reader.gap_closed_table()

	#reader.hplane_point_analysis()
	#tab = reader.inst_info()
	#tab = reader.point_analysis()
	#tab = reader.hplane_analysis()
	#tab = reader.cut_heur_analysis()
	#reader.write_best_params()
	#tab = reader.obj_fails_table()

if __name__ == "__main__":
	main(sys.argv)
