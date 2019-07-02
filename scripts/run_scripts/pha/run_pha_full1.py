####
## The inputted file is either [filename].instances or [file].batch.
## (The extension is not important.)
## The first line of this file will be the directory ``stub'',
## and output will be sent to ${PROJ_DIR}/results/instances/stub if batch mode is off,
## and to ${PROJ_DIR}/results/instances/batches/stub/batchname if it is on.
## Each line contains either a relative input path to an instance (with or without the extension) or a batch name.
## The path is relative to ${PROJ_DIR}/data/instances, e.g., the line will be miplib2/bm23.
## A batch name is distinguished by having the batch end with a '/', e.g., '2/' or 'batch2/'.

## Set up proper path variables
import os
PROJ_DIR = os.path.abspath(os.environ['PHA_DIR'])
EXECUTABLE = PROJ_DIR + "/Release/PHA"

## Solver options
phaActOptions = [1]
numAlg2Rounds = [0,1,2,3,4,-1,-2,-3,-4]
#numRaysToBeCut = [0]
cutLimit = [1000] # negative means divide the limit across splits
useSplitShare = [0,1000]
numCutsIterBilinear = [0]
useUnitVectorsHeur = [0,1]
useCutVertHeur = [0,1]
useTightPointsHeur = [0,1000]
usePresolve = [0,1]

## Set up output and input folders
results_path = PROJ_DIR + '/results'
paramfile = PROJ_DIR + '/data/params/pha_params.txt'
#instances_path = os.getcwd()
instances_path = PROJ_DIR + "/data/instances"
instances_file = instances_path + '/' + "test.instances"

outinfo_stub = 'pha-full' + str(phaActOptions[0])
outinfo_dir = results_path

## Get arguments
from sys import argv
use_batches = False  # set to true/false depending on if mps files are all in one folder or divided up into subfolders
if (len(argv) > 1):
  use_batches = True if argv[1] in ['true', 'True', '1', 't'] else False
  if (use_batches and len(argv) < 2):
    raise ValueError('When using batches, specifying the folder is required')

if (len(argv) > 2):
  instances_file = os.path.abspath(argv[2])

## Where are the instances?
with open(instances_file) as f_in:
  list_to_use = list(filter(None, (line.rstrip() for line in f_in)))

## The first line will be the name of the directory we should use
dir_stub = list_to_use[0]
list_to_use = list_to_use[1:]

#instances_file_name = instances_file.split('/')[-1]
#instances_file_name_split_by_dot = instances_file_name.split('.')
#dir_stub = '.'.join(instances_file_name_split_by_dot[0:len(instances_file_name_split_by_dot)-1])
if use_batches:
  dir_stub = "batches/" + dir_stub

## Finalize outinfo
outinfo_dir = outinfo_dir + '/' + dir_stub 
os.system("mkdir -p " + outinfo_dir)  # make the dir if it does not exist

## Choose order so that deepest for loop are the results you want to see first, fixing all others
batch_name = ''
for usepresolve in usePresolve:
#  for numrays in numRaysToBeCut:
  for cutlimit in cutLimit:
    for numiterblp in numCutsIterBilinear:
      for splitshareopt in useSplitShare:
        for usecutvert in useCutVertHeur:
            for useunitvec in useUnitVectorsHeur:
              for usetight in useTightPointsHeur:
                  for numalg2 in numAlg2Rounds:
                    for actoption in phaActOptions:
                      ## Skip if all zeroes for cut generation
                      if (splitshareopt == 0) and (numiterblp == 0) and (useunitvec == 0) and (usecutvert == 0) and (usetight == 0):
                        continue

                      for inst in list_to_use:
                        ## Check if batch name
                        if (inst[-1] == '/'):
                          batch_name = inst
                          continue

                        ## Check if need to add "mps"
                        inst_name = inst
                        if (inst[-4:] != '.mps') and (inst[-3:] != '.lp') and (inst[-7:] != '.mps.gz') and (inst[-6:] != '.lp.gz'):
                          inst_name = inst_name + '.mps'

                        ## Run on instances_path/inst.mps
                        infile = instances_path + '/' + inst_name
                        curr_out_dir = outinfo_dir + '/' + batch_name
                        outinfo = curr_out_dir + outinfo_stub
              
                        ## In case the out directory does not exist
                        os.system("mkdir -p " + curr_out_dir)

                        ## Arguments
                        extraparams = \
                          ' --opt_file=' + PROJ_DIR + '/data/ip_opt.csv' + \
                          " --hplane_scoring_fn=" + str(actoption) + \
                          " --num_alg2_rounds=" + str(numalg2) + \
                          " --cut_limit=" + str(cutlimit) + \
                          " --use_split_share=" + str(splitshareopt) + \
                          " --num_cuts_iter_bilinear=" + str(numiterblp) + \
                          " --use_unit_vectors_heur=" + str(useunitvec) + \
                          " --use_cut_vert_heur=" + str(usecutvert) + \
                          " --use_tight_points_heur=" + str(usetight) + \
                          " --cut_presolve=" + str(usepresolve) + \
                          " --rounds=" + str(1)
                          
                        print(EXECUTABLE + " -i " + infile + " -o " + curr_out_dir + " --log_file=" + outinfo + " -p " + paramfile + extraparams) 
                          
                        os.system(EXECUTABLE + " -i " + infile + " -o " + curr_out_dir + " --log_file=" + outinfo + " -p " + paramfile + extraparams + " > /dev/null 2>&1")
