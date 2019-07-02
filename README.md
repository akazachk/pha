# Partial hyperplane activation for generalized intersection cuts #
##### Name: README
##### Author: Aleksandr M. Kazachkov, Carnegie Mellon University and Polytechnique Montreal, aleksandr DOT kazachkov AT polymtl DOT ca, http://andrew.cmu.edu/~akazachk/
##### Code for joint paper with Selvaprabu Nadarajah, Egon Balas, and Francois Margot
##### Last edit: July 2019

Implementation of partial hyperplane activation for generalized intersection cuts.
When example code below has square brackets around an argument, this indicates the argument is optional.
Note that the code may yield different results based on the setup in the system being used; e.g., the starting basis chosen by Clp could be affected.

Included in this package is also a set of Python tools for automatically generating LaTeX tables. These are mostly undocumented (though somewhat commented in the code itself), but their sole purpose is to gather statistics from the log files, the core code for which is contained in the function `stat_col` in `scripts/python_scripts/processcsv/processcsv.py`, which can be used to find the average, minimum, maximum for a set of columns from a csv file for a set of instances run with any given set of parameters (useful for parameter analysis).


## Contents ##
1. Quick start
2. Installation requirements
3. Compilation instructions
4. Parameters
5. Algorithms
6. How do I perform more experiments?
7. How do I get tables?
8. Environment variables
9. Notes
10. Package contents
11. Whom do I talk to?


## Quick start ##

0. Ensure that the `COIN_OR` code is properly set up (use the `install_coin.sh` script under `scripts`).

1. These scripts assume that you first compiled the `Release` version of code. Otherwise, it can be done by calling, from `${PHA_DIR}` (assuming all environment variables and dependencies are set up),

        make release

2. The following will run the _best_ set of parameters for each instance (where best is in terms of percent integrality gap closed). The `1` indicates that, to speed things up, the code is run in batches as in the file `data/instances/all.batch`. *NOTE: if you are on a shared computer, running the command below will start 16 processes in parallel.*

        ${PHA_DIR}/scripts/run_scripts/run_experiments.sh best ${PHA_DIR}/data/instances/all.batch 1

	If a non-batched version is preferred, this can be run with

        ${PHA_DIR}/scripts/run_scripts/run_experiments.sh best ${PHA_DIR}/data/instances/all.instances 0
	
	For reducing the number of batches or instances tested via a different file (instead of `all.batch`, one can use `quick.batch`, for example), see instructions later in this README, under the section "How do I perform more experiments?".

3. The result (assuming the batched default command above is used) is that `results/batches/all`, in each batch subfolder, will contain a file `pha-best.csv`. To merge these into one file, which will be copied into `results/batches/all` and `results`, call

        ${PHA_DIR}/scripts/run_scripts/merge_results.sh ${PHA_DIR}/results/batches/all/ 1 best

4. With this file, Table 5 in the paper, containing the results on the best percent integrality gap closed per instance, can be recreated by running

        mkdir -p ${PHA_DIR}/results/tables
        python ${PHA_DIR}/scripts/python_scripts/get_table.py 1 best

5. The output of this goes to `results/tables` and will contain `gap-closed.csv` and `gap-closed.tex`. The first is a file that can be opened as a spreadsheet to compare against Table 5. The second can be put into LaTeX to view as a table; this requires several packages to view exactly as in the paper: `setspace`, `adjustbox`, `booktabs`, `multirow`, `array`, and a definition of a new column type for hiding columns we print to the file but do not wish to display (as they would make the tables too wide):

        \newcolumntype{H}{>{\setbox0=\hbox\bgroup}c<{\egroup}@{}}

6. Tables 3 and 4 can be reproduced using the data in `results/saved_results/pha.csv`. Simply run 

        cp ${PHA_DIR}/results/saved_results/pha.csv ${PHA_DIR}/results/pha-full.csv
        python ${PHA_DIR}/scripts/python_scripts/get_table.py 1
	
  The saved results contained in `results/saved_results/pha.csv` can be reproduced as well, per the more detailed instructions under "How do I perform more experiments". It took us about a month to do all the runs on all the instances. It is possible to choose only a subset of instances for which the reviewer may wish to reproduce the full set of results, e.g., using the file `data/instances/quick.batch`. We also give instructions below on how to write a file specifying a custom set of instances to test on.

7. The output again goes to `results/tables`. Table 3 corresponds to the file named `hplane-analysis-short.csv`, while Table 4 is generated from the file `cut-heur-analysis-short.csv`. Note that the column "Wins" is computed manually for the paper. Alternatively, the user can refer to the file `results/tables/tables.xlsx` if he/she has access to a program for Excel-formatted files.

8. Table 6 requires the generation of a second round of cuts. These can be produced by modifying `scripts/run_scripts/pha/run_pha_best.py` by setting `rounds=2` instead of `rounds=1` on line 77, then repeating step 2 and step 3. We would advise copying the contents of the resulting `pha_best.csv` into `tables.xlsx` (sheet `best_rd2`), because this spreadsheet has certain functions automatically collecting the relevant columns of the table (under the sheet `gap_rd2`).

9. Table 7 is produced manually by comparing the gap closed in our experiments to the relevant rows in Tables 1 and 2 of the paper "On optimizing over lift-and-project closures" by Pierre Bonami published in _Mathematical Programming Computation_ in 2012.


## What do I need? (Installed on review machine by authors) ##

1. A version of Cgl from COIN-OR (latest tested is 0.59), which comes with Clp (latest tested is 1.16). The code may run differently on different versions. Instructions are below (from main directory). The below will work for the `Release` version of the code; the `Debug` version assumes that instead of `build`, a `buildg` directory is created under `Cgl-0.59`, and typically `enable_debug=yes` is set in the `config.site` file. For the review machine, only the `Release` version is prepared by the authors.

        mkdir -p coin-or
        svn co https://projects.coin-or.org/svn/Cgl/stable/0.59 coin-or/Cgl-0.59
        mkdir coin-or/Cgl-0.59/build
        mkdir coin-or/Cgl-0.59/build/share
        cp scripts/coin_scripts/linux_release_config.site coin-or/Cgl-0.59/build/share/config.site
        cd coin-or/Cgl-0.59/build 
        ../configure -C >& last_config.txt
        tail -n 1 last_config.txt
        make
        make install
        make test

2. BLAS and LAPACK (included in Accelerate framework on Mac).
3. For scripts and compiling tables, Python 2 (maybe will work with Python 3) with `matplotlib` (which in turn requires `six`).


## Installation requirements ##

1. Set up all environment variables (see below).
2. Ensure BLAS and LAPACK are working properly. You will need to use the shared library versions. This has been installed locally on the review machine in `lib/compilation_files`. Note that in this case one will need to adjust the environment variable `LD_LIBRARY_PATH` to point to `${PHA_DIR}/lib`.
3. Ensure Cgl code is properly set up and works on the machine. Under `scripts/coin_scripts`, there are config files that have worked on Macs (`mac_*_config.site`) and Linux (`linux_*_config.site`), where * = debug or release. Change the `SYSTEM` variable in the script to `mac` or `linux` as appropriate.
4. Call

        make dir_release
        make release

    If you want to use the `Debug` version, use

        make dir_debug
        make debug

    For a description of the arguments that the executable takes, run
        
        ${PHA_DIR}/[Debug or Release]/PHA --help


## Parameters ##

The options available can be obtained by running

        ${PHA_DIR}/Release/PHA --help

It is required that you specify an instance (with option `-i`) and the output directory (with option `-o`). You may select where the logging is outputted via `--log_file` and you can allow percent gap closed to be calculated by the code by pointing it to a comma-separated file with the optimal IP values for the instances, which we currently keep in `data/ip_opt.csv`, specified to the code via `--opt_file`.

The code includes various other parameters, discussed in Section 6 of the paper (see Table 2 on page 20). The default parameters can be found in `data/params/pha-params.txt` and can be specified via the `-p` option.

The default parameters are (the last eight remain constant throughout our experiments)

        hplane_scoring_fn 1         // select between the hyperplane scoring functions in Section 5.1 on page 17
        num_alg2_rounds 0           // absolute value gives number of hyperplanes to be activated using Algorithm 2 (parameter k_h) on page 12; when this is negative, targeted tilting is not used
        use_unit_vectors_heur 0     // objective functions (R) from Section 5.2 on page 18
        use_split_share 0           // objective functions (S) from Section 5.2 on page 18
        use_tight_points_heur 1000  // objective functions (T) from Section 5.2 on page 18
        use_cut_vert_heur 0         // objective functions (V) from Section 5.2 on page 18
        sics 1                      // number of rounds of standard intersection cuts
        pha 1                       // boolean, whether to enable GICs via partial hyperplane activation
        num_rays_cut 0              // option (unused) to limit number of rays of initial cone to cut
        max_frac_var 0              // limit on number of splits considered
        num_cuts_iter_bilinear 0    // iterative bilinear scheme mentioned in Section 6.4 on the bottom of page 25
        cut_limit 1000              // stop when this many cuts are generated
        cut_presolve 0              // whether to use presolve during cut generation
        rounds 1                    // number of cut rounds


## Algorithms ##

The algorithms in the paper are implemented in the code in `CglPHA.cpp`. They are named as in the paper. Algorithm 1 appears on line 415, Algorithm 2 on line 698, and Algorithm 3 on line 751. 


## How do I perform more experiments? ##

1. There are three things to decide:

    * Which instances do you want to test on? This is given via a file `{name}.instances` (which is intended for running the instances sequentially) or `{name}.batch` (which is for running the instances in batches), that should be put into `data/instances`. Possible choices are:
        * `all`, every instance tested (available as both `all.instances` and `all.batch`)
        * `quick`, a subset of the instances that runs quickly and is good for testing purposes
        * File of your own making. The format of these instance files can be found in `data/instances`. They all start with a first line, which is the name of the folder (e.g., `all`) in which the results will go with `results` (e.g., `results/all` for the non-batch mode, and `results/batches/all` for the batch mode). In the batch mode file, when a line ends with a `/`, it indicates a new "batch" that will go under a subfolder with the correponding name. For both files, the remaining lines indicate where the instances can be found (in `data/instances'). 
    * What setting of PHA do you want to run? The options are just the best set of parameters per instance, a small test set of parameters per instance, or the full set of parameters used in the paper:
        * `best`
        * `test`
        * `full`
		* Run instances sequentially or in batch mode? Batch mode can be used to speed up computations by running subsets of instances in parallel; for those working on shared machines, be aware that, e.g., with `all.batch`, this leads to 16 processes being run in parallel under option `best` and 48 when using option `full`.

2. Having made these choices, call the following sequence of commands

        ${PHA_DIR}/scripts/run_scripts/run_experiments.sh [best, test, or full] [instance directory] [0 or 1]

    After this completes, you create a merged results file (`pha-test.csv` or `pha-best.csv` or `pha-full.csv`) in `results` using

        ${PHA_DIR}/scripts/run_scripts/merge_results.sh ${PHA_DIR}/results/batches/{folder name} [0 or 1] best

3. We first give an example in which one chooses `batches` and `full`. Note that the full set of parameters takes a very long time to run. When the batch argument is specified as `1` and the instance directory is not given, it is by default `data/instances/batches`. This creates three processes per batch, running parameter settings in parallel for the three different hyperplane activation rules.

        ${PHA_DIR}/scripts/run_scripts/run_experiments.sh full 1

    Once this finishes running, merge the results (which, for this case, will be `pha-full0.csv`, `pha-full1.csv`, and `pha-full2.csv` in each of the batch subfolders) by calling

        ${PHA_DIR}/scripts/run_scripts/merge_results.sh ${PHA_DIR}/results/batches/all/ 1 full

    Then copy the file `pha-full.csv` into `results/`:

        cp ${PHA_DIR}/results/batches/all/pha-full.csv ${PHA_DIR}/results/pha-full.csv

5. As an example of testing the best parameters on all the instances without using batches, call
        
        ${PHA_DIR}/scripts/run_scripts/run_experiments.sh best 0 ${PHA_DIR}/data/instances/all

     And when this finishes, merge by calling

        ${PHA_DIR}/scripts/run_scripts/merge_results.sh ${PHA_DIR}/results/all 0 best


## How do I get tables? ##

1. Copy the file you want to use for compiling tables into `results`. This should be automatically done when calling `merge_results.sh`, but the script is still imperfect. Ensure the folder `resutls/tables` exists.

2. To compile table data, use the python scripts in `scripts/python_scripts`. For example, call (again, substitute `best` with `test` or `full` as appropriate)

        python ${PHA_DIR}/scripts/python_scripts/get_table.py 1 best

3. The results are contained in `results/tables`. The tables that end with `short` correspond to the instances for which either SICs or GICs close some of the integrality gap. Note that this means the _averages_ will not match up to what is in Tables 3 and 4 in the paper, because those tables only have the instances in which GICs improve over SICs, but the remaining lines should match.

4. To get the tables in the paper:

        cp ${PHA_DIR}/results/saved_results/pha.csv ${PHA_DIR}/results/instances/pha-full.csv
        python ${PHA_DIR}/scripts/python_scripts/get_table.py 1 full


## Environment variables ##

These variables should be defined for the environment or locally.

1. Define (the Cgl variables will need to be adjusted based on your setup)

        export PHA_DIR="path to pha directory"
        export PHA_DEBUG="${PHA_DIR}/Debug/PHA"
        export PHA_RELEASE="${PHA_DIR}/Release/PHA"
        export COIN_OR_DIR="${PHA_DIR}/coin-or"
        export CGL_VERSION="Cgl-0.59"
        export CGL_DEBUG="${COIN_OR_DIR}/${CGL_VERSION}/buildg"
        export CGL_RELEASE="${COIN_OR_DIR}/${CGL_VERSION}/build"
        export PYTHONPATH="${PYTHONPATH}:${PHA_DIR}/scripts/python_scripts:${PHA_DIR}/scripts/python_scripts/lib64/python:${PHA_DIR}/scripts/python_scripts/lib/python"

2. BLAS, LAPACK: On Mac, you simply need to have the Accelerate framework (invoked by `-framework Accelerate`). On Linux, define (depending on your environment)

        export ENV_LAPACK_LIB="/usr/lib64"
        export ENV_LAPACK_LIB_NAME="liblapack.so.3"
        export ENV_BLAS_LIB="/usr/lib64"
        export ENV_BLAS_LIB_NAME="libblas.so.3"


## Notes ##

The underlying solver used is OsiClp through COIN-OR. We treat slacks differently than does Clp: we assume all slacks/excess variablesare nonnegative, whereas Osi is enforcing that all slacks have a +1 coefficient (so rows in equality form become ax+s, where s is nonpositive for >= rows, nonnegative for <= rows, and = 0 for = rows). We instead have ax + s = b for ax <= b rows and ax - e = b for ax >= b rows, i.e., we add excess variables. Then in the nonbasic space, the LP opt is the origin, and only cuts with right-hand side 1 cut it off.

For example, consider any row ax >= b. Suppose the slack, s, in that row is basic. Thus, the tableau becomes s = b' - a'x (where b' and a' are obtained by multiplying by the basis inverse). We want to replace e = -s. Since the coefficient on e should ultimately by +1 in the final tableau, we need to multiply (-e) = b' - a'x by -1 to get e = -b' + a'x. The rays are thus just a' in that row, not -a' as they are for the <= rows.

Similarly, we need to complement nonbasic slacks from >= rows. This is treated as we would treat any nonbasic variables at their upper-bound in the basis: negate the rays.

We also remark that the results should all be reproducible insofar as the point-ray collections that are generated. However, the cuts themselves may vary. This can depend on your architecture, the version of Clp being used, and shared resources. As one example, there is a time limit imposed on the PRLP of five seconds. If this is encountered more frequently under your setup, you may see fewer cuts from the same collection. In addition, as we are solving Clp iteratively, there is some dependence on the previous iteration's solution. So one timeout can cause a ripple effect on subsequent optimal solutions (either causing more timeouts, or the choosing of different optimal bases, corresponding to different cuts).

Lastly, in our experience, the instance `mod008` is particularly sensitive to the time limit we set and the speed of the review machine.

## Package contents ##
```
.
+-- coin-or/
+-- data/
|   +-- instances/
|       +-- all/
|       +-- all.instances
|       +-- batches.batch
|       +-- best.batch
|       +-- quick.instances
|   +-- ip_opt.csv
|   +-- params/
|       +-- pha-params.txt
+-- makefile
+-- README.md
+-- results/
|   +-- saved_results/
|       +-- pha.csv
|   +-- tables/
|       +-- tables.xlsx
+-- scripts/
|   +-- coin_scripts/
|   +-- python_scripts/
|   +-- run_scripts/
|       +-- get_errors.sh
|       +-- merge_batches_errors.sh
|       +-- merge_batches.sh
|       +-- merge_results.sh
|       +-- run_experiments.sh
|       +-- pha/
|           +-- run_pha_best.py
|           +-- run_pha_full0.py
|           +-- run_pha_full1.py
|           +-- run_pha_full2.py
|           +-- run_pha_sic.py
|           +-- run_pha_test.py
+-- src/
|   +-- classes/
|       +-- AdvCut.cpp
|       +-- AdvCut.hpp
|       +-- Hplane.cpp
|       +-- Hplane.hpp
|       +-- IntersectionInfo.hpp
|       +-- Point.cpp
|       +-- Point.hpp
|       +-- Ray.cpp
|       +-- SolutionInfo.cpp
|       +-- SolutionInfo.hpp
|       +-- Vertex.cpp
|       +-- Vertex.hpp
|   +-- cuts/
|       +-- BBHelper.cpp
|       +-- BBHelper.hpp
|       +-- CglPHA/
|           +-- CglPHA.cpp
|           +-- CglPHA.hpp
|           +-- hplaneActivation.cpp
|           +-- hplaneActivation.hpp
|       +-- CglSIC/
|           +-- CglSIC.cpp
|           +-- CglSIC.hpp
|       +-- CutHelper.cpp
|       +-- CutHelper.hpp
|   +-- PHA.cpp
|   +-- utility/
|       +-- CglGICParam.cpp
|       +-- CglGICParam.hpp
|       +-- Debug.cpp
|       +-- Debug.hpp
|       +-- GlobalConstants.cpp
|       +-- GlobalConstants.hpp
|       +-- optionhelper.cpp
|       +-- optionhelper.hpp
|       +-- optionparser.h
|       +-- Output.cpp
|       +-- Output.hpp
|       +-- Stats.hpp
|       +-- typedefs.hpp
|       +-- Utility.cpp
|       +-- Utility.hpp
```


## Whom do I talk to? ##

Aleksandr M. Kazachkov

aleksandr DOT kazachkov AT polymtl DOT ca

http://andrew.cmu.edu/~akazachk/

Some of this code was also written with Selvaprabu Nadarajah, and some of it is based on code from Francois Margot and Giacomo Nannicini (from the implementation of `CglGMI`).
