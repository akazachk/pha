import sys
import utility
import argparse

__all__ = ["processcsv", "statenum"]
from processcsv import ProcessCSV
from statenum import StatEnum

#try:
#  from matrix2latex import matrix2latex
#except ImportError:
#  # Really ugly hack to please python3 import mechanisms
#  import sys, os
#  SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.expanduser(__file__)))
#
#  sys.path.insert(0, SCRIPT_DIR)
#  from matrix2latex import matrix2latex
#  matrix2latex = matrix2latex.matrix2latex
#  del sys.path[0]             # NOTE: ensure that matrix2latex does not change sys.path

### MAIN ###
def main(argv):
  parser = argparse.ArgumentParser(description="Write the ip opt for each instance in the input file.")
  parser.add_argument("-i", "--infile", default=default_in_fname, dest="in_fname", nargs='?',
                      help="input file name containing instances in the first column and an IP OBJ column")
  parser.add_argument("-f", "--ipfile", default=default_ip_opt_fname, dest="ip_opt_fname", nargs='?',
                      help="name of file containing the IP optimum objective values")
  parser.add_argument("-o", "--outfile", type=str, dest="out_fname",
                      help="output file name equal to input file but with IP objective filled in")
  args = parser.parse_args(argv)
  
  # Process parameters
  in_fname = args.in_fname
  find_dot = utility.index_of('.', in_fname)
  if (find_dot >= 0):
    in_fname_stub = in_fname[:find_dot]
  else:
    in_fname_stub = in_fname

  overwrite = False
  if (args.out_fname is None):
    out_fname = in_fname_stub + "_ip.csv"
  elif (args.out_fname == in_fname):
    overwrite = True
    out_fname = in_fname_stub + "_ip.csv"
  else:
    out_fname = args.out_fname
  
  if __debug__:
    print( "Infile: %s, Outfile: %s, IP file: %s" % (in_fname, out_fname, args.ip_opt_fname) )

  #fill_in_ip_opt(in_fname, out_fname, args.ip_opt_fname, overwrite)

if __name__ == "__main__":
  main(sys.argv)
