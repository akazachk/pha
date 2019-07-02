from __future__ import print_function
import sys

### Useful methods ###
def index_of(val, listOfVals):
  """ 
  Useful method for lists
  (returns -1 instead of ValueError when value not in list)
  """
  try:
    index = listOfVals.index(val)
  except ValueError:
    index = -1

  return index
  
def print_table(tab, sep = None, outfile = None, overwrite = None):
  """ 
  Print table
  tab is the table
  sep is the separator to use between entries
  """
  if (sep is None):
    sep = '\t'
  if (overwrite is None):
    overwrite = False

  # Is an outstream provided?
  close_stream = False
  if (outfile is not None):
    if (type(outfile) is str):
      if (overwrite):
        outfile = open(outfile, "wb")
      else:
        outfile = open(outfile, "ab")
      close_stream = True
  else:
    outfile = sys.stdout

  out_str = ""
  for i in range(len(tab)):
    for j in range(len(tab[i])):
      curr_str = "%s" % (tab[i][j])
      if (j is not len(tab[i]) - 1):
        curr_str = curr_str + sep
      out_str = out_str + curr_str
      print( curr_str, end='', file=outfile )
    out_str = out_str + '\n'
    print( "", file=outfile )

  outfile.flush()
  if (close_stream):
    outfile.close()
  
  return out_str
    
def print_pretty_table(tab):
  """ Print formatted table """
  maxlen = 2 + max([len(str(tab[i][j])) for i in range(len(tab)) for j in range(len(tab[i]))])
  for i in range(len(tab)):
    for j in range(len(tab[i])):
      print( ("{0: <%d}" % maxlen).format(tab[i][j]), end='' )
    print( "" )
    
def format_float(num_dec_places = 2): 
  """ Returns string "{:.xf}" which can be used to format floats """
  return "{:." + str(num_dec_places) + "f}"

def csv2table(f_name):
  import csv
  data = []
  with open(f_name,  'rb') as f:
    reader = csv.reader(f, delimiter=',')
    for line in reader:
      data.append(line)

  return data
