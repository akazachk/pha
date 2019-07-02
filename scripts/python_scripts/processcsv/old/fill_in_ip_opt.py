import processcsv
import csv
import utility
import argparse
import sys
import shutil

# Default values
default_ip_opt_fname = "ip_opt.csv"
default_in_fname = "lg-info.csv"
default_out_fname = "lg-info-ip.csv"
inst_col = 0

def fill_in_ip_opt(in_fname, out_fname, ip_opt_fname, overwrite = None):
  """
  Fills in IP opt for each instance in the relevant column
  Creates a processcsv.ProcessCSV instance for lg_info and ip_opt
  Finds IP opt for each row in lg_info, and creates a new file (out_f) with all info
  """
  if (overwrite is None):
    overwrite = False

  # Read IP opt file in
  ip_opt_reader = processcsv.ProcessCSV(ip_opt_fname, num_header_lines = 1)
  
  # Read lg info
  lg_info_reader = processcsv.ProcessCSV(in_fname, num_header_lines = 2)

  # Open out file
  out_f = open(out_fname, 'w')
  output = csv.writer(out_f)

  # Write file line by line
  lg_ip_obj_col = lg_info_reader.get_col_index("IP OBJ")
  assert lg_ip_obj_col >= 0

  # Write header
  for i in range(len(lg_info_reader._header)):
    output.writerow(lg_info_reader._header[i])

  # Write each row with filled-in value
  for row in lg_info_reader._reader:
    curr_inst = row[inst_col]

    # find_first_val returns a table, with a header row
    # The first row contains all the column information
    val_str = ip_opt_reader.find_first_val(col_info = "IP OBJ", inst_name = curr_inst)[1][1]
    if (len(val_str) > 0):
      curr_inst_ip_obj = float(val_str)
      if __debug__:
        print( "Instance: %s\tIP obj: %f" % (curr_inst, curr_inst_ip_obj) )
      row[lg_ip_obj_col] = curr_inst_ip_obj

    output.writerow(row)

  # Close
  out_f.close()

  if (overwrite):
    # Overwite in_fname
    shutil.move(out_fname, in_fname)
