import sys
import os
from unittest import skip 
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='Process command line inputs.')
parser.add_argument('--x_measurement', metavar='x_measurement', type=float, nargs=1)
parser.add_argument('--x_units', metavar='x_units', type=str, nargs=1)

# ==============
#      Main
# ==============     
def main():
  # Get the command line arguments 
  args = parser.parse_args()
  x_measurement = args.x_measurement[0]
  x_units = args.x_units[0]

  # Check that the command line arguments are valid 
  if x_measurement < 0.0 :
    sys.exit('Error: Expected x_measurement to be positive.')

  if x_units == 'mm':
    x_measurement *= 1e-1
  elif x_units == 'm':
    x_measurement *= 1e+2
  elif x_units == 'cm':
    x_measurement *= 1.00
  else:
    sys.exit('Error: Unrecognized x_units. Expected mm, cm, or m.')

  # Read the necessary files from the reflected shock simulation
  output_folder_path = './output'
  result_flow = read_result_file(output_folder_path,'flow',4)
  result_dist = read_result_file(output_folder_path,'dist',2)
  result_time = read_result_file(output_folder_path,'time',2)
  ts_num = ensure_same_length(result_flow,result_dist,result_time,0)
  result_flow = result_flow[0:ts_num,:]
  result_dist = result_dist[0:ts_num,:]
  result_time = result_time[0:ts_num,:]

  # Extract the relevant variables from the outputted data
  u = result_flow[:,1] * 1.0e+5       # flow velocity in cm/s
  t = result_time[:,0] * 1.0e-6       # time in s 
  x = result_dist[:,0]                # distance in cm
  dt = result_time[:,1] * 1.0e-6      # dt in s 
  dx = result_dist[:,1]               # dx in cm
  du = np.diff(u,prepend=u[0])        # du in cm/s

  # Get the derivatives
  du_dx = du / dx                     # du/dx in cm/(s-cm)
  du_dt = du / dt                     # du/dt in cm/(s^2)

  # Write the file containing x, u_rs, and a_rs for each x position after x_measurement
  output_file = open('shock_motion.inp','w')
  for i in range(ts_num):
    x_from_measurement = x[i] - x_measurement
    output_file.write(f'{x_from_measurement:14.5e}  {u[i]:14.5e}  {du_dt[i]:14.5e} \n')

# ===========================
#      Utility functions
# ===========================
def read_result_file(output_folder_path, result_type, max_num_variables):
  """ Function that reads result-<result_type>.dat from MTCR output. """
  result_file_path = output_folder_path + '/result-' + result_type + '.dat'
  line_count = count_lines(result_file_path)
  result = open(result_file_path,'r')
  result_data = np.empty([line_count, max_num_variables])
  for i, line in enumerate(result):
    for j in range(max_num_variables):
      try:
        result_data[i,j] = float(line.split()[j])
      except: ValueError
  result.close()
  return result_data

def ensure_same_length(*args):
  """ Ensure that a list of input numpy arrays have the same length. The final argument should be the dimension to compare the lengths along. """
  dim = args[-1]
  if not isinstance(dim,int):
    sys.exit(f'Error: Dim should be a scalar integer. Instead, dim={dim}')
  length = args[0].shape[dim]
  for arg in args[1:-1]:
    if arg.shape[dim] != length:
      # sys.exit(f'Error: Inputs were not the same length in ensure_same_length. Should be: {length}, instead got: {arg.shape[dim]}')
      length = min(length,arg.shape[dim])
  return length

def ensure_file_exists(file_path):
  """ Check if a file exists and, if it doesn't, throw an error and exit the program. """
  if not os.path.isfile(file_path):
    sys.exit(f'Error: File {file_path} not found. Exiting program.')

def count_lines(file_path):
  """ Count the total number of lines in a file and return the line_count. Throws an error if line_count is zero. """
  ensure_file_exists(file_path)
  file_being_read = open(file_path,'r')
  num_header_lines = count_header_lines(file_being_read)
  line_count = -num_header_lines
  for _ in file_being_read:
    line_count += 1
  if line_count == 0:
    sys.exit(f'Error: File {file_path} was opened successfully but had zero lines.')
  return line_count

def count_header_lines(file_being_read):
  """ Skip the Tecplot header lines for the result file. """
  last_header_line = "ZONETYPE = ORDERED, DATAPACKING = POINT"
  found_last_header_line = False
  num_header_lines=0
  for line in file_being_read:
    num_header_lines += 1
    if line.strip() == last_header_line.strip():
      found_last_header_line = True
      return num_header_lines
  if found_last_header_line == False:
    sys.exit(f'Error: Last header line not found in file being read.')

if __name__ == '__main__':
  main()