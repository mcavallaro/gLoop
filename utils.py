# This is Free Software - You can use and distribute it under
# the terms of the GNU General Public License, version 3 or later
# (c) Massimo Cavallaro (m.cavallaro@warwick.ac.uk)

import sys
import numpy as np
import glob

size_block = 10
scaling_divisor = 20.


def format_array(data, size_max):
  """
  Chose a number multiple of size_block of elements from the dataset.
  Sort and return the dataset.
  """
  if (size_max % size_block) != 0:
    print('size_max must be multiple of %s'%str(size_block))
    sys.exit()
  l = len(data)
  l_size_block = int(l / size_block) * size_block
  l_size_block = min(size_max, l_size_block)
  # len_data = len_data - len_data % size_block
  tmp = np.random.choice(data, l_size_block, replace=False)
  tmp.sort()
  return tmp / scaling_divisor


def format_files(size_max):
  files = glob.glob('*mRNA.csv*')
  for f in files:
    data = np.genfromtxt(f)
    f = f.strip('Experiment').strip('mRNA.csv')
    data = Process_array(data, size_max)
    np.savetxt(f + str(size_max) + '.csv', data)


def fit_BCK(file_name):
  """
  Fit skew normal to background measures.
  Use R package 'sn' to find the std error of the estimates. 
  """
  from scipy.stats import skewnorm
  data = np.genfromtxt(file_name)
  a, mu, sd = skewnorm.fit(data)
  np.savetxt('BCK.csv', [mu, sd, a])
