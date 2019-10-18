# This is Free Software - You can use and distribute it under
# the terms of the GNU General Public License, version 3 or later
# (c) Massimo Cavallaro (m.cavallaro@warwick.ac.uk)

import matplotlib
matplotlib.use('Agg')
import pymc as pymc
import numpy as np
import datetime as dt
import os
import sys
import utils


def Usage():
    print("Usage: ", sys.argv[0], 'data_file_name')
    sys.exit()


n_steps = 100000
thinning = 5000


def test():
    A = pymc.NegativeBinomial('A', mu=100, alpha=5, size=100)
    synth_data = A.rand()
    index = np.argsort(synth_data)
    data = synth_data[index]  
    main(data, n_steps, './', 0)


def run(data_file_name, n_steps):
    data = np.genfromtxt(data_file_name)
    data = np.random.choice(data, size=500, replace=False)
    folder = data_file_name.strip('.dat')
    status = main(data, n_steps, folder, -1)


def main(data, n_steps, folder, verbose):
    """
    This is used for fitting to synthetic mRNA counts data and smFISH.
    Parameters
    ----------
    data : np.array
        data
    n_steps : int
    folder : str
    verbose : int
    """
    if n_steps < thinning:
        print('n_steps must be larger than %d, as we want %dsamples'%(thinning, thinning))
        return 0

    K = pymc.Gamma('K', alpha=0.001, beta=0.001, value=100)
    Mean = pymc.Gamma('Mean', alpha=0.001, beta=0.001, value=5)

    X = pymc.NegativeBinomial('X', mu=Mean, alpha=K, observed=True, value=data)

    path = folder + '_' + dt.datetime.now().strftime("%d-%H.%M")
    if not os.path.exists(path):
       os.makedirs(path)

    MCMC = pymc.MCMC([K, Mean, X], db='hdf5', dbname=path + '/BetaPoisson.h5')
    MCMC.sample(iter=int(n_steps * 3. / 2.), burn=int(n_steps / 2.), thin=int(n_steps / thinning), verbose=verbose, tune_throughout=True)

    pymc.Matplot.plot(MCMC, path=path)
    MCMC.save_state()

    MCMC.db.close()

    return 1


if __name__ == '__main__':
  if len(sys.argv) != 2:
    Usage()

  data_file_name = sys.argv[1]
  run(data_file_name, n_steps)

