import matplotlib
matplotlib.use('Agg')
import pymc as pymc
import numpy as np
import datetime as dt
import os
import sys
import utils


def Usage():
    print("Usage: ", sys.argv[0], 'data_file_name mean_mu se_mu BCK_params_file_name [db_file_name(.h5)]')
    print('If using Python > 3.3, the launch the script with: >> PYTHONHASHSEED=0 parallel python3 BetaPoissonModel.py ::: ...')
    sys.exit()


size_block = 10
scaling_divisor = 20.
n_steps = 60000000


def test():
  from scipy.stats import skewnorm
  size = int(100)
  poi = np.random.poisson(5, size=size)
  data = poi * 20 + skewnorm.rvs(a=5, scale=30, size=size) - 10
  index = np.argsort(data)
  data = data[index] / scaling_divisor  
  poi = poi[index]

  mean_mu = 5.
  var_mu = 5.

  mu_BKG = -10
  sd_BKG = 30
  alpha_BKG = 5

  main(data, mean_mu, var_mu, mu_BKG, sd_BKG, alpha_BKG, n_steps, './', None)


def run(data_file_name, mean_mu, var_mu, BCK_params_file_name, db_file_name=None):
    """
    Parameters
    ----------
    data_file_name : string
        Name of the file that contains flow-FISH data.
        Must be formatted with utils.format_files().
    mean_mu : float
        Sample mean of smFISH/nanoString result.
    var_mu : float
        Sample variance of smFISH/nanoString result.
    BCK_params_file_name : string
        Name of the file that contains the parameters estimated with utils.fit_BCK() from the background reads.
    db_file_name : string, optional
        Name of the file that contains the previous h5fs database.
    """

    data = np.genfromtxt(data_file_name)
    mu_BKG, sd_BKG, alpha_BKG = np.genfromtxt(BCK_params_file_name)

    folder = data_file_name.strip('.csv')

    return main(data, mean_mu, var_mu, mu_BKG, sd_BKG, alpha_BKG, n_steps, folder, db_file_name)


def main(data, mean_mu, var_mu, mu_BKG, sd_BKG, alpha_BKG, n_steps, folder, db_file_name):
  """
  Parameters
  ----------
  data : np.array
      flow-FISH data
      Must be formatted with utils.format_array()
  mean_mu : float
      Sample mean of smFISH/nanoString result.
  var_mu : float
      Sample variance of smFISH/nanoString result.
  mu_BKG : float
      location parameter of the SkewNormal fit to the background data.
  sd_BKG : float
      scale parameter of the SkewNormal fit to the background data.
  alpha_BKG : float
      skewness parameter of the SkewNormal fit to the background data.
  n_steps : int
      Number of tuning steps.
  folder : str
      Destination folder basename.
  db_file_name : dict, optional
      Name of the file that contains the previous h5fs database.
  """

  thinning = 5000
  if n_steps < thinning:
    print('n_steps must be larger than %d, as we want %dsamples'%(thinning, thinning))
    return 0

  len_data = len(data)
  n_blocks = int(len_data / size_block)
  data = np.reshape(data, [n_blocks, size_block])

  print("n data points", len_data)
  print("size block", size_block)
  print("n of blocks", n_blocks)

  Kon = pymc.Gamma('Kon', alpha=0.001, beta=0.001)
  Koff = pymc.Gamma('Koff', alpha=0.001, beta=0.001)

  # Dummy values just to initialise MB. To avoid numerical error must not be too close to 0 or 1.
  value_init = np.array([0.51246393, 0.28097763, 0.52042068, 0.79960688, 0.67806853, 0.74132544, 0.58796267, 0.47798822, 0.49119662, 0.33277814])
  Mean = pymc.TruncatedNormal('Mean', mu=mean_mu, tau=1./var_mu, a=0.01, b=np.inf, value=mean_mu)
  MB = pymc.Container(
    [pymc.Beta('MB%i'%i, alpha=Kon, beta=Koff, trace=False,
      size=int(size_block), value=value_init) for i in range(n_blocks)]
    )

  mean_kappa = 27.

  trace_ = [True for i in range(n_blocks)]

  init_X = np.empty([n_blocks, size_block], dtype=int)
  for i in range(n_blocks):
    init_X[i] = (data[i] - min(0, data[i, 0])) / mean_kappa * scaling_divisor
  X = pymc.Container(
    [pymc.Poisson('X%i'%i, mu=Mean * (Kon + Koff) / Kon * MB[i], value=init_X[i], trace=trace_[i]) for i in range(n_blocks)]
    )

  Kappa = pymc.TruncatedNormal('Kappa', mu=mean_kappa, tau=0.02, a=1, b=np.inf, value=mean_kappa)
 
  Y = pymc.Container([
      pymc.SkewNormal('YY%i'%i,
        mu=(mu_BKG + Kappa * X[i]) / scaling_divisor,
        tau=(scaling_divisor / sd_BKG) ** 2,
        alpha=alpha_BKG,
        observed=True, value=data[i]) for i in range(n_blocks)]
    )

  Y_Sim = pymc.Container([
      pymc.SkewNormal('Y%i'%i,
        mu=(mu_BKG + Kappa * X[i]) / scaling_divisor,
        tau=(scaling_divisor / sd_BKG) ** 2,
        alpha=alpha_BKG) for i in range(n_blocks)]
    )

  path = folder + '_' + dt.datetime.now().strftime("%d-%H.%M")
  if not os.path.exists(path):
    os.makedirs(path)

  if db_file_name is None:
    MCMC = pymc.MCMC([Kon, Koff, Mean, Kappa, MB, X, Y, Y_Sim], db='hdf5', dbname=path + '/BetaPoissonFACS.h5')
    MCMC.use_step_method(pymc.AdaptiveMetropolis, [Kon, Koff, Mean, Kappa], scales={Kon:0.1, Koff:1, Mean:0.1, Kappa:0.1}, shrink_if_necessary=True) # reduce scales if rejection rate is too high.
    MCMC.sample(iter=int(n_steps * 3. / 2.), burn=int(n_steps / 2.), thin=int(n_steps / thinning), verbose=-1, tune_throughout=True)
  else:
    db = pymc.database.hdf5.load(db_file_name)
    MCMC = pymc.MCMC([Kon, Koff, Mean, Kappa, MB, X, Y, Y_Sim], db=db)

    corr = MCMC.db.getstate()['step_methods']['AdaptiveMetropolis_Koff_Kappa_Mean_Kon']['proposal_sd']
    cov = np.dot(corr, corr.T)
    MCMC.use_step_method(pymc.AdaptiveMetropolis, [Kon, Koff, Mean, Kappa], shrink_if_necessary=False, cov=cov)

    MCMC.assign_step_methods()
    MCMC.restore_sampler_state() #restore sampler and stochastic
    MCMC.restore_sm_state() #restore step methods

    MCMC.sample(iter=int(n_steps * 3. / 2.), burn=int(n_steps / 2.), thin=int(n_steps / thinning), tune_throughout=True, verbose=-1)

  pymc.Matplot.plot(MCMC, path=path)
  MCMC.save_state() #save the rusult of MCMC.get_state()

  if db_file_name is not None:
    db.close()

  MCMC.db.close()

  return 1


if __name__ == '__main__':
  if (len(sys.argv) < 5) or (len(sys.argv) > 6):
    Usage()

  for i in sys.argv:
    print(i)

  data_file_name = sys.argv[1]
  mean_mu = float(sys.argv[2])
  var_mu = np.square(float(sys.argv[3]))
  BCK_params_file_name = sys.argv[4]

  if len(sys.argv) == 6:
    db_file_name = sys.argv[5]
    print('Continuing from previous chain %s'%db_file_name)
    run(data_file_name, mean_mu, var_mu, BCK_params_file_name, db_file_name)
  elif len(sys.argv) == 5:
    print('Sampling from default initial conditions')
    run(data_file_name, mean_mu, var_mu, BCK_params_file_name)
  else:
    Usage() 
