
# Gloop

This repository contains supporting software to reference [1]. Please cite [1] if you find this repository useful. The software is organised as follows.


  * `R` scripts for tidying flow cytometry `.fcs` data and resolving cell-cycle stages (G1/M/G2).

    - `clustering_2.R`
    - `clustering_3.R`
    - `clustering_caller.R`
    - `compensate_caller.R`
    - `CompensateFlowSet.R`

    These require `flowCore` and `flowClust` [2,3].
        
  * `c++` implementation of the Gillespie algorithm for the simulation of gene expression based of the reaction network.
        It requires the GNU Scientific Library (GSL) ver. 2> [4].
      By default it saves the simulation results into a directory named `.\results`).
      A *makefile* is provided for initial setting, compilation, and linking with `gcc`.


    ```{bash}
    make install
    make
    ```

    Simulation parameters are passed from `STDIN`, e.g.:
    ```
     ./main.exe t N $\alpha$ $\beta$ d $\lambda_{on}$ $\lambda_{on}$ l
    ```

  * MCMC samplers implemented in `python` and `pymc` [5] for the three phenomelogical models described in [1]:

    - `BetaPoissonModel.py`
    - `NegativeBinomialModel.py`
    - `PoissonModel.py`

    these can be launched as:
    ```
    python BetaPoissonModel.py data_file_name $\mu_X$ $se_{\mu_X}$ BCK_params_file_name
    ```
    
      `BCK_params_file_name` contains the parameters obtained from the controll cell.
 
      Some diagonistic methods are imported from `pymc3` [6].
      The sampler for the model with no measurement equation and the utils to rescale
      and tidy the calibration data are in separated files:
         
    - `NegativeBinomialModel_no_error.py`
    - `utils.py`

  * `R` and `bash` scripts for the bioinformatic interrogation of *ChIA-PET* data to extract 3'-5' interaction scores *genome wide*:

    - `Download_trim_chia_pet2.sh`
    - `makehicmatrix.sh`
    - `calculate_3_5_interaction_res_2000.r`
    - `loopscore_functions.r`





[1] M. Cavallaro, *et al.*, 3'-5' interactions contribute to transcriptional bursting, bioR$\chi$iv 514174. https://doi.org/10.1101/514174 

[2] F. Hahne, *et al.*, flowCore: a Bioconductor package for high throughput flow cytometry., BMC Bioinformatics. 10 (2009) 106. http://www.ncbi.nlm.nih.gov/pubmed/19358741

[3] K. Lo *et al.*, flowClust: a Bioconductor package for automated gating of flow cytometry data, BMC Bioinformatics. 10 (2009) 145. https://www.ncbi.nlm.nih.gov/pubmed/19442304

[4] M. Galassi *et al.*, GNU Scientific Library Reference Manual (3rd Ed.), ISBN 0954612078. http://www.gnu.org/software/gsl/

[5] A. Patil, D. Huard, C. Fonnesbeck, PyMC: Bayesian Stochastic Modelling in Python, J. Stat. Softw. 35 (2010) 1–81. https://pymc-devs.github.io/pymc/

[6] J. Salvatier, T. Wiecki, C. Fonnesbeck, Probabilistic Programming in Python using PyMC3, PeerJ Comput. Sci. 2 (2015) e55. https://docs.pymc.io/
 



