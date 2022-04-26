# Code for „Addressing COVID-19 Outliers in BVARs with Stochastic Volatility“ by Carriero, Clark, Marcellino and Mertens 

Andrea Carriero (Queen Mary University of London), Todd Clark (Federal Reserve Bank of Cleveland), Massimiliano Marcellino (Bocconi University, IGIER and CEPR) and Elmar Mertens (Deutsche Bundesbank)

Accepted for publication in the *Review of Economics and Statistics*

The usual disclaimers apply; the views and results conveyed in our research are solely those of the authors and do not necessarily reflect the views of the Federal Reserve Bank of Cleveland, the Federal Reserve System, the Eurosystem, or the Deutsche Bundesbank.

Pre-publication versions of paper and supplementary material can be found here: https://www.elmarmertens.com/research/publications

## Overview

This repository provides replication codes and raw input data to replicate (and update) the results shown in our paper and its supplementary appendix. All core scripts are in the main directory. In addition, there are the following subdirectories:
- `data` scripts for data construction based on input files obtained from FRED-MD, available at https://research.stlouisfed.org/econ/mccracken/fred-databases/
- `matlabtoolbox` for general utilities (also available at https://github.com/elmarmertens/em-matlabbox)

The code requires a recent version of Matlab (we used Matlab versions 2019a-2021a) including access to Matlab’s Statistics and Machine Learning Toolbox. The codes employ `parfor` loops that are executed in parallel when a `parpool` has been created in Matlab, which requires availability of the Matlab Parallel Computing Toolbox (otherwise the loops will be executed sequentially).


## General notes 

- All scripts set the MATLAB path to point to toolboxes in `matlabtoolbox`. In addition, most scripts collect output in a temporary directory, which is by default created as subfolder `foo` within the main directory. Edit `localtemp.m` to change the location of this temp directory. Output is collected in a LaTeX file, which is also compiled at the end of each script (provided a LaTeX installation can be found on the path). To control the compilation of output, please edit `finishwrap.m`. To avoid collecting output files, comment out the call to `initwrap` in each script (and make sure to define instead a variable called `wrap` that is set to empty).

- Model names in paper and code:

-- `VARconstvcv` is the constant-variance VAR, labeled CONST in the paper
-- `VARSV` is the VAR with (standard) SV, labeled SV in the paper
-- `VARSVO` is the SVO specification
-- `VARSVOt` combines t and SVO, and is labeled SVO-t in the paper
-- `VARSVnanOutlier` is the VAR that treats (pre-specified) outliers as missing data, labeled SV-OutMiss in the paper

- Additional model variants considered in the supplementary material:

-- `VARSVt` is SV with $t$ errors, called SV-t
-- `VARSVdummy` is the VAR with SV and separate dummies for each month of COVID
-- `VARSVobar` is the common-outlier variant of the SVO model.
-- `VARSV*ar1*` are model variants with AR(1) processes (instead of RW) for SV

## To prepare input data

The folder `data` contains code to transform raw data, obtained 
from FRED-MD, into data files for subsequent use by the analytical routines in the repository's root folder. 

The raw data file that serves as input is `2021-04.csv`, which reflects the 2021-04 vintage of FRED-MD and is available under that name from https://research.stlouisfed.org/econ/mccracken/fred-databases/.  

The script `generate_freddata.m` selects the appropriate contents of `2021-04.csv` and prepares an output files called `fredMD16-2021-04.csv` that contains the 16-variable data inputs  (after all necessary transformations) as described in Table 1 of the paper. The script also prepares LaTeX output to produce Table 1 of the paper.

The supplementary appendix also contains an alternative specification of our VAR models using variables in levels rather than differences. To prepare data for this specification, use `generate_freddataLEVELS.m` which generates `fredMD16levels-2021-04.csv`.

Copies of `fredMD16-2021-04.csv` and  `fredMD16levels-2021-04.csv` have been stored in the root folder of this repository for future use. (In case of data updates, these copies are not created automatically but must be made manually by the user.)

Plots of the transformed input data, as shown in Section I of the supplementary appendix, are produced by calling `chartData16.m` contained in the root folder of this repository.

## To estimate models

- To execute individual models for a given sample, use `doVAR.m` and select the desired model by setting the variable `modeltype` in the preamble to the script. 

- To launch out-of-sample runs for a given model, consider scripts called `goXXX.m`, where `XXX` refers to one of the model labels listed above. After computing the out-of-sample runs each of these `goXXX.m` scripts stores results for further post-processing in a `AAA-BBB-p12.mat` file where `AAA` reflects the name of the input data file (default: `fredMD16-2021-4`) and `BBB` refers to the model name and additional estimation options. By default, these mat files are stored in the current folder, but the user can also chose to move these manually to a different location. Based on these mat files, figures and tables for the paper can be created as described further below. Each of these driver files calls MCMC sampling routines named `mcmcXXX.m` where `XXX` corresponds to the model labels listed above.

- The main directory contains a bash script `gobatch.sh` that can be used to launch a sequence of multiple Matlab scripts from the shell. The Matlab scripts are executed in sequence *and in separate Matlab sessions*. Each Matlab session opens a parallel pool. For example, the shell command `sh gobatch.sh goVARSV.m goVARSVO.m` will launch a command line session of Matlab, start a parallel pool, and then execute `goVARSV.m`; once `goVARSV.m` has been executed, the Matlab sessions closes, a new one is reopened for execution of `goVARSVO.m`. (The shell script supports as many command line arguments as supported by bash and has been written for use on macOS and Linux.) Alternatively, Matlab scripts can, of course, equally be called interactively on the Matlab GUI’s command line.

## To produce figures and tables 
To collect model outputs, there are various scripts called `plotPredicitiveDensity*`, and `oos*` (with the former generating plots while the latter create tables). The preamble to each of these scripts specifies a variable `resultsdir` that is expected to point to a directory of `mat` files as produced by the `goVARXXX.m` routines described above.

Specifically, to collect predictive scores as compiled in Table 2, use `oosMVlogscores.m` (`oosMVlogscoresLevels.m` and `oosMVlogscoresSVar1` produce corresponding tables for VAR specifications in levels and AR(1) processes for SV, respectively, as shown in the supplementary appendix). 

To report relative RMSE and CRPS, as in Table 3, call `oosEvaluationTables.m` (and `oosEvaluationTablesSVAR1.m` for the case of SV-AR(1) models as reported in the supplementary appendix).

To produce the panels of Figure 1 (and additional output for other variables) use `compareSVOt.m`.

Figure 2 is produces by `compareSVpaths.m`.

The panels of Figures 3 and 4 are obtained from `plotPredictiveDensitiesPaper2021.m`. To obtain additional results for other variables, set the variable`doFullYlist` to true in the preamble of the script.
