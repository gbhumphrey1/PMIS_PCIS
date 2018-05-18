PMIS_PCIS
=========

R scripts used to implement the PMIS and PCIS input variable selection (IVS) algorithms as part of the IVS4EM project described in Galelli et. al. (2014). The PMIS algorithm is a filter IVS method developed by Sharma (2000) and later modified by Bowden et al. (2005) and May et al. (2008), where the relevance of potential inputs is evaluated based on the mutual information (MI) between each input variable and the output. The PCIS algorithm (May et al., 2008) is also a filter IVS method, where input relevance is based on partial correlation analysis. Further details of this particular implementation of these algorithms can be found in:

Galelli S., Humphrey G.B., Maier H.R., Castelletti A., Dandy G.C. and Gibbs M.S. (2014)  An evaluation framework for input variable selection algorithms for environmental data-driven models, *Environmental Modelling and Software*, 62, 33-51, DOI: 10.1016/j.envsoft.2014.08.015. ([Link to Paper](http://www.sciencedirect.com/science/article/pii/S1364815214002394))

The purpose of the IVS4EM project is to support a comprehensive framework for the testing and evaluation of IVS algorithms, through the sharing of algorithms (open source code), datasets, and evalution criteria. ([Link to IVS4EM Website](http://ivs4em.deib.polimi.it/))


Contents:
* `PMI_PCIS.R`: code to implement the PMIS and PCIS algorithms.
* `PMIS_run.R`: run the PMIS algorithm to select the most relevent inputs for a given set of input data.
* `PCIS_run.R`: run the PCIS algorithm to select the most relevant inputs for a given set of input data.
* `inp_dat.csv`: an example input data file. Column 1 contains an array of data labels or IDs (e.g. dates on which data were recorded); columns 2 to P+1 contain the P candidate input variables; and column P+2 contains the response variable, while the rows are data points. The first row contains the variable names.

To run the PMIS algorithm, the following command should be used:

`R --args [`*filename*`] [`*out_dir*`] < PMI_run.R`

where *filename* is the name of the name of the input data file (including path) and *out_dir* is the name of the output directory (i.e. the directory to which results will be written. This name should *NOT* include the whole path). The PCIS algorithm is run similarly.


Sharma, A., 2000. Seasonal to interannual rainfall probabilistic forecasts for improved water supply management: Part 1 - a strategy for system predictor identification. *Journal of Hydrology* 239, 232-239.

Bowden, G.J., Maier, H.R., Dandy, G.C., 2005. Input determination for neural network models in water resources applications. Part 1. Background and methodology. *Journal of Hydrology* 301, 75-92.

May, R.J., Maier, H.R., Dandy, G.C., Fernando, T.M.K.G., 2008. Nonlinear variable selection for artificial neural networks using partial mutual information. *Environmental Modelling & Software* 23, 1312-1326.

Copyright 2014 Greer Humphrey.
