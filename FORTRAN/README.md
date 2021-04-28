## Fortran90 Routines for Bayesian Inference

Collection of FORTRAN90 routines to fit the empirical ensemble variance model to observations and infer model parameters in a Bayesian framework. The observational datapoints are the ensemble excess variance measurements in X-ray luminosity bins presented by [Paolillo et al. (2017; their Figure 5)](https://ui.adsabs.harvard.edu/abs/2017MNRAS.471.4398P/abstract). The empirical model components that are allowed to vary are (i) the Black-Hole Mass vs Stellar Mass relation and (ii) the variability Power Spectrum Density.

The [MultiNest library](https://github.com/farhanferoz/MultiNest) is used for Bayesian Inference and should be installed. Additionally the fitsio routines are used to read in the data and the mock AGN catalogue. This library should also be present. For the compilation, the ditrectory of the MultiNest modules and library are defined, e.g.:

`export multimod = /util/MultiNest/MultiNest`

`export multilib = /util/MultiNest/MultiNest`

Simarly for the fitsio library:

`export fitsio = /util/cfisio/lib`

The compilation script using `mpif90` compiler could look like:

`mpif90 -lmpi -O3 -DMPI -fPIC -ffree-line-length-none -c random.f90`

`mpif90 -lmpi -O3 -DMPI -fPIC -ffree-line-length-none -c priors.f90`

`mpif90 -lmpi -O3 -DMPI -fPIC -ffree-line-length-none -c params.f90`

`mpif90 -lmpi -O3 -DMPI -fPIC -ffree-line-length-none -c constants.f90`

`mpif90 -lmpi -O3 -DMPI -fPIC -ffree-line-length-none -c data.f90 -lcfitsio`

`mpif90 -lmpi -O3 -DMPI -fPIC -ffree-line-length-none -c model.f90`

`mpif90 -lmpi -O3 -DMPI -fPIC -ffree-line-length-none -c recon.f90`

`mpif90 -lmpi -O3 -DMPI -fPIC -ffree-line-length-none -c like.f90  -o like.o -I$multimod`

`mpif90 -lmpi -O3 -DMPI -fPIC -ffree-line-length-none -O3 -DMPI -fPIC -ffree-line -length-none -lmpi_f90 -I$multimod -c nestwrap.f90 -o nestwrap.o`

`mpif90 -lmpi -O3 -DMPI -fPIC -ffree-line-length-none -c main.f90 -lcfitsio -I$multimod`

The compilation line below assumes that the MultiNest library is `libnest3` at `$multilib`:

`mpif90 -o a.out params.o priors.o constants.o data.o random.o model.o recon.o like.o nestwrap.o main.o -lmpi -O3 -DMPI -fPIC -ffree-line-length-none -I$multimod -L/usr/lib -L$multilib -lnest3 -llapack -L$fitsio -lcfitsio`

The `main.f90` defines the free parameters of the problem. Currently these are two for the black-hole/stellar mass relation and two for the PSD. The priors for each of these parameters are definedin the lines:

`data low /-1.0, 0.0, -0.8, -2.522/`

`data upper /9.0, 3.0, 0.15, 0.30/`

`data ptype /0, 0, 1, 1/`


where `ptype` is a table that defines the prior type: flat (`0`) within a given range, or gaussian (`1`). In the case of a flat prior the tables `low` and `upper` list the lower and upper limit of the interval. In the case of a gaussian prior the `low` table entry corresponds to the mean and the `upper` table-entry to the sigma.

For each data point of Figure 5 of [Paolillo et al. (2017)](https://ui.adsabs.harvard.edu/abs/2017MNRAS.471.4398P/abstract) there is a separate mock AGN catalogue selected in the same luminosity interval. These mock AGN catalogues are read in as fits binary tables in the lines that look like:

`call init_sample_file(MASTER_SAMPLE(1), "SAMPLE_Z1_0.4_Z2_4.0_DMIN_0.25_DMAX_6205.0_LGLX_42.5.fits")`

The input files can be generated using the python script prepFORT.py

the MultiNest bayesian inference routines are called in the lines:


`  nest_pWrap=0
  call nest_Sample`


the following line reads the MCMC chains produced by MultiNest to determine the median of parameters and their uncertainties, estimate confidence intervals for e.g. the black-hole mass vs stellar mass relation, etc.  

`call readchains()`

The model that is fit to the observations is defined in `model.f90`:

`subroutine Model3XX(mySample, alpha, beta)`

