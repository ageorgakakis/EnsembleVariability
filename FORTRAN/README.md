Collection of FORTRAN90 routines to fit the empirical ensemble variance model to observations and infer model parameters in a Bayesian framework. The observational datapoints are the ensemble excess variance measurements in X-ray luminosity bins presented by Paolillo et al. (2017; their Figure 5). The empirical model components that are allowed to vary are (i) the Black-Hole Mass vs Stellar Mass relation and (ii) the variability Power Spectrum Density.

The [MultiNest library](https://github.com/farhanferoz/MultiNest) is used for Bayesian Inference and should be installed. Additionally the fitsio routines are used to read in the data and the mock AGN catalogue. This library should also be present. For the compilation, it is assumed that the directory of the MultiNest modules and library are defined, e.g.

`export multimod = /util/MultiNest/MultiNest`
`export multilib = /util/MultiNest/MultiNest`

Simarly for the fitsio library

`export fitsio = /util/cfisio/lib`

The compilation script using `mpif90` compiler could look like


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
