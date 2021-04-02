
========================================
AGN Ensemble Variability Empirical Model 
========================================


This project provides routines to model the ensemble excess variance of AGN using as input emprical relations and quantities. These inlcude the stellar mass function of galaxies, the specific accretion-rate distribution of AGN, the black-hole mass vs stellar-mass scaling relation and the Power Spectrum Density of the AGN flux variations. The flow-chart below shows how these components are combined to make predictions on the ensemble excess variance of AGN populations.

.. figure:: model.png

A galaxy sample drawn from the stellar mass function at a given redshift (Panel 1) is seeded with AGN specific accretion-rates using observationally-derived probability distribution functions (Panel 2). This produces a sample of mock AGN (red dots of Panel 3), each of which has been assigned an X-ray luminosity, a host-galaxy stellar mass and a redshift. A parametrisation of the Black-Hole Mass vs Stellar Mass relation (Panel 5) is used to assign black holes to mock AGN and hence, Eddington ratios. The dependence of the AGN variability Power Spectrum Density (PSD; Panel 4) on Black Hole Mass and Eddington ratio is then used to assign an excess variance by integrating the corresponding PSD. The average excess variance of the population binned in luminosity and redshift intervals (Panel 6) can then be directly compared with observational results.
