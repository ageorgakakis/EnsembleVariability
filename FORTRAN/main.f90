program main

  use params
  use constants
  use data
  use like
  use nestwrapper
  use model
  use recon
  
  implicit none
  integer i,j

  ! PRIORS ARE DEFINE HERE FOR
  ! THE MODEL PARAMETERS
  double precision low(sdim)
  data low /-1.0, 0.0, -0.8, -2.522/
  double precision upper(sdim)
  data upper /9.0, 3.0, 0.15, 0.30/
  integer ptype(sdim)

  ! prior type 0: flat between "low" and "upper" range
  ! prior type 0: Gaussian with low==mean, upper==sigma
  data ptype /0, 0, 1, 1/
  double precision cube(sdim)

  
  call init_model(sdim, low, upper, ptype)


  ! each of these files corresponds to one of 
  ! datapoints of the excess variance measurements
  ! in the 7Ms-CDFS presented by Paolillo et al. 2017
  call init_sample_file(MASTER_SAMPLE(1), "SAMPLE_Z1_0.4_Z2_4.0_DMIN_0.25_DMAX_6205.0_LGLX_42.5.fits")
  call init_sample_file(MASTER_SAMPLE(2), "SAMPLE_Z1_0.4_Z2_4.0_DMIN_0.25_DMAX_6205.0_LGLX_43.1.fits")
  call init_sample_file(MASTER_SAMPLE(3), "SAMPLE_Z1_0.4_Z2_4.0_DMIN_0.25_DMAX_6205.0_LGLX_43.5.fits")
  call init_sample_file(MASTER_SAMPLE(4), "SAMPLE_Z1_0.4_Z2_4.0_DMIN_0.25_DMAX_6205.0_LGLX_43.6.fits")
  call init_sample_file(MASTER_SAMPLE(5), "SAMPLE_Z1_0.4_Z2_4.0_DMIN_0.25_DMAX_6205.0_LGLX_43.9.fits")
  call init_sample_file(MASTER_SAMPLE(6), "SAMPLE_Z1_0.4_Z2_4.0_DMIN_0.25_DMAX_6205.0_LGLX_44.1.fits")
  call init_sample_file(MASTER_SAMPLE(7), "SAMPLE_Z1_0.4_Z2_4.0_DMIN_0.25_DMAX_6205.0_LGLX_44.3.fits")
  call init_sample_file(MASTER_SAMPLE(8), "SAMPLE_Z1_0.4_Z2_4.0_DMIN_0.25_DMAX_6205.0_LGLX_44.9.fits")
  N_MASTER_SAMPLE=8
  

  ! runs multinest
  nest_pWrap=0
  call nest_Sample

  ! reads the chains and produces outputs for plotting
  call readchains()
  do i=1,N_MASTER_SAMPLE
     call kill_sample(MASTER_SAMPLE(i))
  enddo
  call kill_model()
   
end program main
