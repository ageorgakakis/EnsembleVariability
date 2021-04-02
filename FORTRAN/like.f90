 module like
	
use params
use utils1
implicit none
      
contains
      
      
!=======================================================================

subroutine slikelihood(Cube,slhood)
  use constants
  use model
  use data
  use priors

  implicit none

  double precision Cube(nest_nPar),slhood  
  integer i,j
  real ymodel, a1
  double precision X, Y, MEAN, SIGMA

  slhood  = 0d0
  
  do i = 1, sdim
     if(param_type(i)==0)then
        param(i) = param_low(i) + ( param_upper(i) - param_low(i) ) * Cube(i)
     else if(param_type(i)==1)then
        MEAN = param_low(i)
        SIGMA = param_upper(i)
        X =  Cube(i)
        call invert_normal_prior(X, Y, MEAN, SIGMA)        
        param(i) = Y
     endif
     !write(*,*)  "PARAM:",i, param(i), param_upper(i),  param_low(i), Cube(i)
  enddo
  !write(30,*) param(3)

  
  do i=1,N_MASTER_SAMPLE
     
     call updateBH_sample(MASTER_SAMPLE(i), real(param(1)), real(param(2)), 0.5)     
     call calc_nu(MASTER_SAMPLE(i))
     call calc_sig(MASTER_SAMPLE(i), real(param(3)), real(param(4)) )

     a1 = SUM(MASTER_SAMPLE(i)%WEIGHT)
     ymodel = SUM(MASTER_SAMPLE(i)%SIG2 * MASTER_SAMPLE(i)%WEIGHT)/a1
     
     slhood = slhood - 0.5 *  dble( (MASTER_SAMPLE(i)%sigma2 - ymodel )**2 / MASTER_SAMPLE(i)%esigma2 / MASTER_SAMPLE(i)%esigma2 )


  enddo

  
end subroutine slikelihood




end module like






