module constants

  implicit none
  
  type :: SAMPLE

     integer :: N
     
     real, dimension (:), allocatable :: Z
     real, dimension (:), allocatable :: LGMBH
     real, dimension (:), allocatable :: LGMS
     real, dimension (:), allocatable :: LGLEDD
     real, dimension (:), allocatable :: LGLX
     real, dimension (:), allocatable :: LGLBOL     
     real, dimension (:), allocatable :: SIG2
     real, dimension (:), allocatable :: PSDNORM 
     real, dimension (:), allocatable :: NUB 
     real, dimension (:), allocatable :: NUMIN
     real, dimension (:), allocatable :: NUMAX
     real, dimension (:), allocatable :: WEIGHT
     real, dimension (:), allocatable :: RANDOM
     real LGLXMIN, LGLXMAX, LGLXMEAN
     real ZMIN, ZMAX     
     real DTMINOBS, DTMAXOBS
     real :: sigma2
     real :: esigma2
    
  end type SAMPLE



  type :: OBS

     real :: ZMIN, ZMAX
     real :: LGLXMIN, LGLXMAX, LGLXMEAN
     real :: DTMIN, DTMAX
     real :: sigma2
     real :: esigma2
         
  end type OBS
  
  type :: BOL

     real, dimension(:), allocatable :: LGLBOL
     real, dimension(:), allocatable :: LGLX 
     real, dimension(:), allocatable :: K
     real LGLMIN
     real LGLMAX
     real DLGL
     integer N
  end type BOL

end module constants
