module model
  use constants
  private
  public :: updateBH_sample, calc_nu, calc_sig, init_model, kill_model
  public :: param, param_low, param_upper, param_type, normdev

  double precision, dimension(:), allocatable :: param, param_low, param_upper
  integer, dimension(:), allocatable :: param_type
  
contains


  subroutine init_model(n, low, upper, ptype)
    implicit none

    integer, intent(in) :: n
    double precision, intent(in), dimension(n) :: low, upper
    integer, intent(in), dimension(n) :: ptype
    integer i
    
    
    allocate ( param(n) )
    allocate ( param_low(n) )
    allocate ( param_upper(n) )
    allocate ( param_type(n) )
    
    do i=1,n
       param_low(i)=low(i)
       param_upper(i)=upper(i)
       param_type(i) = ptype(i)
    enddo

  end subroutine init_model
      

  subroutine kill_model()

    implicit none
    
    if(ALLOCATED(param)) DEALLOCATE (param)
    if(ALLOCATED(param_low)) DEALLOCATE (param_low)
    if(ALLOCATED(param_upper)) DEALLOCATE (param_upper)
    if(ALLOCATED(param_type)) DEALLOCATE (param_type)

  end subroutine kill_model

  
  subroutine normdev(mySAMPLE, SIGMA)

    use random, only : random_normal

    implicit none

    TYPE(SAMPLE) :: mySample
    real, intent(in) :: SIGMA
    integer i
    

    do  i = 1, mySAMPLE%N
       mySAMPLE%RANDOM(i) = random_normal() * SIGMA
       !write(*,*) mySAMPLE%RANDOM(i)
    enddo

  end subroutine normdev

  subroutine updateBH_sample(mySample, A, B, SIGMA)

    implicit none
    
    real, intent(in) :: A, B, SIGMA
    TYPE(SAMPLE) :: mySample
    integer i
    !call normdev(mySAMPLE, SIGMA)
    mySample%LGMBH = A + B * (mySample%LGMS-10.0) + mySAMPLE%RANDOM

    !do  i = 1, mySAMPLE%N
    !   write(*,*) mySAMPLE%RANDOM(i)
    !enddo


    !write(*,*) mySample%LGMBH(1:10)
    !write(*,*) mySample%LGMS(1:10)
    !write(*,*) A, B
    
    
  end subroutine updateBH_sample

  
  subroutine cal_bol(mySample)

    implicit none
    
    TYPE(SAMPLE) :: mySample
    
    
    mySample%LGLBOL = log10(25.0) + mySample%LGLX
    
  end subroutine cal_bol

  
  subroutine cal_edd(mySample)
  
    implicit none
    
    TYPE(SAMPLE) :: mySample
       
    mySample%LGLEDD = mySample%LGLBOL - log10(1.26)- 38.0 - mySample%LGMBH
    
  end subroutine cal_edd

  
  
  subroutine calc_nu(mySample)

    implicit none

    TYPE(SAMPLE) :: mySample
    
    mySample%NUMIN = (1.0+mySample%Z) / ( 86400.0 * mySample%DTMAXOBS )
    mySample%NUMAX = (1.0+mySample%Z) / ( 86400.0 * mySample%DTMINOBS )

  end subroutine calc_nu


  subroutine calc_sig(mySample, alpha, beta)

    implicit none

    TYPE(SAMPLE) :: mySample
    real, optional :: alpha, beta

    call cal_bol(mySample)
    call cal_edd(mySample)
    !call ModelPSD(mySample)
    call ModelPSD(mySample, alpha, beta)
    
    mySample%SIG2 = mySample%PSDNORM * &
         (log(mySample%NUMAX/mySample%NUMIN) &
         - log( (mySample%NUB+mySample%NUMAX) &
         / (mySample%NUB+mySample%NUMIN) ) )
    mySample%SIG2 = mySample%SIG2 / 1.3 / 0.48**(1.1-1)
    
  end subroutine calc_sig


  subroutine ModelPSD(mySample, alpha, beta)
    implicit none
    TYPE(SAMPLE) :: mySample
    real, optional :: alpha, beta
    !if(present(alpha))then
    call Model3XX(mySample, alpha, beta)
    !else 
    !   call Model3(MySample)
    !endif
  end subroutine ModelPSD

  
  subroutine Model1(mySample)

    implicit none
    
    TYPE(SAMPLE) :: mySample
    !  PSD Model 1 of Paolillo+17                
    !  nub=580/(MBH/Msolar)s-1,
    !  nub * PSD(nub) = A / 2 -> A = 2 * nub * PSD(nub) = 0.02

    mySample%PSDNORM = 0.02
    mySample%NUB =log10(580.0) - mySample%LGMBH
    mySample%NUB = 10**mySample%NUB
  end subroutine Model1


  
  subroutine Model2(mySample)

    implicit none
    
    TYPE(SAMPLE) :: mySample

    ! PSD Model 2 of Paolillo+17                
    ! nub is a function of LGBOL and LEDD
    ! nub * PSD(nub) = A / 2 -> A = 2 * nub * PSD(nub) = 0.02
    
    mySample%PSDNORM = 0.02
    mySample%LGLBOL = mySample%LGLEDD + log10(1.26)+ 38.0 + mySample%LGMBH
    mySample%NUB = log10(200./86400.) + (mySample%LGLBOL-44.0) - 2.0*( mySample%LGMBH-6.0)
    mySample%NUB = 10**mySample%NUB
  end subroutine Model2


  subroutine Model3(mySample)

    implicit none
    
    TYPE(SAMPLE) :: mySample

    ! PSD Model 3 of Paolillo+17   
    ! nub=580/(MBH/Msolar)s-1,
    ! PSD normalisation is 3e-3 x LEdd^-0.8
    ! PSD(nub) = A/nub * (1+nub/nub)^-1->
    ! nub * PSD(nub) = A / 2 -> A = 2 * nub * PSD(nub) =  3e-3 x LEdd^-0.8
    
    mySample%PSDNORM = log10(3e-3)-0.8*mySample%LGLEDD
    mySample%PSDNORM =10**mySample%PSDNORM
    mySample%NUB =log10(580.0) - mySample%LGMBH
    mySample%NUB = 10**mySample%NUB
    
  end subroutine Model3


  subroutine Model3X(mySample, alpha)

    implicit none
    
    TYPE(SAMPLE) :: mySample
    real alpha

    ! variant of PSD Model 3 of Paolillo+17   
    ! nub=580/(MBH/Msolar)s-1,
    ! PSD normalisation is 3e-3 x LEdd^alpha
    ! where alpha a free parameter

    mySample%PSDNORM = log10(3e-3)+alpha*mySample%LGLEDD
    mySample%PSDNORM =10**mySample%PSDNORM
    mySample%NUB =log10(580.0) - mySample%LGMBH
    mySample%NUB = 10**mySample%NUB
    
  end subroutine Model3X


  subroutine Model3XX(mySample, alpha, beta)

    implicit none
    
    TYPE(SAMPLE) :: mySample
    real alpha, beta

    ! variant of PSD Model 3 of Paolillo+17   
    ! nub=580/(MBH/Msolar)s-1,
    ! PSD normalisation is beta x LEdd^alpha
    ! where alpha, beta free parameters

    mySample%PSDNORM = beta+alpha*mySample%LGLEDD
    mySample%PSDNORM =10**mySample%PSDNORM
    mySample%NUB =log10(580.0) - mySample%LGMBH
    mySample%NUB = 10**mySample%NUB
    
  end subroutine Model3XX


  subroutine Model4(mySample)

    implicit none
    
    TYPE(SAMPLE) :: mySample

    ! Model 4 of  Paolillo+17   

    mySample%PSDNORM = log10(3e-3)-0.8*mySample%LGLEDD
    mySample%PSDNORM =10**mySample%PSDNORM
    mySample%LGLBOL = mySample%LGLEDD + log10(1.26)+ 38.0 + mySample%LGMBH
    mySample%NUB = log10(200./86400.) + (mySample%LGLBOL-44.0) - 2.0*( mySample%LGMBH-6.0)
    mySample%NUB = 10**mySample%NUB
  end subroutine Model4



  
end module model
  
