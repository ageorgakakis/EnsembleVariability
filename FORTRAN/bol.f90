module bolometric
  
  use constants
  
  implicit none

  private
  public :: GLBOL, findBOL, initBOL, killBOL
  

  ! global bolometric correction structure
  ! that strores the bolometric corections
  ! parameters and lookup tables
  ! as defined in initBOL
  TYPE(BOL) GLBOL

contains
    
  
  subroutine findBOL(LGLX, myBOL, LGLBOL)

    implicit none

    TYPE(BOL), intent(in) :: myBOL
    real, intent(in) :: LGLX
    real, intent(out) :: LGLBOL
    INTEGER N, FLAG


    N=0
    FLAG=0
    DO WHILE (FLAG==0 .AND. N<SIZE(myBOL%LGLBOL) )
       N=N+1
       if(LGLX<myBOL%LGLX(N))FLAG=1
       !write(*,*) "H",LGLX, N, myBOL%LGLX(N), FLAG

    ENDDO
        
    N=N-1
    IF(N<=0)N=1

    
    LGLBOL = myBOL%LGLBOL(N) +(myBOL%LGLBOL(N+1) - myBOL%LGLBOL(N)) &
         / (myBOL%LGLX(N+1) - myBOL%LGLX(N) ) * ( LGLX - myBOL%LGLX(N))
    

  end subroutine findBOL


  subroutine initBOL(myBOL)

    implicit  none
    
    TYPE(BOL), intent(out) :: myBOL
    
    myBOL%LGLMIN = 6
    myBOL%LGLMAX = 15
    myBOL%DLGL = 0.1
    myBOL%N = ( myBOL%LGLMAX - myBOL%LGLMIN ) /  myBOL%DLGL
    
    allocate ( myBOL%LGLBOL(myBOL%N) )
    allocate ( myBOL%LGLX(myBOL%N ) )
    allocate ( myBOL%K(myBOL%N) )
    
    ! using the Duras+20 corrections
    call calcduras(myBOL)
    !call callconstant(myBOL)
  end subroutine initBOL
  
  
  subroutine calcduras(myBOL)

    ! bolometric corrections from Duras et al. 2020
    ! https://ui.adsabs.harvard.edu/abs/2020A%26A...636A..73D/abstract
    ! for their total AGN sample
    ! The bolometric correction is for the 2-10keV X-ray band

    implicit none
    TYPE(BOL) myBOL
    real a, b, c, LGLBOL
    integer I

    a = 10.96;b=11.93;c=17.79
    !write(*,*) myBOL%LGLMIN,  myBOL%DLGL, SIZE(myBOL%K),SIZE(myBOL%LGLX)
    do I=1,myBOL%N
       LGLBOL = myBOL%LGLMIN + (i-1) * myBOL%DLGL
       myBOL%K(I)  = a * (1+ (LGLBOL/b)**c)
       myBOL%LGLX(I) = LGLBOL - LOG10(myBOL%K(I)) + LOG10(3.826)+33.0
       myBOL%LGLBOL(I) = LGLBOL +  LOG10(3.826)+33.0
       !write(*,*) myBOL%LGLBOL(I), myBOL%LGLX(I), myBOL%K(I)
    enddo
       
  end subroutine calcduras



  subroutine calcconst(myBOL)

    ! simple bolometric correction
    ! using a constant scaling factor of 25
    ! the correction is for the 2-10keV band
    ! this is implemented as a lookup table
    ! for consistency

    implicit none
    TYPE(BOL) myBOL
    real c, LGLBOL
    integer I

    c=25.0

    !write(*,*) myBOL%LGLMIN,  myBOL%DLGL, SIZE(myBOL%K),SIZE(myBOL%LGLX)
    do I=1,myBOL%N
       LGLBOL = myBOL%LGLMIN + (i-1) * myBOL%DLGL
       myBOL%K(I)  = c
       myBOL%LGLX(I) = LGLBOL - LOG10(myBOL%K(I)) + LOG10(3.826)+33.0
       myBOL%LGLBOL(I) = LGLBOL +  LOG10(3.826)+33.0

    enddo
       
  end subroutine calcconst



  subroutine killBOL(myBOL)
    implicit none
    
    TYPE(BOL), intent(out) :: myBOL
    
    if(ALLOCATED(myBOL%LGLBOL)) DEALLOCATE (myBOL%LGLBOL)
    if(ALLOCATED(myBOL%K)) DEALLOCATE (myBOL%K)
    if(ALLOCATED(myBOL%LGLX)) DEALLOCATE (myBOL%LGLX)
    
  end subroutine killBOL


end module bolometric
