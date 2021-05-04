module data

  use constants

  implicit none
    
  
  private
  public :: init_sample, init_sample_file, kill_sample, init_obs
  public :: MASTER_SAMPLE, N_MASTER_SAMPLE
  
  TYPE(SAMPLE),  dimension (11) :: MASTER_SAMPLE
  INTEGER N_MASTER_SAMPLE

contains


  subroutine init_sample(mySample)
    implicit none
    TYPE(SAMPLE), intent(out) :: mySample
    call kill_sample(mySample)
  end subroutine init_sample

  subroutine init_sample_file(mySample, filename)
    implicit none
    TYPE(SAMPLE), intent(out) :: mySample
    character(len=*), intent(in) :: filename
    call kill_sample(mySample)
    call read_sample(mySample, filename)        
  end subroutine init_sample_file


  subroutine init_arrays_sample(mySample, N)
    implicit none
    TYPE(SAMPLE) :: mySample
    integer, intent(in) :: N

    ALLOCATE (mySample%Z(N))
    ALLOCATE (mySample%LGMBH(N))
    ALLOCATE (mySample%LGMS(N))
    ALLOCATE (mySample%LGLEDD(N))
    ALLOCATE (mySample%LGLX(N))
    ALLOCATE (mySample%LGLBOL(N))
    ALLOCATE (mySample%PSDNORM(N))
    ALLOCATE (mySample%NUB(N))
    ALLOCATE (mySample%NUMIN(N))
    ALLOCATE (mySample%NUMAX(N))
    ALLOCATE (mySample%SIG2(N))
    ALLOCATE (mySample%WEIGHT(N))
    ALLOCATE (mySample%RANDOM(N))
    
  end subroutine init_arrays_sample
  
  
  subroutine kill_sample(mySample)

    implicit none

    TYPE(SAMPLE) :: mySample

    if(ALLOCATED(mySample%Z)) DEALLOCATE (mySample%Z)
    if(ALLOCATED(mySample%LGMBH)) DEALLOCATE (mySample%LGMBH)
    if(ALLOCATED(mySample%LGMS)) DEALLOCATE (mySample%LGMS)
    if(ALLOCATED(mySample%LGLEDD)) DEALLOCATE (mySample%LGLEDD)
    if(ALLOCATED(mySample%LGLX)) DEALLOCATE (mySample%LGLX)
    if(ALLOCATED(mySample%SIG2)) DEALLOCATE (mySample%SIG2)      
    if(ALLOCATED(mySample%LGLBOL)) DEALLOCATE (mySample%LGLBOL)      
    if(ALLOCATED(mySample%PSDNORM)) DEALLOCATE (mySample%PSDNORM)      
    if(ALLOCATED(mySample%NUB)) DEALLOCATE (mySample%NUB)      
    if(ALLOCATED(mySample%NUMIN)) DEALLOCATE (mySample%NUMIN)      
    if(ALLOCATED(mySample%NUMAX)) DEALLOCATE (mySample%NUMAX)      
    if(ALLOCATED(mySample%WEIGHT)) DEALLOCATE (mySample%WEIGHT)      
    if(ALLOCATED(mySample%RANDOM)) DEALLOCATE (mySample%RANDOM)      
    
  end subroutine kill_sample


  

  subroutine read_sample(mySample, filename)

    !use fitstools
   
    
    implicit none
    TYPE(SAMPLE), intent(out) :: mySample
    character(len=*), intent(in) :: filename

    
    integer :: IROW, NROWS
    integer :: colnumZ, colnumLX, colnumMS, colnumW
    integer :: unit, readwrite, blocksize, nelems, felems, anyf, status
    integer :: hdutype, extver
    real    :: nullval
    logical casesen
    character :: coltemplate*10, extname*10, comment*10
    real a1

    write(*,*) filename
    
    ! open the FITS file, with read-only access.  The returned BLOCKSIZE
    ! parameter is obsolete and should be ignored. 
    unit=0
    readwrite=0
    status=0
    call ftopen(unit,filename,readwrite,blocksize,status)

    hdutype=2
    extver=0    
    extname="SAMPLE"
    call FTMNHD(unit, hdutype, extname, extver, status)
    
    ! read the header
    !call FTGKYJ(unit,"NAXIS2",mySample%N,comment, status)
    !Get the number of rows or columns in the current FITS table.
    call FTGNRW(unit, NROWS, status)
    mySample%N = NROWS

    call FTGKYE(unit,"ZMIN",mySample%ZMIN,comment, status)
    call FTGKYE(unit,"ZMAX",mySample%ZMAX,comment, status)
    call FTGKYE(unit,"LGLXMIN",mySample%LGLXMIN,comment, status)
    call FTGKYE(unit,"LGLXMAX",mySample%LGLXMAX,comment, status)
    call FTGKYE(unit,"DTMINOBS",mySample%DTMINOBS,comment, status)
    call FTGKYE(unit,"DTMAXOBS",mySample%DTMAXOBS,comment, status)
    call FTGKYE(unit,"SIG2",mySample%sigma2,comment, status)
    call FTGKYE(unit,"ESIG2",mySample%esigma2,comment, status)

    !write(*,*) mySample%DTMINOBS, mySample%DTMAXOBS, mySample%sigma2, mySample%esigma2
    
    ! init arrays
    call init_arrays_sample(mySample, NROWS)
    !write(*,*) SIZE(mySample%Z), SIZE(mySample%LGLX),SIZE(mySample%LGMS)
    
    !Get the table column number (and name) of the column whose name matches an input template name.
    casesen = .false.
    coltemplate = "Z"
    call FTGCNO(unit,casesen,coltemplate, colnumZ,status)
    coltemplate = "LGLX"
    call FTGCNO(unit,casesen,coltemplate, colnumLX,status)
    coltemplate = "LGM"
    call FTGCNO(unit,casesen,coltemplate, colnumMS,status)
    coltemplate = "WEIGHT"
    call FTGCNO(unit,casesen,coltemplate, colnumW,status)

    felems=1
    nelems=1
    do irow=1, mySample%N     
       call ftgcve(unit,colnumZ,irow,felems,nelems,nullval,a1, &
            anyf,status)
       mySample%Z(irow)=a1
       call ftgcve(unit,colnumLX,irow,felems,nelems,nullval,a1, &
            anyf,status)
       mySample%LGLX(irow)=a1
       call ftgcve(unit,colnumMS,irow,felems,nelems,nullval,a1, &
            anyf,status)
       mySample%LGMS(irow)=a1
       call ftgcve(unit,colnumW,irow,felems,nelems,nullval,a1, &
            anyf,status)
       mySample%WEIGHT(irow)=a1
 
       !write(*,*) irow,mySample%Z(irow), mySample%LGLX(irow),mySample%LGMS(irow), mySample%WEIGHT(irow)
    enddo
    
    mySAMPLE%RANDOM=0.0
    
    call cal_bol(mySample)

    call ftclos(unit, status)
    call ftfiou(unit, status) 
  end subroutine read_sample


  subroutine cal_bol(mySample)

    use bolometric

    implicit none
 
       
    TYPE(SAMPLE) :: mySample
    INTEGER I 

    do I=1,SIZE(mySample%LGLX)
       call findBOL(mySample%LGLX(I), GLBOL, mySample%LGLBOL(I))
       !write(*,*) I, mySample%LGLX(I), mySample%LGLBOL(I)
    enddo
    !mySample%LGLBOL = log10(25.0) + mySample%LGLX
    
  end subroutine cal_bol

 

  
  subroutine init_obs(myObs, ZMIN, ZMAX, LGLXMIN, LGLXMAX, LGLXMEAN, DTMIN, DTMAX, SIGMA2, ESIGMA2)

    implicit none
    
    TYPE(OBS), intent(out) :: myObs
    real, intent(in) :: ZMIN, ZMAX
    real, intent(in) :: LGLXMIN, LGLXMAX, LGLXMEAN
    real, intent(in) :: DTMIN, DTMAX
    real, intent(in) :: SIGMA2
    real, intent(in) :: ESIGMA2

    myObs%ZMIN=ZMIN
    myObs%ZMAX=ZMAX
    myObs%DTMIN=DTMIN
    myObs%DTMAX=DTMAX
    myObs%LGLXMIN=LGLXMIN
    myObs%LGLXMAX=LGLXMAX
    myObs%LGLXMEAN=LGLXMEAN
    myObs%SIGMA2=SIGMA2
    myObs%ESIGMA2=ESIGMA2
    
  end subroutine init_obs

end module data
  
