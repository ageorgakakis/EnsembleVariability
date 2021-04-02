module recon

  use constants
  use params
  use model
  use data
  use priors

  implicit none

contains 


subroutine readchains()
  
  implicit none
  
  character (len=100) filename
  real, dimension(:,:), allocatable :: cube
  real, dimension(:), allocatable :: LGBH
  real, dimension(:), allocatable :: XSV
  real, dimension(:,:), allocatable :: XSVQ
  real, dimension(:), allocatable :: temp1
  real, allocatable :: masked(:)
  integer, dimension(:), allocatable :: ORDER
  integer i, j, nchain, stat, i1,i2
  real redshift, lgl, lglmin, lglmax, lgm, dlgm
  double precision X, Y, MEAN, SIGMA
  TYPE(SAMPLE) FSAMPLE


  filename = trim(nest_root)//"post_equal_weights.dat"

  open(1,file=filename,status='old')
  nchain = 0
  do j=1,100000   
     READ(1, *, iostat=stat)lgl
     if(stat/=0)exit
     !write(*,*) lgl
     nchain=nchain+1
  enddo
  close(1)

  ! prepare arrays
  allocate(cube(nchain, sdim))
  allocate(temp1(nchain))
  allocate(LGBH(nchain))
  allocate(order(nchain))
  allocate(XSV(nchain))
  allocate(XSVQ(5,N_MASTER_SAMPLE))
  


  !write(*,*) nchain
  open(1,file=filename,status='old')

  ! read the chains in the MultiNest 0-1 space
  stat=0
  do i=1,nchain  
     READ(1, *)(cube(i,j), j=1,sdim)
  enddo
  close(1)

  ! converts parameters from chains from MultiNest 0-1 space
  ! to physical space and estimates median parameters and 
  ! and confidence intervals (10, 16, 84, 90%)
  ! OUTPUT fort.89  
  do i = 1, sdim
     if(param_type(i)==0)then
        cube(:,i) = param_low(i) + ( param_upper(i) - param_low(i) ) * cube(:,i)
     else if(param_type(i)==1)then        
        MEAN = param_low(i)
        SIGMA = param_upper(i)
        do j=1,SIZE(cube(:,i))
           X =  DBLE(Cube(j,i))
           call invert_normal_prior(X, Y, MEAN, SIGMA)        
           cube(j,i) = REAL(Y)
        enddo
     endif

     temp1 = cube(:,i)
     call quick_sort(temp1, ORDER) 

     !write(*,*)temp1(1:10)
     write(89,*)i,temp1(int(nchain*0.5)), &
          temp1(int(nchain*0.1)), temp1(int(nchain*0.16)), &
          temp1(int(nchain*0.84)), temp1(int(nchain*0.90))
  enddo

  ! produce full chain with physical parameters
  ! OUTPUT fort.90
  do i=1,nchain  
     WRITE(90, *)(cube(i,j), j=1,sdim)
  enddo
  write(*,*) "Finished writting out physical params"
  !write out physical params
  
  
  ! M-MBH relations based on the chains
  ! for each luminosity bin estimates median and
  ! percentiles: 5, 10, 16, 84, 90, 95%
  ! OUTPUT FILE is fort.88
  dlgm=0.1
  i1=int(8./dlgm)
  i2=int(14./dlgm)
  do i=i1,i2,1
     lgm=(i-1)*dlgm
     
     LGBH = cube(:,1)+cube(:,2)*(lgm-10)
     temp1 = LGBH
     call quick_sort(temp1, ORDER)      
     write(88,*) lgm,temp1(int(nchain*0.5)),&
          temp1(int(nchain*0.16)),temp1(int(nchain*0.84)),&
          temp1(int(nchain*0.1)),temp1(int(nchain*0.90)),&
          temp1(int(nchain*0.05)),temp1(int(nchain*0.95))
  enddo
  write(*,*) "Finished writting out Mstar-MBH relation"
  ! M-MBH relations

!  do j=1,N_MASTER_SAMPLE  
!     do i=1,nchain     
!        call updateBH_sample(MASTER_SAMPLE(j), real(cube(i,1)), real(cube(i,2)), 0.5)
!        call calc_nu(MASTER_SAMPLE(j))
!        call calc_sig(MASTER_SAMPLE(j), real(cube(i,3)), real(cube(i,4)) )
!        
!        XSV(i) = SUM(MASTER_SAMPLE(j)%SIG2 * MASTER_SAMPLE(j)%WEIGHT) & 
!             /  SUM(MASTER_SAMPLE(j)%WEIGHT)
!        !write(*,*) i, XSV(i),  real(cube(i,3)), real(cube(i,4)), MASTER_SAMPLE(j)%sigma2
!     enddo
!     temp1 = XSV
!     call quick_sort(temp1, ORDER) 
!     XSVQ(1,j)=temp1(int(nchain*0.5))
!     XSVQ(2,j)=temp1(int(nchain*0.1))
!     XSVQ(3,j)=temp1(int(nchain*0.16))
!     XSVQ(4,j)=temp1(int(nchain*0.84))
!     XSVQ(5,j)=temp1(int(nchain*0.90))
!     !write(*,*) (XSVQ(i,j), i=1,5)
!  enddo


!  do j=1,N_MASTER_SAMPLE  
!     write(91,*)MASTER_SAMPLE(j)%LGLXMIN,MASTER_SAMPLE(j)%LGLXMAX, MASTER_SAMPLE(j)%sigma2,  MASTER_SAMPLE(j)%esigma2
!     write(92,*)MASTER_SAMPLE(j)%LGLXMIN,MASTER_SAMPLE(j)%LGLXMAX, (XSVQ(i,j), i=1,5)
!  enddo
!  deallocate(XSVQ)


  ! produces ensemble excess variance as a function of lumonsity
  ! takes long time to run!!!!
  ! output file fort.93
  call init_sample_file(FSAMPLE, "SAMPLE_Z1_0.4_Z2_4.0_DMIN_0.25_DMAX_6205.0_LGLX_40.0.fits")
  call normdev(FSAMPLE, 0.5)
  !write(*,*) "FSAMPLE SIZE:",FSAMPLE%N
  
  !stop
  dlgm=0.5
  i1=int(42/dlgm)
  i2=int(47./dlgm)
  
  allocate(XSVQ(5,i2+1))
  write(*,*) SHAPE(XSVQ)
  allocate(masked(FSAMPLE%N))
  do j=i1,i2,1
     lgm=(j-1)*dlgm
     write(*,*) lgm
     do i=1,nchain
        !write(*,*) i,SHAPE(XSV), shape(temp1)
        call updateBH_sample(FSAMPLE, real(cube(i,1)), real(cube(i,2)), 0.5)
        call calc_nu(FSAMPLE)
        call calc_sig(FSAMPLE, real(cube(i,3)), real(cube(i,4)) )
        
        !write(*,*) FSAMPLE%SIG2(1:10) 
        !write(*,*) FSAMPLE%LGLX(1:10)
        

        masked = MERGE(FSAMPLE%WEIGHT, 0., FSAMPLE%LGLX<lgm+dlgm/2.0.and.FSAMPLE%LGLX>=lgm-dlgm/2)
        XSV(i)=0.0
        !write(*,*) sum(masked)
        if(sum(masked)>0)then
           XSV(i) = SUM(FSAMPLE%SIG2 * masked) & 
                /  SUM(masked)
           endif
        !write(*,*) i, XSV(i),SHAPE(masked)
     enddo

     temp1 = XSV
     !write(*,*) temp1
     call quick_sort(temp1, ORDER) 
     !write(*,*) "HERE"
     XSVQ(1,j)=temp1(int(nchain*0.5))
     XSVQ(2,j)=temp1(int(nchain*0.05))
     XSVQ(3,j)=temp1(int(nchain*0.16))
     XSVQ(4,j)=temp1(int(nchain*0.84))
     XSVQ(5,j)=temp1(int(nchain*0.95))
     write(93,*) lgm+log10(1.503),(XSVQ(i,j), i=1,5)
     !write(93,*) lgm,(XSV(1000), i=1,5)
  enddo



  call kill_sample(FSAMPLE)  
  deallocate(masked)
  deallocate(cube)
  deallocate(temp1)
  deallocate(LGBH)
  deallocate(order)
  deallocate(XSV)
  deallocate(XSVQ)

end subroutine readchains




RECURSIVE SUBROUTINE quick_sort(list, order)

  ! Quick sort routine from:
  ! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
  ! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
  ! Modified by Alan Miller to include an associated integer array which gives
  ! the positions of the elements in the original order.
  
  IMPLICIT NONE
  REAL, DIMENSION (:), INTENT(IN OUT)  :: list
  INTEGER, DIMENSION (:), INTENT(OUT)  :: order

  ! Local variable
  INTEGER :: i

  DO i = 1, SIZE(list)
     order(i) = i
  END DO
  !write(*,*)"HE", SIZE(order), SIZE(list)
  CALL quick_sort_1(1, SIZE(list))

CONTAINS

  RECURSIVE SUBROUTINE quick_sort_1(left_end, right_end)

    INTEGER, INTENT(IN) :: left_end, right_end

    !     Local variables
    INTEGER             :: i, j, itemp
    REAL   :: reference, temp
    INTEGER, PARAMETER  :: max_simple_sort_size = 6
    
    IF (right_end < left_end + max_simple_sort_size) THEN
       ! Use interchange sort for small lists
       CALL interchange_sort(left_end, right_end)

    ELSE
       ! Use partition ("quick") sort
       reference = list((left_end + right_end)/2)
       i = left_end - 1; j = right_end + 1
       
       DO
    ! Scan list from left end until element >= reference is found
          DO
             i = i + 1
             IF (list(i) >= reference) EXIT
          END DO
          ! Scan list from right end until element <= reference is found
          DO
             j = j - 1
             IF (list(j) <= reference) EXIT
          END DO
          
          
          IF (i < j) THEN
             ! Swap two out-of-order elements
             temp = list(i); list(i) = list(j); list(j) = temp
             itemp = order(i); order(i) = order(j); order(j) = itemp
          ELSE IF (i == j) THEN
             i = i + 1
             EXIT
          ELSE
             EXIT
          END IF
       END DO
       
       IF (left_end < j) CALL quick_sort_1(left_end, j)
       IF (i < right_end) CALL quick_sort_1(i, right_end)
    END IF
    
  END SUBROUTINE quick_sort_1


  SUBROUTINE interchange_sort(left_end, right_end)
    
    INTEGER, INTENT(IN) :: left_end, right_end
    
    !     Local variables
    INTEGER             :: i, j, itemp
    REAL                :: temp
    
    DO i = left_end, right_end - 1
       DO j = i+1, right_end
          IF (list(i) > list(j)) THEN
             temp = list(i); list(i) = list(j); list(j) = temp
             itemp = order(i); order(i) = order(j); order(j) = itemp
          END IF
       END DO
    END DO
    
  END SUBROUTINE interchange_sort
  
END SUBROUTINE quick_sort






end module recon
