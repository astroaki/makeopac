!=======================================================================
!            MAKING OPACITY TABLE USING BOHREN-HUFFMAN PROGRAM
!                        ADAPTED BY B.T. DRAINE
!           MADE INTO F90 AND INTO CURRENT FORM BY C.P. DULLEMOND
!                  ADDED USER-DEFINED SCATTERING ANGLES
!                        BY CPD AND AKIMASA KATAOKA
!
!***********************************************************************
! COMMENTS FROM ORIGINAL CODE:
! Program to interactively call Bohren-Huffman Mie theory program
!
! CALLBHMIE will interactively prompt for:
! 1. refractive index of surrounding medium
! 2. either refractive index or dielectric constant of sphere
! 3. radius of sphere
! 4. wavelength (in vacuo)
! 5. number of angles at which to calculate scattering intensities
!
! CALLBHMIE will return:
! 1. Q_ext, Q_abs, Q_sca, g, Q_back
! 2. If NANG>0, then will also return scattering matrix elements
!    S_11, S_33, S_34, and POL
!
! Adapted by B.T.Draine, Princeton Univ. Obs.
!***********************************************************************
!=======================================================================
program bhmakeopac
  implicit none
  integer, parameter :: MXNANG=1000
  integer :: IREADEP,J,NAN,NANG,JJ
  doubleprecision :: DANG,PI,sum,error,errmax
  real :: QABS,QBACK,QEXT,QSCA,RAD,REFMED,GSCA,POL
  real :: S11,S12,S33,S34,WAVEL,X
  complex :: REFREL,CXEPS,S1(2*MXNANG-1),S2(2*MXNANG-1)
  doubleprecision :: anglerad(MXNANG)
  doubleprecision, allocatable :: optcnst_lambda_cm(:),optcnst_n(:),optcnst_k(:)
  doubleprecision, allocatable :: kappa_lambda_cm(:)
  doubleprecision, allocatable :: kappa_abs(:),kappa_sca(:),kappa_g(:)
  doubleprecision, allocatable :: zscat(:,:,:),angle(:),mu(:),scalefact(:)
  doubleprecision, allocatable :: agrain_cm(:),dagr(:),weight(:),mgrain(:)
  integer :: nlamoc,ilamoc,ilam00,iformat
  integer :: nlamkap,ilamkap,nang180,nagr,ia,leng
  doubleprecision :: agrain_cm_min,agrain_cm_max,xigrain,dum(1:3)
  doubleprecision :: agrain_cm_0,agrain_cm_width,agrain_cm_da
  doubleprecision :: siggeom,factor,plindex,epsl
  character*160 :: filename,material,str0,str1,dustmodel
  logical :: notfinished,fex,found
  PI=4.D0*ATAN(1.D0)
  !
  ! Defaults
  !
  REFMED = 1.d0
  errmax = 0.01
  !
  ! Open and read parameter file
  !
  open(unit=1,file='dust_param.inp')
  read(1,*) iformat
  read(1,*) dustmodel
  dustmodel = trim(dustmodel)
  select case(dustmodel)
  case('powerlawdistr')
     read(1,*) material
     read(1,*) nagr
     read(1,*) agrain_cm_min
     read(1,*) agrain_cm_max
     read(1,*) plindex           ! Should be -3.5 for standard MRN
     read(1,*) xigrain
  case('gaussdistr')
     read(1,*) material
     read(1,*) nagr
     read(1,*) agrain_cm_0       ! Center of the Gaussian
     read(1,*) agrain_cm_width   ! Width of the Gaussian
     read(1,*) agrain_cm_da      ! Half-width of the size grid from
     read(1,*) xigrain
  case default
     write(*,*) 'Error: Dust model ',dustmodel,' not known.'
     close(1)
     stop
  end select
  read(1,*) nang180           ! Nr of angles between 0 and 180 degrees
  read(1,*,end=209) errmax
209 continue
  close(1)
  !
  ! Do a check
  !
  inquire(file='scatangle.inp',exist=fex)
  if(nang180.eq.0) then
     !
     ! This means that we want to read in an angular grid
     !
     if(.not.fex) then
        write(*,*) 'ERROR: In the dust_param.inp file nang180 is set to 0,'
        write(*,*) 'meaning there MUST be a file scatangle.inp instead.'
        write(*,*) 'It does not seem to be there...'
        stop
     endif
     !
     ! Read the scattering angle theta grid
     ! NOTE: Only the part from 0 to 90 degrees. The other part (180 to 90)
     !       is a mirror copy.
     !
     open(unit=1,file='scatangle.inp')
     read(1,*) NANG    ! If you have one bin per angle (e.g. 0, 1, ..., 90) this should be 91
     NANG180 = 2*NANG-1
     IF(nang180.GT.2*MXNANG-1)STOP '***Error: NANG > MXNANG'
     NAN=NANG180
     allocate(mu(nan),angle(nan))
     do j=1,nang
        read(1,*) angle(j)
        mu(j)     = cos(angle(j)*PI/180.)
     enddo
     close(1)
     do j=1,nang-1
        JJ=2*NANG-J
        angle(jj) = 180.d0-angle(j)
        mu(jj)    = cos(angle(jj)*PI/180.)
     enddo
     if(angle(1).ne.0.d0) then
        write(*,*) 'Error: in scatangle.inp the first angle MUST be 0'
        stop
     endif
     if(angle(nang).ne.90.d0) then
        write(*,*) 'Error: in scatangle.inp the last angle MUST be 90'
        stop
     endif
     mu(1) = 1.0d0
     mu(nan) = -1.0d0
  else
     !
     ! Use a simple linear grid in scattering angle theta
     !
     if(fex) then
        write(*,*) 'Warning: There is a file scatangle.inp, but in '
        write(*,*) 'dust_param.inp you specify nang180 (not zero),'
        write(*,*) 'meaning the scatangle.inp file is ignored.'
        write(*,*) 'If you want to use the scatangle.inp file to '
        write(*,*) 'specify the scattering angle grid, please set '
        write(*,*) 'nang180 (line 7 of dust_param.inp) to 0.'
     endif
     if(nang180.lt.60) then
        write(*,*) 'Warning: You have less than 60 angles between '
        write(*,*) '0 and 180 deg.'
        write(*,*) 'Are you sure that this is what you want?'
     endif
     !
     ! Set up the linear spaced grid in scattering angle theta
     !
     ! NANG=number of angles between 0 and 90 degrees (incl. 0 and 90)
     ! Scattering matrix elements are calculated for 2*NANG-1 angles
     ! including 0, 90, and 180 degrees.
     !
     IF(nang180.GT.2*MXNANG-1)STOP '***Error: NANG > MXNANG'
     NANG=(NANG180+1)/2
     DANG=0.5d0*PI/(NANG-1.d0)
     NAN=NANG180
     allocate(mu(nan),angle(nan))
     do j=1,nan
        angle(j)=DANG*((1.d0*J)-1.E0)*180.E0/PI
        mu(j)=cos(angle(j)*PI/180.)
     enddo
     mu(1) = 1.0d0
     mu(nan) = -1.0d0
  endif
  !
  ! For BHMIE_ANGGRID later, create the angle grid in radian.
  !
  anglerad(1:nang) = angle(1:nang)*PI/180.
  !
  ! Now handle the different grain size distribution models
  !
  select case(dustmodel)
  case('powerlawdistr')
     !
     ! --------------------------------------------
     ! DUST MODEL: POWERLAW GRAIN SIZE DISTRIBUTION
     ! --------------------------------------------
     !
     ! Make the grain size grid
     !
     if(nagr.lt.4) then
        write(*,*) 'ERROR: Must have at least 4 grain size sampling points.'
        stop
     endif
     allocate(agrain_cm(nagr),dagr(nagr),weight(nagr),mgrain(nagr))
     do ia=1,nagr
        agrain_cm(ia) = agrain_cm_min * (agrain_cm_max/agrain_cm_min)**((ia-1.d0)/(nagr-1.d0))
     enddo
     dagr(1)    = 0.5d0*(agrain_cm(2)-agrain_cm(1))
     dagr(nagr) = 0.5d0*(agrain_cm(nagr)-agrain_cm(nagr-1))
     do ia=2,nagr-1
        dagr(ia) = 0.5d0*(agrain_cm(ia+1)-agrain_cm(ia-1))
     enddo
     !
     ! Compute mass of grains
     !
     do ia=1,nagr
        mgrain(ia) = (4.d0*pi/3.d0)*xigrain*agrain_cm(ia)**3
     enddo
     !
     ! Make the grain size distribution according to the
     !
     !   N(a)da ~ a^plindex
     !
     ! recipe. Note: An MRN distribution has plindex = -3.5
     !
     sum = 0.d0
     do ia=1,nagr
        weight(ia) = agrain_cm(ia)**plindex * mgrain(ia) * dagr(ia)
        sum = sum + weight(ia)
     enddo
     weight(:) = weight(:) / sum
  case('gaussdistr')
     !
     ! --------------------------------------------
     ! DUST MODEL: GAUSSIAN GRAIN SIZE DISTRIBUTION
     ! --------------------------------------------
     !
     ! Make the grain size grid
     !
     if(nagr.lt.4) then
        write(*,*) 'ERROR: Must have at least 4 grain size sampling points.'
        stop
     endif
     if(agrain_cm_da.lt.0.5*agrain_cm_width) then
        write(*,*) 'ERROR: The grain size window must be at least as large as the half width of the Gauss'
        stop
     endif
     allocate(agrain_cm(nagr),dagr(nagr),weight(nagr),mgrain(nagr))
     agrain_cm_min = agrain_cm_0 - agrain_cm_da
     agrain_cm_max = agrain_cm_0 + agrain_cm_da
     do ia=1,nagr
        agrain_cm(ia) = agrain_cm_min * (agrain_cm_max/agrain_cm_min)**((ia-1.d0)/(nagr-1.d0))
     enddo
     dagr(1)    = 0.5d0*(agrain_cm(2)-agrain_cm(1))
     dagr(nagr) = 0.5d0*(agrain_cm(nagr)-agrain_cm(nagr-1))
     do ia=2,nagr-1
        dagr(ia) = 0.5d0*(agrain_cm(ia+1)-agrain_cm(ia-1))
     enddo
     !
     ! Compute mass of grains
     !
     do ia=1,nagr
        mgrain(ia) = (4.d0*pi/3.d0)*xigrain*agrain_cm(ia)**3
     enddo
     !
     ! Make the grain size distribution according to the
     !
     !   N(a)da ~ exp(-0.5*(a-a0)^2/width^2) da
     !
     ! recipe. The width is mainly to get rid of the strong resonances.
     !
     !
     sum = 0.d0
     do ia=1,nagr
        weight(ia) = exp(-0.5d0*(agrain_cm(ia)-agrain_cm_0)**2/agrain_cm_width**2) * dagr(ia)
        sum = sum + weight(ia)
     enddo
     weight(:) = weight(:) / sum
  case default
     write(*,*) 'Error: Dust model ',dustmodel,' not known.'
     stop
  end select
  !
  ! Open optical constants file
  !
  filename = trim(material)//".lnk"
  nlamoc = 0  ! nlamoc stands for nr of lambda points in the optical constants file
  notfinished = .true.
  open(unit=1,file=filename)
  do while(notfinished)
     read(1,*,end=20) dum
     nlamoc = nlamoc + 1
  enddo
20 continue
  close(1)
  !
  ! Read the optical constants file
  !
  allocate(optcnst_lambda_cm(nlamoc),optcnst_n(nlamoc),optcnst_k(nlamoc))
  open(unit=1,file=filename)
  do ilamoc=1,nlamoc
     read(1,*) optcnst_lambda_cm(ilamoc),optcnst_n(ilamoc),optcnst_k(ilamoc)
  enddo
  close(1)
  optcnst_lambda_cm = optcnst_lambda_cm * 1d-4
  !
  ! Check if the file wavelength_micron_bh.inp is available
  ! If so, we will use those wavelengths for the output, otherwise
  ! we will use the same wavelengths sample as in the .lnk file.
  !
  inquire(file='wavelength_micron_bh.inp',exist=fex)
  if(fex) then
     !
     ! Read the wavelength_micron_bh.inp file: this will now be the
     ! wavelength grid used for the output
     !
     open(unit=1,file='wavelength_micron_bh.inp')
     read(1,*) nlamkap
     allocate(kappa_lambda_cm(nlamkap))
     do ilamkap=1,nlamkap
        read(1,*) kappa_lambda_cm(ilamkap)
     enddo
     close(1)
     kappa_lambda_cm = kappa_lambda_cm * 1d-4
  else
     !
     ! Use the same wavelengths sample as in the .lnk file.
     !
     nlamkap = nlamoc
     allocate(kappa_lambda_cm(nlamkap))
     kappa_lambda_cm(:) = optcnst_lambda_cm(:)
  endif
  !
  ! Allocate the rest of the output arrays
  !
  allocate(kappa_abs(nlamkap),kappa_sca(nlamkap),kappa_g(nlamkap),scalefact(nlamkap))
  allocate(zscat(6,2*MXNANG-1,nlamkap))
  !
  ! Reset things
  !
  zscat(:,:,:) = 0.d0
  kappa_abs(:) = 0.d0
  kappa_sca(:) = 0.d0
  kappa_g(:)   = 0.d0
  !
  ! Start the loop over grain sizes
  !
  do ia=1,nagr
     !
     ! Compute geometric cross section
     !
     siggeom = pi*agrain_cm(ia)**2
     !
     ! Now do the loop over wavelengths
     !
     do ilamkap=1,nlamkap
        !
        ! Prepare the parameters for BHMie
        !
        ! The complex index of refraction
        !
        found = .false.
        do ilamoc=nlamoc-1,1,-1
           if((kappa_lambda_cm(ilamkap)-optcnst_lambda_cm(ilamoc))*   &
              (kappa_lambda_cm(ilamkap)-optcnst_lambda_cm(ilamoc+1)).le.0.d0) then
              ilam00 = ilamoc
              found = .true.
           endif
        enddo
        if(.not.found) then
           write(*,*) 'Error: Wavelength was not found in the range given by the optical constants file.'
           stop
        endif
        ilamoc = ilam00
        epsl   = (kappa_lambda_cm(ilamkap)-optcnst_lambda_cm(ilamoc))/   &
                 (optcnst_lambda_cm(ilamoc+1)-optcnst_lambda_cm(ilamoc))
        refrel = cmplx((1.d0-epsl)*optcnst_n(ilamoc)+epsl*optcnst_n(ilamoc+1), &
                       (1.d0-epsl)*optcnst_k(ilamoc)+epsl*optcnst_k(ilamoc+1))/refmed
        !
        ! Radius of the grain in cm
        !
        rad = agrain_cm(ia)
        !
        ! Wavelength in cm
        !
        wavel = kappa_lambda_cm(ilamkap)
        !
        ! Compute the dimensionless grain size size
        !
        X=2.E0*PI*RAD*REFMED/WAVEL
        !
        ! Call BHMie
        !
        CALL BHMIE_ANGGRID(X,REFREL,NANG,ANGLERAD,S1,S2,QEXT,QSCA,QBACK,GSCA)
        QABS=QEXT-QSCA
        !
        ! Put results into array
        !
        ! Note: The averaging of g has to be done multiplied by kappa_scat,
        !       otherwise the answer is wrong.
        !
        kappa_abs(ilamkap) = kappa_abs(ilamkap) + weight(ia)*qabs*siggeom/mgrain(ia)
        kappa_sca(ilamkap) = kappa_sca(ilamkap) + weight(ia)*qsca*siggeom/mgrain(ia)
        kappa_g(ilamkap)   = kappa_g(ilamkap)   + weight(ia)*GSCA*qsca*siggeom/mgrain(ia)
        !
        ! Compute conversion factor from the Sxx matrix elements
        ! from the Bohren & Huffman code to the Zxx matrix elements we
        ! use (such that 2*pi*int_{-1}^{+1}Z11(mu)dmu=kappa_scat).
        ! This includes the factor k^2 (wavenumber squared) to get
        ! the actual cross section in units of cm^2 / ster, and there
        ! is the mass of the grain to get the cross section per gram.
        !
        factor = (kappa_lambda_cm(ilamkap)/(2*PI))**2/mgrain(ia)
        !
        ! Also store the Z matrix elements
        !
        NAN=2*NANG-1
        DO J=1,NAN
           S11=0.5E0*CABS(S2(J))*CABS(S2(J))
           S11=S11+0.5E0*CABS(S1(J))*CABS(S1(J))
           S12=0.5E0*CABS(S2(J))*CABS(S2(J))
           S12=S12-0.5E0*CABS(S1(J))*CABS(S1(J))
           POL=-S12/S11
           S33=REAL(S2(J)*CONJG(S1(J)))
           S34=AIMAG(S2(J)*CONJG(S1(J)))
           zscat(1,j,ilamkap) = zscat(1,j,ilamkap) + weight(ia) * S11 * factor
           zscat(2,j,ilamkap) = zscat(2,j,ilamkap) + weight(ia) * S12 * factor
           zscat(3,j,ilamkap) = zscat(3,j,ilamkap) + weight(ia) * S11 * factor
           zscat(4,j,ilamkap) = zscat(4,j,ilamkap) + weight(ia) * S33 * factor
           zscat(5,j,ilamkap) = zscat(5,j,ilamkap) + weight(ia) * S34 * factor
           zscat(6,j,ilamkap) = zscat(6,j,ilamkap) + weight(ia) * S33 * factor
        enddo
        !
        ! End loop over all wavelengths
        !
     enddo
     !
     ! End loop over all grain sizes
     !
  enddo
  !
  ! Now renormalize the g factor
  !
  kappa_g(:) = kappa_g(:) / kappa_sca(:)
  !
  ! Check if the sum of the S11 over all angles is indeed kappa_scat
  !
  do ilamkap=1,nlamkap
     sum = 0.d0
     do j=2,nan
        sum = sum + 0.25d0*(zscat(1,j-1,ilamkap)+zscat(1,j,ilamkap))* &
              abs(mu(j)-mu(j-1))
     enddo
     sum = sum * 4*PI
     error = abs(sum/kappa_sca(ilamkap)-1.d0)
     if(error.gt.errmax) then
        write(*,*) 'ERROR: At lambda=',kappa_lambda_cm(ilamkap)*1d4,'micron the error in the',&
             ' scattering integral is ',error,'which is larger than the error limit',errmax
        write(*,*) '    kappa_scat                     = ',kappa_sca(ilamkap)
        write(*,*) '    2*pi*int_{-1}^{+1} Z_11(mu)dmu = ',sum
        write(*,*) '    Please use a larger number of angle points or take a weaker error limit (5th line in param.inp).'
        write(*,*) '    For now: we do not rescale and just hope for the best...'
        close(2)
        scalefact(ilamkap) = 1.0
     else
        scalefact(ilamkap) = kappa_sca(ilamkap) / sum
        zscat(1:6,1:nan,ilamkap) = zscat(1:6,1:nan,ilamkap) * scalefact(ilamkap)
     endif
  enddo
  !
  ! Write the results
  !
  filename = 'dustkapscatmat_'//trim(material)//'.inp'
  open(unit=2,file=filename)
  leng = len_trim(material)
  if(leng.lt.10) then
     write(str0,'(I1)')
  elseif(leng.lt.100) then
     write(str0,'(I2)')
  elseif(leng.lt.1000) then
     write(str0,'(I3)')
  else
     write(*,*) 'Dust opacity name too long'
     stop
  endif
  str1 = '(A41,A'//trim(str0)//')'
  write(2,str1) '# Opacity and scattering matrix file for ',trim(material)
  write(2,'(A109)') '# Please do not forget to cite in your publications the original ' &
       //'paper of these optical constant measurements'
  write(2,'(A44)') '# Made with the make_scatmat_distr.f90 code,'
  write(2,'(A70)') '# using the bhmie.f Mie code of Bohren and Huffman (version by Draine)'
  write(2,'(A26)') '# Grain size distribution:'
  write(2,'(A23,E13.6,A3)') '#   agrain_min       = ',agrain_cm_min,' cm'
  write(2,'(A23,E13.6,A3)') '#   agrain_max       = ',agrain_cm_max,' cm'
  write(2,'(A23,F13.6)') '#   powerlaw index   = ',plindex
  write(2,'(A23,I4,A3)') '#   nr of sizes used = ',nagr
  write(2,'(A20,F9.6,A7)') '# Material density =',xigrain,' g/cm^3'
  write(2,*) 1     ! Format number
  write(2,*) nlamkap
  write(2,*) nan
  write(2,*)
  do ilamkap=1,nlamkap
     write(2,'(4(E13.6,1X))') kappa_lambda_cm(ilamkap)*1e4,kappa_abs(ilamkap), &
                              kappa_sca(ilamkap),kappa_g(ilamkap)
  enddo
  write(2,*)
  do j=1,nan
     write(2,'(F13.6)') angle(j)
  enddo
  write(2,*)
  do ilamkap=1,nlamkap
     do j=1,nan
        write(2,'(6(E13.6,1X))') zscat(1:6,j,ilamkap)
     enddo
     write(2,*)
  enddo
  close(2)
  !
  ! Just for information to the user: write out the scaling factor used
  !
  open(unit=1,file='scalefactor.out')
  write(1,*) nlamkap
  do ilamkap=1,nlamkap
     write(1,'(2(E13.6,1X))') kappa_lambda_cm(ilamkap)*1e4,scalefact(ilamkap)
  enddo
  close(1)
  !
  ! Deallocate stuff
  !
  deallocate(optcnst_lambda_cm,optcnst_n,optcnst_k)
  deallocate(kappa_abs,kappa_sca,kappa_g)
  deallocate(zscat,angle,mu,scalefact)
  deallocate(agrain_cm,dagr,weight,mgrain)
end program bhmakeopac
