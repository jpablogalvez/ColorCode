!==========================================================
!
       module colour
!
       implicit none
!
       include 'colour.h'
!
       contains
!
!======================================================================!
!
! References
! ----------
!
! - Beck, M. E. (2005), Estimation of physiologically perceived color
!    from TDDFT-derived excitation spectra. Int. J. Quantum Chem.,
!    101: 683–689. https://doi.org/10.1002/qua.20326
! - Beyond λ[lambda]max: Transforming Visible Spectra into 24-Bit Color
!    Values Darren L. Williams, Thomas J. Flaherty, Casie L. Jupe,
!    Stephanie A. Coleman, Kara A. Marquez, and Jamie J. Stanton Journal
!    of Chemical Education 2007 84 (11), 1873.
!    https://doi.org/10.1021/ed084p1873
! - J. A. Stephen Viggiano, Nanette Salvaggio, Nitin Sampat, 
!    "Chromaticity Matrix to Tristimulus Matrix Conversion for RGB Color
!    Spaces–Even In the Dark"  in Proc. IS&T Int’l. Symp. on Electronic 
!    Imaging: Color Imaging xpIII: Displaying, Processing, Hardcopy, and
!    Applications,  2018,  pp 324-1 - 324-7. 
!    https://doi.org/10.2352/ISSN.2470-1173.2018.16.COLOR-324
! - SMPTE RP 177-1993, Derivation of Basic Television Color Equations
! - ASTM E308-01. Standard Practice for Computing the Colors of  
!    Objects by Using the CIE System; ASTM International: West 
!    Conshohoken, PA, 2001
! - http://www.brucelindbloom.com
!
       subroutine colorate(illu,iadapt,itype,nref,xref,yref,           &
                           Xw,Yw,Zw,X,Y,Z,R,G,B)
!
       implicit none
!
! Input/Output variables
!
       character(len=10),intent(in)             ::  illu      !
       real(kind=4),dimension(nref),intent(in)  ::  xref      ! 
       real(kind=4),dimension(nref),intent(in)  ::  yref      !
       real(kind=4),intent(out)                 ::  X,Y,Z     !  XYZ tristimulus values
       real(kind=4),intent(inout)               ::  Xw,Yw,Zw  !  Illuminant white point XYZ coordinates
       real(kind=4),intent(out)                 ::  R,G,B     !
       integer,intent(in)                       ::  nref      !
       integer,intent(in)                       ::  iadapt    !
       integer,intent(in)                       ::  itype     !
!
! Local variables
! 
       character(len=10)                        ::  refillu   !  Reference illuminant
       real(kind=4)                             ::  xr,yr,zr  !  RGB system red point chromatic coordinates
       real(kind=4)                             ::  xg,yg,zg  !  RGB system green point chromatic coordinates
       real(kind=4)                             ::  xb,yb,zb  !  RGB system blue point chromatic coordinates
       real(kind=4)                             ::  rx,ry,rz  !
       real(kind=4)                             ::  gx,gy,gz  !
       real(kind=4)                             ::  bx,by,bz  !
       real(kind=4)                             ::  sr,sg,sb  !  Tristimulus sums
       real(kind=4)                             ::  Xwref     !  Reference illuminant white
       real(kind=4)                             ::  Ywref     !  Reference illuminant white
       real(kind=4)                             ::  Zwref     !  Reference illuminant white
       real(kind=4)                             ::  gam       !  Gamma value
       real(kind=4)                             ::  saux      !
!
! Converting a spectrum to a colour
! ---------------------------------
!
! Computing XYZ tristimulus values from spectral data
! ...................................................
!
       if ( itype .eq. 1 ) then
!
         call xyzexp(nref,xref,yref,X,Y,Z)
!
       else
         select case (trim(illu))
!
           case('inverse')
!
             call xyzinverse(X,Y,Z)
!
             Xw = 1.0d0
             Yw = 1.0d0
             Zw = 1.0d0
!
           case('a','c','d50','d55','d65','d75','e')
!
             call xyzstandard(X,Y,Z,Xw,Yw,Zw)
!
           case('black-body')
!
             call xyzbb(X,Y,Z,Xw,Yw,Zw)
!
         end select
       end if
!
! Setting up RGB system information
! .................................
!
       call rgbmodel(refillu,gam,Xwref,Ywref,Zwref,xr,yr,zr,           &
                     xg,yg,zg,xb,yb,zb)
!
! Chromatic adaptation
! ....................
!
       if ( (trim(illu).ne.trim(refillu)) .and. (iadapt.ne.4) )        &
                call adaptation(iadapt,X,Y,Z,Xw,Yw,Zw,Xwref,Ywref,Zwref)
!
! Building the transformation matrix
! ..................................
!
       rx = yg*zb - yb*zg
       ry = xb*zg - xg*zb
       rz = xg*yb - xb*yg
!
       gx = yb*zr - yr*zb
       gy = xr*zb - xb*zr
       gz = xb*yr - xr*yb
!
       bx = yr*zg - yg*zr
       by = xg*zr - xr*zg
       bz = xr*yg - xg*yr
!
!
!
       sr = rx*Xwref + ry*Ywref + rz*Zwref
       sg = gx*Xwref + gy*Ywref + gz*Zwref
       sb = bx*Xwref + by*Ywref + bz*Zwref
!
!
!
       rx = rx/sr
       ry = ry/sr
       rz = rz/sr
!
       gx = gx/sg
       gy = gy/sg
       gz = gz/sg
!
       bx = bx/sb
       by = by/sb
       bz = bz/sb       
!
! Conversion from XYZ to RGB
! ..........................
!
! XYZ tristimulus values must be in the nominal range [0, 1]
! 
       saux = max(X,Y,Z)
!
       X = X/saux
       Y = Y/saux
       Z = Z/saux
!
       R = rx*X + ry*Y + rz*Z
       G = gx*X + gy*Y + gz*Z  
       B = bx*X + by*Y + bz*Z
!
! Approximating out-of-gamut colours
! ----------------------------------
!

!
! Applying Gamma correction (companding)
! .........................
!
       if ( trim(illu) .eq. 'srgb' ) then
!
         if ( R .le. 0.0031308 ) then
           R = 12.92*R
         else
           R = 1.055*R**(1.0/2.4) - 0.055
         end if
!
         if ( G .le. 0.0031308 ) then
           G = 12.92*G
         else
           G = 1.055*G**(1.0/2.4) - 0.055
         end if
!
         if ( B .le. 0.0031308 ) then
           B = 12.92*B
         else
           B = 1.055*B**(1.0/2.4) - 0.055
         end if
!
         saux = max(R,G,B)
!
         R = R/saux
         G = G/saux
         B = B/saux
!
       else if ( trim(illu) .eq. 'eci' ) then
! 
         if ( R .le. 0.008856 ) then
           R = R*0.9033
         else
           R = 1.16*R**(1.0/3.0) - 0.16
         end if
!
         if ( G .le. 0.008856 ) then
           G = G*0.9033
         else
           G = 1.16*G**(1.0/3.0) - 0.16
         end if
!
         if ( B .le. 0.008856 ) then
           B = B*0.9033
         else
           B = 1.16*B**(1.0/3.0) - 0.16
         end if
!
         saux = max(R,G,B)
!
         R = R/saux
         G = G/saux
         B = B/saux
!
       else
!
         saux = max(R,G,B)
!
         R = (R/saux)**(1.0d0/gam)
         G = (G/saux)**(1.0d0/gam)
         B = (B/saux)**(1.0d0/gam)
!
       end if
!
       return
       end subroutine colorate
!
!======================================================================!
!
       subroutine xyzexp(n,xdat,ydat,X1,Y1,Z1)
!
       use variables
!
       implicit none
!
! Input/Output variables
!
       real(kind=4),intent(out)              ::  X1,Y1,Z1  !  XYZ tristimulus values
       real(kind=4),dimension(n),intent(in)  ::  xdat      ! 
       real(kind=4),dimension(n),intent(in)  ::  ydat      !
       integer,intent(in)                    ::  n         !

!
! Local variables
!
       real(kind=4),dimension(n)             ::  xspec     ! 
       real(kind=4),dimension(n)             ::  yspec     !
!
! Computing XYZ tristimulus values from the experimental spectrum
! ---------------------------------------------------------------
!
       xspec(:) = xdat(:)
       yspec(:) = ydat(:)
!
       if ( ixunit .ne. 2 ) then   !  Spectrum must be in wavelength units
         xspec(:) = 1.0E7/xspec(:)
       end if
!
       if ( iyunit .ne. 3 ) then   !  BB method needs transmittance spectrum

       end if
!
       if ( iyunit .ne. 2 ) then   !  Inverse method needs absorbance spectrum

       end if
!

!
       X = simpuneven(n,xspec,yspec,Xcie)
       Y = simpuneven(n,xspec,yspec,Ycie)
       Z = simpuneven(n,xspec,yspec,Zcie)
!
       return
       end subroutine xyzexp
!
!======================================================================!
!
! References
! ----------
!
! - Beck, M. E. (2005), Estimation of physiologically perceived color
!    from TDDFT-derived excitation spectra. Int. J. Quantum Chem.,
!    101: 683–689. doi: 10.1002/qua.20326
!
       subroutine xyzinverse(X1,Y1,Z1)
!
       use variables
!
       implicit none
!
! Input/Output variables
!
       real(kind=4),intent(out)               ::  X1,Y1,Z1  !  XYZ tristimulus values
!
! Local variables
!
       real(kind=4),dimension(:),allocatable  ::  xspec    !
       real(kind=4),dimension(:),allocatable  ::  yspec    !
       real(kind=4)                           ::  xa,xb    !  Integration limits
       real(kind=4)                           ::  thr      !  Threshold value
       real(kind=4)                           ::  dif      !
       real(kind=4)                           ::  h        !  Step length
       real(kind=4)                           ::  smax     !
       real(kind=4)                           ::  smin     !
       real(kind=4)                           ::  X0       !
       real(kind=4)                           ::  Y0       !
       real(kind=4)                           ::  Z0       !
       integer                                ::  n        !
!
! Computing XYZ tristimulus values through the inverse spectrum method
! --------------------------------------------------------------------
!
       thr = 1.0E-8
!
       xa = 300.0d0
       xb = 800.0d0
!
       X1 = 0.0
       Y1 = 0.0
       Z1 = 0.0
!
       dif = 100.0
!
! Setting the initial number of subintervals
!
       h = dxnm
!
       n = (xb - xa)/h
       if ( mod(n,2) .ne. 0 ) n = n + 1
!
! Integrating iteratively until the RMSD of the XYZ values between two
!  successive iterations is lower than a predefined threshold
!
       do while ( dif .gt. thr )
!
         X0 = X1
         Y0 = Y1
         Z0 = Z1
!
! Precomputing an approximation of the spectrum's shape in the energy
!  domain folding the oscillator strengths with broadening functions
!
         h = (xb - xa)/n
!
         allocate(xspec(n+1),yspec(n+1))
!
         call precompinv(n+1,xspec,yspec,xa,h)
!
! Estimating the inverse spectrum by taking the mirror image of the
!  spectra
!
         smax = maxval(yspec)
         smin = minval(yspec)
!
         yspec(:) = smax - smin  - yspec(:)
!
! Integrating inverse spectrum (calculation of XYZ tristimulus values)
!
         X1 = simpeven(n+1,xspec,yspec,h,Xcie)
         Y1 = simpeven(n+1,xspec,yspec,h,Ycie)
         Z1 = simpeven(n+1,xspec,yspec,h,Zcie)
!
! Increasing the number of subintervals for the next iteration
!
         n = 2*n
!
         dif = sqrt((X1-X0**2)+(Y1-Y0)**2+(Z1-Z0)**2)
!
         deallocate(xspec,yspec)
!
       end do
!
       return
       end subroutine xyzinverse
!
!======================================================================!
!
! References
! ----------
!
! - Beyond λ[lambda]max: Transforming Visible Spectra into 24-Bit Color
!    Values Darren L. Williams, Thomas J. Flaherty, Casie L. Jupe,
!    Stephanie A. Coleman, Kara A. Marquez, and Jamie J. Stanton Journal
!    of Chemical Education 2007 84 (11), 1873.
!    https://doi.org/10.1021/ed084p1873
!
       subroutine xyzstandard(X1,Y1,Z1,Xw,Yw,Zw)
!
       use variables
       use illuminants
!
       implicit none
!
! Input/Output variables
!
       real(kind=4),intent(out)    ::  X1,Y1,Z1  !  XYZ tristimulus values
       real(kind=4),intent(out)    ::  Xw,Yw,Zw  !  White point coordinates
!
! Local variables
!
       real(kind=4),dimension(81)  ::  xspec     !
       real(kind=4),dimension(81)  ::  yspec     !
       real(kind=4),dimension(81)  ::  yillu     !
       real(kind=4),dimension(81)  ::  yprod     !
       real(kind=4)                ::  n         !
       real(kind=4)                ::  xa,xb     !  Integration limits
       real(kind=4)                ::  h         !  Step length
!
! Computing XYZ tristimulus values using an standard illuminant
! -------------------------------------------------------------
!
       xa = 380.0d0
       xb = 780.0d0
!
       h = 5.0d0
!
! Choosing illuminant
!
       select case (trim(illu))
         case('a')
           call illuminantA(xspec,yillu,Xw,Yw,Zw)
         case('c')
           call illuminantC(xspec,yillu,Xw,Yw,Zw)
         case('d50')
           call illuminantD50(xspec,yillu,Xw,Yw,Zw)
         case('d55')
           call illuminantD55(xspec,yillu,Xw,Yw,Zw)
         case('d65')
           call illuminantD65(xspec,yillu,Xw,Yw,Zw)
         case('d75')
           call illuminantD75(xspec,yillu,Xw,Yw,Zw)
         case('e')
           call illuminantE(xspec,yillu,Xw,Yw,Zw)
       end select
!
! Computing molar absorption coefficent at every wavelength
!
       call precompillu(81,xspec,yspec)
!
! Converting to transmittance spectrum
!
       yspec(:) = 10.0d0**(-yspec(:)*conc*path)
!
! Calculating XYZ tristimulus values
!
       N = simpeven(81,xspec,yillu,h,Ycie)
!
       yprod(:) = yillu(:)*yspec(:)
!
       X1 = simpeven(81,xspec,yprod,h,Xcie)/N
       Y1 = simpeven(81,xspec,yprod,h,Ycie)/N
       Z1 = simpeven(81,xspec,yprod,h,Zcie)/N
!
       return
       end subroutine xyzstandard
!
!======================================================================!
!
! References
! ----------
!
! - Beyond λ[lambda]max: Transforming Visible Spectra into 24-Bit Color
!    Values Darren L. Williams, Thomas J. Flaherty, Casie L. Jupe,
!    Stephanie A. Coleman, Kara A. Marquez, and Jamie J. Stanton Journal
!    of Chemical Education 2007 84 (11), 1873.
!    https://doi.org/10.1021/ed084p1873
!
       subroutine xyzbb(X1,Y1,Z1,Xw,Yw,Zw)
!
       use lineshape
!
       implicit none
!
! Input/Output variables
!
       real(kind=4),intent(out)  ::  X1,Y1,Z1   !  XYZ tristimulus values
       real(kind=4),intent(out)  ::  Xw,Yw,Zw   !  White point coordinates
!
! Local variables
!
       real(kind=4)              ::  xa,xb      !  Integration limits
       real(kind=4)              ::  fx         ! 
       real(kind=4)              ::  h          !  Step length
       real(kind=4)              ::  BBxa,BBxb  ! 
       real(kind=4)              ::  BBya,BByb  ! 
       real(kind=4)              ::  BBza,BBzb  ! 
       real(kind=4)              ::  BBx        ! 
       real(kind=4)              ::  fxa,fxb    !  Value of the function at the integration limits
       real(kind=4)              ::  fya,fyb    !  Value of the function at the integration limits
       real(kind=4)              ::  fza,fzb    !  Value of the function at the integration limits
       real(kind=4)              ::  X0         !
       real(kind=4)              ::  Y0         !
       real(kind=4)              ::  Z0         !
       real(kind=4)              ::  Xw0        !
       real(kind=4)              ::  Yw0        !
       real(kind=4)              ::  Zw0        !
       real(kind=4)              ::  thr        !  Threshold value
       real(kind=4)              ::  dif        !
       real(kind=4)              ::  fxodd      !  Auxliary variable
       real(kind=4)              ::  fyodd      !  Auxliary variable
       real(kind=4)              ::  fzodd      !  Auxliary variable
       real(kind=4)              ::  fxeven     !  Auxliary variable
       real(kind=4)              ::  fyeven     !  Auxliary variable
       real(kind=4)              ::  fzeven     !  Auxliary variable
       real(kind=4)              ::  BBxodd     !  Auxliary variable
       real(kind=4)              ::  BByodd     !  Auxliary variable
       real(kind=4)              ::  BBzodd     !  Auxliary variable
       real(kind=4)              ::  BBxeven    !  Auxliary variable
       real(kind=4)              ::  BByeven    !  Auxliary variable
       real(kind=4)              ::  BBzeven    !  Auxliary variable
       real(kind=4)              ::  width      !
       real(kind=4)              ::  sig        !
       real(kind=4)              ::  fg         !
       real(kind=4)              ::  fl         !
       real(kind=4)              ::  f          !
       real(kind=4)              ::  eta        !
       integer                   ::  m          !  Number of points
       integer                   ::  n          !  Number of subintervals
       integer                   ::  i          !
!
! Computing XYZ tristimulus values using a Black Body as illuminant
! -----------------------------------------------------------------
!
       thr = 1.0E-8
!
       h = dxnm
!
       xa = 300.0d0
       xb = 800.0d0
!
       X1 = 0.0
       Y1 = 0.0
       Z1 = 0.0
!
       Xw = 0.0
       Yw = 0.0
       Zw = 0.0
!
       m = (xb - xa)/h + 1
       if ( mod(n,2) .eq. 0 ) m = m + 1
!
! Calculating the value of the funtion at the lower limit of integration
!
       BBx = blackbody(xa)
!
       BBxa = BBx*Xcie(xa)
       BBya = BBx*Ycie(xa)
       BBza = BBx*Zcie(xa)
!
       if ( ixunit .eq. 1 ) then
         width = 1.0E7/hwhm
         fg    = 1.0E7/gauhwhm
         fl    = 1.0E7/lorhwhm
       else
         width = hwhm
         fg    = gauhwhm
         fl    = lorhwhm
       end if
!
       f = (fg**5 + 2.69269*fg**4*fl + 2.42843*fg**3*fl*2              &
            + 4.47163*fg**2*fl**3 + 0.07842*fg*fl**4 + fl**5)**(1.0/5.0)
!
       eta = 1.36603*(fl/f) - 0.47719*(fl/f)**2.0 + 0.11116*(fl/f)**3.0
!
       if ( iline .eq. 1 ) then
!
         sig = width*sqrt(2.0*log(2.0))
!
         fx = 0.0
         do  i = 1, nband
           fx = fx + real(Aabs)*inten(i)*(1.0E7/xa)                    &
                              *gauss((1.0E7/xa-1.0E7/freq(i)),1.0E7/sig)
         end do
!
       else if ( iline .eq. 2 ) then
!
         fx = 0.0
         do  i = 1, nband
           fx = fx + real(Aabs)*inten(i)*(1.0E7/xa)                    &
                              *lor((1.0E7/xa-1.0E7/freq(i)),1.0E7/width)
         end do
!
       else
!
         sig = f*sqrt(2.0*log(2.0))
!
         fx = 0.0
         do  i = 1, nband
           fx = fx + real(Aabs)*inten(i)*(1.0E7/xa)                    &
                  *voigt((1.0E7/xa-1.0E7/freq(i)),1.0E7/f,1.0E7/sig,eta)
         end do
!
       end if
!
       fx = 10.0**(-fx*conc*path)
!
       fxa = BBxa*fx
       fya = BBya*fx
       fza = BBza*fx
!
       n = m - 1
!
! Integrating iteratively until the RMSD of the XYZ values between two
!  successive iterations is lower than a predefined threshold
!
       dif = 100
!
       do while ( dif .gt. thr )
!
         X0 = X1
         Y0 = Y1
         Z0 = Z1
!
         Xw0 = Xw
         Yw0 = Yw
         Zw0 = Zw
!
         h = (xb - xa)/n
!
! Calculating Simpson's rule sums for the selected lineshape
!
         if ( iline .eq. 1 ) then
           call simpitergauss(fxodd,fyodd,fzodd,fxeven,fyeven,         &
                              fzeven,BBxodd,BByodd,BBzodd,BBxeven,     &
                              BByeven,BBzeven,n,xa,h,sig)
         else if ( iline .eq. 2 ) then
           call simpitergauss(fxodd,fyodd,fzodd,fxeven,fyeven,         &
                              fzeven,BBxodd,BByodd,BBzodd,BBxeven,     &
                              BByeven,BBzeven,n,xa,h,width)
         else
           call simpitervoigt(fxodd,fyodd,fzodd,fxeven,fyeven,         &
                              fzeven,BBxodd,BByodd,BBzodd,BBxeven,     &
                              BByeven,BBzeven,n,xa,h,f,sig,eta)
         end if
!
! Calculating the value of the funtion at the upper limit of integration
!
         xb = xa + (m-1)*h
!
         BBx = blackbody(xb)
!
         BBxb = BBx*Xcie(xb)
         BByb = BBx*Ycie(xb)
         BBzb = BBx*Zcie(xb)
!
         if ( iline .eq. 1 ) then
!
           fx = 0.0
           do  i = 1, nband
             fx = fx + real(Aabs)*inten(i)*(1.0E7/xb)                  &
                              *gauss((1.0E7/xb-1.0E7/freq(i)),1.0E7/sig)
           end do
!
         else if ( iline .eq. 2 ) then
!
           fx = 0.0
           do  i = 1, nband
             fx = fx + real(Aabs)*inten(i)*(1.0E7/xb)                  &
                              *lor((1.0E7/xb-1.0E7/freq(i)),1.0E7/width)
           end do
!
         else
!
           fx = 0.0
           do  i = 1, nband
             fx = fx + real(Aabs)*inten(i)*(1.0E7/xb)                  &
                  *voigt((1.0E7/xb-1.0E7/freq(i)),1.0E7/f,1.0E7/sig,eta)
           end do
!
         end if
!
         fx = 10.0**(-fx*conc*path)
!
         fxb = BBxb*fx
         fyb = BByb*fx
         fzb = BBzb*fx
!
! Calculating the approximated area (integral)
!
         X1 = h/3.0d0*(fxa + fxb + 2*fxeven + 4*fxodd)
         Y1 = h/3.0d0*(fya + fyb + 2*fyeven + 4*fyodd)
         Z1 = h/3.0d0*(fza + fzb + 2*fzeven + 4*fzodd)
!
         Xw = h/3.0d0*(BBxa + BBxb + 2*BBxeven + 4*BBxodd)
         Yw = h/3.0d0*(BBya + BByb + 2*BByeven + 4*BByodd)
         Zw = h/3.0d0*(BBza + BBzb + 2*BBzeven + 4*BBzodd)
!
! Increasing the number of subintervals for the next iteration
!
         n = 2*n
!
         dif = sqrt((X1-X0**2) + (Y1-Y0)**2 + (Z1-Z0)**2               &
                    + (Xw-Xw0**2) + (Yw-Yw0)**2 + (Zw-Zw0)**2)
!
       end do
!
       Xw = Xw/Yw
       Zw = Zw/Yw
       Yw = 1.0d0
!
       return
       end subroutine xyzbb
!
!======================================================================!
!
       subroutine simpitergauss(fxodd,fyodd,fzodd,fxeven,fyeven,       &
                                fzeven,BBxodd,BByodd,BBzodd,BBxeven,   &
                                BByeven,BBzeven,n,xa,h,sig)
!
       use lineshape
!
       implicit none
!
! Input/Output variables
!
       real(kind=4),intent(out)  ::  fxodd     !
       real(kind=4),intent(out)  ::  fyodd     !
       real(kind=4),intent(out)  ::  fzodd     !
       real(kind=4),intent(out)  ::  fxeven    !
       real(kind=4),intent(out)  ::  fyeven    !
       real(kind=4),intent(out)  ::  fzeven    !
       real(kind=4),intent(out)  ::  BBxodd    !
       real(kind=4),intent(out)  ::  BByodd    !
       real(kind=4),intent(out)  ::  BBzodd    !
       real(kind=4),intent(out)  ::  BBxeven   !
       real(kind=4),intent(out)  ::  BByeven   !
       real(kind=4),intent(out)  ::  BBzeven   !
       real(kind=4),intent(in)   ::  xa        !
       real(kind=4),intent(in)   ::  h         !
       real(kind=4),intent(in)   ::  sig       !
       integer,intent(in)        ::  n         !
!
! Local variables
!
       real(kind=4)              ::  xp        !
       real(kind=4)              ::  fx        !
       real(kind=4)              ::  BBx       !
       integer                   ::  i,j       !
!
! Computing the Simpson's rule sums for a Gaussian broadening
! -----------------------------------------------------------
!
! Adding the area of the function at each odd-subinterval
!
       fxodd  = 0.0
       fyodd  = 0.0
       fzodd  = 0.0
!
       BBxodd = 0.0
       BByodd = 0.0
       BBzodd = 0.0
!
       xp = xa
!
       do i = 2, n, 2
!
         xp = xp + (i-1)*h
!
         fx = 0.0
         do  j = 1, nband
           fx = fx + real(Aabs)*inten(j)*(1.0E7/xp)                    &
                              *gauss((1.0E7/xp-1.0E7/freq(j)),1.0E7/sig)
         end do
         fx = 10.0**(-fx*conc*path)
!
         BBx = blackbody(xp)
!
         BBxodd = BBxodd + BBx*Xcie(xp)
         BByodd = BByodd + BBx*Ycie(xp)
         BBzodd = BBzodd + BBx*Zcie(xp)

         fxodd = fxodd + BBxodd*fx
         fyodd = fyodd + BByodd*fx
         fzodd = fzodd + BBzodd*fx
!
       end do
!
! Adding the area of the function at each even-subinterval
!
       fxeven  = 0.0
       fyeven  = 0.0
       fzeven  = 0.0
!
       BBxeven = 0.0
       BByeven = 0.0
       BBzeven = 0.0
!
       xp = xa
!
       do i = 3, n-1, 2
!
         xp = xp + (i-1)*h
!
         fx = 0.0
         do  j = 1, nband
           fx = fx + real(Aabs)*inten(j)*(1.0E7/xp)                    &
                              *gauss((1.0E7/xp-1.0E7/freq(j)),1.0E7/sig)
         end do
         fx = 10.0**(-fx*conc*path)
!
         BBx = blackbody(xp)
!
         BBxeven = BBxeven + BBx*Xcie(xp)
         BByeven = BByeven + BBx*Ycie(xp)
         BBzeven = BBzeven + BBx*Zcie(xp)
!
         fxeven = fxeven + BBxeven*fx
         fyeven = fyeven + BByeven*fx
         fzeven = fzeven + BBzeven*fx
!
       end do
!
       return
       end subroutine simpitergauss
!
!======================================================================!
!
       subroutine simpiterlor(fxodd,fyodd,fzodd,fxeven,fyeven,fzeven,  &
                              BBxodd,BByodd,BBzodd,BBxeven,BByeven,    &
                              BBzeven,n,xa,h,width)
!
       use lineshape
!
       implicit none
!
! Input/Output variables
!
       real(kind=4),intent(out)  ::  fxodd     !
       real(kind=4),intent(out)  ::  fyodd     !
       real(kind=4),intent(out)  ::  fzodd     !
       real(kind=4),intent(out)  ::  fxeven    !
       real(kind=4),intent(out)  ::  fyeven    !
       real(kind=4),intent(out)  ::  fzeven    !
       real(kind=4),intent(out)  ::  BBxodd    !
       real(kind=4),intent(out)  ::  BByodd    !
       real(kind=4),intent(out)  ::  BBzodd    !
       real(kind=4),intent(out)  ::  BBxeven   !
       real(kind=4),intent(out)  ::  BByeven   !
       real(kind=4),intent(out)  ::  BBzeven   !
       real(kind=4),intent(in)   ::  xa        !
       real(kind=4),intent(in)   ::  h         !
       real(kind=4),intent(in)   ::  width     !
       integer,intent(in)        ::  n         !
!
! Local variables
!
       real(kind=4)              ::  xp        !
       real(kind=4)              ::  fx        !
       real(kind=4)              ::  BBx       !
       integer                   ::  i,j       !
!
! Computing the Simpson's rule sums for a Gaussian broadening
! -----------------------------------------------------------
!
! Adding the area of the function at each odd-subinterval
!
       fxodd  = 0.0
       fyodd  = 0.0
       fzodd  = 0.0
!
       BBxodd = 0.0
       BByodd = 0.0
       BBzodd = 0.0
!
       xp = xa
!
       do i = 2, n, 2
!
         xp = xp + (i-1)*h
!
         fx = 0.0
         do  j = 1, nband
           fx = fx + real(Aabs)*inten(j)*(1.0E7/xp)                    &
                               *lor((1.0E7/xp-1.0E7/freq(j)),1.0E7/width)
         end do
         fx = 10.0**(-fx*conc*path)
!
         BBx = blackbody(xp)
!
         BBxodd = BBxodd + BBx*Xcie(xp)
         BByodd = BByodd + BBx*Ycie(xp)
         BBzodd = BBzodd + BBx*Zcie(xp)

         fxodd = fxodd + BBxodd*fx
         fyodd = fyodd + BByodd*fx
         fzodd = fzodd + BBzodd*fx
!
       end do
!
! Adding the area of the function at each even-subinterval
!
       fxeven  = 0.0
       fyeven  = 0.0
       fzeven  = 0.0
!
       BBxeven = 0.0
       BByeven = 0.0
       BBzeven = 0.0
!
       xp = xa
!
       do i = 3, n-1, 2
!
         xp = xp + (i-1)*h
!
         fx = 0.0
         do  j = 1, nband
           fx = fx + real(Aabs)*inten(j)*(1.0E7/xp)                    &
                               *lor((1.0E7/xp-1.0E7/freq(j)),1.0E7/width)
         end do
         fx = 10.0**(-fx*conc*path)
!
         BBx = blackbody(xp)
!
         BBxeven = BBxeven + BBx*Xcie(xp)
         BByeven = BByeven + BBx*Ycie(xp)
         BBzeven = BBzeven + BBx*Zcie(xp)
!
         fxeven = fxeven + BBxeven*fx
         fyeven = fyeven + BByeven*fx
         fzeven = fzeven + BBzeven*fx
!
       end do
!
       return
       end subroutine simpiterlor
!
!======================================================================!
!
       subroutine simpitervoigt(fxodd,fyodd,fzodd,fxeven,fyeven,       &
                                fzeven,BBxodd,BByodd,BBzodd,BBxeven,   &
                                BByeven,BBzeven,n,xa,h,f,sig,eta)
!
       use lineshape
!
       implicit none
!
! Input/Output variables
!
       real(kind=4),intent(out)  ::  fxodd     !
       real(kind=4),intent(out)  ::  fyodd     !
       real(kind=4),intent(out)  ::  fzodd     !
       real(kind=4),intent(out)  ::  fxeven    !
       real(kind=4),intent(out)  ::  fyeven    !
       real(kind=4),intent(out)  ::  fzeven    !
       real(kind=4),intent(out)  ::  BBxodd    !
       real(kind=4),intent(out)  ::  BByodd    !
       real(kind=4),intent(out)  ::  BBzodd    !
       real(kind=4),intent(out)  ::  BBxeven   !
       real(kind=4),intent(out)  ::  BByeven   !
       real(kind=4),intent(out)  ::  BBzeven   !
       real(kind=4),intent(in)   ::  xa        !
       real(kind=4),intent(in)   ::  h         !
       real(kind=4),intent(in)   ::  f         !
       real(kind=4),intent(in)   ::  sig       !
       real(kind=4),intent(in)   ::  eta       !
       integer,intent(in)        ::  n         !
!
! Local variables
!
       real(kind=4)              ::  xp        !
       real(kind=4)              ::  fx        !
       real(kind=4)              ::  BBx       !
       integer                   ::  i,j       !
!
! Computing the Simpson's rule sums for a Gaussian broadening
! -----------------------------------------------------------
!
! Adding the area of the function at each odd-subinterval
!
       fxodd  = 0.0
       fyodd  = 0.0
       fzodd  = 0.0
!
       BBxodd = 0.0
       BByodd = 0.0
       BBzodd = 0.0
!
       xp = xa
!
       do i = 2, n, 2
!
         xp = xp + (i-1)*h
!
         fx = 0.0
         do  j = 1, nband
           fx = fx + real(Aabs)*inten(j)*(1.0E7/xp)                    &
                *voigt(1.0E7*(1.0/xp-1.0/freq(j)),1.0E7/f,1.0E7/sig,eta)
         end do
         fx = 10.0**(-fx*conc*path)
!
         BBx = blackbody(xp)
!
         BBxodd = BBxodd + BBx*Xcie(xp)
         BByodd = BByodd + BBx*Ycie(xp)
         BBzodd = BBzodd + BBx*Zcie(xp)

         fxodd = fxodd + BBxodd*fx
         fyodd = fyodd + BByodd*fx
         fzodd = fzodd + BBzodd*fx
!
       end do
!
! Adding the area of the function at each even-subinterval
!
       fxeven  = 0.0
       fyeven  = 0.0
       fzeven  = 0.0
!
       BBxeven = 0.0
       BByeven = 0.0
       BBzeven = 0.0
!
       xp = xa
!
       do i = 3, n-1, 2
!
         xp = xp + (i-1)*h
!
         fx = 0.0
         do  j = 1, nband
           fx = fx + real(Aabs)*inten(j)*(1.0E7/xp)                    &
                *voigt(1.0E7*(1.0/xp-1.0/freq(j)),1.0E7/f,1.0E7/sig,eta)
         end do
         fx = 10.0**(-fx*conc*path)
!
         BBx = blackbody(xp)
!
         BBxeven = BBxeven + BBx*Xcie(xp)
         BByeven = BByeven + BBx*Ycie(xp)
         BBzeven = BBzeven + BBx*Zcie(xp)
!
         fxeven = fxeven + BBxeven*fx
         fyeven = fyeven + BByeven*fx
         fzeven = fzeven + BBzeven*fx
!
       end do
!
       return
       end subroutine simpitervoigt
!
!======================================================================!
!
       subroutine precompinv(nspec,xspec,yspec,xa,h)
!
       use lineshape
!
       implicit none
!
! Input/Output variables
!
       real(kind=4),dimension(:),intent(out)  ::  xspec  !
       real(kind=4),dimension(:),intent(out)  ::  yspec  !
       real(kind=4),intent(in)                ::  xa     !  Lower integration limit
       real(kind=4),intent(in)                ::  h      !  Step length
       integer,intent(in)                     ::  nspec  !
!
! Local variables
!
       real(kind=4)                           ::  width  !
       real(kind=4)                           ::  ghwhm  !
       real(kind=4)                           ::  lhwhm  !
!
! Precomputing the spectrum for the given integration limits and step
! -------------------------------------------------------------------
!
       if ( ixunit .eq. 1 ) then
         width = 1.0E7/hwhm
         ghwhm = 1.0E7/gauhwhm
         lhwhm = 1.0E7/lorhwhm
       else
         width = hwhm
         ghwhm = gauhwhm
         lhwhm = lorhwhm
       end if
!
! Folding the oscillator strengths with broadening functions
!  (inverse spectrum method)
!
       if ( iline .eq. 1 ) then
         call gaussnm(nspec,xspec,yspec,xa,h,width,1)
       else if ( iline .eq. 2 ) then
         call lornm(nspec,xspec,yspec,xa,h,width,1)
       else
         call voigtnm(nspec,xspec,yspec,xa,h,ghwhm,lhwhm,1)
       end if
!
!~        if ( iline .eq. 1 ) then               !  Gaussian broadening
!~          call fgaussnm(nspec,xspec,yline,h,a,width,1)
!~        else if ( iline .eq. 2 ) then          !  Lorentzian broadening
!~          call flornm(nspec,xspec,yline,h,a,width,1)
!~        else                                   !  Voigt broadening
!~          call fvoigtnm(nspec,xspec,yline,h,a,ghwhm,lhwhm,1)
!~        end if
!
       return
       end subroutine precompinv
!
!======================================================================!
!
       subroutine precompillu(nspec,xspec,yspec)
!
       use lineshape
!
       implicit none
!
! Input/Output variables
!
       real(kind=4),dimension(:),intent(in)   ::  xspec  !
       real(kind=4),dimension(:),intent(out)  ::  yspec  !
       integer,intent(in)                     ::  nspec  !
!
! Local variables
!
       real(kind=4)                           ::  width  !
       real(kind=4)                           ::  ghwhm  !
       real(kind=4)                           ::  lhwhm  !
!
! Precomputing the spectrum for the given integration limits and step
! -------------------------------------------------------------------
!
       if ( ixunit .eq. 1 ) then
         width = 1.0E7/hwhm
         ghwhm = 1.0E7/gauhwhm
         lhwhm = 1.0E7/lorhwhm
       else
         width = hwhm
         ghwhm = gauhwhm
         lhwhm = lorhwhm
       end if
!
! Calculating the molar absorption coefficients
!  (standard illuminant method)
!
       if ( iline .eq. 1 ) then                     !  Gaussian broadening
         call pregaussnm(nspec,xspec,yspec,width,1)
       else if ( iline .eq. 2 ) then                !  Lorentzian broadening
         call prelornm(nspec,xspec,yspec,width,1)
       else                                         !  Voigt broadening
         call prevoigtnm(nspec,xspec,yspec,ghwhm,lhwhm,1)
       end if
!
       return
       end subroutine precompillu
!
!======================================================================!
!
! References
! ----------
!
! - http://www.brucelindbloom.com/
! - Michael S Tooms. The Bradford Colour Adaptation Transform Colour,
!    Reproduction in Electronic Imaging Systems, John Wiley & Sons, Ltd,
!    2015, Part V, 645-647
!
       subroutine adaptation(iadapt,X1,Y1,Z1,Xw,Yw,Zw,Xwref,Ywref,Zwref)
!
       implicit none
!
! Input/Output variables
!
       real(kind=4),intent(inout)   ::  X1,Y1,Z1     !  XYZ tristimulus values
       real(kind=4),intent(in)      ::  Xw,Yw,Zw     !  Source white point coordinates
       real(kind=4),intent(in)      ::  Xwref        !  Reference white point coordinate
       real(kind=4),intent(in)      ::  Ywref        !  Reference white point coordinate
       real(kind=4),intent(in)      ::  Zwref        !  Reference white point coordinate
       integer,intent(in)           ::  iadapt       !  Adaptation method
!
! Local variables
!
       real(kind=4)                 ::  a11,a12,a13  !  Adapatation matrix elements
       real(kind=4)                 ::  a21,a22,a23  !  Adapatation matrix elements
       real(kind=4)                 ::  a31,a32,a33  !  Adapatation matrix elements
!
       real(kind=4)                 ::  detA         !  Determinant of the matrix A
!
       real(kind=4)                 ::  b11,b12,b13  !  Inverse adaptation matrix elements
       real(kind=4)                 ::  b21,b22,b23  !  Inverse adaptation matrix elements
       real(kind=4)                 ::  b31,b32,b33  !  Inverse adaptation matrix elements
!
       real(kind=4)                 ::  m11,m12,m13  !  Transformation matrix elements
       real(kind=4)                 ::  m21,m22,m23  !  Transformation matrix elements
       real(kind=4)                 ::  m31,m32,m33  !  Transformation matrix elements
!
       real(kind=4)                 ::  rs,gs,bs     !  Source cone response domain
       real(kind=4)                 ::  rd,gd,bd     !  Destination cone response domain
!
! Transforming a source color into a destination color
! ----------------------------------------------------
!
! Going from a source color to a destination color by a linear
!  transformation which is dependent on the source reference white
!  and the destination reference white
!
       if ( iadapt .eq. 1 ) then
!
         a11 =  0.8951000
         a12 =  0.2664000
         a13 = -0.1614000
!
         a21 = -0.7502000
         a22 =  1.7135000
         a23 =  0.0367000

!
         a31 =  0.0389000
         a32 = -0.0685000
         a33 =  1.0296000
!
       else if ( iadapt .eq. 2 ) then
!
         a11 =  0.4002400
         a12 =  0.7076000
         a13 = -0.0808100
!
         a21 = -0.2263000
         a22 =  1.1653200
         a23 =  0.0457000
!
         a31 =  0.0000000
         a32 =  0.0000000
         a33 =  0.9182200
!
       else
!
         X1 = (Xwref/Xw)*X1
         Y1 = (Ywref/Yw)*Y1
         Z1 = (Zwref/Zw)*Z1
         return
!
       end if
!
! Computing the inverse of the cone response domains matrix
!
       b11 = a11*a33 - a23*a32
       b12 = a13*a32 - a12*a33
       b13 = a12*a23 - a13*a22
!
       b21 = a23*a31 - a21*a33
       b22 = a11*a33 - a13*a31
       b23 = a13*a21 - a11*a23
!
       b31 = a21*a32 - a22*a31
       b32 = a12*a31 - a11*a32
       b33 = a11*a22 - a12*a21
!
       detA = a11*b11 + a12*b21 + a13*b31
!
       b11 = b11/detA
       b12 = b12/detA
       b13 = b13/detA
!
       b21 = b21/detA
       b22 = b22/detA
       b23 = b23/detA
!
       b31 = b31/detA
       b32 = b32/detA
       b33 = b33/detA
!
! Transforming the reference whites from XYZ into a cone response domain
!
       rs = a11*Xw + a12*Yw + a13*Zw
       gs = a21*Xw + a22*Yw + a23*Zw
       bs = a31*Xw + a32*Yw + a33*Zw
!
       rd = a11*Xwref + a12*Ywref + a13*Zwref
       gd = a21*Xwref + a22*Ywref + a23*Zwref
       bd = a31*Xwref + a32*Ywref + a33*Zwref
!
! Scaling  the vector components by factors dependent upon both the 
!  source and destination reference whites
!
       rd = rd/rs
       gd = gd/gs
       bd = bd/bs
!
       rs = X1
       gs = Y1
       bs = Z1
!
! Scaling elements of the Ma matrix
!
       a11 = a11*rd
       a12 = a12*rd
       a13 = a13*rd
!
       a21 = a21*gd
       a22 = a22*gd
       a23 = a23*gd
!
       a31 = a31*bd
       a32 = a32*bd
       a33 = a33*bd
!
! Computing the transformation matrix (Ma^-1 * Ma')
!
       m11 = b11*a11 + b12*a21 + b13*a31
       m12 = b11*a12 + b12*a22 + b13*a32
       m13 = b11*a13 + b12*a23 + b13*a33
!
       m21 = b21*a11 + b22*a21 + b23*a31
       m22 = b21*a12 + b22*a22 + b23*a32
       m23 = b21*a13 + b22*a23 + b23*a33
!
       m31 = b31*a11 + b32*a21 + b33*a31
       m32 = b31*a12 + b32*a22 + b33*a32
       m33 = b31*a13 + b32*a23 + b33*a33
!
! Transform from cone response domains back to XYZ
! 
       X1 = m11*rs + m12*gs + m13*bs
       Y1 = m21*rs + m22*gs + m23*bs
       Z1 = m31*rs + m32*gs + m33*bs
!
       return
       end subroutine adaptation
!
!======================================================================!
!
! References
! ----------
!
! - https://en.wikipedia.org/wiki/Simpson%27s_rule
! - Shklov, N. (December 1960). "Simpson's Rule for Unequally Spaced
!    Ordinates". The American Mathematical Monthly. 67 (10): 1022.
!    doi:10.2307/2309244
!
       real(kind=4) function simpuneven(m,x,f,cie) result(Ics)
!
       implicit none
!
! Input/Output variables
!
       real(kind=4),dimension(m),intent(in)  ::  x      !
       real(kind=4),dimension(m),intent(in)  ::  f      !
       real(kind=4),external                 ::  cie    !
       integer,intent(in)                    ::  m      !
!
! Local variables
!
       real(kind=4),dimension(m)             ::  dx     ! Step length
       real(kind=4)                          ::  h0     !
       real(kind=4)                          ::  h1     !
       real(kind=4)                          ::  hsum   !
       real(kind=4)                          ::  hprod  !
       real(kind=4)                          ::  hdiv   !
       integer                               ::  n      !
       integer                               ::  i      !
!
! Integrating experimental spectrum using Composite Simpson's rule for
!  unevenly spaced intervals (irregulary spaced data)
!
       n = m - 1
!
       do i = 1, n
         dx(i) = f(i+1) - f(i)
       end do
!
       Ics = 0.0
!
       do i = 2, n, 2
!
         h0 = dx(i-1)
         h1 = dx(i)
!
         hsum  = h0 + h1
         hprod = h0*h1
         hdiv  = h1/h0
!
         Ics = Ics + hsum/6.0*( (2 - hdiv)*f(i-1)*cie(x(i-1))          &
                               + (hsum**2/hprod)*f(i)*cie(x(i))        &
                               + (2 - 1.0/hdiv)*f(i+1)*cie(x(i+1)) )
!
       end do
!
       if ( mod(n,2) .ne. 0 ) then
!
         h0 = dx(n-1)
         h1 = dx(n-2)
!
         Ics = Ics + f(n)*cie(x(n))*(2.0*h0**2 + 3.0*h0*h1)/           &
                                                      (6.0*(h1 + h0))  &
                   + f(n-1)*cie(x(n-1))*(h0**2 + 3.0*h0*h1)/(6.0*h1)   &
                   - f(n-2)*cie(x(n-2))*h0**3/(6.0*h1*(h1 + h0))
!
       end if
!
       return
       end function simpuneven
!
!======================================================================!
!
       real(kind=4) function simpeven(m,x,f,h,cie) result(Ics)
!
       implicit none
!
! Input/Output variables
!
       real(kind=4),dimension(m),intent(in)  ::  x      !
       real(kind=4),dimension(m),intent(in)  ::  f      !
       real(kind=4)                          ::  h      !  Step length
       real(kind=4),external                 ::  cie    !
       integer,intent(in)                    ::  m      !
!
! Local variables
!
       real(kind=4)                          ::  fa,fb  !  Value of the function at the integration limits
       real(kind=4)                          ::  seven  !  Auxliary variable
       real(kind=4)                          ::  sodd   !  Auxliary variable
       integer                               ::  n      !
       integer                               ::  i      !
!
! Integrating spectrum using Composite Simpson's method
!
       if ( mod(m,2) .eq. 0 ) then
         stop 'Simpson method needs an even number of subintervals'
       end if
!
       fa = f(1)*cie(x(1))
       fb = f(m)*cie(x(m))
!
       n = m - 1
!
! Adding the area of the function at each odd-subinterval
!
       sodd  = 0.0
!
       do i = 2, n, 2
         sodd = sodd + f(i)*cie(x(i))
       end do
!
! Adding the area of the function at each even-subinterval
!
       seven = 0.0
!
       do i = 3, n-1, 2
         seven = seven + f(i)*cie(x(i))
       end do
!
! Calculating the approximated area (integral)
!
       Ics = h/3.0*(fa + fb + 2*seven + 4*sodd)
!
       return
       end function simpeven
!
!======================================================================!
!
       real(kind=4) function blackbody(xp) result(f)
!
       use variables
       use parameters 
!
       implicit none
!
! Input/Output variables
!
       real(kind=4),intent(in)  ::  xp   !
!
!  Computing relative spectral density of electromagnetic radiation 
!   emitted by a black body according to Planck's law for a reference
!   wavelength of 560 nm
!
       f = (560.d0/xp)**5*(exp(cc2/CCT) - 1)/(exp(c2*1.0E9/xp/CCT) - 1)
!
       return
       end function blackbody
!
!======================================================================!
!
! References
! ----------
!
! - Chris Wyman, Peter-Pike Sloan, and Peter Shirley, Simple Analytic
!    Approximations to the CIE XYZ Color Matching Functions, Journal
!    of Computer Graphics Techniques (JCGT), vol. 2, no. 2, 1-11, 2013.
!    http://jcgt.org/published/0002/02/01/
!
       real(kind=4) function Xcie(x)
!
       implicit none
!
! Input/Output variables
!
       real(kind=4),intent(in)  ::  x  !
!
! Computing the X CIE1931 color-matching function using the multi-lobe
!  Gaussian fit
!
       Xcie = Gcie(x,x0a,x0b,x0g,x0d) + Gcie(x,x1a,x1b,x1g,x1d)        &
                                               + Gcie(x,x2a,x2b,x2g,x2d)
!
       return
       end function Xcie
!
!======================================================================!
!
! References
! ----------
!
! - Chris Wyman, Peter-Pike Sloan, and Peter Shirley, Simple Analytic
!    Approximations to the CIE XYZ Color Matching Functions, Journal
!    of Computer Graphics Techniques (JCGT), vol. 2, no. 2, 1-11, 2013.
!    http://jcgt.org/published/0002/02/01/
!
       real(kind=4) function Ycie(x)
!
       implicit none
!
! Input/Output variables
!
       real(kind=4),intent(in)  ::  x  !
!
! Computing the Y CIE1931 color-matching function using the multi-lobe
!  Gaussian fit
!
       Ycie = Gcie(x,y0a,y0b,y0g,y0d) + Gcie(x,y1a,y1b,y1g,y1d)
!
       return
       end function Ycie
!
!======================================================================!
!
! References
! ----------
!
! - Chris Wyman, Peter-Pike Sloan, and Peter Shirley, Simple Analytic
!    Approximations to the CIE XYZ Color Matching Functions, Journal
!    of Computer Graphics Techniques (JCGT), vol. 2, no. 2, 1-11, 2013.
!    http://jcgt.org/published/0002/02/01/
!
       real(kind=4) function Zcie(x)
!
       implicit none

!
! Input/Output variables
!
       real(kind=4),intent(in)  ::  x  !
!
! Computing the Z CIE1931 color-matching function using the multi-lobe
!  Gaussian fit
!
       Zcie = Gcie(x,z0a,z0b,z0g,z0d) + Gcie(x,z1a,z1b,z1g,z1d)
!
       return
       end function Zcie
!
!======================================================================!
!
       real(kind=4) function Gcie(x,xa,xb,xg,xd)
!
       implicit none
!
! Input/Output variables
!
       real(kind=4),intent(in)  ::  x   !
       real(kind=4),intent(in)  ::  xa  !
       real(kind=4),intent(in)  ::  xb  !
       real(kind=4),intent(in)  ::  xg  !
       real(kind=4),intent(in)  ::  xd  !
!
!  Computing multi-lobe, piecewise Gaussian function
!
       if ( x .lt. xb ) then
         Gcie = xa*exp(-1.0/2.0*((x-xb)*xg)**2)
       else
         Gcie = xa*exp(-1.0/2.0*((x-xb)*xd)**2)
       end if
!
       return
       end function Gcie
!
!======================================================================!
!
! References
! ----------
!
! - http://www.brucelindbloom.com/
!
       subroutine rgbmodel(refillu,gam,Xwref,Ywref,Zwref,xr,yr,zr,     &
                           xg,yg,zg,xb,yb,zb)
!
       use variables
!
       implicit none
!
! Input/Output variables
!
       character(len=10),intent(out)  ::  refillu   !
       real(kind=4),intent(out)       ::  Xwref     !  Reference illuminant white
       real(kind=4),intent(out)       ::  Ywref     !  Reference illuminant white
       real(kind=4),intent(out)       ::  Zwref     !  Reference illuminant white       real(kind=4)                              ::  xr,yr,zr  !  RGB system red point chromatic coordinates
       real(kind=4),intent(out)       ::  xr,yr,zr  !  RGB system red point chromatic coordinates
       real(kind=4),intent(out)       ::  xg,yg,zg  !  RGB system green point chromatic coordinates
       real(kind=4),intent(out)       ::  xb,yb,zb  !  RGB system blue point chromatic coordinates
       real(kind=4),intent(out)       ::  gam       !
!
! Setting up RGB system information
! .................................
!
       select case(trim(model))
!
! Reference white D65
!
         case ('adobe','bruce','hdtv','ebu','smpte','srgb')
!
           Xwref = 0.95047
           Ywref = 1.00000 
           Zwref = 1.08883
!
           gam = 2.2d0
!
           refillu = 'd65'
!
         case ('apple')
!
           Xwref = 0.95047
           Ywref = 1.00000 
           Zwref = 1.08883
!
           gam = 1.8d0
!
           refillu = 'd65'
!
! Reference white D50
!
         case ('best','beta','don','ekta','widegamut')
!
           Xwref = 0.96422
           Ywref = 1.00000 
           Zwref = 0.82521
!
           gam = 2.2d0
!
           refillu = 'd50'
!
         case ('colormatch','eci','prophoto')
!
           Xwref = 0.96422
           Ywref = 1.00000 
           Zwref = 0.82521
!
           gam = 1.8d0
!
           refillu = 'd50'
!
! Reference white E
!
         case ('cie')
!
           xr = 0.7350
           yr = 0.2650
!
           xg = 0.2740
           yg = 0.7170
!
           xb = 0.1670
           yb = 0.0090
!
           Xwref = 1.00000
           Ywref = 1.00000 
           Zwref = 1.00000
!
           gam = 2.2d0
!
           refillu = 'e'
!
! Reference white C
!
         case ('ntsc')
!
           xr = 0.6700
           yr = 0.3300
!
           xg = 0.2100
           yg = 0.7100
!
           xb = 0.1400
           yb = 0.0800
!
           Xwref = 0.98074 
           Ywref = 1.00000 
           Zwref = 1.18232
!
           gam = 2.2d0
!
           refillu = 'c'
!
       end select
!
       select case(trim(model))
!
         case ('adobe')
!
           xr = 0.6400
           yr = 0.3300
!
           xg = 0.2100
           yg = 0.7100
!
           xb = 0.1500
           yb = 0.0600
!
         case ('apple')
!
           xr = 0.6250
           yr = 0.3400
!
           xg = 0.2800
           yg = 0.5950
!
           xb = 0.1550
           yb = 0.0700
!
         case ('best')
!
           xr = 0.7347
           yr = 0.2653
!
           xg = 0.2150
           yg = 0.7750
!
           xb = 0.1300
           yb = 0.0350
!
         case ('beta')
!
           xr = 0.6888
           yr = 0.3112
!
           xg = 0.1986
           yg = 0.7551
!
           xb = 0.1265
           yb = 0.0352
!
         case ('bruce')
!
           xr = 0.6400
           yr = 0.3300
!
           xg = 0.2800
           yg = 0.6500
!
           xb = 0.1500
           yb = 0.0600
!
         case ('colormatch')
!
           xr = 0.6300
           yr = 0.3400
!
           xg = 0.2950
           yg = 0.6050
!
           xb = 0.1500
           yb = 0.0750
!
         case ('don')
!
           xr = 0.6960
           yr = 0.3000
!
           xg = 0.2150
           yg = 0.7650
!
           xb = 0.1300
           yb = 0.0350
!
         case ('eci')
!
           xr = 0.6700
           yr = 0.3300
!
           xg = 0.2100
           yg = 0.7100
!
           xb = 0.1400
           yb = 0.0800
!
         case ('ekta')
!
           xr = 0.6950
           yr = 0.3050
!
           xg = 0.2600
           yg = 0.7000
!
           xb = 0.1100
           yb = 0.0050
!
         case ('hdtv')
!
           xr = 0.6700
           yr = 0.3300
!
           xg = 0.2100 
           yg = 0.7100
!
           xb = 0.1500
           yb = 0.0600
!
         case ('ebu')
!
           xr = 0.6400
           yr = 0.3300
!
           xg = 0.2900
           yg = 0.6000
!
           xb = 0.1500
           yb = 0.0600
!
         case ('prophoto')
!
           xr = 0.7347
           yr = 0.2653
!
           xg = 0.1596
           yg = 0.8404
!
           xb = 0.0366
           yb = 0.0001
!
         case ('smpte')
!
           xr = 0.6300
           yr = 0.3400
!
           xg = 0.3100
           yg = 0.5950
!
           xb = 0.1550
           yb = 0.0700
!
         case ('srgb')
!
           xr = 0.6400
           yr = 0.3300
!
           xg = 0.3000
           yg = 0.6000
!
           xb = 0.1500
           yb = 0.0600
!
         case ('widegamut')
!
           xr = 0.7350
           yr = 0.2650
!
           xg = 0.1150
           yg = 0.8260
!
           xb = 0.1570
           yb = 0.0180
!
       end select
!
       zr = 1.0000 - xr - yr
       zg = 1.0000 - xg - yg
       zb = 1.0000 - xb - yb
!
       return
       end subroutine rgbmodel
!
!======================================================================!
!
       end module colour
!
!======================================================================!
