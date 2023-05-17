!======================================================================!
!
       module lineshape
!
       use variables
       use parameters
!
       implicit none
!
       contains
!
!======================================================================!
!
       subroutine gausscm(nspec,xspec,yspec,x0,xstep,hwhm,n)
!
       implicit none
!
! Input/Output variables
!
       real(kind=4),dimension(nspec),intent(out)  ::  xspec  !
       real(kind=4),dimension(nspec),intent(out)  ::  yspec  !
       real(kind=4),intent(in)                    ::  xstep  !
       real(kind=4),intent(in)                    ::  x0     !
       real(kind=4),intent(in)                    ::  hwhm   !
       integer,intent(in)                         ::  nspec  !
       integer,intent(in)                         ::  n      !
!
! Local variables
!
       real(kind=4)                               ::  x      !
       real(kind=4)                               ::  sig    !
       integer                                    ::  i,j    !  Indexes
!
! Computing spectral lineshape
!
       sig = hwhm/sqrt(2.0*log(2.0))
!
       x = x0
!
       do i = 1, nspec
!
         xspec(i) = x
!
         yspec(i) = 0.0
!
         do  j = 1, nband
           yspec(i) = yspec(i) + real(Aabs)*inten(j)*x**n              &
                                                   *gauss(x-freq(j),sig)
!~            yline(i) = yline(i) + eps*inten(j)*x**n                     &
!~                                                  *gauss(x-freq(j),sig)
!~            yline(i) = yline(i) + Aabs*inten(j)*(1.0E7/x)*gauss((1.0E7/x)-(1.0E7/freq(j)),1.0E7/sig)
         end do
!
         x = x + xstep
!
       end do
!
       return
       end subroutine gausscm
!
!======================================================================!
!
       subroutine lorcm(nspec,xspec,yspec,x0,xstep,hwhm,n)
!
       implicit none
!
! Input/Output variables
!
       real(kind=4),dimension(nspec),intent(out)  ::  xspec  !
       real(kind=4),dimension(nspec),intent(out)  ::  yspec  !
       real(kind=4),intent(in)                    ::  xstep  !
       real(kind=4),intent(in)                    ::  x0     !
       real(kind=4),intent(in)                    ::  hwhm   !
       integer,intent(in)                         ::  nspec  !
       integer,intent(in)                         ::  n      !
!
! Local variables
!
       real(kind=4)                               ::  x      !
       integer                                    ::  i,j    !  Indexes
!
! Computing spectral lineshape
!
       x = x0
!
       do i = 1, nspec
!
         xspec(i) = x
!
         yspec(i) = 0.0
!
         do  j = 1, nband
!~            yline(i) = yline(i) + eps*inten(j)*x**n*lor(x-freq(j),hwhm)
           yspec(i) = yspec(i) + real(Aabs)*inten(j)*x**n              &
                                                    *lor(x-freq(j),hwhm)
!~            yline(i) = yline(i) + Aabs*inten(j)*(1.0E7/x)*lor((1.0E7/x)-(1.0E7/freq(j)),1.0E7/hwhm)
         end do
!
         x = x + xstep
!
       end do
!
       return
       end subroutine lorcm
!
!======================================================================!
!
       subroutine voigtcm(nspec,xspec,yspec,x0,xstep,fg,fl,n)
!
       implicit none
!
! Input/Output variables
!
       real(kind=4),dimension(nspec),intent(out)  ::  xspec  !
       real(kind=4),dimension(nspec),intent(out)  ::  yspec  !
       real(kind=4),intent(in)                    ::  xstep  !
       real(kind=4),intent(in)                    ::  x0     !
       real(kind=4),intent(in)                    ::  fg     !
       real(kind=4),intent(in)                    ::  fl     !
       integer,intent(in)                         ::  nspec  !
       integer,intent(in)                         ::  n      !
!
! Local variables
!
       real(kind=4)                               ::  x      !
       real(kind=4)                               ::  f      !  Total HWHM
       real(kind=4)                               ::  sig    !
       real(kind=4)                               ::  eta    !
       integer                                    ::  i,j    !  Indexes
!
! Computing spectral lineshape
!
       f = (fg**5 + 2.69269*fg**4*fl + 2.42843*fg**3*fl*2              &
            + 4.47163*fg**2*fl**3 + 0.07842*fg*fl**4 + fl**5)**(1.0/5.0)
!
       eta = 1.36603*(fl/f) - 0.47719*(fl/f)**2.0 + 0.11116*(fl/f)**3.0
!
       sig = f/sqrt(2.0*log(2.0))
!
       x = x0
!
       do i = 1, nspec
!
         xspec(i) = x
!
         yspec(i) = 0.0
!
         do  j = 1, nband
           yspec(i) = yspec(i) + real(Aabs)*inten(j)*x**n              &
                                             *voigt(x-freq(j),f,sig,eta)
         end do
!
         x = x + xstep
!
       end do
!
       return
       end subroutine voigtcm
!
!======================================================================!
!
       subroutine gaussnm(nspec,xspec,yspec,x0,xstep,hwhm,n)
!
       implicit none
!
! Input/Output variables
!
       real(kind=4),dimension(nspec),intent(out)  ::  xspec  !
       real(kind=4),dimension(nspec),intent(out)  ::  yspec  !
       real(kind=4),intent(in)                    ::  xstep  !
       real(kind=4),intent(in)                    ::  x0     !
       real(kind=4),intent(in)                    ::  hwhm   !
       integer,intent(in)                         ::  nspec  !
       integer,intent(in)                         ::  n      !
!
! Local variables
!
       real(kind=4)                               ::  x      !
       real(kind=4)                               ::  sig    !
       integer                                    ::  i,j    !  Indexes
!
! Computing spectral lineshape
!
       sig = hwhm*sqrt(2.0*log(2.0))
!
       x = x0
!
       do i = 1, nspec
!
         xspec(i) = x
!
         yspec(i) = 0.0
!
         do  j = 1, nband
           yspec(i) = yspec(i) + real(Aabs)*inten(j)*(1.0E7/x)**n      &
                             *gauss(1.0E7*(1.0/x-1.0/freq(j)),1.0E7/sig)
         end do
!
         x = x + xstep
!
       end do
!
       return
       end subroutine gaussnm
!
!======================================================================!
!
       subroutine lornm(nspec,xspec,yspec,x0,xstep,hwhm,n)
!
       implicit none
!
! Input/Output variables
!
       real(kind=4),dimension(nspec),intent(out)  ::  xspec  !
       real(kind=4),dimension(nspec),intent(out)  ::  yspec  !
       real(kind=4),intent(in)                    ::  xstep  !
       real(kind=4),intent(in)                    ::  x0     !
       real(kind=4),intent(in)                    ::  hwhm   !
       integer,intent(in)                         ::  nspec  !
       integer,intent(in)                         ::  n      !
!
! Local variables
!
       real(kind=4)                               ::  x      !
       integer                                    ::  i,j    !  Indexes
!
! Computing spectral lineshape
!
       x = x0
!
       do i = 1, nspec
!
         xspec(i) = x
!
         yspec(i) = 0.0
!
         do  j = 1, nband
           yspec(i) = yspec(i) + real(Aabs)*inten(j)*(1.0E7/x)**n      &
                              *lor(1.0E7*(1.0/x-1.0/freq(j)),1.0E7/hwhm)
         end do
!
         x = x + xstep
!
       end do
!
       return
       end subroutine lornm
!
!======================================================================!
!
       subroutine voigtnm(nspec,xspec,yspec,x0,xstep,fg,fl,n)
!
       implicit none
!
! Input/Output variables
!
       real(kind=4),dimension(nspec),intent(out)  ::  xspec  !
       real(kind=4),dimension(nspec),intent(out)  ::  yspec  !
       real(kind=4),intent(in)                    ::  xstep  !
       real(kind=4),intent(in)                    ::  x0     !
       real(kind=4),intent(in)                    ::  fg     !
       real(kind=4),intent(in)                    ::  fl     !
       integer,intent(in)                         ::  nspec  !
       integer,intent(in)                         ::  n      !
!
! Local variables
!
       real(kind=4)                               ::  x      !
       real(kind=4)                               ::  f      ! Total HWHM
       real(kind=4)                               ::  sig    !
       real(kind=4)                               ::  eta    !
       integer                                    ::  i,j    !  Indexes
!
! Computing spectral lineshape
!
       f = (fg**5 + 2.69269*fg**4*fl + 2.42843*fg**3*fl*2              &
            + 4.47163*fg**2*fl**3 + 0.07842*fg*fl**4 + fl**5)**(1.0/5.0)
!
       eta = 1.36603*(fl/f) - 0.47719*(fl/f)**2.0 + 0.11116*(fl/f)**3.0
!
       sig = f*sqrt(2.0*log(2.0))
!
       x = x0
!
       do i = 1, nspec
!
         xspec(i) = x
!
         yspec(i) = 0.0
!
         do  j = 1, nband
           yspec(i) = yspec(i) + real(Aabs)*inten(j)*(1.0E7/x)**n      &
                 *voigt(1.0E7*(1.0/x-1.0/freq(j)),1.0E7/f,1.0E7/sig,eta)
         end do
!
         x = x + xstep
!
       end do
!
       return
       end subroutine voigtnm
!
!======================================================================!
!
       subroutine pregaussnm(nspec,xspec,yspec,hwhm,n)
!
       implicit none
!
! Input/Output variables
!
       real(kind=4),dimension(nspec),intent(in)   ::  xspec  !
       real(kind=4),dimension(nspec),intent(out)  ::  yspec  !
       real(kind=4),intent(in)                    ::  hwhm   !
       integer,intent(in)                         ::  nspec  !
       integer,intent(in)                         ::  n      !
!
! Local variables
!
       real(kind=4)                               ::  x      !
       real(kind=4)                               ::  sig    !
       integer                                    ::  i,j    !  Indexes
!
! Computing spectral lineshape
!
       sig = hwhm*sqrt(2.0*log(2.0))
!
       do i = 1, nspec
!
         x = xspec(i)
!
         yspec(i) = 0.0
!
         do  j = 1, nband
           yspec(i) = yspec(i) + real(Aabs)*inten(j)*(1.0E7/x)**n      &
                             *gauss(1.0E7*(1.0/x-1.0/freq(j)),1.0E7/sig)
         end do
!
       end do
!
       return
       end subroutine pregaussnm
!
!======================================================================!
!
       subroutine prelornm(nspec,xspec,yspec,hwhm,n)
!
       implicit none
!
! Input/Output variables
!
       real(kind=4),dimension(nspec),intent(in)   ::  xspec  !
       real(kind=4),dimension(nspec),intent(out)  ::  yspec  !
       real(kind=4),intent(in)                    ::  hwhm   !
       integer,intent(in)                         ::  nspec  !
       integer,intent(in)                         ::  n      !
!
! Local variables
!
       real(kind=4)                               ::  x      !
       integer                                    ::  i,j    !  Indexes
!
! Computing spectral lineshape
!
       do i = 1, nspec
!
         x = xspec(i)
!
         yspec(i) = 0.0
!
         do  j = 1, nband
           yspec(i) = yspec(i) + real(Aabs)*inten(j)*(1.0E7/x)**n      &
                              *lor(1.0E7*(1.0/x-1.0/freq(j)),1.0E7/hwhm)
         end do
!
       end do
!
       return
       end subroutine prelornm
!
!======================================================================!
!
       subroutine prevoigtnm(nspec,xspec,yspec,fg,fl,n)
!
       implicit none
!
! Input/Output variables
!
       real(kind=4),dimension(nspec),intent(in)   ::  xspec  !
       real(kind=4),dimension(nspec),intent(out)  ::  yspec  !
       real(kind=4),intent(in)                    ::  fg     !
       real(kind=4),intent(in)                    ::  fl     !
       integer,intent(in)                         ::  nspec  !
       integer,intent(in)                         ::  n      !
!
! Local variables
!
       real(kind=4)                               ::  x      !
       real(kind=4)                               ::  f      ! Total HWHM
       real(kind=4)                               ::  sig    !
       real(kind=4)                               ::  eta    !
       integer                                    ::  i,j    !  Indexes
!
! Computing spectral lineshape
!
       f = (fg**5 + 2.69269*fg**4*fl + 2.42843*fg**3*fl*2              &
            + 4.47163*fg**2*fl**3 + 0.07842*fg*fl**4 + fl**5)**(1.0/5.0)
!
       eta = 1.36603*(fl/f) - 0.47719*(fl/f)**2.0 + 0.11116*(fl/f)**3.0
!
       sig = f*sqrt(2.0*log(2.0))
!
       do i = 1, nspec
!
         x = xspec(i)
!
         yspec(i) = 0.0
!
         do  j = 1, nband
           yspec(i) = yspec(i) + real(Aabs)*inten(j)*(1.0E7/x)**n      &
                 *voigt(1.0E7*(1.0/x-1.0/freq(j)),1.0E7/f,1.0E7/sig,eta)
         end do
!
       end do
!
       return
       end subroutine prevoigtnm
!
!======================================================================!
!
       subroutine fgaussnm(nspec,xspec,yspec,x0,xstep,hwhm,n)
!
       implicit none
!
! Input/Output variables
!
       real(kind=4),dimension(nspec),intent(out)  ::  xspec  !
       real(kind=4),dimension(nspec),intent(out)  ::  yspec  !
       real(kind=4),intent(in)                    ::  xstep  !
       real(kind=4),intent(in)                    ::  x0     !
       real(kind=4),intent(in)                    ::  hwhm   !
       integer,intent(in)                         ::  nspec  !
       integer,intent(in)                         ::  n      !
!
! Local variables
!
       real(kind=4)                               ::  x      !
       real(kind=4)                               ::  sig    !
       integer                                    ::  i,j    !  Indexes
!
! Computing spectral lineshape
!
       sig = hwhm*sqrt(2.0*log(2.0))
!
       x = x0
!
       do i = 1, nspec
!
         xspec(i) = x
!
         yspec(i) = 0.0
!
         do  j = 1, nband
           yspec(i) = yspec(i) + real(Aabs)*inten(j)*(1.0E7/x)**n      &
                             *gauss(1.0E7*(1.0/x-1.0/freq(j)),1.0E7/sig)
         end do
!
         x = x + xstep
!
       end do
!
       return
       end subroutine fgaussnm
!
!======================================================================!
!
       real(kind=4) function gauss(dw,sig)
!
       implicit none
!
! Input/Output variables
!
       real(kind=4),intent(in)  ::  dw   !
       real(kind=4),intent(in)  ::  sig  !
!
       gauss = exp(-dw**2/2/sig**2)/sig/sqrt(2.0*real(pi))
!
       end function gauss
!
!======================================================================!
!
       real(kind=4) function lor(dw,hwhm)
!
       implicit none
!
! Input/Output variables
!
       real(kind=4),intent(in)  ::  dw    !
       real(kind=4),intent(in)  ::  hwhm  !
!
       lor = 1.0/( real(pi)*hwhm*(1.0 + (dw/hwhm)**2) )
!
       return
       end function lor
!
!======================================================================!
!
       real(kind=4) function voigt(dw,hwhm,sig,eta)
!
       implicit none
!
! Input/Output variables
!
       real(kind=4),intent(in)  ::  dw    !
       real(kind=4),intent(in)  ::  hwhm  !
       real(kind=4),intent(in)  ::  sig   !
       real(kind=4),intent(in)  ::  eta   !
!
       voigt = eta*lor(dw,hwhm) + (1.0 - eta)*gauss(dw,sig)
!
       return
       end function voigt
!
!======================================================================!
!
       subroutine deltacm()
!
       use sorting
!
       implicit none
!
! Local variables
!
       integer             ::  i
!
! Computing stick spectrum
!
       do i = 1, nstick
         xstick(i) = xmin + (i-1)*dx
         ystick(i) = 0.0
       end do
!
       xstick(npts+1:npts+nband) = freq(:nband)
       xstick(npts+nband+1:)     = freq(:nband)!~        ystick(npts+1:) = Aabs*freq(:nband)*inten(:nband)  !  Absorption coefficient
!~        ystick(npts+1:) = eps*freq(:nband)*inten(:nband)   !  Absorption coefficient (Requena)
!~        ystick(npts+1:) = fnm*freq(:nband)*inten(:nband)   !  Oscillator strength (Requena)
       ystick(npts+1:npts+nband) = 1054.94*fosc(:nband)                !  Integrated absorption coefficient
       ystick(npts+nband+1:)     = 0.0
!
!~ write(*,*) 'BUENO',1054.94*fosc(:nband)
!~ write(*,*) ' TEST',Bnm*Na*hplanck/log(10.0)/clight*freq(:nband)*inten(:nband)
!~ write(*,*) 'CONTR',Aabs*freq(:nband)*inten(:nband)
!
       call sinsort(nstick,xstick,ystick)
!
       return
       end subroutine deltacm
!
!======================================================================!
!
       subroutine deltanm()
!
       use sorting
!
       implicit none
!
! Local variables
!
       integer             ::  i
!
! Computing stick spectrum
!
       do i = 1, nstick
         xstick(i) = xmin + (i-1)*dx
         ystick(i) = 0.0
       end do
!
       xstick(npts+1:npts+nband) = freq(:nband)
       xstick(npts+nband+1:)     = freq(:nband)
!~        ystick(npts+1:) = Aabs*freq(:nband)*inten(:nband)  !  Absorption coefficient
!~        ystick(npts+1:) = eps*freq(:nband)*inten(:nband)   !  Absorption coefficient (Requena)
!~        ystick(npts+1:) = fnm*freq(:nband)*inten(:nband)   !  Oscillator strength (Requena)
       ystick(npts+1:npts+nband) = 1054.94*fosc(:nband)      !  Integrated absorption coefficient
       ystick(npts+nband+1:)     = 0.0
!
!~ write(*,*) 'BUENO',1054.94*fosc(:nband)
!~ write(*,*) ' TEST',Bnm*Na*hplanck/log(10.0)/clight*freq(:nband)*inten(:nband)
!~ write(*,*) 'CONTR',Aabs*freq(:nband)*inten(:nband)
!
       call sinsort(nstick,xstick,ystick)
!
       return
       end subroutine deltanm
!
!======================================================================!
!
       end module lineshape
!
!======================================================================!
