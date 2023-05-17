!==========================================================
!
       module print_plot
!
       use DISLIN
!
       implicit none
!
       contains
!
!======================================================================!
!
       subroutine initialize()
!
       use variables
       use lineshape
!
       implicit none
!
! Allocating memory
!
       if ( ixunit .eq. 1 ) then
         npts = int((xmaxcm-xmincm)/dxcm) + 1
       else
         npts = int((xmaxnm-xminnm)/dxnm) + 1
       end if
!
       nstick = npts + nband*2
!
       allocate(xstick(nstick),ystick(nstick))
       allocate(xline(npts),yline(npts))
!
! Computing stick spectra
!
       if ( ixunit .eq. 1 ) then
         call deltacm()
       else
         call deltanm()
       end if
!
       if ( iyunit .eq. 2 ) then
         ystick(:) = ystick(:)*conc*path
       else if ( iyunit .eq. 3 ) then
         ystick(:) = 10.0**(-ystick(:)*cold*lold+2)
       end if
!
       if ( (ilog.eq.2) .and. (inorm.eq.1) ) then
         ystick(:) = log10(ystick)/ylinemaxl
       else if ( ilog .eq. 2 ) then
         ystick(:) = log10(ystick)
       else if ( inorm .eq. 1 ) then
         ystick(:) = ystick(:)/ylinemax
       end if
!
! Computing convoluted spectra
!
       call convolute()
!
! Calculating initial limits of the spectrum
!
       call oldmargin()
!
       return
       end subroutine initialize
!
!======================================================================!
!
       subroutine iniplot(id)
!
       use variables
!
       implicit none
!
! Input/Output variables
!
       integer,intent(in) :: id
!
! Local variables
!
       integer            ::  ic
!
! Defines an external graphics window for X11 and Windows displays
!
       CALL SETXID(idplot,'WIDGET')
!
! Returns the selected element of a list widget
!
       if ( iunits .eq. 1 ) then
         CALL GWGLIS(idlog,ilog)
         CALL GWGLIS(idxunit,ixunit)
         CALL GWGLIS(idyunit,iyunit)
!  Transmittance units do not admit normalization/logarithmic scale
!~          if ( iyunit .eq. 3 ) then
!~            CALL SWGATT(idnorm,'INACTIVE','STATUS')
!~            CALL SWGATT(idlog,'INACTIVE','STATUS')
!~          else
!~            CALL SWGATT(idnorm,'ACTIVE','STATUS')
!~            CALL SWGATT(idlog,'ACTIVE','STATUS')
!~          end if
       end if
!
! Select metafile format (`XWIN`, `PDF`, ...).
!
       if ( id .eq. idexport ) then
         CALL METAFL('PNG')
       else
         CALL METAFL('XWIN')
       end if
!
! Set background to white, foreground to black
!
       CALL SCRMOD('REVERS')
!
! Initialise DISLIN
! -----------------
!
       CALL DISINI()
!
! The parameters in GRAF will be calculated automatically by DISLIN
!
!~        CALL SETSCL(xscl,2,'X')
!~        CALL SETSCL(yscl,2,'Y')
!
! Calculates parameters for GRAF from a minimum and maximum of data
!  values
!  Allows to extend the calculated axis limits to a full axis step
!
       CALL GAXPAR(xscl(1),xscl(2),'EXTEND','X',                       &  ! FLAG: if GAXPAR called from NEWMARGIN use 'NOEXTEND'
                   xminaux,xmaxaux,xor,xstp,nxdig)
       CALL GAXPAR(yscl(1),yscl(2),'EXTEND','Y',                       &
                   yminaux,ymaxaux,yor,ystp,nydig)
!
       if ( int((xmaxaux-xminaux)/xstp) .gt. 5 )                       &
                                              xstp = (xmaxaux-xminaux)/5
!
! Clears the screen, a graphics window or the page of a raster format
!  such as TIFF, PNG, PPM and BMP
!
       CALL ERASE()
!
! Add border around page
!
       CALL PAGERA()
!
! Select single stroke font (SIMPLX).
!
!~        CALL SIMPLX()
!
! Sets a complex font
!
       CALL COMPLX()
!
! Suppresses listing of data points that lie outside of the axis scaling
!
       CALL NOCHEK()
!
! Defines the thickness of curves
!
       CALL THKCRV(2)
!
! Defines the position of ticks
!
       CALL TICPOS('LABELS','XY')
!
! Sets the number of ticks
!
       CALL TICKS(2,'XY')
!
! Removes a part of an axis or a complete axis from an axis system
!  Used to locate ticks only with labels
!
       CALL SETGRF('NAME','NAME','LINE','LINE')
!
! Determines the position of axis systems
!
       CALL AXSPOS(400,1900)
!
! Defines axis lengths for a 2-D axis system (sets axis size)
!
!~        CALL AXSLEN(2000,1700)
       CALL AXSLEN(1900,1700)
!
! Defines an alternate set of control characters for plotting indices
!  and exponents. The default characters '[', ']' and '$' are replaced
!  by '^', '_' and '%'.
!
!~        CALL NEWMIX()
!
! Used to enable TeX mode (latex) in DISLIN
!  In TeX mode, all character strings passed to DISLIN routines can
!   contain TeX instructions for plotting mathematical formulas
!
       CALL TEXMOD('ON')
!
! Add axis titles
!
       if ( (ilog.eq.1) .and. (iyunit.eq.1) ) then
         CALL NAME('$\epsilon$ (M$^{-1}$cm$^{-1}$)','Y')
       else if ( (ilog.eq.2) .and. (iyunit.eq.1) ) then
         CALL NAME('log $\epsilon$ (M$^{-1}$cm$^{-1}$)','Y')
       else if ( (ilog.eq.1) .and. (iyunit.eq.2) ) then
         CALL NAME('Absorbance','Y')
       else if ( (ilog.eq.2) .and. (iyunit.eq.2) ) then
         CALL NAME('log [Absorbance])','Y')
       else if ( iyunit .eq. 3 ) then
         CALL NAME('Transmittance (%)','Y')
       end if
!
       if ( ixunit .eq. 1 ) then
         CALL NAME('Wavenumber (cm$^{-1}$)','X')
       else if ( ixunit .eq. 2 ) then
         CALL NAME('Wavelength (nm)','X')
       else if ( ixunit .eq. 3 ) then
         CALL NAME('Energy (eV)','X')
       end if
!
! Changes the value of a text widget
!
       if ( imenu .eq. 1 ) then
         CALL SWGFLT(idxmin,xminaux,-2)
         CALL SWGFLT(idxmax,xmaxaux,-2)
         CALL SWGFLT(idymin,yminaux,-2)
         CALL SWGFLT(idymax,ymaxaux,-2)
       end if
!
! Defines the number of decimal places in bar labels
!
       CALL LABDIG(nxdig,'X')                    ! FLAG: Use with GAXPAR
       CALL LABDIG(nydig,'Y')                    ! FLAG: Use with GAXPAR
!
! Add text lines for titles
!
       CALL TITLIN('Absorption Spectrum',4)
!
! Calculates an explicit colour value
!
       ic = INTRGB(1.0,1.0,1.0)
!
! Defines the background colour
!
       CALL AXSBGD(ic)
!
! Plots a two-dimensional axis system
!
!~        CALL GRAF(xscl(1),xscl(2),xor,xstp,yscl(1),yscl(2),yor,ystp)
       CALL GRAF(xminaux,xmaxaux,xor,xstp,yminaux,ymaxaux,yor,ystp)
!
! Defines the colours used for plotting text and lines
!
       CALL COLOR('FORE')
!
! Plots a title over an axis system
!
       CALL TITLE()
!
! Printing stick spectra if requested
! -----------------------------------
!
! Returns the status of a button or push button widget
!
       CALL GWGBUT(idstick,istick)
!
       if ( istick .ne. 0 ) then
         CALL CURVE(xstick,ystick,nstick)
       end if
!
! Printing convoluted spectra if requested
! ----------------------------------------
!
       CALL GWGBUT(idconv,iconv)
!
       if ( iconv .ne. 0 ) then
!
!~ write(*,*) 'CONVOLUTION ACTIVATED'
!~ write(*,*) icolor
!~          if ( icolor .eq. 1 ) then
!~ write(*,*) 'COLOR ACTIVATED', R, G , B
!~            ic = INTRGB(R,G,B)
!~            CALL SETRGB(R,G,B)
!~            CALL COLOR('FORE')
!~          end if
!
         CALL THKCRV(4)
         CALL CURVE(xline,yline,npts)
!
       end if
!
! Printing reference spectra if requested
! ---------------------------------------
!
       CALL GWGBUT(idhidden,ihidden)
!
       if ( (iref.ne.0) .and. (ihidden.eq.0) ) then
         CALL THKCRV(4)
         CALL COLOR('RED')
         CALL CURVE(xref,yref,nref)
       end if
!
! Terminates DISLIN
! -----------------
!
       CALL DISFIN()
!
       return
       end subroutine iniplot
!
!======================================================================!
!
       subroutine plot(id)
!
       use variables
!
       implicit none
!
! Input/Output variables
!
       integer,intent(in)  ::  id
!
! Printing a beautiful spectrum
!
       call oldmargin()
!
       call iniplot(1)
!
       return
       end subroutine plot
!
!======================================================================!
!
       subroutine convolute()
!
       use variables
       use lineshape
!
       implicit none
!
! Convoluting stick spectrum
!
       if ( ixunit .eq. 1 ) then
!
         if ( iline .eq. 1 ) then
           call gausscm(npts,xline,yline,xmin,dx,hwhm,1)
         else if ( iline .eq. 2 ) then
           call lorcm(npts,xline,yline,xmin,dx,hwhm,1)
         else
           call voigtcm(npts,xline,yline,xmin,dx,gauhwhm,lorhwhm,1)
         end if
!
       else
!
         if ( iline .eq. 1 ) then
           call gaussnm(npts,xline,yline,xmin,dx,hwhm,1)
         else if ( iline .eq. 2 ) then
           call lornm(npts,xline,yline,xmin,dx,hwhm,1)
         else
           call voigtnm(npts,xline,yline,xmin,dx,gauhwhm,lorhwhm,1)
         end if
!
       end if
!
       ylinemax  = maxval(yline)
       ylinemaxl = log10(ylinemax)
!
! Printing spectra with selected modifications
!
       if ( iyunit .eq. 2 ) then
!
         yline(:)  = yline(:)*conc*path
         ylinemax  = maxval(yline)
         ylinemaxl = log10(ylinemax)
!
       else if ( iyunit .eq. 3 ) then
!
         ylinemax  = maxval(yline*cold*lold)
         ylinemaxl = log10(ylinemax)
         yline(:)  = 10.0**(-yline(:)*cold*lold+2)
         return
!
       end if
!
       if ( (ilog.eq.2) .and. (inorm.eq.1) ) then
!
         yline(:)  = log10(yline)/ylinemaxl
!
       else if ( ilog .eq. 2 ) then
!
         yline(:)  = log10(yline)
!
       else if ( inorm .eq. 1 ) then
!
         yline(:)  = yline(:)/ylinemax
!
       end if
!
!~        if ( iyunit .eq. 2 ) then
!~ !
!~          yline(:)  = yline(:)*conc*path
!~          ystick(:) = ystick(:)*conc*path
!~          ylinemax  = maxval(yline)
!~          ylinemaxl = log10(ylinemax)
!~ !
!~        else if ( iyunit .eq. 3 ) then
!~ !
!~          ylinemax  = maxval(yline*cold*lold)
!~          ylinemaxl = log10(ylinemax)
!~          yline(:)  = 10.0**(-yline(:)*cold*lold+2)
!~          ystick(:) = 10.0**(-ystick(:)*cold*lold+2)
!~          return
!~ !
!~        end if
!~ !
!~        if ( (ilog.eq.2) .and. (inorm.eq.1) ) then
!~ !
!~          yline(:)  = log10(yline)/ylinemaxl
!~          ystick(:) = log10(ystick)/ylinemaxl
!~ !
!~        else if ( ilog .eq. 2 ) then
!~ !
!~          yline(:)  = log10(yline)
!~          ystick(:) = log10(ystick)
!~ !
!~        else if ( inorm .eq. 1 ) then
!~ !
!~          yline(:)  = yline(:)/ylinemax
!~          ystick(:) = ystick(:)/ylinemax
!~ !
!~        end if
!
       return
       end subroutine convolute
!
!======================================================================!
!
       subroutine newplot(id)
!
       use variables
       use parameters
       use colour
!
       implicit none
!
       integer,intent(in)  ::  id
!
! Modifying the convoluted spectrum
! ---------------------------------
!
! Requests the value of a box widget (returns the selected element)
!  Reading lineshape function type selected (gau, lor or voigt)
!
       CALL GWGBOX(idline,iline)
!
! Sets widget attributes
!  Enables and disables HWHM parameters for different profiles
!
       if ( iline .eq. 3 ) then
         CALL SWGATT(idhwhm,'INACTIVE','STATUS')
         CALL SWGATT(idgauhwhm,'ACTIVE','STATUS')
         CALL SWGATT(idlorhwhm,'ACTIVE','STATUS')
       else
         CALL SWGATT(idhwhm,'ACTIVE','STATUS')
         CALL SWGATT(idgauhwhm,'INACTIVE','STATUS')
         CALL SWGATT(idlorhwhm,'INACTIVE','STATUS')
       end if
!
! Requests the value of a text widget as real number (returns the input)
!  Reading the Half Width at Half Maximum
!
       CALL GWGFLT(idhwhm,hwhm)
       CALL GWGFLT(idgauhwhm,gauhwhm)
       CALL GWGFLT(idlorhwhm,lorhwhm)
!
       if ( ixunit .eq. 1 ) then
         hwhm    = hwhm*ev2cm
         gauhwhm = gauhwhm*ev2cm
         lorhwhm = lorhwhm*ev2cm
       else
         hwhm    = 1.0E7/(hwhm*ev2cm)
         gauhwhm = 1.0E7/(gauhwhm*ev2cm)
         lorhwhm = 1.0E7/(lorhwhm*ev2cm)
       end if
!
! Recomputing the RGB coordinates of the theoretical spectrum           ! FLAG: TODO
!
       call colorate(illu,iadapt,2,nref,xref,yref,Xwp,Ywp,Zwp,         &
                     X,Y,Z,R,G,B)
!
! Recomputing the spectral lineshape
!
       call convolute()
!
! Printing spectra with the new modifications
!
       call iniplot(1)
!
       return
       end subroutine newplot
!
!======================================================================!
!
       subroutine newmargin(id)
!
       use variables
!
       implicit none
!
       integer,intent(in)  ::  id
!
! Changing the lower and upper limits of the X- and Y-axis
!
       CALL GWGFLT(idxmin,xscl(1))
       CALL GWGFLT(idxmax,xscl(2))
       CALL GWGFLT(idymin,yscl(1))
       CALL GWGFLT(idymax,yscl(2))
!
! OBSOLETE AFTER USING GAXPAR/SETSCL         !  FLAG: TODO prevent GAXPAR from changing selected limits
!                                                     when calling iniplot
!
! Print 5+1 labels in each axis
       xstp = (xscl(2)-xscl(1))/5
       ystp = (yscl(2)-yscl(1))/5
! Start printing labels at the minimum value of each axis
       xor = xscl(1)
       yor = yscl(1)
!
! Printing spectra with the new modifications
!
       call iniplot(1)
!
       return
       end subroutine newmargin
!
!======================================================================!
!
       subroutine resetplot(id)
!
       use variables
       use parameters
       use colour
!
       implicit none
!
       integer,intent(in)  ::  id
!
! Changing scaling factor to default value
!
!~        factor = 1.0
!
!~        CALL SWGFLT(idfactor,factor,4)
!
! Changing HWHM to default values
!
!~        CALL GWGLIS(idxunit,ixunit)
!
       CALL SWGATT(idhwhm,'ACTIVE','STATUS')
       CALL SWGATT(idgauhwhm,'ACTIVE','STATUS')
       CALL SWGATT(idlorhwhm,'ACTIVE','STATUS')
!
       hwhm    = hwhmev
       gauhwhm = hwhmev
       lorhwhm = hwhmev
!
       CALL SWGFLT(idhwhm,hwhm,4)
       CALL SWGFLT(idgauhwhm,gauhwhm,4)
       CALL SWGFLT(idlorhwhm,lorhwhm,4)
!
       if ( iline .eq. 3 ) then
         CALL SWGATT(idhwhm,'INACTIVE','STATUS')
         CALL SWGATT(idgauhwhm,'ACTIVE','STATUS')
         CALL SWGATT(idlorhwhm,'ACTIVE','STATUS')
       else
         CALL SWGATT(idhwhm,'ACTIVE','STATUS')
         CALL SWGATT(idgauhwhm,'INACTIVE','STATUS')
         CALL SWGATT(idlorhwhm,'INACTIVE','STATUS')
       end if
!
       if ( ixunit .eq. 1 ) then
         hwhm    = hwhm*ev2cm
         gauhwhm = gauhwhm*ev2cm
         lorhwhm = lorhwhm*ev2cm
       else
         hwhm    = 1.0E7/(hwhm*ev2cm)
         gauhwhm = 1.0E7/(gauhwhm*ev2cm)
         lorhwhm = 1.0E7/(lorhwhm*ev2cm)
       end if
!
! Recomputing the RGB coordinates of the theoretical spectrum           ! FLAG: TODO
!
       call colorate(illu,iadapt,2,nref,xref,yref,Xwp,Ywp,Zwp,         &
                     X,Y,Z,R,G,B)
!
! Recomputing the spectral lineshape
!
       call convolute()
!
! Changing the lower and upper limits of the axes to default values
!
       call oldmargin()
!
! Printing spectra with the new modifications
!
       call iniplot(1)
!
       return
       end subroutine resetplot
!
!======================================================================!
!
       subroutine importref(id)
!
       use variables
       use utils,     only: print_end,                                 &
                            uniinp
!
       implicit none
!
! Input/Output variables
!
       integer,intent(in)  ::  id
!
! Local variables
!
       character(len=256)  ::  line
       character(len=8)    ::  ext
       character(len=1)    ::  aux
       integer             ::  io
       integer             ::  i
!
! Requesting user a reference file
!
       if ( iref .eq. 0 ) CALL DWGFIL('Please, choose a reference '//  &
                                                'spectrum file',ref,'*')
!
       open(unit=uniinp,file=trim(ref),action='read',                  &
                                                 status='old',iostat=io)
!
       if ( io .ne. 0 ) then
         write(*,'(2X,68("="))')
         write(*,'(3X,A)') 'ERROR:  Missing reference spectrum file'
         write(*,*)
         write(*,'(3X,A)') '  Reference file "'//trim(ref)//'" not'//  &
                                                                ' found'
         write(*,'(2X,68("="))')
         write(*,*)
         call print_end()
       end if
!
       CALL SWGATT(idsepref,'ACTIVE','STATUS')
       CALL SWGATT(idlabref,'ACTIVE','STATUS')
       CALL SWGATT(idnoref,'ACTIVE','STATUS')
       CALL SWGATT(idhidden,'ACTIVE','STATUS')
       CALL SWGATT(idcolref,'ACTIVE','STATUS')
!
! Checking extension of the reference file
!
       ext = ref(len_trim(ref)-3:)
!
       select case ( trim(ext) )
!
         case ('.dat','.csv','.jdx')
!
           iref    = 1
           ihidden = 0
!
         case default
!
           write(*,*)
           write(*,'(2X,68("="))')
           write(*,'(3X,A)')    'ERROR:  Invalid reference spectru'//  &
                                                      'm file extension'
           write(*,*)
           write(*,'(3X,A)') 'Unrecognised extension : '//trim(ext)
           write(*,*)
           write(*,'(3X,A)') 'Please, choose between : ".dat", ".c'//  &
                                                         'sv" or ".jdx"'
           write(*,'(2X,68("="))')
           write(*,*)
           call print_end()
!
       end select
!
! Counting number of lines in the file
!
       nref = 0
!
       select case ( trim(ext) )
!
         case ('.dat','.csv')
!
           do
             read(uniinp,*,iostat=io)
             if ( io .ne. 0 ) exit
             nref = nref + 1
           end do
!
         case ('.jdx')
!
           do
             read(uniinp,*,iostat=io) line
             if ( io .ne. 0 ) exit
             line = adjustl(line)
             if ( line(1:1) .ne. '#' ) nref = nref + 1
           end do
!
       end select
!
! Allocating memory for the reference spectrum
!
       allocate(xref(nref),yref(nref))
!
! Reading reference spectrum
!
       rewind(uniinp)
!
       select case ( trim(ext) )
!
         case ('.dat')
!
           do i = 1, nref
!
             read(uniinp,'(A)') line
             read(line,*,iostat=io) xref(i),yref(i)
!
             if ( io .ne. 0 ) call print_errimp(line)
!
           end do
!
         case ('.csv')
!
           do i = 1, nref
!
             read(uniinp,'(A)') line
             read(line,*,iostat=io) xref(i),yref(i)
!
             if ( io .ne. 0 ) call print_errimp(line)
!
           end do
!
         case ('.jdx')
!
           i = 1
           do
!
             read(uniinp,'(A)',iostat=io) line
             if ( io .ne. 0 ) exit
!
             line = adjustl(line)
!
             if ( line(1:1) .ne. '#' ) then
               read(line,*,iostat=io) xref(i),yref(i)
               if ( io .ne. 0 ) call print_errimp(line)
               i = i + 1
             end if
!
           end do
!
       end select
!
       yrefmax  = maxval(yref)
       yrefmaxl = log10(yrefmax)
!
       if ( (ilog.eq.2) .and. (inorm.eq.1) ) then
!
!~          yref(:) = log10(yref)/yrefmaxl
         yscl(2) = 1.1
!
       else if ( ilog .eq. 2 ) then
!
!~          yref(:) = log10(yref)
         yscl(2) = max(ylinemaxl,yrefmaxl)
!
       else if ( inorm .eq. 1 ) then
!
!~          yref(:) = yref(:)/yrefmax
         yscl(2) = max(ylinemax,yrefmax)
!
       end if
!
! Closing reference spectrum file
!
       close(uniinp)
!
! Printing spectra with the new modifications
!
       if ( id .ne. idimport ) return
!
       call iniplot(1)
!
       return
       end subroutine importref
!
!======================================================================!
!
       subroutine print_errimp(line)
!
       use variables
       use utils,     only: print_end
!
       implicit none
!
! Input/Output variables
!
       character(len=256)  ::  line
!
       write(*,'(2X,68("="))')
       write(*,'(3X,A)') 'ERROR:  Error while reading reference sp'//  &
                                                           'ectrum file'
       write(*,*)
       write(*,'(3X,A)') 'Unrecognized line  :'
       write(*,'(5X,A)') trim(adjustl(line))
       write(*,*)
       write(*,'(3X,A)') 'Please, check file "'//trim(ref)//'"'
       write(*,'(2X,68("="))')
       write(*,*)
       call print_end()
!
       return
       end subroutine print_errimp
!
!======================================================================!
!
       subroutine axismenu(id)
!
       use variables
!
       implicit none
!
! Input/Output variables
!
       integer,intent(in)  ::  id
!
! Local variables
!
       integer             ::  ipaux
!
! Creating widget for graph manipulation options
! ----------------------------------------------
!
       if ( imenu .eq. 1 ) return
!
       imenu = 1
!
! Modifies the appearance of the popup menu bar
!  'NOQUIT' suppresses the 'QUIT' entry in the 'EXIT' menu
!  With this user cannot terminate a program by accident
!
       CALL SWGPOP('NOQUIT')
!
       CALL SWGWTH(-15)
!
       CALL SWGTIT('Axis manipulation options')
!
       CALL SWGHLP('With this menu you can change the limits of th'//  &
                                                      'e X- and Y-axis')
!
       CALL WGINI('VERT',idaxsmenu)
!
! Creating graph manipulation options
! ...................................
!
       CALL WGLAB(idaxsmenu,'X-axis',idlab)
!
       CALL WGBAS(idaxsmenu,'HORI',ipaux)
!
       CALL SWGWTH(-10)
       CALL WGLTXT(ipaux,'Start:','0.00',60,idxmin)
       CALL WGLTXT(ipaux,'Stop:','0.00',60,idxmax)
!
       CALL WGLAB(idaxsmenu,'Y-axis',idlab)
!
       CALL WGBAS(idaxsmenu,'HORI',ipaux)
!
       CALL SWGWTH(-10)
       CALL WGLTXT(ipaux,'Start:','0.00',60,idymin)
       CALL WGLTXT(ipaux,'Stop:','0.00',60,idymax)
!
! Changes the value of a text widget
!
       CALL REAWGT()
!
       CALL SWGFLT(idxmin,xminaux,-2)
       CALL SWGFLT(idxmax,xmaxaux,-2)
       CALL SWGFLT(idymin,yminaux,-2)
       CALL SWGFLT(idymax,ymaxaux,-2)
!
       CALL SWGCBK(idxmin,newmargin)
       CALL SWGCBK(idxmax,newmargin)
       CALL SWGCBK(idymin,newmargin)
       CALL SWGCBK(idymax,newmargin)
!
       CALL SWGWTH(-15)
!
       CALL WGFIN()
!
       imenu = 0
!
       return
       end subroutine axismenu
!
!======================================================================!
!
       subroutine unitsmenu(id)
!
       use variables
!
       implicit none
!
! Input/Output variables
!
       integer,intent(in)  ::  id
!
! Local variables
!
       integer             ::  ipaux
       integer             ::  ipmenu
!
! Creating widget for axis options
! --------------------------------
!
       if ( iunits .eq. 1 ) return
!
       iunits = 1
!
! Modifies the appearance of the popup menu bar
!  'NOQUIT' suppresses the 'QUIT' entry in the 'EXIT' menu
!
       CALL SWGPOP('NOQUIT')
!
       CALL SWGWTH(-22)
!
       CALL SWGTIT('Units/scale options')
!
       CALL SWGHLP('With this menu you can change the units and/or'//  &
                                               ' scale of the spectrum')
!
       CALL WGINI('VERT',ipmenu)
!
! Creating spectrum units options
! ...............................
!
       CALL WGBAS(ipmenu,'HORI',ipaux)
!
       CALL SWGWTH(-5)
       CALL WGLAB(ipaux,'X Units',idlab)
!
! Creates a dropping list widget.
!  This list widget can be used to save space in the parent widget.
!  This widget is used whenever an application must present a list of
!   names from which the user can choose.
!
       CALL SWGWTH(-15)
       CALL WGDLIS(ipaux,'Wavenumber|Wavelength',ixunit,idxunit)
!
       CALL WGBAS(ipmenu,'HORI',ipaux)
!
       CALL SWGWTH(-5)
       CALL WGLAB(ipaux,'Y Units',idlab)
       CALL SWGWTH(-15)
       CALL WGDLIS(ipaux,'Absorption coefficient|Absorbance|Transm'//  &
                   'ittance',iyunit,idyunit)
!
! Creating scaling options
! ........................
!
       CALL WGSEP (ipmenu,idsep)
!
       CALL WGBAS(ipmenu,'HORI',ipaux)
!
       CALL SWGWTH(-5)
       CALL WGLAB(ipaux,'Y Scale',idlab)
       CALL SWGWTH(-15)
       CALL WGDLIS(ipaux,'Linear|Logarithmic',ilog,idlog)
!
! Creating experimental conditions options
! ........................................
!
       CALL WGSEP (ipmenu,idsep)
       CALL WGLTXT(ipmenu,'Concentration (M)','1.00',40,idconc)
!
       if ( iyunit .eq. 1 ) then
         CALL SWGFLT(idconc,1.0,2)
       else if ( iyunit .eq. 2 ) then
         CALL SWGFLT(idconc,conc,-2)
       end if
!
       CALL WGLTXT(ipmenu,'Path-length (cm)','1.00',40,idpath)
!
! Updating widgets according to the units
! .......................................
!
       CALL REAWGT()
!
       if ( iyunit .eq. 1 ) then
         CALL SWGFLT(idconc,1.0,2)
         CALL SWGFLT(idpath,1.0,2)
         CALL SWGATT(idconc,'INACTIVE','STATUS')
         CALL SWGATT(idpath,'INACTIVE','STATUS')
       else
         CALL SWGFLT(idconc,conc,-2)
         CALL SWGFLT(idpath,path,-2)
         CALL SWGATT(idlog,'ACTIVE','STATUS')
         CALL SWGATT(idpath,'ACTIVE','STATUS')
       end if
!
! Changing units/scale when new list item is selected
! ...................................................
!
       CALL SWGCBK(idlog,newlog)
!
       CALL SWGCBK(idxunit,newxunits)
       CALL SWGCBK(idyunit,newyunits)
!
       CALL SWGCBK(idconc,newabs)
       CALL SWGCBK(idpath,newabs)
!
       CALL SWGWTH(-15)
!
       CALL WGFIN()
!
       iunits = 0
!
       return
       end subroutine unitsmenu
!
!======================================================================!
!
       subroutine colormenu(id)
!
       use variables
!
       implicit none
!
! Input/Output variables
!
       integer,intent(in)  ::  id
!
! Local variables
!
       integer             ::  ipmenu
       integer             ::  ipaux
!
! Creating widget for colour options
! ----------------------------------
!
       if ( icolmenu .eq. 1 ) return
!
       icolmenu = 1      
!
       CALL SWGPOP('NOQUIT')
!
       CALL SWGWTH(-20)
!
       CALL SWGTIT('Colour calculation help')
!
       CALL SWGHLP('With this menu you can change the method to co'//  &
                                           'mpute the predicted colour')
!
       CALL SWGTIT('Colour calculation options')
!
       CALL WGINI('VERT',ipmenu)
!
       CALL WGBAS(ipmenu,'HORI',ipaux)
!
! Creating colour method options
! ..............................
!
       CALL SWGSPC(-1.0,-0.0)
       CALL WGBAS(ipmenu,'HORI',ipaux)
!
       CALL SWGWTH(-10)
       CALL WGLAB(ipaux,'Illuminant',idlab)
!
       CALL SWGWTH(-16)
       CALL WGDLIS(ipaux,'A (Tungsten filament)|C (Average dayligh'//  &
                         't)|D50 (Horizon daylight)|D55 (Vertical '//  &
                         'daylight)|D65 (Daylight, overcast)|D75 ('//  &
                         'Phase of daylight)|E (Equal-energy Radia'//  &
                         'tor)|Inverse method|Black Body',iillu,idillu)
!
       CALL SWGSPC(4.0,0.5)
       CALL WGBAS(ipmenu,'VERT',ipaux)
       CALL SWGSPC(-1.0,-0.0)
       CALL WGBAS(ipmenu,'HORI',ipaux)
!
       CALL SWGWTH(-10)
       CALL WGLAB(ipaux,'RGB Model',idlab)
!
!
       CALL SWGWTH(-16)
       CALL WGDLIS(ipaux,'Adobe RGB|Apple RGB|Best RGB|Beta RGB|Br'//  &
                         'uce RGB|CIE RGB|Color Match RGB|Don RGB '//  &
                         ' 4|ECI RGB v2|Ekta space PS5|HDTV RGB|NT'//  &
                         'SC RGB|EBU (PAL/SECAM) RGB|ProPhoto RGB|'//  &
                         'SMPTE RGB|sRGB|Wide Gamut RGB',imodel,idmodel)
!
       CALL SWGSPC(4.0,0.5)
       CALL WGBAS(ipmenu,'VERT',ipaux)
       CALL SWGSPC(-1.0,-0.0)
       CALL WGBAS(ipmenu,'HORI',ipaux)
!
       CALL SWGWTH(-10)
       CALL WGLAB(ipaux,'Adaptation method',idlab)
!
       CALL SWGWTH(-16)
       CALL WGDLIS(ipaux,'Bradford|von Kries|XYZ Scaling|None',        &
                   iadapt,idadapt)
!
       CALL SWGSPC(4.0,0.5)
       CALL SWGWTH(-27)
       CALL WGBAS(ipmenu,'VERT',ipaux)
!
       CALL WGLTXT(ipaux,'Black body T (K)','0.00',59,idcct)
!
       CALL REAWGT()
       CALL SWGFLT(idcct,cct,2)
!
! Experimental spectrum does not admidt standard illuminant method
! ................................................................
!



!
! Updating widgets according to the units
! .......................................
!
!  Enables correlated color temperature option if the illuminant is 
!   a Black Body (iillu = 9)
!
       if ( iillu .ne. 9 ) then      
         CALL SWGATT(idcct,'INACTIVE','STATUS')
       else
         CALL SWGATT(idcct,'ACTIVE','STATUS')
       end if
!
! Changing color when a new option is selected
! ............................................
!
!~        CALL SWGCBK(idillu,newcolor)
!~        CALL SWGCBK(idmodel,newcolor)
!~        CALL SWGCBK(idadapt,newcolor)
!~        CALL SWGCBK(idcct,newcolor)
!
       CALL WGFIN()
!
       icolmenu = 0
!
       CALL SWGSPC(4.0,0.5)
!
       return
       end subroutine colormenu
!
!======================================================================!
!
       subroutine newxunits(id)
!
       use variables
!
       implicit none
!
! Input/Output variables

       integer,intent(in)  ::  id
!
! Local variables
!
       integer             ::  iaux
!
! Modifying units of the X-axis
! -----------------------------
!
       iaux = ixunit
!
       CALL GWGLIS(idxunit,ixunit)
!
       if ( ixunit .eq. iaux ) return
!
       deallocate(xstick,ystick)
       deallocate(xline,yline)
!
       freq(:nband) = 1.0E7/freq(:nband)
!
       hwhm    = 1.0E7/hwhm
       gauhwhm = 1.0E7/gauhwhm
       lorhwhm = 1.0E7/lorhwhm
!
       if ( ixunit .eq. 1 ) then       ! FLAG: TODO allow to reuse values input by user
         dx   = dxcm
         xmin = xmincm
         xmax = xmaxcm
       else
         dx   = dxnm
         xmin = xminnm
         xmax = xmaxnm
       end if
!
! Changes the value of a text widget and the text string of label and
!  push button widgets.
!
! https://groups.google.com/g/dislin-users/c/wpvKp8r1gtY/m/CsZZyEo6BQAJ
!
!   "swgtxt changes the contents of text widgets and also the text of
!    label widgets. If you are using the widget routine wgltxt, which
!    contains a label and a text field, swgtxt changes the text field
!    if you pass the widget id returned by wgltxt to swgtxt. If you pass
!    the value  id - 1 to swgtxt, the label field will be changed"
!
!~        if ( ixunit .eq. 1 ) then
!~          CALL SWGTXT(idhwhm-1,'HWHM (cm**-1)')
!~        else if ( ixunit .eq. 2 ) then
!~          CALL SWGTXT(idhwhm-1,'HWHM (nm)')
!~        else if ( ixunit .eq. 3 ) then
!~          CALL SWGTXT(idhwhm-1,'HWHM (eV)')
!~        end if
!
! Recomputing the spectral lineshape
!
       call initialize()
!
! Printing spectra with the new modifications
!
       call iniplot(1)
!
       return
       end subroutine newxunits
!
!======================================================================!
!
       subroutine newyunits(id)
!
       use variables
       use lineshape
!
       implicit none
!
! Input/Output variables

       integer,intent(in)  ::  id
!
! Local variables
!
       integer             ::  iaux
!
! Modifying units of the Y-axis
! -----------------------------
!
       iaux = iyunit
!
       CALL GWGLIS(idyunit,iyunit)
!
       if ( iyunit .eq. iaux ) return
!
! Transforming the actual scaling to non-normalized linear scaling
! ................................................................
!
!  Transmittance units do not admit normalization/logarithmic scale
!~ !
!~        if ( iaux .ne. 3 ) then
!~          if ( (ilog.eq.2) .and. (inorm.eq.1) ) then
!~ !
!~ ! Renormalize and go back to the linear scaling
!~ !
!~            yline(:)  = 10.0**(yline(:)*ylinemaxl)
!~            ystick(:) = 10.0**(ystick(:)*ylinemaxl)
!~            if ( iref .eq. 1 ) yref(:) = 10.0**(yref(:)*yrefmaxl)
!~ !
!~          else if ( ilog .eq. 2 ) then
!~ !
!~ ! Go to the linear scaling
!~ !
!~            yline(:)  = 10.0**yline(:)
!~            ystick(:) = 10.0**ystick(:)
!~            if ( iref .eq. 1 ) yref(:) = 10.0**yref(:)
!~ !
!~          else if ( inorm .eq. 1 ) then
!~ !
!~ ! Renormalize from linear scaling
!~ !
!~            yline(:)  = yline(:)*ylinemax
!~            ystick(:) = ystick(:)*ylinemax
!~            if ( iref .eq. 1 ) yref(:) = yref(:)*yrefmax
!~ !
!~          end if
!~        end if
!~ !
!~ ! Changing units consistently
!~ ! ...........................
!~ !
!~        if ( (iaux.eq.1) .and. (iyunit.eq.2) ) then
!~ !
!~ !  Going from molar absorption coefficient to absorbance
!~ !
!~          yline(:)  = yline(:)*conc*path
!~          ystick(:) = ystick(:)*conc*path
!~ ! Changing concentration and path-length labels to selected values
!~          CALL SWGFLT(idconc,conc,-2)
!~          CALL SWGFLT(idpath,path,-2)
!~ !
!~        else if ( (iaux.eq.2) .and. (iyunit.eq.1) ) then
!~ !
!~ !  Going from absorbance to molar absorption coefficient
!~ !
!~          yline(:)  = yline(:)/conc/path
!~          ystick(:) = ystick(:)/conc/path
!~ ! Changing concentration and path-length labels to 1 M and 1 cm
!~          CALL SWGFLT(idconc,1.0,2)
!~          CALL SWGFLT(idpath,1.0,2)
!~ !
!~        else if ( (iaux.eq.2) .and. (iyunit.eq.3) ) then
!~ !
!~ ! Going from absorbance to transmittance
!~ !
!~          yline(:)  = 10.0**(-yline(:))
!~          ystick(:) = 10.0**(-ystick(:))
!~          if ( iref .eq. 1 ) yref(:) = 10.0**(-yref(:))
!~ !
!~          CALL SWGFLT(idconc,conc,-2)
!~          CALL SWGFLT(idpath,path,-2)
!~ !
!~          cold = conc
!~          lold = path
!~ !
!~        else if ( (iaux.eq.1) .and. (iyunit.eq.3) ) then
!~ !
!~ ! Going from molar absorption coefficient to transmittance
!~ !
!~          yline(:)  = 10.0**(-yline(:))
!~          ystick(:) = 10.0**(-ystick(:))
!~          if ( iref .eq. 1 ) yref(:) = 10.0**(-yref(:))
!~ !
!~          cold = 1.0d0
!~          lold = 1.0d0
!~ !
!~          CALL SWGFLT(idconc,1.0,2)
!~          CALL SWGFLT(idpath,1.0,2)
!~ !
!~        else if ( (iaux.eq.3) .and. (iyunit.eq.2) ) then
!~ !
!~ ! Going from transmittance to absorbance
!~ !
!~          yline(:)  = -log10(yline(:))!/cold/lold
!~          ystick(:) = -log10(ystick(:))!/cold/lold
!~          if ( iref .eq. 1 ) yref(:) = -log10(yref(:))
!~ !
!~          CALL SWGFLT(idconc,conc,-2)
!~          CALL SWGFLT(idpath,path,-2)
!~ !
!~        else if ( (iaux.eq.3) .and. (iyunit.eq.1) ) then
!~ !
!~ ! Going from transmittance to molar absorption coefficient
!~ !
!~ WRITE(*,*) 'BEFORE',MAXVAL(YLINE),MINVAL(YLINE)
!~          yline(:)  = -log10(yline(:)**(cold*lold))
!~          ystick(:) = -log10(ystick(:)**(cold*lold))
!~ WRITE(*,*) 'AFTER',MAXVAL(YLINE),MINVAL(YLINE)
!~          if ( iref .eq. 1 ) yref(:) = -log10(yref(:))
!~ !
!~ !         call convolute()
!~ !
!~ !         if ( ixunit .eq. 1 ) then
!~ !          call deltacm()
!~ !         else
!~ !          call deltanm()
!~ !         end if
!~ !
!~          CALL SWGFLT(idconc,1.0,2)
!~          CALL SWGFLT(idpath,1.0,2)
!~ !
!~        end if
!~ !
!~ ! Transforming back to the previous scaling
!~ ! .........................................
!~ !
!~        yscl(2) = 1.1
!~ !
!~        if ( iyunit .ne. 3 ) then
!~ !
!~          ylinemax  = maxval(yline(:))
!~          ylinemaxl = log10(ylinemax)
!~ !
!~          if ( (ilog.eq.2) .and. (inorm.eq.1) ) then
!~ !
!~ ! Go back to the logarithmic scaling and normalize
!~ !
!~            yline(:)  = log10(yline(:))/ylinemaxl
!~            ystick(:) = log10(ystick(:))/ylinemaxl
!~            if ( iref .eq. 1 ) yref(:) = log10(yref(:))/yrefmaxl
!~ !
!~          else if ( ilog .eq. 2 ) then
!~ !
!~ ! Go back to the logarithmic scaling
!~ !
!~            yline(:)  = log10(yline(:))
!~            ystick(:) = log10(ystick(:))
!~            if ( iref .eq. 1 ) yref = log10(yref(:))
!~ !
!~            yscl(2) = max(ylinemaxl,yrefmaxl)
!~ !
!~          else if ( inorm .eq. 1 ) then
!~ !
!~ ! Normalize from linear scaling
!~ !
!~            yline(:)  = yline(:)/ylinemax
!~            ystick(:) = ystick(:)/ylinemax
!~            if ( iref .eq. 1 ) yref(:) = yref(:)/yrefmax
!~ !
!~          else
!~ !
!~            yscl(2) = max(ylinemax,yrefmax)
!~ !
!~          end if
!~ !
!~        end if
!
! Updating widgets according to the new units
! ...........................................
!
       CALL SWGATT(idconc,'ACTIVE','STATUS')
       CALL SWGATT(idpath,'ACTIVE','STATUS')
!
       if ( (iaux.eq.1) .and. (iyunit.eq.2) ) then
!
!  Going from molar absorption coefficient to absorbance
!
         CALL SWGFLT(idconc,conc,-2)
         CALL SWGFLT(idpath,path,-2)
!
       else if ( (iaux.eq.2) .and. (iyunit.eq.1) ) then
!
!  Going from absorbance to molar absorption coefficient
!
         CALL SWGFLT(idconc,1.0,2)
         CALL SWGFLT(idpath,1.0,2)
!
       else if ( (iaux.eq.2) .and. (iyunit.eq.3) ) then
!
! Going from absorbance to transmittance
!
         CALL SWGFLT(idconc,conc,-2)
         CALL SWGFLT(idpath,path,-2)
!
         cold = conc
         lold = path
!
       else if ( (iaux.eq.1) .and. (iyunit.eq.3) ) then
!
! Going from molar absorption coefficient to transmittance
!
         CALL SWGFLT(idconc,1.0,2)
         CALL SWGFLT(idpath,1.0,2)
!
         cold = 1.0d0
         lold = 1.0d0
!
       else if ( (iaux.eq.3) .and. (iyunit.eq.1) ) then
!
! Going from transmittance to molar absorption coefficient
!
         CALL SWGFLT(idconc,1.0,2)
         CALL SWGFLT(idpath,1.0,2)
!
       else if ( (iaux.eq.3) .and. (iyunit.eq.2) ) then
!
! Going from transmittance to absorbance
!
         CALL SWGFLT(idconc,conc,-2)
         CALL SWGFLT(idpath,path,-2)
!
         cold = conc
         lold = path
!
       end if
!
       if ( iyunit .ne. 3 ) then               ! FLAG :  USE SWG after activate to set user value
         CALL SWGATT(idnorm,'ACTIVE','STATUS')
         CALL SWGATT(idlog,'ACTIVE','STATUS')
!~          CALL SWGBUT(idnorm,0)
!~          CALL SWGLIS(idxunit,3)
       else                                    ! FLAG : USE SWG before deactivate to set user value
!~          CALL SWGBUT(idnorm,inorm)
!~          CALL SWGLIS(idxunit,ixunit)
         CALL SWGATT(idnorm,'INACTIVE','STATUS')
         CALL SWGATT(idlog,'INACTIVE','STATUS')
       end if
!
       if ( iyunit .eq. 1 ) then    
         CALL SWGATT(idconc,'INACTIVE','STATUS')
         CALL SWGATT(idpath,'INACTIVE','STATUS')
       else
         CALL SWGATT(idconc,'ACTIVE','STATUS')
         CALL SWGATT(idpath,'ACTIVE','STATUS')
       end if
!
! Computing spectra with the new units
! ....................................
!
       if ( ixunit .eq. 1 ) then
         call deltacm()
       else
         call deltanm()
       end if
!
       if ( iyunit .eq. 2 ) then
         ystick(:) = ystick(:)*conc*path
       else if ( iyunit .eq. 3 ) then
         ystick(:) = 10.0**(-ystick(:)*cold*lold+2)
       end if
!
       if ( (ilog.eq.2) .and. (inorm.eq.1) ) then
         ystick(:) = log10(ystick)/ylinemaxl
       else if ( ilog .eq. 2 ) then
         ystick(:) = log10(ystick)
       else if ( inorm .eq. 1 ) then
         ystick(:) = ystick(:)/ylinemax
       end if
!
       call convolute()
!
! Selecting axis limits depending on the units
! ............................................
!
       if ( iyunit .eq. 3 ) then
         yscl(2) = 110.0
       else if ( inorm .eq. 1 ) then
         yscl(2) = 1.1
       else if ( ilog .eq. 2 ) then
         yscl(2) = max(ylinemaxl,yrefmaxl)
       else
         yscl(2) = max(ylinemax,yrefmax)
       end if
!
! Printing spectra with the new modifications
!
       call iniplot(1)
!
       return
       end subroutine newyunits
!
!======================================================================!
!
       subroutine newlog(id)
!
       use variables
!
       implicit none
!
! Input/Output variables

       integer,intent(in)  ::  id
!
! Local variables
!
       integer             ::  iaux
!
! Modifying units of the X-axis
! -----------------------------
!
       iaux = ilog
!
       CALL GWGLIS(idlog,ilog)
       CALL GWGBUT(idnorm,inorm)
!
       if ( ilog .eq. iaux ) return
!
       if ( ilog .eq. 1 ) then                !  Changing to linear scaling
         if ( inorm .eq. 1 ) then             !  Keeping function normalized
!
! Renormalize from logarithmic scaling and change to linear one
!
           yline(:)  = 10.0**(yline(:)*ylinemaxl)
           ystick(:) = 10.0**(ystick(:)*ylinemaxl)
           if ( iref .eq. 1 ) yref(:) = 10**(yref(:)*yrefmaxl)
!
! Normalize in the linear scale
!
           yline(:)  = yline(:)/ylinemax
           ystick(:) = ystick(:)/ylinemax
           if ( iref .eq. 1 ) yref(:) = yref(:)/yrefmax
!
           yscl(2) = 1.1
!
         else
!
! Change from logarithmic scaling to linear
!
           yline     = 10.0**yline
           ystick    = 10.0**ystick
           if ( iref .eq. 1 ) yref = 10**yref
!
           yscl(2) = max(ylinemax,yrefmax)
!
         end if
       else                                   !  Changing to logarithmic scaling
         if ( inorm .eq. 1 ) then             !  Keeping function normalized
!
! Renormalize from linear scaling and change to logarithmic one
!
           yline(:)  = log10(yline(:)*ylinemax)
           ystick(:) = log10(ystick(:)*ylinemax)
           if ( iref .eq. 1 ) yref(:) = log10(yref(:)*yrefmax)
!
! Normalize in the logartihmic scale
!
           yline(:)  = yline(:)/ylinemaxl
           ystick(:) = ystick(:)/ylinemaxl
           if ( iref .eq. 1 ) yref(:) = yref(:)/yrefmaxl
!
           yscl(2) = 1.1
!
         else
!
! Change from linear scaling to logarithmic
!
           yline     = log10(yline)
           ystick    = log10(ystick)
           if ( iref .eq. 1 ) yref = log10(yref)
!
           yscl(2) = max(ylinemaxl,yrefmaxl)
!
         end if
       end if
!
! Printing spectra with the new modifications
!
       call iniplot(1)
!
       return
       end subroutine newlog
!
!======================================================================!
!
       subroutine newabs(id)
!
       use lineshape
       use colour
!
       implicit none
!
! Input/Output variables

       integer,intent(in)  ::  id
!
! Modifying experimental settings
! -------------------------------
!
       CALL GWGFLT(idconc,conc)
       CALL GWGFLT(idpath,path)
!
       cold = conc
       lold = path
!
! Computing spectra with the new settings
! .......................................
!
       if ( ixunit .eq. 1 ) then
         call deltacm()
       else
         call deltanm()
       end if
!
       if ( iyunit .eq. 2 ) then
         ystick(:) = ystick(:)*conc*path
       else if ( iyunit .eq. 3 ) then
         ystick(:) = 10.0**(-ystick(:)*cold*lold+2)
       end if
!
       if ( (ilog.eq.2) .and. (inorm.eq.1) ) then
         ystick(:) = log10(ystick)/ylinemaxl
       else if ( ilog .eq. 2 ) then
         ystick(:) = log10(ystick)
       else if ( inorm .eq. 1 ) then
         ystick(:) = ystick(:)/ylinemax
       end if
!
       call convolute() 
!
       if ( iyunit .eq. 3 ) then     ! FLAG: copied from OLDMARGIN
!
         yscl(2) = 110.0
!
       else if ( inorm .eq. 0 ) then
!
         if ( ilog .eq. 1 ) then
           yscl(2) = max(ylinemax,yrefmax)
         else
           yscl(2) = max(ylinemaxl,yrefmaxl)
         end if
!
       else
!
         yscl(2) = 1.1
!
       end if
!
! Recomputing the RGB coordinates of the theoretical spectrum           ! FLAG: TODO
!
       call colorate(illu,iadapt,2,nref,xref,yref,Xwp,Ywp,Zwp,         &
                     X,Y,Z,R,G,B)
!
! Printing spectra with the new modifications
!
       call iniplot(1)
!
       return
       end subroutine newabs
!
!======================================================================!
!
       subroutine removeref(id)
!
       use variables
!
       implicit none
!
! Input/Output variables
!
       integer,intent(in)  ::  id
!
! Removing experimental spectra
! -----------------------------
!
       iref    = 0
       ihidden = 0
!
       CALL SWGBUT(idhidden,ihidden)
!
       yrefmax = 0.0
!
! Removing experimental spectra options
!
       CALL SWGATT(idsepref,'INVISIBLE','STATUS')
       CALL SWGATT(idlabref,'INVISIBLE','STATUS')
       CALL SWGATT(idnoref,'INVISIBLE','STATUS')
       CALL SWGATT(idhidden,'INVISIBLE','STATUS')
       CALL SWGATT(idcolref,'INVISIBLE','STATUS')
!
! Deallocating memory
!
       deallocate(xref,yref)
!
! Plotting graph without experimental spectra
!
       call iniplot(1)
!
       return
       end subroutine removeref
!
!======================================================================!
!
       subroutine oldmargin()
!
       use variables
       use parameters
!
       implicit none
!
! Local variables
!
       real(kind=4)  ::  w0  !
       integer       ::  i   !
!
! Changing the lower and upper limits of the axes to default values
!
       if ( ixunit .eq. 1 ) then
!
         i       = minloc(freq(:nband),dim=1)
         w0      = freq(i)
         xtail   = fscl*real(Aabs)*inten(i)!*w0
         xtail   = hwhm**2/(2.0*w0*xtail) - hwhm/(2.0*xtail)           &
                            *sqrt((hwhm/w0)**2  - 4.0*xtail*(xtail - 1))
         xscl(1) = max(w0+xtail,xsclmincm)
!
         i       = maxloc(freq(:nband),dim=1)
         w0      = freq(i)
         xtail   = fscl*real(Aabs)*inten(i)!*w0
         xtail   = hwhm**2/(2.0*w0*xtail) + hwhm/(2.0*xtail)           &
                            *sqrt((hwhm/w0)**2  - 4.0*xtail*(xtail - 1))
         xscl(2) = min(w0+xtail,xsclmaxcm)
!
       else
!
         i       = minloc(freq(:nband),dim=1)
         w0      = freq(i)
         xtail   = fscl*real(Aabs)*inten(i)!*(1.0E7/w0)
         xtail   = w0*1.0E7/(2.0*hwhm**2*xtail)                        &
                          - 1.0E7/(2.0*xtail*hwhm)                     &
                            *sqrt((w0/hwhm)**2  - 4.0*xtail*(xtail - 1))
         xscl(1) = max(w0+xtail,xsclminnm)
!
         i       = maxloc(freq(:nband),dim=1)
         w0      = freq(i)
         xtail   = fscl*real(Aabs)*inten(i)!*(1.0E7/w0)
         xtail   = w0*1.0E7/(2.0*hwhm**2*xtail)                        &
                          + 1.0E7/(2.0*xtail*hwhm)                     &
                            *sqrt((w0/hwhm)**2  - 4.0*xtail*(xtail - 1))
         xscl(2) = min(w0+xtail,xsclmaxnm)
!
       end if
!
       yscl(1) = 0.0
!
       if ( iyunit .eq. 3 ) then
!
         yscl(2) = 110.0
!
       else if ( inorm .eq. 0 ) then
!
         if ( ilog .eq. 1 ) then
           yscl(2) = max(ylinemax,yrefmax)
         else
           yscl(2) = max(ylinemaxl,yrefmaxl)
         end if
!
       else
!
         yscl(2) = 1.1
!
       end if
!
       return
       end subroutine oldmargin
!
!======================================================================!
!
       subroutine normalize(id)
!
       use variables
!
       implicit none
!
! Input/Output variables
!
       integer,intent(in)  ::  id
!
! Normalizing or renormalizing spectrums
!
       CALL GWGBUT(idnorm,inorm)
!
       if ( inorm .eq. 1 ) then
!
         if ( ilog .eq. 1 ) then
           yline(:)  = yline(:)/ylinemax
           ystick(:) = ystick(:)/ylinemax
           if ( iref .eq. 1 ) yref(:) = yref(:)/yrefmax
         else
           yline(:)  = yline(:)/ylinemaxl
           ystick(:) = ystick(:)/ylinemaxl
           if ( iref .eq. 1 ) yref = yref/yrefmaxl
         end if
!
         yscl(2) = 1.1
!
       else
!
         if ( ilog .eq. 1 ) then
!
           yline(:)  = yline(:)*ylinemax
           ystick(:) = ystick(:)*ylinemax
           if ( iref .eq. 1 ) yref(:) = yref(:)*yrefmax
!
           yscl(2) = max(ylinemax,yrefmax)
!
         else
!
           yline(:)  = yline(:)*ylinemaxl
           ystick(:) = ystick(:)*ylinemaxl
           if ( iref .eq. 1 ) yref(:) = yref(:)*yrefmaxl
!
           yscl(2) = max(ylinemaxl,yrefmaxl)
!
         end if
!
       end if
!
! Ploting graph with the new scale
!
       call iniplot(1)
!
       return
       end subroutine normalize
!
!======================================================================!
!
       subroutine plotcolor(id)
!
       use variables
       use parameters
       use colour
!
       implicit none
!
! Input/Output variables
!
       integer,intent(in)  ::  id
!
! 
!
       if ( id .eq. idcolor ) then
         CALL GWGBUT(idcolor,icolor)
       else if ( id .eq. idcolref ) then
         CALL GWGBUT(idcolref,icolref)
       end if
!
       call iniplot(0)
!
       return
       end subroutine plotcolor
!
!======================================================================!
!
       end module print_plot
!
!======================================================================!
