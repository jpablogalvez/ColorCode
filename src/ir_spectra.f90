!======================================================================!
!
       program ir_spectra
!
       use variables
       use lineshape
       use colour
       use print_plot
       use utils,      only: print_start,                              &
                             print_end,                                &
                             leninp,                                   &
                             lenout,                                   &
                             uniout
!
       implicit none
!
       character(len=leninp)  ::  inp     !  Input file name
       character(len=lenout)  ::  outp    !  Output file name
       integer                ::  io      !  Status
       integer                ::  i       !  Index
!
! Printing header
!
       write(*,*)
!~        write(*,'(5X,82("#"))')
!~        write(*,'(5X,10("#"),3X,A,3X,13("#"))') 'CompuTherm - Comput'//  &
!~                                     'ational Thermochemistry Calculator'
!~        write(*,'(5X,82("#"))')
!~        write(*,*)
!~        write(*,'(9X,A)') 'Welcome to CompuTherm, a very simple pr'//  &
!~                                           'ogram for the calculation of'
!~        write(*,'(9X,A)') ' thermodynamic properties from electroni'//  &
!~                                               'c structure calculations'
!~        write(*,*)
!~        write(*,'(1X,90("-"))')
!~        write(*,*)
!
       call print_start()
!
! Reading command line options
!
       call command_line(inp,outp)
!
       freq(:)  = 0.0d0
       fosc(:)  = 0.0d0
       inten(:) = 0.0d0
!
! Requesting user an input file if needed via GUI
!
       if ( len_trim(inp) .eq. 0 ) then
!
         CALL DWGFIL('Please, choose an input file',inp,'*.out')
!
!  The routine DWGERR returns a status for the routines DWGFIL, DWGTXT
!   and DWGLIS. The routine can be used to check directly after the
!   routines above if the OK button is pressed in the routines.
!
         CALL DWGERR(io)
!
         if ( io .ne. 0 ) then
           write(*,'(2X,68("="))')
           write(*,'(3X,A)') 'ERROR:  Error while reading input fi'//  &
                                                               'le name'
           write(*,*)
           write(*,'(3X,A)') 'An error occurred while reading inpu'//  &
                                                           't file name'
           write(*,'(3X,A)') 'Remember that you can choose an inpu'//  &
                                         't file name directly from the'
           write(*,'(3X,A)') ' command-line using the flag "--file"'
           write(*,'(2X,68("="))')
           write(*,*)
           call print_end()
         end if
!
         if ( len_trim(outp) .ne. 0 ) then
           if ( outp(len_trim(outp)-3:) .ne. '.dat' )                  &
                                               outp = trim(outp)//'.dat'
         else
           outp = inp(:len_trim(inp)-4)//'.dat'
         end if
!
       end if
!
! By default, the plot file name consists of the keyword 'dislin' and an
!  extension that depends on the file format. An alternate filename can 
!  be set with SETFIL.
!
       CALL SETFIL(outp(:len_trim(inp)-4)//'.png')
!
! Determines if a new plot file name is created for existing files.
!  'DELETE'	means that the existing file will be overwritten.
!
       CALL FILMOD('DELETE')
!
! Checking electronic structure program used for the calculation
!

!
! Reading input file
!
       call read_inp(inp)
!
       if ( ixunit .eq. 2 ) freq(:nband) = 1.0E7/freq(:nband)
!
! Reding reference file if requested
!
       if ( iref .eq. 1 ) call importref(0)
!
! Initializing theoretical spectrums
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
       call initialize()
!
       call colorate(illu,iadapt,2,nref,xref,yref,Xwp,Ywp,Zwp,         &
                     X,Y,Z,R,G,B)
!
! Experimental spectrum does not admidt standard illuminant method      ! FLAG: TODO
!



!
! Computing RGB coordinates of the theoretical spectrum if needed       ! FLAG: TODO
!



!
! Computing RGB coordinates associated with the experimental spectrum   ! FLAG: TODO
!



!
! Generating a GUI usin DISLIN
! ----------------------------
!
! Setting a title for the main window
!
       CALL SWGTIT('UV-vis Spectroscopy')
!
! Setting a character string that will be displayed if the Help menu
!  is clicked
!
       CALL SWGHLP('This program is made to plot beautiful spectrums')
!
! Sets widget options
!  'CLOSE' changes the behaviour of the close button of the main widget
!
       CALL SWGOPT('OK','CLOSE')
!
! Creating parent window (IP is the top level widget id)
!
       CALL WGINI('HORI',ip)
!
! Setting the default width of horizontal and parent/base widgets
!
       CALL SWGWTH(-15)
!
! Creating a container widget
!
       CALL WGBAS(ip,'VERT',ip1)
!
! Setting the default width of horizontal and parent/base widgets
!
       CALL SWGWTH(-70)
!
! Creating a container widget
!
       CALL WGBAS(ip,'VERT',ip2)
!
! Modifying the height of draw widgets
!  Using same height/width ratio as DINA4
!
       CALL SWGDRW(2100.0/2970.0)
!
! Creating a draw widget
!
       CALL WGDRAW(ip2,idplot)
!
! Generating popup menu entries
! -----------------------------
!
! Creates a popup menu in the menu bar of the main widget, or a popup
!  submenu of a popup menu
!
       CALL WGPOP(ip,'Data',iddata)
       CALL WGPOP(ip,'Axis',idaxis)
!
! Modifies the appearance of certain widgets
!  If CLASS = 'POPUP', CTYPE can have the values 'STRING' and 'MENU'.
!  STRING means that a popup menu can be directly connected with a
!  callback routine. Normally, menu entries in a popup menu can be
!  connected with callback routines
!
       CALL SWGTYP('STRING','POPUP')
       CALL WGPOP(ip,'Colour',idcolopt)
!
! Creating Data popup menu entries
! ................................
!
! Creates an entry in a popup menu
!
       CALL WGAPP(iddata,'Import',idimport)
!
       CALL WGAPP(iddata,'Export',idexport)
!
! Creating Axis popup menu entries
! ................................
!
       CALL WGAPP(idaxis,'Limits',idscl)
       CALL WGAPP(idaxis,'Units',idunits)
!
! Creating spectral lineshape options
! ...................................
!
! Creates a label widget to present the section
!
       CALL WGLAB(ip1,' ',idlab)
       CALL WGLAB(ip1,'Spectral lineshape',idlab)
!
! Select lineshape by pressing a radio button, just a small number to
!  choose from
!
       CALL WGBOX(ip1,'Gaussian|Lorentzian|Voigt',iline,idline)
!
! Select scaling factor and HWHM with the keyword
!  Labelled text box for input
!
       CALL WGSEP (ip1,idsep)
!
!~       CALL WGLTXT(ip1,'Scaling factor','1.0000',35,idfactor)
!
       CALL WGLTXT(ip1,'HWHM (eV)','0.0000',35,idhwhm)
       CALL WGLTXT(ip1,'Gaussian   HWHM','0.0000',35,idgauhwhm)
       CALL WGLTXT(ip1,'Lorentzian HWHM','0.0000',35,idlorhwhm)
!
       CALL REAWGT()
!
! Changes the value of a text widget as real number (returns the input)
!  Setting the Half Width at Half Maximum values
!
       if ( ixunit .eq. 1 ) then
         CALL SWGFLT(idhwhm,real(hwhm*cm2ev),4)
         CALL SWGFLT(idgauhwhm,real(gauhwhm*cm2ev),4)
         CALL SWGFLT(idlorhwhm,real(lorhwhm*cm2ev),4)
       else
         CALL SWGFLT(idhwhm,real(1.0E7/hwhm*cm2ev),4)
         CALL SWGFLT(idgauhwhm,real(1.0E7/gauhwhm*cm2ev),4)
         CALL SWGFLT(idlorhwhm,real(1.0E7/lorhwhm*cm2ev),4)
       end if
!
!~        CALL SWGFLT(idhwhm,hwhm,4) 
!~        CALL SWGFLT(idgauhwhm,gauhwhm,4)
!~        CALL SWGFLT(idlorhwhm,lorhwhm,4)
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
! Separates widgets by drawing horizontal or vertical lines, or menu
!  entries by drawing horizontal lines
!
       CALL WGSEP (ip1,idsep)
!
! Single button with label - click to turn on or off
!
       CALL WGBUT(ip1,'Sticks',istick,idstick)
!
       CALL WGBUT(ip1,'Convolute',iconv,idconv)
!
       CALL WGBUT(ip1,'Normalize',inorm,idnorm)
!
       CALL WGBUT(ip1,'Color',icolor,idcolor)
!
       CALL REAWGT()
!
       if ( iyunit .ne. 3 ) then
         CALL SWGATT(idnorm,'ACTIVE','STATUS')
       else
         CALL SWGATT(idnorm,'INACTIVE','STATUS')
       end if
!
! Creating data transformation options
! ....................................
!
!~       CALL WGLAB(ip1,' ',idlab)
!~       CALL WGLAB(ip1,'Data transformations',idlab)
!~ !
!~       yshift = 0.0
!~       CALL WGLTXT(ip1,'Vertical shift','0.00',40,idxshift)
!~ !
!~       xshift = 0.0
!~       CALL WGLTXT(ip1,'Horizontal shift','0.00',40,idyshift)
!
! Plotting spectra
! ................
!
       CALL WGSEP (ip1,idsep)
!
       CALL WGPBUT(ip1,'Plot',idbut)
!
! Callback routine - called when new list item selected or PLOT button pressed
!
       CALL SWGCBK(idbut,plot)
!
! Resetting the axes values
! .........................
!
       CALL WGPBUT(ip1,'Reset',idreset)
       CALL SWGCBK(idreset,resetplot)
!
! Add an OK button for user convenience
!
       CALL WGOK(ip1,idok)
!
!
! Creating experimental spectra options
! .....................................
!
       CALL WGSEP(ip1,idsepref)
       CALL WGLAB(ip1,'Experimental spectrum',idlabref)
!
! Removing experimental spectra
!
       CALL WGPBUT(ip1,'Remove',idnoref)
!
! Hidden experimental spectra if requested
!
       CALL WGBUT(ip1,'Hide',ihidden,idhidden)
!
! Colouring experimental spectra
!
       CALL WGBUT(ip1,'Color',icolref,idcolref)
!
! Realizes a widget tree. Since the windows ID of a widget can only be 
!  calculated for X11 if the widget is already realized, this routine is 
!  useful if the windows ID of a widget is needed before WGFIN.
!  Normally, the widget tree is realized in WGFIN.
!
       CALL REAWGT()
!
       if ( iref .eq. 0 ) then
         CALL SWGATT(idsepref,'INVISIBLE','STATUS')
         CALL SWGATT(idlabref,'INVISIBLE','STATUS')
         CALL SWGATT(idnoref,'INVISIBLE','STATUS')
         CALL SWGATT(idhidden,'INVISIBLE','STATUS')
         CALL SWGATT(idcolref,'INVISIBLE','STATUS')
       else
         CALL SWGATT(idsepref,'ACTIVE','STATUS')
         CALL SWGATT(idlabref,'ACTIVE','STATUS')
         CALL SWGATT(idnoref,'ACTIVE','STATUS')
         CALL SWGATT(idhidden,'ACTIVE','STATUS')
         CALL SWGATT(idcolref,'ACTIVE','STATUS')
       end if
!
! Convoluted spectra settings
! ---------------------------
!
! Sets the status of a button widget
!  Plot the spectrum without pushing the buttom "Plot"
!
       CALL SWGBUT(idbut,1)
!
! Changing line style
!
       CALL SWGCBK(idline,newplot)
!
! Changing half width at half maximum   ! FLAG: doesnt work when units menu is open
!
       CALL SWGCBK(idhwhm,newplot)
       CALL SWGCBK(idgauhwhm,newplot)
       CALL SWGCBK(idlorhwhm,newplot)
!
! Plotting stick spectra
!
       CALL SWGCBK(idstick,iniplot)
!
! Plotting convoluted spectra
!
       CALL SWGCBK(idconv,iniplot)
!
! Normalize experimental and theoretical spectrums
!
       CALL SWGCBK(idnorm,normalize)
!
! Plotting coloured spectra
!
!~        CALL SWGCBK(idcolor,plotcolor)
!
! Color manipulation options
! --------------------------
!
       CALL SWGCBK(idcolopt,colormenu)
!
! Graph manipulation settings
! ---------------------------
!
       CALL SWGCBK(idscl,axismenu)
!
       CALL SWGCBK(idunits,unitsmenu)
!
! Experimental spectra options
! ----------------------------
!
! Importing experimental spectra
!
       CALL SWGCBK(idimport,importref)
!
! Removing experimental spectra
!
       CALL SWGCBK(idnoref,removeref)
!
! Hiding experimental spectra
!
       CALL SWGCBK(idhidden,iniplot)
!
! Colouring experimental spectra
!
!~        CALL SWGCBK(idcolref,plotcolor)
!
! Terminating DISLIN
! ------------------
!
       CALL SWGCBK(idexport,iniplot)
!
       CALL WGFIN()
!
! Saving calculated spectrum in .dat format
!
       open(unit=uniout,file=trim(outp),action='write')
!
       write(*,'(2X,A)') 'Printing convoluted spectrum on file '//     &
                                                              trim(outp)
       write(*,*)
!
       do i = 1, npts
         if ( xline(i) .gt. xmaxaux ) exit
         if ( xline(i) .ge. xminaux ) write(uniout,*) xline(i),yline(i)
       end do
!
       close(uniout)
!
! Deallocating variables
!
       deallocate(xline,yline)
       deallocate(xstick,ystick)
       if ( allocated(xref) ) deallocate(xref,yref)
!
! Printing end message
!
       call print_end()
!
       end program ir_spectra
!
!======================================================================!
!
       subroutine command_line(inp,outp)
!
       use variables
       use parameters
       use utils
!
       implicit none
!
! Input/Output variables
!
       character(len=leninp),intent(out)  ::  inp      !  Input file name
       character(len=lenout),intent(out)  ::  outp     !  Output file name
!
! Local variables
!
       character(len=lencmd)              ::  cmd      !  Command executed
       character(len=lenarg)              ::  code     !  Executable name
       character(len=lenarg)              ::  arg      !  Argument read
       character(len=lenarg)              ::  next     !  Next argument to be read
       logical                            ::  flgdx    !
       logical                            ::  flgxmin  !
       logical                            ::  flgxmax  !
       integer                            ::  io       !  Status
       integer                            ::  i        !  Index
!
! Setting defaults
!
       inp   = ''
       outp  = ''
       ref   = ''
!
       iref     = 0     ! = 0 (No reference), = 1 (Reference on memory)
       imenu    = 0     ! = 0 (No axes menu), = 1 (Axes menu window on)
       iunits   = 0     ! = 0 (No units menu), = 1 (Units menu window on)
       icolmenu = 0     ! = 0 (No color menu), = 1 (Color menu window on)
!
       ihidden = 0      ! = 0 (Show reference), = 1 (Hide reference)
       icolref = 0      ! = 0 (Reference color off), = 1 (Color on)
       inorm   = 0      ! = 0 (Original spectrums), = 1 (Normalized spectrums)
!
       ixunit = 2       ! = 1 (Wavenumber), = 2 (Wavelength)
       iyunit = 1       ! = 1 (Epsilon), = 2 (Absorbance), = 3 (Transmittance)
       ilog   = 1       ! = 1 (Linear), = 2 (Logarithmic)
!
       istick = 1       ! = 0 (Sticks off), = 1 (Sticks on)
       iconv  = 1       ! = 0 (Convolution off), = 1 (Convolution on)
       iline  = 1       ! = 1 (Gaussian), = 2 (Lorentzian), = 3 (Voigt)
       icolor = 0       ! = 0 (Color off), = 1 (Color on)
!
       iillu  = 5       ! = 1 (A)  , = 2 (C)  , = 3 (D50), = 4( D55),
                        ! = 5 (D65), = 6 (D75), = 7 (E)  , = 8 (Inverse)
                        ! = 9 (Black Body)
!
       imodel = 16      ! = 1  (Adobe RGB)   , = 2  (Apple RGB)
                        ! = 3  (Best RGB)    , = 4  (Beta RGB) 
                        ! = 5  (Bruce RGB)   , = 6  (CIE RGB)
                        ! = 7  (Color match) , = 8  (Don RGB 4)
                        ! = 9  (ECI RGB v2)  , = 10 (Ekta Space PS5)
                        ! = 11 (HDTV RGB)    , = 12 (NTSC RGB)    
                        ! = 13 (EBU [PAL/SECAM] RGB) 
                        ! = 14 (ProPhoto RGB), = 15 (SMPTE RGB) 
                        ! = 16 (sRGB)        , = 17 (Wide Gamut RGB)
!
       iadapt = 1       ! = 1 (Bradford)   , = 2 (von Kries)
                        ! = 3 (XYZ Scaling), = 4 (None)
!
       factor = 1.0d0
!
       yrefmax = 0.0d0
!
       illu  = 'd65'
       model = 'srgb'
       cct   = 6504.0d0
!
       conc = 1.0E-3
       path = 1.0d0
!
       cold = 1.0d0
       lold = 1.0d0
!
       hwhm    = hwhmev
       gauhwhm = hwhmev
       lorhwhm = hwhmev
!
       dx     = dxnm
       xmin   = xminnm
       xmax   = xmaxnm
!
       ymin   = 0.0d0
       ymax   = 1.0d0
!
       xshift = 0.0d0
       yshift = 0.0d0
!
       flgdx   = .TRUE.
       flgxmin = .TRUE.
       flgxmax = .TRUE.
!
! Reading command line
!
       call get_command_argument(0,code)
       call get_command(cmd)
!
! Reading command line options
!
       i = 1
       do
         call get_command_argument(i,arg)
         if ( len_trim(arg) == 0 ) exit
         i = i+1
         select case ( arg )
           case ('-f','-file','--file')
!
             call get_command_argument(i,inp,status=io)
             call check_arg(inp,io,arg,cmd)
             i = i + 1
!
           case ('-o','-out','-outp','-output','--output')
!
             call get_command_argument(i,outp,status=io)
             call check_arg(outp,io,arg,cmd)
             i = i + 1
!
           case ('-r','-ref','--ref','-reference','--reference')
!
             call get_command_argument(i,ref,status=io)
             call check_arg(ref,io,arg,cmd)
             iref = 1
             i = i + 1
!
           case ('-s','--scaling','--factor','--scaling-factor')
!
             call get_command_argument(i,next,status=io)
             read(next,*) factor
             i = i + 1
!
           case ('-w','-hwhm','--hwhm')
!
             call get_command_argument(i,next,status=io)
             read(next,*) hwhm
             i = i + 1
!
           case ('-stick','--stick')
!
             istick = 1
!
           case ('-nostick','--nostick')
!
             istick = 0
!
           case ('-conv','--conv','--convolute')
!
             iconv = 1
!
           case ('-noconv','--noconv','--noconvolute')
!
             iconv = 0
!
           case ('-col','--color','--colour')
!
             icolor  = 1
             icolref = 1
!
           case ('-nocol','--nocolor','--nocolour')
!
             icolor  = 0
             icolref = 0
!
           case ('-i','-iluminant','-illuminant','--illuminant',       &
                                                          '--iluminant')
!
             call get_command_argument(i,next,status=io)
!
             next = lowercase(next)
             select case (trim(next))
               case('a')
                 illu  = 'a'
                 iillu = 1
               case('c')
                 illu  = 'c'
                 iillu = 2
               case('d50')
                 illu  = 'd50'
                 iillu = 3
               case('d55')
                 illu  = 'd55'
                 iillu = 4
               case('s','sunlight','standard','d65')
                 illu  = 'd65'
                 iillu = 5
               case('d75')
                 illu  = 'd75'
                 iillu = 6
               case('e')
                 illu  = 'e'
                 iillu = 7
               case('w','white','inverse','invert','beck','beck2005')
                 illu  = 'inverse'
                 iillu = 8
               case('bb','blackbody','black-body')
                 illu = 'black-body'
                 iillu = 9
               case default
                 write(*,*)
                 write(*,'(2X,68("="))')
                 write(*,'(3X,A)')    'ERROR:  Invalid value intro'//  &
                                         'duced for --illuminant option'
                 write(*,*)
                 write(*,'(3X,2(A))') 'Unrecognised value     : ',     &
                                                              trim(next)
                 write(*,*)
                 write(*,'(3X,A)')    'Please, to know the possibl'//  &
                                                     'e options execute'
                 write(*,*)
                 write(*,'(4X,2(A))')  code(1:len_trim(code)), ' -h'
                 write(*,'(2X,68("="))')
                 write(*,*)
                 call print_end()
             end select
!
             i = i + 1
!
           case ('-cct','--correlated-color-temperature',              &
                 '--temperature','--colour-temperature',               &
                 '--color-temperature')
!
             call get_command_argument(i,next,status=io)
             read(next,*) cct
             i = i + 1
!
           case ('-rgb','--rgb','--rgb-model','--rgb-space')
!
             call get_command_argument(i,next,status=io)
!
             next = lowercase(next)
             select case (trim(next))
               case('adobe','adobe-rgb')
                 model  = 'adobe'
                 imodel = 1
               case('apple','apple-rgb')
                 model  = 'apple'
                 imodel = 2
               case('best','best-rgb')
                 model  = 'best'
                 imodel = 3
               case('beta','beta-rgb')
                 model  = 'beta'
                 imodel = 4
               case('bruce','bruce-rgb')
                 model  = 'bruce'
                 imodel = 5
               case('cie',' cie-rgb')
                 model  = 'cie'
                 imodel = 6
               case('colormatch','colormatch-rgb')
                 model  = 'colormatch'
                 imodel = 7
               case('don','don-rgb','don-rgb-4')
                 model  = 'don'
                 imodel = 8
               case('eci','eci-rgb','eci-rgb-v2')
                 model  = 'eci'
                 imodel = 9
               case('ekta','ekta-rgb','ekta-space','ekta-space-ps5')
                 model  = 'ekta'
                 imodel = 10
               case('hdtv','hdtv-rgb')
                 model  = 'hdtv'
                 imodel = 11
               case('ntsc','ntsc-rgb')
                 model  = 'ntsc'
                 imodel = 12
               case('ebu','ebu-rgb','pal/secam','pal/secam-rgb')
                 model  = 'ebu'
                 imodel = 13
               case('prophoto','prophoto-rgb')
                 model  = 'prophoto'
                 imodel = 14
               case('smpte','smpte-rgb')
                 model  = 'smpte'
                 imodel = 15
               case('srgb','rec.709')
                 model  = 'srgb'
                 imodel = 16
               case('widegamut','wide-gamut','wide-rgb',               &
                    'widegamut-rgb','wide-gamut-rgb')
                 model  = 'widegamut'
                 imodel = 17
               case default
                 write(*,*)
                 write(*,'(2X,68("="))')
                 write(*,'(3X,A)')    'ERROR:  Invalid value intro'//  &
                                          'duced for --rgb-space option'
                 write(*,*)
                 write(*,'(3X,2(A))') 'Unrecognised value     : ',     &
                                                              trim(next)
                 write(*,*)
                 write(*,'(3X,A)')    'Please, to know the possibl'//  &
                                                     'e options execute'
                 write(*,*)
                 write(*,'(4X,2(A))')  code(1:len_trim(code)), ' -h'
                 write(*,'(2X,68("="))')
                 write(*,*)
                 call print_end()
             end select
!
             i = i + 1
           case ('-adapt','--adaptation')
!
             call get_command_argument(i,next,status=io)
!
             next = lowercase(next)
             select case (trim(next))
               case('bradford')
                 iadapt = 1
               case('vonkries','von-kries')
                 iadapt = 2
               case('xyz','xyzscaling','xyz-scaling')
                 iadapt = 3
               case('none')
                 iadapt = 4
               case default
                 write(*,*)
                 write(*,'(2X,68("="))')
                 write(*,'(3X,A)')    'ERROR:  Invalid value intro'//  &
                                         'duced for --adaptation option'
                 write(*,*)
                 write(*,'(3X,2(A))') 'Unrecognised value     : ',     &
                                                              trim(next)
                 write(*,*)
                 write(*,'(3X,A)') 'Please, choose between : "brad'//  &
                            'ford", "vonkries", "xyz-scaling" or "none"'
                 write(*,*)
                 write(*,'(4X,2(A))')  code(1:len_trim(code)), ' -h'
                 write(*,'(2X,68("="))')
                 write(*,*)
                 call print_end()
             end select
!
             i = i + 1
!
           case ('-c','-conc','--conc','--concentration')
!
             call get_command_argument(i,next,status=io)
             read(next,*) conc
             i = i + 1
!
           case ('-b','-pl','-path','--pathlength','--path-length ')
!
             call get_command_argument(i,next,status=io)
             read(next,*) path
             i = i + 1
!
           case ('-l','-lineshape','--lineshape')
!
             call get_command_argument(i,next,status=io)
!
             next = lowercase(next)
             select case (trim(next))
               case('g','gau','gauss','gaussian')
                 iline = 1
               case('l','lor','lorentz','lorentzian')
                 iline = 2
               case('v','voigt')
                 iline = 3
               case default
                 write(*,*)
                 write(*,'(2X,68("="))')
                 write(*,'(3X,A)')    'ERROR:  Invalid value intro'//  &
                                          'duced for --lineshape option'
                 write(*,*)
                 write(*,'(3X,2(A))') 'Unrecognised value     : ',     &
                                                              trim(next)
                 write(*,*)
                 write(*,'(3X,A)') 'Please, choose between : "gaus'//  &
                                        'sian", "lorentzian" or "voigt"'
                 write(*,'(2X,68("="))')
                 write(*,*)
                 call print_end()
             end select
!
             i = i + 1
!
           case ('-dx','-reso','--reso','--resolution')
!
             call get_command_argument(i,next,status=io)
             read(next,*) dx
             flgdx = .FALSE.
             i = i + 1
!
           case ('-xmin','--xmin','--minimum-x')
!
             call get_command_argument(i,next,status=io)
             read(next,*) xmin
             flgxmin = .FALSE.
             i = i + 1
!
           case ('-xmax','--xmax','--maximum-x')
!
             call get_command_argument(i,next,status=io)
             read(next,*) xmax
             flgxmax = .FALSE.
             i = i + 1
!
           case ('-ymin','--ymin','--minimum-y')
!
             call get_command_argument(i,next,status=io)
             read(next,*) ymin
             i = i + 1
!
           case ('-ymax','--ymax','--maximum-y')
!
             call get_command_argument(i,next,status=io)
             read(next,*) ymax
             i = i + 1
!
           case ('-xshift','--xshift','-hshift','--hshift',            &
                                                   '--horizontal-shift')
!
             call get_command_argument(i,next,status=io)
             read(next,*) xshift
             i = i + 1
!
           case ('-yshift','--yshift','-vshift','--vshift',            &
                                                     '--vertical-shift')
!
             call get_command_argument(i,next,status=io)
             read(next,*) yshift
             i = i + 1
!
           case ('-scl','-scale','--scl','--scale')
!
             call get_command_argument(i,next,status=io)
!
             next = lowercase(next)
             select case (trim(next))
               case('lin','linear')
                 ilog = 1
               case('log','logarithm','logarithmic')
                 ilog = 2
               case default
                 write(*,*)
                 write(*,'(2X,68("="))')
                 write(*,'(3X,A)')    'ERROR:  Invalid value intro'//  &
                                              'duced for --scale option'
                 write(*,*)
                 write(*,'(3X,2(A))') 'Unrecognised value     : ',     &
                                                              trim(next)
                 write(*,*)
                 write(*,'(3X,A)') 'Please, choose between : "line'//  &
                                                  'ar" or "logarithmic"'
                 write(*,'(2X,68("="))')
                 write(*,*)
                 call print_end()
             end select
!
             i = i + 1
!
           case ('-u','-unit','-units','--unit','--units')
!
             call get_command_argument(i,next,status=io)
!
             next = lowercase(next)
!
             select case (trim(next))
!
               case('cm-1','cm**-1','wavenumber','wvn')
!
                 ixunit = 1
!
                 if ( flgdx )   dx   = dxcm
                 if ( flgxmin ) xmin = xmincm
                 if ( flgxmax ) xmax = xmaxcm
!
               case('nm','wavelength','wvl')
!
                 ixunit = 2
!
                 if ( flgdx )   dx   = dxnm
                 if ( flgxmin ) xmin = xminnm
                 if ( flgxmax ) xmax = xmaxnm
!
               case default
!
                 write(*,*)
                 write(*,'(2X,68("="))')
                 write(*,'(3X,A)')    'ERROR:  Invalid value intro'//  &
                                          'duced for --units option'
                 write(*,*)
                 write(*,'(3X,2(A))') 'Unrecognised value     : ',     &
                                                              trim(next)
                 write(*,*)
                 write(*,'(3X,A)') 'Please, choose between : "wave'//  &
                                               'number" or "wavelength"'
                 write(*,'(2X,68("="))')
                 write(*,*)
                 call print_end()
!
             end select
!
             i = i + 1
!
           case ('-h','-help','--help')
!
             call print_help()
             call print_end()
!
           case default
!
             write(*,*)
             write(*,'(2X,68("="))')
             write(*,'(3X,A)')    'ERROR:  Unknown statements from'//  &
                                  ' command line'
             write(*,*)
             write(*,'(4X,A)')     trim(cmd)
             write(*,*)
             write(*,'(3X,2(A))') 'Unrecognised command-line option'// &
                                  '  :  ', arg
             write(*,'(3X,A)')    'Please, to request help execute'
             write(*,*)
             write(*,'(4X,2(A))')  code(1:len_trim(code)), ' -h'
             write(*,'(2X,68("="))')
             write(*,*)
             call print_end()
!
         end select
       end do
!
! General settings and fatal errors checks
!
       if ( len_trim(outp) .ne. 0 ) then
         if ( outp(len_trim(outp)-3:) .ne. '.dat' )                    &
                                               outp = trim(outp)//'.dat'
       else if ( len_trim(inp) .ne. 0 ) then
         outp = inp(:len_trim(inp)-4)//'.dat'
       end if
!
       if ( iyunit .ne. 1 ) then
         cold = conc
         lold = path
       end if
!
!~        if ( ixunit .eq. 1 ) then
!~          hwhm    = hwhm*ev2cm
!~          gauhwhm = hwhm*ev2cm
!~          lorhwhm = hwhm*ev2cm
!~        else
!~          hwhm    = 1.0E7/(hwhm*ev2cm)
!~          gauhwhm = 1.0E7/(hwhm*ev2cm)
!~          lorhwhm = 1.0E7/(hwhm*ev2cm)
!~        end if
!
       return
       end subroutine command_line
!
!======================================================================!
!
       subroutine read_inp(inp)
!
       use variables
       use parameters
       use utils
!
       implicit none
!
! Input/Output variables
!
       character(len=leninp),intent(in)  ::  inp     !  Input file name
!
! Local variables
!
       character(len=lenline)            ::  line    !  Line read
       character(len=lenline)            ::  key
       real(kind=4)                      ::  s1
       integer                           ::  i1
       integer                           ::  keylen  !  Keyword length
       integer                           ::  io      !  Input/Output status
!
!
!
       open(unit=uniinp,file=trim(inp),action='read',                  &
                                                 status='old',iostat=io)
!
       if ( io .ne. 0 ) then
         write(*,'(2X,68("="))')
         write(*,'(3X,A)') 'ERROR:  Missing input file'
         write(*,*)
         write(*,'(3X,A)') '  Input file "'//trim(inp)//'" not found'
         write(*,'(2X,68("="))')
         write(*,*)
         call print_end()
       end if
!
       key = 'ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE M'//  &
                                                                'OMENTS'
!
       keylen = len(key)
!
       do while ( io .eq. 0 )
!
         do
           read(uniinp,'(A)',iostat=io) line!
!~          if ( io .ne. 0 ) then
!~            write(*,'(2X,68("="))')
!~            write(*,'(3X,A)') 'ERROR:  Missing information in the i'//  &
!~                                                              'nput file'
!~            write(*,*)
!~            write(*,'(3X,A)') 'Excited states information was not f'//  &
!~                                                 'ound in the input file'
!~            write(*,'(3X,A)') 'Please, check that the input contain'//  &
!~                                              's the correct calculation'
!~            write(*,'(2X,68("="))')
!~            write(*,*)
!~            call print_end()
!~          end if
!
           line = adjustl(line)
           if ( io .ne. 0 ) exit
           if ( line(:keylen) .eq. trim(key) ) exit
         end do
!
         if ( io .ne. 0 ) exit
!
         read(uniinp,'(A)') line
         read(uniinp,'(A)') line
         read(uniinp,'(A)') line
         read(uniinp,'(A)') line
!
         nband = 1
         read(uniinp,'(A)') line
!
         do
           read(line,*) i1,freq(nband),s1,fosc(nband),inten(nband)
!~          read(line,*) i1,s1,freq(nband),fosc(nband),inten(nband)
           read(uniinp,'(A)') line
           if ( len_trim(line) .eq. 0 ) exit
           nband = nband + 1
         end do
!
       end do
!
write(*,*) freq(:nband)
!
!~ write(*,*) 'BEFORE',inten(:nband)
!~        inten(:) = (qe*a0)**2*inten(:)           ! FLAG: change of units
!~ write(*,*) 'AFTER ',inten(:nband)
!
!~ write(*,*) 'DIPOLE MOMENT AU2SI ',qe*a0
!~ write(*,*) 'PREFACTOR SI UNITS  ',Aabs
!~ write(*,*) 'PREFACTOR AU        ',Aabs*(qe*a0)**2
!~ write(*,*) 'PREAFACTOR JAVIER   ',2.5115E22*a0**2*10
!
       close(uniinp)
!
       return
       end subroutine read_inp
!
!======================================================================!
!
       subroutine print_help()
!
       implicit none
!
       write(*,'(2X,68("="))')
       write(*,'(3X,A)') 'Command-line options:'
       write(*,*)
       write(*,'(5X,A)') '-h,--help             Print usage inform'//  &
                                                         'tion and exit'
       write(*,'(5X,A)') '-f,--file             Input file name'
       write(*,'(5X,A)') '-o,--output           Output file name'
       write(*,'(5X,A)') '-r,--reference        Reference file name'
       write(*,*)
       write(*,'(5X,A)') '--[no]stick           Print stick spectrum'
       write(*,'(5X,A)') '--[no]convolute       Print convoluted s'//  &
                                                               'pectrum'
       write(*,*)
       write(*,'(5X,A)') '-l,--lineshape        Lineshape function type'
       write(*,'(5X,A)') '                       (gauss | lor | voigt)'
       write(*,'(5X,A)') '-s,--scaling-factor   Frequencies scalin'//  &
                                                              'g factor'
       write(*,'(5X,A)') '-w,--hwhm             Half width at half'//  &
                                                         ' maximum (eV)'
       write(*,*)
       write(*,'(5X,A)') '--[no]colour          Print coloured spectrum'
       write(*,'(5X,A)') '-i,--illuminant       Standard illuminan'//  &
                                                                't type'
       write(*,'(5X,A)') '                       ( ... | ... )'
       write(*,'(5X,A)') '-rgb,--rgb-space      RGB working space'
       write(*,'(5X,A)') '                       ( ... | ... )'
       write(*,'(5X,A)') '-adapt,--adaptation   Chromatic adaptati'//  &
                                                             'on method'
       write(*,'(5X,A)') '                       ( ... | ... )'
       write(*,'(5X,A)') '-cct,--temperature    Correlated color t'//  &
                                                        'emperature (K)'
       write(*,*)
       write(*,'(5X,A)') '-c,--concentration    Concentration (mol'//  &
                                                               '*L**-1)'
       write(*,'(5X,A)') '-b,--path-length      Path-length (cm)'
       write(*,*)
       write(*,'(5X,A)') '-u,--units            X-axis units'
       write(*,'(5X,A)') '                       (wavenumber | wav'//  &
                                                              'elength)'
       write(*,'(5X,A)') '-scl,--scale          Y-axis scale'
       write(*,'(5X,A)') '                       (linear | logarithmic)'
       write(*,*)
       write(*,'(5X,A)') '-reso,--resolution    Spectrum resolutio'//  &
                                                                'n (nm)'
       write(*,'(5X,A)') '-xmin,--minimum-x     Lower limit of the'//  &
                                                               ' X-axis'
       write(*,'(5X,A)') '-xmax,--maximum-x     Upper limit of the'//  &
                                                               ' X-axis'
       write(*,'(5X,A)') '-ymin,--minimum-y     Lower limit of the'//  &
                                                               ' Y-axis'
       write(*,'(5X,A)') '-ymax,--maximum-y     Upper limit of the'//  &
                                                               ' Y-axis'
       write(*,*)
       write(*,'(5X,A)') '-xshift,--hshift      Horizontal shift'
       write(*,'(5X,A)') '-yshift,--vshift      Vertical shift'
       write(*,'(2X,68("="))')
       write(*,*)
!
       return
       end subroutine print_help
!
!======================================================================!
