!======================================================================!
!
       module variables
!
       implicit none
!
! Molecular system information
!
       integer,parameter                      ::  nmax = 200
       real(kind=4),dimension(nmax)           ::  freq      !  Transition energies
       real(kind=4),dimension(nmax)           ::  inten     !  Absorbance
       real(kind=4),dimension(nmax)           ::  fosc      !
       integer                                ::  nat       !  Number of atoms
       integer                                ::  nband     !
!
! Spectral lineshape variables
!
       real(kind=4),parameter                 ::  hwhmev = 0.3
!
       character(len=256)                     ::  ref        !  Reference file name
       real(kind=4)                           ::  factor     !  Frequencies scaling factor
       real(kind=4)                           ::  hwhm       !  Half width at half maximum
       real(kind=4)                           ::  gauhwhm    !  HWHM for Gaussian broadening
       real(kind=4)                           ::  lorhwhm    !  HWHM for Lorentzian broadening
       real(kind=4)                           ::  CCT        !  Correlated colour temperature
       integer                                ::  iline      !  Lineshape identifier
       integer                                ::  istick     !  Stick spectra option
       integer                                ::  iconv      !  Convolution option
       integer                                ::  icolor     !  Color popup menu status
       integer                                ::  icolmenu   !  Color option
       integer                                ::  iillu      !  Illuminant option
       integer                                ::  imodel     !  RGB model option (working space)
       integer                                ::  iadapt     !  Adapation method option
       integer                                ::  ilog       !  Scaling option
       integer                                ::  iunits     !  Units popup menu status
       integer                                ::  ixunit     !  X-axis units option
       integer                                ::  iyunit     !  Y-axis units option
       integer                                ::  imenu      !  Axis popup menu status
!
! Experimental spectra options
!
       integer                                ::  iref       !
       integer                                ::  ihidden    !
       integer                                ::  icolref    !
       integer                                ::  inorm      !
!
! Curves information
!
       real(kind=4),dimension(:),allocatable  ::  xline      !  Frequencies values
       real(kind=4),dimension(:),allocatable  ::  yline      !  Convoluted spectral lineshape
       real(kind=4),dimension(:),allocatable  ::  xstick     !  Stick spectra frequencies
       real(kind=4),dimension(:),allocatable  ::  ystick     !  Stick spectra intensities
       real(kind=4),dimension(:),allocatable  ::  xref       !  Reference spectra frequencies
       real(kind=4),dimension(:),allocatable  ::  yref       !  Reference spectra intensities
       integer                                ::  npts       !  Number of convoluted spectrum points
       integer                                ::  nstick     !  Number of stick spectrum points
       integer                                ::  nref       !  Number of reference specutrum points
!
! Graph manipulation variables
!
       real(kind=4),parameter                 ::  fscl = 5.0E-4
!
       real(kind=4),parameter                 ::  dxcm   = 0.5d0
       real(kind=4),parameter                 ::  xmincm = 10000.0d0
       real(kind=4),parameter                 ::  xmaxcm = 1000000.0d0
!
       real(kind=4),parameter                 ::  dxnm   = 0.01d0
!~        real(kind=4),parameter                 ::  dxnm   = 5.0d0
       real(kind=4),parameter                 ::  xminnm = 0.0d0
       real(kind=4),parameter                 ::  xmaxnm = 1000.0d0
!
       real(kind=4),parameter                 ::  xsclmincm = 0.0d0
       real(kind=4),parameter                 ::  xsclmaxcm = 1000000.0d0
!
       real(kind=4),parameter                 ::  xsclminnm = 0.0d0
       real(kind=4),parameter                 ::  xsclmaxnm = 800.0d0
!
       real(kind=4),dimension(2)              ::  xscl       !  Plotting limits of the X-axis
       real(kind=4),dimension(2)              ::  yscl       !  Plotting limits of the Y-axis
       real(kind=4),dimension(2)              ::  xold       !
       real(kind=4),dimension(2)              ::  yold       !
       real(kind=4)                           ::  dx         !  Spacing between points of the X-axis
       real(kind=4)                           ::  xmin       !  Lower limit of the X-axis
       real(kind=4)                           ::  xmax       !  Upper limit of the X-axis
       real(kind=4)                           ::  xor        !  First X-axis label
       real(kind=4)                           ::  xstp       !  Step between labels of the X-axis
       real(kind=4)                           ::  ymin       !  Lower limit of the Y-axis
       real(kind=4)                           ::  ymax       !  Upper limit of the Y-axis
       real(kind=4)                           ::  yor        !  First Y-axis label
       real(kind=4)                           ::  ystp       !  Step between labels of the Y-axis
       real(kind=4)                           ::  xtail      !
       integer                                ::  nxdig      !
       integer                                ::  nydig      !
!
       real(kind=4)                           ::  xminaux    !
       real(kind=4)                           ::  xmaxaux    !
       real(kind=4)                           ::  yminaux    !
       real(kind=4)                           ::  ymaxaux    !
!
! Colouring spectrum variables
!
       character(len=10)                      ::  illu       !  Standard illuminant type
       character(len=10)                      ::  model      !  RGB working space
       real(kind=4)                           ::  X,Y,Z      !  Theoretical tristimulus values
       real(kind=4)                           ::  XX,YY,ZZ   !  Experimental tristimulus values
       real(kind=4)                           ::  Xwp        !  Reference white point tristimulus values
       real(kind=4)                           ::  Ywp        !  Reference white point tristimulus values
       real(kind=4)                           ::  Zwp        !  Reference white point tristimulus values
       real(kind=4)                           ::  R,G,B      !  Theoretical RGB values
       real(kind=4)                           ::  RR,GG,BB   !  Experimental RGB values
       real(kind=4)                           ::  conc       !  Molar concentration
       real(kind=4)                           ::  path       !  Path-length
       real(kind=4)                           ::  cold       !
       real(kind=4)                           ::  lold       !
!
! Data transformation variables
!  
       real(kind=4)                           ::  xshift     !  Horizontal shift
       real(kind=4)                           ::  yshift     !  Vertical shift
!
!
!
       integer                                ::  ip         !
       integer                                ::  ip1        !
       integer                                ::  ip2        !
       integer                                ::  ip3        !
       integer                                ::  idaxsmenu
!~        integer                                ::  ipmenu     !
!
! Plotting information identifiers
!
       integer                                ::  idline     !
       integer                                ::  idfactor   !
       integer                                ::  idhwhm     !
       integer                                ::  idgauhwhm  !
       integer                                ::  idlorhwhm  !
       integer                                ::  idstick    !
       integer                                ::  idconv     !
       integer                                ::  idcolor    !
       integer                                ::  idxshift   !
       integer                                ::  idyshift   !
       integer                                ::  idlog      !
       integer                                ::  idconc     !
       integer                                ::  idpath     !
       integer                                ::  idunits    !
       integer                                ::  idxunit    !
       integer                                ::  idyunit    !
       integer                                ::  idillu     !  Illuminant
       integer                                ::  idmodel    !  RGB model
       integer                                ::  idadapt    !  Adaptation method
       integer                                ::  idcct      !  Correlated color temperature
!
! Axis properties identifiers
!
       integer                                ::  idplot     !
       integer                                ::  idref      !
       integer                                ::  idreso     !
       integer                                ::  idxmin     !
       integer                                ::  idxmax     !
       integer                                ::  idymin     !
       integer                                ::  idymax     !
       integer                                ::  idstep     !
!
! Popup menu identifiers
!
       integer                                ::  iddata     !
       integer                                ::  idaxis     !
       integer                                ::  idcolopt   !
       integer                                ::  idimport   !
       integer                                ::  idexport   !
       integer                                ::  idscl      !
!
! Experimental spectra identifiers
!
       real(kind=4)                           ::  ylinemax   !
       real(kind=4)                           ::  ylinemaxl  !
       real(kind=4)                           ::  yrefmax    !
       real(kind=4)                           ::  yrefmaxl   !
       integer                                ::  idlabref   !
       integer                                ::  idsepref   !
       integer                                ::  idhidden   !
       integer                                ::  idnoref    !
       integer                                ::  idcolref   !
       integer                                ::  idnorm     !
!
! Widget identifiers
!
       integer                                ::  idbut      !
       integer                                ::  idreset    !
       integer                                ::  idlab      !
       integer                                ::  idsep      !
       integer                                ::  idok       !
       integer                                ::  iwidth     !
       integer                                ::  iheight    !
!
       end module variables
!
!======================================================================!
