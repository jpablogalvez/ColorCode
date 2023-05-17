!=======================================================================
!
! References
! ----------
!
! - Chris Wyman, Peter-Pike Sloan, and Peter Shirley, Simple Analytic
!    Approximations to the CIE XYZ Color Matching Functions, Journal 
!    of Computer Graphics Techniques (JCGT), vol. 2, no. 2, 1-11, 2013.
!    http://jcgt.org/published/0002/02/01/
!
! Fitting coefficients for multi-lobe Gaussian fit to the
!  X CIE1931 color-matching function
!
       real(kind=4),parameter  ::  x0a = 0.362
       real(kind=4),parameter  ::  x0b = 442.0
       real(kind=4),parameter  ::  x0g = 0.0624
       real(kind=4),parameter  ::  x0d = 0.0374
!
       real(kind=4),parameter  ::  x1a = 1.056
       real(kind=4),parameter  ::  x1b = 599.8
       real(kind=4),parameter  ::  x1g = 0.0264
       real(kind=4),parameter  ::  x1d = 0.0323
!
       real(kind=4),parameter  ::  x2a = -0.065
       real(kind=4),parameter  ::  x2b = 501.1
       real(kind=4),parameter  ::  x2g = 0.0490
       real(kind=4),parameter  ::  x2d = 0.0382
!
! Fitting coefficients for multi-lobe Gaussian fit to the
!  Y CIE1931 color-matching function
!
       real(kind=4),parameter  ::  y0a = 0.821
       real(kind=4),parameter  ::  y0b = 568.8
       real(kind=4),parameter  ::  y0g = 0.0213
       real(kind=4),parameter  ::  y0d = 0.0247
!
       real(kind=4),parameter  ::  y1a = 0.286
       real(kind=4),parameter  ::  y1b = 530.9
       real(kind=4),parameter  ::  y1g = 0.0613
       real(kind=4),parameter  ::  y1d = 0.0322
!
! Fitting coefficients for multi-lobe Gaussian fit to the
!  Z CIE1931 color-matching function
!
       real(kind=4),parameter  ::  z0a = 1.217
       real(kind=4),parameter  ::  z0b = 437.0
       real(kind=4),parameter  ::  z0g = 0.0845
       real(kind=4),parameter  ::  z0d = 0.0278
!
       real(kind=4),parameter  ::  z1a = 0.681
       real(kind=4),parameter  ::  z1b = 459.0
       real(kind=4),parameter  ::  z1g = 0.0385
       real(kind=4),parameter  ::  z1d = 0.0725
!
!=======================================================================
