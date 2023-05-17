!==========================================================
!
! References
! ----------
!
! - ASTM E308-01. Standard Practice for Computing the Colors of  
!    Objects by Using the CIE System; ASTM International: West 
!    Conshohoken, PA, 2001
! - http://www.brucelindbloom.com/
!
       module illuminants
!
       implicit none
!       
       contains
!
!======================================================================!
!
       subroutine illuminantA(xspec,yspec,Xw,Yw,Zw)
!
       implicit none    
!
! Input/Output variables
!
       real(kind=4),dimension(81),intent(out)  ::  xspec     !
       real(kind=4),dimension(81),intent(out)  ::  yspec     !
       real(kind=4),intent(out)                ::  Xw,Yw,Zw  !  White point coordinates
!
! Illuminant white point tristimulus values
!
       Xw = 1.09850 
       Yw = 1.00000 
       Zw = 0.35585
!
! Relative Spectral Power Distributions of CIE Standard Illuminant A
!  at 5-nm Intervals from 380 to 780 nm
!
       call setlamb(xspec)
!
       yspec(1)  = 0.0980
       yspec(2)  = 0.1090
       yspec(3)  = 0.1209
       yspec(4)  = 0.1335
       yspec(5)  = 0.1471
       yspec(6)  = 0.1615
       yspec(7)  = 0.1768
       yspec(8)  = 0.1929
       yspec(9)  = 0.2099
       yspec(10) = 0.2279
       yspec(11) = 0.2467
       yspec(12) = 0.2664
       yspec(13) = 0.2870
       yspec(14) = 0.3085
       yspec(15) = 0.3309
       yspec(16) = 0.3541
       yspec(17) = 0.3781
       yspec(18) = 0.4030
       yspec(19) = 0.4287
       yspec(20) = 0.4552
       yspec(21) = 0.4824
       yspec(22) = 0.5104
       yspec(23) = 0.5391
       yspec(24) = 0.5685
       yspec(25) = 0.5986
       yspec(26) = 0.6293
       yspec(27) = 0.6606
       yspec(28) = 0.6925
       yspec(29) = 0.7250
       yspec(30) = 0.7579
       yspec(31) = 0.7913
       yspec(32) = 0.8252
       yspec(33) = 0.8595
       yspec(34) = 0.8941
       yspec(35) = 0.9291
       yspec(36) = 0.9644
       yspec(37) = 1.0000
       yspec(38) = 1.0358
       yspec(39) = 1.0718
       yspec(40) = 1.1080
       yspec(41) = 1.1444
       yspec(42) = 1.1808
       yspec(43) = 1.2173
       yspec(44) = 1.2539
       yspec(45) = 1.2904
       yspec(46) = 1.3270
       yspec(47) = 1.3635
       yspec(48) = 1.3999
       yspec(49) = 1.4362
       yspec(50) = 1.4724
       yspec(51) = 1.5084
       yspec(52) = 1.5442
       yspec(53) = 1.5798
       yspec(54) = 1.6152
       yspec(55) = 1.6503
       yspec(56) = 1.6851
       yspec(57) = 1.7196
       yspec(58) = 1.7538
       yspec(59) = 1.7877
       yspec(60) = 1.8212
       yspec(61) = 1.8543
       yspec(62) = 1.8870
       yspec(63) = 1.9193
       yspec(64) = 1.9512
       yspec(65) = 1.9826
       yspec(66) = 2.0136
       yspec(67) = 2.0441
       yspec(68) = 2.0741
       yspec(69) = 2.1036
       yspec(70) = 2.1327
       yspec(71) = 2.1612
       yspec(72) = 2.1892
       yspec(73) = 2.2167
       yspec(74) = 2.2436
       yspec(75) = 2.2700
       yspec(76) = 2.2959
       yspec(77) = 2.3212
       yspec(78) = 2.3459
       yspec(79) = 2.3701
       yspec(80) = 2.3937
       yspec(81) = 2.4168
!
       return       
       end subroutine illuminantA
!
!======================================================================!
!
       subroutine illuminantC(xspec,yspec,Xw,Yw,Zw)
!
       implicit none    
!
! Input/Output variables
!
       real(kind=4),dimension(81),intent(out)  ::  xspec     !
       real(kind=4),dimension(81),intent(out)  ::  yspec     !
       real(kind=4),intent(out)                ::  Xw,Yw,Zw  !  White point coordinates
!
! Illuminant white point tristimulus values
!
       Xw = 0.98074 
       Yw = 1.00000 
       Zw = 1.18232
!
! Relative Spectral Power Distributions of CIE Standard Illuminant C
!  at 5-nm Intervals from 380 to 780 nm
!
       call setlamb(xspec)
!
       yspec(1)  = 0.3300
       yspec(2)  = 0.3992
       yspec(3)  = 0.4740
       yspec(4)  = 0.5517
       yspec(5)  = 0.6330
       yspec(6)  = 0.7181
       yspec(7)  = 0.8060
       yspec(8)  = 0.8953
       yspec(9)  = 0.9810
       yspec(10) = 1.0580
       yspec(11) = 1.1240
       yspec(12) = 1.1775
       yspec(13) = 1.2150
       yspec(14) = 1.2345
       yspec(15) = 1.2400
       yspec(16) = 1.2360
       yspec(17) = 1.2310
       yspec(18) = 1.2330
       yspec(19) = 1.2380
       yspec(20) = 1.2409
       yspec(21) = 1.2390
       yspec(22) = 1.2292
       yspec(23) = 1.2070
       yspec(24) = 1.1690
       yspec(25) = 1.1210
       yspec(26) = 1.0698
       yspec(27) = 1.0230
       yspec(28) = 0.9881
       yspec(29) = 0.9690
       yspec(30) = 0.9678
       yspec(31) = 0.9800
       yspec(32) = 0.9994
       yspec(33) = 1.0210
       yspec(34) = 1.0395
       yspec(35) = 1.0520
       yspec(36) = 1.0567
       yspec(37) = 1.0530
       yspec(38) = 1.0411
       yspec(39) = 1.0230
       yspec(40) = 1.0015
       yspec(41) = 0.9780
       yspec(42) = 0.9543
       yspec(43) = 0.9320
       yspec(44) = 0.9122
       yspec(45) = 0.8970
       yspec(46) = 0.8883
       yspec(47) = 0.8840
       yspec(48) = 0.8819
       yspec(49) = 0.8810
       yspec(50) = 0.8806
       yspec(51) = 0.8800
       yspec(52) = 0.8786
       yspec(53) = 0.8780
       yspec(54) = 0.8799
       yspec(55) = 0.8820
       yspec(56) = 0.8820
       yspec(57) = 0.8790
       yspec(58) = 0.8722
       yspec(59) = 0.8630
       yspec(60) = 0.8530
       yspec(61) = 0.8400
       yspec(62) = 0.8221
       yspec(63) = 0.8020
       yspec(64) = 0.7824
       yspec(65) = 0.7630
       yspec(66) = 0.7436
       yspec(67) = 0.7240
       yspec(68) = 0.7040
       yspec(69) = 0.6830
       yspec(70) = 0.6630
       yspec(71) = 0.6440
       yspec(72) = 0.6280
       yspec(73) = 0.6150
       yspec(74) = 0.6020
       yspec(75) = 0.5920
       yspec(76) = 0.5850
       yspec(77) = 0.5810
       yspec(78) = 0.5800
       yspec(79) = 0.5820
       yspec(80) = 0.5850
       yspec(81) = 0.5910
!
       return       
       end subroutine illuminantC
!
!======================================================================!
!
       subroutine illuminantD50(xspec,yspec,Xw,Yw,Zw)
!
       implicit none    
!
! Input/Output variables
!
       real(kind=4),dimension(81),intent(out)  ::  xspec     !
       real(kind=4),dimension(81),intent(out)  ::  yspec     !
       real(kind=4),intent(out)                ::  Xw,Yw,Zw  !  White point coordinates
!
! Illuminant white point tristimulus values
!
       Xw = 0.96422
       Yw = 1.00000 
       Zw = 0.82521
!
! Relative Spectral Power Distributions of CIE Standard Illuminant D50
!  at 5-nm Intervals from 380 to 780 nm
!
       call setlamb(xspec)
!
       yspec(1)  = 0.2449
       yspec(2)  = 0.2718
       yspec(3)  = 0.2987
       yspec(4)  = 0.3959
       yspec(5)  = 0.4931
       yspec(6)  = 0.5291
       yspec(7)  = 0.5651
       yspec(8)  = 0.5827
       yspec(9)  = 0.6003
       yspec(10) = 0.5893
       yspec(11) = 0.5782
       yspec(12) = 0.6632
       yspec(13) = 0.7482
       yspec(14) = 0.8104
       yspec(15) = 0.8725
       yspec(16) = 0.8893
       yspec(17) = 0.9061
       yspec(18) = 0.9099
       yspec(19) = 0.9137
       yspec(20) = 0.9324
       yspec(21) = 0.9511
       yspec(22) = 0.9354
       yspec(23) = 0.9196
       yspec(24) = 0.9384
       yspec(25) = 0.9572
       yspec(26) = 0.9617
       yspec(27) = 0.9661
       yspec(28) = 0.9687
       yspec(29) = 0.9713
       yspec(30) = 0.9961
       yspec(31) = 1.0210
       yspec(32) = 1.0143
       yspec(33) = 1.0075
       yspec(34) = 1.0154
       yspec(35) = 1.0232
       yspec(36) = 1.0116
       yspec(37) = 1.0000
       yspec(38) = 0.9887
       yspec(39) = 0.9774
       yspec(40) = 0.9833
       yspec(41) = 0.9892
       yspec(42) = 0.9621
       yspec(43) = 0.9350
       yspec(44) = 0.9559
       yspec(45) = 0.9769
       yspec(46) = 0.9848
       yspec(47) = 0.9927
       yspec(48) = 0.9916
       yspec(49) = 0.9904
       yspec(50) = 0.9738
       yspec(51) = 0.9572
       yspec(52) = 0.9729
       yspec(53) = 0.9886
       yspec(54) = 0.9726
       yspec(55) = 0.9567
       yspec(56) = 0.9693
       yspec(57) = 0.9819
       yspec(58) = 1.0060
       yspec(59) = 1.0300
       yspec(60) = 1.0107
       yspec(61) = 0.9913
       yspec(62) = 0.9326
       yspec(63) = 0.8738
       yspec(64) = 0.8949
       yspec(65) = 0.9160
       yspec(66) = 0.9225
       yspec(67) = 0.9289
       yspec(68) = 0.8487
       yspec(69) = 0.7685
       yspec(70) = 0.8168
       yspec(71) = 0.8651
       yspec(72) = 0.8955
       yspec(73) = 0.9258
       yspec(74) = 0.8540
       yspec(75) = 0.7823
       yspec(76) = 0.6796
       yspec(77) = 0.5769
       yspec(78) = 0.7031
       yspec(79) = 0.8292
       yspec(80) = 0.8060
       yspec(81) = 0.7827
!
       return       
       end subroutine illuminantD50
!
!======================================================================!
!
       subroutine illuminantD55(xspec,yspec,Xw,Yw,Zw)
!
       implicit none    
!
! Input/Output variables
!
       real(kind=4),dimension(81),intent(out)  ::  xspec     !
       real(kind=4),dimension(81),intent(out)  ::  yspec     !
       real(kind=4),intent(out)                ::  Xw,Yw,Zw  !  White point coordinates
!
! Illuminant white point tristimulus values
!
       Xw = 0.95682
       Yw = 1.00000 
       Zw = 0.92149
!
! Relative Spectral Power Distributions of CIE Standard Illuminant D55
!  at 5-nm Intervals from 380 to 780 nm
!
       call setlamb(xspec)
!
       yspec(1)  = 0.3258
       yspec(2)  = 0.3534
       yspec(3)  = 0.3809
       yspec(4)  = 0.4952
       yspec(5)  = 0.6095
       yspec(6)  = 0.6475
       yspec(7)  = 0.6855
       yspec(8)  = 0.7007
       yspec(9)  = 0.7158
       yspec(10) = 0.6975
       yspec(11) = 0.6791
       yspec(12) = 0.7676
       yspec(13) = 0.8561
       yspec(14) = 0.9180
       yspec(15) = 0.9799
       yspec(16) = 0.9923
       yspec(17) = 1.0046
       yspec(18) = 1.0019
       yspec(19) = 0.9991
       yspec(20) = 1.0133
       yspec(21) = 1.0274
       yspec(22) = 1.0041
       yspec(23) = 0.9808
       yspec(24) = 0.9938
       yspec(25) = 1.0068
       yspec(26) = 1.0069
       yspec(27) = 1.0070
       yspec(28) = 1.0034
       yspec(29) = 0.9999
       yspec(30) = 1.0210
       yspec(31) = 1.0421
       yspec(32) = 1.0316
       yspec(33) = 1.0210
       yspec(34) = 1.0253
       yspec(35) = 1.0297
       yspec(36) = 1.0148
       yspec(37) = 1.0000
       yspec(38) = 0.9861
       yspec(39) = 0.9722
       yspec(40) = 0.9748
       yspec(41) = 0.9775
       yspec(42) = 0.9459
       yspec(43) = 0.9143
       yspec(44) = 0.9293
       yspec(45) = 0.9442
       yspec(46) = 0.9478
       yspec(47) = 0.9514
       yspec(48) = 0.9468
       yspec(49) = 0.9422
       yspec(50) = 0.9233
       yspec(51) = 0.9045
       yspec(52) = 0.9139
       yspec(53) = 0.9233
       yspec(54) = 0.9059
       yspec(55) = 0.8885
       yspec(56) = 0.8959
       yspec(57) = 0.9032
       yspec(58) = 0.9213
       yspec(59) = 0.9395
       yspec(60) = 0.9195
       yspec(61) = 0.8996
       yspec(62) = 0.8482
       yspec(63) = 0.7968
       yspec(64) = 0.8126
       yspec(65) = 0.8284
       yspec(66) = 0.8384
       yspec(67) = 0.8484
       yspec(68) = 0.7754
       yspec(69) = 0.7024
       yspec(70) = 0.7477
       yspec(71) = 0.7930
       yspec(72) = 0.8215
       yspec(73) = 0.8499
       yspec(74) = 0.7844
       yspec(75) = 0.7188
       yspec(76) = 0.6234
       yspec(77) = 0.5279
       yspec(78) = 0.6436
       yspec(79) = 0.7593
       yspec(80) = 0.7387
       yspec(81) = 0.7182
!
       return       
       end subroutine illuminantD55
!
!======================================================================!
!
       subroutine illuminantD65(xspec,yspec,Xw,Yw,Zw)
!
       implicit none    
!
! Input/Output variables
!
       real(kind=4),dimension(81),intent(out)  ::  xspec     !
       real(kind=4),dimension(81),intent(out)  ::  yspec     !
       real(kind=4),intent(out)                ::  Xw,Yw,Zw  !  White point coordinates
!
! Illuminant white point tristimulus values
!
       Xw = 0.95047
       Yw = 1.00000 
       Zw = 1.08883
!
! Relative Spectral Power Distributions of CIE Standard Illuminant D65
!  at 5-nm Intervals from 380 to 780 nm
!
       call setlamb(xspec)
!
       yspec(1)  = 0.4998
       yspec(2)  = 0.5231
       yspec(3)  = 0.5465
       yspec(4)  = 0.6870
       yspec(5)  = 0.8275
       yspec(6)  = 0.8712
       yspec(7)  = 0.9149
       yspec(8)  = 0.9246
       yspec(9)  = 0.9343
       yspec(10) = 0.9006
       yspec(11) = 0.8668
       yspec(12) = 0.9577
       yspec(13) = 1.0486
       yspec(14) = 1.1094
       yspec(15) = 1.1701
       yspec(16) = 1.1741
       yspec(17) = 1.1781
       yspec(18) = 1.1634
       yspec(19) = 1.1486
       yspec(20) = 1.1539
       yspec(21) = 1.1592
       yspec(22) = 1.1237
       yspec(23) = 1.0881
       yspec(24) = 1.0908
       yspec(25) = 1.0935
       yspec(26) = 1.0858
       yspec(27) = 1.0780
       yspec(28) = 1.0630
       yspec(29) = 1.0479
       yspec(30) = 1.0624
       yspec(31) = 1.0769
       yspec(32) = 1.0605
       yspec(33) = 1.0441
       yspec(34) = 1.0423
       yspec(35) = 1.0405
       yspec(36) = 1.0202
       yspec(37) = 1.0000
       yspec(38) = 0.9817
       yspec(39) = 0.9633
       yspec(40) = 0.9606
       yspec(41) = 0.9579
       yspec(42) = 0.9224
       yspec(43) = 0.8869
       yspec(44) = 0.8935
       yspec(45) = 0.9001
       yspec(46) = 0.8980
       yspec(47) = 0.8960
       yspec(48) = 0.8865
       yspec(49) = 0.8770
       yspec(50) = 0.8549
       yspec(51) = 0.8329
       yspec(52) = 0.8349
       yspec(53) = 0.8370
       yspec(54) = 0.8186
       yspec(55) = 0.8003
       yspec(56) = 0.8012
       yspec(57) = 0.8021
       yspec(58) = 0.8125
       yspec(59) = 0.8228
       yspec(60) = 0.8028
       yspec(61) = 0.7828
       yspec(62) = 0.7400
       yspec(63) = 0.6972
       yspec(64) = 0.7067
       yspec(65) = 0.7161
       yspec(66) = 0.7298
       yspec(67) = 0.7435
       yspec(68) = 0.6798
       yspec(69) = 0.6160
       yspec(70) = 0.6574
       yspec(71) = 0.6989
       yspec(72) = 0.7249
       yspec(73) = 0.7509
       yspec(74) = 0.6934
       yspec(75) = 0.6359
       yspec(76) = 0.5501
       yspec(77) = 0.4642
       yspec(78) = 0.5661
       yspec(79) = 0.6681
       yspec(80) = 0.6509
       yspec(81) = 0.6338
!
       return       
       end subroutine illuminantD65
!
!======================================================================!
!
       subroutine illuminantD75(xspec,yspec,Xw,Yw,Zw)
!
       implicit none    
!
! Input/Output variables
!
       real(kind=4),dimension(81),intent(out)  ::  xspec     !
       real(kind=4),dimension(81),intent(out)  ::  yspec     !
       real(kind=4),intent(out)                ::  Xw,Yw,Zw  !  White point coordinates
!
! Illuminant white point tristimulus values
!
       Xw = 0.94972
       Yw = 1.00000 
       Zw = 1.22638
!
! Relative Spectral Power Distributions of CIE Standard Illuminant D75
!  at 5-nm Intervals from 380 to 780 nm
!
       call setlamb(xspec)
!
       yspec(1)  = 0.6670
       yspec(2)  = 0.6833
       yspec(3)  = 0.6996
       yspec(4)  = 0.8595
       yspec(5)  = 1.0193
       yspec(6)  = 1.0691
       yspec(7)  = 1.1189
       yspec(8)  = 1.1235
       yspec(9)  = 1.1280
       yspec(10) = 1.0794
       yspec(11) = 1.0309
       yspec(12) = 1.1214
       yspec(13) = 1.2120
       yspec(14) = 1.2710
       yspec(15) = 1.3301
       yspec(16) = 1.3268
       yspec(17) = 1.3236
       yspec(18) = 1.2984
       yspec(19) = 1.2732
       yspec(20) = 1.2706
       yspec(21) = 1.2680
       yspec(22) = 1.2229
       yspec(23) = 1.1778
       yspec(24) = 1.1719
       yspec(25) = 1.1659
       yspec(26) = 1.1515
       yspec(27) = 1.1370
       yspec(28) = 1.1118
       yspec(29) = 1.0856
       yspec(30) = 1.0955
       yspec(31) = 1.1044
       yspec(32) = 1.0837
       yspec(33) = 1.0629
       yspec(34) = 1.0560
       yspec(35) = 1.0490
       yspec(36) = 1.0245
       yspec(37) = 1.0000
       yspec(38) = 0.9781
       yspec(39) = 0.9562
       yspec(40) = 0.9491
       yspec(41) = 0.9421
       yspec(42) = 0.9060
       yspec(43) = 0.8700
       yspec(44) = 0.8711
       yspec(45) = 0.8723
       yspec(46) = 0.8668
       yspec(47) = 0.8614
       yspec(48) = 0.8486
       yspec(49) = 0.8358
       yspec(50) = 0.8116
       yspec(51) = 0.7875
       yspec(52) = 0.7859
       yspec(53) = 0.7843
       yspec(54) = 0.7661
       yspec(55) = 0.7480
       yspec(56) = 0.7456
       yspec(57) = 0.7432
       yspec(58) = 0.7487
       yspec(59) = 0.7542
       yspec(60) = 0.7350
       yspec(61) = 0.7158
       yspec(62) = 0.6771
       yspec(63) = 0.6385
       yspec(64) = 0.6446
       yspec(65) = 0.6508
       yspec(66) = 0.6657
       yspec(67) = 0.6807
       yspec(68) = 0.6226
       yspec(69) = 0.5644
       yspec(70) = 0.6034
       yspec(71) = 0.6424
       yspec(72) = 0.6670
       yspec(73) = 0.6915
       yspec(74) = 0.6389
       yspec(75) = 0.5863
       yspec(76) = 0.5062
       yspec(77) = 0.4262
       yspec(78) = 0.5198
       yspec(79) = 0.6135
       yspec(80) = 0.5984
       yspec(81) = 0.5832
!
       return       
       end subroutine illuminantD75
!
!======================================================================!
!
       subroutine illuminantE(xspec,yspec,Xw,Yw,Zw)
!
       implicit none    
!
! Input/Output variables
!
       real(kind=4),dimension(81),intent(out)  ::  xspec     !
       real(kind=4),dimension(81),intent(out)  ::  yspec     !
       real(kind=4),intent(out)                ::  Xw,Yw,Zw  !  White point coordinates
!
! Illuminant white point tristimulus values
!
       Xw = 1.00000
       Yw = 1.00000 
       Zw = 1.00000
!
! Relative Spectral Power Distributions of CIE Standard Illuminant D75
!  at 5-nm Intervals from 380 to 780 nm
!
       call setlamb(xspec)
!
       yspec(:) = 1.0d0
!
       return       
       end subroutine illuminantE
!
!======================================================================!
!
       subroutine setlamb(xspec)
!
       implicit none    
!
! Input/Output variables
!
       real(kind=4),dimension(81),intent(out)  ::  xspec  !
!
! Table of wavelengths at 5-nm Intervals from 380 to 780 nm
!
       xspec(1)  = 380.0
       xspec(2)  = 385.0
       xspec(3)  = 390.0
       xspec(4)  = 395.0
       xspec(5)  = 400.0
       xspec(6)  = 405.0
       xspec(7)  = 410.0
       xspec(8)  = 415.0
       xspec(9)  = 420.0
       xspec(10) = 425.0
       xspec(11) = 430.0
       xspec(12) = 435.0
       xspec(13) = 440.0
       xspec(14) = 445.0
       xspec(15) = 450.0
       xspec(16) = 455.0
       xspec(17) = 460.0
       xspec(18) = 465.0
       xspec(19) = 470.0
       xspec(20) = 475.0
       xspec(21) = 480.0
       xspec(22) = 485.0
       xspec(23) = 490.0
       xspec(24) = 495.0
       xspec(25) = 500.0
       xspec(26) = 505.0
       xspec(27) = 510.0
       xspec(28) = 515.0
       xspec(29) = 520.0
       xspec(30) = 525.0
       xspec(31) = 530.0
       xspec(32) = 535.0
       xspec(33) = 540.0
       xspec(34) = 545.0
       xspec(35) = 550.0
       xspec(36) = 555.0
       xspec(37) = 560.0
       xspec(38) = 565.0
       xspec(39) = 570.0
       xspec(40) = 575.0
       xspec(41) = 580.0
       xspec(42) = 585.0
       xspec(43) = 590.0
       xspec(44) = 595.0
       xspec(45) = 600.0
       xspec(46) = 605.0
       xspec(47) = 610.0
       xspec(48) = 615.0
       xspec(49) = 620.0
       xspec(50) = 625.0
       xspec(51) = 630.0
       xspec(52) = 635.0
       xspec(53) = 640.0
       xspec(54) = 645.0
       xspec(55) = 650.0
       xspec(56) = 655.0
       xspec(57) = 660.0
       xspec(58) = 665.0
       xspec(59) = 670.0
       xspec(60) = 675.0
       xspec(61) = 680.0
       xspec(62) = 685.0
       xspec(63) = 690.0
       xspec(64) = 695.0
       xspec(65) = 700.0
       xspec(66) = 705.0
       xspec(67) = 710.0
       xspec(68) = 715.0
       xspec(69) = 720.0
       xspec(70) = 725.0
       xspec(71) = 730.0
       xspec(72) = 735.0
       xspec(73) = 740.0
       xspec(74) = 745.0
       xspec(75) = 750.0
       xspec(76) = 755.0
       xspec(77) = 760.0
       xspec(78) = 765.0
       xspec(79) = 770.0
       xspec(80) = 775.0
       xspec(81) = 780.0
!
       return       
       end subroutine setlamb
!
!======================================================================!
!              
       end module illuminants
!
!======================================================================!
