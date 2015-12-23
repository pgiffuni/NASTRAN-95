BLOCK DATA ofp5bd
!OFP5BD
 INTEGER :: esingl, e1, e21, e41, e61, e81
 COMMON /ofpbd5/ esingl(64),e1(100),e21(100),e41(100),e61(100) ,e81(100)
!*****
!     SPACING ARRAY - ESINGL
!*****
 DATA    esingl( 1) / 4H/    /
 DATA    esingl( 2) / 4H15X  /
 DATA    esingl( 3) / 4H10X  /
 DATA    esingl( 4) / 4H5X   /
 DATA    esingl( 5) / 4H1X   /
 DATA    esingl( 6) / 4H/10X /
 DATA    esingl( 7) / 4H16X  /
 DATA    esingl( 8) / 4H2H1  /
 DATA    esingl( 9) / 4H2H2  /
 DATA    esingl(10) / 4H2H3  /
 DATA    esingl(11) / 4H2H4  /
 DATA    esingl(12) / 4H2H5  /
 DATA    esingl(13) / 4H7X   /
 DATA    esingl(14) / 4H/16X /
 DATA    esingl(15) / 4H/13X /
 DATA    esingl(16) / 4H4X   /
 DATA    esingl(17) / 4H/14X /
 DATA    esingl(18) / 4H11X  /
 DATA    esingl(19) / 4H/24X /
 DATA    esingl(20) / 4H1H0  /
 DATA    esingl(21) / 4H2H / /
 DATA    esingl(22) / 4H2HEN /
 DATA    esingl(23) / 4H2HDA /
 DATA    esingl(24) / 4H2HDB /
 DATA    esingl(25) / 4H/1H0 /
 DATA    esingl(26) / 4H23X  /
 DATA    esingl(27) / 4H/26X /
 DATA    esingl(28) / 4H/9X  /
 DATA    esingl(29) / 4H/12X /
 DATA    esingl(30) / 4H/1H  /
 DATA    esingl(31) / 4H/20X /
 DATA    esingl(32) / 4H/32X /
 DATA    esingl(33) / 4H/28X /
 DATA    esingl(34) / 4H/15X /
 DATA    esingl(35) / 4H/19X /
 DATA    esingl(36) / 4H/21X /
 DATA    esingl(37) / 4H/11X /
 DATA    esingl(38) / 4H/17X /
 DATA    esingl(39) / 4H2X   /
 DATA    esingl(40) / 4H1HX  /
 DATA    esingl(41) / 4H2HXY /
 DATA    esingl(42) / 4H1HA  /
 DATA    esingl(43) / 4H2HLX /
 DATA    esingl(44) / 4H1HY  /
 DATA    esingl(45) / 4H2HYZ /
 DATA    esingl(46) / 4H1HB  /
 DATA    esingl(47) / 4H2HLY /
 DATA    esingl(48) / 4H1HZ  /
 DATA    esingl(49) / 4H2HZX /
 DATA    esingl(50) / 4H1HC  /
 DATA    esingl(51) / 4H2HLZ /
 DATA    esingl(52) / 4H2HCP /
 DATA    esingl(53) / 4H2HMP /
 DATA    esingl(54) / 4H2HC  /
 DATA    esingl(55) / 4H3X   /
 DATA    esingl(56) / 4H/30X /
 DATA    esingl(57) / 4H9X   /
 DATA    esingl(58) / 4H/23X /
 DATA    esingl(59) / 4H6X   /
 DATA    esingl(60) / 4H39X  /
 DATA    esingl(61) / 4H24X  /
 DATA    esingl(62) / 4H     /
 DATA    esingl(63) / 4H     /
 DATA    esingl(64) / 4H     /
 
 
!     FORMAT BUILDING BLOCK E-ARRAY
 
 
!                    -STANDARD-                 -ALTERNATES-
!                 ****************        ***********************
 DATA    e1  / 4H1P,E ,4H15.6        , 4H0P,F ,4H6.1  ,4H,9X  &
     , 4H1P,E ,4H16.6        , 4H0P,F ,4H7.1  ,4H,9X  &
     , 4H1P,E ,4H17.6        , 4H0P,F ,4H8.1  ,4H,9X  &
     , 4H1P,E ,4H18.6        , 4H0P,F ,4H9.1  ,4H,9X  &
     , 4H1P,E ,4H19.6        , 4H0P,F ,4H10.1 ,4H,9X  &
     , 4H1P,E ,4H20.6        , 4H0P,F ,4H11.1 ,4H,9X  &
     , 4H1P,E ,4H21.6        , 4H0P,F ,4H12.1 ,4H,9X  &
     , 4H1P,E ,4H30.6        , 4H0P,F ,4H21.1 ,4H,9X  &
     , 4H1P,E ,4H26.6        , 4H0P,F ,4H17.1 ,4H,9X  &
     , 4H1P,E ,4H24.6        , 4H0P,F ,4H15.1 ,4H,9X  &
     , 4H0P,F ,4H11.4        , 4H0P,F ,4H8.1  ,4H,3X  &
     , 4H0P,F ,4H14.4        , 4H0P,F ,4H11.1 ,4H,3X  &
     , 4H1P,E ,4H28.6        , 4H0P,F ,4H19.1 ,4H,9X  &
     , 4H1P,E ,4H37.6        , 4H0P,F ,4H28.1 ,4H,9X  &
     , 4H1P,E ,4H22.6        , 4H0P,F ,4H17.1 ,4H,5X  &
     , 4H1P,E ,4H14.6        , 4H0P,F ,4H5.1  ,4H,9X  &
     , 4H0P,F ,4H15.4        , 4H0P,F ,4H12.1 ,4H,3X  &
     , 4H0P,F ,4H9.4         , 4H0P,F ,4H6.1  ,4H,3X  &
     , 4H0P,F ,4H15.3        , 4H0P,F ,4H12.1 ,4H 3X  &
     , 4H1P,E ,4H23.6        , 4H0P,F ,4H14.1 ,4H,9X    /

 DATA    e21 / 4H1P,E ,4H35.6        , 4H0P,F ,4H26.1 ,4H,9X  &
     , 4H1P,E ,4H25.6        , 4H0P,F ,4H16.1 ,4H,9X  &
     , 4H1P,E ,4H50.6        , 4H0P,F ,4H41.1 ,4H,9X  &
     , 4H0P,F ,4H46.4        , 4H0P,F ,4H43.1 ,4H,3X  &
     , 4H     ,4H            , 4H0P,F ,4H12.1 ,4H,3X  &
     , 4H0P,F ,4H20.4        , 4H0P,F ,4H17.1 ,4H,3X  &
     , 4H0P,F ,4H16.4        , 4H0P,F ,4H13.1 ,4H,3X  &
     , 4H0P,F ,4H22.4        , 4H0P,F ,4H19.1 ,4H,3X  &
     , 4H1P,E ,4H27.6        , 4H0P,F ,4H18.1 ,4H,9X  &
     , 4H0P,F ,4H12.5        , 4H0P,F ,4H11.1 ,4H,3X  &
     , 4H1P,E ,4H13.5        , 4H0P,F ,4H5.1  ,4H,8X  &
     , 4H0P,F ,4H13.3        , 4H0P,F ,4H9.1  ,4H,4X  &
     , 4H0P,F ,4H18.4        , 4H0P,F ,4H15.1 ,4H,3X  &
     , 4H0P,F ,4H26.4        , 4H0P,F ,4H23.1 ,4H,3X  &
     , 4H1P,E ,4H14.5        , 4H0P,F ,4H6.1  ,4H,8X  &
     , 4H0P,F ,4H14.3        , 4H0P,F ,4H10.1 ,4H,4X  &
     , 4H0P,F ,4H5.2         , 4H0P,F ,4H4.1  ,4H,1X  &
     , 4H1P,E ,4H13.6        , 4H0P,F ,4H4.1  ,4H,9X  &
     , 4H     ,4H            , 4H     ,4H     ,4H     &
     , 4H1P,E ,4H9.1         , 4HA1   ,4H,8X  ,4H      /

 DATA    e41 / 4H6X,A ,4H1,3X        , 4HI7   ,4H,3X  ,4H     &
     , 4HI15  ,4H            , 4H     ,4H     ,4H     &
     , 4HI9,1 ,4HX           , 4H     ,4H     ,4H     &
     , 4H1H0, ,4HI8          , 4H     ,4H     ,4H     &
     , 4H1X,I ,4H13          , 4H     ,4H     ,4H     &
     , 4H1X,I ,4H8           , 4H     ,4H     ,4H     &
     , 4H1H0, ,4HI7          , 4H     ,4H     ,4H     &
     , 4H6X,I ,4H8           , 4H     ,4H     ,4H     &
     , 4H1X,I ,4H15          , 4H     ,4H     ,4H     &
     , 4H1X,I ,4H12          , 4H     ,4H     ,4H     &
     , 4HI10  ,4H            , 4H     ,4H     ,4H     &
     , 4HI7,1 ,4HX           , 4H     ,4H     ,4H     &
     , 4H3X,A ,4H4           , 4H     ,4H     ,4H     &
     , 4H1H0, ,4HI13         , 4H     ,4H     ,4H     &
     , 4H1X,I ,4H20          , 4H     ,4H     ,4H     &
     , 4H5X,A ,4H1,3X        , 4HI5   ,4H,4X  ,4H     &
     , 4H1X,I ,4H22          , 4H     ,4H     ,4H     &
     , 4HI12  ,4H            , 4H     ,4H     ,4H     &
     , 4H1X,I ,4H19          , 4H     ,4H     ,4H     &
     , 4HI16  ,4H            , 4H     ,4H     ,4H      /
 
 DATA    e61 / 4HI8   ,4H            , 4HA4   ,4H,4X  ,4H    &
     , 4HI9   ,4H            , 4HA4   ,4H,5X  ,4H     &
     , 4HI11  ,4H            , 4HA4   ,4H,7X  ,4H     &
     , 4HI20  ,4H            , 4HA4   ,4H,16X ,4H     &
     , 4HI19  ,4H            , 4HA4   ,4H,15X ,4H     &
     , 4H1X,I ,4H23          , 4H     ,4H     ,4H     &
     , 4HI23  ,4H            , 4H     ,4H     ,4H     &
     , 4HI28  ,4H            , 4H     ,4H     ,4H     &
     , 4H/1H  ,4H,I18        , 4H     ,4H     ,4H     &
     , 4H1H0, ,4HI15         , 4H     ,4H     ,4H     &
     , 4H1H0, ,4HI14         , 4H     ,4H     ,4H     &
     , 4H0P,F ,4H22.4        , 4HI9   ,4H,13X ,4H     &
     , 4H0P,F ,4H16.4        , 4HI5   ,4H,11X ,4H     &
     , 4H0P,F ,4H10.4        , 4H     ,4H     ,4H     &
     , 4H1H0, ,4HI19         , 4H     ,4H     ,4H     &
     , 4H1H0, ,4HI20         , 4H     ,4H     ,4H     &
     , 4HI10, ,4H5X          , 4H     ,4H     ,4H     &
     , 4H     ,4H            , 4H     ,4H     ,4H     &
     , 4HI8,  ,4H2X          , 4H3X,3 ,4HHCEN ,4H,A4  &
     , 4HF8.3 ,4H            , 4H     ,4H     ,4H      /

 DATA    e81 / 4H1H0, ,4HI27         , 4H     ,4H     ,4H    &
     , 4H1H0, ,4HI5          , 4H     ,4H     ,4H     &
     , 4H1H0, ,4HI3          , 4H     ,4H     ,4H     &
     , 4HI4   ,4H            , 4H     ,4H     ,4H     &
     , 4H1P,E ,4H11.4        , 4H0P,F ,4H4.1  ,4H,7X  &
     , 4HA4   ,4H            , 4H     ,4H     ,4H     &
     , 4H1P9E ,4H11.3        , 4H0P9( ,4HF9.3 ,4H,2X) &
!AIX , 4H  9E ,4H11.3        , 4H  9( ,4HF9.3 ,4H,2X) &
     , 4H0P,F ,4H22.3        , 4H0P,F ,4H20.1 ,4H,2X  &
     , 4H/1PE ,4H11.3        , 4H/0PF ,4H7.1  ,4H,4X  &
!AIX , 4H/, E ,4H11.3        , 4H/0P, ,4HF7.1 ,4H,4X  &
     , 4H0P,F ,4H19.4        , 4HI6   ,4H,13X ,4H     &
     , 4HF8.2 ,4H            , 4H     ,4H     ,4H     &
     , 4H1P,E ,4H12.5        , 4H     ,4H     ,4H     &
     , 4H1H0, ,4HI12         , 4H     ,4H     ,4H     &
     , 4H4X,I ,4H8           , 4H4X,A ,4H4,4X ,4H     &
     , 4H     ,4H            , 4H     ,4H     ,4H     &
     , 4H     ,4H            , 4H     ,4H     ,4H     &
     , 4H     ,4H            , 4H     ,4H     ,4H     &
     , 4H     ,4H            , 4H     ,4H     ,4H     &
     , 4H     ,4H            , 4H     ,4H     ,4H     &
     , 4H     ,4H            , 4H     ,4H     ,4H      /

END BLOCK DATA
