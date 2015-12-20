SUBROUTINE ifp1pc (i81,icont,pocard,org,porg)
     
!     SUBROUTINE TO PERFORM FIRST-LEVEL CHECKING OF STRUCTURE PLOTTER
!     CONTROL CARD FORMAT.
 
 
 INTEGER, INTENT(IN)                      :: i81
 INTEGER, INTENT(OUT)                     :: icont
 INTEGER, INTENT(IN)                      :: pocard(1)
 INTEGER, INTENT(OUT)                     :: org
 INTEGER, INTENT(OUT)                     :: porg
 EXTERNAL        rshift,complf
 LOGICAL :: flag(3),bit64
 INTEGER :: case(400),ctype(21),idvpr(3),camera(5),origin(11),  &
     axes(3),maxes(3),cntur(20),setpr(33),setp2(12),  &
     coord(25),lblpr(5),pltpr(28),nast(2), core(1),corey(401), &
     blank
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm
 COMMON /system/ isys,nout,nogo,skp(16),pltopt,sys21,ilink, skp63(63),intra
 COMMON /xifp1 / BLANK,bit64
 COMMON /zzzzzz/ corex(1)
 EQUIVALENCE     (proj,ctype(11)), (defo,idvpr( 1)),  &
     (symm,pltpr(13)), (anti,pltpr(14)), (magn,cntur(13)), (thru,pltpr(22)),  &
     (poin,lblpr( 2)), (core(1),corey(401)),  &
     (corex(1),corey(1),case(1)), (hidd,pltpr(24))
 DATA    ctype / 4HPLOT, 4HORTH, 4HPERS, 4HSTER, 4HAXES, 4HVIEW,  &
     4HMAXI, 4HCSCA, 4HFIND, 4HCONT, 4HPROJ, 4HOCUL,  &
     4HCAME, 4HPAPE, 4HPEN , 4HPTIT, 4HSCAL, 4HORIG, 4HVANT, 4HSET , 4HREGI/
 DATA    camera/ 4HFILM, 4HPAPE, 4HBOTH, 4HBLAN, 4HFRAM/
 DATA    axes  / 4HX   , 4HY   , 4HZ   /
 DATA    maxes / 4HMX  , 4HMY  , 4HMZ  /
 DATA    cntur / 4HMAJP, 4HMINP, 4HMAXS, 4HXNOR, 4HYNOR, 4HZNOR,  &
     4HXYSH, 4HXZSH, 4HYZSH, 4HXDIS, 4HYDIS, 4HZDIS,  &
     4HMAGN, 4HNRM1, 4HNRM2, 4HSH12, 4HSH1Z, 4HSH2Z, 4HBDSH, 4HSTRA/
 DATA    setpr / 4HINCL, 4HEXCL, 4HEXCE, 4HELEM, 4HGRID, 4HALL ,  &
     4HAERO, 4HAXIF, 4HBAR , 4HCONE, 4HCONR, 4HHEXA,  &
     4HFLUI, 4HIHEX, 4HPLOT, 4HQDME, 4HQDPL, 4HQUAD,  &
     4HROD , 4HSHEA, 4HSLOT, 4HTETR, 4HTORD, 4HTRAP,  &
     4HTRBS, 4HTRIA, 4HTRME, 4HTRPL, 4HTUBE, 4HTWIS, 4HVISC, 4HWEDG, 4HHBDY/
 DATA    setp2 / 4HAX  , 4HRG  , 4H1   , 4H2   , 4H3   , 4H4   ,  &
     4HD2  , 4HD3  , 4HD4  , 4HM   , 4HM1  , 4HM2  /
 DATA    pltpr / 4HSET , 4HSTAT, 4HMODA, 4HCMOD, 4HFREQ, 4HTRAN,  &
     4HCONT, 4HRANG, 4HTIME, 4HPHAS, 4HMAGN, 4HORIG,  &
     4HSYMM, 4HANTI, 4HPEN , 4HDENS, 4HSYMB, 4HLABE,  &
     4HSHAP, 4HVECT, 4HOUTL, 4HTHRU, 4HMAXI, 4HHIDD,  &
     4HSHRI, 4HNOFI, 4HFILL, 4HOFFS/
 DATA    idvpr / 4HDEFO, 4HVELO, 4HACCE/
 DATA    coord / 4HYX  , 4HZX  , 4HZY  , 4HXY  , 4HXZ  , 4HYZ  ,  &
     4HX   , 4HY   , 4HZ   , 4HXYZ , 4HRXY , 4HRXZ , 4HRYZ , 4HR   , 4HRN  ,  &
     4HXN  , 4HYN  , 4HZN  , 4HXYN , 4HXZN , 4HYZN ,  &
     4HXYZN, 4HRXYN, 4HRXZN, 4HRYZN  /
 DATA    lblpr / 4HGRID, 4HPOIN, 4HELEM, 4HBOTH, 4HEPID/
 DATA    ter   / 4HTER /, plan / 4HPLAN/, sepa / 4HSEPA/
 DATA    lag   / 4HLAG /, nast / 4HSC  , 4HCALC/,ilnk  / 4HNS01/
 
 
!     INITIALIZE
 
 IF (intra <= 1 .AND. ilink == ilnk) GO TO 15
 DO  i = 1,200
   core(i)= pocard(i)
 END DO
 15 allon  = complf(0)
 eor    = rshift(allon,1)
 isplot = 0
 iwrd   = i81
 
!     BRANCH FOR CONTINUATION CARD
!                                  SET   PLOT  FIND
 IF (icont /= 0) THEN
    SELECT CASE ( icont )
     CASE (    1)
       GO TO 10
     CASE (    2)
       GO TO  2111
     CASE (    3)
       GO TO  2210
     CASE (    4)
       GO TO  1067
   END SELECT
 END IF
 
 IF (core(iwrd) < 0.0) THEN
   GO TO  9800
 ELSE IF (core(iwrd) == 0.0) THEN
   GO TO   350
 ELSE
   GO TO    20
 END IF
 10 IF (core(iwrd) <=   0) GO TO 320
 20 IF (core(iwrd) == eor) GO TO 350
 mode = core(iwrd)
 iwrd = iwrd + 1
 
!     BRANCH FOR CARD TYPE
 
 100 iword = core(iwrd)
 IF (bit64) CALL mvbits (BLANK,0,32,iword,0)
 DO  i = 1,20
   IF (iword == ctype(i))  &
       GO TO (400,  500,  500,  500,  600,  700,  800,  900, 1000,  &
       1100, 1200, 1300, 1400,  320,  320,  320, 1800, 1900, 2000, 2100), i
   
!    1         PLOT  ORTH  PERS  STER  AXES  VIEW  MAXI  CSCA  FIND
!    2         CONT  PROJ  OCUL  CAME  PAPE   PEN  PTIT  SCAL  ORIG
!    3         VANT   SET
   
 END DO
 GO TO 9802
 320 IF (mode <= 0) GO TO 330
 iwrd = iwrd + 2
 mode = mode - 1
 GO TO 320
 330 IF (core(iwrd) < 0.0) THEN
   GO TO   335
 ELSE
   GO TO   340
 END IF
 335 IF (core(iwrd) == -4) iwrd = iwrd + 1
 iwrd = iwrd + 2
 GO TO 330
 340 IF (core(iwrd) == 0 .OR. core(iwrd) == eor) GO TO 350
 mode = core(iwrd)
 iwrd = iwrd + 1
 GO TO 320
 350 icont = 0
 IF (core(iwrd) == 0) icont = 1
 GO TO 9998
 
!     BRANCH TO PLOT OR PLOTTER
 
 400 iword = core(iwrd+1)
 IF (bit64) CALL mvbits (BLANK,0,32,iword,0)
 IF (iword == ter) GO TO 410
 isplot = 1
 GO TO 2200
 
!     PLOTTER CARD
 
 410 iword = core(iwrd+2)
 IF (bit64) CALL mvbits (BLANK,0,32,iword,0)
 IF (iword == nast(1) .OR. iword == nast(2)) GO TO 9804
 GO TO 320
 
!     PROJECTION CARD
 
 500 iwrd  = iwrd + 2
 mode  = mode - 1
 iword = core(iwrd)
 IF (bit64) CALL mvbits (BLANK,0,32,iword,0)
 IF (iword == proj) GO TO 510
 ASSIGN 510 TO irtn
 iprm = proj
 GO TO 9806
 510 iwrd = iwrd + 2
 mode = mode - 1
 IF (mode > 0) THEN
   GO TO   100
 ELSE
   GO TO   330
 END IF
 
!     AXES CARD
 
 600 iwrd = iwrd + 2
 mode = mode - 1
 iword = core(iwrd)
 IF (bit64) CALL mvbits (BLANK,0,32,iword,0)
 IF (core(iwrd) == allon .OR. iword == BLANK) GO TO 600
 DO  j = 1,3
   flag(j) = .false.
 END DO
 i = 0
 GO TO 607
 606 iword = core(iwrd)
 IF (bit64) CALL mvbits (BLANK,0,32,iword,0)
 607 IF (core(iwrd) == allon .OR. iword == BLANK) GO TO 610
 DO  j = 1,3
   IF (iword == axes(j) .OR. iword == maxes(j)) flag(j) = .true.
 END DO
 i = i + 1
 610 iwrd = iwrd + 2
 mode = mode - 1
 IF (i < 3) GO TO 606
 
 ASSIGN 320 TO irtn
 IF (.NOT.flag(1) .OR. .NOT.flag(2) .OR. .NOT.flag(3)) GO TO 9810
 620 iword = core(iwrd)
 IF (bit64) CALL mvbits (BLANK,0,32,iword,0)
 IF (iword == symm .OR. iword == anti) GO TO 630
 IF (core(iwrd) == 0 .OR. core(iwrd) == eor) GO TO 350
 IF (core(iwrd) /= allon .AND. iword /= BLANK) GO TO 100
 iwrd = iwrd + 2
 mode = mode - 1
 IF (mode > 0) THEN
   GO TO   620
 ELSE
   GO TO  9812
 END IF
 630 iwrd = iwrd + 2
 mode = mode - 1
 IF (mode > 0) THEN
   GO TO   100
 ELSE
   GO TO   330
 END IF
 
!     VIEW COMMAND
 
 700 nreal = 3
 nopt  = 0
 GO TO 1310
 
!     MAXIMUM DEFORMATION CARD
 
 800 nreal = 1
 nopt  = 0
 iword = core(iwrd+2)
 IF (bit64) CALL mvbits (BLANK,0,32,iword,0)
 IF (iword == defo) GO TO 1310
 ASSIGN 320 TO irtn
 iprm = core(iwrd+2)
 GO TO 9808
 
!     CSCALE CARD
 
 900 ASSIGN 320 TO irtn
 910 iwrd = iwrd + 2
 mode = mode - 1
 IF (mode > 0) THEN
   GO TO   920
 ELSE
   GO TO   930
 END IF
 920 iword = core(iwrd)
 IF (bit64) CALL mvbits (BLANK,0,32,iword,0)
 IF (core(iwrd) == allon .OR. iword == BLANK) GO TO 910
 GO TO 9812
 930 IF (core(iwrd)+1 < 0.0) THEN
   GO TO   960
 ELSE IF (core(iwrd)+1 == 0.0) THEN
   GO TO   940
 ELSE
   GO TO  9816
 END IF
 940 WRITE  (nout,950)
 950 FORMAT (/5X,'REAL VALUE, NOT INTEGER, IS NOW USED FOR CSCALE')
 GO TO 9816
 
 960 nreal = 1
 nopt  = 0
 GO TO 1700
 
!     FIND COMMAND
 
 1000 iwrd = iwrd + 2
 mode = mode - 1
 1005 IF (core(iwrd) == 0 .OR. core(iwrd) == eor) GO TO 1080
 ASSIGN 1070 TO irtn
 IF (mode > 0) THEN
   GO TO  1006
 ELSE
   GO TO  9812
 END IF
 1006 iword = core(iwrd)
 IF (bit64) CALL mvbits (BLANK,0,32,iword,0)
 IF (core(iwrd) == allon .OR. iword == BLANK) GO TO 1000
 DO  i = 17,21
   itype = i - 16
   IF (iword == ctype(i)) THEN 
     SELECT CASE ( itype )
       CASE (    1)
         GO TO 1020
       CASE (    2)
         GO TO  1030
       CASE (    3)
         GO TO  1040
       CASE (    4)
         GO TO  1030
       CASE (    5)
         GO TO  1050
     END SELECT
   ENDIF
!                SCAL  ORIG  VANT   SET  REGI
   
 END DO
 iprm = core(iwrd)
 GO TO 9808
 
 1020 nreal = 1
 1021 iwrd  = iwrd + 2
 mode  = mode - 1
 IF (mode <= 0) GO TO 1061
 iword = core(iwrd)
 IF (bit64) CALL mvbits (BLANK,0,32,iword,0)
 IF (core(iwrd) == allon .OR. iword == BLANK) GO TO 1021
 GO TO 1005
 
 1030 iprm = core(iwrd)
 ASSIGN 1005 TO irtn
 1031 iwrd = iwrd + 2
 mode = mode - 1
 IF (mode > 0) THEN
   GO TO  1032
 ELSE
   GO TO  1033
 END IF
 1032 iword = core(iwrd)
 IF (bit64) CALL mvbits (BLANK,0,32,iword,0)
 IF (core(iwrd) == allon .OR. iword == BLANK) GO TO 1031
 GO TO 9814
 1033 integ = 1
 IF (core(iwrd) == eor) GO TO 9814
 IF (core(iwrd) ==  -1) integ = 0
 IF (core(iwrd) ==  -4) iwrd  = iwrd + 1
 IF (itype /= 2) GO TO 1034
 forg = core(iwrd+1)
 org  = org + 1
 origin(org) = forg
 1034 iwrd = iwrd + 2
 IF (porg >= 0) GO TO 1066
 porg  = 0
 porg1 = forg
 GO TO 1066
 
 1040 iwrd = iwrd + 2
 mode = mode - 1
 ASSIGN 1070 TO irtn
 IF (mode > 0) THEN
   GO TO  1042
 END IF
 1041 iprm = poin
 GO TO 9806
 1042 iword = core(iwrd)
 IF (bit64) CALL mvbits (BLANK,0,32,iword,0)
 IF (iword == poin) GO TO 1000
 iprm = core(iwrd)
 GO TO 9808
 
 1050 nreal = 4
 1060 iwrd  = iwrd + 2
 mode  = mode - 1
 IF (mode <= 0) GO TO 1062
 iword = core(iwrd)
 IF (bit64) CALL mvbits (BLANK,0,32,iword,0)
 IF (core(iwrd) == allon .OR. iword == BLANK) GO TO 1060
 GO TO 9818
 1061 IF (core(iwrd) == 0 .OR. core(iwrd) == eor) GO TO 1080
 1062 integ = 0
 ASSIGN 1005 TO irtn
 DO  i = 1,nreal
   IF (core(iwrd) == -1) integ = 1
   IF (core(iwrd) == -4) iwrd  = iwrd + 1
   iwrd = iwrd + 2
 END DO
 1066 IF (integ > 0) THEN
   GO TO  9816
 END IF
 1067 IF (core(iwrd) == 0 .OR. core(iwrd) == eor) GO TO 1080
 mode = core(iwrd)
 iwrd = iwrd + 1
 GO TO 1005
 
 1070 IF (core(iwrd) == 0 .OR. core(iwrd) == eor) GO TO 1080
 iwrd = iwrd + 1
 GO TO 1070
 1080 icont = 0
 IF (core(iwrd) == 0) icont = 4
 GO TO 9998
 
!     CONTOUR
 
 1100 iwrd = iwrd + 2
 mode = mode - 1
 ASSIGN 320 TO irtn
 1105 IF (core(iwrd) == 0 .OR. core(iwrd) == eor) GO TO 350
 iword = core(iwrd)
 IF (bit64) CALL mvbits (BLANK,0,32,iword,0)
 IF (core(iwrd) == allon .OR. iword == BLANK) GO TO 1110
 DO  i = 1,20
   IF (iword == cntur(i)) GO TO 320
 END DO
 iprm = core(iwrd)
 GO TO 9808
 1110 iwrd = iwrd + 2
 mode = mode - 1
 IF (mode > 0) THEN
   GO TO  1105
 ELSE
   GO TO  9812
 END IF
 
!     PROJECTION PLANE SEPARATION
 
 1200 iwrd = iwrd + 2
 mode = mode - 1
 ASSIGN 320 TO irtn
 IF (mode > 0) THEN
   GO TO  1220
 END IF
 1210 iprm = plan
 GO TO 9806
 1220 iword = core(iwrd)
 IF (bit64) CALL mvbits (BLANK,0,32,iword,0)
 IF (iword /= plan) GO TO 1231
 iword = core(iwrd+2)
 IF (bit64) CALL mvbits (BLANK,0,32,iword,0)
 IF (iword == sepa) GO TO 1240
 1231 iprm = core(iwrd)
 GO TO 9808
 1240 nreal = 1
 nopt  = 0
 GO TO 1310
 
!     OCULAR SEPARATION
 
 1300 nreal = 1
 nopt  = 0
 iword = core(iwrd+2)
 IF (bit64) CALL mvbits (BLANK,0,32,iword,0)
 IF (iword == sepa) GO TO 1310
 ASSIGN 320 TO irtn
 iprm = core(iwrd+2)
 GO TO 9808
 
 1310 iwrd = iwrd + 2
 mode = mode - 1
 IF (mode > 0) THEN
   GO TO  1310
 ELSE
   GO TO  1700
 END IF
 
!     CAMERA
 
 1400 iwrd = iwrd + 2
 mode = mode - 1
 IF (mode <= 0) GO TO 1420
 iword = core(iwrd)
 IF (bit64) CALL mvbits (BLANK,0,32,iword,0)
 IF (core(iwrd) == allon .OR. iword == BLANK ) GO TO 1400
 IF (core(iwrd) == eor   .OR. core(iwrd) == 0) GO TO 9820
 DO  i = 1,4
   IF (iword == camera(i)) GO TO 1415
 END DO
 iprm = core(iwrd)
 ASSIGN 320 TO irtn
 GO TO 9808
 1415 iwrd = iwrd + 2
 mode = mode - 1
 IF (mode <= 0) GO TO 1420
 iword = core(iwrd)
 IF (bit64) CALL mvbits (BLANK,0,32,iword,0)
 IF (core(iwrd) == allon .OR. iword == BLANK) GO TO 1415
 i = i + 1
 IF (iword == camera(4) .OR. iword == camera(5)) GO TO 1415
 ASSIGN 320 TO irtn
 IF (i-4 == 0) THEN
   GO TO  9812
 ELSE
   GO TO   100
 END IF
 1420 IF (core(iwrd) == eor .OR. core(iwrd) == 0) IF (i-3) 350,350,9820
 ASSIGN 320 TO irtn
 IF (core(iwrd)+1 == 0.0) THEN
   GO TO  1430
 ELSE
   GO TO  9816
 END IF
 1430 iwrd = iwrd + 2
 GO TO 10
 
!     TEST FOR REAL VALUES
 
 1700 iro = 0
 nro = nreal
 1710 integ = 0
 ASSIGN 320 TO irtn
 DO  i = 1,nro
   IF (core(iwrd) >= 0 .OR. core(iwrd) < -4) IF (iro) 9818,9818,1712
   1712 IF (core(iwrd) == -1) integ = 1
   IF (core(iwrd) == -4) iwrd  = iwrd + 1
   iwrd = iwrd + 2
 END DO
 IF (integ == 0) GO TO 1730
 ASSIGN 1730 TO irtn
 GO TO 9816
 1730 IF (core(iwrd) < 0.0) THEN
   GO TO  1740
 ELSE IF (core(iwrd) == 0.0) THEN
   GO TO   350
 ELSE
   GO TO    20
 END IF
 1740 IF (iro == 1 .OR. nopt == 0) GO TO 9812
 iro = 1
 nro = nopt
 GO TO 1710
 
!     SCALE
 
 1800 nreal = 1
 nopt  = 1
 GO TO 1310
 
!     ORIGIN
 
 1900 nreal = 3
 nopt  = 0
 1905 iwrd  = iwrd + 2
 mode  = mode - 1
 IF (mode > 0) THEN
   GO TO  1906
 ELSE
   GO TO  1907
 END IF
 1906 iword = core(iwrd)
 IF (bit64) CALL mvbits (BLANK,0,32,iword,0)
 IF (core(iwrd) == allon .OR. iword == BLANK) GO TO 1905
 1907 IF (core(iwrd) == -1) GO TO 1910
 iprm = ctype(18)
 ASSIGN 320 TO irtn
 GO TO 9814
 1910 IF (core(iwrd) == -4) iwrd = iwrd + 1
 iwrd = iwrd + 2
 ASSIGN 320 TO irtn
 IF (core(iwrd) == eor) GO TO 9818
 IF (core(iwrd) < 0) GO TO 1700
 mode  = core(iwrd)
 iwrd  = iwrd + 1
 iword = core(iwrd)
 IF (bit64) CALL mvbits (BLANK,0,32,iword,0)
 IF (core(iwrd) == allon .OR. iword == BLANK) GO TO 1310
 GO TO 9812
 
!     VANTAGE POINT
 
 2000 nreal = 3
 nopt  = 1
 iword = core(iwrd+2)
 IF (bit64) CALL mvbits (BLANK,0,32,iword,0)
 IF (iword == poin) GO TO 1310
 ASSIGN 320 TO irtn
 iprm = core(iwrd+2)
 GO TO 9808
 
!     SET DEFINITION CARD
 
 2100 nint = 0
 nthru= 0
 2105 iwrd = iwrd + 2
 mode = mode - 1
 IF (mode > 0) THEN
   GO TO  2108
 END IF
 2106 IF (core(iwrd) == -1) GO TO 2110
 ASSIGN 2107 TO irtn
 GO TO 9816
 2107 IF (core(iwrd) == -4) iwrd = iwrd + 1
 GO TO 2110
 2108 iword = core(iwrd)
 IF (bit64) CALL mvbits (BLANK,0,32,iword,0)
 IF (core(iwrd) == allon .OR. iword == BLANK) GO TO 2105
 iprm = ctype(20)
 ASSIGN 2120 TO irtn
 GO TO 9814
 
 2110 iwrd  = iwrd + 2
 nreal = 0
 2111 IF (core(iwrd) < 0.0) THEN
   GO TO  2112
 ELSE IF (core(iwrd) == 0.0) THEN
   GO TO  2113
 ELSE
   GO TO  2114
 END IF
 2112 nint  = nint + 1
 IF (core(iwrd) == -1 .OR. nreal /= 0) GO TO 2110
 ASSIGN 2110 TO irtn
 GO TO 9816
 2113 icont = 2
 nthru = 0
 GO TO 9998
 2114 IF (core(iwrd) /= eor) GO TO 2115
 icont = 0
 GO TO 9998
 2115 mode  = core(iwrd)
 iwrd  = iwrd + 1
 2120 iword = core(iwrd)
 IF (bit64) CALL mvbits (BLANK,0,32,iword,0)
 IF (core(iwrd) /= allon .AND. iword /= BLANK) GO TO 2121
 iwrd  = iwrd + 2
 mode  = mode - 1
 IF (mode > 0) THEN
   GO TO  2120
 ELSE
   GO TO  2111
 END IF
 2121 iword = core(iwrd)
 IF (bit64) CALL mvbits (BLANK,0,32,iword,0)
 IF (iword /= thru) GO TO 2130
 nthru = nthru + 1
 IF (core(iwrd-3) == -1 .AND. core(iwrd+2) == -1) GO TO 2122
 ASSIGN 2123 TO irtn
 nreal = 1
 GO TO 9822
 2122 IF (nthru == 1) GO TO 2123
 IF (nint >= 2 .AND. core(iwrd-2) > core(iwrd-4)) GO TO 2123
 ASSIGN 2123 TO irtn
 GO TO 9824
 2123 nint = 0
 iwrd = iwrd + 2
 mode = mode - 1
 IF (mode > 0) THEN
   GO TO  2130
 ELSE
   GO TO  2111
 END IF
 2130 IF (core(iwrd) ==   0) GO TO 2113
 IF (core(iwrd) == eor) GO TO 2114
 iword = core(iwrd)
 IF (bit64) CALL mvbits (BLANK,0,32,iword,0)
 IF (core(iwrd) == allon .OR. iword == BLANK) GO TO 2135
 DO  i = 1,33
   IF (iword == setpr(i))  &
       GO TO (2135, 2135, 2135, 2135, 2138, 2135, 2135, 2142, 2135,  &
       2135, 2135, 2143, 2144, 2145, 2135, 2146, 2135, 2143,  &
       2135, 2135, 2147, 2135, 2135, 2148, 2135, 2149, 2135,  &
       2135, 2135, 2135, 2135, 2135, 2135), i
   
!    1          INCL  EXCL  EXCE  ELEM  GRID   ALL  AERO  AXIF   BAR
!    2          CONE  CONR  HEXA  FLUI  IHEX  PLOT  QDME  QDPL  QUAD
!    3           ROD  SHEA  SLOT  TETR  TORD  TRAP  TRBS  TRIA  TRME
!    4          TRPL  TUBE  TWIS VISCX  WEDG  HBDY
   
 END DO
 ASSIGN 2135 TO irtn
 iprm = core(iwrd)
 GO TO 9808
 2135 iwrd = iwrd + 2
 mode = mode - 1
 IF (mode <= 0) GO TO 2136
 iword = core(iwrd)
 IF (bit64) CALL mvbits (BLANK,0,32,iword,0)
 IF (core(iwrd) == allon .OR. iword == BLANK) GO TO 2135
 GO TO 2130
 2136 nthru = 0
 GO TO 2111
 
 2138 iwrd = iwrd + 2
 mode = mode - 1
 IF (mode > 0) THEN
   GO TO  2140
 END IF
 2139 ASSIGN 2136 TO irtn
 iprm = poin
 GO TO 9806
 2140 iword = core(iwrd)
 IF (bit64) CALL mvbits (BLANK,0,32,iword,0)
 IF (iword == poin) GO TO 2135
 ASSIGN 2130 TO irtn
 iprm = core(iwrd)
 GO TO 9808
 
 2142 istt = 4
 istb = 6
 GO TO 2150
 
 2143 istt = 3
 istb = 6
 GO TO 2150
 
 2144 istt = 7
 istb = 9
 GO TO 2150
 
 2145 istt = 3
 istb = 5
 GO TO 2150
 
 2146 istt = 10
 istb = 12
 GO TO 2150
 
 2147 istt = 5
 istb = 6
 GO TO 2150
 
 2148 istt = 1
 istb = 2
 GO TO 2150
 
 2149 istt = 1
 istb = 5
 
 2150 iword = core(iwrd+1)
 IF (bit64) CALL mvbits (BLANK,0,32,iword,0)
 DO  i = istt,istb
   IF (iword == setp2(i)) GO TO 2135
 END DO
 ASSIGN 2135 TO irtn
 iprm = core(iwrd)
 GO TO 9808
 
!     PLOT COMMAND CARD
 
 2200 iwrd = iwrd + 2
 mode = mode - 1
 2202 IF (core(iwrd) == 0 .OR. core(iwrd) == eor) GO TO 2215
 iword = core(iwrd)
 IF (bit64) CALL mvbits (BLANK,0,32,iword,0)
 IF (core(iwrd) == allon .OR. iword == BLANK) GO TO 2207
 DO  i = 1,28
   IF (iword == pltpr(i))  &
       GO TO (2208, 2220, 2220, 2220, 2230, 2230, 2207, 2250, 2250,  &
       2260, 2207, 2208, 2280, 2280, 2208, 2208, 2208, 2290,  &
       2207, 2281, 2207, 2248, 2240, 2207, 2245, 2207, 2207, 2208), i
   
!    1            SET  STAT  MODA  CMOD  FREQ  TRAN  CONT  RANG  TIME
!    2           PHAS  MAGN  ORIG  SYMM  ANTI   PEN  DENS  SYMB  LABE
!    3           SHAP  VECT  OUTL  THRU  MAXI  HIDD  SHRI  NOFI  FILL
!    4           OFFS
   
 END DO
 ASSIGN 2207 TO irtn
 iprm = core(iwrd)
 GO TO 9808
 
 2207 iwrd = iwrd + 2
 mode = mode - 1
 IF (mode > 0) THEN
   GO TO  2202
 ELSE
   GO TO  2210
 END IF
 
 2208 iprm = core (iwrd)
 ASSIGN 2202 TO irtn
 2209 iwrd = iwrd + 2
 mode = mode - 1
 IF (mode <= 0) GO TO 2210
 iword = core(iwrd)
 IF (bit64) CALL mvbits (BLANK,0,32,iword,0)
 IF (core(iwrd) == allon .OR. iword == BLANK) GO TO 2209
 GO TO 9814
 
 2210 IF (core(iwrd) >= 0) GO TO 2215
 IF (i /= 12) GO TO 2214
 porg = core(iwrd+1)
 IF (org <= 0) GO TO 9830
 DO  i = 1,org
   IF (porg == origin(i)) GO TO 2214
 END DO
 GO TO 9830
 2214 iwrd = iwrd + 2
 GO TO 2210
 
 2215 IF (core(iwrd) /= 0) GO TO 2216
 icont = 3
 GO TO 9998
 2216 IF (core(iwrd) /= eor) GO TO 2217
 icont = 0
 GO TO 9998
 2217 mode = core(iwrd)
 iwrd = iwrd + 1
 GO TO 2202
 
 2220 ipr1 = core(iwrd  )
 ipr2 = core(iwrd+1)
 iwrd = iwrd + 2
 mode = mode - 1
 IF (mode > 0) THEN
   GO TO  2223
 END IF
 2222 iprm = defo
 ASSIGN 2210 TO irtn
 GO TO 9806
 2223 iword = core(iwrd)
 IF (bit64) CALL mvbits (BLANK,0,32,iword,0)
 DO  i = 1,3
!                                     DEFO  VELO  ACCE
   IF (iword == idvpr(i)) THEN
      SELECT CASE ( i )
       CASE (    1)
         GO TO 2207
       CASE (    2)
         GO TO  9826
       CASE (    3)
         GO TO  9826
     END SELECT
   END IF
 END DO
 ASSIGN 2207 TO irtn
 iprm = core(iwrd)
 GO TO 9808
 
 2230 iwrd = iwrd + 2
 mode = mode - 1
 IF (mode > 0) THEN
   GO TO  2232
 END IF
 2231 ASSIGN 2210 TO irtn
 iprm = defo
 GO TO 9806
 
 2232 iword = core(iwrd)
 IF (bit64) CALL mvbits (BLANK,0,32,iword,0)
 DO  i = 1,3
   IF (iword == idvpr(i)) GO TO 2207
 END DO
 ASSIGN 2207 TO irtn
 iprm = core(iwrd)
 GO TO 9808
 
 2250 nreal = 2
 ASSIGN 2202 TO irtn
 2251 iwrd = iwrd + 2
 mode = mode - 1
 IF (mode <= 0) GO TO 2252
 iword = core(iwrd)
 IF (bit64) CALL mvbits (BLANK,0,32,iword,0)
 IF (core(iwrd) == allon .OR. iword == BLANK) GO TO 2251
 GO TO 9818
 2252 integ = 0
 DO  i = 1,nreal
   IF (core(iwrd) >=  0) GO TO 2257
   IF (core(iwrd) == -1) integ = 1
   IF (core(iwrd) == -4) iwrd  = iwrd + 1
   iwrd = iwrd + 2
 END DO
 IF (integ > 0) THEN
   GO TO  2256
 ELSE
   GO TO  2210
 END IF
 2256 ASSIGN 2210 TO irtn
 GO TO 9816
 2257 ASSIGN 2215 TO irtn
 GO TO 9818
 
 2260 iwrd = iwrd + 2
 mode = mode - 1
 nreal= 1
 IF (mode > 0) THEN
   GO TO  2262
 END IF
 2261 ASSIGN 2210 TO irtn
 iprm = lag
 GO TO 9806
 2262 iword = core(iwrd)
 IF (bit64) CALL mvbits (BLANK,0,32,iword,0)
 IF (iword == lag) GO TO 2251
 ASSIGN 2251 TO irtn
 iprm = core(iwrd)
 GO TO 9808
 
 2280 ncrd = 9
 icrd = 1
 ivc  = 0
 GO TO 2282
 2281 ncrd = 25
 icrd = 4
 ivc  = 1
 2282 ASSIGN 2210 TO irtn
 iax  = 0
 2283 iwrd = iwrd + 2
 mode = mode - 1
 IF (mode > 0) THEN
   GO TO  2284
 ELSE
   GO TO  9810
 END IF
 2284 iword = core(iwrd)
 IF (bit64) CALL mvbits (BLANK,0,32,iword,0)
 IF (core(iwrd) == allon .OR. iword == BLANK) GO TO 2283
 DO  i = icrd,ncrd
   IF (iword == coord(i)) GO TO 2286
 END DO
 IF (iax > 0) THEN
   GO TO  2202
 ELSE
   GO TO  9810
 END IF
 2286 iwrd = iwrd + 2
 mode = mode - 1
 IF (mode > 0) THEN
   GO TO  2287
 ELSE
   GO TO  2215
 END IF
 2287 IF (ivc > 0) THEN
   GO TO  2202
 END IF
 2288 IF (iax > 0) THEN
   GO TO  2202
 END IF
 2289 iax = 1
 GO TO 2284
 
 2290 iwrd = iwrd + 2
 mode = mode - 1
 IF (mode > 0) THEN
   GO TO  2291
 ELSE
   GO TO  2210
 END IF
 2291 iword = core(iwrd)
 IF (bit64) CALL mvbits (BLANK,0,32,iword,0)
 IF (core(iwrd) == allon .OR. iword == BLANK) GO TO 2290
 DO  i = 1,5
   IF (iword == lblpr(i)) THEN
      SELECT CASE ( i )
       CASE (    1)
         GO TO 2290
       CASE (    2)
         GO TO  2207
       CASE (    3)
         GO TO  2207
       CASE (    4)
         GO TO  2207
       CASE (    5)
         GO TO  2207
     END SELECT
   END IF
!                                     GRID  POIN  ELEM  BOTH  EPID
 END DO
 GO TO 2202
 
 2240 iwrd  = iwrd + 2
 mode  = mode - 1
 IF (mode > 0) THEN
   GO TO  2242
 END IF
 2241 ASSIGN 2210 TO irtn
 GO TO 9812
 2242 iword = core(iwrd)
 IF (bit64) CALL mvbits (BLANK,0,32,iword,0)
 IF (iword == defo) GO TO 2243
 ASSIGN 2243 TO irtn
 iprm  = core(iwrd)
 GO TO 9808
 2243 nreal = 1
 GO TO 2251
 
 2245 iwrd  = iwrd + 2
 iword = core(iwrd)
 IF (bit64) CALL mvbits (BLANK,0,32,iword,0)
 IF (iword == hidd) GO TO 2207
 mode  = mode - 1
 IF (mode < 0) THEN
   GO TO  2241
 ELSE
   GO TO  2210
 END IF
 
 2248 IF (core(iwrd-3) == -1 .AND. core(iwrd+2) == -1) GO TO 2207
 ASSIGN 2207 TO irtn
 GO TO 9822
 
!     SET UP ERROR MESSAGE
 
 9800 ASSIGN 9900 TO ierr
 msgno = 348
 GO TO 9890
 9802 ASSIGN 9902 TO ierr
 msgno = 349
 GO TO 9890
 9804 ASSIGN 9904 TO ierr
 msgno = 350
 GO TO 9895
 9806 ASSIGN 9906 TO ierr
 msgno = 351
 GO TO 9890
 9808 ASSIGN 9908 TO ierr
 msgno = 351
 GO TO 9890
 9810 ASSIGN 9910 TO ierr
 msgno = 352
 GO TO 9890
 9812 ASSIGN 9912 TO ierr
 msgno = 353
 GO TO 9890
 9814 ASSIGN 9914 TO ierr
 msgno = 354
 GO TO 9895
 9816 ASSIGN 9916 TO ierr
 msgno = 355
 GO TO 9890
 9818 ASSIGN 9918 TO ierr
 msgno = 356
 GO TO 9890
 9820 ASSIGN 9920 TO ierr
 msgno = 357
 GO TO 9895
 9822 ASSIGN 9922 TO ierr
 msgno = 358
 GO TO 9890
 9824 ASSIGN 9924 TO ierr
 msgno = 359
 GO TO 9890
 9826 ASSIGN 9926 TO ierr
 msgno = 360
 GO TO 9890
 9828 ASSIGN 9928 TO ierr
 msgno = 361
 GO TO 9895
 9830 ASSIGN 9930 TO ierr
 msgno = 362
 GO TO 9895
 
 9890 CALL page2 (2)
 WRITE  (nout,9891) ufm,msgno
 9891 FORMAT (a23,i4)
 IF (pltopt <= 2) nogo = 1
 GO TO 9898
 9895 CALL page2 (2)
 WRITE  (nout,9896) uwm,msgno
 9896 FORMAT (a25,i4)
 
 9898 GO TO ierr, (9900,9902,9904,9906,9908,9910,9912,9914,9916,9918,  &
     9920,9922,9924,9926,9928,9930)
 
 9900 WRITE  (nout,9901)
 9901 FORMAT (5X,'FIRST CHARACTER ON CARD IS NUMERIC. INCORRECT FORMAT',  &
     ' OR INCORRECT CONTINUATION ON PREVIOUS CARD')
 GO TO 320
 
 9902 WRITE  (nout,9903) core(iwrd)
 9903 FORMAT (5X,'PLOT COMMAND ',a4,' NOT RECOGNIZED.  CHECK SPELLING ',  &
     'AND FORMAT ON THIS CARD AND CONTINUATION ON PREVIOUS ONE')
 GO TO 320
 9904 WRITE  (nout,9905)
 9905 FORMAT (1H+,30X,' - ONLY NASTRAN GENERAL PURPOSE PLOTTER IS ',  &
     'SUPPORTED')
 GO TO 320
 
 9906 WRITE  (nout,9907) iprm
 9907 FORMAT (1H+,30X,' - KEYWORD ',a4,' NOT FOUND')
 GO TO irtn, (320,1070,2110,2136,2210,510)
 
 9908 WRITE  (nout,9909) iprm
 9909 FORMAT (1H+,30X,' - KEYWORD ',a4,' NOT RECOGNIZED')
 GO TO irtn, (320,1070,2130,2135,2202,2207,2243,2251)
 
 9910 WRITE  (nout,9911)
 9911 FORMAT (1H+,30X,' - COORDINATE AXES INCORRECTLY DEFINED')
 GO TO irtn, (320,2210)
 
 9912 WRITE  (nout,9913)
 9913 FORMAT (1H+,30X,' - INCORRECT FORMAT')
 GO TO irtn, (320,1070,2210)
 
 9914 WRITE  (nout,9915) iprm
 9915 FORMAT (1H+,30X,3H - ,a4,' IDENTIFICATION NUMBER NOT DEFINED')
 GO TO irtn, (320,1005,1910,2120,2202)
 
 9916 WRITE  (nout,9917)
 9917 FORMAT (1H+,30X,' - DATA TYPE IS INCORRECT')
 GO TO irtn, (1005,1730,2107,2110,2210,320)
 
 9918 WRITE  (nout,9919)
 9919 FORMAT (1H+,30X,' - ONE OR MORE REQUIRED REAL VALUES MISSING')
 GO TO irtn, (320,1005,2202,2215)
 
 9920 WRITE  (nout,9921)
 9921 FORMAT (1H+,30X,' - CAMERA OPTION NOT SPECIFIED')
 GO TO 320
 
 9922 WRITE  (nout,9923)
 9923 FORMAT (1H+,30X,' - THRU MUST BE PRECEDED AND FOLLOWED BY INTEGER'  &
     ,      ' VALUES')
 GO TO irtn, (2123,2207)
 
 9924 WRITE  (nout,9925)
 9925 FORMAT (1H+,30X,' - THRU RANGE OVERLAPS RANGE OF PREVIOUS THRU')
 GO TO 2123
 
 9926 WRITE  (nout,9927) ipr1,ipr2
 9927 FORMAT (1H+,30X,' - ONLY DEFORMATION VALID WITH ',2A4)
 GO TO 2207
 
 9928 WRITE  (nout,9929) forg,porg
 9929 FORMAT (1H+,30X,' - A NEW ORIGIN',i8,' WAS DEFINED IN A FIND ',  &
     'CARD, BUT IT IS NOT USED BY THE IMMEDIATE PLOT CARD',  &
     /5X,'(ORIGIN',i8,' WILL BE USED FOR THIS PLOT)',/)
 GO TO 9999
 
 9930 WRITE  (nout,9931) porg
 9931 FORMAT (1H+,30X,' - ORIGIN',i8,' IS UNDEFINED')
 GO TO 2207
 
 9998 IF (isplot == 0 .OR. porg == -1) RETURN
 IF (porg == 0) porg = porg1
 IF (forg /= 0 .AND. forg /= porg) GO TO 9828
 9999 forg = 0
 porg = 0
 RETURN
END SUBROUTINE ifp1pc
