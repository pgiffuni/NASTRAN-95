SUBROUTINE trbscd
     
!     THIS SUBROUTINE CALCULATES THE STIFFNESS AND MASS MATRICES FOR
!     THE BASIC BENDING TRIANGLE.  THE MASS MATRIX MAY BE CALCULATED
!     EITHER BY THE CONVENTIONAL OR THE CONSISTENT MASS METHODS (USING
!     EMASTQ OR INCLUDED CODE) ACCORDING TO THE PARAMETER ICMBAR.
!     THIS ELEMENT MAY NOT BE USED IN A HEAT PROBLEM.
 
!     DOUBLE PRECISION VERSION
 
!     ECPT FOR THIS ELEMENT
 
!     INDEX  NAME      TYPE      DESCRIPTION
!     ----- -------    ----    ------------------
!      1    IELID        I     ELEMENT ID
!      2    NGRID(1)     I     FIRST GRID POINT
!      3    NGRID(2)     I     SECOND GRID POINT
!      4    NGRID(3)     I     THIRD GRID POINT
!      5    ANGLE        R     ANGLE OF MATERIAL
!      6    MATID1       I     MATERIAL ID 1
!      7    EYE          R     MOMENT OF INERTIA
!      8    MATID2       I     MATERIAL ID 2
!      9    T2           R     T2
!     10    FMU          R     NON-STRUCTURAL MASS
!     11    Z11          R     Z1
!     12    Z22          R     Z2
!     13    NECPT(13)    I     COORD  SYSTEM ID 1
!     14    X1           R
!     15    Y1           R     COORDINATES
!     16    Z1           R
!     17    NECPT(17)    I     COORD SYSTEM ID 2
!     18    X2           R
!     19    Y2           R     COORDINATES
!     20    Z2           R
!     21    NECPT(21)    I     COORD SYSTEM ID 3
!     22    X3           R
!     23    Y3           R     COORDINATES
!     24    Z3           R
!     25    ELTEMP       R     ELEMENT TEMPERATURE
 
 LOGICAL :: iheat,nogo
 INTEGER :: elid,estid,dict(9),ipart(3),necpt(25)
 REAL :: ecpt(25)
 DOUBLE PRECISION :: a,prod,temp9,xsubb,sxubc,ysubc,bfact,e,kout  &
     ,                kk,ksav,m(324),mout(324)
 COMMON /emgprm/  ixtr,jcore,ncore,dm(12),ismb(3),iprec,nogo,heat, icmbar
 COMMON /emgdic/  qq,ldict,ngrids,elid,estid
 COMMON /emgest/  ielid,ngrid(3)
 COMMON /emgtrx/  a(225),prod(9),temp9(9),xsubb,sxubc,ysubc,bfact,  &
     e(18),kout(324),kk(324),ksav(81)
 COMMON /system/  ksystm(60)
 EQUIVALENCE      (ksystm(2),ioutpt),(ksystm(56),iheat),  &
     (ecpt(1),necpt(1),ielid),(dict5,dict(5)), (kk(1),mout(1)),(kout(1),m(1))
 DATA    ipart /  1,2,3 /
 
 ip = iprec
 
!     IF THIS IS A HEAT PROBLEM THIS SHOULD NOT CALL US, SO RETURN
 
 IF (iheat) RETURN
 
!     CREATE AN ARRAY POINTING TO THE GRID POINTS IN INCREASING  SIL
!     ORDER
 
 100 DO  i = 1,2
   ip1 = i + 1
   ii = ipart(i)
   DO  j = ip1,3
     jj = ipart(j)
     IF (ngrid(ii) <= ngrid(jj)) CYCLE
     ipart(i) = jj
     ipart(j) = ii
     ii = jj
     GO TO 100
   END DO
 END DO
 
!     IF STIFFNESS MATRIX IS DESIRED CALL ETRBKD, OTHERWISE ONLY MASS
!     MATRIX IS DESIRED
 
 IF (ismb(1) == 0) GO TO 200
 
 CALL etrbkd (0)
 IF (nogo) RETURN
 dict5 = bfact
 
!     RE ORDER THE MATRIX BY INCREASING SIL VALUE.    NOTE THAT
 
!     KK  = KK(1 TO  9)     KK   = KK(10 TO 18)     KK   = KK(19 TO  27)
!       AA                    AB                      AC
 
!     KK  = KK(28 TO  36)  KK   = KK(37 TO  45)   KK   =  KK(46 TO  54)
!       BA                   BB                     BC
 
!     KK  = KK(55 TO  63)  KK   = KK(64 TO  72)   KK  =  KK(73 TO  81)
!       CA                   CB                     CC
 
!     AND
 
!     KOUT  = KOUT(1 - 36) KOUT  = KOUT( 4 - 6)   KOUT  = KOUT( 7 -  9)
!         I I    (10 - 12)     I I     (13 - 15)      I I     (16 - 18)
!          1 1   (19 - 21)      1 2    (22- 24)        1 3    (25 - 27)
 
!     ETC
 
 
 DO  i = 1,3
   ii = ipart(i)
   DO  j = 1,3
     jj = ipart(j)
     DO  k = 1,3
       DO  l = 1,3
         ik   = (ii-1)*27 + (jj-1)*9 + (k-1)*3 + l
         iout = (i -1)*27 + (j -1)*3 + (k-1)*9 + l
         kout(iout) = kk(ik)
       END DO
     END DO
   END DO
 END DO
 
!     NOW OUTPUT THE MATRIX
 
 dict(1) = estid
 dict(2) = 1
 dict(3) = 9
 dict(4) = 4 + 8 + 16
 
 CALL emgout (kout,kout,81,1,dict,1,ip)
 
!     NOW CALCULATE THE MASS MATRIX IF NEEDED
 
 200 IF (ismb(2) == 0) RETURN
 
!     WHICH MASS METHOD TO BE USED (CONVENTIONAL OR CONSISTENT)
 
 IF (icmbar >= 0) GO TO 300
 
 CALL emadtq (3,m)
 
!     REORDER THE DIAGONAL MASS MATRIX
 
 DO  i = 1,3
   ii = (i-1)*3 + 1
   ij = ipart(i)
   jj = (ij-1)*3 + 1
   DO  j = 1,3
     iout = ii + j - 1
     ik   = jj + j - 1
     mout(iout) = m(ik)
   END DO
 END DO
 
!     NOW OUTPUT THE MATRIX
 
 dict(1) = estid
 dict(2) = 2
 dict(3) = 9
 dict(4) = 7
 
 CALL emgout (mout,mout,9,1,dict,2,ip)
 
 RETURN
 
!     THE COUPLED MASS MATRIX CALCULATIONS ARE MADE HERE VIA ETRBMD
 
 300 CALL etrbmd
 IF (nogo) RETURN
 
!     INSERT THE MATRICES INTO THE OUTPUT MATRIX IN INCREASING SIL ORDER
 
 DO  i = 1,3
   ii = ipart(i)
   DO  j = 1,3
     jj = ipart(j)
     DO  k = 1,3
       DO  l = 1,3
         ia   = (ii-1)*36 + (jj-ii)*9 + (k-1)*3 + l
         iout = (i -1)*27 + (j - 1)*3 + (k-1)*9 + l
         mout(iout) = m(ia)
       END DO
     END DO
   END DO
 END DO
 
!     NOW OUTPUT THE MASS MATRIX
 
 dict(1) = estid
 dict(2) = 1
 dict(3) = 9
 dict(4) = 4 + 8 + 16
 
 CALL emgout (mout,mout,81,1,dict,2,ip)
 RETURN
 
END SUBROUTINE trbscd
