SUBROUTINE ifs4p (*,*,*)
     
    LOGICAL :: abort,baddat,badfor,lharm,lflsym,fphys1,ifpdco
    INTEGER :: t1,t4,thru,SAVE(24),nm(2),ty1,ty2,  &
        ret,bcdyes,bcdno,bcds,bcda,bcdnon,bcdaxi
    REAL :: z(100)
    CHARACTER (LEN=25) :: sfm
    CHARACTER (LEN=29) :: uim
    CHARACTER (LEN=25) :: uwm
    CHARACTER (LEN=23) :: ufm
    COMMON /machin/ mach
    COMMON /xmssg / ufm,uwm,uim,sfm
    COMMON /system/ nbuf,nout,abort
    COMMON /ifpx1 / ncds,t1(2,310)
    COMMON /ifpx3 / t4(2,314)
    COMMON /ifpdta/ id,n,k,kx,ky,i(100),m(100),mf(100),m1(100),  &
        m1f(100),kn,baddat,badfor,nopen,nparam,iax,nax,  &
        iaxf,naxf,lharm,knt,slotdf(5),gc(7),ll(6)
    COMMON /cifs4p/ j(20),km,lflsym,fphys1
    EQUIVALENCE     (z(1),m(1)),(kout,j(2))
    DATA    thru  , bcdyes,bcdno /4HTHRU,4HYES ,4HNO  /
    DATA    bcds  , bcda,bcdnon  /4HS   ,4HA   ,4HNONE/
    DATA    bcdaxi/ 4HAXIS/
 
    IF (k > 100) GO TO 81
    GO TO (     5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
        5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
        5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
        5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
        5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
        5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
        5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
        5,   5,   5,   5,   5,   5,   5,   5, 790, 800,  &
        5,   5,   5,   5,   5,   5,   5,   5,   5, 900,  &
        900,   5,   5,   5,   5,   5,   5, 980,   5,   5 ), k
81  IF (kx > 100) GO TO 82
    GO TO (     5,1020,   5,1040,1050,   5,   5,   5,   5,   5,  &
        5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
        5,1220,   5,1050,   5,   5,   5,   5,   5,   5,  &
        5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
        5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
        5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
        5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
        5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
        5,   5,   5,   5,   5,   5,   5,   5,1020,   5,  &
        5,   5,   5,   5,1950,1960,   5,   5,1990,   5 ), kx
82  IF (ky > 100) GO TO 83
    GO TO (  2100,2200,2300,2400,3100,3200,3300,3400,3500,3600,  &
        3700,3800,3900,4000,   5,   5,4100,4200,4300,4400,  &
        4500,4600,   5,   5,   5,   5,   5,   5,   5,   5,  &
        5,   5,   5,   5,   5,   5,   5,   5,1990,   5,  &
        5,1990,   5,   5,   5,   5,   5,   5,   5,   5,  &
        5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
        5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
        5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
        5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
        5,   5,   5,   5,6501,6601,   5,   5,   5,   5 ), ky
83  kz = k - 300
    IF (kz > 39) GO TO 5
    GO TO (     5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
        5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
        8000,9000,9100,9200,9300,9000,9400,9500,   5,   5,  &
        5,   5,4300,4400,4100,4200,   5,   5,   5      ), kz
5   CALL page2 (2)
    WRITE  (nout,6) sfm
6   FORMAT (a25,' 322, ILLEGAL ENTRY TO IFS4P.')
    abort =.true.
    RETURN 1
7   badfor =.true.
    RETURN 1
8   baddat =.true.
    RETURN 1
    3 DO  l = 1,n
        i(l) = m(l)
    END DO
2   RETURN
9   RETURN 3
 
    !******              79-CTRIARG,80-CTRAPRG             ****************
 
790 i1 = 4
    GO TO 791
800 i1 = 5
791 IF (m(1) <= 0 .OR. m(i1+2) <= 0) GO TO 8
    DO  l = 2,i1
        IF (m(l) <= 0) GO TO 8
        IF (l == 2) CYCLE
        DO  l1 = l,i1
            IF (m(l-1) == m(l1)) GO TO 8
        END DO
    END DO
    n = i1 + 2
    GO TO 3
 
    !*******       MATS1,MATT1        **************************************
 
    900 DO  l = 1,11
        IF (m(l) < 0) GO TO 8
        i(l) = m(l)
    END DO
    n = 11
    GO TO 2
 
    !*******       TEMPD              **************************************
 
    980 DO  l = 1,7,2
        IF (m(l) == 0 .AND. m(l+1) == 0) CYCLE
        IF (m(l) <= 0) GO TO 8
        n = n + 2
        i(n-1) = m(l  )
        i(n  ) = m(l+1)
        IF (n <= 2) CYCLE
        DO  l1 = 4,n,2
            IF (i(n-1) == i(l1-3)) GO TO 8
        END DO
    END DO
    IF (n > 0) THEN
        GO TO     2
    ELSE
        GO TO     8
    END IF
 
    !**************    MATT2,189-MATT3     *********************************
 
    1020 DO  l = 1,16
        IF (m(l) < 0) GO TO 8
        i(l) = m(l)
    END DO
    IF (m(1) == 0) GO TO 8
    n = 16
    GO TO 2
 
    !******           104-CTORDRG           ************************
 
1040 IF (m(1) <= 0 .OR. m(3) <= 0 .OR. m(4) <= 0 .OR. m(3) == m(4) .OR.  &
        z(5) < 0.0 .OR. z(5) > 180.0 .OR. z(6) < 0.0 .OR. z(6) > 180.0) GO TO 8
    IF (mf(2) == 0) m(2) = m(1)
    IF (m(2)  <= 0) GO TO 8
    n = 7
    GO TO 3
 
    !*******       SPOINT,124-EPOINT    ************************************
 
1050 IF (mf(2) == 3) GO TO 1056
    DO  l = 1,8
        IF (mf(l) /= 1 .AND. mf(l) /= 0) GO TO 7
        IF (m(l) < 0) THEN
            GO TO     8
        ELSE IF (m(l) == 0) THEN
            GO TO  1055
        END IF
1052    IF (m(l) > 999999) GO TO 8
        n = n + 1
        i(n) = m(l)
        IF (n <= 1) CYCLE
        DO  l1 = 2,n
            IF (i(n) == i(l1-1)) GO TO 8
        END DO
1055 CONTINUE
     END DO
     IF (n > 0) THEN
         GO TO     2
     ELSE
         GO TO     8
     END IF
1056 IF (m(2) /= thru) GO TO 8
     IF (mf(1) /= 1 .OR. mf(3) /= 1) GO TO 7
     k2078 = 208
     IF (k == 124) k2078 = 207
     l1 = 1
     l2 = 4
     DO  l = l2,8
         IF (mf(l) /= 0) GO TO 7
     END DO
     IF (m(l2) > 9999999) GO TO 8
     ii = m(l1) - 1
     l2 = m(l2) - m(l1)
     IF (ii < 0 .OR. l2 <= 0) GO TO 8
     l1 = 1
     DO  l = 1,l2
         ii = ii + 1
         CALL WRITE (k2078,ii,1,0)
     END DO
     i(1) = ii + 1
     n = 1
     GO TO 2
 
     !*******         122-MAT3        *****************************
 
1220 IF (m(1) <= 0 .OR. z(2) < 0. .OR. z(3) < 0. .OR. z(4) < 0. .OR.  &
         z(9) < 0. .OR. z(10) < 0. .OR. z(11) < 0.) GO TO 8
     IF (ABS(z(5)) <= 1. .AND. ABS(z(6)) <= 1. .AND. ABS(z(7)) <= 1.) GO TO 1222
     CALL page2 (2)
     WRITE  (nout,1221) uwm,t1(1,k),t1(2,k),knt
1221 FORMAT (a25,' 301, BULK DATA CARD ',2A4,' CONTAINS INCONSISTENT',  &
         ' DATA.',10X,'SORTED CARD COUNT =',i7)
1222 n = 16
     GO TO 3
 
 
     !*******       195-RANDPS       ****************************************
 
1950 IF (m(1) <= 0 .OR. m(2) <= 0 .OR. m(3) < m(2) .OR. m(6) < 0) GO TO 8
     IF (m(2) == m(3) .AND. z(5) /= 0.0) GO TO 8
     n = 6
     IF (kout <= 2) GO TO 1955
     IF (m(1) == j(kout)) GO TO 3
     IF (kout == j(   1)) GO TO 8
1955 kout = kout + 1
     j(kout) = m(1)
     GO TO 3
 
     !*******       196-RANDT1       ****************************************
 
1960 IF (kout <= 2) GO TO 8
     DO  in = 3,kout
         IF (m(1) == j(in)) GO TO 1962
     END DO
     GO TO 8
1962 IF (m(1) <= 0 .OR. m(2) <= 0 .OR. z(3) < 0.0 .OR. z(4) <= z(3)) GO TO 8
     n = 4
     GO TO 3
 
     !*****         199-PLOAD2,239-QBDY1,242-QVOL   *************************
 
1990 IF (km   /= 0) GO TO 1991
     IF (mf(1) /= 1 .OR. mf(2) /= 2 .AND. mf(2) /= 0) GO TO 7
     IF (m(1) <= 0) GO TO 8
     l = 3
     isid = m(1)
     iqvl = m(2)
     GO TO 1992
1991 l = 1
1992 IF (mf(8) == 3) GO TO 7
     ntot  = 0
     k2078 = 209
1993 IF (m(l)  == 0) GO TO 1998
     IF (m(l)  < 0) GO TO 8
     IF (mf(l) == 3) GO TO 7
     IF (mf(l+1) == 3) GO TO 1995
     IF (mf(l) /= 1 .AND. mf(l) /= 0) GO TO 7
     n = n + 3
     i(n-2) = isid
     i(n-1) = iqvl
     i(n) = m(l)
     l = l + 1
     IF (n < 48) GO TO 1997
     CALL WRITE (k2078,i,n,0)
     ntot = ntot + n
     n = 0
     GO TO 1997
1995 IF (m(l+1) /= thru) GO TO 8
     IF (mf(l+3) /= 1 .AND. mf(l+3) /= 0) GO TO 7
     l1 = m(l  ) - 1
     l2 = m(l+3) - l1
     IF (l2 <= 1 .OR. l1 < 0) GO TO 8
     DO  ii = 1,l2
         n = n + 3
         i(n-2) = isid
         i(n-1) = iqvl
         i(n) = ii + l1
         IF (n < 48) CYCLE
         CALL WRITE (k2078,i,n,0)
         ntot = ntot + n
         n = 0
     END DO
     l = l + 4
1997 IF (l <= 8) GO TO 1993
1998 t4(2,k) = t4(2,k) + ntot
     IF (m1(1) == 0 .AND. m1(2) == 0) GO TO 1999
     km = 0
     GO TO 2
1999 km = 1
     GO TO 2
 
     !**********          201-TEMPP1          *******************************
 
2100 IF (km /= 0) GO TO 2120
     nn = 6
     n  = 6
     id = m(1)
     IF (mf(5) == -32767) GO TO 2106
     IF (mf(7) /= 0 .OR. mf(8) /= 0) badfor =.true.
     2101 DO  l = 3,6
         IF (mf(l) == 0 .OR. mf(l) == 2) CYCLE
         badfor =.true.
     END DO
2103 CONTINUE
     IF (mf(1) /= 1 .OR. mf(2) /= 1) badfor =.true.
     IF (m(1) <= 0 .OR. m(2) <= 0) baddat =.true.
     DO  l = 1,n
         i(l) = m(l)
         SAVE(l) = m(l)
     END DO
     GO TO 2110
     2106 DO  l = 3,4
         IF (mf(l) == 0 .OR. mf(l) == 2) CYCLE
         badfor =.true.
     END DO
     GO TO 2103
2110 IF (m1(1) == 0 .AND. m1(2) == 0) GO TO 2115
     km = 0
     kn = 0
     GO TO 9
2115 kn = 1
     km = km + 1
     GO TO 9
2120 IF (mf(2) == 3 .OR. mf(5) == 3) GO TO 2150
     n  = 0
     DO  l = 1,8
         IF (mf(l) == 0) CYCLE
         IF (mf(l) == 1) GO TO 2125
         IF (mf(l) == -32767) EXIT
         badfor =.true.
         CYCLE
2125     IF (m(l) > 0) GO TO 2130
         baddat =.true.
         CYCLE
2130     SAVE(2) = m(l)
         CALL WRITE (209,SAVE,nn,0)
     END DO
2145 CONTINUE
     GO TO 2110
2150 n = 0
     IF (mf(7) == 0 .AND. mf(8) == 0) GO TO 2155
     IF (mf(4) == 0 .AND. mf(5) == -32767) GO TO 2155
     badfor =.true.
     GO TO 2110
2155 l1 =-1
     DO  l = 1,4,3
         IF (mf(l) == 0 .AND. mf(l+1) == 0 .AND. mf(l+2) == 0) CYCLE
         IF (mf(l) == 1 .AND. mf(l+1) == 3 .AND. mf(l+2) == 1) GO TO 2160
         IF (mf(l+1) == -32767) EXIT
         badfor =.true.
         CYCLE
2160     l1 = l1 + 1
         l2 = l1 + l
         IF (m(l2) > 0 .AND. m(l2+1) == thru .AND. m(l2+3) > m(l2)) GO TO 2165
         baddat =.true.
         CYCLE
2165     l3 = m(l2  )
         l4 = m(l2+3)
         DO  l5 = l3,l4
             SAVE(2) = l5
             CALL WRITE (209,SAVE,nn,0)
         END DO
     END DO
2185 CONTINUE
     GO TO 2110
 
     !*******       202-TEMPP2         **************************************
 
2200 IF (km /= 0) GO TO 2120
     nn = 8
     n  = 8
     id = m(1)
     IF (mf(5) == -32767) GO TO 2106
     IF (mf(7) /= 0 .AND. mf(7) /= 2 .OR.  &
         mf(8) /= 0 .AND. mf(8) /= 2) badfor =.true.
     GO TO 2101
 
     !*******       203-TEMPP3         **************************************
 
2300 IF (km /= 0) GO TO 2330
     nn = 24
     n  = 0
     id = m(1)
     l1 = 1
     IF (mf(1) /= 1 .OR. mf(2) /= 1) badfor =.true.
     DO  l = 3,8
         IF (mf(l) == 0 .OR. mf(l) == 2) CYCLE
         IF (mf(l) == -32767) GO TO 2302
         badfor =.true.
         CYCLE
         2302 DO  l5 = l,8
             m(l5) = 0
         END DO
         mf(7) = 0
         mf(8) = 0
     END DO
     IF (m(1) <= 0 .OR. m(2) <= 0) baddat =.true.
     IF (z(3) >= z(5)) baddat =.true.
     zz = z(5)
     IF (mf(7) == 0 .AND. mf(8) == 0) GO TO 2310
     IF (zz >= z(7)) baddat =.true.
2310 zz = z(7)
     DO  l = 1,8
         i(l) = m(l)
         SAVE(l) = m(l)
     END DO
2321 l1 = l1 + 8
     IF (m1(1) == 0 .AND. m1(2) == 0) GO TO 2328
     km = 0
     kn = 0
2322 IF (l1 > nn) GO TO 2326
     DO  l = l1,nn
         i(l) = 0
         SAVE(l) = 0
     END DO
2326 n  = nn
     GO TO 9
2328 km = km + 1
     kn = 1
     IF (km-3 < 0) THEN
         GO TO     9
     ELSE
         GO TO  2322
     END IF
2330 IF (km > 2) GO TO 2120
     n  = 0
     l3 = 8*km
     DO  l = 1,7,2
         IF (mf(l) == 0 .AND. mf(l+1) == 0) GO TO 2340
         IF (mf(l) /= -32767) GO TO 2335
         mf(7) = 0
         mf(8) = 0
         GO TO 2340
2335 CONTINUE
     IF (mf(l  ) /= 0 .AND. mf(l  ) /= 2  .OR.  &
         mf(l+1) /= 0 .AND. mf(l+1) /= 2) badfor =.true.
     IF (zz >= z(l)) baddat =.true.
2340 zz = z(l)
     l5 = l3 + l
     i(l5) = m(l)
     SAVE(l5) = m(l)
     i(l5+1)  = m(l+1)
     SAVE(l5+1) = m(l+1)
 END DO
 GO TO 2321
 
 !*******       204-TEMPRB         **************************************
 
2400 IF (km /= 0) GO TO 2430
 nn = 16
 n  = 0
 id = m(1)
 l1 = 1
 IF (mf(1) /= 1 .OR. mf(2) /= 1) badfor =.true.
 DO  l = 3,8
     IF (mf(l) == 0 .OR. mf(l) == 2) CYCLE
     IF (mf(l) == -32767) GO TO 2402
     badfor =.true.
     CYCLE
     2402 DO  l5 = l,8
         m(l5) = 0
     END DO
 END DO
 IF (m(1) <= 0 .OR. m(2) <= 0) baddat =.true.
 DO  l = 1,8
     i(l) = m(l)
     SAVE(l) = m(l)
 END DO
2421 l1 = l1 + 8
 IF (m1(1) == 0 .AND. m1(2) == 0) GO TO 2428
 km = 0
 kn = 0
2422 IF (l1 > nn) GO TO 2426
 DO  l = l1,nn
     i(l) = 0
     SAVE(l) = 0
 END DO
2426 n  = nn
 GO TO 9
2428 km = km + 1
 kn = 1
 IF (km-2 < 0) THEN
     GO TO     9
 ELSE
     GO TO  2422
 END IF
2430 IF (km > 1) GO TO 2120
 n  = 0
 DO  l = 1,8
     IF (mf(l) == -32767) GO TO 2455
     IF (mf(l) /= 0 .AND. mf(l) /= 2) badfor =.true.
     i(l+8) = m(l)
     SAVE(l+8) = m(l)
 END DO
 GO TO 2421
 2455 DO  l = 5,8
     i(l+8) = 0
     SAVE(l+8) = 0
 END DO
 GO TO 2421
 
 !    TEMPG IS MODELLED AFTER TEMPP3
 !    TEMPP4 IS MODELLED AFTER TEMPP1,EXCEPT THAT TEMPP1 HAS ONE LESS C
 
 
 !*******      295-TEMPG     ********************************************
 
6501 CONTINUE
     GO TO 2300
 
 !*******      296-TEMPP4    ********************************************
 
6601 CONTINUE
     IF (km /= 0) GO TO 6630
     nn = 14
     n  = 0
     id = m(1)
     l1 = 1
     IF (mf(1) /= 1 .OR. mf(2) /= 1) badfor =.true.
     DO  l = 3,8
         IF (mf(l) == 0 .OR. mf(l) == 2) CYCLE
         IF (mf(l) == -32767) GO TO 6602
         badfor =.true.
         CYCLE
         6602 DO  l5 = l,8
             m(l5) = 0
         END DO
     END DO
     IF (m(1) <= 0 .OR. m(2) <= 0) baddat =.true.
     DO  l = 1,8
         i(l) = m(l)
         SAVE(l) = m(l)
     END DO
6621 l1 = l1 + 8
     IF (m1(1) == 0 .AND. m1(2) == 0) GO TO 6628
     km = 0
     kn = 0
6622 IF (l1 > nn) GO TO 6626
     DO  l = l1,nn
         i(l) = 0
         SAVE(l) = 0
     END DO
6626 n  = nn
     GO TO 9
6628 km = km + 1
     kn = 1
     IF (km-2 < 0) THEN
         GO TO     9
     ELSE
         GO TO  6622
     END IF
6630 IF (km > 1) GO TO 2120
     n  = 0
     l3 = 8*km
     IF (mf(7) /= 0 .AND. mf(8) /= 0) badfor =.true.
     DO  l = 1,6
         IF (mf(l) == 0 .OR. mf(l) == 2) CYCLE
         IF (mf(l) == -32767) GO TO 6632
         badfor =.true.
         CYCLE
         6632 DO  l6 = l,6
             m(l6) = 0
         END DO
     END DO
     DO  l = 1,6
         l5 = l3 + l
         i(l5) = m(l)
         SAVE(l5) = m(l)
     END DO
     GO TO 6621
 
     !*******       205-GRIDB          **************************************
 
3100 ASSIGN 3105 TO ret
     GO TO 3890
3105 IF (m(1) <= 0 .OR. m(6) < 0 .OR. m(8) <= 0) GO TO 8
     IF (ifpdco(m(7))) GO TO 8
     n    = 5
     i(1) = m(1)
     i(2) = m(4)
     i(3) = m(6)
     i(4) = m(7)
     i(5) = m(8)
     GO TO 2
 
     !*******       206-FSLIST         **************************************
 
3200 IF (km /= 0) GO TO 3270
     ASSIGN 3205 TO ret
     GO TO 3890
3205 IF (mf(1) == 0 .OR. mf(1) == 2) GO TO 3207
3206 badfor =.true.
     GO TO 3250
3207 IF (mf(1) == 0 .OR. (mf(1) == 2 .AND. z(1) > 0.0)) GO TO 3209
3208 baddat =.true.
     GO TO 3250
3209 IF (mf(1) == 0) m(1) = 1
     i(1) = m(1)
     n  = 1
     l1 = 2
     l2 = 0
     IF (mf(2) /= 3) GO TO 3220
     IF (m(2)  /= bcdaxi) GO TO 3208
     n  = n + 1
     i(n) = 0
     l1 = l1 + 1
     l2 = 1
     3220 DO  l = l1,8
         l3 = l + l2
         IF (mf(l) == 3) GO TO 3230
         IF (mf(l) == 0) GO TO 3235
         IF (mf(l) /= 1) GO TO 3206
         IF (m(l3) <= 0) GO TO 3208
         n  = n + 1
         i(n) = m(l3)
     END DO
3227 IF (n > 0) THEN
         GO TO  3250
     ELSE
         GO TO  3208
     END IF
3230 IF (m(l3) /= bcdaxi) GO TO 3208
     n  = n + 1
     i(n) = 0
3235 IF (l == 8) GO TO 3245
     l  = l + 1
     DO  l2 = l,8
         IF (mf(l2) /= 0) GO TO 3206
     END DO
3245 IF (m1(1) == 0 .AND. m1(2) == 0) GO TO 3208
     GO TO 3227
3250 IF (m1(1) == 0 .AND. m1(2) == 0) GO TO 3255
     km = 0
     kn = 0
     n  = n + 1
     i(n) =-1
     GO TO 9
3255 kn = 1
     km = km + 1
     GO TO 9
3270 l1 = 1
     l2 = 0
     GO TO 3220
 
     !*******       207-RINGFL         **************************************
 
3300 ASSIGN 3310 TO ret
     GO TO 3890
     3310 DO  l = 1,5,4
         IF (m(l) == 0.AND.m(l+1) == 0.AND.m(l+2) == 0.AND.m(l+3) == 0) CYCLE
         IF (m(l) <= 0 .OR. z(l+1) <= 0.0) GO TO 8
         n = n + 4
         IF (n > 4 .AND. m(l) == m(l-4)) GO TO 8
         IF (m(l) <= 99999) GO TO 3320
         CALL page2 (2)
         WRITE (nout,3512) ufm
         GO TO 8
3320     i(n-3) = m(l  )
         i(n-2) = m(l+1)
         i(n-1) = m(l+2)
         i(n  ) = m(l+3)
     END DO
     IF (n > 0) THEN
         GO TO     2
     ELSE
         GO TO     8
     END IF
 
     !*******       208-PRESPT         **************************************
 
3400 ASSIGN 3410 TO ret
     GO TO 3890
3410 IF (m(1) <= 0) GO TO 8
     DO  l = 3,7,2
         IF (m(l) == 0 .AND. m(l+1) == 0) CYCLE
         IF (m(l) <= 0) GO TO 8
         n      = n + 3
         i(n-2) = m(1)
         i(n-1) = m(l)
         i(n  ) = m(l+1)
     END DO
     IF (n > 0) THEN
         GO TO     2
     ELSE
         GO TO     8
     END IF
 
     !*******       209-CFLUID2        **************************************
 
3500 kfl = 2
3505 ASSIGN 3510 TO ret
     GO TO 3890
3510 IF (m(1) <= 0) GO TO 8
     IF (m(1) <= 99999) GO TO 3513
     CALL page2 (2)
     WRITE  (nout,3512) ufm
3512 FORMAT (a23,' 5004, FLUID POINT ID ON CFLUID OR RINGFL CARD ',  &
         'EXCEEDS 999999 LIMIT')
     GO TO 8
     3513 DO  l = 2,kfl
         IF (m(l) <= 0) GO TO 8
         IF (l  == kfl) CYCLE
         l2 = l + 1
         DO  l1 = l2,kfl
             IF (m(l) == m(l1)) GO TO 8
         END DO
     END DO
     i(1) = m(1)
     n    = kfl + 3
     IF (mf(6) == 0) m(6) = 1
     IF (mf(7) == 0) m(7) = 1
     i(kfl+2) = m(6)
     i(kfl+3) = m(7)
     DO  l = 1,kfl
         i(l+1) = m(l+1)
     END DO
     GO TO 2
 
     !*******       210-CFLUID3        **************************************
 
3600 kfl = 3
     GO TO 3505
 
     !*******       211-CFLUID4        **************************************
 
3700 kfl = 4
     GO TO 3505
 
     !*******       212-AXIF           **************************************
 
3800 n = 0
     IF (km   /= 0) GO TO 3850
     IF (iaxf > 0) GO TO 3840
     iaxf = iaxf + 1
     IF (mf(1) /= 1 .OR. mf(2) /= 0 .AND. mf(2) /= 2 .OR. mf(3) /= 0  &
         .AND. mf(3) /= 2 .OR. mf(4) /= 0 .AND. mf(4) /= 2 .OR.  &
         mf(5) /= 3) badfor =.true.
     IF (mf(7) /= 0 .OR. mf(8) /= 0 .OR. mf(6) /= 0 .AND. mf(6) /= 3)  &
         badfor =.true.
     IF (mf(3) == 0) m(3) = 1
     IF (m(5) /= bcdyes .AND. m(5) /= bcdno) baddat =.true.
     IF (m(5) == bcdyes) m(5) = 1
     IF (m(5) == bcdno ) m(5) = 0
     CALL WRITE (215,m,5,0)
     IF (mf(6) == 3) GO TO 3820
     IF (m1(1) /= 0 .OR. m1(2) /= 0) baddat =.true.
     GO TO 3875
3820 IF (m(7) /= bcdnon) baddat =.true.
     IF (m1(1) == 0 .AND. m1(2) == 0) baddat =.true.
     GO TO 3875
3840 CALL page2 (2)
     WRITE  (nout,3841) ufm
3841 FORMAT (a23,' 4121, ONLY ONE (1) AXIF CARD ALLOWED IN BULK DATA.')
     abort =.true.
     GO TO 3875
3850 IF (mf(2) == 3) GO TO 3860
     DO  l = 1,8
         IF (mf(l) == 0) CYCLE
         IF (mf(l) == 1) GO TO 3853
         badfor =.true.
         n = 0
         GO TO 3856
3853     IF (m(l) <= naxf) baddat =.true.
         n = n + 1
         naxf = m(l)
         i(n) = m(l)
     END DO
     IF (n <= 0) baddat =.true.
3856 GO TO 3875
3860 IF (mf(4) == 3) GO TO 3870
     l1 = 1
     l2 = 1
     IF (mf(1) == 1 .AND. mf(3) == 1) GO TO 3862
3861 badfor =.true.
     GO TO 3875
     3862 DO  l = 4,8
         IF (mf(l) /= 0) GO TO 3861
     END DO
     IF (m(1) < m(4) .AND. m(1) >= 0) GO TO 3864
     baddat =.true.
     GO TO 3875
3864 IF (m(1) <= naxf) baddat =.true.
     IF (m(1) > 0) GO TO 3866
     CALL WRITE (215,0,1,0)
     GO TO 3867
3866 l2 = m(1)
3867 l3 = m(4)
     DO  l = l2,l3,l1
         CALL WRITE (215,l,1,0)
     END DO
     naxf = l3
     GO TO 3875
3870 l1 = m(7)
     l2 = l1
     IF (mf(1) == 1 .AND. mf(3) == 1 .AND. mf(5) == 1 .AND. mf(6) == 0  &
         .AND. mf(7) == 0 .AND. mf(8) == 0) GO TO 3872
     GO TO 3861
3872 IF (m(1) < m(4) .AND. m(7) > 0 .AND. m(7) <= m(4) .AND.  &
         MOD(m(4)-m(1) , m(7)) == 0) GO TO 3874
     baddat =.true.
     GO TO 3875
3874 GO TO 3864
3875 IF (m1(1) == 0 .AND. m1(2) == 0) GO TO 3878
     km = 0
     kn = 0
     IF (naxf < 100) GO TO 3877
     CALL page2 (2)
     WRITE  (nout,3876) ufm,naxf
3876 FORMAT (a23,' 4125, MAXIMUM ALLOWABLE HARMONIC ID IS 99.  DATA ',  &
         'CONTAINS MAXIMUM =',i20)
     abort =.true.
3877 CONTINUE
     n  = n + 1
     i(n) =-1
     GO TO 9
3878 kn = 1
     km = km + 1
     GO TO 9
3890 IF (iaxf > 0) GO TO 3892
     IF (lharm) CALL page2 (2)
     IF (lharm) WRITE (nout,3891) ufm
3891 FORMAT (a23,' 4122, AXIF CARD REQUIRED.')
     lharm =.false.
     abort =.true.
3892 GO TO ret, (3105,3205,3310,3410,3510,4504,4610)
 
     !*******       213-BDYLIST        **************************************
 
3900 GO TO 3200
 
     !*******       214-FREEPT         **************************************
 
4000 GO TO 3400
 
     !*******       217-CTETRA,  335-CFTETRA  *******************************
 
4100 n = 6
     4105 DO  l = 1,n
         IF (m(l) <= 0) GO TO 8
     END DO
     n1 = n - 1
     DO  l = 3,n1
         l2 = l + 1
         DO  l1 = l2,n
             IF (m(l) == m(l1)) GO TO 8
         END DO
     END DO
     GO TO 3
 
     !*******       218-CWEDGE,  336-CFWEDGE  *******************************
 
4200 n = 8
     GO TO 4105
 
     !*******       219-CHEXA1,  333-CFHEX1   *******************************
 
4300 IF (mf(15) /= 0 .OR. mf(16) /= 0) GO TO 7
     n = 10
     GO TO 4105
 
     !*******       220-CHEXA2,  334-CFHEX2   *******************************
 
4400 IF (mf(15) /= 0 .OR. mf(16) /= 0) GO TO 7
     n = 10
     GO TO 4105
 
     !*******       221-DMIAX          **************************************
 
4500 IF (.NOT.fphys1) GO TO 4501
     fphys1 =.false.
     nm(1) = 0
     nm(2) = 0
4501 IF (km   /= 0) GO TO 4505
     IF (m(3) == 0) GO TO 4503
     IF (m(1) /= nm(1) .OR. m(2) /= nm(2)) GO TO 4510
     IF (mf(2) /= 1 .OR. mf(3) /= 1 .AND. mf(3) /= 0) GO TO 4511
     IF (mf(4) /= 1 .AND. mf(4) /= 0) GO TO 4511
     IF (m(3) <= 0 .OR. m(4) < 0 .OR. m(4) > 6) GO TO 4511
     IF (IABS(m(5)) > naxf) GO TO 4511
     IF (mf(5) /= 0 .OR. mf(6) /= 0 .OR. mf(7) /= 0 .OR. mf(8) /= 0) GO TO 4511
     IF (m1(1) /= 0 .OR. m1(2) /= 0) GO TO 4511
     n  = 2
     i(2) = m(4)
     l1 = 4
     l2 = 5
     ASSIGN 4595 TO ret
     GO TO 4520
4503 ASSIGN 4504 TO ret
     GO TO 3890
4504 IF (mf(1) /= 3 .OR. m(1) == nm(1) .AND. m(2) == nm(2)) GO TO 4510
     ifo = m(4)
     ty1 = m(5)
     ity1= 2*MOD(ty1,2)
     ty2 = m(6)
     IF (mach /= 12) GO TO 45045
     IF (ty2 == 2.OR.ty2 == 4) ty2 = ty2 - 1
45045 CONTINUE
      IF (ifo /= 1 .AND. ifo /= 2 .AND. ifo /= 6) GO TO 4510
      IF (ty1 <= 0 .OR. ty1 > 4 .OR. ty2 <= 0 .OR. ty2 > 4) GO TO 4510
      IF (ty2 == 1 .AND. ty1 == 3) GO TO 4510
      nm(1) = m(1)
      nm(2) = m(2)
      IF (mf(6) /= 0 .OR. mf(7) /= 0 .OR. mf(8) /= 0) GO TO 4511
      IF (m1f(2) /= 3 .OR. m1(3) /= nm(1) .OR. m1(4) /= nm(2)) GO TO 4511
      n = 9
      GO TO 3
4505  IF (m(1) <= 0 .OR. m(2) < 0 .OR. m(2) > 6) GO TO 4511
      IF (mf(1) /= 1 .OR. mf(2) /= 1 .AND. mf(2) /= 0) GO TO 4511
      IF (mf(4) /= 0 .AND. mf(4)+ity1 /= 4) GO TO 4511
      IF (mf(5) /= 0 .AND. ty1 /= 3 .AND. ty1 /= 4) GO TO 4511
      IF (IABS(m(3)) > naxf) GO TO 4511
      IF (mf(3) /= 1 .AND. mf(3) /= 0) GO TO 4511
      n  = 3
      i(2) = m(2)
      l1 = 3
      l2 = 3
      ASSIGN 4506 TO ret
      GO TO 4520
4506  i(3) = m(4)
      IF (ty1 == 1) GO TO 4508
      n = 4
      i(4) = m(5)
      IF (ty1 == 2 .OR. ty1 == 3) GO TO 4508
      n = 6
      i(5) = m(6)
      i(6) = m(7)
4508  IF (m1(1) == 0 .AND. m1(2) == 0) GO TO 4595
      n = n + 2
      i(n-1) =-1
      i(n  ) =-1
      IF (m1(1) == t1(1,k) .AND. m1(2) == t1(2,k) .AND. m1(3) == nm(1)  &
          .AND. m1(4) == nm(2)) GO TO 4592
      n = n + 2
      i(n-1) =-1
      i(n  ) =-1
      GO TO 4592
4510  nm(1) = m(1)
      nm(2) = m(2)
4511  abort =.true.
      CALL page2 (2)
      WRITE  (nout,4512) ufm,nm(1),nm(2)
4512  FORMAT (a23,' 4126, BAD DATA OR FORMAT OR NON-UNIQUE NAME, DMIAX',  &
          1X ,2A4)
      GO TO 4590
4520  IF (mf(l1) == 1) GO TO 4521
      i(1) = m(l2-2)
      GO TO 4525
4521  IF (m(l2) < 0) GO TO 4522
      i(1) = 1000000*(1+m(l2)) + m(l2-2)
      GO TO 4525
4522  i(1) = 500000*(1-m(l2)*2) + m(l2-2)
4525  GO TO ret, (4506,4595)
4590  IF (m1(1) == 0 .AND. m1(2) == 0) GO TO 4595
4592  km = 0
      kn = 0
      GO TO 2
4595  kn = 1
      km = km + 1
      GO TO 2
 
      !*******       222-FLSYM          **************************************
 
4600  IF (lflsym) GO TO 4690
      lflsym =.true.
      ASSIGN 4610 TO ret
      GO TO 3890
4610 CONTINUE
     IF (mf(1) /= 1 .OR. mf(2) /= 3 .OR. mf(3) /= 3) badfor =.true.
     DO  l = 4,8
         IF (mf(l) /= 0) badfor =.true.
     END DO
     IF (m(1) < 2 .OR. m(2) /= bcds .AND. m(2) /= bcda .OR.  &
         m(4) /= bcds .AND. m(4) /= bcda) baddat =.true.
     IF (MOD(m(1),2) /= 0) baddat =.true.
     IF (m(2) == bcds) m(2) = +1
     IF (m(2) == bcda) m(2) = -1
     IF (m(4) == bcds) m(3) = +1
     IF (m(4) == bcda) m(3) = -1
     n = 3
     GO TO 3
4690 CALL page2 (2)
     WRITE  (nout,4691) ufm
4691 FORMAT (a23,' 4123, ONLY ONE (1) FLSYM CARD ALLOWED IN BULK DATA')
     abort =.true.
     GO TO 2
 
     !*******     321-CEMLOOP     *******************************************
 
8000 IF (m(1) <= 0 .OR. m(13) < 0) GO TO 8
     IF (m(3) == 0) GO TO 8002
     IF (m(5) /= 0) GO TO 8
     DO  iem = 7,13
         IF (m(iem) /= 0) GO TO 8
     END DO
     GO TO 8003
8002 dx1 = z(4) - z(10)
     dy1 = z(5) - z(11)
     dz1 = z(6) - z(12)
     dx2 = z(7) - z(10)
     dy2 = z(8) - z(11)
     dz2 = z(9) - z(12)
     dl1 = dx1**2 + dy1**2 + dz1**2
     dl2 = dx2**2 + dy2**2 + dz2**2
     IF (ABS(dl1-dl2) > 1.e-4) GO TO 8
     dc1 = dy1*dz2 - dy2*dz1
     dc2 = dx2*dz1 - dx1*dz2
     dc3 = dx1*dy2 - dy1*dx2
     dlc = SQRT(dc1**2 + dc2**2 + dc3**2)
     IF (dlc/SQRT(dl2) < .0001) GO TO 8
8003 n = 13
     GO TO 3
 
     !*******    322-SPCFLD,   326-REMFLUX      *****************************
 
9000 IF (m(1) <= 0) GO TO 8
     IF (m(2) < 0) GO TO 8
     IF (m(6) /= -1) GO TO 9003
     DO  l = 7,8
         IF (mf(l) /= 0) GO TO 7
     END DO
     n = 6
     GO TO 3
9003 IF (mf(7) == 3) GO TO 9005
     DO  l = 6,8
         IF (mf(l) /= 1 .AND. mf(l) /= 0) GO TO 7
         IF (m(l) < 0) GO TO 8
         IF (m(l) == 0) CYCLE
         n = n + 6
         i(n-5) = m(1)
         i(n-4) = m(2)
         i(n-3) = m(3)
         i(n-2) = m(4)
         i(n-1) = m(5)
         i(n  ) = m(l)
     END DO
     IF (n > 0) THEN
         GO TO     2
     ELSE
         GO TO     8
     END IF
9005 IF (m(7) /= thru) GO TO 8
     IF (mf(6) /= 1 .OR. mf(8) /= 1) GO TO 7
     l1 = 6
     l2 = 9
     ii = m(l1) - 1
     l2 = m(l2) - m(l1)
     IF (ii < 0 .OR. l2 <= 0) GO TO 8
     l1 = 1
     DO  l = 1,5
         i(l) = m(l)
     END DO
     n = 6
     DO  l = l1,l2
         i(6) = l + ii
         CALL WRITE (209,i,n,0)
     END DO
     i(6) = ii + l2 + 1
     GO TO 2
 
     !*****   323-CIS2D8   **************************************************
 
9100 IF (m( 1) <= 0 .OR. m( 2) <= 0 ) GO TO 8
     IF (m(11) < 0 .OR. z(12) < 0.) GO TO 8
     IF (m(11) == 0) m(11) = 2
     IF (m(11) /= 2 .AND. m(11) /= 3) GO TO 8
     DO  l = 3,10
         IF (m(l) <= 0) GO TO 8
     END DO
     DO  l = 3,9
         lp1 = l + 1
         DO  lll = lp1,10
             IF (m(l) == m(lll)) GO TO 8
         END DO
     END DO
     n = 12
     GO TO 3
 
     !*****   324-PIS2D8   **************************************************
 
9200 IF (z(3) <= 0.) GO TO 8
     IF (m(1) <= 0 .OR. m(2) <= 0) GO TO 8
     n = 3
     GO TO 3
 
     !*****   325-GEMLOOP   *************************************************
 
9300 IF (mf(1) /= 1) GO TO 7
     IF (mf(2) /= 2 .AND. mf(2) /= 0) GO TO 7
     IF (mf(3) /= 1 .AND. mf(3) /= 0) GO TO 7
     IF (m(1) <= 0 .OR. m(3) < 0) GO TO 8
 
     !     FOR NOW, CID MUST BE 0
 
     IF (m(3) /= 0) GO TO 7
     npts = 0
     DO  l = 4,49,3
         IF (mf(l) == 3) GO TO 9320
         npts = npts + 1
         IF (mf(l  ) /= 2 .AND. mf(l  ) /= 0) GO TO 7
         IF (mf(l+1) /= 2 .AND. mf(l+1) /= 0) GO TO 7
         IF (mf(l+2) /= 2 .AND. mf(l+2) /= 0) GO TO 7
     END DO
     GO TO 8
9320 IF (npts < 2) GO TO 8
     DO  lll = l,49
         m(lll) = 0
     END DO
     DO  l = 1,3
         i(l) = m(l)
     END DO
     i(4) = npts
     DO  l = 4,48
         i(l+1) = m(l)
     END DO
     n = 49
     GO TO 2
 
     !*****   327-BFIELD   **************************************************
 
9400 IF (m(1) < 0) GO TO 8
     IF (m(2) /= -1) GO TO 9405
     DO  l = 3,8
         IF (mf(l) /= 0) GO TO 7
     END DO
     n = 2
     GO TO 3
9405 IF (mf(3) == 3) GO TO 9420
     DO  l = 2,8
         IF (mf(l) /= 1 .AND. mf(l) /= 0) GO TO 7
         IF (m(l) < 0) GO TO 8
         IF (m(l) == 0) CYCLE
         n = n + 2
         i(n-1) = m(1)
         i(n  ) = m(l)
     END DO
     IF (n > 0) THEN
         GO TO     2
     ELSE
         GO TO     8
     END IF
9420 IF (m(3) /= thru) GO TO 7
     IF (mf(2) /= 1 .OR. mf(4) /= 1) GO TO 7
     l1 = 2
     l2 = 5
     ii = m(l1) - 1
     l2 = m(l2) - m(l1)
     IF (ii < 0 .OR. l2 <= 0) GO TO 8
     l1 = 1
     i(1) = m(1)
     n  = 2
     DO  l = l1,l2
         i(2) = l + ii
         CALL WRITE (201,i,n,0)
     END DO
     i(2) = ii + l2 + 1
     GO TO 2
 
     !*****   328-MDIPOLE     ***********************************************
 
9500 IF (m(1) <= 0  .OR. m( 2) < 0 ) GO TO 8
     IF (z(9) < 0. .OR. z(10) < 0.) GO TO 8
     n = 10
     GO TO 3
 
 END SUBROUTINE ifs4p
