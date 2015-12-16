SUBROUTINE ifs3p (*,*,*)

 
 LOGICAL :: noud,nos,baddat,badfor,abort,lh,iax,idfreq,lharm,  &
     grdmsg,ifpdco,perm,prol,rbe,first,prt
!HURNB 11/93
 LOGICAL :: oneh,blankh
!HURNE
 INTEGER :: r,r1,g1,t3,t4,thru,arigid,brigid,crigid,drigid,  &
     erigid,frigid,blnk,endt,ia(6),ib(6),ic(6),ja(6),  &
     jb(6),jc(6),crtr,crba,crbe,q(92)
 DIMENSION       rm(50)
!HURNB 11/93
 DIMENSION nam(2),iones(4)
!HURNE
 CHARACTER (LEN=19) :: gcc
 CHARACTER (LEN=19) :: scc
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim,sfm
 COMMON /system/ nbuf,nout,abort,idummy(52),ithrml,dum21(21),ipiez
 COMMON /ifpdta/ id,n,k,kx,ky,i(100),m(100),mf(100),m1(100),  &
     m1f(100),kn,baddat,badfor,nopen,nparam,iax,nn,  &
     iaxf,naxf,lharm,knt,slotdf(5),gc(7),ll(6), &
!HURNB 11/93  &
     nns,oneh,blankh,iaxg
!HURNE
 COMMON /zzzzzz/ ibuff(1)
 COMMON /ifpx2 / t3(2,270)
 COMMON /ifpx3 / t4(2,270)
 COMMON /cifs3p/ grdmsg,la1,l7,km,l0,g1,lh,  &
     igdst2,igdst6,igdst7,igdst8,iddsf, idfreq,idrad,nvar,ids,jms,kms,lplf
 EQUIVALENCE     (m(1),rm(1)),(line,idummy(9)),(nbpw,idummy(37))
!HURNB 11/93
 EQUIVALENCE (xin,ixin)
!HURNE
 DATA prol,endt    / .false.,4HENDT/, perm  /.false.             /
 DATA first,prt    /  2*.true.     /
 DATA lud,lz,kk,ls /  4HUD  ,4HZ    , 4HK   ,4HS   /, nt1 /250   /
 DATA arigid/4HCRIG/, brigid/4HD1  /, crigid/4HD2  /, irigid /1  /
 DATA drigid/4HD3  /, mset  /4HMSET/, blnk  /4H    /, thru/4HTHRU/
 DATA erigid/4H1   /, frigid/4H2   /, ind /4HIN  /
 DATA crtr  /4HCRTR/, crba  /4HCRBA/, crbe  /4HCRBE/, ium /4HUM  /
 
 DATA scc   /'SORTED CARD COUNT =' /, gcc  /'GENERATED CARD    -'/
!HURNB 11/93
 DATA iscr1 /301/, iones/4*-1/, nam/4HIFS3,4HP   /
!HURNE
 
 IF (k > 100) GO TO 81
 GO TO (   100, 200,   5,   5,   5,   5,   5,   5,   5,   5,  &
     5,3980,4020,   5,   5,   5,1700,   5,   5,   5,  &
     5,   5,   5,   5,   5,   5,   5,2800,   5,   5,  &
     5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
     5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
     5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
     5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
     5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
     5,3960,4020,4060,   5,   5,   5,   5,   5,   5,  &
     5,3981,   5,   5,   5,   5,   5,   5,   5,   5  ),k
 81 IF (kx > 100) GO TO 82
 GO TO (     5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
     5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
     5,   5,4060,   5,   5,1260,   5,   5,   5,   5,  &
     1310,1310,   5,   5,   5,   5,   5,1380,1390,   5,  &
     5,   5,1430,1440,1450,1460,1470,1480,1490,1500,  &
     1500,1520,1530,1540,1550,1560,1560,   5,   5,   5,  &
     5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
     5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
     5,1820,1820,1820,1420,   5,   5,   5,   5,   5,  &
     5,   5,   5,   5,   5,   5,   5,   5,   5,   5  ),kx
 82 IF (ky > 100) GO TO 83
 GO TO (     5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
     5,   5,   5,   5,   5,3981,   5,   5,   5,   5,  &
     5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
     5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
     5,   5,   5,1400,   5,   5,   5,   5,   5,   5,  &
     5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
     5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
     5,   5,1415,1415,   5,   5,   5,   5,2010,   5,  &
     5,   5,   5,2060,2111,2030,2040,2030,   5,1410,  &
     5,   5,   5,   5,   5,   5,7300,7000,   5,   5  ),ky
 83 kz = ky - 100
 IF (kz > 53) GO TO 5
 IF (kz < 47 .OR. kz > 51 .OR. .NOT.first) GO TO 90
 first = .false.
 IF (.NOT.prt) GO TO 90
 CALL page1
 WRITE  (nout,85) uim
 85 FORMAT (a29,', CONVERSIONS OF RIGID ELEMENTS, CRROD, CRBAR, ',  &
     'CRTRPLT, CRBE1, AND CRBE2, TO CRIGDR, CRIGD2, OR CRIGD3',  &
     /5X,'ARE AS FOLLOWS (BLANK FIELDS MAY BE PRINTED AS ZEROS',  &
     '. CONTINUATION FIELDS ARE NOT PRINTED) -',/)
 line = 8
 90 GO TO (     5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
     5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
     5,   5,   5,   5,   5,   5,   5,   5,5100,5200,  &
     5,3980,   5,   5,   5,   5,   5,   5,   5,   5,  &
     5,   5,   5,   5,2920,3010,6000,6100,6300,7000,  &
     6400,6500,6600                                     ),kz
 5 CALL page2 (2)
 WRITE  (nout,6) sfm
 6 FORMAT (a25,' 322, ILLEGAL ENTRY TO IFS3P.')
 abort  = .true.
 RETURN 1
 7 badfor = .true.
 RETURN 1
 8 baddat = .true.
 RETURN 1
 3 DO  l = 1,n
   i(l) = m(l)
 END DO
 9 RETURN 3
 
!*******              1-GRID            ********************************
 
 100 IF (mf(2) == 0) m(2) = igdst2
 IF (mf(6) == 0) m(6) = igdst6
 IF (mf(7) == 0) m(7) = igdst7
 IF (mf(8) == 0) m(8) = igdst8
 IF (m(1) <= 0 .OR. m(2) < 0 .OR. m(6) < -1) GO TO 8
 IF (m(6) >= 0 .OR. grdmsg) GO TO 105
 CALL page2 (2)
 WRITE  (nout,103) uwm
 103 FORMAT (a23,' 302, ONE OR MORE GRID CARDS HAVE DISPLACEMENT ',  &
     'COORDINATE SYSTEM ID OF -1')
 grdmsg = .true.
 105 IF (ifpdco(m(7))) GO TO 8
 IF (ifpdco(m(8))) GO TO 8
 IF (mf(8) /= 0) GO TO 7
 n = 8
 GO TO 3
 
!*******        2-GRDSET       ****************************************
 
 200 IF (g1 == 0) GO TO 8
 g1 = 0
 IF (m(2) == 0 .AND. m(6) == 0 .AND. m(7) == 0 .AND. m(8) == 0) GO TO 8
 IF (m(2) < 0 .OR. m(6) < -1 .OR. m(7) < 0 .OR. m(8) < 0) GO TO 8
 IF (ifpdco(m(7)) .OR. ifpdco(m(8))) GO TO 8
 IF (mf(8) /= 0) GO TO 7
 igdst2 = m(2)
 igdst6 = m(6)
 igdst7 = m(7)
 igdst8 = m(8)
 RETURN 2
 
!*****         126-FREQ       ******************************************
 
 1260 IF (idfreq) iddsf = 0
 idfreq = .false.
 GO TO 1430
 
!******     131-RLOAD1, 132-RLOAD2    **********************************
 
 1310 IF (m(5) == 0 .AND. m(6) == 0) GO TO 8
 IF (m(1) <= 0 .OR. m(2) <= 0 .OR. m(3) < 0 .OR. m(4) < 0) GO TO 8
 IF (m(5) < 0 .OR. m(6) < 0) GO TO 8
 n = 6
 GO TO 3
 
!*******       138-TLOAD1      *****************************************
 
 1380 IF (m(1) <= 0 .OR. m(2) <= 0 .OR. m(3) < 0 .OR. m(5) <= 0) GO TO 8
 IF (m(4) < 0 .OR. m(4) > 4) GO TO 8
 n = 5
 GO TO 3
 
!*******       139-TLOAD2      *****************************************
 
 1390 IF (m(1) <= 0 .OR. m(2) <= 0 .OR. m(3) < 0) GO TO 8
 IF (rm(5) < 0. .OR. rm(6) <= rm(5) .OR. rm(7) < 0.) GO TO 8
 IF (m(4) < 0 .OR. m(4) > 4) GO TO 8
 n = 10
 GO TO 3
 
!******        244-RADMTX     *****************************************
 
 1400 IF (km == 1) GO TO 1431
 km = 1
 IF (mf(1) /= 1) badfor = .true.
 id = m(1)
 IF (id <= idrad) baddat = .true.
 idrad = id
 i(1)  = id
 n  = 1
 l1 = 2
 GO TO 1432
 
!******        290-VARIAN        **************************************
 
 1410 IF (km == 1) GO TO 1431
 km = 1
 IF (nvar /= 0) GO TO 8
 nvar = 1
 GO TO 1431
 
!*****      273-AEFACT , 274-FLFACT    ********************************
 
 1415 IF (km == 1) GO TO 1431
 km = 1
 IF (mf(1) /= 1) badfor = .true.
 IF (m(1)  <= 0) baddat = .true.
 i(1) = m(1)
 n  = 1
 l1 = 2
 IF (mf(3) /= 3) GO TO 1432
 IF (m(3) /= thru .OR. m(4) /= blnk) baddat = .true.
 IF (mf(2) /= 2 .OR. mf(4) /= 2 .OR. mf(5) /= 1 .OR. mf(6) /= 2)  &
     badfor = .true.
 IF (m(6) <=    1) baddat = .true.
 IF (m(5) == m(2)) baddat = .true.
 imid = 0
 IF (rm(5)-rm(7) >= 0. .AND. rm(7)-rm(2) < 0.) imid = 1
 IF (rm(5)-rm(7) <= 0. .AND. rm(7)-rm(2) > 0.) imid = 1
 IF (m1(1) /= 0 .OR. m1(2) /= 0) GO TO 1416
 badfor = .true.
 GO TO 1438
 1416 IF (badfor .OR. baddat) GO TO 1435
 IF (imid == 0) GO TO 1418
 rm(7) = 0.5*(rm(2) + rm(5))
 CALL page2 (3)
 WRITE  (nout,1417) uwm,i(1)
 1417 FORMAT (a25,' 528, FACTOR FMID IN FLFACT SET',i9,' DOES NOT LIE ',  &
     'BETWEEN F1 AND FNF.', /5X,'IT IS BEING RESET TO (F1 + ', 'FNF)/2.0')
 1418 t4(2,k) = t4(2,k) + 1
 CALL WRITE (204,i,1,0)
 l = 1
 1419 term1 = (m(6)-l)*(rm(5)-rm(7))
 term2 = (l-1)*(rm(7)-rm(2))
 anum  = rm(2)*term1 + rm(5)*term2
 den   = term1 + term2
 factor= anum/ den
 t4(2,k) = t4(2,k) + 1
 CALL WRITE (204,factor,1,0)
 l = l + 1
 IF (l <= m(6)) GO TO 1419
 i(1) = -1
 t4(2,k) = t4(2,k) + 1
 CALL WRITE (204,i,1,0)
 n  = 0
 km = 0
 kn = 0
 GO TO 9
 
!*****      143-DSFACT(1430), 185-PLFACT(1420)     ********************
 
 1420 IF (lplf < 0) THEN
   GO TO     8
 ELSE IF (lplf == 0) THEN
   GO TO  1425
 ELSE
   GO TO  1430
 END IF
 1425 lplf  = 1
 iddsf = 0
 1430 IF (km == 1) GO TO 1431
 km = 1
 IF (mf(1) /= 1) badfor = .true.
 id = m(1)
 IF (id <= iddsf) baddat = .true.
 iddsf = id
 i(1)  = id
 IF (mf(2) /= 2) badfor = .true.
 n  = 2
 l1 = 3
 i(n) = m(2)
 GO TO 1432
 1431 l1 = 1
 1432 DO  l = l1,8
   IF (mf(l) == 0) GO TO 1436
   IF (mf(l) /= 2) badfor = .true.
   n = n + 1
   i(n) = m(l)
 END DO
 IF (m1(1) == 0 .AND. m1(2) == 0) GO TO 1438
 1435 km = 0
 n  = n + 1
 i(n) =-1
 kn = 0
 GO TO 9
 1436 IF (l == 1) badfor = .true.
 DO  l2 = l,8
   IF (mf(l2) /= 0) badfor = .true.
 END DO
 IF (m1(1) /= 0 .OR. m1(2) /= 0) GO TO 1435
 badfor = .true.
 1438 kn = 1
 GO TO 9
 
!******        144-AXIC           **************************************
 
 1440 IF (iax) GO TO 1445
 iax = .true.
 nn  = 998
 DO  l = 1,nt1
   IF (t4(1,l) > 0) t3(1,l) = t3(1,k)
 END DO
!HURD2 11/93
!     IF (M(1).LT.0 .OR. M(1).GT.998 .OR. M(2).NE.0) GO TO 8
!     NN = M(1)
!HURNB 11/93
 
! M.LT.0 CHECK IS REMOVED TO ALLOW FOR SINGLE HARMONIC
 
!     IF(M(1).LT.0.OR.M(1).GT.998.OR.M(2).NE.0)GO TO 8
 IF(             m(1) > 998.OR.m(2) /= 0)GO TO 8
 nns = m(1)
 nn = IABS(m(1))
 oneh = .false.
 IF(nns < 0)oneh = .true.
!HURNE
 n  = 2
 IF (nn > 15 .AND. nbpw <= 32) GO TO 1448
 GO TO 3
 1445 CALL page2 (2)
 WRITE  (nout,1446) ufm
 1446 FORMAT (a23,' 329, ONLY ONE(1) AXIC CARD ALLOWED.')
 abort = .true.
 GO TO 2
 1448 WRITE  (nout,1449) uwm
 1449 FORMAT (a25,', POTENTIAL SYSTEM FATAL ERROR DUE TO LARGE HARMONIC'  &
     ,      ' (LARGER THAN 15) ON 32-BIT WORD MACHINE')
 GO TO 3
!  OR GO TO 1447
 
!******        145-RINGAX         **************************************
 
 1450 IF (m(1) <= 0 .OR. rm(3) <= 0.) GO TO 8
 ih = nn
 ASSIGN 1451 TO r
 ASSIGN    8 TO r1
 GO TO 21
 1451 IF (ifpdco(m(7))) GO TO 8
 n = 4
 i(1) = m(1)
 i(2) = m(3)
 i(3) = m(4)
 i(4) = m(7)
 GO TO 2
 
!******        146-CCONEAX        **************************************
 
 1460 IF (m(1) <= 0 .OR. m(3) <= 0 .OR. m(4) <= 0) GO TO 8
 IF (mf(2) == 0) m(2) = m(1)
 IF (m(2) <= 0 .OR. m(4) == m(3)) GO TO 8
 ih = nn
 ASSIGN 1461 TO r
 ASSIGN    8 TO r1
 GO TO 21
 1461 n = 4
 GO TO 3
 
!******        147-PCONEAX        **************************************
 
 1470 IF (m(1) <= 0) GO TO 8
 IF (m(2) == 0 .AND. m(3) /= 0 .OR. m(2) < 0) GO TO 8
 IF (m(4) == 0 .AND. m(5) /= 0 .OR. m(4) < 0) GO TO 8
 IF (m(6) == 0 .AND. m(7) /= 0 .OR. m(6) < 0) GO TO 8
 IF (m(2) /= 0 .AND. m(3) == 0) GO TO 8
 IF (m(6) /= 0 .AND. m(7) == 0) GO TO 8
 ih = nn
 ASSIGN 1471 TO r
 ASSIGN    8 TO r1
 GO TO 21
 1471 n = 24
 GO TO 3
 
!******        148-SPCAX       *****************************************
 
 1480 IF (m(1) <= 0 .OR. m(2) <= 0 .OR. m(3) < 0) GO TO 8
 IF (ifpdco(m(4))) GO TO 8
!HURNB 11/93
 IF(mf(3) == 0)GO TO 1481
!HURNE
 ASSIGN    8 TO r1
 ASSIGN 1489 TO r
 ih = m(3)
 GO TO 21
 1489 n = 5
 GO TO 3
!HURNB 11/93
 
! HID IS BLANK - GENERATE HID FOR THIS SPCAX FOR ALL HARMONICS
 
 1481 nharms=nns+1
 IF(oneh)nharms=1
 DO  il=1,nharms
   n=n+5
   i(n-4)=m(1)
   i(n-3)=m(2)
   i(n-1)=m(4)
   i( n )=m(5)
   i(n-2)=il-1
   IF(oneh)i(n-2)=nn
 END DO
 GO TO 2
!HURNE
 
!******        149-MPCAX       *****************************************
 
 1490 IF (m(7) > 6) baddat = .true.
 IF (ithrml == 1 .AND. m(7) > 1) baddat = .true.
 IF (km /= 0) GO TO 1492
 km = 1
 nt = 0
!HURNB 11/93
 blankh=.false.
!HURNE
 IF (mf(1) /= 1 .OR. mf(2) /= 0 .OR. mf(3) /= 0 .OR. mf(4) /= 0)  &
     badfor = .true.
 l1 = 5
!HURNB 11/93
 IF(mf(6) == 0)blankh=.true.
 IF(blankh)CALL gopen(iscr1,ibuff(2*nbuf+1),1)
!HURNE
 ASSIGN 1491 TO r
 GO TO 1493
 1491 IF (m(1) <= 0) baddat = .true.
 id = m(1)
 n  = 1
 i(n) = id
 ih = nn
 ASSIGN 1497 TO r
 ASSIGN    8 TO r1
 GO TO 21
 1492 l1 = 1
 IF (m(3) > 6) baddat = .true.
 IF (ithrml == 1 .AND. m(3) > 1) baddat = .true.
 ASSIGN 1496 TO r
 1493 DO  l = l1,8
   IF (mf(l) == 0) CYCLE
   IF (l == 4 .OR. l == 8) GO TO 1494
   IF (mf(l) /= 1) badfor = .true.
   CYCLE
   1494 IF (mf(l) /= 2) badfor = .true.
 END DO
 GO TO r, (1491,1496)
 1496 n = 0
 1497 DO  l = l1,5,4
   IF (m(l  ) == 0 .AND. m(l+1) == 0 .AND. m(l+2) == 0 .AND.  &
       m(l+3) == 0) CYCLE
   IF (m(l) <= 0 .OR. m(l+1) < 0 .OR. m(l+2) < 0 .OR. m(l+3) == 0  &
       .AND. l1 == 5) baddat = .true.
!HURNB 11/93
   IF(blankh.AND.l1 == 1.AND.mf(l+1) /= 0)badfor=.true.
   IF(.NOT.blankh.AND.l1 == 1.AND.mf(l+1) == 0)badfor=.true.
!HURNE
   n = n + 4
   i(n-3) = m(l  )
   i(n-2) = m(l+1)
   i(n-1) = m(l+2)
   i(n  ) = m(l+3)
 END DO
 nt = nt + n
 IF (n < 4) baddat = .true.
 kn = 1
!HURNB 11/93
 IF(m1(1) /= 0.OR.m1(2) /= 0)GO TO 1499
 IF(.NOT.blankh)GO TO 9
 CALL WRITE(iscr1,i,n,0)
!      WRITE(6,10005)N,(I(IL),IL=1,N)
!10005 FORMAT(6H MPCAX,6I5)
 n=0
 GO TO 9
!HURNE
 IF (m1(1) == 0 .AND. m1(2) == 0) GO TO 9
!HURNB 11/93
 1499 CONTINUE
!HURNE
 n  = n + 4
 i(n-3) = -1
 i(n-2) = -1
 i(n-1) = -1
 i(n  ) = -1
 kn = 0
 km = 0
 IF (nt < 9) baddat = .true.
!HURNB 11/93
 IF(.NOT.blankh)GO TO 9
 
! MPCAX CARD DONE - GENERATE CARDS FOR ALL HARMONICS ASSUMING THE ONE JU
! STORED (WITH BLANK HARMONIC) IS FOR THE ZERO HARMONIC
 
 IF(nt > nopen)CALL mesage(-8,0,nam)
 CALL WRITE(iscr1,i,n-4,1)
!     WRITE(6,10006)N,(I(IL),IL=1,N)
!0006 FORMAT(7H MPCAX1,10I5)
 CALL CLOSE(iscr1,1)
 CALL gopen(iscr1,ibuff(2*nbuf+1),0)
 CALL READ(*14990,*14991,iscr1,ibuff(3*nbuf+1),nopen,0,nnt)
 14990 CALL mesage(-8,0,nam)
 14991 CALL CLOSE(iscr1,1)
!     WRITE(6,10007)NT,NNT,(IBUFF(3*NBUF+IL),IL=1,NNT)
!0007 FORMAT(7H MPCAX2,10I5)
 IF(nt /= nnt)CALL mesage(-61,0,0)
 
! ALL MPCAX CARD INFO FOR THIS CARD IS READ IN. GENERATE FOR ALL HARMONI
 
 nharms=nns+1
 IF(oneh)nharms=1
 DO  l=1,nharms
   ill=l-1
   IF(oneh)ill=IABS(nns)
   DO  il=3,nt,4
     ibuff(3*nbuf+il)=ill
   END DO
   t4(2,k)=t4(2,k)+nt
   CALL WRITE(215,ibuff(3*nbuf+1),nt,0)
   t4(2,k)=t4(2,k)+4
   CALL WRITE(215,iones,4,0)
 END DO
 n=0
!HURNE
 GO TO 9
 
!******        151-SUPAX, 150-OMITAX     *******************************
 
 1500 l = 1
 1501 IF (m(l) == 0 .AND. m(l+1) == 0 .AND. m(l+2) == 0) GO TO 1510
 IF (m(l) <= 0 .OR.  m(l+1) < 0) GO TO 8
 IF (ifpdco(m(l+2))) GO TO 8
 ASSIGN 1507 TO r
 ASSIGN    8 TO r1
 ih = m(l+1)
 GO TO 21
 1507 n  = n + 3
 IF (n > 3 .AND. m(l) == m(l-3) .AND. m(l-1) == m(l-4) .AND.  &
     m(l-2) == m(l-5)) GO TO 8
 i(n-2) = m(l  )
 i(n-1) = m(l+1)
 i(n  ) = m(l+2)
 1510 l = l + 3
 IF (l == 4) GO TO 1501
 IF (n > 0) THEN
   GO TO     2
 ELSE
   GO TO     8
 END IF
 
!******        152-POINTAX        **************************************
 
 1520 n = 3
 1521 IF (m(1) <= 0 .OR. m(2) <= 0) GO TO 8
 1522 ASSIGN 3 TO r
 1523 ih = nn
 ASSIGN 8 TO r1
 GO TO 21
 1524 ASSIGN 2 TO r
 GO TO 1523
 
!******        153-SECTAX         **************************************
 
 1530 n = 5
 IF (rm(3) > 0.0) THEN
   GO TO  1521
 ELSE
   GO TO     8
 END IF
 
!******        154-PRESAX         **************************************
 
 1540 n = 6
 IF (m(1) <= 0 .OR. m(4) <= 0 .OR. m(4) == m(3)) GO TO 8
 IF (ipiez == 1) GO TO 1522
 IF (m(3)  <= 0) GO TO 8
 IF (ABS(rm(5)) >= ABS(rm(6)) .AND. SIGN(1.,rm(5)) == SIGN(1.,rm(6  &
     ))) GO TO 8
 GO TO 1522
 
!******        155-TEMPAX         **************************************
 
 1550 DO  l = 1,5,4
   IF (m(l) == 0 .AND. m(l+1) == 0 .AND. m(l+2) == 0) CYCLE
   IF (m(l) <= 0 .OR.  m(l+1) <= 0) GO TO 8
   n = n + 4
   i(n-3) = m(l  )
   i(n-2) = m(l+1)
   i(n-1) = m(l+2)
   i(n  ) = m(l+3)
 END DO
 IF (n > 0) THEN
   GO TO  1524
 ELSE
   GO TO     8
 END IF
 
!******     156-FORCEAX, 157-MOMAX    *******************************
 
 1560 IF (m(1) <= 0 .OR. m(2) <= 0) GO TO 8
 IF (mf(3) == 2.OR.mf(3) == 4) GO TO 8
 IF (mf(3) /= 3 .AND. m(3) < 0) GO TO 8
 n = 8
 l = 4
 i(1) = m(1)
 i(2) = m(2)
 i(3) = m(3)
 i(4) = 0
 IF (mf(3) == 3) i(4) = m(4)
 IF (mf(3) == 3) l = 5
 i(5) = m(l)
 i(6) = m(l+1)
 i(7) = m(l+2)
 i(8) = m(l+3)
 GO TO 2
 
!******        17-MPC       ******************************************
 
 1700 IF (m(3) > 6 .OR. m(6) > 6) baddat = .true.
 IF (ithrml /= 1) GO TO 1710
 IF (m(3) > 1 .OR. m(6) > 1) baddat = .true.
 1710 IF (km /= 0) GO TO 1724
 km = 1
 nt = 0
 IF (mf(1) /= 1 .OR. mf(8) /= 0) badfor = .true.
 ASSIGN 1712 TO r
 GO TO 1725
 1712 IF (m(1) <= 0) baddat = .true.
 id = m(1)
 IF (m(2) <= 0 .OR. m(3) < 0 .OR. m(4) == 0) baddat = .true.
 IF (ids == id .AND. jms == m(2) .AND. kms == m(3)) baddat = .true.
 ids = id
 jms = m(2)
 kms = m(3)
 n = 4
 DO  l = 1,4
   i(l) = m(l)
 END DO
 GO TO 1745
 1724 IF (mf(1) /= 0 .OR. mf(8) /= 0) badfor = .true.
 ASSIGN 1737 TO r
 1725 DO  l = 2,7
   IF (mf(l) == 0) CYCLE
   IF (l == 4 .OR. l == 7) GO TO 1733
   IF (mf(l) /= 1) badfor = .true.
   CYCLE
   1733 IF (mf(l) /= 2) badfor = .true.
 END DO
 GO TO r, (1712,1737)
 1737 n = 0
 IF (m(2) == 0 .AND. m(3) == 0 .AND. m(4) == 0) GO TO 1745
 IF (m(2) <= 0 .OR.  m(3) < 0) baddat = .true.
 n = 3
 DO  l = 2,4
   i(l-1) = m(l)
 END DO
 1745 IF (m(5) == 0 .AND. m(6) == 0 .AND. m(7) == 0) GO TO 1751
 IF (m(5) <= 0 .OR.  m(6) < 0) baddat = .true.
 n = n + 3
 i(n-2) = m(5)
 i(n-1) = m(6)
 i(n  ) = m(7)
 1751 IF (n <= 0) baddat = .true.
 nt = nt + n
 DO  l = 1,8
   m(l) = 0
 END DO
 kn = 1
 IF (m1(1) == 0 .AND. m1(2) == 0) GO TO 9
 n = n + 3
 i(n-2) = -1
 i(n-1) = -1
 i(n  ) = -1
 kn = 0
 km = 0
 IF (nt < 7) baddat = .true.
 GO TO 9
 
!******     182-DAREA, 183-DELAY, 184-DPHASE      *******************
 
 1820 IF (m(1) <= 0) GO TO 8
 DO  l = 2,5,3
!HURNB 11/93
!     WRITE(6,10003)L,M(L),M(L+1),M(L+2),N,NNS,(I(IL),IL=1,N)
!0003 FORMAT(7H DAREA0,6I10/(1X,24I5))
!HURNE
   IF (m(l) == 0 .AND. m(l+1) == 0 .AND. m(l+2) == 0) CYCLE
   IF (m(l) <= 0 .OR.  m(l+1) < 0 .OR.  m(l+1) > 6) GO TO 8
   n = n + 4
   i(n-3) = m(1)
   i(n-2) = m(l)
   i(n-1) = m(l+1)
   i(n  ) = m(l+2)
!HURNB 11/93
   IF(.NOT.iax)CYCLE
   IF(m(l) >= 1000000)CYCLE
   
! FOR AXIC PROBLEMS AND GRID ID ON DAREA .LT. 10**6, GENERATE DAREAS FOR
! HARMONICS, COMPUTING THE GRID ID.  ASSUME  PRESSURE VALUE IS GIVEN FOR
! ZERO HARMONIC; FOR HIGHER HARMONICS, HALVE IT.
   
   nharms=nns+1
   IF(oneh)nharms=1
   DO  il=1,nharms
     ill=il
     IF(nns >= 0 .AND. il == 1)GO TO 1823
     IF(il > 1)GO TO 1821
     
! NNS.LT.0 .AND. IL.EQ.1
     
     ill=nn+1
     GO TO 1822
     1821 n=n+4
     i(n-3)=m(1)
     i(n-1)=m(l+1)
     1822 xin=0.5*rm(l+2)
     i(n)=ixin
     1823 i(n-2)=m(l)+1000000*ill
   END DO
!HURNE
 END DO
!HURNB 11/93
!      WRITE(6,10001)NHARMS,NNS,N,(I(IL),IL=1,N)
!10001 FORMAT(6H DAREA,3I10/(1X,24I5))
!HURNE
 IF (n > 0) THEN
   GO TO     2
 ELSE
   GO TO     8
 END IF
 
!******        279-CRIGD1         **********************************
 
 2010 CONTINUE
 kn = 1
 SELECT CASE ( irigid )
   CASE (    1)
     GO TO 2011
   CASE (    2)
     GO TO 2012
 END SELECT
 2011 CONTINUE
 irigid = irigid + 1
 IF (mf(1) /= 1) badfor = .true.
 IF (m(1) <= 0) baddat = .true.
 i(1) = m(1)
 n = 2
 IF (mf(2) /= 1) badfor = .true.
 IF (m(2) < 1) baddat = .true.
 i(2) = m(2)
 IF (mf(4) == 3) GO TO 2020
 irg = 3
 GO TO 2013
 2012 CONTINUE
 n = 0
 irg = 1
 2013 CONTINUE
 DO  l = irg,8
   l1 = l
   IF (m(l)  <= 0) GO TO 2018
   IF (mf(l) /= 1) badfor = .true.
   i(n+1) = m(l)
   DO  j = 1,6
     i(n+1+j) = j
   END DO
   n = n + 7
 END DO
 IF (m1(1) == 0 .AND. m1(2) == 0) GO TO 9
 2016 irigid = 1
 DO  j = 1,7
   i(n+j) = -1
 END DO
 IF (m1(1) == arigid .AND. m1(2) == brigid) i(n+2) = 0
 n  = n + 7
 kn = 0
 GO TO 9
 2018 CONTINUE
 DO  lk = l1,8
   IF (m(lk) /= 0) baddat = .true.
   IF (mf(lk) /= 0) badfor = .true.
 END DO
 IF (m1(1) == 0 .AND. m1(2) == 0) GO TO 2022
 GO TO 2016
 2020 IF (m(4) == thru .AND. m(5) == blnk) GO TO 2024
 baddat = .true.
 IF (m1(1) /= 0 .OR. m1(2) /= 0) GO TO 2016
 2022 badfor = .true.
 GO TO 9
 2024 IF (m1(1) == 0 .AND. m1(2) == 0) GO TO 2022
 IF (mf(3) /= 1 .OR.mf(5) /= 1) badfor = .true.
 IF (m(3) <= 0  .OR. m(6) <= 0) baddat = .true.
 IF (m(6) <= m(3)) baddat = .true.
 DO  l = 6,8
   IF (mf(l) /= 0) badfor = .true.
 END DO
 IF (badfor .OR. baddat) GO TO 2016
 t4(2,k) = t4(2,k) + 2
 CALL WRITE (210,m,2,0)
 l = m(3)
 2026 i(1) = l
 DO  j = 1,6
   i(j+1) = j
 END DO
 t4(2,k) = t4(2,k) + 7
 CALL WRITE (210,i,7,0)
 l = l + 1
 IF (l <= m(6)) GO TO 2026
 irigid = 1
 DO  j = 1,7
   i(j) = -1
 END DO
 IF (m1(1) == arigid .AND. m1(2) == brigid) i(2) = 0
 n  = 0
 kn = 0
 t4(2,k) = t4(2,k) + 7
 CALL WRITE (210,i,7,0)
 GO TO 9
 
!******        284-CRIGD2        **********************************
 
 2060 CONTINUE
 kn = 1
 SELECT CASE ( irigid )
   CASE (    1)
     GO TO 2061
   CASE (    2)
     GO TO 2062
 END SELECT
 2061 irigid = irigid + 1
 IF (mf(1) /= 1) badfor = .true.
 IF (m(1)  <= 0) baddat = .true.
 i(1) = m(1)
 n = 2
 IF (mf(2) /= 1) badfor = .true.
 IF (m(2)  < 1) baddat = .true.
 i(2) = m(2)
 irg  = 3
 GO TO 2063
 2062 CONTINUE
 n   = 0
 irg = 1
 2063 CONTINUE
 DO   l = irg,8,2
   l1 = l
   IF (m(l  ) <= 0) GO TO 2068
   IF (m(l+1) <= 0) baddat = .true.
   IF (mf(l) /= 1 .OR. mf(l+1) /= 1) badfor = .true.
   i(n+1) = m(l)
   IF (ifpdco(m(l+1))) baddat = .true.
   DO  j = 1,6
     i(n+1+j) = ll(j)
   END DO
   n = n + 7
 END DO
 IF (m1(1) == 0 .AND. m1(2) == 0) GO TO 9
 2066 irigid = 1
 DO  j = 1,7
   i(n+j) = -1
 END DO
 IF (m1(1) == arigid .AND. m1(2) == crigid) i(n+2) = 0
 n  = n + 7
 kn = 0
 GO TO 9
 2068 CONTINUE
 IF (m1(1) == 0 .AND. m1(2) == 0) baddat = .true.
 DO  lk = l1,8
   IF (m(lk)  /= 0) baddat = .true.
   IF (mf(lk) /= 0) badfor = .true.
 END DO
 GO TO 2066
 
!******      298-CRIGD3, 350-CRBE1       ******************************
 
 7000 kn = 1
 SELECT CASE ( irigid )
   CASE (    1)
     GO TO 7020
   CASE (    2)
     GO TO 7160
   CASE (    3)
     GO TO 7200
   CASE (    4)
     GO TO 7240
 END SELECT
 7020 irigid = 2
 jrigid = 1
 knt1 = knt
 l1 = 2
 l2 = 6
 l6 = 0
 IF (mf(1) /= 1 .OR. mf(2) /= 1 .OR. mf(3) /= 1) badfor = .true.
 IF (m(1) < 1 .OR. m(2) < 1 .OR. m(3) < 1) baddat = .true.
 n  = 1
 i(1) = m(1)
 q(1) = m(1)
 l8   = 1
 ncomp= 0
 7040 l5 = l2 + 2
 DO  l = l1,l2,2
   l3 = l + 1
   IF (mf(l-l6) == 0) GO TO 7120
   IF (mf(l-l6) /= 1 .OR. mf(l-l6+1) /= 1) badfor = .true.
   IF (m(l) < 1 .OR. m(l+1) < 1) baddat = .true.
   IF (.NOT.prt) GO TO 7050
   q(l8+1) = m(l )
   q(l8+2) = m(l3)
   l8 = l8 + 2
   7050 i(n+1) = m(l)
   IF (ifpdco(m(l+1))) baddat = .true.
   DO  j = 1,6
     i(n+j+1) = ll(j)
     IF (irigid == 4) CYCLE
     IF (ll(j) /= 0) ncomp = ncomp + 1
   END DO
   n = n + 7
   IF (irigid == 4) CYCLE
   IF (ncomp  > 6) baddat = .true.
 END DO
 IF (mf(l5-l6) /= 0) badfor = .true.
 IF (m(l5) /= 0) baddat = .true.
 SELECT CASE ( jrigid )
   CASE (    1)
     GO TO 7100
   CASE (    2)
     GO TO 7100
   CASE (    3)
     GO TO 7220
 END SELECT
 7100 IF (m1(1) /= 0  .OR.  m1(2) /= 0) GO TO 7110
 IF (m1f(2) /= 0 .AND. ncomp < 6) baddat = .true.
 GO TO 9
 7110 badfor = .true.
 irigid = 1
 GO TO 9
 7120 DO  lk = l3,l5
   IF (mf(lk-l6) /= 0) badfor = .true.
   IF (m(lk) /= 0) baddat = .true.
 END DO
 SELECT CASE ( jrigid )
   CASE (    1)
     GO TO 7100
   CASE (    2)
     GO TO 7100
   CASE (    3)
     GO TO 7220
 END SELECT
 
 7160 IF (mf(1) /= 0) GO TO 7200
 irigid = 3
 jrigid = 2
 l1 = 2
 l2 = 6
 l6 = 0
 IF (mf(2) /= 1 .OR. mf(3) /= 1) badfor = .true.
 IF (m(1) /= 0 .OR. m(2) < 1 .OR. m(3) < 1) baddat = .true.
 n = 0
 GO TO 7040
 
 7200 irigid = 4
 jrigid = 3
 l1 = 3
 l2 = 7
 l6 = 1
 l7 = l8
 IF (mf(1) /= 3 .OR. mf(2) /= 1 .OR. mf(3) /= 1) badfor = .true.
 IF ((m(1) /= mset .AND. m(1) /= ium) .OR. m(2) /= blnk .OR.  &
     m(3) < 1 .OR. m(4) < 1) GO TO 7250
 n = 1
 i(1) = mset
 GO TO 7040
 7220 IF (m1(1) == 0 .AND. m1(2) == 0) GO TO 9
 irigid = 1
 DO  j = 1,7
   i(n+j) = -1
 END DO
 IF (m1(1) == arigid .AND. m1(2) == drigid) i(n+2) = 0
 IF (m1(1) == crbe   .AND. m1(2) == erigid) i(n+2) = 0
 n  = n + 7
 kn = 0
 IF (kz /= 50 .OR. .NOT.prt) GO TO 9
 lk = (l8+4)/3 + 2
 CALL page2 (lk)
 WRITE  (nout,7232) scc,knt1,(q(j),j=1,l7)
 7232 FORMAT (/25X,a19,i7,1H-,5X,'CRBE1 ',7I8, /,(71X,6I8))
 lk = l7 + 1
 WRITE  (nout,7234) (q(j),j=lk,l8)
 7234 FORMAT (69X,'UM',6I8, /,(71X,6I8))
 WRITE  (nout,7236) gcc,(q(j),j=1,l7)
 7236 FORMAT (25X,a19,13X,'CRIGD3',7I8, /,(71X,6I8))
 WRITE  (nout,7238) (q(j),j=lk,l8)
 7238 FORMAT (67X,'MSET',6I8, /,(71X,6I8))
 GO TO 9
 
 7240 l1 = 2
 l2 = 6
 l6 = 0
 IF (mf(1) /= 0 .OR. mf(2) /= 1 .OR. mf(3) /= 1) badfor = .true.
 IF (m(1) /= 0 .OR. m(2) < 1 .OR. m(3) < 1) baddat = .true.
 n = 0
 GO TO 7040
 7250 WRITE (nout,6475) ufm,blnk,q(1),knt1
 GO TO 8
 
!******         297-CRIGDR      *************************************
 
 7300 DO  l = 1,5,4
   IF (m(l  ) == 0 .AND. m(l+1) == 0 .AND. m(l+2) == 0 .AND.  &
       m(l+3) == 0) CYCLE
   IF (m(l) <= 0 .OR. m(l+1) <= 0 .OR. m(l+2) <= 0 .OR. m(l+3) <= 0) GO TO 8
   IF (m(l+1) == m(l+2)) GO TO 8
   IF (m(l+3) > 3) GO TO 7310
   n = n + 4
   IF (n > 4 .AND. m(l) == m(l-4)) GO TO 8
   i(n-3) = m(l  )
   i(n-2) = m(l+1)
   i(n-1) = m(l+2)
   i(n  ) = m(l+3)
   CYCLE
   7310 WRITE (nout,6475) ufm,blnk,m(l),knt
   baddat = .true.
 END DO
 IF (n > 0) THEN
   GO TO     2
 ELSE
   GO TO     8
 END IF
 
!******       347-CRROD        *****************************************
 
!     MAP THIS RIGID ELEMENT INTO CRIGID3 FORM
 
 6000 IF (mf(1)+mf(2)+mf(3) /= 3) GO TO 7
 IF (m(1) <= 0 .OR. m(2) <= 0 .OR. m(3) <= 0) GO TO 8
 IF (m(2) == m(3)) GO TO 8
 IF (m(4) < 0 .OR. m(5) < 0) GO TO 6480
 l = m(4) + m(5)
 IF (l < 1 .OR. l > 3) GO TO 6480
 IF (m(4) /= 0 .AND. m(5) /= 0) GO TO 6480
 IF (.NOT.prt) GO TO 6004
 CALL page2 (3)
 IF (m(4) /= 0) WRITE (nout,6002) scc,knt,(m(j),j=1,4),gcc,m(1),  &
     m(3),m(2),m(4)
 IF (m(4) == 0) WRITE (nout,6003) scc,knt,(m(j),j=1,3),m(5),gcc,  &
     (m(j),j=1,3),m(5)
 6002 FORMAT (/25X,a19,i7,1H-,5X,'CRROD ',4I8, /25X,a19,13X,'CRIGDR',4I8)
 6003 FORMAT (/25X,a19,i7,1H-,5X,'CRROD ',3I8,8X,i8,  &
     /25X,a19,13X,'CRIGDR',4I8)
 6004 l = m(3)
 IF (m(4) == 0) GO TO 6005
 l = m(2)
 m(2) = m(3)
 m(3) = l
 m(5) = m(4)
 6005 m(4) = m(5)
 n = 4
 GO TO 3
 
!******       348-CRBAR        *****************************************
 
!     MAP THIS RIGID ELEMENT INTO CRIGD3 FORM
 
 6100 IF (mf(1)+mf(2)+mf(3) /= 3) GO TO 7
 IF (m(1) <= 0 .OR. m(2) <= 0 .OR. m(3) <= 0) GO TO 8
 IF (m(2) == m(3)) GO TO 8
 rbe = .false.
 IF (m(6) == 0 .AND. m(7) == 0) rbe = .true.
 IF (m(4) == 0 .AND. m(5) == 0) GO TO 6470
 IF (ifpdco(m(4))) GO TO 6470
 lk = 1
 DO  l = 1,6
   lll = ll(l)
   IF (rbe .AND. lll == 0) m(6) = m(6) + l*lk
   IF (lll == 0) lk = lk*10
   ia(l) = lll
 END DO
 IF (ifpdco(m(5))) GO TO 6470
 lk = 1
 DO  l = 1,6
   lll = ll(l)
   IF (rbe .AND. lll == 0) m(7) = m(7) + l*lk
   IF (lll == 0) lk = lk*10
   ib(l) = lll
 END DO
 IF (rbe) GO TO 6130
 IF (ifpdco(m(6))) GO TO 6480
 DO  l = 1,6
   IF (ia(l) ==     0) GO TO 6120
   IF (ia(l) == ll(l)) GO TO 6480
   6120 ja(l) = ll(l)
 END DO
 IF (ifpdco(m(7))) GO TO 6480
 DO  l = 1,6
   IF (ib(l) ==     0) GO TO 6125
   IF (ib(l) == ll(l)) GO TO 6480
   6125 jb(l) = ll(l)
 END DO
 
 6130 IF (.NOT.prt) GO TO 6133
 CALL page2 (4)
 WRITE  (nout,6131) scc,knt,(m(l),l=1,7),gcc,m(1),m(2),m(4),m(3),  &
     m(5),m(2),m(6),m(3),m(7)
 6131 FORMAT (/25X,a19,i7,1H-,5X,'CRBAR ',7I8,  &
     /25X,a19,13X,'CRIGD3',5I8, /67X,'MSET',4I8)
 
!     KZ=48 (CRBAR),   KZ=49 (CRTRPLT)
 
 6133 ncomp = 0
 DO  l = 1,6
   IF (ia(l) /= 0) ncomp = ncomp + 1
   IF (ib(l) /= 0) ncomp = ncomp + 1
   IF (kz   /= 49) CYCLE
   IF (ic(l) /= 0) ncomp = ncomp + 1
 END DO
 IF (ncomp /= 6) GO TO 6470
 lk = 0
 IF (kz == 49) lk = 1
 i(1) = m(1)
 n  = 2
 IF (m(4+lk) == 0) GO TO 6143
 i(n) = m(2)
 DO  j = 1,6
   i(n+j) = ia(j)
 END DO
 n = n + 7
 6143 IF (m(5+lk) == 0) GO TO 6147
 i(n) = m(3)
 DO  j = 1,6
   i(n+j) = ib(j)
 END DO
 n = n + 7
 6147 IF (kz /= 49 .OR. m(6+lk) == 0) GO TO 6160
 i(n) = m(4)
 DO  j = 1,6
   i(j+n) = ic(j)
 END DO
 n = n + 7
 
 6160 i(n) = mset
 n = n + 1
 IF (.NOT.rbe) GO TO 6170
 DO  j = 1,6
   IF (ia(j) == 0) ia(j) =-j
   IF (ia(j) > 0) ia(j) = 0
   IF (ib(j) == 0) ib(j) =-j
   IF (ib(j) > 0) ib(j) = 0
   IF (kz   /= 49) CYCLE
   IF (ic(j) == 0) ic(j) =-j
   IF (ic(j) > 0) ic(j) = 0
 END DO
 6170 IF (kz == 49) lk = 3
 IF (m(6+lk) == 0) GO TO 6177
 i(n) = m(2)
 DO  j = 1,6
   IF (     rbe) i(n+j) =-ia(j)
   IF (.NOT.rbe) i(n+j) = ja(j)
 END DO
 n = n + 7
 6177 IF (m(7+lk) == 0) GO TO 6182
 i(n) = m(3)
 DO  j = 1,6
   IF (     rbe) i(n+j) =-ib(j)
   IF (.NOT.rbe) i(n+j) = jb(j)
 END DO
 n = n + 7
 6182 IF (kz /= 49 .OR. m(8+lk) == 0) GO TO 6190
 i(n) = m(4)
 DO  j = 1,6
   IF (     rbe) i(n+j) =-ic(j)
   IF (.NOT.rbe) i(n+j) = jc(j)
 END DO
 n = n + 7
 6190 n = n - 1
 DO  j = 1,7
   i(n+j) = -1
 END DO
 IF (m1(1) == crtr .OR. m1(1) == crba) i(n+2) = 0
 n = n + 7
 GO TO 9
 
!******      349-CRTRPLT      ******************************************
 
!     MAP THIS RIGID ELEMENT INTO CRIGD3 FORM
 
 6300 IF (mf(1)+mf(2)+mf(3)+mf(4) /= 4) GO TO 7
 IF (m(1) <= 0 .OR. m(2) <= 0 .OR. m(3) <= 0 .OR. m(4) <= 0) GO TO 8
 IF (m(2) == m(3) .OR. m(2) == m(4) .OR. m(3) == m(4)) GO TO 8
 IF (m(5) == 0 .AND. m(6) == 0 .AND. m(7) == 0) GO TO 6470
 rbe  = .false.
 IF (m(9) == 0 .AND. m(10) == 0 .AND. m(11) == 0) rbe = .true.
 IF (ifpdco(m(5))) GO TO 6470
 lk = 1
 DO  l = 1,6
   lll = ll(l)
   IF (rbe .AND. lll == 0) m(9) = m(9) + l*lk
   IF (lll == 0) lk = lk*10
   ia(l) = lll
 END DO
 IF (ifpdco(m(6))) GO TO 6470
 lk = 1
 DO  l = 1,6
   lll = ll(l)
   IF (rbe .AND. lll == 0) m(10) = m(10) + l*lk
   IF (lll == 0) lk = lk*10
   ib(l) = lll
 END DO
 IF (ifpdco(m(7))) GO TO 6470
 lk = 1
 DO  l = 1,6
   lll = ll(l)
   IF (rbe .AND. lll == 0) m(11) = m(11) + l*lk
   IF (lll == 0) lk = lk*10
   ic(l) = lll
 END DO
 IF (rbe) GO TO 6365
 IF (ifpdco(m(9))) GO TO 6480
 DO  l = 1,6
   IF (ia(l) ==     0) GO TO 6340
   IF (ia(l) == ll(l)) GO TO 6480
   6340 ja(l) = ll(l)
 END DO
 IF (ifpdco(m(10))) GO TO 6480
 DO  l = 1,6
   IF (ib(l) ==     0) GO TO 6350
   IF (ib(l) == ll(l)) GO TO 6480
   6350 jb(l) = ll(l)
 END DO
 IF (ifpdco(m(11))) GO TO 6480
 DO  l = 1,6
   IF (ic(l) ==     0) GO TO 6360
   IF (ic(l) == ll(l)) GO TO 6480
   6360 jc(l) = ll(l)
 END DO
 6365 IF (.NOT.prt) GO TO 6133
 knt1 = knt
 IF (.NOT.rbe) knt1 = knt - 1
 CALL page2 (5)
 WRITE (nout,6370) scc,knt1,(m(l),l=1,7),(m(l),l=9,11), gcc,m(1),  &
     m(2),m(5),m(3),m(6),m(4),m(7),m(2),m(9),m(3),m(10),m(4),m(11)
 6370 FORMAT (/25X,a19,i7,1H-,5X,'CRTRPLT',i7,6I8, /63X,3I8,  &
     /25X,a19,13X,'CRIGD3',7I8, /67X,'MSET',6I8)
 GO TO 6133
 
!******    351-CRBE2         *******************************************
 
!     MAP THIS RIGID ELEMENT INTO CRIGD2 FORM
 
 6400 kn = 1
 SELECT CASE ( irigid )
   CASE (    1)
     GO TO 6405
   CASE (    2)
     GO TO 6410
 END SELECT
 6405 irigid = irigid + 1
 knt1 = knt
 l6 = 60
 l7 = l6
 l8 = 0
 IF (mf(1)+mf(2)+mf(3) /=   3) GO TO 7
 IF (m(1) <= 0 .OR. m(2) <= 0) GO TO 8
 i(1) = m(1)
 i(2) = m(2)
 q(1) = m(1)
 q(2) = m(2)
 m3   = m(3)
 l8   = l8+2
 q(l7+1) = m(1)
 q(l7+2) = m(2)
 q(l7+3) = m3
 l7  = l7 + 3
 n   = 2
 irg = 4
 IF (ifpdco(m3)) baddat = .true.
 IF (m3 ==  0) baddat = .true.
 GO TO 6420
 6410 n   = 0
 irg = 1
 6420 DO  l = irg,8
   IF (mf(l) == 0) GO TO 6450
   IF (mf(l) /= 1) badfor = .true.
   IF (m(l)  <= 0) baddat = .true.
   IF (l8 >= l6) GO TO 6422
   q(l8+1) = m(l)
   q(l8+2) = m3
   6422 l8 = l8 + 2
   IF (l7 < 92) q(l7+1) = m(l)
   l7 = l7 + 1
   i(n+1) = m(l)
   DO  j = 1,6
     i(n+1+j) = ll(j)
   END DO
   n = n + 7
 END DO
 IF (m1(1) == 0 .AND. m1(2) == 0) GO TO 9
 6440 irigid = 1
 DO  j = 1,7
   i(n+j) = -1
 END DO
 IF (m1(1) == crbe .AND. m1(2) == frigid) i(n+2) = 0
 n  = n + 7
 kn = 0
 IF (.NOT.prt) GO TO 9
 l3 = l7
 l5 = l8
 IF (l3 > 92) l3 = 92
 IF (l5 > l6) l5 = l6
 j = (l5+2)/8 + (l3-l6+2)/8 + 2
 CALL page2 (j)
 l6 = l6 + 1
 WRITE (nout,6447) scc,knt1,(q(j),j=l6,l3)
 WRITE (nout,6448) gcc,(q(j),j=1,l5)
 6447 FORMAT (/25X,a19,i7,1H-,5X,'CRBE2 ',8I8, /,(63X,8I8))
 6448 FORMAT ( 25X,a19,13X,'CRIGD2',8I8, /,(63X,8I8))
 IF (l8 > l6 .OR. l7 > 102) WRITE (nout,6449)
 6449 FORMAT (57X,'*** ABOVE PRINTOUT MAY BE IMCOMPLETE.  DATA IS OK')
 GO TO 9
 6450 l1 = l
 IF (l1 > 8) GO TO 6460
 DO  l = l1,8
   IF (m(l)  /= 0) baddat = .true.
   IF (mf(l) /= 0) badfor = .true.
 END DO
 6460 IF (m1(1) == 0 .AND. m1(2) == 0) baddat = .true.
 GO TO 6440
 
 6470 WRITE  (nout,6475) ufm,ind,m(1),knt1
 6475 FORMAT (a23,', ILLEGAL ',a2,'DEPENDENT D.O.F.',  &
     ' FOR RIGID ELEMENT',i9,' SORTED COUNT',i8)
 GO TO 8
 6480 WRITE (nout,6475) ufm,blnk,m(1),knt1
 GO TO 8
 
!******    352-CRBE3         *******************************************
 
!     CARD 3, OR CARDS 2 AND 3, CAN BE OMITTED IF THE CARD(S) CONTAINS
!     ALL BLANKS.
!     CARD 5, OR CARDS 4 AND 5, CAN BE OMITTED IF THE CARD(S) CONTIANS
!     ALL BLANKS, OR DEFAULT FOR THE 'UM' OPTION IS USED
 
!     ACTUALLY THIS CRBE3 INPUT CARD IS NOT WHAT SHOWN IN THE USER'S
!     MANUAL. THE LIST OF G(I,J) CAN BE AS LONG AS NEEDED. THEREFORE
!     CARDS 2 AND 3 CAN BE EXPANDED BEYOND THE 3 GRID POINTS AS SHOWN.
!     THE 4TH AND 5TH CARDS CAN BE EXPANDED TOO. THE WI AND CI FIELDS
!     NEED NOT BE IN THE FIELDS AS SHOWN IN THE EXAMPLE OF THE MANUAL
 
!     CHANGES DONE IN 92 VERSION WERE REMOVED AND REPLACED BY 91 CODE
!     SEE 93 CODE FOR THESE CHANGES
 
!     IM HERE IS CARD NUMBER COUNT
 
 6500 CONTINUE
 IF (km /= 0) GO TO 6510
 km = 1
 im = 1
 IF (mf(1)+mf(3)+mf(4) /= 3) badfor = .true.
 IF (m(1) <= 0 .OR. m(3) <= 0 .OR. m(4) <= 0) baddat = .true.
 IF (ifpdco(m(4))) baddat = .true.
 IF (mf(5) /= 2) baddat = .true.
 i(1) = m(1)
 i(2) = m(3)
 i(3) = m(4)
 
! ... NOTE - COMPONENTS IN LL NOT SENT OUT IN CRBE3
 
 n  = 3
 l1 = 5
 GO TO 6520
 
 
 6510 IF (mf(1) == 3) GO TO 6560
 IF (im    == 0) GO TO 6565
 l1 = 1
 6520 DO  l = l1,8
   IF (mf(l) /= 2) GO TO 6530
   IF (l1 == 5) GO TO 6525
   n   = n + 1
   i(n)=-1
   6525 im  = 1
!WKBI 11/93 SPR93018
   l1  = 1
   n   = n + 1
   i(n)= m(l)
   CYCLE
   6530 IF (mf(l) == 0) CYCLE
   IF (mf(l) /= 1 .OR. m(l) <= 0) baddat =.true.
   IF (im  ==  -1) GO TO 6535
   IF (ifpdco(m(l))) baddat =.true.
   6535 im =-1
   n  = n + 1
   i(n) = m(l)
 END DO
 IF (m1(1) /= 0) GO TO 6550
 kn = 1
 GO TO 9
 6550 n  = n + 1
 i(n) = -1
 6555 kn = 0
 km = 0
 n  = n + 1
 i(n) = -3
 GO TO 9
 6560 IF (m(1) /= ium) baddat =.true.
 i(n+1) = -1
 i(n+2) = -2
 n  = n + 2
 im = 0
 l1 = 3
 GO TO 6570
 6565 l1 = 2
 6570 DO  l = 2,6,2
   IF (mf(l) == 0) GO TO 6575
   IF (mf(l  ) /= 1 .OR. m(l1  ) <= 0) baddat =.true.
   IF (mf(l+1) /= 1 .OR. m(l1+1) <= 0) baddat =.true.
   IF (ifpdco(m(l1+1))) baddat =.true.
   i(n+1) = m(l1  )
   i(n+2) = m(l1+1)
   n  = n  + 2
   6575 l1 = l1 + 2
 END DO
 IF (m1(1) /= 0) GO TO 6555
 GO TO 9
 
!******    353-CRSPLINE      *******************************************
 
 6600 CONTINUE
 IF (km /= 0) GO TO 6610
 km = 1
 im = -1
 IF (mf(1) /= 1 .OR. m(1) <= 0) GO TO 6680
 IF (mf(2) == 0 ) rm(2) = .1
 IF (rm(2) <= 0.) GO TO 6680
 IF (mf(3) /= 1 .OR. m(3) <= 0) GO TO 6680
 i(1) = m(1)
 i(2) = m(2)
 i(3) = m(3)
 n  = 3
 l1 = 4
 GO TO 6620
 6610 l1 = 1
 IF (im == -9) GO TO 6680
 6620 DO  l = l1,8
   IF (mf(l) /= 0 .AND. mf(l) /= 1) GO TO 6680
   IF (im == -1 .AND. m(l) < 0) GO TO 6680
   IF (im == -1 .AND. m(l) == 0) GO TO 6650
   IF (im  ==  -1) GO TO 6630
   IF (ifpdco(m(l))) GO TO 6680
   
! ... NOTE - COMPONENTS IN LL NOT SENT OUT IN CRSPLINE
   
   6630 im = -im
   n  = n + 1
   i(n) = m(l)
 END DO
 IF (m1(1) == 0 .AND. m1(2) == 0) GO TO 6670
 n   = n + 1
 i(n)= 0
 6650 im  = -9
 n   = n + 1
 i(n)= -1
 IF (l == 8) GO TO 6670
 l1 = l
 DO  l = l1,8
   IF (mf(l) /= 0) GO TO 6680
 END DO
 
 6670 kn = 1
 IF (m1(1) == 0 .AND. m1(2) == 0) GO TO 9
 kn = 0
 km = 0
 n  = n + 1
 i(n) = -1
 GO TO 9
 6680 baddat = .true.
 GO TO 6670
 
!******        285-CTRIAAX       ***************************************
 
 2111 IF (m(1) <= 0 .OR. m(2) <= 0) GO TO 8
 IF (m(3) <= 0 .OR. m(4) <= 0) GO TO 8
 IF (m(3) == m(4)) GO TO 8
 IF (m(3) == m(5)) GO TO 8
 ih = nn
 ASSIGN    8 TO r1
 ASSIGN 2172 TO r
 GO TO 21
 2172 n = 6
 GO TO 3
 
!******       286-PTRIAX, 288-PTRAPAX   *******************************
 
 2030 IF (m(1) <= 0) GO TO 8
 ih = nn
 ASSIGN 8 TO r1
 ASSIGN 2031 TO r
 GO TO 21
 2031 n = 17
 GO TO 3
 
!*******       287-CTRAPAX             ********************************
 
 2040 IF (m(1) <= 0 .OR. m(2) <= 0) GO TO 8
 IF (m(3) == m(4)) GO TO 8
 IF (m(3) == m(5)) GO TO 8
 ih = nn
 ASSIGN 8 TO r1
 ASSIGN 2041 TO r
 GO TO 21
 2041 n = 7
 GO TO 3
 
!******         28-GENEL         **************************************
 
 2800 SELECT CASE ( l0 )
   CASE (    1)
     GO TO 2802
   CASE (    2)
     GO TO 2810
   CASE (    3)
     GO TO 2830
   CASE (    4)
     GO TO 2810
   CASE (    5)
     GO TO 2836
   CASE (    6)
     GO TO 2844
   CASE (    7)
     GO TO 2858
   CASE (    8)
     GO TO 2844
 END SELECT
 2802 l0 = l0 + 1
 kzflag = 0
 l8 = 0
 noud = .true.
 nos  = .true.
 IF (mf(1) /= 1 .OR. mf(2) /= 0) badfor = .true.
 IF (m(1) <= 0) baddat = .true.
 id = m(1)
 i(1) = id
 n  = 1
 l3 = 3
 GO TO 2812
 2810 l3 = 1
 2811 n = 0
 2812 DO  l = l3,8
   IF (mf(l) /= 0 .AND. mf(l) /= 1) badfor = .true.
 END DO
 l5 = 1
 DO  l = l3,7,2
   IF (m(l) == 0) GO TO 2824
   l5 = l  + 2
   l8 = l8 + 1
   n  = n  + 2
   i(n-1) = m(l  )
   i(n  ) = m(l+1)
   IF (m(l) <= 0) GO TO 2816
   IF (m(l+1) >= 0 .AND. m(l+1) <= 6) CYCLE
   2816 baddat = .true.
 END DO
 IF (m1f(2) /= 3) GO TO 2864
 2820 n = n + 2
 i(n-1) = -1
 i(n) = l8
 l0 = l0 + 1
 IF (l0 == 5) GO TO 2822
 l6 = l8
 GO TO 2864
 2822 l7 = l8
 GO TO 2864
 2824 DO  l = l5,7,2
   IF (m(l) /= 0 .OR. m(l+1) /= 0) baddat = .true.
 END DO
 IF (l5 <= 1) baddat = .true.
 IF (m1f(2) == 3) GO TO 2820
 baddat = .true.
 GO TO 2864
 2830 l0 = l0 + 1
 lb = 0
 IF (mf(1) /= 3 .OR. (m(1) /= lz .AND. m(1) /= kk)) GO TO 2831
 l0 = l0 + 1
 lb = 2
 i(1) =-1
 i(2) = 0
 GO TO 2838
 2831 l8 = 0
 IF (mf(1) /= 3 .OR. mf(2) /= 0) badfor = .true.
 IF (m(1) /= lud) baddat = .true.
 l3 = 3
 noud = .false.
 DO  l = 2,8
   m(l) = m(l+1)
 END DO
 GO TO 2811
 2836 IF (m(1) /= lz .AND. m(1) /= kk) baddat = .true.
 2838 l9 = (l6*(l6+1))/2
 lb = lb + 1
 IF (m(1) == lz) kzflag = 1
 IF (m(1) == kk) kzflag = 2
 i(lb) = kzflag
 2840 l0 = l0 + 1
 l8 = 0
 IF (mf(1) /= 3) badfor = .true.
 l3 = 2
 DO  l = 2,8
   m(l) = m(l+1)
 END DO
 GO TO 2846
 2844 l3 = 1
 lb = 0
 2846 DO  l = l3,8
   IF (mf(l) /= 2 .AND. mf(l) /= 0) badfor = .true.
 END DO
 n  = lb
 l5 = l9 - l8 + l3 - 1
 IF (l5 <= 8) GO TO 2850
 l5 = 8
 2850 DO  l = l3,l5
   n = n + 1
   i(n) = m(l)
 END DO
 l5 = l9 - l8 + l3
 l8 = l8 + n - lb
 IF (l9 > l8) GO TO 2864
 IF (l9 == l8) GO TO 2855
 DO  l = l5,8
   IF (m(l) /= 0) baddat = .true.
 END DO
 2855 IF (l0 == 8) GO TO 2856
 l0 = l0 + 1
 GO TO 2864
 2856 l0 = 1
 GO TO 2864
 2858 IF (m(1) /= ls) baddat = .true.
 l9 = l6*l7
 lb = 1
 i(1) = l7
 nos  = .false.
 GO TO 2840
 2864 DO  l = 1,8
   m(l) = 0
 END DO
 kn = 1
 IF (m1(1) == 0 .AND. m1(2) == 0) GO TO 9
 kn = 0
 IF (id <= la1) GO TO 2868
 la1 = id
 IF (.NOT.noud .AND. l7 /= 6 .AND. nos .AND. kzflag == 1) GO TO 2868
 IF (l7 == 0 .AND. .NOT. nos) GO TO 2868
 l7 = 0
 IF (l0 == 1 .AND. .NOT.noud) GO TO 9
 n = n + 1
 i(n) = 0
 l0 = 1
 GO TO 9
 2868 baddat = .true.
 l0 = 1
 l7 = 0
 la1= id
 GO TO 9
 
!******        345-STREAML1      **************************************
 
 2920 IF (km == 1) GO TO 2921
 km = 1
 IF (mf(1) /= 1) badfor = .true.
 IF (m(1) <= 0) baddat = .true.
 IF (m(1) <= 0) baddat = .true.
 i(1) = m(1)
 n = 1
 IF (mf(3) == 3 .AND. m(3) == thru) GO TO 2928
 l1 = 2
 GO TO 2922
 2921 l1 = 1
 2922 DO  l = l1,8
   IF (mf(l) /= 0 .AND. mf(l) /= 1) badfor = .true.
 END DO
 DO  l = l1,8
   IF (m(l) < 0) THEN
     GO TO  2925
   ELSE IF (m(l) == 0) THEN
     GO TO  2926
   END IF
   2924 n = n + 1
   i(n) = m(l)
   CYCLE
   2925 baddat = .true.
   2926 CONTINUE
 END DO
 IF (n < l1) baddat = .true.
 kn = 1
 IF (m1(1) == 0 .AND. m1(2) == 0) GO TO 9
 2927 km = 0
 n = n + 1
 i(n) = -1
 kn = 0
 GO TO 9
 2928 IF (m1(1) /= 0 .OR. m1(2) /= 0) GO TO 2929
 kn = 1
 badfor = .true.
 GO TO 9
 2929 IF (mf(2) /= 1 .OR. mf(4) /= 1) badfor = .true.
 IF (m(2) <= 0 .OR. m(5) <= 0 .OR. (m(2) > m(5))) baddat = .true.
 IF (badfor .OR. baddat) GO TO 2927
 CALL WRITE (204,i,n,0)
 l1 = m(2)
 l2 = m(5)
 DO  l = l1,l2
   CALL WRITE (204,l,1,0)
 END DO
 n = 0
 GO TO 2927
 
!******        346-STREAML2      **************************************
 
!     THEORY DEPENDENT RESTRICTION -  (3.GE. NSTNS .LE.10)
 
 3010 IF (m(1) <= 0) GO TO 8
 IF (m(2) < 3 .OR. m(2) > 10) GO TO 8
 IF (rm(4) <= 0.0) GO TO 8
 DO  l = 6,9
   IF (rm(l) <= 0.0) GO TO 8
 END DO
 IF (rm( 3) <= -90.0 .OR. rm( 3) >= 90.0) GO TO 8
 IF (rm(10) <= -90.0 .OR. rm(10) >= 90.0) GO TO 8
 n = 10
 GO TO 3
 
!******         82-PARAM         ***********************************
 
 3960 IF (mf(1) /= 3 .OR. mf(2) <= 0 .OR. mf(3) /= 0 .AND.  &
     mf(3) /= mf(2)) GO TO 3968
 IF (mf(3) /= 0 .AND. mf(3) /= 2 .AND. mf(3) /= 4) GO TO 3968
 DO  l = 4,8
   IF (mf(l) /= 0) GO TO 3968
 END DO
 IF (nparam+7 <= nopen) GO TO 3964
 CALL page2 (2)
 WRITE  (nout,3962) sfm
 3962 FORMAT (a25,' 330, NO ROOM IN CORE FOR PARAM CARDS.')
 3963 abort = .true.
 GO TO 2
 3964 ip = 2*nbuf + nparam
 ibuff(ip+1) = m(1)
 ibuff(ip+2) = m(2)
 ibuff(ip+3) = mf(2)
 ibuff(ip+4) = m(3)
 nparam = nparam + 4
 IF (mf(2) <= 2 .AND. mf(3) == 0) GO TO 2
 ibuff(ip+5) = m(4)
 nparam = nparam + 1
 IF (mf(2) <= 4 .AND. mf(3) == 0) GO TO 2
 IF (mf(3) == 4) GO TO 3965
 ibuff(ip+3) = 5
 GO TO 2
 3965 ibuff(ip+3) = 6
 ibuff(ip+6) = m(5)
 ibuff(ip+7) = m(6)
 nparam = nparam + 2
 GO TO 2
 3968 WRITE  (nout,3969) ufm,m(1),m(2),knt
 3969 FORMAT (a23,' 331, IMPROPER PARAM CARD ',2A4,10X,  &
     'SORTED CARD COUNT =',i7)
 CALL page2 (2)
 GO TO 3963
 
!*******    12-SPC1(3980), 92-OMIT1(3981), 216-ASET1(3981)   ***********
!          332-CFLSTR(3980)
 
 3980 iz = 2
 ifile = 210
 IF (k == 332) ifile = 208
 GO TO 3983
 3981 iz = 1
 ifile = 210
 3983 IF (km /= 0) GO TO 3990
 km = 1
 IF (mf(iz) /= 0 .AND. mf(iz) /= 1) badfor = .true.
 IF (k == 332) GO TO 3984
 IF (ifpdco(m(iz))) baddat = .true.
 IF (iz /= 2) GO TO 3986
 3984 IF (mf(1) /= 1) badfor = .true.
 IF (m(1)  <= 0) baddat = .true.
 3986 id = m(1)
 i(1) = m(1)
 IF (iz == 2) i(2) = m(2)
 n  = iz
 l1 = iz + 1
 IF (mf(iz+2) == 3 .AND. m(iz+2) == thru) GO TO 4000
 GO TO 3991
 3990 l1 = 1
 3991 DO  l = l1,8
   IF (mf(l) /= 0 .AND. mf(l) /= 1) badfor = .true.
 END DO
 DO  l = l1,8
   IF (mf(l) == 1) GO TO 3994
 END DO
 baddat = .true.
 3994 DO  l = l1,8
   IF (m(l) < 0) THEN
     GO TO  3996
   ELSE IF (m(l) == 0) THEN
     GO TO  3998
   END IF
   3995 n = n + 1
   i(n) = m(l)
   CYCLE
   3996 baddat = .true.
   3998 CONTINUE
 END DO
 kn = 1
 IF (m1(1) == 0 .AND. m1(2) == 0) GO TO 9
 3999 km = 0
 n  = n + 1
 i(n) =-1
 kn = 0
 GO TO 9
 4000 IF (m1(1) /= 0 .OR. m1(2) /= 0) GO TO 4001
 kn = 1
 badfor = .true.
 GO TO 9
 4001 IF (mf(iz+1) /= 1 .OR. mf(iz+3) /= 1) badfor = .true.
 IF (m(iz+1) <= 0 .OR. m(iz+4) <= m(iz+1)) baddat = .true.
 DO  l = iz,4
   IF (mf(l+4) /= 0) badfor = .true.
 END DO
 IF (badfor .OR. baddat) GO TO 3999
 CALL WRITE (ifile,m,iz,0)
 l1 = m(iz+1)
 l2 = m(iz+4)
 l = l1
 4010 CALL WRITE (ifile,l,1,0)
 l = l + 1
 IF (l <= l2) GO TO 4010
 n = 0
 GO TO 3999
 
!******      13-SPCADD, 83-MPCADD    **********************************
 
 4020 IF (km == 1) GO TO 4990
 km = 1
 IF (m(1) <= 0) baddat = .true.
 id = m(1)
 i(1) = id
 IF (m(2) <= 0 .OR. m(3) < 0) baddat = .true.
 IF (m(3) == 0) CALL page2 (2)
 IF (m(3) == 0) WRITE (nout,4024) uwm
 4024 FORMAT (a25,' 4124, THE SPCADD OR MPCADD UNION CONSISTS OF A ',  &
     'SINGLE SET.')
 n = 1
 GO TO 4992
 
!******        84-LOAD, 123-DLOAD      *******************************
 
 4060 IF (km == 1) GO TO 4068
 km = 1
 IF (mf(1) /= 0 .AND. mf(1) /= 1 .OR. mf(2) /= 0 .AND. mf(2) /= 2)  &
     badfor = .true.
 IF (m(1) <= 0) baddat = .true.
 id = m(1)
 i(1) = id
 i(2) = m(2)
 IF (m(4) <= 0) baddat = .true.
 n = 2
 GO TO 4070
 4068 n = 0
 4070 l8 = n + 1
 DO  l = l8,7,2
   IF (mf(l  ) /= 0 .AND. mf(l) /= 2 .OR. mf(l+1) /= 0 .AND.  &
       mf(l+1) /= 1) badfor = .true.
 END DO
 4076 n = n + 2
 IF (m(n) < 0) THEN
   GO TO  4078
 ELSE IF (m(n) == 0) THEN
   GO TO  4080
 ELSE
   GO TO  4084
 END IF
 4078 baddat = .true.
 4080 n  = n - 2
 l7 = 1
 l8 = n + 1
 DO  l = l8,8
   IF (mf(l) /= 0) baddat = .true.
 END DO
 IF (n <= 0) baddat = .true.
 GO TO 4086
 4084 i(n-1) = m(n-1)
 i(n  ) = m(n  )
 IF (n < 8) GO TO 4076
 4086 kn = 1
 IF (m1(1) == 0 .AND. m1(2) == 0) GO TO 4088
 km = 0
 n  = n + 2
 i(n-1) =-1
 i(n  ) =-1
 kn = 0
 GO TO 4090
 4088 IF (l7 /= 1) GO TO 9
 baddat = .true.
 4090 l7 = 0
 GO TO 9
 
!     ******************************************************************
 
 21 IF (.NOT.iax) GO TO 22
 IF (ih > nn .OR. ih < 0) GO TO 25
 GO TO 24
 22 IF (lh) WRITE (nout,23) ufm
 23 FORMAT (a23,' 332, AXIC CARD REQUIRED.')
 IF (lh) CALL page2 (2)
 lh = .false.
 abort = .true.
 24 GO TO  r, (1489,1507,1521,1451,1461,1471,1497,2031,2041,2172,3,2)
 25 GO TO r1, (8)
 
!***** TEMPORARY UNFIX FOR SPCADD AND MPCADD ***************************
 
 4990 n = 0
 4992 DO  l = 1,8
   IF (mf(l) /= 0 .AND. mf(l) /= 1) badfor = .true.
 END DO
 4995 n = n + 1
 IF (m(n) < 0) THEN
   GO TO  4996
 ELSE IF (m(n) == 0) THEN
   GO TO  4998
 ELSE
   GO TO  5002
 END IF
 4996 baddat = .true.
 4998 n  = n - 1
 l7 = 1
 l8 = n + 1
 DO  l = l8,8
   IF (mf(l) /= 0) baddat = .true.
 END DO
 IF (n <= 0) baddat = .true.
 GO TO 5004
 5002 i(n) = m(n)
 IF (n < 8) GO TO 4995
 5004 kn = 1
 IF (m1(1) == 0 .AND. m1(2) == 0) GO TO 5006
 km = 0
 n  = n + 1
 i(n) =-1
 kn = 0
 GO TO 5008
 5006 IF (l7 /= 1) GO TO 9
 baddat = .true.
 5008 l7 = 0
 GO TO 9
 
!******     329-PROLATE     ********************************************
 
 5100 IF (km /= 0) GO TO 5115
 IF (prol) baddat = .true.
 prol = .true.
 km = 1
 IF (mf(1) /= 2 .OR. mf(2) /= 2) badfor = .true.
 IF (rm(1) <= rm(2)) baddat = .true.
 DO  l = 3,6
   IF (mf(l) /= 1) badfor = .true.
   IF (m(l)  < 0) baddat = .true.
 END DO
 IF (m(3) <  2) baddat = .true.
 IF (m(4) <  2) baddat = .true.
 IF (m(5) > 30) baddat = .true.
 IF (m(6) > m(5)) m(6) = m(5)
 id = m(1)
 nsegs = m(3)
 msegs = m(4)
 itot1 = (nsegs-1)*(msegs+1) + 2
 itot2 = (nsegs-1)*msegs + 2
 DO  l = 1,6
   i(l) = m(l)
 END DO
 n  = 6
 l1 = 7
 items = 0
 GO TO 5120
 5115 l1 = 1
 5120 DO  l = l1,8
   IF (mf(l) /= 1 .AND. mf(l) /= 3) badfor = .true.
   IF (mf(l) == 3) GO TO 5140
   items = items + 1
   IF (m(l) <= 0) baddat = .true.
   n = n + 1
   i(n) = m(l)
 END DO
 kn = 1
 5135 IF (m1(1) == 0 .AND. m1(2) == 0) GO TO 9
 baddat = .true.
 GO TO 5150
 5140 IF (m(l) /= endt) GO TO 5145
 IF (m1(1) == 0 .AND. m1(2) == 0) baddat = .true.
 GO TO 5150
 5145 baddat = .true.
 GO TO 5135
 5150 km = 0
 kn = 0
 IF (k == 330) GO TO 9
 IF (items /= itot1 .AND. items /= itot2) baddat = .true.
 GO TO 9
 
!******      330-PERMBDY       *****************************************
 
 5200 IF (km /= 0) GO TO 5210
 IF (perm) baddat = .true.
 perm = .true.
 km = 1
 5210 DO  l = 1,8
   IF (mf(l) /= 1 .AND. mf(l) /= 3) badfor = .true.
   IF (mf(l) == 3) GO TO 5140
   IF (m(l)  <= 0) baddat = .true.
   n = n + 1
   i(n) = m(l)
 END DO
 kn = 1
 GO TO 5135
 
 2 RETURN
END SUBROUTINE ifs3p
