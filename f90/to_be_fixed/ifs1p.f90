SUBROUTINE ifs1p (*,*,*)
     
 
 , INTENT(IN OUT)                         :: *
 , INTENT(IN OUT)                         :: *
 , INTENT(IN OUT)                         :: *
 LOGICAL :: abort,baddat,badfor,iax,lharm,slot,ifpdco
 INTEGER :: m(100),klotdf(5),b1,bardf2,bardf5,bardf6,bardf7,  &
     bardf8,hbdynm(2,7),hbdyix(7),thru,blk,bcdc,bcdr, bcds,e(40)
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim,sfm
 COMMON /system/ ksystm(80)
 COMMON /BLANK / e
 COMMON /ifpdta/ id,n,k,kx,ky,i(100),rm(100),mf(100),m1(100),  &
     m1f(100),kn,baddat,badfor,nopen,nparam,iax,nax,  &
     iaxf,naxf,lharm,knt,slotdf(5),gc(7),ll(6)
 COMMON /cifs1p/ b1,bardf2,bardf5,bardf6,bardf7,bardf8,km,slot, idrdl
 EQUIVALENCE     (ksystm(2),nout),(ksystm(3),abort),(m(1),rm(1)),  &
     (slotdf(1) ,klotdf(1))
 DATA   hbdynm / 4HPOIN , 4HT , 4HLINE , 4H  &
     , 4HREV  , 4H , 4HAREA , 4H3  &
     , 4HAREA , 4H4 , 4HELCY , 4HL  &
     , 4HFTUB , 4HE   /
 DATA   hbdyix / 1,2,2,3,4,2,2  /
 DATA   thru   / 4HTHRU         /
 DATA   blk    , bcdc,bcdr,bcds /  1H ,1HC,1HR,1HS/
 DATA   it1,it2, it3 / 2HT1, 2HT2, 2HT3           /
 
 IF (k > 100) GO TO 81
 GO TO (     5,   5,   5,  40, 500, 600, 700, 800, 900,1000,  &
     1111,   5,   5,1400,1400,1600,   5,1800,1800,2000,  &
     2000,2200,2200,2400,2500,2600,2500,   5,2900,2920,  &
     318,   5,2980,3011,3020,3020,3012,2980,3013,3020,  &
     3014,3015,3016,3210,3220,3255,3260,3281,3282,3283,  &
     5,3360,3360,3360,3360,3360,3460,3460,3460,3460,  &
     3540,3540,3580,3600,3620,3623,3674,3697,3620,3623,  &
     3675,3698,3620,3800,3676,3699,3860,3880,   5,   5,  &
     2500,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
     5,   5,   5,   5,   5,   5,   5,   5,   5,   5), k
 81 IF (kx > 100) GO TO 82
 GO TO (     5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
     5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
     3017,   5,   5,   5,1250,   5,1270,1280,1290,1290,  &
     5,   5,   5,   5,  40,1360,1370,   5,   5,   5,  &
     5,1420,   5,   5,   5,   5,   5,   5,   5,   5,  &
     5,   5,   5,   5,   5,   5,   5,1580,   5,   5,  &
     5,   5,   5,   5,   5,1660,   5,   5,   5,   5,  &
     5,   5,   5,   5,   5,   5,   5,   5, 100, 200,  &
     300,   5,   5,   5,   5,   5,   5,   5,   5,1900,  &
     5,   5,   5,   5,   5,   5,   5,   5,   5,   5), kx
 82 IF (ky > 100)  GO TO 83
 GO TO (     5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
     5,   5,   5,   5,1400,   5,   5,   5,   5,   5,  &
     5,   5,4100,4200,4300,4400,4500,4600,4700,4800,  &
     4900,5000,5050,3900,4000,5100,3950,4050,   5,5150,  &
     5200,   5,5250,   5,   5,   5,   5,   5,3460,3018,  &
     5,   5,   5,   5,   5,1600,5240,5245,3460,3019,  &
     5,   5,   5,   5,   5,   5,   5,5300,   5,   5,  &
     5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
     5,   5,   5,   5,   5,   5,   5,   5,5175,   5,  &
     6101,6201,6301,6401,   5,   5,   5,   5,7501,7601), ky
 83 kz = ky - 100
 IF (kz > 59) GO TO 5
 GO TO (     5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
     5,   5,   5,   5,4060,4070,4080,4090,4060,4080,  &
     5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
     1111,   5,   5,   5,   5,   5,3900,   5,   5,1260,  &
     1270,1230,1235,1240,   5,   5,   5,   5,   5,   5,  &
     5,   5,   5,   5,   5,3280,3360,3460,7700     ), kz
 5 CALL page2 (2)
 WRITE  (nout,6) sfm
 6 FORMAT (a25,' 322, ILLEGAL ENTRY TO IFS1P.')
 abort  = .true.
 RETURN 1
 7 badfor = .true.
 RETURN 1
 8 baddat = .true.
 RETURN 1
 3 DO  l = 1,n
   i(l) = m(l)
 END DO
 2 RETURN
 9 RETURN 3
 
!*****         4-SEQGP,135-SEQEP    ************************************
 
 40 DO  l = 1,7,2
   IF (m(l) == 0 .AND. m(l+1) == 0) CYCLE
   IF (m(l) <= 0 .OR.  m(l+1) <= 0) GO TO 8
   n = n + 2
   i(n-1) = m(l  )
   i(n  ) = m(l+1)
   IF (n <= 2) CYCLE
   DO  l1 = 4,n,2
     IF (i(n-1) == i(l1-3) .OR. i(n) == i(l1-2)) GO TO 8
   END DO
 END DO
 IF (n > 0) THEN
   GO TO     2
 ELSE
   GO TO     8
 END IF
 
!*****         179-BAROR       *****************************************
 
 100 IF (b1 == 0) GO TO 8
 b1 = 0
 IF (m(2) == 0 .AND. m(5) == 0 .AND. m(6) == 0 .AND. m(7) == 0  &
     .AND. m(8) == 0) GO TO 8
 IF (m(2) < 0 .OR. m(8) < 0 .OR. m(8) > 2) GO TO 8
 IF (mf(8) /= 0) GO TO 110
 IF (mf(5) == 1 .AND. mf(6) /= 0 .AND. mf(7) /= 0) GO TO 8
 IF (mf(5) == 1 .AND. mf(6) == 0 .AND. mf(7) == 0) m(8) = 2
 IF (mf(5) == 2 .OR.  mf(6) == 2 .OR.  mf(7) == 2) m(8) = 1
 110 bardf2 = m(2)
 bardf5 = m(5)
 bardf6 = m(6)
 bardf7 = m(7)
 bardf8 = m(8)
 RETURN 2
 
!*****         180-CBAR        *****************************************
 
 200 IF (mf(2)  /= 0) GO TO 201
 IF (bardf2 == 0) GO TO 203
 m(2) = bardf2
 GO TO 201
 203 m(2) = m(1)
 201 CONTINUE
 IF (mf(5) == 0) m(5) = bardf5
 IF (mf(8) == 0) m(8) = bardf8
 IF (mf(5) >= 3 .OR.  mf(6) >= 3 .OR. mf(7) >= 3) GO TO 8
 IF (m(8) == 0 .AND. (mf(5) == 2 .OR. mf(6) == 2 .OR. mf(7) == 2)) m(8) = 1
 IF (m(8) == 0 .AND. mf(5) == 1 .AND. mf(6)+mf(7) == 0) m(8) = 2
 IF (m(8) <= 0 .OR.  m(8) > 2) GO TO 8
 IF (m(8)  == 2) GO TO 205
 IF (mf(6) == 0) m(6) = bardf6
 IF (mf(7) == 0) m(7) = bardf7
 205 CONTINUE
 IF (m(1) <= 0 .OR. m(2) <= 0 .OR. m(3) <= 0 .OR. m(4) <= 0) GO TO 8
 IF (m(8) == 1 .AND. (mf(5) /= 2 .AND. mf(5) /= 0 .OR.  &
     m(5) == 0 .AND. m(6) == 0 .AND. m(7) == 0)) GO TO 8
 IF ((m(8) == 2 .OR. m(8) == 3) .AND. (mf(5) /= 1 .AND.mf(5) /= 0  &
     .OR. m(5) <= 0 .OR. m(6) /= 0 .OR. m(7) /= 0)) GO TO 8
 IF (ifpdco(m(9))) GO TO 8
 IF (m(9) > 65432) GO TO 8
 IF (ifpdco(m(10))) GO TO 8
 IF (m(10) > 65432) GO TO 8
 IF (m(3) == m(4) .OR. m(3) == m(5) .AND. m(8) == 2) GO TO 8
 IF (m(8) == 2 .AND. m(4) == m(5)) GO TO 8
 n = 16
 GO TO 3
 
!*****         181-PBAR        *****************************************
 
 300 n = 19
 IF (rm(4) < 0. .OR. rm(5) < 0. .OR. rm(4)*rm(5) < rm(19)**2) GO TO 8
 GO TO 2903
 
!*****         31-PVISC        *****************************************
 
 310 DO  l = 1,5,4
   IF (m(l) == 0 .AND. m(l+1) == 0 .AND. m(l+2) == 0) CYCLE
   IF (m(l) <= 0) GO TO 8
   n = n + 3
   i(n-2) = m(l  )
   i(n-1) = m(l+1)
   i(n  ) = m(l+2)
   IF (e(kl) < 0) CYCLE
   IF (m(l) > e(kl)) GO TO 314
   e(kl) = -m(l)
   CYCLE
   314 e(kl) = m(l)
 END DO
 IF (n > 0) THEN
   GO TO     2
 ELSE
   GO TO     8
 END IF
 318 kl = 33
 GO TO 310
 
!*****         5-CORD1R        *****************************************
 
 500 l50 = 1
 510 DO  l = 1,5,4
   IF (m(l) == 0 .AND. m(l+1) == 0 .AND. m(l+2) == 0 .AND. m(l+3) == 0) CYCLE
   IF (m(l) <= 0 .OR. m(l+1) <= 0 .OR. m(l+2) <= 0 .OR. m(l+3) <= 0) GO TO 8
   IF (m(l+1) == m(l+2) .OR. m(l+1) == m(l+3) .OR. m(l+3) == m(l+2)) GO TO 8
   n = n + 6
   IF (n > 6 .AND. m(l) == m(l-4)) GO TO 8
   i(n-5) = m(l  )
   i(n-4) = l50
   i(n-3) = 1
   i(n-2) = m(l+1)
   i(n-1) = m(l+2)
   i(n  ) = m(l+3)
 END DO
 IF (n > 0) THEN
   GO TO     2
 ELSE
   GO TO     8
 END IF
 
!*****         6-CORD1C        *****************************************
 
 600 l50 = 2
 GO TO 510
 
!*****         7-CORD1S        *****************************************
 
 700 l50 = 3
 GO TO 510
 
!*****         8-CORD2R        *****************************************
 
 800 i(2) = 1
 810 i(1) = m(1)
 IF (m(1) <= 0 .OR. m(2) < 0) GO TO 8
 IF (m(3) == m(6) .AND. m(4) == m( 7) .AND. m(5) == m( 8)) GO TO 8
 IF (m(3) == m(9) .AND. m(4) == m(10) .AND. m(5) == m(11)) GO TO 8
 IF (m(6) == m(9) .AND. m(7) == m(10) .AND. m(8) == m(11)) GO TO 8
 i(3) = 2
 DO  l = 2,11
   i(l+2) = m(l)
 END DO
 n = 13
 GO TO 2
 
!*****         9-CORD2C        *****************************************
 
 900 i(2) = 2
 GO TO 810
 
!*****         10-CORD2S       *****************************************
 
 1000 i(2) = 3
 GO TO 810
 
!*****   11-PLOTEL,   331-CFFREE   *************************************
 
 1100 DO  l = 1,5,4
   IF (m(l) == 0 .AND. m(l+1) == 0 .AND. m(l+2) == 0) CYCLE
   IF (m(l) <= 0 .OR.  m(l+1) <= 0 .OR.  m(l+2) <= 0) GO TO 8
   IF (m(l+1) == m(l+2)) GO TO 8
   n = n + 3
   i(n-2) = m(l  )
   i(n-1) = m(l+1)
   i(n  ) = m(l+2)
   IF (e(kl) < 0) CYCLE
   IF (m(l) > e(kl)) GO TO 1107
   e(kl) = -m(l)
   CYCLE
   1107 e(kl) = m(l)
 END DO
 IF (n > 0) THEN
   GO TO     2
 ELSE
   GO TO     8
 END IF
 1111 kl = 10
 GO TO 1100
 
!*********       342-CFTUBE    *****************************************
 
 1230 n = 4
 IF (m(1) <= 0  .OR. m(3) <= 0 .OR. m(4) <= 0) GO TO 8
 IF (mf(2) == 0) m(2) = m(1)
 IF (m(2) <= 3 .OR. m(3) == m(4)) GO TO 8
 GO TO 3
 
!*********       343-PFTUBE    ****************************************
 
 1235 n = 5
 IF (m(1)  <= 0) GO TO 8
 IF (rm(2) <= 0. .OR. rm(3) < 0. .OR. rm(4) <= 0.) GO TO 8
 IF (rm(5) == 0.) rm(5) = rm(4)
 IF (rm(5) < 0.) GO TO 8
 GO TO 3
 
!*********       344-NFTUBE    *****************************************
 
 1240 n = 5
 IF (mf(2) /= 1 .OR. m(1) <= 0) GO TO 8
 IF (mf(2) /= 1 .OR. m(2) <= 0) GO TO 8
 IF (mf(3) /= 1 .OR. m(3) <= 0) GO TO 8
 IF (m(2) == m(3)) GO TO 8
 IF (mf(4) /= 0 .AND. mf(4) /= 2) GO TO 8
 IF (mf(5) == 1 .AND.  m(5) < 0) GO TO 8
 IF (mf(5) > 2) GO TO 7
 GO TO 3
 
!***********       125-FREQ1      **************************************
 
 1250 IF (m(1) <= 0 .OR. rm(2) < 0. .OR. rm(3) <= 0. .OR. m(4) <= 0) GO TO 8
 n = 4
 GO TO 3
 
!*****             340-NOLIN5         **********************************
 
 1260 IF (km /= 0) GO TO 1262
 km  = 1
 kn  = 1
 nmo = 8
 IF (mf(1) /= 1 .OR.  m(1) <= 0 ) baddat =.true.
 IF (mf(2) /= 2 .OR. rm(2) <= 0.) baddat =.true.
 IF (mf(3) /= 2 .OR. rm(3) <= 0.) baddat =.true.
 IF (mf(4) /= 2 .OR. rm(4) <= 0.) baddat =.true.
 IF (mf(5) == 1 .AND. m(5) < 0 ) baddat =.true.
 IF (mf(6) == 1 .AND. m(6) < 0 ) baddat =.true.
 IF (mf(7) == 1 .AND. m(7) < 0 ) baddat =.true.
 IF (mf(8) == 1 .AND. m(8) < 0 ) baddat =.true.
 IF (mf(5) == 2 .AND.rm(5) < 0.) baddat =.true.
 IF (mf(6) == 2 .AND.rm(6) < 0.) baddat =.true.
 IF (mf(7) == 2 .AND.rm(7) < 0.) baddat =.true.
 IF (mf(8) == 2 .AND.rm(8) < 0.) baddat =.true.
 n = 8
 DO  l = 1,8
   i(l) = m(l)
 END DO
 GO TO 1265
 1262 n = 8
 nmo = nmo + 8
 DO  l = 1,8
   i(l) = m(l)
   IF (mf(l) == 0) CYCLE
   IF (mf(l) /= 1 .OR. m(l) <= 0) baddat =.true.
 END DO
 1265 IF (m1(1) == 0 .AND. m1(2) == 0) GO TO 9
 km = 0
 kn = 0
 IF (nmo == 16) GO TO 9
 IF (nmo > 16) baddat =.true.
 DO  l = 1,8
   n = n + 1
   i(n) = 0
 END DO
 GO TO 9
 
!*****         127-NOLIN1,341-NOLIN6    ********************************
 
 1270 IF (mf(8) == 0) THEN
   GO TO  1282
 ELSE
   GO TO     8
 END IF
 1272 IF (m(1) <= 0 .OR. m(2) <= 0 .OR. m(3) < 0 .OR. m(5) <= 0 .OR.  &
     m(6) < 0) GO TO 8
 IF (m(3) > 6) GO TO 8
 IF ((m(6) > 6 .AND. m(6) < 10) .OR. m(6) > 16) GO TO 8
 n = 8
 GO TO 3
 
!*****         128-NOLIN2        ***************************************
 
 1280 IF (m(8) < 0 .OR. mf(8) /= 1 .AND. mf(8) /= 0) GO TO 8
 IF ((m(8) > 6 .AND. m(8) < 10) .OR. m(8) > 16) GO TO 8
 1282 IF (mf(7) /= 1 .OR. m(7) <= 0) GO TO 8
 GO TO 1272
 
!*****         129-NOLIN3,130-NOLIN4        ****************************
 
 1290 IF (mf(8) /= 0 .OR. mf(7) /= 2 .AND. mf(7) /= 0) GO TO 8
 GO TO 1272
 
!*****         136-TF         ******************************************
 
 1360 IF (km /= 0) GO TO 1363
 nmo = 5
 id  = m(1)
 IF (id <= 0 .OR. m(2) <= 0 .OR. m(3) < 0) GO TO 1427
 IF (mf(1) /= 1 .OR. mf(2) /= 1 .OR. mf(3) > 1) badfor =.true.
 IF ((mf(4) /= 2 .AND. mf(4) /= 0) .OR. (mf(5) /= 2 .AND.  &
     mf(5) /= 0) .OR. (mf(6) /= 2 .AND. mf(6) /= 0)) badfor =.true.
 n = 6
 1361 DO  l = 1,n
   i(l) = m(l)
 END DO
 GO TO 1428
 1363 IF (m(1) <= 0 .OR. m(2) < 0) GO TO 1427
 IF (mf(1) /= 1 .OR. mf(2) > 1) badfor =.true.
 IF ((mf(3) /= 2 .AND. mf(3) /= 0) .OR. (mf(4) /= 2 .AND.  &
     mf(4) /= 0) .OR. (mf(5) /= 2 .AND. mf(5) /= 0)) badfor =.true.
 n = 5
 GO TO 1361
 
!*****         137-TIC        ******************************************
 
 1370 IF (m(1) <= 0 .OR. m(2) <= 0 .OR. m(3) < 0 .OR. m(3) > 6) GO TO 8
 n = 5
 GO TO 3
 
!*****         14-SUPORT,15-OMIT,215-ASET         **********************
 
 1400 l = 1
 1402 IF (m(l) == 0 .AND. m(l+1) == 0) GO TO 1409
 IF (m(l) <= 0) GO TO 8
 IF (ifpdco(m(l+1))) GO TO 8
 iz = 6
 IF (m(l+1) == 0) iz = 1
 DO  l2 = 1,iz
   IF (iz /= 1 .AND. ll(l2) == 0) CYCLE
   n = n + 2
   i(n-1) = m(l  )
   i(n  ) = ll(l2)
   IF (n <= 2) CYCLE
   DO  l1 = 4,n,2
     IF (i(n-1) == i(l1-3) .AND. i(n) == i(l1-2)) GO TO 8
   END DO
 END DO
 1409 l = l + 2
 IF (l <= 7) GO TO 1402
 IF (n > 0) THEN
   GO TO     2
 ELSE
   GO TO     8
 END IF
 
!*****         142-TSTEP        ****************************************
 
 1420 IF (mf(5) /= 0 .OR. mf(6) /= 0 .OR. mf(7) /= 0 .OR. mf(8) /= 0) GO TO 7
 IF (km /= 0) GO TO 1422
 nmo = 3
 id  = m(1)
 IF (id <= 0 .OR. mf(1) /= 1) GO TO 1427
 n   = 1
 i(n) = m(1)
 GO TO 1425
 1422 IF (mf(1) /= 0) GO TO 1427
 1425 IF (mf(2) /= 1 .OR. mf(4) /= 1 .OR. mf(3) /= 2) GO TO 1427
 IF (m(4) <= 0 .OR. rm(3) <= 0. .OR. m(2) < m(4)) GO TO 1427
 n = n + 3
 i(n-2) = m(2)
 i(n-1) = m(3)
 i(n  ) = m(4)
 GO TO 1428
 1427 baddat =.true.
 1428 IF (m1(1) == 0 .AND. m1(2) == 0) GO TO 1429
 km = 0
 kn = 0
 IF (nmo <= 0) GO TO 9
 DO  l = 1,nmo
   n = n + 1
   i(n) =-1
 END DO
 GO TO 9
 1429 km = 1
 kn = 1
 GO TO 9
 
!*****         158-EIGP        *****************************************
 
 1580 IF (m(1) <= 0) GO TO 8
 DO  l = 2,5,3
   IF (m(l) == 0 .AND. m(l+1) == 0 .AND. m(l+2) == 0) CYCLE
   IF (m(l+2) <= 0) GO TO 8
   n = n + 4
   i(n-3) = m(  1)
   i(n-2) = m(  l)
   i(n-1) = m(l+1)
   i(n  ) = m(l+2)
 END DO
 IF (n > 0) THEN
   GO TO     2
 ELSE
   GO TO     8
 END IF
 
!*****         16-SPC , 256-SPCD ***********************************
 
 1600 IF (m(1) <= 0) GO TO 8
 l = 2
 1601 IF (m(l) == 0 .AND. m(l+1) == 0 .AND. m(l+2) == 0) GO TO 1609
 IF (m(l) <= 0 .OR. m(l+1) < 0) GO TO 8
 IF (ifpdco(m(l+1))) GO TO 8
 n = n + 4
 IF (n > 4 .AND. m(l) == m(l-3) .AND. m(l+1) == m(l-2)) GO TO 8
 i(n-3) = m(1  )
 i(n-2) = m(l  )
 i(n-1) = m(l+1)
 i(n  ) = m(l+2)
 1609 l = l + 3
 IF (l == 5) GO TO 1601
 IF (n > 0) THEN
   GO TO     2
 ELSE
   GO TO     8
 END IF
 
!***********       166-FREQ2      **************************************
 
 1660 IF (rm(2) > 0.0) THEN
   GO TO  1250
 ELSE
   GO TO     8
 END IF
 
!*******       18-FORCE,19-MOMENT   **************************
 
 1800 IF (m(2) > 0) THEN
   GO TO  1900
 ELSE
   GO TO     8
 END IF
 
!***************        190-RFORCE    *****************************
 
 1900 IF (mf(3) /= 0 .AND. mf(3) /= 1) GO TO 8
 IF (m(1) <= 0 .OR. m(2) < 0 .OR. m(3) < 0) GO TO 8
 IF (m(5) /= 0 .OR. m(6) /= 0 .OR. m(7) /= 0) GO TO 1905
 IF (m(4) /= 0) GO TO 8
 rm(5) = 1.0
 1905 n = 7
!WKBDB 2/95 SPR94015
!      IF (K .NE. 190) GO TO 3
!      IF (M(8) .EQ. 0) M(8) = 1
!      IF (M(8) .LT.0 .OR. M(8).GT.2) GO TO 8
!      N = 8
!WKBDE 2/95 SPR94015
 GO TO 3
 
!*****         20-FORCE1,21-MOMENT1   **********************************
 
 2000 IF (m(1) <= 0 .OR. m(2) <= 0 .OR. m(4) <= 0 .OR. m(5) <= 0) GO TO 8
 IF (m(4) == m(5)) GO TO 8
 n = 5
 GO TO 3
 
!*****         22-FORCE2,23-MOMENT2   **********************************
 
 2200 IF (m(1) <= 0 .OR. m(2) <= 0 .OR. m(4) <= 0) GO TO 8
 IF (m(5) <= 0 .OR. m(6) <= 0 .OR. m(7) <= 0) GO TO 8
 IF (m(4) == m(5) .OR. m(6) == m(7) .OR.  m(4) == m(6) .AND.  &
     m(5) == m(7) .OR. m(4) == m(7) .AND. m(5) == m(6)) GO TO 8
 n = 7
 GO TO 3
 
!*****         24-PLOAD        *****************************************
 
 2400 IF (m(1) <= 0 .OR. m(3) <= 0 .OR. m(4) <= 0 .OR. m(5) <= 0) GO TO 8
 IF (m(6) < 0 .OR. m(6) == 0 .AND. mf(6) /= 0) GO TO 8
 DO  l  = 4,6
   DO  l1 = l,6
     IF (m(l-1) == m(l1)) GO TO 8
   END DO
 END DO
 n = 6
 GO TO 3
 
!*****         25-SLOAD,27-TEMP,81-DEFORM    ***************************
 
 2500 IF (m(1) <= 0) GO TO 8
 DO  l = 2,6,2
   IF (m(l) == 0 .AND. m(l+1) == 0) CYCLE
   IF (m(l) <= 0) GO TO 8
   n = n + 3
   i(n-2) = m(1  )
   i(n-1) = m(l  )
   i(n  ) = m(l+1)
   IF (n <= 3) CYCLE
   DO  l1 = 6,n,3
     IF (i(n-1) == i(l1-4)) GO TO 8
   END DO
 END DO
 IF (n > 0) THEN
   GO TO     2
 ELSE
   GO TO     8
 END IF
 
!*****         26-GRAV         *****************************************
 
 2600 IF (m(1) <= 0 .OR. m(2) < 0) GO TO 8
 IF (m(4) /= 0 .OR. m(5) /= 0 .OR. m(6) /= 0) GO TO 2605
 IF (m(3) /= 0) GO TO 8
 rm(4) = 1.0
 2605 n = 6
 GO TO 3
 
!*****         29-PROD         *****************************************
 
 2900 n = 6
 2903 IF (m(2) <= 0) GO TO 8
 2906 IF (m(1) <= 0) GO TO 8
 GO TO 3
 
!*****         30-PTUBE        *****************************************
 
 2920 n = 5
 IF (rm(3) <= 0.0 .OR. rm(4) < 0.0 .OR. rm(4) > 0.5*rm(3)) GO TO 8
 IF (rm(4) == 0.0) rm(4) = 0.5*rm(3)
 GO TO 2903
 
!*****         33-PTRIA1,38-PQUAD1    **********************************
 
 2980 IF (m(2) < 0 .OR.  m(4) < 0 .OR.  m(6) < 0) GO TO 8
 IF (m(2) == 0 .AND. m(4) == 0 .AND. m(6) == 0) GO TO 8
 DO  l = 2,6,2
   IF (m(l) == 0 .AND. m(l+1) /= 0) GO TO 8
 END DO
 n = 10
 GO TO 2906
 
!*****      34-PTRIA2,37-PTRMEM,39-PQUAD2,41-PQDMEM  *******************
!*****      42-PSHEAR,43-PTWIST,121-PTORDRG,250-PQDMEM1    *************
!*****      260-PQDMEM2
 
 3000 DO  l = 1,5,4
   IF (m(l) == 0 .AND. m(l+1) == 0 .AND. m(l+2) == 0) CYCLE
   IF (m(l) <= 0 .OR.  m(l+1) <= 0) GO TO 8
   IF (rm(l+2) <= 0.0) GO TO 8
   n = n + 4
   i(n-3) = m(l  )
   i(n-2) = m(l+1)
   i(n-1) = m(l+2)
   i(n  ) = m(l+3)
   IF (e(kl) < 0) CYCLE
   IF (m(l) > e(kl)) GO TO 3004
   e(kl) = -m(l)
   CYCLE
   3004 e(kl) = m(l)
 END DO
 IF (n > 0) THEN
   GO TO     2
 ELSE
   GO TO     8
 END IF
 3011 kl = 30
 GO TO 3000
 3012 kl = 31
 GO TO 3000
 3013 kl = 27
 GO TO 3000
 3014 kl = 24
 GO TO 3000
 3015 kl = 28
 GO TO 3000
 3016 kl = 32
 GO TO 3000
 3017 kl = 29
 GO TO 3000
 3018 kl = 25
 GO TO 3000
 3019 kl = 26
 GO TO 3000
 
!*****         35-PTRBSC,36-PTRPLT,40-PQDPLT     ***********************
 
 3020 IF (m(2) < 0 .OR. m(4) < 0 .OR. m(2) == 0 .AND. m(4) == 0) GO TO 8
 DO  l = 2,4,2
   IF (m(l) == 0 .AND. m(l+1) /= 0) GO TO 8
 END DO
 n = 8
 GO TO 2906
 
!*****         44-PMASS,45-PDAMP    ***********************************
 
 3200 DO  l = 1,7,2
   IF (m(l) == 0 .AND. m(l+1) == 0) CYCLE
   IF (m(l) <= 0) GO TO 8
   n = n + 2
   i(n-1) = m(l  )
   i(n  ) = m(l+1)
   IF (e(kl) < 0) CYCLE
   IF (m(l) > e(kl)) GO TO 3204
   e(kl) = -m(l)
   CYCLE
   3204 e(kl) = m(l)
 END DO
 IF (n > 0) THEN
   GO TO     2
 ELSE
   GO TO     8
 END IF
 3210 kl = 23
 GO TO 3200
 3220 kl = 21
 GO TO 3200
 
!*****         46-PELAS             ************************************
 
 3240 DO  l = 1,5,4
   IF (m(l) == 0 .AND. m(l+1) == 0 .AND. m(l+2) == 0) CYCLE
   IF (m(l) <= 0) GO TO 8
   n = n + 4
   i(n-3) = m(l  )
   i(n-2) = m(l+1)
   i(n-1) = m(l+2)
   i(n  ) = m(l+3)
   IF (e(kl) < 0) CYCLE
   IF (m(l) > e(kl)) GO TO 3244
   e(kl) = -m(l)
   CYCLE
   3244 e(kl) = m(l)
 END DO
 IF (n > 0) THEN
   GO TO     2
 ELSE
   GO TO     8
 END IF
 3255 kl = 22
 GO TO 3240
 
!*****         47-CONROD       *****************************************
 
 3260 IF (m(1) <= 0 .OR. m(2) <= 0 .OR. m(3) <= 0 .OR. m(4) <= 0) GO TO 8
 IF (m(2) == m(3)) GO TO 8
 n = 8
 GO TO 3
 
!*****         48-CROD,49-CTUBE,50-CVISC,356-CPSE2    *****************
 
 3280 DO  l = 1,5,4
   IF (m(l) == 0 .AND. m(l+1) == 0 .AND. m(l+2) == 0 .AND. m(l+3) == 0) CYCLE
   IF (m(l) <= 0 .OR. m(l+2) <= 0 .OR. m(l+3) <= 0) GO TO 8
   IF (mf(l+1) == 0) m(l+1) = m(l)
   IF (m(l+1) <= 0 .OR. m(l+2) == m(l+3)) GO TO 8
   n = n + 4
   i(n-3) = m(l  )
   i(n-2) = m(l+1)
   i(n-1) = m(l+2)
   i(n  ) = m(l+3)
   IF (e(kl) < 0) CYCLE
   IF (m(l) > e(kl)) GO TO 3287
   e(kl) = -m(l)
   CYCLE
   3287 e(kl) = m(l)
 END DO
 IF (n > 0) THEN
   GO TO     2
 ELSE
   GO TO     8
 END IF
 3281 kl = 1
 GO TO 3280
 3282 kl = 2
 GO TO 3280
 3283 kl = 3
 GO TO 3280
 
!*****       52-CTRIA1,53-CTRIA2,54-CTRBSC,55-CTRPLT,56-CTRMEM    ****
!            357-CPSE3
 
 3360 IF (m(3) <= 0 .OR. m(4) <= 0 .OR. m(5) <= 0) GO TO 8
 IF (m(3) == m(4) .OR. m(4) == m(5) .OR. m(3) == m(5)) GO TO 8
 n = 6
 IF (k == 357) n = 5
 3370 IF (mf(2) == 0) m(2) = m(1)
 GO TO 2903
 
!*****       57-CQUAD1,58-CQUAD2,59-CQDPLT,60-CQDMEM,249-CQDMEM1    ****
!*****       259-CQDMEM2,358-CPSE4
 
 3460 IF (m(3) <= 0 .OR. m(4) <= 0 .OR. m(5) <= 0 .OR. m(6) <= 0) GO TO 8
 IF (m(3) == m(4) .OR. m(4) == m(5) .OR. m(5) == m(6) .OR.  &
     m(3) == m(5) .OR. m(4) == m(6) .OR. m(3) == m(6)) GO TO 8
 n = 7
 IF (k == 358) n = 6
 GO TO 3370
 
!*****         61-CSHEAR,62-CTWIST    **********************************
 
 3540 IF (m(3) <= 0 .OR. m(4) <= 0 .OR. m(5) <= 0 .OR. m(6) <= 0) GO TO 8
 IF (m(3) == m(4) .OR. m(4) == m(5) .OR. m(5) == m(6) .OR.  &
     m(3) == m(5) .OR. m(4) == m(6) .OR. m(3) == m(6)) GO TO 8
 n = 6
 GO TO 3370
 
!*****         63-CONM1        *****************************************
 
 3580 IF (m(1) < 0 .OR. m(2) <= 0 .OR. m(3) < 0) GO TO 8
 n = 24
 GO TO 3
 
!*****         64-CONM2        *****************************************
 
 3600 IF (m(1) < 0 .OR. m(2) <= 0) GO TO 8
 DO  l = 1,7
   i(l) = m(l)
 END DO
 DO  l = 8,13
   i(l) = m(l+1)
 END DO
 n = 13
 GO TO 2
 
!*****         65-CMASS1,69-CDAMP1,73-CELAS1,70-CDAMP2,66-CMASS2    ****
 
 3620 IF (mf(2) == 0) m(2) = m(1)
 IF (m(2)  <= 0) GO TO 8
 3623 n = 6
 3626 IF (m(1) <= 0) GO TO 8
 IF (m(3) < 0 .OR.  m(4) < 0 .OR. m(5) < 0 .OR. m(6) < 0) GO TO 8
 IF (m(4) > 6 .OR.  m(6) > 6 .OR. m(3) == 0 .AND. m(5) == 0) GO TO 8
 IF (m(3) == 0 .AND. m(4) /= 0 .OR. m(5) == 0 .AND. m(6) /= 0) GO TO 8
 IF (m(3) == m(5) .AND. m(4) == m(6)) GO TO 8
 icell = m(4)
 m(4)  = m(5)
 m(5)  = icell
 GO TO 3
 
!*****         67-CMASS3,75-CELAS3,71-CDAMP3    ************************
 
 3660 DO  l = 1,5,4
   IF (m(l) == 0 .AND. m(l+1) == 0 .AND. m(l+2) == 0 .AND. m(l+3) == 0) CYCLE
   IF (m(l) <= 0 .OR. m(l+2) < 0 .OR. m(l+3) < 0) GO TO 8
   IF (mf(l+1) == 0) m(l+1) = m(l)
   IF (m(l+1) <= 0 .OR. m(l+2) == m(l+3)) GO TO 8
   n = n + 4
   i(n-3) = m(l  )
   i(n-2) = m(l+1)
   i(n-1) = m(l+2)
   i(n  ) = m(l+3)
   IF (e(kl) < 0) CYCLE
   IF (m(l) > e(kl)) GO TO 3667
   e(kl) = -m(l)
   CYCLE
   3667 e(kl) = m(l)
 END DO
 IF (n > 0) THEN
   GO TO     2
 ELSE
   GO TO     8
 END IF
 3674 kl = 4
 GO TO 3660
 3675 kl = 5
 GO TO 3660
 3676 kl = 6
 GO TO 3660
 
!*****         68-CMASS4,76-CELAS4,72-CDAMP4    ************************
 
 3680 DO  l = 1,5,4
   IF (m(l) == 0 .AND. m(l+1) == 0 .AND. m(l+2) == 0 .AND. m(l+3) == 0) CYCLE
   IF (m(l) <= 0 .OR. m(l+2) < 0 .OR. m(l+3) < 0) GO TO 8
   IF (m(l+2) == m(l+3)) GO TO 8
   n = n + 4
   i(n-3) = m(l  )
   i(n-2) = m(l+1)
   i(n-1) = m(l+2)
   i(n  ) = m(l+3)
   IF (e(kl) < 0) CYCLE
   IF (m(l) > e(kl)) GO TO 3687
   e(kl) = -m(l)
   CYCLE
   3687 e(kl) = m(l)
 END DO
 IF (n > 0) THEN
   GO TO     2
 ELSE
   GO TO     8
 END IF
 3697 kl = 7
 GO TO 3680
 3698 kl = 8
 GO TO 3680
 3699 kl = 9
 GO TO 3680
 
!*****         74-CELAS2       *****************************************
 
 3800 n = 8
 GO TO 3626
 
!*****         77-MAT1         *****************************************
 
 3860 IF (m(1) <= 0 .OR. (rm(2) == 0 .AND. rm(3) == 0)) GO TO 8
 IF ((rm(2) < 0. .OR. rm(3) < 0.) .AND. ksystm(78) >= 0) GO TO 8
 n = 12
 IF (m(12) < 0) GO TO 8
 l = 3
 IF (mf(2) == 0 .OR. rm(2) == 0.) l = l - 1
 IF (mf(3) == 0 .OR. rm(3) == 0.) l = l - 1
 IF (mf(4) == 0 .OR. rm(4) == 0.) l = l - 1
 IF (l >= 2) GO TO 3865
 CALL page2 (3)
 WRITE  (nout,3862) uwm,m(1)
 3862 FORMAT (a25,' 2251, TWO OF THE E, G AND NU ON MAT1 CARD ',i8,  &
     ' ARE ZEROS OR BLANKS.', /5X, 'POTENTIAL ERROR MAY OCCUR LATER')
 3865 IF (mf(2) == 2 .AND. mf(3) == 2 .AND. mf(4) == 2) GO TO 3
 IF (mf(2) == 0) rm(2) = 2.0*rm(3)*(1.0+rm(4))
 IF (mf(3) == 0) rm(3) = rm(2)/(2.0*(1.0+rm(4)))
 IF (mf(4) == 0) rm(4) = rm(2)/(2.0*rm(3)) - 1.0
 IF (rm(4) >= -1.0 .AND. rm(4) <= 0.5) GO TO 3
 CALL page2 (2)
 WRITE  (nout,3870) uwm,m(1),rm(4)
 3870 FORMAT (a25,' 2251, PHYSICALLY UNREALISTIC VALUE FOR NU ON MAT1 ',  &
     'CARD ',i8,'.  VALUE = ',1P,e16.4)
 GO TO 3
 
!*****         78-MAT2         *****************************************
 
 3880 n = 17
 IF (m(17) < 0) GO TO 8
 IF (m(1) > 0) THEN
   GO TO     3
 ELSE
   GO TO     8
 END IF
 
!*****     234-MAT4,   337-MATF   **************************************
 
 3900 IF (m(1) <= 0) GO TO 8
 IF (rm(2) <= 0.0) GO TO 8
 IF (rm(3) <= 0.0 .AND. mf(3) == 2) GO TO 8
 n = 3
 GO TO 3
 
!*****         237-MATT4           *************************************
 
 3950 IF (m(1) <= 0) GO TO 8
 IF (m(2) < 0) GO TO 8
 n = 2
 GO TO 3
 
!*****         235-MAT5            *************************************
 
 4000 IF (m(1) <= 0) GO TO 8
 IF (rm(8) <= 0.0 .AND. mf(8) == 2) GO TO 8
 n = 8
 GO TO 3
 
!*****         238-MATT5           *************************************
 
 4050 IF (m(1)  <= 0) GO TO 8
 IF (mf(8) /= 0) GO TO 7
 n = 7
 GO TO 3
 
!*****  315-MATPZ1,   319-MAT6  ****************************************
 
 4060 IF (m(1) <= 0) GO TO 8
 n = 15
 IF (k == 319) n = 31
 GO TO 3
 
!*****    316-MATPZ2        ********************************************
 
 4070 IF (m(1) <= 0) GO TO 8
 n = 52
 GO TO 3
 
!*****     317-MTTPZ1,   320-MATT6   ***********************************
 
 4080 n = 15
 IF (k == 320) n = 31
 DO  l = 1,n
   IF (m(l) < 0) GO TO 8
 END DO
 IF (m(1) == 0) GO TO 8
 GO TO 3
 
!*****   318-MTTPZ2    *************************************************
 
 4090 DO  l = 1,52
   IF (m(l) < 0) GO TO 8
 END DO
 IF (m(1) == 0) GO TO 8
 n = 52
 GO TO 3
 
!*****         223-AXSLOT         **************************************
 
 4100 IF (slot) GO TO 8
 slot = .true.
 iaxf = iaxf + 2
 slotdf(1) = rm(1)
 slotdf(2) = rm(2)
 IF (m(3) < 0) baddat =.true.
 klotdf(3) = m(3)
 slotdf(4) = rm(4)
 IF (m(5) < 0) baddat =.true.
 klotdf(5) = m(5)
 n = 5
 GO TO 3
 
!*****         224-CAXIF2         **************************************
 
 4200 IF (mf(4) /= 0 .OR. mf(5) /= 0) GO TO 7
 n = 3
 4250 IF (m(1)  <= 0) GO TO 8
 IF (mf(6) == 0) rm(6) = slotdf(1)
 IF (mf(7) == 0) rm(7) = slotdf(2)
 IF (mf(8) == 0) m(8)  = klotdf(3)
 DO  l = 2,n
   IF (m(l) <= 0) GO TO 8
   IF (l == 2) CYCLE
   DO  l1 = 3,l
     IF (m(l1-1) == m(l)) GO TO 8
   END DO
 END DO
!     CHECK FOR RHO .GE. 0.0
!     CHECK FOR B .GE. 0.0
!     CHECK FOR N .GE. 0
 DO  l = 6,8
   l1 = l + n - 5
   i(l1) = m(l)
 END DO
 DO  l = 1,n
   i(l) = m(l)
 END DO
 n = n + 3
 GO TO 2
 
!*****         225-CAXIF3         **************************************
 
 4300 IF (mf(5) /= 0) GO TO 7
 n = 4
 GO TO 4250
 
!*****         226-CAXIF4         **************************************
 
 4400 n = 5
 GO TO 4250
 
!*****         227-CSLOT3         **************************************
 
 4500 IF (mf(5) /= 0) GO TO 7
 n = 4
 4550 IF (mf(6) == 0) rm(6) = slotdf(1)
 IF (mf(7) == 0) rm(7) = slotdf(2)
 IF (mf(8) == 0) m(8)  = klotdf(5)
!     CHECK FOR ALL KINDS OF THINGS
 DO  l = 6,8
   l1 = l + n - 5
   i(l1) = m(l)
 END DO
 DO  l = 1,n
   i(l) = m(l)
 END DO
 n = n + 4
 i(n) = klotdf(3)
 GO TO 2
 
!*****         228-CSLOT4         **************************************
 
 4600 n = 5
 GO TO 4550
 
!*****         229-GRIDF          **************************************
 
 4700 IF (m(1) <= 0) GO TO 8
 IF (rm(2) <= 0.0) GO TO 8
 n = 3
 GO TO 3
 
!*****         230-GRIDS          **************************************
 
 4800 IF (m(1) <= 0) GO TO 8
 IF (m(5) < 0) GO TO 8
 IF (mf(4) == 0) rm(4) = slotdf(4)
 n = 5
 GO TO 3
 
!*****         231-SLBDY          **************************************
 
 4900 IF (km /= 0) GO TO 4905
 km = 1
 IF (mf(1) /= 2 .AND. mf(1) /= 0) badfor =.true.
 IF (mf(1) == 0) m(1) = klotdf(1)
 IF (mf(2) /= 1 .AND. mf(2) /= 0) badfor =.true.
 IF (mf(2) == 0) m(2) = klotdf(5)
 IF (m(2) < 0) baddat =.true.
 i(1) = m(1)
 i(2) = m(2)
 n  = 2
 iz = 3
 GO TO 4906
 4905 iz = 1
 4906 DO  l = iz,8
   IF (mf(l) == 0) GO TO 4940
   IF (m(l)  <= 0) baddat =.true.
   n = n + 1
   i(n) = m(l)
 END DO
 4910 IF (m1(1) == 0 .AND. m1(2) == 0) GO TO 4920
 km = 0
 n  = n + 1
 i(n) =-1
 kn = 0
 4920 GO TO 9
 4940 iz = l + 1
 DO   l = iz,8
   IF (mf(l) /= 0) badfor =.true.
 END DO
 IF (m1(1) == 0 .AND. m1(2) == 0) badfor =.true.
 GO TO 4910
 
!*****         232-CHBDY           *************************************
 
 5000 IF (m(1) <= 0) GO TO 8
 i(1) = m(1)
 IF (m(2) < 0) GO TO 8
 i(2) = m(2)
 DO  l = 1,7
   IF (m(3) == hbdynm(1,l) .AND. m(4) == hbdynm(2,l)) GO TO 5020
 END DO
 GO TO 8
 5020 i(3) = l
 l1 = hbdyix(l)
 DO  l2 = 1,l1
   IF (m(l2+4) <= 0 .OR. m(l2+9) < 0) GO TO 8
   i(l2+3) = m(l2+4)
   i(l2+7) = m(l2+9)
 END DO
 IF (l1 == 4) GO TO 5035
 DO  l2 = l1,3
   IF (m(l2+5) /= 0 .OR. m(l2+10) /= 0) GO TO 8
   i(l2+4) = 0
   i(l2+8) = 0
 END DO
 5035 DO  l2 = 12,14
   i(l2) = m(l2+2)
 END DO
 n = 15
 i(15) = m(9)
 GO TO 2
 
!*****         233-QHBDY           *************************************
 
 5050 IF (m(1) <= 0) GO TO 8
 i(1) = m(1)
 DO  l = 1,5
   IF (m(2) == hbdynm(1,l) .AND. m(3) == hbdynm(2,l)) GO TO 5060
 END DO
 GO TO 8
 5060 i(2) = l
 l1 = hbdyix(l)
 DO  l2 = 1,l1
   IF (m(l2+5) <= 0) GO TO 8
   i(l2+4) = m(l2+5)
 END DO
 IF (l1 == 4) GO TO 5075
 DO  l2 = l1,3
   IF (m(l2+6) /= 0) GO TO 8
   i(l2+5) = 0
 END DO
 5075 i(3) = m(4)
 IF (l >= 3 .AND. mf(4) /= 0) GO TO 7
 IF (l < 3 .AND. rm(5) <= 0.0) GO TO 8
 i(4) = m(5)
 n = 8
 GO TO 2
 
!*****         236-PHBDY           *************************************
 
 5100 IF (m(1) <= 0) GO TO 8
 IF (m(2) < 0) GO TO 8
 IF (rm(3) < 0.0) GO TO 8
 IF (rm(4) < 0.0 .OR. rm(4) > 1.0) GO TO 8
 IF (rm(5) < 0.0 .OR. rm(5) > 1.0) GO TO 8
 IF (mf(5) == 0) rm(5) = rm(4)
 n = 7
 GO TO 3
 
!*****         240-QBDY2           *************************************
 
 5150 IF (m(1) <= 0) GO TO 8
 IF (m(2) <= 0) GO TO 8
 n = 6
 GO TO 3
 
!*****                  289-VIEW                 ***************
 
 5175 n = 6
 IF (m(1) > 0) GO TO 3
 GO TO 8
 
!*****         241-QVECT           *************************************
 
 5200 IF (km /= 0) GO TO 5215
 IF (m(1) <= 0) baddat =.true.
 IF (mf(2) /= 2 .AND. mf(2) /= 0) badfor =.true.
 i(1) = m(1)
 i(2) = m(2)
 DO  l = 3,6
   IF (mf(l) == 1) GO TO 5205
   IF (mf(l) /= 2 .AND. mf(l) /= 0) badfor =.true.
   i(l) = m(l)
   CYCLE
   5205 IF (m(l) < 0) baddat =.true.
   i(l) = m(l)
 END DO
 l = 6
 k914 = 209
 GO TO 5216
 5215 l  = 1
 5216 km = 1
 kn = 1
 n  = 6
 l4 = l
 IF (mf(l) /= 1) badfor =.true.
 IF (m(l4) <= 0) baddat =.true.
 5220 IF (l == 8) GO TO 5235
 IF (mf(l)   == 3) GO TO 5225
 IF (mf(l+1) == 0) GO TO 5234
 IF (m(l4) <= 0) baddat =.true.
 i(n) = m(l4)
 l  = l  + 1
 l4 = l4 + 1
 CALL WRITE (k914,i,n,0)
 GO TO 5220
 5225 IF (mf(l+1) /= 1 .OR. m(l4) /= thru) GO TO 5232
 IF (m(l4-1) >= m(l4+2)) GO TO 5232
 l1 = m(l4-1) + 1
 l2 = m(l4+2) - 1
 IF (l2 <= l1) GO TO 5230
 5227 l3 = l1
 i(n) = l3
 CALL WRITE (k914,i,n,0)
 l1 = l1 + 1
 IF (l1 <= l2) GO TO 5227
 5230 l  = l  + 1
 l4 = l4 + 2
 GO TO 5220
 5232 baddat =.true.
 l  = l  + 1
 l4 = l4 + 2
 GO TO 5220
 5234 IF (m1(1) == 0 .AND. m1(2) == 0) badfor =.true.
 5235 IF (mf(l) /= 1) badfor =.true.
 i(n) = m(l4)
 IF (m1(1) == 0 .AND. m1(2) == 0) GO TO 9
 kn = 0
 km = 0
 GO TO 9
 
!*****         257-CYJOIN       ***************************************
 
 5240 IF (km /= 0) GO TO 5253
 IF (m(1) /= 1 .AND. m(1) /= 2) baddat =.true.
 i(1) = m(1)
 IF (mf(2) == 3) GO TO 5242
 IF (mf(2) /= 0) badfor =.true.
 i(2) = blk
 l4 = 3
 GO TO 5244
 5242 IF (m(2) /= bcdc .AND. m(2) /= bcdr .AND. m(2) /= bcds .AND.  &
     m(2) /= it1  .AND. m(2) /= it2  .AND. m(2) /= it3) baddat=.true.
 i(2) = m(2)
 l4 = 4
 5244 km = 1
 i(3) = blk
 n = 3
 l = 3
 k914 = 210
 IF (mf(l) /= 1) badfor =.true.
 IF (m(l4) <= 0) baddat =.true.
 GO TO 5252
 
!*****         258-CNGRNT          *************************************
 
 5245 IF (km /= 0) GO TO 5253
 k914 = 208
 GO TO 5253
 
!*****         243-RADLST          *************************************
 
 5250 IF (km /= 0) GO TO 5253
 IF (idrdl == 1) badfor =.true.
 idrdl = 1
 k914  = 214
 5253 l  = 1
 n  = 0
 5251 km = 1
 l4 = l
 IF (mf(l) /= 1) badfor =.true.
 IF (m(l4) <= 0) baddat =.true.
 5252 IF (l > 8) GO TO 5260
 IF (mf(l) == 0) GO TO 5262
 IF (mf(l) == 3) GO TO 5254
 IF (m(l4) <= 0) baddat =.true.
 IF (n < 49) GO TO 5255
 CALL WRITE (k914,i,n,0)
 n = 0
 5255 n = n + 1
 i(n) = m(l4)
 l  = l  + 1
 l4 = l4 + 1
 GO TO 5252
 5254 IF (l == 8) GO TO 5258
 IF (mf(l+1) /= 1 .OR. m(l4) /= thru) GO TO 5258
 IF (m(l4-1) >= m(l4+2)) GO TO 5258
 l1 = m(l4-1) + 1
 l2 = m(l4+2)
 5256 l3 = l1
 IF (n < 49) GO TO 5257
 CALL WRITE (k914,i,n,0)
 n = 0
 5257 n = n + 1
 i(n) = l3
 l1 = l1 + 1
 IF (l1 <= l2) GO TO 5256
 l  = l  + 2
 l4 = l4 + 3
 GO TO 5252
 5258 baddat =.true.
 l  = l  + 1
 l4 = l4 + 2
 GO TO 5252
 5260 IF (m1(1) == 0 .AND. m1(2) == 0) GO TO 5266
 5261 km = 0
 n  = n + 1
 i(n) =-1
 kn = 0
 GO TO 9
 5262 DO  l2 = l,8
   IF (mf(l2) /= 0) badfor =.true.
 END DO
 IF (m1(1) /= 0 .OR. m1(2) /= 0) GO TO 5261
 badfor =.true.
 5266 kn = 1
 GO TO 9
 
!*****         268-SET1       ******************************************
 
 5300 IF (km /= 0) GO TO 5253
 IF (mf(1) /= 1) badfor =.true.
 i(1) = m(1)
 n = 1
 l = 2
 k914 = 204
 GO TO 5251
 
!*****       291-CTRIM6         ****************************************
 
 6101 IF (m(3) <= 0 .OR. m(4) <= 0 .OR. m(5) <= 0 .OR. m(6) <= 0 .OR.  &
     m(7) <= 0 .OR. m(8) <= 0) GO TO 8
 IF (m(3) == m(4) .OR. m(3) == m(5) .OR. m(3) == m(6) .OR.  &
     m(3) == m(7) .OR. m(3) == m(8) .OR. m(4) == m(5) .OR.  &
     m(4) == m(6) .OR. m(4) == m(7) .OR. m(4) == m(8)) GO TO  8
 IF (m(5) == m(6) .OR. m(5) == m(7) .OR. m(5) == m(8) .OR.  &
     m(6) == m(7) .OR. m(6) == m(8) .OR. m(7) == m(8)) GO TO 8
 DO  l = 1,8
   IF (mf(l) /= 1) GO TO 7
 END DO
 IF (mf(9) /= 0 .AND. mf(9) /= 2) GO TO 7
 IF (mf(2) == 0) m(2) = m(1)
 n = 9
 GO TO 2903
 
!*****       292-PTRIM6         ****************************************
 
 6201 IF (m(2) < 0 .OR. rm(3) < 0.0 .OR. rm(4) < 0.0 .OR.  &
     rm(5) < 0.0) GO TO 8
 IF (rm(3) == 0.0)  GO TO  8
 IF (mf(1) /= 1 .AND. mf(2) /= 1) GO TO 7
 IF (mf(3) /= 2) GO TO 7
 DO  l = 4,6
   IF (mf(l) /= 0 .AND. mf(l) /= 2) GO TO 7
 END DO
 n = 6
 GO TO 2906
 
!*****       293-CTRPLT1        ****************************************
 
 6301 GO TO 6101
 
!*****       294-PTRPLT1        ****************************************
 
 6401 IF (m(2) < 0 .OR. m(6) < 0 .OR. m(2) == 0 .AND. m(6) == 0) GO TO 8
 IF (m(2) == 0  .AND. m(3) /= 0 ) GO TO 8
 IF (m(6) == 0  .AND. m(7) /= 0 ) GO TO 8
 IF (mf(1) /= 1 .AND. mf(2) /= 1) GO TO 7
 IF (mf(6) /= 0 .AND. mf(6) /= 1) GO TO 7
 IF (mf(3) /= 2) GO TO 7
 IF (mf(4) /= 0 .AND. mf(4) /= 2) GO TO 7
 IF (mf(5) /= 0 .AND. mf(5) /= 2) GO TO 7
 DO  l = 7,16
   IF (mf(l) /= 0 .AND. mf(l) /= 2) GO TO 7
 END DO
 n = 16
 GO TO 2906
 
!*****       295-CTRSHL         ****************************************
 
 7501 GO TO 6101
 
!*****       296-PTRSHL         ****************************************
 
 7601 CONTINUE
 IF (m(2) < 0 .OR. m(6) < 0 .OR. m(10) < 0 .OR. m(2) == 0 .AND.  &
     m(6) == 0 .AND. m(10) == 0) GO TO 8
 IF (m(2) == 0 .AND. rm(3) /= 0.0) GO TO 8
 IF (m(6) == 0 .AND. rm(7) /= 0.0) GO TO 8
 IF (m(10) == 0 .AND. rm(11) /= 0.0) GO TO 8
 IF (rm(3) < 0.0 .OR. rm(4) < 0.0 .OR. rm(5) < 0.0) GO TO 8
 IF (rm(7) < 0.0 .OR. rm(8) < 0.0 .OR. rm(9) < 0.0) GO TO 8
 IF (rm(11) < 0.0 .OR. rm(12) < 0.0 .OR. rm(13) < 0.0) GO TO 8
 IF (mf(10) /= 0 .AND. mf(10) /= 1) GO TO 7
 IF (mf(1) /= 1) GO TO 7
 IF (mf(2) /= 0 .AND. mf(2) /= 1) GO TO 7
 IF (mf(6) /= 0 .AND. mf(6) /= 1) GO TO 7
 DO  l = 3,11,4
   IF (mf(l) /= 0 .AND. mf(l) /= 2) GO TO 7
   IF (mf(l+1) /= 0 .AND. mf(l+1) /= 2) GO TO 7
   IF (mf(l+2) /= 0 .AND. mf(l+2) /= 2) GO TO 7
 END DO
 n = 20
 GO TO 2906
 
!*********     359-PPSE       ******************************************
 
 7700 n = 5
 IF (m(1)  <=  0) GO TO 8
 IF (rm(2) == 0.) GO TO 8
 rm(3) = 0.0
 rm(4) = 0.0
 rm(5) = 0.0
 GO TO 3
 
END SUBROUTINE ifs1p
