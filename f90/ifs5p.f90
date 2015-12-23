SUBROUTINE ifs5p (*,*,*)
     
    EXTERNAL        lshift,rshift,orf
    LOGICAL :: abort,baddat,badfor,ifpdco
    INTEGER :: m(100),ret,thru,nfdh(10),itype(12),iscal(4),  &
               orf,rshift,lshift,c,p,t1,blank,met(4),mot(3),gc
    CHARACTER (LEN=25) :: sfm
    CHARACTER (LEN=29) :: uim
    CHARACTER (LEN=25) :: uwm
    CHARACTER (LEN=23) :: ufm

    COMMON /xmssg / ufm,uwm,uim,sfm
    COMMON /ifpx1 / ncds,t1(2,310)
    COMMON /system/ nbuf,nout,abort,junk(42),kdumel(9)
    COMMON /bitpos/ kb(32,2)
    COMMON /ifpdta/ id,n,k,kx,ky,i(100),rm(100),mf(100),m1(100), m1f(100),  &
            kn,baddat,badfor,nopen,nparam,iax,nax,iaxf,naxf,  &
            lharm,knt,slotdf(5),gc(7),ll(6)
    COMMON /cifs5p/ km,c,p,icont,iaero,ipopt

    EQUIVALENCE     (m(1),rm(1)), (BLANK,iblank)
    DATA    thru  / 4HTHRU/
    DATA    BLANK / 1H    /
    DATA    iyes  , ino   /    4HYES , 4HNO   /
    DATA    ms,ml / 4HS   ,    4HL   /
    DATA    mot   / 1HZ,  1HY, 2HZY  /
    DATA    met   / 1HK,  2HPK,2HKE,   3HINV  /
    DATA    nmt   / 4     /
    DATA    itype, iscal  /  &
        4HFX  ,4HFY  ,4HFZ  ,4HFXE ,4HFYE ,4HFZE ,4HMX  ,4HMY  ,  &
        4HMZ  ,4HMXE ,4HMYE ,4HMZE ,4HLE  ,4HFR  ,4HLEPR,4HFRPR/
 
    IF (k > 100) GO TO 81
    GO TO (     5,   5, 100,   5,   5,   5,   5,   5,   5,   5,  &
        5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
        5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
        5, 200,   5,   5,   5,   5,   5,   5,   5,   5,  &
        5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
        300,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
        5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
        5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
        5,   5,   5,   5,   5,   5,   5, 400,   5,   5,  &
        5,   5,   5,   5,   5,   5,   5,   5, 500, 600 ), k
81  IF (kx > 100) GO TO 82
    GO TO (   700,   5, 800,   5,   5, 900,1000,1100,1200,1300,  &
        1400,1500,1600,1700,1800,1900,2000,2100,   5,   5,  &
        5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
        5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
        5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
        5,   5,   5,   5,   5,   5,   5,   5,2200,2300,  &
        2400,   5,2500,2600,2700,   5,2800,2900,3000,3100,  &
        3200,3300,3400,3500,3600,3700,3800,3900,   5,   5,  &
        5,   5,   5,   5,   5,4000,4100,   5,   5,   5,  &
        5,   5,4400,4500,   5,   5,   5,6000,   5,   5 ), kx
82  IF (ky > 100) GO TO 83
    GO TO (     5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
        5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
        5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
        5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
        5,   5,   5,   5,4600,4600,   5,   5,   7,   8,  &
        5000,5100,5200,5300,5400,   5,   5,   5,   5,   5,  &
        5,   5,6400,6500,6600,6700,6800,   5,5600,5700,  &
        5800,5900,   5,   5,6100,6200,6300,6900,   5,   5,  &
        5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
        5,   5,   5,   5,   5,   5,   5,   5,   5,   5 ), ky
83  kz = ky-100
    IF (kz > 39) GO TO 5
    GO TO (  6400,6400,6400,6510,6520,6530,6850,7600,6400,7700,  &
        3300,3300,3300,3350,   5,   5,   5,   5,   5,   5,  &
        5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
        5,   5,   5,   5,   5,   5,   5,4700,4710      ), kz
5   CALL page2 (2)
    WRITE  (nout,6) sfm
6   FORMAT (a25,' 322, ILLEGAL ENTRY TO IFS5P.')
    abort =.true.
    IF (k == 0) GO TO 9999
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
 
!*****         3-ADUM1        ******************************************
 
100 CONTINUE
    idumel = 1
    GO TO 8100
 
!*****         32-ADUM2       ******************************************
 
200 CONTINUE
    idumel = 2
    GO TO 8100
 
!*****         51-ADUM3       ******************************************
 
300 CONTINUE
    idumel = 3
    GO TO 8100
 
!*****         88-ADUM4       ******************************************
 
400 CONTINUE
    idumel = 4
    GO TO 8100
 
!*****         99-ADUM5       ******************************************
 
500 CONTINUE
    idumel = 5
    GO TO 8100
 
!*****         100-ADUM6      ******************************************
 
600 CONTINUE
    idumel = 6
    GO TO 8100
 
!*****         101-ADUM7      ******************************************
 
700 CONTINUE
    idumel = 7
    GO TO 8100
 
!*****         103-ADUM8      ******************************************
 
800 CONTINUE
    idumel = 8
    GO TO 8100
 
!*****         106-ADUM9      ******************************************
 
900 CONTINUE
    idumel = 9
    GO TO 8100
 
!*****         107-CDUM1      ******************************************
 
1000 CONTINUE
     idumel = 1
     GO TO 8200
 
 !*****         108-CDUM2      ******************************************
 
1100 CONTINUE
     idumel = 2
     GO TO 8200
 
 !*****         109-CDUM3      ******************************************
 
1200 CONTINUE
     idumel = 3
     GO TO 8200
 
 !*****         110-CDUM4      ******************************************
 
1300 CONTINUE
     idumel = 4
     GO TO 8200
 
 !*****         111-CDUM5      ******************************************
 
1400 CONTINUE
     idumel = 5
     GO TO 8200
 
 !*****         112-CDUM6      ******************************************
 
1500 CONTINUE
     idumel = 6
     GO TO 8200
 
 !*****         113-CDUM7      ******************************************
 
1600 CONTINUE
     idumel = 7
     GO TO 8200
 
 !*****         114-CDUM8      ******************************************
 
1700 CONTINUE
     idumel = 8
     GO TO 8200
 
 !*****         115-CDUM9      ******************************************
 
1800 CONTINUE
     idumel = 9
     GO TO 8200
 
 !*****         116-PDUM1      ******************************************
 
1900 CONTINUE
     idumel = 1
     GO TO 8300
 
 !*****         117-PDUM2      ******************************************
 
2000 CONTINUE
     idumel = 2
     GO TO 8300
 
 !*****         118-PDUM3      ******************************************
 
2100 CONTINUE
     idumel = 3
     GO TO 8300
 
 !*****         159-PDUM4      ******************************************
 
2200 CONTINUE
     idumel = 4
     GO TO 8300
 
 !*****         160-PDUM5      ******************************************
 
2300 CONTINUE
     idumel = 5
     GO TO 8300
 
 !*****         161-PDUM6      ******************************************
 
2400 CONTINUE
     idumel = 6
     GO TO 8300
 
 !*****         163-PDUM7      ******************************************
 
2500 CONTINUE
     idumel = 7
     GO TO 8300
 
 !*****         164-PDUM8      ******************************************
 
2600 CONTINUE
     idumel = 8
     GO TO 8300
 
 !*****         165-PDUM9      ******************************************
 
2700 CONTINUE
     idumel = 9
     GO TO 8300
 
     !*****         167-CONCT1     ******************************************
 
2800 IF (km == 1) GO TO 2850
     nss = 0
     IF (mf(1) /= 1) GO TO 7
     DO  l = 2,8
         IF (mf(l) /= 3 .AND. mf(l) /= 0) GO TO 7
         IF (mf(l) == 3) nss = nss + 1
     END DO
     IF (m1(1) /= 0 .OR. m1(2) /= 0) GO TO 7
     IF (m(1) <= 0) GO TO 8
     IF (nss  == 1) GO TO 8
     i(1) = nss
     i(2) = m(1)
     n  = 2
     nb = 0
     DO  l = 2,8
         IF (mf(l) == 0) GO TO 2809
         n = n + 2
         nfdh(l-1) = 1
         i(n-1) = m(n-2+nb)
         i(n  ) = m(n-1+nb)
         CYCLE
2809     nb = nb + 1
         nfdh(l-1) = 0
     END DO
     km = 1
     GO TO 2
2850 km = 0
     DO  l = 1,8
         IF (mf(l) > 1) GO TO 7
         IF (m(l) <= 0 .AND. mf(l) == 1) GO TO 8
     END DO
     DO  l = 2,8
         IF (mf(l) == 1 .AND. nfdh(l-1) == 0) GO TO 8
     END DO
     i(1) = m(1)
     n = 1
     DO  l = 2,8
         IF (nfdh(l-1) == 0) CYCLE
         n = n + 1
         i(n) = m(l)
     END DO
     kn = 1
     km = 1
     IF (m1(1) == 0 .AND. m1(2) == 0) GO TO 2
     kn = 0
     km = 0
     n  = n + 1
     i(n) = -1
     GO TO 2
 
     !*****         168-CONCT      ******************************************
 
2900 IF (km == 1) GO TO 2950
     km = 1
     DO  l = 1,2
         IF (mf(l  ) /= 1) GO TO 7
         IF (mf(l+2) /= 3) GO TO 7
         IF (m(l) <= 0) GO TO 8
     END DO
     DO  l = 1,6
         i(l) = m(l)
     END DO
     n = 6
     IF (m1(1) /= 0 .AND. m1(2) /= 0) GO TO 7
     GO TO 2
     2950 DO  l = 1,8
         IF (mf(l) /= 0 .AND. mf(l) /= 1) GO TO 7
     END DO
     DO  l = 1,8
         IF (mf(l) == 1 .AND. m(l) <= 0) GO TO 8
     END DO
     n = 0
     DO  l = 1,8,2
         kdlh = mf(l) + mf(l+1)
         IF (kdlh /= 0 .AND. kdlh /= 2) GO TO 8
         IF (kdlh == 0) CYCLE
         n = n + 2
         i(n-1) = m(n-1)
         i(n  ) = m(n  )
     END DO
     IF (m1(1) == 0 .AND. m1(2) == 0) GO TO 2
     n = n + 2
     i(n-1) = -1
     i(n  ) = -1
     km = 0
     GO TO 2
 
     !*****         169-TRANS      ******************************************
 
3000 IF (mf(1) /= 1 .OR. mf(2) /= 0) GO TO 7
     DO  l = 3,11
         IF (mf(l) /= 2 .AND. mf(l) /= 0) GO TO 7
     END DO
     IF (m(1) <= 0) GO TO 8
     v11   = rm( 6) - rm(3)
     v12   = rm( 7) - rm(4)
     v13   = rm( 8) - rm(5)
     v21   = rm( 9) - rm(3)
     v22   = rm(10) - rm(4)
     v23   = rm(11) - rm(5)
     tr1   = v12*v23 - v13*v22
     tr2   = v11*v23 - v13*v21
     tr3   = v11*v22 - v12*v21
     tmag  = SQRT(tr1**2 + tr2**2 + tr3**2)
     v1mag = SQRT(v11**2 + v12**2 + v13**2)
     v2mag = SQRT(v21**2 + v22**2 + v23**2)
     IF (v1mag == 0.0) GO TO 8
     IF (v2mag == 0.0) GO TO 8
     angsin = tmag/v1mag/v2mag
     IF (angsin < 0.087) GO TO 8
     i(1) = m(1)
     DO  l = 3,11
         i(l-1) = m(l)
     END DO
     n = 10
     GO TO 2
 
     !*****         170-RELES      ******************************************
 
3100 IF (km == 1) GO TO 3170
     km = 1
     IF (mf(1) /= 1) GO TO 7
     IF (mf(2) /= 3) GO TO 7
     IF (m(1)  <= 0) GO TO 8
     i(1) = m(1)
     i(2) = m(2)
     i(3) = m(3)
     l1 = 3
     n  = 3
     3180 DO  l = l1,8,2
         kdlh = mf(l) + mf(l+1)
         IF (kdlh /= 0 .AND. kdlh /= 2) GO TO 8
         IF (kdlh == 0) CYCLE
         n = n + 2
         i(n-1) = m(n-1)
         i(n  ) = m(n  )
     END DO
     IF (m1(1) == 0 .AND. m1(2) == 0) GO TO 2
     n = n + 2
     i(n-1) = -1
     i(n  ) = -1
     km = 0
     GO TO 2
3170 n  = 0
     l1 = 1
     GO TO 3180
 
     !*****         171-LOADC      ******************************************
 
3200 IF (km == 1) GO TO 3250
     km = 1
     IF ((mf(1) /= 0 .AND. mf(1) /= 1) .OR.  &
         (mf(2) /= 0 .AND. mf(2) /= 2)) GO TO 7
     IF (m(1) <= 0 .OR. m(2) == 0) GO TO 8
     IF (mf(3) /= 3 .OR. (mf(6) /= 3 .AND. mf(6) /= 0)) GO TO 7
     i(1) = m(1)
     i(2) = m(2)
     n   = 2
     ldh = 0
     3260 DO  l = 3,8,3
         kdlh = mf(l) + mf(l+1) + mf(l+2)
         IF (kdlh /= 0 .AND. kdlh /= 6) GO TO 8
         IF (kdlh == 0) CYCLE
         n = n + 4
         i(n-3) = m(n-3+ldh)
         i(n-2) = m(n-2+ldh)
         i(n-1) = m(n-1+ldh)
         i(n)   = m(n  +ldh)
     END DO
     IF (m1(1) == 0 .AND. m1(2) == 0) GO TO 2
     n = n + 4
     i(n-3) = iblank
     i(n-2) = iblank
     i(n-1) = -1
     i(n  ) = -1
     km = 0
     GO TO 2
3250 n  = 0
     ldh = 2
     GO TO 3260
 
     !*****         172-SPCSD ,  311-DAREAS          **********************
     !              312-DELAYS,  313-DPHASES
 
3300 IF (m(1) <= 0) GO TO 8
     IF (m(4) <= 0) GO TO 8
     IF (m(5) < 0) GO TO 8
     IF (m(7) < 0) GO TO 8
     IF (m(8) < 0) GO TO 8
     n = 12
     IF (m(7) == 0 ) n = 9
     m(n-2) = -1
     m(n-1) = -1
     m(n  ) = -1
     GO TO 3
 
     !*****         314-TICS         **************************************
 
3350 IF (m(1) <= 0) GO TO 8
     IF (m(4) <= 0) GO TO 8
     IF (m(5) < 0) GO TO 8
     DO  l = 8,11
         m(l) = -1
     END DO
     n = 11
     GO TO 3
 
     !*****         173-SPCS1      ******************************************
 
3400 IF (km == 1) GO TO 3410
     km = 1
     IF (mf(1) /= 1) badfor =.true.
     IF (mf(2) /= 3) badfor =.true.
     IF (m(4)  < 0) baddat =.true.
     CALL WRITE (210,m,4,0)
     j1 = 4
     l1 = 5
     GO TO 3920
3410 l1 = 1
     j1 = 1
     GO TO 3920
 
     !*****         174-SPCS       ******************************************
 
 
     !     SAME AS RELES DATA CARD
 
3500 GO TO 3100
 
     !*****         175-BDYC       ******************************************
3600 IF (km == 1) GO TO 3650
 
     IF (mf(8) /= 0 .OR. mf(1) /= 1) GO TO 7
     IF (m(1) <= 0) GO TO 8
     3660 DO  l = 2,7,2
         IF (mf(  l) /= 0 .AND. mf(l  ) /= 3) GO TO 7
         IF (mf(l+1) /= 0 .AND. mf(l+1) /= 1) GO TO 7
     END DO
     i(1) = m(1)
     n  = 1
     j1 = 1
     IF (km == 1) j1 = 0
     DO  l = 2,7,2
         kdlh = mf(l) + mf(l+1)
         IF (kdlh /= 0 .AND. kdlh /= 4) GO TO 8
         IF (kdlh == 0) CYCLE
         n  = n  + 3
         j1 = j1 + 3
         i(j1-2) = m(n-2)
         i(j1-1) = m(n-1)
         i(j1  ) = m(n  )
     END DO
     n  = j1
     km = 1
     IF (m1(1) == 0 .AND. m1(2) == 0) GO TO 2
     km = 0
     n  = n  + 3
     j1 = j1 + 3
     i(j1-2) = iblank
     i(j1-1) = iblank
     i(j1  ) = -1
     GO TO 2
3650 IF (mf(1) /= 0 .OR. mf(8) /= 0) GO TO 7
     GO TO 3660
 
     !*****         176-MPCS       ******************************************
 
3700 IF (km == 1) GO TO 3750
     km = 1
     IF (mf(1) /= 1) GO TO 7
     IF (mf(2) /= 3) GO TO 7
     IF (mf(3) /= 1) GO TO 7
     IF (mf(4) /= 1) GO TO 7
     IF (mf(5) /= 2) GO TO 7
     IF (m(1)  <= 0) GO TO 8
     IF (m(4)  <= 0) GO TO 8
     IF (m(5)  < 0) GO TO 8
     IF (m(6)  == 0) GO TO 8
     DO  l = 1,6
         i(l) = m(l)
     END DO
     n = 6
     IF (m1(1) == 0 .AND. m1(2) == 0) GO TO 2
     GO TO 7
3750 IF (mf(1) /= 0) GO TO 7
     IF (mf(2) /= 3) GO TO 7
     DO  l = 3,6,3
         IF (mf(l)+mf(l+2)+mf(l+1) == 0) CYCLE
         IF (mf(l) /= 1 .OR. mf(l+1) /= 1) GO TO 7
         IF (mf(l+2) /= 2) GO TO 7
         IF (m(l+1) <= 0 .AND. mf(l+2) <= 0) GO TO 8
     END DO
     n = 0
     DO  l = 3,8,3
         kdlh = mf(l) + mf(l+1) + mf(l+2)
         IF (kdlh /= 0 .AND. kdlh /= 4) GO TO 8
         IF (kdlh == 0) CYCLE
         i(n+1) = m(2)
         i(n+2) = m(3)
         n = n + 5
         i(n-2) = m(l+1)
         i(n-1) = m(l+2)
         i(n  ) = m(l+3)
     END DO
     IF (m1(1) == 0 .AND. m1(2) == 0) GO TO 2
     i(n+1) = iblank
     i(n+2) = iblank
     n = n + 5
     i(n-2) = -1
     i(n-1) = -1
     i(n  ) = -1
     km = 0
     GO TO 2
 
     !*****         177-BDYS       ******************************************
 
     3800 DO  l = 1,7
         IF (mf(l) /= 1 .AND. mf(l) /= 0) GO TO 7
         IF (mf(l) == 1 .AND.  m(l) <= 0) GO TO 8
     END DO
     IF (mf(1) == 0) GO TO 7
     n = 1
     i(n) = m(1)
     DO  l = 2,7,2
         kdlh = mf(l) + mf(l+1)
         IF (kdlh /= 2 .AND. kdlh /= 0) GO TO 8
         IF (kdlh == 0) CYCLE
         n = n + 2
         i(n-1) = m(n-1)
         i(n  ) = m(n  )
     END DO
     n = n + 2
     i(n-1) = -1
     i(n  ) = -1
     GO TO 2
 
     !*****         178-BDYS1      ******************************************
 
3900 IF (km == 1) GO TO 3910
     km = 1
     IF (mf(1) /= 1 .OR. mf(2) > 1) badfor =.true.
     IF (m(1) < 1 .OR. m(2) < 0) baddat =.true.
     CALL WRITE (210,m,2,0)
     j1 = 3
     l1 = 3
     GO TO 3920
3910 j1 = 1
     l1 = 1
 
     !     COMMON PROCESSING FOR SPCS1 AND BDYS1 CARDS
 
3920 IF (mf(j1) /= 0) GO TO 3925
     j1 = j1 + 1
     l1 = l1 + 1
     GO TO 3960
3925 IF (mf(j1) == 1) GO TO 3930
     badfor =.true.
     GO TO 3965
3930 IF (j1 > 6) GO TO 3955
     IF (mf(j1+1) /= 3) GO TO 3955
     IF (m(l1+1) == thru) GO TO 3935
     baddat =.true.
     GO TO 3965
3935 IF (mf(j1+2) == 1) GO TO 3940
     badfor =.true.
     GO TO 3965
3940 IF (m(l1+3) > m(l1)) GO TO 3945
     baddat =.true.
     GO TO 3965
3945 ig1 = m(l1  )
     ig2 = m(l1+3)
     DO  j = ig1,ig2
         CALL WRITE (210,j,1,0)
     END DO
     j1 = j1 + 3
     l1 = l1 + 4
     GO TO 3960
3955 CALL WRITE (210,m(l1),1,0)
     j1 = j1 + 1
     l1 = l1 + 1
3960 IF (j1 <= 8) GO TO 3920
3965 IF (m1(1) /= 0 .OR. m1(2) /= 0) GO TO 3970
     kn = 1
     n  = 0
     GO TO 9
3970 km = 0
     kn = 0
     n  = 1
     i(1) = -1
     GO TO 9
 
     !*****         186-GNEW       ******************************************
 
4000 IF (mf(1) /= 1) GO TO 7
     IF (mf(2) /= 3) GO TO 7
     IF (mf(3) /= 1 .AND. mf(3) /= 0) GO TO 7
     IF (mf(4) /= 1) GO TO 7
     IF (mf(5) /= 1) GO TO 7
     IF (m(1)  <= 0) GO TO 8
     IF (m(4)  < 0) GO TO 8
     IF (m(5)  <= 0) GO TO 8
     IF (m(6)  <= 0) GO TO 8
     n = 6
     GO TO 3
 
     !*****         187-GTRAN      ******************************************
 
4100 IF (mf(1) /= 1) GO TO 7
     IF (mf(2) /= 3) GO TO 7
     IF (mf(3) /= 1) GO TO 7
     IF (mf(4) /= 1 .AND. mf(4) /= 0) GO TO 7
     IF (m(1)  <= 0) GO TO 8
     IF (m(4)  <= 0) GO TO 8
     IF (m(5)  < 0) GO TO 8
     n = 5
     GO TO 3
 
     !*****         193-USET       ******************************************
 
4400 ASSIGN 4405 TO ret
4401 n = 0
     IF (m(2) /= BLANK) GO TO 8
     DO  l = 1,32
         IF (m(1) == kb(l,2)) GO TO 4404
     END DO
     GO TO 8
4404 id = kb(l,1)
     GO TO ret, (4405,4505)
     4405 DO  l = 3,7,2
         IF (m(l) == 0 .AND. m(l+1) == 0) CYCLE
         IF (m(l) <= 0) GO TO 8
         IF (ifpdco(m(l+1))) GO TO 8
         lz = 6
         IF (m(l+1) == 0) lz = 1
         DO  l2 = 1,lz
             IF (lz /= 1 .AND. ll(l2) == 0) CYCLE
             n = n + 3
             i(n-2) = id
             i(n-1) = m(l  )
             i(n  ) = ll(l2)
             IF (n <= 3) CYCLE
             DO  l1 = 6,n,3
                 IF (i(n-1) == i(l1-4) .AND. i(n) == i(l1-3)) GO TO 8
             END DO
         END DO
     END DO
     IF (n > 0) THEN
         GO TO     2
     ELSE
         GO TO     8
     END IF
 
     !*****         194-USET1      ******************************************
 
4500 IF (km /= 0) GO TO 4510
     km = 1
     ASSIGN 4505 TO ret
     GO TO 4401
4505 n = 2
     i(1) = id
     IF (mf(2) /= 0 .AND. mf(2) /= 1) badfor =.true.
     IF (ifpdco(m(3))) baddat =.true.
     i(2) = m(3)
     IF (mf(4) == 3 .AND. m(5) == thru) GO TO 4550
     l1 = 4
     l3 =-1
     l2 = 9
     GO TO 4511
4510 l1 = 1
     l3 = 0
     l2 = 8
     4511 DO  l = l1,l2
         IF (mf(l+l3) /= 0 .AND. mf(l+l3) /= 1) badfor =.true.
     END DO
     DO  l = l1,l2
         IF (mf(l+l3) == 1) GO TO 4525
     END DO
     baddat =.true.
     4525 DO  l = l1,l2
         IF (m(l) < 0) THEN
             GO TO  4535
         ELSE IF (m(l) == 0) THEN
             GO TO  4540
         END IF
4530     n = n + 1
         i(n) = m(l)
         CYCLE
4535     baddat =.true.
4540 CONTINUE
     END DO
     kn = 1
     IF (m1(1) == 0 .AND. m1(2) == 0) GO TO 9
4545 km = 0
     n  = n + 1
     i(n) = -1
     kn = 0
     GO TO 9
4550 IF (m1(1) /= 0 .OR. m1(2) /= 0) GO TO 4555
     kn = 1
     badfor =.true.
     GO TO 9
4555 IF (mf(3) /= 1 .OR. mf(5) /= 1 ) badfor =.true.
     IF (m(4) <= 0 .OR. m(7) <= m(4)) baddat =.true.
     DO  l = 1,3
         IF (mf(l+5) /= 0) badfor =.true.
     END DO
     IF (badfor .OR. baddat) GO TO 4545
     CALL WRITE (210,i,2,0)
     l1 = m(4)
     l2 = m(7)
     DO  l = l1,l2
         CALL WRITE (210,l,1,0)
     END DO
     n = 0
     GO TO 4545
 
 !*****         245-SAME 246-NOSAME          ****************************
 
4600 CONTINUE
     ialt = 1
     IF (m(3) == thru) ialt = 3
     kdx  = ialt + icont
     SELECT CASE ( kdx )
         CASE (    1)
             GO TO 4620
         CASE (    2)
             GO TO 4620
         CASE (    3)
             GO TO 4630
         CASE (    4)
             GO TO 4640
     END SELECT
     4620 DO  in1 = 1,8,2
         in2 = in1 + 1
         IF (mf(in1) == 0 .AND. mf(in2) == 0) CYCLE
         IF (mf(in1) /= 1 .OR.  mf(in2) /= 1) badfor =.true.
         IF (m(in1) <= 0 .OR.  m(in2) <= 0) baddat =.true.
   
         n = n + 2
         i(n-1) = m(in1)
         i(n  ) = m(in2)
     END DO
     GO TO 4680
 
4630 IF (mf(1) /= 1 .OR. mf(2) /= 1 .OR. mf(3) /= 3 .OR. mf(4) /= 1 .OR.  &
         mf(5) /= 1 .OR. mf(6) /= 1 .OR. mf(7) /= 3 .OR. mf(8) /= 1)  &
         badfor =.true.
     IF (m(1) <= 0 .OR. m(2) <= 0 .OR. m(5) <= 0 .OR. m(6) <= 0 .OR.  &
         m(7) <= 0 .OR. m(8) /= thru .OR. m(10) <= 0) baddat =.true.
     IF (m(5) <= m(2) .OR. m(10) <= m(7)) baddat =.true.
     irange = m(5) - m(2)
     IF ((m(10)-m(7)) /= irange)  baddat =.true.
     i(1) = -1
     i(2) = irange + 1
     i(3) = m(1)
     i(4) = m(2)
     i(5) = m(6)
     i(6) = m(7)
     n = 6
     GO TO 4680
     4640 DO  in1 = 1,6,5
         in2 = in1 + 1
         in3 = in2 + 1
         in4 = in3 + 1
         in5 = in4 + 1
         IF (mf(in1) == 0 .AND. mf(in2) == 0 .AND. mf(in3) == 0 .AND.  &
             mf(in4) == 0) CYCLE
         IF (mf(in1) /= 1 .OR. mf(in2) /= 1 .OR. mf(in3) /= 3 .OR.  &
             mf(in4) /= 1) badfor = .true.
         IF (m(in1) <= 0 .OR. m(in2) <= 0 .OR. m(in3) /= thru .OR.  &
             m(in5) <= 0) baddat =.true.
         IF (m(in5) <= m(in2) .OR. m(in5)-m(in2) /= irange) baddat =.true.
         i(n+1) = m(in1)
         i(n+2) = m(in2)
         n = n + 2
     END DO
 
4680 IF (m1f(1) == 0 .AND. m1f(2) == 0) GO TO 4685
     icont  = 0
     i(n+1) =-1
     i(n+2) =-1
     n = n + 2
     GO TO 9
4685 icont = 1
     GO TO 9
 
     !*****         338-CELBOW     ******************************************
 
4700 IF (m(2) == 0) m(2) = m(1)
     n = 8
     GO TO 3
 
     !*****         339-PELBOW     ******************************************
 
4710 n = 24
     GO TO 3
 
     !*****         251-CIHEX1     ******************************************
 
5000 n = 10
     5010 DO  l = 1,n
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
 
     !*****         252-CIHEX2     ******************************************
 
5100 n = 22
     GO TO 5010
 
     !*****         253-CIHEX3     ******************************************
 
5200 n = 34
     GO TO 5010
 
     !*****         254-PIHEX      ******************************************
 
5300 n = 7
     IF (m(1) <= 0 .OR. m(2) <= 0) GO TO 8
     IF (m(3) < 0) GO TO 8
     IF ((m(4) < 2 .OR. m(4) > 4) .AND. m(4) /= 0) GO TO 8
     DO  l = 5,7
         IF (mf(l) ==   0) GO TO 5310
         IF (mf(l) /=   2) GO TO 7
         IF (rm(l) < 0.0) GO TO 8
         CYCLE
5310     rm(l) = -1.0
     END DO
     IF (rm(5) >= 0.0  .AND. rm(5) < 1.0  ) GO TO 8
     IF (rm(6) > 180.0 .OR. rm(7) > 180.0) GO TO 8
     GO TO 3
 
     !*****         255-PLOAD3     ******************************************
 
5400 IF (m(1) <= 0) GO TO 8
     DO  l = 3,6,3
         IF (m(l) == 0 .AND. m(l+1) == 0 .AND. m(l+2) == 0) CYCLE
         IF (m(l) < 0 .OR.  m(l+1) < 0 .OR.  m(l+2) < 0) GO TO 8
         n = n + 5
         i(n-4) = m(1)
         i(n-3) = m(2)
         i(n-2) = m(l)
         i(n-1) = m(l+1)
         i(n  ) = m(l+2)
     END DO
     IF (n > 0) THEN
         GO TO     2
     ELSE
         GO TO     8
     END IF
 
     !*****    263-CAERO1, 301-CAERO2, 302-CAERO3, 303-CAERO4  *******
     !         309-CAERO5
 
6400 IF (m(1) <= 0) GO TO 8
     IF (m(2) <= 0) GO TO 8
     DO  l = 3,8
         IF (m(l) < 0) GO TO 8
     END DO
     IF (k == 302) GO TO 6410
     IF (k == 303) GO TO 6420
     IF (k == 309) GO TO 6420
     IF (m(4) == 0 .AND. m(6) == 0) GO TO 8
     IF (m(5) == 0 .AND. m(7) == 0) GO TO 8
     IF (m(8) <= 0) GO TO 8
6405 IF (rm(12) < 0.0) GO TO 8
     IF (rm(16) < 0.0) GO TO 8
     IF (rm(12) == 0.0 .AND. rm(16) == 0.0) GO TO 8
     n = 16
     GO TO 3
 
     !*****     CAERO3      ************************************************
 
6410 IF (m(4)   == 0 ) GO TO 8
     IF (rm(12) == 0.) GO TO 8
     GO TO 6405
 
     !*****     CAERO4   CAERO5    ******************************************
 
6420 IF (m(4) == 0 .AND. m(5) == 0) GO TO 8
     IF (m(6) > 2) GO TO 8
     GO TO 6405
 
     !*****         264-PAERO1     ******************************************
 
6500 IF (m(1) <= 0) GO TO 8
     DO  l = 2,8
         IF (m(l) < 0) GO TO 8
     END DO
     n = 8
     GO TO 3
 
     !*****      304 - PAERO2    ***************
 
6510 IF (m(1) <= 0) GO TO 8
     DO  l = 1,3
         IF (m(2) == mot(l)) GO TO 6512
     END DO
     GO TO 8
6512 IF (rm(4) <= 0.0) GO TO 8
     IF (rm(5) <= 0.0) GO TO 8
     DO  l = 6,15
         IF (m(l) < 0) GO TO 8
     END DO
     n = 15
     GO TO 3
 
     !*****      305 - PAERO3    ****************
 
6520 IF (m(1) <= 0 .OR. m(2) <= 0 .OR. m(3) < 0) GO TO 8
     IF (m(2) > 50) GO TO 8
     n = 0
     IF (m(3) == 0) n = 4
     IF (m(3) == 1) n = 12
     IF (m(3) == 2) n = 16
     IF (n == 0) GO TO 8
     m(4) = n
     n = n + 4
     IF (n == 8) GO TO 6522
     DO  l = 9,n
         IF (mf(l) == -32767) GO TO 8
     END DO
     IF (rm(12) < rm(10)) GO TO 8
     IF (rm(16) < rm(14)) GO TO 8
     IF (n ==  16) GO TO 6522
     IF (rm(20) < rm(18)) GO TO 8
6522 GO TO 3
 
     !*****      306 - PAERO4   **********************
 
6530 IF (km /= 0) GO TO 6535
     km = 1
     IF (mf(1) /= 1 .OR. m(1) <= 0) GO TO 6540
     DO  l = 2,5
         IF (mf(1) < 0 .OR. mf(l) > 1) GO TO 6540
     END DO
     IF (m(3) < 0) GO TO 6540
     IF (m(2) == 0 .AND. m(3) /= 0) GO TO 6540
     IF (m(2) > 0 .AND. m(3) == 0) GO TO 6540
     IF (m(2) /= 0 .AND. m(4) /= 0) GO TO 6540
     IF (m(4) < 0 .OR.  m(4) > 3) GO TO 6540
     IF (m(4) == 0 .AND. m(5) /= 0) GO TO 6540
     IF (m(4) > 0 .AND. m(5) == 0) GO TO 6540
     DO  l = 1,5
         i(l) = m(l)
     END DO
     n  = 5
     l1 = 6
     GO TO 6533
6535 l1 = 1
     6533 DO  l = l1,8
         IF (mf(l) ==  0) GO TO 6550
         IF (mf(l) /=  2) GO TO 6540
         IF (rm(l) < 0.) GO TO 6540
         n = n + 1
         i(n) = m(l)
     END DO
6539 kn = 1
     IF (m1(1) == 0 .AND. m1(2) == 0) GO TO 9
     kn = 0
     km = 0
     n  = n + 1
     i(n) = -1
     GO TO 9
6540 baddat = .true.
     GO TO 6539
6550 IF (m1(1) == 0 .AND. m1(2) == 0) baddat = .true.
     GO TO 6539
 
     !*****   310 - PAERO5  ************
 
7700 IF (km /= 0) GO TO 6535
     km = 1
     DO  l = 1,3
         IF (mf(l) /= 1 .OR. m(l) <= 0) GO TO 6540
     END DO
     DO  l = 4,7
         IF (mf(l) < 0 .OR. mf(l) > 1) GO TO 6540
     END DO
     IF (m(4) /= 0 .AND. m(5) == 0) GO TO 6540
     IF (m(6) /= 0 .AND. m(7) == 0) GO TO 6540
     DO  l = 1,7
         i(l) = m(l)
     END DO
     n = 7
     GO TO 6539
 
     !*****         265-AERO       ******************************************
 
6600 IF (iaero /= 0) GO TO 8
     iaero = 1
     IF (m(1) < 0) GO TO 8
     n = 6
     GO TO 3
 
     !*****         266-SPLINE1    ******************************************
 
6700 IF (m(2) <= 0 .OR. m(3) <= 0 .OR. m(4) <= 0 .OR. m(5) <= 0 .OR.  &
         m(1) <= 0 .OR. rm(6) < 0.0) GO TO 8
     n = 6
     GO TO 3
 
     !*****         267-SPLINE2    ******************************************
 
6800 IF (m(1) <= 0 .OR. m(2) <= 0 .OR. m(3) <= 0 .OR. m(4) <= 0 .OR.  &
         m(5) <= 0 .OR.m(8) < 0) GO TO 8
     n = 10
     GO TO 3
 
     !*****    307 - SPLINE3       *********************
 
6850 IF (km /= 0) GO TO 6852
     km = 1
     IF (mf(1) /= 1 .OR. mf(2) /= 1 .OR. mf(3) /= 1 .OR. mf(4) /= 1) GO TO 6540
     IF (m(2) <= 0 .OR. m(3) < 0) GO TO 6540
     IF (ifpdco(m(4))) GO TO 6540
     IF (gc(2) /= 0) GO TO 6540
     DO  l = 1,4
         i(l) = m(l)
     END DO
     n  = 4
     l1 = 5
     GO TO 6853
6852 l1 = 1
     6853 DO  l = l1,8,4
         IF (mf(l  ) == 0) CYCLE
         IF (mf(l  ) /= 1) GO TO 6540
         IF (mf(l+1) /= 1) GO TO 6540
         IF (ifpdco(m(l+1))) GO TO 6540
         IF (gc(2  ) /= 0) GO TO 6540
         IF (mf(l+2) /= 2) GO TO 6540
         IF (m(l)    <= 0) GO TO 6540
         n = n + 3
         i(n  ) = m(l+2)
         i(n-1) = m(l+1)
         i(n-2) = m(l  )
     END DO
     GO TO 6539
 
     !*****         269-SET2       ******************************************
 
5600 IF (m(1) <= 0 .OR. m(2) <= 0) GO TO 8
     n = 8
     GO TO 3
 
     !*****         270-MKAERO2    ******************************************
 
5700 n = 0
     DO   l = 2,8,2
         IF (mf(l) == 0 .AND. mf(l-1) == 0) CYCLE
         IF (mf(l) == 0 .OR.  mf(l-1) == 0) GO TO 7
         n = n + 2
         i(n-1) = m(l-1)
         IF (rm(l) <= 0.0) GO TO 8
         i(n) = m(l)
     END DO
     IF (n == 0) GO TO 8
     GO TO 2
 
     !*****         271-MKAERO1    ******************************************
 
5800 IF (mf(1) /= 2 .OR. mf(9) /= 2) GO TO 7
     IF (rm(9) <= 0.0) GO TO 8
     DO  l = 2,8
         IF (mf(l) == 0) m(l) = -1
         IF (mf(l+8) /= 0 .AND. rm(l+8) <= 0.0) GO TO 8
         IF (mf(l+8) == 0) m(l+8) = -1
     END DO
     n = 16
     GO TO 3
 
     !*****         257-FLUTTER    ******************************************
 
5900 IF (m(1) <= 0 .OR. m(4) < 0 .OR. m(5) < 0 .OR. m(6) < 0) GO TO 8
     DO  l = 1,nmt
         IF (m(2) == met(l)) GO TO 5920
     END DO
     GO TO 8
5920 CONTINUE
     IF (m(7) /= ms .AND. m(7) /= ml) GO TO 8
     n = 10
     GO TO 3
 
     !******    308 - GUST
 
7600 IF (m(1) <= 0 .OR. m(2) <= 0) GO TO 8
     IF (rm(3) == 0.0 .OR. rm(5) == 0.0) GO TO 8
     n = 5
     GO TO 3
 
     !*****         198-PLOAD1       ****************************************
 
6000 IF (m(1) <= 0 .OR. m(2) <= 0) GO TO 8
     i(1) = m(1)
     i(2) = m(2)
     DO  l = 1,12
         IF (m(3) == itype(l)) GO TO 6020
     END DO
     GO TO 8
6020 i(3) = l
     DO  l = 1,4
         IF (m(5) == iscal(l)) GO TO 6040
     END DO
     GO TO 8
6040 i(4) = l
     IF (rm(9) == 0.0) rm(9) = rm(7)
     IF (rm(9) < rm(7)) GO TO 8
     DO  l = 7,10
         i(l-2) = m(l)
     END DO
     n = 8
     GO TO 2
 
     !*****         275-CBARAO       ****************************************
 
6100 IF (m(1) <= 0) GO TO 8
     i(1) = m(1)
     DO  l = 1,2
         IF (m(2) == iscal(l)) GO TO 6120
     END DO
     GO TO 8
6120 i(2) = l
     DO  l = 4,9
         i(l-1) = m(l)
     END DO
     n = 9
     IF (mf(3) == 2) GO TO 6140
     IF (mf(3) /= 1) GO TO 7
     IF (i(3)  <= 0) GO TO 8
     IF (i(3) > 20) i(3) = 20
     IF (rm(5) <= 0.0 .OR. rm(6) <= 0.0) GO TO 8
     i(9) = -1
     GO TO 2
6140 i(9) = 1
     DO  l = 4,9
         IF (rm(l) < 0.0) GO TO 8
     END DO
     GO TO 2
 
     !*****         276-PLIMIT       ****************************************
 
6200 IF (mf(1) /= 3) GO TO 7
     IF (mf(2) /= 2 .AND. mf(2) /= 0) GO TO 7
     IF (rm(3) < 0.0) GO TO 8
     IF (rm(3) == 0.0 .AND. rm(4) == 0.0) GO TO 8
     IF (rm(4) == 0.0) GO TO 6210
     IF (mf(3) /= 2 .OR. rm(4) <= rm(3)) GO TO 8
6210 IF (mf(5) == 3 ) GO TO 6230
     DO  l = 4,8
         IF (mf(l) /= 0 .AND. mf(l) /= 1) GO TO 7
         IF (m(l+1) < 0) GO TO 8
     END DO
     GO TO 6240
6230 IF (m(6) /= thru) GO TO 8
     IF (mf(4) /= 1 .OR. mf(6) /= 1) GO TO 7
     IF (m(8) <= m(5)) GO TO 8
6240 n = 9
     GO TO 3
 
     !*****         277-POPT         ****************************************
 
6300 IF (m(1) <= 0 .OR. m(4) == 0) GO TO 8
     IF (ipopt /= 0) GO TO 8
     ipopt = 1
     IF (rm(2) < 0.0) GO TO 8
     IF (rm(3) <= 0.0) GO TO 8
     IF (m(5) /= iyes .AND. m(5) /= ino) GO TO 8
     n = 6
     GO TO 3
 
     !******       278  PLOADX   ******************************************
 
6900 IF (m(1) <= 0) GO TO 8
     IF (m(4) <= 0 .OR. m(5) <= 0 .OR. m(6) <= 0) GO TO 8
     n = 6
     GO TO 3
 
 
 
 !     ******************************************************************
 
 !     PROCESS ADUM-I CARDS.
 
8100 CONTINUE
     IF (m(1) <= 0) GO TO 8
     IF (m(2) < 0) GO TO 8
     IF (m(3) < 0) GO TO 8
     IF (m(4) /= 3 .AND. m(4) /= 6) GO TO 8
     IF (mf(5) /= 0 .OR. mf(6) /= 0 .OR. mf(7) /= 0 .OR. mf(8) /= 0) GO TO 7
     kdumel(idumel) = m(4) + 10*(m(3) + 1000*(m(2) + 1000*m(1)))
 
     !     PUT IN CONNECTION AND PROPERTY CARD NAME IF SUPPLIED BY USER
 
     IF (mf(5) /= 3) GO TO 8150
     nbpc = junk(36)
     ncpw = junk(38)
     nsht = nbpc*(ncpw-1)
     nm1  = t1(1,k)
     nm2  = t1(2,k)
     nm1  = rshift(lshift(nm1,nbpc),nbpc)
     c    = lshift(rshift(c,nsht),nsht)
     nm1  = orf(nm1,c)
     p    = lshift(rshift(p,nsht),nsht)
     DO  l = 1,ncds
         IF (nm1 == t1(1,l) .AND. nm2 == t1(2,l)) GO TO 8120
     END DO
     GO TO 8150
8120 t1(1,l) = m(5)
     t1(2,l) = m(6)
     nm1 =  orf(p,rshift(lshift(nm1,nbpc),nbpc))
     DO  l = 1,ncds
         IF (nm1 == t1(1,l) .AND. nm2 == t1(2,l)) GO TO 8140
     END DO
     GO TO 8150
8140 m(5) = orf(p,rshift(lshift(m(5),nbpc),nbpc))
     t1(1,l) = m(5)
     t1(2,l) = m(6)
8150 CONTINUE
     RETURN 3
 
 !     ******************************************************************
 
 !     PROCESS CDUM-I CARDS.
 
8200 CONTINUE
 
     !     ==============
     !     ONLY DO THIS FOR FIRST ONE IF I CAN FIGURE OUT HOW
 
     ASSIGN 8210 TO ret
     GO TO 9010
8210 CONTINUE
     !     ==============
 
     IF (mf(1) /= 1 .OR. mf(2) /= 1) GO TO 7
     IF (m(1) <= 0  .OR. m(2) <= 0) GO TO 8
     l1 = ndumg + 2
     DO  l = 3,l1
         IF (mf(l) /= 1) GO TO 7
         IF (m(l)  <= 0) GO TO 8
         IF (l == 3) CYCLE
         l3 = l - 1
         DO  l2 = 3,l3
             IF (m(l2)-m(l) == 0) THEN
                 GO TO     8
             ELSE
                 GO TO  8215
             END IF
8215     CONTINUE
         END DO
     END DO
     n = ndumc
     GO TO 3
 
 !     ******************************************************************
 
 !     PROCESS PDUM-I CARDS.
 
8300 CONTINUE
 
     !     ==============
     !     ONLY DO THIS FOR FIRST ONE IF I CAN FIGURE OUT HOW
 
     ASSIGN 8310 TO ret
     GO TO 9010
8310 CONTINUE
     !     ==============
 
     IF (mf(1) /= 1 .OR. mf(2) /= 1) GO TO 7
     IF (m(1) <= 0  .OR. m(2) <= 0) GO TO 8
     n = ndump
     GO TO 3
 
 !     ******************************************************************
 
 !     DECODE ADUM-I CARD CONTENTS AS PACKED INTO /SYSTEM/
 
9010 CONTINUE
     ndumg = kdumel(idumel)/10000000
     ndumd = kdumel(idumel) - 10000000*ndumg
     ndumc = ndumd/10000
     ndump = (ndumd - ndumc*10000)/10
     ndumd = kdumel(idumel) - (kdumel(idumel)/10)*10
     ndumc = ndumg + ndumc + 2
     ndump = ndump + 2
     IF (ndumc > 24) GO TO 8
     IF (ndump > 24) GO TO 8
     GO TO ret, (8210,8310)

9999 RETURN
 
 END SUBROUTINE ifs5p
