SUBROUTINE cyct1
     
!     GENERATE CYCLIC TRANSFORMATION MATRIX, TRANSFORM VECTORS
 
!     DMAP CALLING SEQUENCE
 
!     CYCT1   VIN/VOUT,GCYC/V,Y,CTYPE/V,Y,CDIR/V,Y,N/V,Y,KMAX/
!             V,Y,NLOAD/V,N,NOGO $
 
 LOGICAL :: lback, lcos, ldrl, ldsa, lnmult, lvin, lvout
 INTEGER :: buf, cdir, ctype, gcyc, hback, hdrl, hdsa, hrot,  &
     iz, mcb(7), pkin, pkincr, pkirow, pknrow, pkout,  &
     precis, scrt, subr(2), sysbuf, vin, vout, outpt
 REAL :: rz(1)
 DOUBLE PRECISION :: dc, dc1, dfac, dfak, ds, ds1, dz(1)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg /  ufm
 COMMON /system/  ksystm(65)
 COMMON /packx /  pkin,pkout,pkirow,pknrow,pkincr
 COMMON /zzzzzz/  iz(1)
 COMMON /BLANK /  ctype(2),cdir(2),nn,kmaxi,nload,nogo
 EQUIVALENCE      (ksystm( 1),sysbuf), (ksystm( 2) ,outpt),  &
     (ksystm(55),iprec ), (iz(1),rz(1),dz(1))
 DATA    subr  /  4HCYCT, 4H1   /,   hback  /  4HBACK    /
 DATA    vin   /  101  /, vout,gcyc /201,202/, scrt  /301/
 DATA    hrot  ,  hdrl  , hdsa  /4HROT ,4HDRL ,4HDSA     /
 
 
!     FIND NECESSARY PARAMETERS
 
 nogo   = 1
 precis = 2
 IF (iprec /= 2) precis = 1
 ldrl = ctype(1) == hdrl
 ldsa = ctype(1) == hdsa
 IF (.NOT.((ctype(1) == hrot) .OR. ldrl.OR.ldsa)) GO TO 310
 10 lback = cdir(1) == hback
 
!     CURRENT DOCUMENTED USAGE DOES NOT USE NEGATIVE VALUES OF KMAXI
!     OTHER THAN THE DEFAULT OF -1     10/02/73
!     LOGIC IS INCLUDED IN THE ROUTINE TO USE NEGATIVE KMAXI BUT IS NOT
!     FULLY CHECKED OUT.  THE FOLLOWING STATEMENT NEGATES ALL THIS LOGIC
 
 IF (kmaxi < 0) kmaxi = nn/2
 kmax = kmaxi
 kmin = 0
 IF (kmax >= 0) GO TO 20
 kmax =-kmax
 kmin = kmax
 20 IF (2*kmax > nn .OR. nn <= 0) GO TO 330
 30 IF (nload <= 0) GO TO 350
 40 nloads = nload
 IF (ldsa) nloads = 2*nload
 nloadt = nload
 IF (ldrl .OR. ldsa) nloadt = 2*nload
 numrow = nn
 IF (.NOT.lback) GO TO 50
 numrow = 2*(kmax-kmin+1)
 IF (kmin   ==  0) numrow = numrow - 1
 IF (2*kmax == nn) numrow = numrow - 1
 50 numrow = nloadt*numrow
 
!     DEFINE OPEN CORE POINTERS AND GINO BUFFER
!                POINTERS                BEGIN    END
!          TABLE OF COS (2.0*PI*N/NN)     ICOS    NCOS
!          TABLE OF SIN (2.0*PI*N/NN)     ISIN    NSIN
!          AREA TO ASSEMBLE COLUMNS       ICOL    NCOL
!          (NOTE  N = LITTLE N, NN = CAPITAL N)   N = 0,(NN-1)
!          (ALLOW FOR (NLOADS-1) ZEROS BEFORE FIRST ENTRY IN COL)
 
 buf  = korsz(iz) - sysbuf + 1
 icos = 1
 ncos = icos + nn - 1
 isin = ncos + 1
 nsin = isin + nn - 1
 icol = nsin + 1
 jcol = icol + nloads - 1
 ncol = jcol + numrow - 1
 IF (2*ncol >= buf) CALL mesage (-8,0,subr)
 
!     CHECK DATA BLOCK TRAILERS
 
 mcb(1) = gcyc
 CALL rdtrl (mcb(1))
 IF (mcb(1) <= 0) GO TO 370
 60 mcb(1) = vout
 CALL rdtrl (mcb(1))
 lvout  = mcb(1) > 0
 mcb(1) = vin
 CALL rdtrl (mcb(1))
 lvin = mcb(1) > 0
 IF (.NOT.lvin) mcb(2) = 0
 lnmult = mcb(2) /= numrow
 IF (lvin .AND. lvout .AND. lnmult) GO TO 390
 IF (nogo > 0) THEN
   GO TO    70
 ELSE
   GO TO   410
 END IF
 
!     THE PARAMETERS ARE OK
!     PREPARE TRIGONOMETRIC TABLES,  DC1=COS(2*PI/NN), PI = 4*ATAN(1)
!     MOVABLE POINTERS  JXXX=N , KXXX= NN-N
 
 70 rn   = FLOAT(nn)
 dfac = (8.0D0*DATAN(1.0D0))/DBLE(rn)
 dc1  = DCOS(dfac)
 ds1  = DSIN(dfac)
 jcos = icos
 kcos = ncos + 1
 jsin = isin
 ksin = nsin + 1
 dz(jcos) = 1.0D0
 dz(jsin) = 0.0D0
 80 IF (kcos-jcos-2 < 0) THEN
   GO TO   120
 ELSE IF (kcos-jcos-2 == 0) THEN
   GO TO    90
 ELSE
   GO TO   100
 END IF
 90 dc   =-1.0D0
 ds   = 0.0D0
 GO TO 110
 100 dc   = dc1*dz(jcos) - ds1*dz(jsin)
 ds   = ds1*dz(jcos) + dc1*dz(jsin)
 110 jcos = jcos + 1
 jsin = jsin + 1
 kcos = kcos - 1
 ksin = ksin - 1
 dz(jcos) = dc
 dz(jsin) = ds
 dz(kcos) = dc
 dz(ksin) =-ds
 GO TO 80
 
!     ZERO THE AREA FOR FORMING THE COLUMN
 
 120 DO  j = icol,ncol
   dz(j) = 0.0D0
 END DO
 
!     OPEN GCYC MATRIX,  GET READY TO USE PACK
 
 CALL gopen (gcyc,iz(buf),1)
 CALL makmcb (mcb,gcyc,numrow,2,precis)
 pkin   = 2
 pkout  = precis
 pkirow = 1
 pknrow = numrow
 pkincr = 1
 IF (lback) GO TO 240
 
!     START LOOPING ON COLUMNS OF MATRIX OF TYPE FORE.
!     FORM A COLUMN AND PACK IT OUT
!          K = KMIN,KMAX  ALTERNATE COSINE AND SINE COLUMNS
 
 dfac = 2.0D0/DBLE(rn)
 IF (ldrl) dfac = 0.5D0*dfac
 k    = kmin
 140 dfak = dfac
 IF (k == 0 .OR. 2*k == nn) dfak = 0.5D0*dfak
 lcos  = .true.
 ktrig = icos
 ntrig = ncos
 GO TO 160
 150 lcos  = .false.
 ktrig = isin
 ntrig = nsin
 160 DO  kcol = jcol,ncol,nloadt
   dz(kcol) = dfak*dz(ktrig)
   ktrig = ktrig + k
   IF (ktrig > ntrig) ktrig = ktrig - nn
 END DO
 
!     PACK OUT NLOADT COLUMNS  (FOR EITHER FORE OR BACK)
!      IF  ROT OR DSA   WE ARE READY
!      IF     DRL       PRODUCE INTERMEDIATE TERMS FIRST (EXPAND)
 
 180 nxcol = 1
 IF (.NOT.ldrl) GO TO 220
 DO  kcol = jcol,ncol,nloadt
   kcol2 = kcol + nloads
   dz(kcol2) = dz(kcol)
 END DO
 GO TO 220
 200 nxcol = 2
 DO  kcol = jcol,ncol,nloadt
   kcol2 = kcol + nloads
   dz(kcol2) = -dz(kcol)
 END DO
 220 kcol = jcol
 230 CALL pack (dz(kcol),gcyc,mcb)
 kcol = kcol - 1
 IF (kcol >= icol) GO TO 230
 IF (ldrl .AND. nxcol == 1) GO TO 200
 IF (lback) GO TO 280
 
!     BOTTOM OF LOOP FOR TYPE FORE
 
 IF (k /= 0 .AND. 2*k /= nn .AND. lcos) GO TO 150
 k = k + 1
 IF (k-kmax > 0) THEN
   GO TO   290
 ELSE
   GO TO   140
 END IF
 
!     START LOOPING ON COLUMNS OF MATRIX OF TYPE BACK
!         N = 1,NN
 
 240 n = 1
 250 k = 0
 kcos = icos
 kcol = jcol
 260 IF (k < kmin) GO TO 270
 dz(kcol) = dz(kcos)
 kcol = kcol + nloadt
 IF (k == 0 .OR. 2*k == nn) GO TO 270
 dz(kcol) = dz(kcos+nn)
 kcol = kcol + nloadt
 270 kcos = kcos + n - 1
 IF (kcos > ncos) kcos = kcos - nn
 k = k + 1
 IF (k-kmax > 0) THEN
   GO TO   180
 ELSE
   GO TO   260
 END IF
 
!     BOTTOM OF LOOP FOR TYPE BACK
 
 280 n = n + 1
 IF (n-nn > 0) THEN
   GO TO   290
 ELSE
   GO TO   250
 END IF
 
!     THE GCYC MATRIX IS NOW COMPLETE
 
 290 CALL CLOSE  (gcyc,1)
 CALL wrttrl (mcb(1))
 
!     IF WE HAVE TO FORM VOUT, USE SSG2B.  (VOUT = VIN*GCYC)
 
 IF (lnmult) GO TO 300
 CALL ssg2b (vin,gcyc,0,vout,0,precis,1,scrt)
 300 CONTINUE
 RETURN
 
!     FATAL MESSAGES
 
 310 nogo = -1
 WRITE  (outpt,320) ufm,ctype(1)
 320 FORMAT (a23,' 4063, ILLEGAL VALUE (',a4,') FOR PARAMETER CTYPE.')
 GO TO 10
 330 nogo = -1
 WRITE  (outpt,340) ufm,nn,kmaxi
 340 FORMAT (a23,' 4064, ILLEGAL VALUES (',i8,1H,,i8,  &
     ') FOR PARAMETERS (NSEGS,KMAX).')
 GO TO 30
 350 nogo = -1
 WRITE  (outpt,360) ufm,nload
 360 FORMAT (a23,' 4065, ILLEGAL VALUE (',i8,') FOR PARAMETER NLOAD.')
 GO TO 40
 370 nogo = -1
 WRITE  (outpt,380) ufm
 380 FORMAT (a23,' 4066, SECOND OUTPUT DATA BLOCK MUST NOT BE PURGED.')
 GO TO 60
 390 nogo = -1
 WRITE  (outpt,400) ufm,mcb(2),numrow
 400 FORMAT (a23,' 4067, VIN HAS',i9,' COLS, GCYC HAS',i9,6H rows.)
 410 CALL mesage (-61,0,subr)
 RETURN
END SUBROUTINE cyct1
