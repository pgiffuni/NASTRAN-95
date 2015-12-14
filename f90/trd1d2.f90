SUBROUTINE trd1d2
     
!     THIS ROUTINE COMPUTES NON-LINEAR LOADS FOR TRANSIENT ANALYSIS
 
!     THIS ROUTINE IS SUITABLE FOR DOUBLE PRECISION OPERATION
 
 LOGICAL :: dec
 INTEGER :: iz(1),pnl,dit,FILE,sysbuf,itlist(13),NAME(2), nmtd(2)
 DIMENSION        z(1)
 DOUBLE PRECISION :: x,y,dz,h,fx,fy
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg /  ufm,uwm
 COMMON /system/  sysbuf,iout
 COMMON /machin/  mach
 COMMON /zzzzzz/  dz(1)
 COMMON /packx /  it1,it2,ii,nrow,incr
 COMMON /trdd1 /  nlft,dit,nlftp,nout,icount,iloop,modal,lcore,  &
     icore,iu,ip,ipnl(7),nmodes,nstep,pnl,ist,iu1, deltat,ifrst,tabs,sigma,tim
 EQUIVALENCE      (z(1),iz(1),dz(1))
 DATA    itlist/  4,1105,11,1,1205,12,2,1305,13,3,1405,14,4/
 DATA    NAME  /  4HNLFT,4HTRDD /
 DATA    nmtd  /  4HTRD1,4HD2   /
 DATA    kount /  0             /
 
!     IDENTIFICATION OF VARIABLES
 
!     NLFT    NON-LINEAR FUNCTION TABLE
!     PNL     NON-LINEAR FORCES --MATRIX
!     DIT     DIRECT INPUT TABLES
!     NLFTP   NON-LINEAR FUNCTION SET SELECTION
!     NOUT    OUT PUT  EVERY NOUT TIME STEPS( PLUS 1 AND NSTEP)
!     ICOUNT  CURRENT INTERATION COUNTER
!     ILOOP   LOOP ON NUMBER OF TIME STEP CHANGES
!     MODAL   LESS THAN ZERO IMPLIES THIS IS A DIRECT FORMULATION
!     LCORE   AMOUNT OF CORE FOR TRD1D
!     ICORE   POINTER TO FIRST CELL OF OPEN CORE
!     IU      POINTER TO LATEST DISPLACEMENT VECTOR
!     IU1     POINTER TO DISPLACEMENT VECTOR -- ONE TIME STEP BACK
!     IP      POINTER TO LOAD VECTOR
!     NMODES  NUMBER OF MODES IN PROBLEM
!     NSTEP   NUMBER OF TIME STEPS
!     ITLIST  LIST OF CARD TYPES FOR DYNAMIC TABLES
!     NROW    SIZE OF SOLUTION SET
!     IBUF1   POINTER TO BUFFER
!     NCARDS  NUMBER OF LOAD CARDS IN SELECTED SET
!     ICARDS  POINTER  TO FIRST CARD
!     NTABL   NUMBER OF TABLES
!     ITABL   POINTER TO FIRST TABLE
!     IPNL    MATRIX CONTROL BLOCK FOR PNL
 
!     DESCRIPTION OF TYPES OF NON-LINEAR LOADING
 
!     TYPE    DESCRIPTION
!     ----    -----------
 
!       1     DISPLACEMENT-DEPENDENT NOLIN1 LOAD
!       2     DISPLACEMENT-DEPENDENT/DISPLACEMENT-DEPENDENT NOLIN2 LOAD
!       3     DISPLACEMENT-DEPENDENT NOLIN3 LOAD
!       4     DISPLACEMENT-DEPENDENT NOLIN4 LOAD
!       5     VELOCITY-DEPENDENT NOLIN1 LOAD
!       6     VELOCITY-DEPENDENT/DISPLACEMENT-DEPENDENT NOLIN2 LOAD
!       7     VELOCITY-DEPENDENT NOLIN3 LOAD
!       8     VELOCITY-DEPENDENT NOLIN4 LOAD
!       9     VELOCITY-DEPENDENT/VELOCITY-DEPENDENT NOLIN2 LOAD
!      10     DISPLACEMENT-DEPENDENT/VELOCITY-DEPENDENT NOLIN2 LOAD
!      11     TEMPERATURE-DEPENDENT CONVECTION NON-LINEAR LOAD (FTUBE)
!      12     TEMPERATURE-DEPENDENT EMISSIVITIES-ABSORPTIVITIES, NOLIN5
!      13     DISPLACEMENT-DEPENDENT/VELOCITY-DEPENDENT NOLIN6 LOAD
!      14     VELOCITY-DEPENDENT/DISPLACEMENT-DEPENDENT NOLIN6 LOAD
 
!     DETERMINE ENTRY NUMBER
 
 dec = mach == 5 .OR. mach == 6 .OR. mach == 21
 ipx = ip
 
 IF ((iloop == 1 .AND. icount > 1) .OR. (iloop > 1 .AND.  &
     icount > 0)) GO TO 170
 IF (ifrst  /= 0) GO TO 170
 
!     FIRST TIME FOR TIME STEP
 
 CALL sswtch (10,ialg)
 ibuf1 = lcore + icore - sysbuf
 FILE  = nlft
 lcore = lcore - sysbuf - 1
 icrq  =-lcore
 IF (lcore <= 0) GO TO 430
 CALL OPEN (*400,nlft,iz(ibuf1),0)
 
!     FIND SELECTED SET ID
 
 CALL READ (*420,*10,nlft,iz(icore+1),lcore,0,iflag)
 icrq = lcore
 GO TO 430
 10 DO  i = 3,iflag
   k = i + icore
   IF (iz(k) == nlftp) GO TO 30
 END DO
 CALL mesage (-31,nlftp,NAME)
 
!     FOUND SET ID -- POSITION TO RECORD IN NLFT
 
 30 k = i-3
 IF (k == 0) GO TO 50
 DO  i = 1,k
   CALL fwdrec (*420,nlft)
 END DO
 
!     BRING IN  8 WORDS PER CARD
!     FORMAT =    TYPE,SILD,SILE,A,SILD,SILE,A OR SILD,SILE
!     CONVERT TO  TYPE,ROWP,ROWP,A,ROWP OR A
!     COUNT NUMBER OF CARDS
 
 50 ncards = 0
 icards = icore + 1
 k      = icards
 60 icrq   = 8 - lcore
 IF (icrq  > 0) GO TO 430
 CALL READ (*420,*80,nlft,iz(k),8,0,iflag)
 IF (modal < 0) GO TO 70
 
!     MODAL FORM -- CONVERT SILE TO ROW POSITIONS AND STORE IN SILD
 
 IF (iz(k+2) == 0) GO TO 440
 iz(k+1) = iz(k+2) + nmodes
 IF (iz(k+5) == 0) GO TO 440
 iz(k+4) = iz(k+5) +  nmodes
 IF (iz(k) /= 2 .AND. iz(k) /= 6 .AND. iz(k) /= 9 .AND. iz(k) /= 10) GO TO 70
 IF (iz(k+7) == 0) GO TO 440
 iz(k+6) = iz(k+7) + nmodes
 70 CONTINUE
 
!     MOVE UP
 
 iz(k+2) = iz(k+4)
 iz(k+4) = iz(k+6)
 k       = k + 5
 lcore   = lcore  - 5
 ncards  = ncards + 1
 GO TO 60
 
!     END OF RECORD-- DONE
 
 80 CALL CLOSE (nlft,1)
 
!     EXTRACT LIST OF  UNIQUE TABLES FROM CARD TYPES 1,5,11 AND 14
 
 l     = icards
 ntabl = 0
 itabl = k
 DO  i = 1,ncards
   izl = iz(l)
   IF (izl /= 1 .AND. izl /= 5 .AND. (izl < 11 .OR. izl > 14)) GO TO 110
   IF (izl /= 11 .AND. izl /= 12) GO TO 85
   izl = iz(l+4)
   IF (iz(l) /= 11) GO TO 83
   
!     NFTUBE CARD
   
   81 nxx = numtyp(izl)
   IF (dec .AND. izl > 16000 .AND. izl <= 99999999) nxx = 1
   IF (nxx-1 == 0) THEN
     GO TO    85
   ELSE
     GO TO   110
   END IF
   
!     NOLIN5 CARD
   
   83 nxx   = numtyp(iz(l+3))
   IF (dec .AND. iz(l+3) > 16000 .AND. iz(l+3) <= 99999999) nxx = 1
   IF (nxx /= 1) GO TO 81
   itid1 = iz(l+3)
   nxx   = numtyp(izl)
   IF (dec .AND. izl > 16000 .AND. izl <= 99999999) nxx = 1
   IF (nxx /= 1) GO TO 87
   itid2 = iz(l+4)
   numtb = 2
   GO TO 89
   85 itid1 = iz(l+4)
   87 numtb = 1
   89 CONTINUE
   
!     FIND OUT IF UNIQUE TABLE
   
   IF (ntabl == 0) GO TO 100
   DO  m = 1,ntabl
     k = itabl + m
     IF (iz(k) == itid1) GO TO 110
   END DO
   
!     NEW TABLE
   
   100 ntabl = ntabl + 1
   k     = itabl + ntabl
   iz(k) = itid1
   110 CONTINUE
   IF (numtb == 1) GO TO 115
   numtb = 1
   itid1 = itid2
   GO TO 89
   115 l     = l + 5
 END DO
 
 iz(itabl) = ntabl
 lcore = lcore - ntabl - 1
 icrq  =-lcore
 IF (lcore <= 0) GO TO 430
 IF (ntabl == 0) GO TO 150
 
!     INITIALIZE TABLES
 
 k     = itabl + ntabl + 1
 CALL pretab (dit,iz(k),iz(k),iz(ibuf1),lcore,l,iz(itabl),itlist)
 lcore = lcore - l
 IF (ialg == 0) GO TO 140
 in1   = (k + l)/2
 in2   = in1 + nrow
 in3   = in2 + nrow
 lcore = lcore - 6*nrow
 icrq  =-lcore
 IF (lcore < 0) GO TO 430
 
!     ZERO LOAD VECTORS
 
 DO  i = 1,nrow
   k     = in1 + i
   dz(k) = 0.0D0
   k     = in2 + i
   dz(k) = 0.0D0
   k     = in3 + i
   dz(k) = 0.0D0
 END DO
 140 CONTINUE
 150 RETURN
 
!     COMPUTE LOADS
 
 170 k   = icards + ncards*5 - 1
 IF (ialg == 0) GO TO 180
 ipx = in1
 DO  i = 1,nrow
   l   = in1 + i
   dz(l) = 0.0D0
 END DO
 
!     LOOP THRU EACH LOAD CARD OR COLLECTION (NOLIN5, NOLIN6)
 
 180 h  = 1.0D0/deltat
 i  = icards
 190 CONTINUE
 fx = 0.0D0
 fy = 1.0D0
 m  = iu + iz(i+2)
 mm = iu + iz(i+4)
 n  = iu1+ iz(i+2)
 nn = iu1+ iz(i+4)
 x  = dz(m)
 y  = (x-dz(n))*h
 l  = iz(i)
!     L  =     1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14
 SELECT CASE ( l )
   CASE (    1)
     GO TO 200
   CASE (    2)
     GO TO 210
   CASE (    3)
     GO TO 220
   CASE (    4)
     GO TO 230
   CASE (    5)
     GO TO 205
   CASE (    6)
     GO TO 213
   CASE (    7)
     GO TO 225
   CASE (    8)
     GO TO 235
   CASE (    9)
     GO TO 215
   CASE (   10)
     GO TO 217
   CASE (   11)
     GO TO 250
   CASE (   12)
     GO TO 260
   CASE (   13)
     GO TO 240
   CASE (   14)
     GO TO 245
 END SELECT
 
!     NOLIN 1
 
 200 xsp = x
 CALL tab (iz(i+4),xsp,fxsp)
 fx = fxsp
 GO TO 290
 205 x  = y
 GO TO 200
 
!     NOLIN 2
 
 210 y  = dz(mm)
 fx = x*y
 GO TO 290
 213 x  = y
 GO TO 210
 215 x  = y
 217 fx = x*(dz(mm) - dz(nn))*h
 GO TO 290
 
!     NOLIN 3
 
 220 IF (x <= 0.0D0) GO TO 290
 fx = x**z(i+4)
 GO TO 290
 225 x  = y
 GO TO 220
 
!     NOLIN 4
 
 230 IF (x >= 0.0D0) GO TO 290
 fx =-DABS(x)**z(i+4)
 GO TO 290
 235 x  = y
 GO TO 230
 
!     NOLIN6
 
 240 x  = y
 fy = x*DABS(x)
 x  = dz(m)
 GO TO 200
 245 y  = dz(mm)
 fy = y*DABS(y)
 GO TO 200
 
!     NFTUBE.  LOOKUP VDOT IF NEEDED
 
 250 fxsp =  z(i+4)
 izl  = iz(i+4)
 nxx  = numtyp(izl)
 IF (dec .AND. izl > 16000 .AND. izl <= 99999999) nxx = 1
 IF (nxx  ==   1) CALL tab (iz(i+4),tim,fxsp)
 IF (fxsp >= 0.0) m = iu + iz(i+1)
 fx   = fxsp*dz(m)
 l    = ipx + iz(i+2)
 dz(l)= dz(l) + fx*z(i+3)
 fy   =-1.0D0
 GO TO 290
 
!     NOLIN5
 
!     A. COMPUTE SURFACE AVERAGE TEMPERATURES
 
 260 mm = 0
 nn = 0
 tavga = 0.0
 tavgb = 0.0
 j  = 1
 DO  l = 1,4
   IF (l == 3) j = 6
   m  = iz(i+j)
   IF (m == 0) GO TO 265
   m  = iu + m
   tavga = tavga + z(m)
   mm = mm + 1
   265 m  = iz(i+j+10)
   IF (m == 0) GO TO 270
   m  = iu + m
   tavgb = tavgb + z(m)
   nn = nn + 1
   270 j  = j + 1
 END DO
 tavga = tavga/FLOAT(mm)
 tavgb = tavgb/FLOAT(nn)
 aa    = z(i+3)
 ab    = z(i+4)
 fab   = z(i+8)
 fabsq = fab*fab
 eta   = z(i+13)
 etb   = z(i+14)
 nxx   = numtyp(iz(i+13))
 IF (dec .AND. iz(i+13) > 16000 .AND. iz(i+13) <= 99999999) nxx = 1
 IF (nxx == 1) CALL tab (iz(i+13),tavga,eta)
 nxx   = numtyp(iz(i+14))
 IF (dec .AND. iz(i+14) > 16000 .AND. iz(i+14) <= 99999999) nxx = 1
 IF (nxx == 1) CALL tab (iz(i+14),tavgb,etb)
 alpha = z(i+18)
 alphb = z(i+19)
 nxx   = numtyp(iz(i+18))
 IF (dec .AND. iz(i+18) > 16000 .AND. iz(i+18) <= 99999999) nxx = 1
 IF (nxx == 1) CALL tab (iz(i+18),tavga,alpha)
 nxx   = numtyp(iz(i+19))
 IF (dec .AND. iz(i+19) > 16000 .AND. iz(i+19) <= 99999999) nxx = 1
 IF (nxx == 1) CALL tab (iz(i+19),tavgb,alphb)
 alpha = alpha - 1.0
 alphb = alphb - 1.0
 
!     B. COMPUTE DENOMINATOR
 
 xh  = sigma*eta*(tavga+tabs)**4
 xk  = sigma*etb*(tavgb+tabs)**4
 fxsp= alpha*fab*xk - aa*xh + fab*xk - (alphb*fabsq*xh)/ab
 fysp= alphb*fab*xh - ab*xk + fab*xh - (alpha*fabsq*xk)/aa
 fab = 1.0 - (alpha*alphb/aa)*(fabsq/ab)
 fx  = fxsp/(fab*FLOAT(mm))
 fy  = fysp/(fab*FLOAT(nn))
 
!     C. APPLY FORCES ON AREAS A AND  B
 
 j = 1
 DO  l = 1,4
   IF (l == 3) j = 6
   m = iz(i+j)
   IF (m == 0) GO TO 275
   m = ipx + m
   dz(m) = dz(m) + fx
   275 m = iz(i+j+10)
   IF (m == 0) GO TO 280
   m = ipx + m
   dz(m) = dz(m) + fy
   280 j = j + 1
 END DO
 i = i + 20
 GO TO 320
 
!     FINISH APPLYING SCALE FACTOR AND ADD
 
 290 l     = ipx + iz(i+1)
 dz(l) = dz(l) + fx*fy*z(i+3)
 IF (DABS(dz(l)) < 1.0D-36) dz(l) = 0.0D0
 IF (DABS(dz(l)) < 1.0D+36) GO TO 310
 kount = kount + 1
 IF (kount == 1 .OR. kount == 4) WRITE (iout,295)
 IF (kount <= 3) WRITE (iout,300) uwm,dz(l)
 295 FORMAT (/1X,28(4H****),/)
 300 FORMAT (a25,' 3309, UNUSUALLY LARGE VALUE COMPUTED FOR NONLINEAR',  &
     ' FORCING FUNCTION',5X,d15.5)
 310 i = i + 5
 320 IF (i < k) GO TO 190
 
!     END OF LOAD LOOP
 
!     DONE
 
 IF (ialg == 0) GO TO 380
 DO  i = 1,nrow
   
!     SUM OVER LAST THREE LOADS
   
   l  = ip  + i
   k  = in1 + i
   m  = in2 + i
   kk = in3 + i
   dz(l) = dz(l) + (dz(k)+dz(m)+dz(kk))/3.0D0
 END DO
 
!     SWITCH POINTERS
 
 k   = in1
 in1 = in2
 in2 = in3
 in3 = k
 380 RETURN
 
!     ERROR MESSAGES
 
 400 WRITE  (iout,405) ufm
 405 FORMAT (a23,', NON-LINEAR FORCING LOAD (NLFT) WAS NOT GENERATED',  &
     ' PREVIOUSLY')
 ip1 =-37
 410 CALL mesage (ip1,FILE,nmtd)
 RETURN
 420 ip1 =-2
 GO TO 410
 430 ip1  =-8
 FILE = icrq
 GO TO 410
 
!     LOADED POINT  NOT E-POINT IN MODAL FORMULATION
 
 440 CALL mesage (-44,nlftp,iz(k))
 RETURN
END SUBROUTINE trd1d2
