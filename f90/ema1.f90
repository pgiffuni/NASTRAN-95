SUBROUTINE ema1
     
!     EMA1 ASSEMBLES A STRUCTURAL MATRIX FOR THE MODEL FROM EACH OF
!     THE INDIVIDUAL ELEMENT STRUCTURAL MATRICES.
 
!     EMA1   GPECT,KDICT,KELEM,SIL,ECT / KGG / C,N,NOK4/ C,N,WTMASS
 
!     NOK4 .NE. -1 MEANS MULTIPLY BY DAMPING FACTOR (GE)
!     ABS(WTMASS-1.0) .GT. 1.E-6 MEANS MULTIPLY BY WTMASS
 
!     EMA1 USES 2 SCRATCH FILES
 
 LOGICAL :: last
 INTEGER :: system  ,sysbuf ,zi(1)  ,rd     ,rdrew  ,wrt   ,  &
     wrtrew  ,clsrew ,cls    ,gpect  ,sil    ,ect   ,  &
     subnam(2)       ,scr1   ,scr2   ,elem   ,dofg  ,  &
     mcbkgg(7)       ,typin1 ,typou1 ,typin2 ,prec  ,  &
     trlsil(7)       ,even   ,buf(10),buf1   ,buf2  ,  &
     buf3    ,openw  ,openr  ,silnbr ,opcls  ,mcb(7),  &
     tt(3)   ,oldcod ,scalas(32)     ,dof           , BLOCK(20)
 REAL :: zs(1)   ,xs(4)
 DOUBLE PRECISION :: zd      ,xd     ,d
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg /  ufm     ,uwm    ,uim    ,sfm
 COMMON /system/  system(80)
 COMMON /names /  rd      ,rdrew  ,wrt    ,wrtrew ,clsrew ,cls
 COMMON /BLANK /  nok4    ,wtmass
 COMMON /packx /  typin1  ,typou1 ,ii1    ,jj1    ,incr1
 COMMON /unpakx/  typin2  ,ii2    ,jj2    ,incr2
 COMMON /zblpkx/  xd(2)   ,ix
 COMMON /gpta1 /  nelem   ,jlast  ,incre  ,elem(1)
 COMMON /ma1xx /  d(18)
 COMMON /zzzzzz/  zd(1)
 EQUIVALENCE      (system(1),sysbuf), (system(2),nout),  &
     (trlsil(2),nbrsil), (trlsil(3),luset),  &
     (system(22),mach ), (zd(1),zs(1),zi(1)), (xd(1),xs(1))
 
!     DEFINITION OF INPUT DATA BLOCKS
 
 DATA    gpect ,  kdict, kelem, sil, ect  / 101   ,  102  , 103  , 104, 105  /
 
!     DEFINITION OF OUTPUT DATA BLOCKS
 
 DATA    kgg   /  201  /
 
!     DEFINITION OF SCRATCH FILES
 
 DATA    scr1  ,  scr2 / 301, 302 /
 
!     MISCELANEOUS DATA
 
 DATA    subnam/  4HEMA1,4H    /, nhema1/ 4HEMA1/,  &
     large /  2147483647   /, lpcb  / 8     /
 
!     DATA    TERMS /  1, 0, 9, 0, 0, 18 /,
!    1        SCL   /  1, 1, 0           /
 
!     STATEMENT FUNCTION
 
 even(n) = 2*((n+1)/2)
 
!     PERFORM GENERALIZATION
 
 lcore = korsz(zd)
 trlsil(1) = sil
 CALL rdtrl (trlsil)
 WRITE  (nout,999) (trlsil(i),i=1,7)
 999 FORMAT (1H ,7I10)
 isil0 = lcore - nbrsil - 1
 lcore = isil0
 buf1  = lcore - sysbuf
 buf2  = buf1 - sysbuf
 buf3  = buf2 - sysbuf
 buf(1)= kelem
 CALL rdtrl (buf)
 WRITE (nout,999) (buf(i),i=1,7)
 prec  = buf(2)
 CALL makmcb (mcbkgg,kgg,luset,6,prec)
 openw = wrtrew
 openr = rdrew
 last  = .false.
 silnbr= 0
 opcls = cls
 maxdct= 0
 maxvec= 0
 oldcod= 0
 
!     SET SWITCH FOR MULTIPLICATION BY DAMPING AND/OR WEIGHT MASS FACTOR
 
 eps = ABS(wtmass-1.0)
 IF (eps < 1.e-6 .AND. nok4 < 0) ASSIGN 244 TO kfact
 IF (eps < 1.e-6 .AND. nok4 >= 0) ASSIGN 245 TO kfact
 IF (eps >= 1.e-6 .AND. nok4 < 0) ASSIGN 246 TO kfact
 IF (eps >= 1.e-6 .AND. nok4 >= 0) ASSIGN 247 TO kfact
 
!     READ THE CONTENTS OF THE SIL DATA BLOCK INTO CORE
 
 CALL gopen (sil,zi(buf1),rdrew)
 CALL fread (sil,zi(isil0+1),nbrsil,1)
 CALL CLOSE (sil,clsrew)
 zi(isil0+nbrsil+1) = luset + 1
 CALL cdcbug (nhema1,100,zi(isil0+1),nbrsil+1)
 
!     READ THE KDICT AND ECT DATA BLOCKS. WRITE A MODIFIED KDICT ON SCR2
!     WHICH INCLUDES THE INTERNAL GRID NUMBERS FOR EACH ELEMENT.
!     THE FORMAT FOR EACH RECORD ON SCR2 IS...
!     3-WORD RECORD HEADER
!        1  ELEMENT TYPE
!        2  NBR OF WORDS PER ENTRY( N )
!        3  NBR OF GRID POINTS PER ENTRY
!     N-WORD ELEMENT ENTRY
!        1  ELEMENT ID( INTERNAL NUMBER )
!        2  FORM OF COLUMN PARTITIONS( 1=RECT, 2=DIAG )
!        3  NUMBER OF TERMS PER COLUMN PARTITION
!        4  SCALAR CODE DEFINING DOF PER GRID POINT
!        5  GE
!        6  INTERNAL INDEX OF 1ST GRID POINT
!        7  GINO ADDRESS OF 1ST COLUMN PARTITION
!       ...
!       N-1 INTERNAL INDEX OF LAST GRID POINT
!        N  GINO ADDRESS OF LAST COLUMN PARTITION
 
!     NOTE...
!     GRID POINTS ARE IN SORT BY INTERNAL INDEX. ZERO INDICATES
!     MISSING GRID POINT. ANY ZERO-S ARE LAST IN LIST.
 
 CALL gopen (kdict,zi(buf1),rdrew )
 CALL gopen (ect  ,zi(buf2),rdrew )
 CALL gopen (scr2 ,zi(buf3),wrtrew)
 111 CALL READ  (*124,*111,kdict,buf(4),3,0,j)
 CALL cdcbug (nhema1,111,buf(4),3)
 112 CALL ectloc (*900,ect,buf,i)
 CALL cdcbug (nhema1,112,buf,3)
 IF (elem(i+2) == buf(4)) GO TO 114
 CALL skprec (ect,1)
 GO TO 112
 114 buf(5) = buf(5) + buf(6)
 CALL WRITE (scr2,buf(4),3,0)
 igrid  = elem(i+12)
 nbrgrd = elem(i+ 9)
 nwdect = elem(i+ 5)
 idict  = nwdect + 1
 nwddct = buf(5) - buf(6)
 ngrid  = igrid + nbrgrd - 1
 maxdct = MAX0(maxdct,buf(5))
 IF (nbrgrd /= buf(6)) GO TO 901
 115 CALL READ (*122,*122,ect,zi,nwdect,0,j)
 CALL cdcbug (nhema1,115,zi,nwdect)
 CALL fread (kdict,zi(idict),nwddct,0)
 CALL cdcbug (nhema1,116,zi(idict),nwddct)
 DO  j = igrid,ngrid
   IF (zi(j) == 0) zi(j) = large
 END DO
 CALL sort (0,0,1,1,zi(igrid),nbrgrd)
 DO  j = igrid,ngrid
   IF (zi(j) == large) zi(j) = 0
 END DO
 CALL cdcbug (nhema1,118,zi(igrid),nbrgrd)
 CALL WRITE  (scr2,zi(idict),nwddct-nbrgrd,0)
 iloc = idict + nwddct - nbrgrd
 DO  j = 1,nbrgrd
   CALL WRITE (scr2,zi(igrid+j-1),1,0)
   CALL WRITE (scr2,zi(iloc +j-1),1,0)
 END DO
 maxvec = MAX0(maxvec,zi(idict+2)*prec)
 GO TO 115
 122 CALL skprec (kdict,1)
 CALL WRITE (scr2,0,0,1)
 GO TO 111
 124 CALL CLOSE (kdict,clsrew)
 CALL CLOSE (  ect,clsrew)
 CALL CLOSE ( scr2,clsrew)
 tt(1) = maxdct
 tt(2) = maxvec
 CALL cdcbug (nhema1,125,tt,2)
 
!     READ GPECT AND PREPARE THE SCR1 DATA BLOCK. FOR EACH GRID/SCALAR
!     POINT, TWO RECORDS ARE WRITTEN. THE 1ST CONTAINS 6 WORDS...
!       1  INTERNAL INDEX OF GRID/SCALAR POINT
!       2  DOF OF POINT (1=SCALAR, 6=GRID)
!       3  DOF OF EACH CONNECTED POINT (0 IF NO CONNECTED POINTS)
!       4  NUMBER OF CONNECTED POINTS
!       5  INDEX OF  1ST CONNECTED POINT
!       6  INDEX OF LAST CONNECTED POINT
 
!     THE 2ND RECORD IS A PACKED COLUMN WHICH CONTAINS A NON-ZERO TERM
!     FOR EACH CONNECTED POINT.
 
 typin1 = 1
 typou1 = 1
 incr1  = 1
 incr2  = 1
 CALL makmcb (mcb,scr1,nbrsil,1,1)
 ilook0 = nbrsil + 1
 IF (ilook0+luset+1 >= buf3) CALL mesage (-8,0,subnam)
 DO  i = 1,nbrsil
   j = zi(isil0+i)
   zi(ilook0+j) = i
 END DO
 CALL cdcbug (nhema1,131,zi(ilook0+1),luset)
 CALL gopen (gpect,zi(buf1),rdrew )
 CALL gopen (scr1 ,zi(buf2),wrtrew)
 DO  ii = 1,nbrsil
   nbrcon = buf(4)
   minnbr = buf(5)
   maxnbr = buf(6)
   IF (ii /= 1) GO TO 130
   nbrcon = nbrsil
   minnbr = 1
   maxnbr = nbrsil
   130 CALL fread (gpect,buf,2,0)
   buf(1) = ii
   buf(3) = 0
   buf(4) = 0
   buf(5) = large
   buf(6) = 0
   IF (nbrcon == 0) GO TO 134
   DO  i = minnbr,maxnbr
     zi(i) = 0
   END DO
   134 CALL READ (*138,*138,gpect,tt,3,0,i)
   CALL cdcbug (nhema1,134,tt,3)
   nbrgrd = IABS(tt(1)) - 2
   DO  i = 1,nbrgrd
     CALL fread (gpect,silnbr,1,0)
     j = zi(ilook0+silnbr )
     IF (zs(j) /= 0) CYCLE
     buf(3) = MAX0(buf(3),zi(isil0+j+1)-zi(isil0+j))
     buf(4) = buf(4) + 1
     buf(5) = MIN0(buf(5),j)
     buf(6) = MAX0(buf(6),j)
     zs(j)  = 1.0
   END DO
   GO TO 134
   138 CALL WRITE (scr1,buf,6,1)
   CALL cdcbug (nhema1,138,buf,6)
   IF (buf(4) == 0) GO TO 142
   
!     PACK COLUMN FOR POINT WITH CONNECTED POINTS
   
   ii1 = buf(5)
   jj1 = buf(6)
   CALL cdcbug (nhema1,139,zi(ii1),jj1-ii1+1)
   CALL pack (zs(ii1),scr1,mcb)
   CYCLE
   
!     HERE IF PIVOT HAS NO CONNECTED POINTS
   
   142 CONTINUE
   mcb(2) = mcb(2) + 1
   
!     CLOSE FILES
   
 END DO
 CALL CLOSE (gpect,clsrew)
 CALL CLOSE ( scr1,clsrew)
 CALL wrttrl (mcb)
 
!     ALLOCATE STORAGE FOR MAXIMUM COLUMN OF ELEMENT MATRIX
!     AND MAXIMUM ENTRY FROM MODIFIED KDICT( SCR2 )
 
 idict = maxvec + 1
 igrid = idict + 5
 ipvt  = idict + maxdct
 lcore = even( buf2 ) - 1
 
 
!     BEGIN A PASS BY OPENING SCR1 AND SETTING ALLOCATION POINTERS
 
 
 150 CALL gopen (scr1,zi(buf1),openr)
 ii = ipvt
 jj = lcore
 
!     BEGIN A PIVOT ALLOCATION BY READING PIVOT CONTROL BLOCK FROM SCR1
 
 160 CONTINUE
 tt(1) = ii
 tt(2) = jj
 CALL cdcbug (nhema1,160,tt,2)
 IF (ii+lpcb >= jj) GO TO 202
 CALL fread (scr1,zi(ii),6,1)
 silnbr = zi(ii)
 zi(ii+6) = 0
 zi(ii+7) = 0
 IF (zi(ii+3) == 0) GO TO 195
 
!     ATTEMPT TO ALLOCATE SPACE FOR CONNECTED GRID VECTOR
!     AND FOR MATRICES CONNECTED TO THE PIVOT
 
 nwdcgv = zi(ii+5) - zi(ii+4) + 1
 nwdmat = prec*zi(ii+1)*zi(ii+2)*zi(ii+3)
 IF (ii+lpcb >= jj-nwdcgv-nwdmat) GO TO 200
 imat = jj - nwdmat
 zi(ii+6) = imat - even(nwdcgv)
 zi(ii+7) = imat
 jj   = zi(ii+6)
 nmat = imat + nwdmat - 1
 DO  i = imat,nmat
   zs(i) = 0
 END DO
 icgvec = jj
 ncgvec = icgvec + nwdcgv - 1
 
!     UNPACK CONNECTED GRID VECTOR. CONVERT NON-ZERO POSITIONS TO
!     RELATIVE POINTERS (IN PRECISION OF PROBLEM) TO THE CORRESPONDING
!     1ST TERM OF THE ELEMENT MATRIX
 
 ii2 = zi(ii+4)
 jj2 = zi(ii+5)
 ntrmec = zi(ii+2)
 kk = 1
 typin2 = 1
 CALL unpack (*902,scr1,zs(icgvec))
 DO  i = icgvec,ncgvec
   IF (zi(i) == 0) CYCLE
   zi(i) = kk
   kk = kk + ntrmec
 END DO
 CALL cdcbug (nhema1,174,zi(ii),8)
 CALL cdcbug (nhema1,175,zi(icgvec),nwdcgv)
 IF (kk-1 /= zi(ii+2)*zi(ii+3)) GO TO 903
 
!     TEST FOR LAST PIVOT. IF NOT, TRY TO ALLOCATE ANOTHER PIVOT
 
 195 IF (silnbr == nbrsil) GO TO 210
 ii = ii + lpcb
 GO TO 160
 
!     HERE IF CURRENT PIVOT CANNOT BE ALLOCATED -- MAKE SURE AT LEAST
!     ONE PIVOT HAS BEEN ALLOCATED.
 
 200 CALL bckrec (scr1)
 202 IF (ii == ipvt) CALL mesage (-8,0,subnam)
 npvt = ii - lpcb
 GO TO 220
 
!     HERE WHEN LAST PIVOT HAS BEEN READ AND ALLOCATED
 
 210 last = .true.
 opcls= clsrew
 npvt = ii
 
 
!     CLOSE SCR1, OPEN SCR2 AND KELEM. PREPARE TO ASSEMBLE
!     STRUCTURAL MATRIX FOR THOSE PIVOTS CURRENTLY ALLOCATED.
 
 
 220 CONTINUE
 CALL CLOSE (scr1,opcls)
 CALL gopen (scr2, zi(buf1),rdrew)
 CALL gopen (kelem,zi(buf2),rdrew)
 
!     READ HEADER FOR CURRENT ELEMENT TYPE FROM SCR2
 
 230 CONTINUE
 CALL READ (*260,*230,scr2,tt,3,0,i)
 CALL cdcbug (nhema1,230,tt,3)
 nwddct = tt(2)
 ngrid  = igrid + 2*(tt(3)-1)
 
!     READ AN ELEMENT DEFINITION. IF ANY GRID POINT IS IN CURRENT
!     ALLOCATION, PREPARE TO PROCESS IT.
 
 240 CALL READ (*230,*230,scr2,zi(idict),nwddct,0,i)
 CALL cdcbug (nhema1,240,zi(idict),nwddct)
 DO  i = igrid,ngrid,2
   IF (zi(i) >= zi(ipvt) .AND. zi(i) <= zi(npvt))  &
       GO TO kfact, (244,245,246,247)
 END DO
 GO TO 240
 244 factor = 1.0
 GO TO 248
 245 factor = zs(idict+4)
 GO TO 248
 246 factor = wtmass
 GO TO 248
 247 factor = wtmass*zs(idict+4)
 
!     DECODE RELATIVE COLUMN NUMBERS
 
 248 IF (oldcod == zi(idict+3)) GO TO 250
 icode = zi(idict+3)
 CALL decode (icode,scalas,nsca)
 oldcod = zi(idict+3)
 
!     READ EACH COLUMN OF THE ELEMENT MATRIX.
!     ADD IT TO THE STRUCTURAL MATRIX.
 
 250 nwdcol = prec*zi(idict+2)
 IF (zi(idict+1) == 2) nwdcol = prec
 252 ii = ipvt + (zi(i)-zi(ipvt))*lpcb
 tt(1) = i
 tt(2) = zi(i)
 tt(3) = nsca
 CALL cdcbug (nhema1,252,tt,3)
 CALL filpos (kelem,zi(i+1))
 icgvec = zi(ii+6)
 imat   = zi(ii+7)
 DO  j = 1,nsca
   CALL fread (kelem,zi,nwdcol,0)
   CALL cdcbug(nhema1,254,zi,nwdcol)
   IF (prec == 1) CALL ema1s (j,nsca,scalas,zi(ii),zi(idict),  &
       zi(icgvec),zi(imat),zi,factor)
   IF (prec == 2) CALL ema1d (j,nsca,scalas,zi(ii),zi(idict),  &
       zi(icgvec),zi(imat),zi,factor)
 END DO
 255 IF (i == ngrid) GO TO 240
 i = i + 2
 IF (zi(i) >= zi(ipvt) .AND. zi(i) <= zi(npvt)) GO TO 252
 GO TO 255
 
!     ALL COLUMNS OF STRUCTURAL MATRIX NOW ALLOCATED ARE COMPLETE.
!     OPEN KGG AND PACK COLUMNS.
 
 260 CALL CLOSE (scr2,clsrew)
 CALL CLOSE (kelem,clsrew)
 CALL gopen (kgg,zi(buf1),openw)
 DO  ii = ipvt,npvt,lpcb
   dof    = zi(ii+1)
   dofg   = zi(ii+2)
   nbrcon = zi(ii+3)
   icgvec = zi(ii+6)
   imat   = zi(ii+7)
   ii1    = zi(ii+4)
   ii2    = zi(ii+5)
   kk     = imat
   CALL cdcbug (nhema1,260,zi(imat),((ii2-ii1+1)*(dof*dofg)))
   
!     PACK COLUMNS WITH BLDPK
   
   DO  jj = 1,dof
     CALL bldpk (prec,prec,kgg,BLOCK,1)
     IF (nbrcon == 0) GO TO 266
     i = icgvec
     DO  j = ii1,ii2
       IF (zi(i) == 0) GO TO 263
       k  = zi(isil0+j)
       n  = k + MIN0(dofg,zi(isil0+j+1)-zi(isil0+j)) - 1
       ll = kk
       DO  silnbr = k,n
         CALL bldpki (zs(ll),silnbr,kgg,BLOCK)
         ll = ll + prec
       END DO
       kk = kk + dofg*prec
       263 i  = i + 1
     END DO
     266 CALL bldpkn (kgg,BLOCK,mcbkgg)
   END DO
 END DO
 CALL CLOSE (kgg,opcls)
 
!     TEST FOR COMPLETION OF LAST PASS
 
 IF (last) GO TO 310
 openr = rd
 openw = wrt
 GO TO 150
 
!     KGG NOW COMPLETE -- WRITE ITS TRAILER.
 
 310 CONTINUE
 CALL wrttrl (mcbkgg)
 RETURN
 
!     FATAL ERRORS
 
 900 kerr = 112
 GO TO 990
 901 kerr = 114
 GO TO 990
 902 kerr = 172
 GO TO 990
 903 kerr = 174
 GO TO 990
 
!     PROCESS LOGIC ERROR
 
 990 WRITE  (nout,991) sfm,kerr
 991 FORMAT (a25,' 3102, EMA1 LOGIC ERROR',i4)
 IF (mach == 2 .OR. mach == 5 .OR. mach == 21) kerr = -kerr
 CALL gperr (subnam,kerr)
 RETURN
END SUBROUTINE ema1
