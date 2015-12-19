SUBROUTINE merge (irp,icp,core)
     
!     MERGE WILL PUT UP TO 4 MATRICES, IA11,IA21,IA12,IA22, TOGETHER
!     INTO NAMEA -- THIS ROUTINE IS THE INVERSE OF PARTN
 
!     THE ARGUMENTS ARE EXACTLY THE SAME IN MEANING AND OPTION AS FOR
!     PARTITION
 
 
 INTEGER, INTENT(IN)                      :: irp(1)
 INTEGER, INTENT(IN)                      :: icp(1)
 INTEGER, INTENT(OUT)                     :: core(1)
 EXTERNAL        rshift,andf
 DIMENSION  a11(4),b11(4),  BLOCK(40),NAME(2)
 COMMON /parmeg/ namea,ncola,nrowa,iforma,itypa,ia(2), ia11(7,4),lcare,rule
 COMMON /system/ sysbuf,nout
 COMMON /two   / two1(32)
 COMMON /zblpkx/ ic11(4),ii
 DATA    NAME  / 4HMERG,4HE    /
 
!     CHECK FILES
 
 lcore = IABS(lcare)
 k = namea
 DO  i = 1,4
   IF (k == 0) GO TO 15
   DO  j = i,4
     IF (ia11(1,j) == k) GO TO 440
   END DO
   15 k = ia11(1,i)
 END DO
 
!     PICK UP PARAMETERS AND INITIALIZE
 
 irew  = 0
 IF (lcare < 0) irew = 2
 ncola1= ncola
 ncola = 0
 ia(1) = 0
 ia(2) = 0
 istor = 0
 iotp  = itypa
 nmat  = 0
 DO  i = 1,4
   IF (ia11(1,i) <= 0) CYCLE
!WKBD 2/94 SPR93025 IF (IA11(5,I) .NE. ITYPA) IOTP = 4
   nmat = nmat + 1
   DO  j = 2,5
     IF (ia11(j,i) == 0) GO TO 460
   END DO
 END DO
 ntypa = iotp
 IF (ntypa == 3) ntypa = 2
 ibuf   = lcore - sysbuf + 1
 ibufcp = ibuf - nrowa
 IF (ibufcp > 0) THEN
   GO TO    40
 ELSE
   GO TO   420
 END IF
 40 lcore = ibufcp - 1
 CALL ruler (rule,icp,zcpct,ocpct,core(ibufcp),nrowa,core(ibuf),1)
 IF (irp(1) == icp(1) .AND. irp(1) /= 0) GO TO 60
 ibufrp = ibufcp - (ncola1+31)/32
 IF (ibufrp > 0) THEN
   GO TO    50
 ELSE
   GO TO   420
 END IF
 50 CALL ruler (rule,irp,zrpct,orpct,core(ibufrp),ncola1,core(ibuf),0)
 lcore = ibufrp - 1
 GO TO 70
 60 istor = 1
 
!     OPEN INPUT FILES
 
 70 IF (lcore-nmat*sysbuf < 0) GO TO 420
 DO  i = 1,4
   IF (ia11(1,i) < 0) THEN
     GO TO    90
   ELSE IF (ia11(1,i) == 0) THEN
     GO TO   100
   END IF
   80 lcore = lcore - sysbuf
   CALL OPEN (*90,ia11(1,i),core(lcore+1),irew)
   CALL skprec (ia11(1,i),1)
   CYCLE
   90 ia11(1,i) = 0
   100 CONTINUE
 END DO
 
!     OPEN OUTPUT FILE
 
 CALL gopen (namea,core(ibuf),1)
 
!     FIX POINTERS -- SORT ON ABS VALUE
 
 k = ibufcp - 1
 l = ibufcp
 DO  i = 1,nrowa
   k = k + 1
   IF (core(k) < 0.0) THEN
     GO TO   110
   ELSE
     GO TO   120
   END IF
   110 core(l) = i
   l = l + 1
   120 CONTINUE
 END DO
 m = l - 1
 k = ibufcp
 DO  i = 1,nrowa
   130 IF (core(k)-i < 0.0) THEN
     GO TO   150
   ELSE IF (core(k)-i == 0.0) THEN
     GO TO   160
   END IF
   140 core(l) = i
   l = l + 1
   CYCLE
   150 IF (k == m) GO TO 140
   k = k + 1
   GO TO 130
   160 CONTINUE
 END DO
 
!     LOOP ON COLUMNS OF OUTPUT
 
 km = 0
 l2 = ibufcp
 l3 = ibufcp + zcpct
 DO  loop = 1,ncola1
   CALL bldpk (iotp,itypa,namea,0,0)
   IF (istor == 1) GO TO 190
   j  = (loop-1)/32 + ibufrp
   km = km + 1
   IF (km > 32) km = 1
   itemp = andf(core(j),two1(km))
   IF (km == 1) itemp = rshift(andf(core(j),two1(km)),1)
   IF (itemp /= 0) GO TO 180
   
!     IA11 AND IA21 BEING USED
   
   170 l1 = 0
   IF (l2 == l3-1) GO TO 200
   l2 = l2 + 1
   GO TO 200
   
!     IA12 AND IA22 BEING USED
   
   180 l1 = 2
   l3 = l3 + 1
   GO TO 200
   
!     USE ROW STORE
   
   190 IF (core(l2) == loop) GO TO 170
   IF (core(l3) == loop) GO TO 180
   GO TO 460
   
!     BEGIN ON SUBMATRICES
   
   200 io = 0
   DO  j = 1,2
     k = l1 + j
     IF (ia11(1,k) == 0) THEN
       GO TO   220
     END IF
     210 m = 20*j - 19
     CALL intpk (*220,ia11(1,k),BLOCK(m),iotp,1)
     io = io + j
     220 CONTINUE
   END DO
   IF (io == 0) THEN
     GO TO   380
   END IF
   
!     PICK UP NON ZERO
   
   230 ieol = 0
   jeol = 0
   ipos = 9999999
   jpos = 9999999
   iaz  = 1
   ibz  = 1
   nam1 = ia11(1,l1+1)
   nam2 = ia11(1,l1+2)
   IF (io-2 == 0) THEN
     GO TO   280
   END IF
   240 iaz  = 0
   250 IF (ieol == 0) THEN
     GO TO   260
   ELSE
     GO TO   370
   END IF
   260 CALL intpki (a11(1),i,nam1,BLOCK(1),ieol)
   k    = ibufcp + i - 1
   ipos = core(k)
   IF (io == 1) GO TO 310
   io   = 1
   280 ibz  = 0
   290 IF (jeol == 0) THEN
     GO TO   300
   ELSE
     GO TO   340
   END IF
   300 CALL intpki (b11(1),j,nam2,BLOCK(21),jeol)
   k = ibufcp + zcpct + j - 1
   jpos = core(k)
   310 IF (ipos-jpos < 0) THEN
     GO TO   350
   END IF
   
!     PUT IN B11
   
   320 DO  m = 1,ntypa
     ic11(m) = b11(m)
   END DO
   ii = jpos
   CALL zblpki
   GO TO 290
   340 jpos = 9999999
   ibz  = 1
   IF (iaz+ibz == 2) GO TO 380
   350 DO  m = 1,ntypa
     ic11(m) = a11(m)
   END DO
   ii = ipos
   CALL zblpki
   GO TO 250
   370 iaz  = 1
   ipos = 9999999
   IF (iaz+ibz /= 2) GO TO 320
   
!     OUTPUT COLUMN
   
   380 CALL bldpkn (namea,0,namea)
   
 END DO
 
!     DONE -- CLOSE OPEN MATRICES
 
 DO  i = 1,4
   IF (ia11(1,i) > 0) CALL CLOSE (ia11(1,i),1)
 END DO
 CALL CLOSE (namea,1)
 GO TO 500
 
 420 mn = -8
 GO TO 480
 440 WRITE  (nout,450) k
 450 FORMAT ('0*** SYSTEM OR USER ERROR, DUPLICATE GINO FILES AS ',  &
     'DETECTED BY MERGE ROUTINE - ',i5)
 nm = -37
 GO  TO 480
 460 mn = -7
 480 CALL mesage (mn,0,NAME)
 
 500 RETURN
END SUBROUTINE merge
