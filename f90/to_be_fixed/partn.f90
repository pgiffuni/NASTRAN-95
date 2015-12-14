SUBROUTINE partn (irp,icp,core)
     
 
 INTEGER, INTENT(IN)                      :: irp(1)
 INTEGER, INTENT(IN)                      :: icp(1)
 INTEGER, INTENT(IN OUT)                  :: core(1)
 EXTERNAL         rshift,andf
 INTEGER :: andf,two1, sysbuf,rshift
 
 DIMENSION  ias(7,4),head(2),block1(40),NAME(2)
 COMMON /parmeg/  namea,ncola,nrowa,iforma,itypa,ia(2),  &
     ia11(7),ia21(7),ia12(7),ia22(7),lcare,rule
 COMMON /system/  sysbuf
 COMMON /two   /  two1(32)
 COMMON /zntpkx/  a11(4),ii,ieol,IEOR
 EQUIVALENCE     (ias(1,1),ia11(1))
 DATA    iln   /  20 /, NAME / 4HPART,4HN   /
 
!     ZERO 6 AND 7 OF OUTPUT BLOCKS
 
 iotp  = itypa
 iopen = 0
 DO  i = 1,4
   DO  j = 6,7
     ias(j,i) = 0
   END DO
   IF (ias(1,i) == 0) THEN
     GO TO    40
   END IF
   20 IF (ias(5,i) /= itypa) iotp = 4
   iopen = iopen + 1
   DO  j = 2,5
     IF (ias(j,i) > 0) THEN
       GO TO    30
     ELSE
       GO TO   340
     END IF
   END DO
   ias(2,i) = 0
 END DO
 lcore  = lcare
 ibuf   = lcore- sysbuf + 1
 ibufcp = ibuf - nrowa
 ibufrp = ibufcp - (ncola+31)/32
 IF (ibufrp > 0) THEN
   GO TO    50
 ELSE
   GO TO   300
 END IF
 50 lcore = ibufrp - 1
 inorp = 0
 CALL ruler (rule,icp,zcpct,ocpct,core(ibufcp),nrowa,core(ibuf),1)
 IF (irp(1) == icp(1) .AND. irp(1) /= 0 .AND. nrowa == ncola) GO TO 60
 CALL ruler (rule,irp,zrpct,orpct,core(ibufrp),ncola,core(ibuf),0)
 GO TO 70
 60 inorp = 1
 lcore = ibufcp - 1
 
!     OPEN OUTPUT MATRICES
 
 70 IF (iopen*sysbuf > lcore) GO TO 300
 DO  i = 1,4
   IF (ias(1,i) == 0) THEN
     GO TO   100
   END IF
   80 lcore = lcore - sysbuf
   CALL OPEN  (*90,ias(1,i),core(lcore+1),1)
   CALL fname (ias(1,i),head)
   CALL WRITE (ias(1,i),head,2,1)
   CYCLE
   90 ias(1,i) = 0
 END DO
 
!     OPEN INPUT MATRIX
 
 CALL gopen (namea,core(ibuf),0)
 
!     LOOP FOR EACH COLUMN
 
 km = 0
 DO  loop = 1,ncola
   IF (inorp /= 0) GO TO 110
   
!     COLUMN PARTITION A SEQ. OF ZEROS AND ONES
   
   km = km + 1
   IF (km > 32) km = 1
   l = ibufrp + (loop-1)/32
   itemp = andf(core(l),two1(km))
   IF (km == 1) itemp = rshift(andf(core(l),two1(km)),1)
   IF (itemp == 0) THEN
     GO TO   130
   ELSE
     GO TO   120
   END IF
   110 l  = ibufcp + loop - 1
   IF (core(l) < 0.0) THEN
     GO TO   130
   END IF
   120 l1 = 2
   GO TO 140
   130 l1 = 0
   
!     BEGIN BLDPK ON TWO SUBS
   
   140 j = 0
   DO  l = 1,2
     k = l1 + l
     m = iln*(l-1) + 1
     IF (ias(1,k) == 0) THEN
       GO TO   160
     END IF
     150 CALL bldpk (iotp,ias(5,k),ias(1,k),block1(m),1)
     j = j + 1
   END DO
   IF (j == 0) THEN
     GO TO   260
   END IF
   
!     SEARCH COLUMN FOR NON-ZERO ELEMENTS
   
   170 CALL intpk (*230,namea,0,iotp,0)
   
!     LOOP FOR ROWS WITHIN COLUMN
   
   180 IF (ieol == 0) THEN
     GO TO   190
   ELSE
     GO TO   230
   END IF
   190 CALL zntpki
   
!     COMPUTE ROW POSITION AND OUTPUT MATRIX
   
   l = ibufcp + ii - 1
   ipos = IABS(core(l))
   IF (core(l) < 0.0) THEN
     GO TO   200
   ELSE
     GO TO   210
   END IF
   200 m1 = l1 + 1
   m  = 1
   GO TO 220
   210 m1 = l1 + 2
   m  = iln+ 1
   220 IF (ias(1,m1) == 0) GO TO 180
   CALL bldpki (a11(1),ipos,ias(1,m1),block1(m))
   GO TO 180
   230 DO  l = 1,2
     k = l + l1
     m = iln*(l-1) + 1
     IF (ias(1,k) == 0) THEN
       GO TO   250
     END IF
     240 CALL bldpkn (ias(1,k),block1(m),ias(1,k))
   END DO
   CYCLE
   260 CALL skprec (namea,1)
 END DO
 
!     ALL DONE - CLOSE OPEN MATRICES
 
 CALL CLOSE (namea,1)
 DO  i = 1,4
   IF (ias(1,i) == 0) THEN
     GO TO   290
   END IF
   280 CALL CLOSE (ias(1,i),1)
 END DO
 RETURN
 
 300 ipm1 =-8
 310 CALL mesage (ipm1,ipm2,NAME)
 340 ipm1 =-7
 GO TO 310
END SUBROUTINE partn
