SUBROUTINE lamx
     
!     LAMX MAKES OR EDITS THE LAMA DATA BLOCK
 
!     LAMX  EDIT,LAMA/LAMB/C,Y,NLAM=0 $
!     IF NLAM LT 0 MAKE LAMB A MATRIX OF 5 COLUMNS
!     LAMA  OMEGA FREQ GM GS
!     UNTIL GM = 0.0
 
 
 INTEGER :: sysbuf,ist(10),trl(7),bufa,bufb,bufe,edit
 
 DIMENSION d(3),z(7)
 
 COMMON /BLANK / nlam
 COMMON /zzzzzz/ iz(1)
 COMMON /system/ sysbuf,nout
 COMMON /condas/ pi,twopi
 COMMON /unpakx/ ito,ii,ie,incr
 COMMON /packx / ityin,ityout,iii,nnn,incr1
 COMMON /output/ hdg(96)
 
 EQUIVALENCE     (d(1),a),(d(2),b),(d(3),c)
 EQUIVALENCE     (z(1),iz(1))
 
 DATA    edit  , lama,lamb /101,102,201/
 DATA    ist   / 21,6,7*0,7/
 DATA    lma   / 1/,   ied /1/,  iz2 /2/
 DATA    nam   / 4HLAMX    /
 
!     INITILIZE AND DECIDE MODE OF OPERATIONS
 
 icore  = korsz(z)
 trl(1) = lama
 CALL rdtrl (trl)
 IF (trl(1) < 0) lma = 0
 trl(1) = edit
 CALL rdtrl (trl)
 IF (trl(1) < 0) ied = 0
 ncol = trl(2)
 IF (ncol == 0) ied = 0
 IF (lma == 0 .AND. ied == 0) GO TO 1000
 ito = 1
 ii  = 1
 incr= 1
 ie  = trl(3)
 IF (ie > 3) ie = 3
 b = 0.0
 c = 0.0
 bufb = icore - sysbuf
 CALL gopen (lamb,z(bufb),1)
 IF (lma == 0) GO TO 200
 bufa = bufb - sysbuf
 CALL gopen (lama,z(bufa),0)
 IF (nlam < 0) GO TO 500
 bufe = bufa - sysbuf
 IF (ied == 0) GO TO 5
 
!      EDITING LAMA FROM EDIT
 
 CALL gopen (edit,z(bufe),0)
 
!     WRITE HEADER
 
 5 CALL READ (*10,*10,lama,z,bufe,1,nwr)
 10 CALL WRITE (lamb,z,nwr,1)
 IF (ied == 0) GO TO 100
 
!     MAKE RECORDS
 
 j = 0
 DO  i = 1,ncol
   CALL READ (*60,*60,lama,z,7,0,nwr)
   CALL unpack (*40,edit,d)
   IF (a == 0.0 .AND. b == 0.0 .AND. c == 0.0) GO TO 40
   IF (c < 0.0) CYCLE
   z(5) = z(5)*(1.0+b) + a
   z(4) = z(5)*twopi
   z(3) = z(4)*z(4)
   IF (c /= 0.0) z(6) = c
   z(7) = z(6)*z(3)
   40 j    = j+1
   iz(1)= j
   IF (nlam <= 0) GO TO 45
   IF (j > nlam) GO TO 180
   45 CALL WRITE (lamb,z,7,0)
 END DO
 60 GO TO 180
 
!     COPY LAMA TO LAMB FOR NLAM RECORDS
 
 100 IF (nlam == 0) GO TO 190
 j = nlam
 m = 7*nlam
 CALL READ (*180,*110,lama,z,m,0,nwr)
 CALL WRITE (lamb,z,7*nlam,0)
 GO TO 180
 110 CALL WRITE (lamb,z,nwr,0)
 180 trl(1) = lamb
 trl(2) = j
 CALL wrttrl (trl)
 190 CALL CLOSE (lama,1)
 CALL CLOSE (lamb,1)
 CALL CLOSE (edit,1)
 GO TO 1000
 
!      MAKE A NEW LAMB
 
 200 bufe = bufb - sysbuf
 CALL gopen (edit,z(bufe),0)
 IF (nlam > 0) ncol = MIN0(ncol,nlam)
 
!     WRITE HEADER
 
 CALL WRITE (lamb,ist,50,0)
 CALL WRITE (lamb,hdg,96,1)
 
!     MAKE RECORDS
 
 DO  i = 1,ncol
   CALL unpack (*310,edit,d)
   GO TO 210
   310 d(1) = 0.0
   d(2) = 0.0
   d(3) = 0.0
   210 iz(  1) = i
   iz(iz2) = i
   z(5) = a
   z(4) = twopi*a
   z(3) = z(4)*z(4)
   z(6) = c
   z(7) = c*z(3)
   CALL WRITE (lamb,z,7,0)
 END DO
 j = ncol
 GO TO 180
 
!     BUILD LAMB AS A MATRIX
 
 500 trl(1) = lamb
 trl(2) = 0
 trl(4) = 1
 trl(5) = 1
 trl(6) = 0
 trl(7) = 0
 ityin  = 1
 ityout = 1
 iii    = 1
 incr1  = 7
 CALL fwdrec (*190,lama)
 CALL READ (*190,*510,lama,z,bufa,0,nwr)
 CALL mesage (8,0,nam)
 GO TO 190
 510 nloop = 0
 DO  i = 1,nwr,7
   IF (z(i+5) == 0.0) EXIT
   nloop = nloop +1
 END DO
 530 IF (nloop == 0) GO TO 190
 trl(3) = nloop
 nnn = nloop
 l   = 3
 DO  i = 1,5
   CALL pack (z(l),lamb,trl)
   l = l + 1
 END DO
 CALL wrttrl (trl)
 GO TO 190
 1000 RETURN
END SUBROUTINE lamx
