SUBROUTINE gencos
     
!     GENCOS  GENERATES DIRECTION COSINE MATRIX, UP TO NX3, FOR DDAM.
!     THE SHOCK DIRECTIONS ARE GIVEN BY A COORDINATE SYSTEM (PROBABLY
!     RECTANGULAR, BUT NOT NECESSARILY) DEFINED ON A CORDIJ CARD.
!     THE ID OF THAT SYSTEM MUST BE SPECIFIED BY PARAM SHOCK ID.
!     THE DIRECTIONS OF INTEREST MUST BE SPECIFIED ON A PARAM DIRECT DIR
!     CARD WHERE DIR=1,2,3,12,13,23,OR 123 GIVING THE SHOCK DIRECTIONS
!     DESIRED IN THE SHOCK COORDINATE SYSTEM.  (DEFAULT IS 123)  WE WILL
!     BE CONVERTING A ROW VECTOR IN THE GLOBAL SYSTEM TO A ROW VECTOR IN
!     THE SHOCK SYSTEM.  TO CONVERT A COLUMN VECTOR FROM GLOBAL TO SHOCK
!     FIRST CONVERT TO BASIC.  THEN TRANSFORM FROM BASIC TO SHOCK, I.E.
!     (VECTOR-SHOCK) = (TRANSPOSE(T-SHOCK TO BASIC))*
!                      (T-GLOBAL TO BASIC)*(VECTOR-GLOBAL)
!     BUT BECAUSE WE ARE TRANSFORMING ROW VECTORS, THE EQUATION IS
!     TRANSPOSED . NSCALE =1 MEANS THERE ARE SCALAR POINTS,=0 MEANS NO
 
!     GENCOS    BGPDT,CSTM/DIRCOS/C,Y,SHOCK=0/C,Y,DIRECT=123/
!               V,N,LUSET/V,N,NSCALE $
 
 LOGICAL :: REC,all
 INTEGER :: bgpdt,cstm,dircos,buf1,FILE,shock,DIRECT,otpe
 DIMENSION       nam(2),mcb(7),iz(1),tshock(9),coord(4),icoord(4),  &
     tpoint(9),tfinal(9),idir(3),isub(3)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /BLANK / shock,DIRECT,luset,nscale
 COMMON /system/ ibuf,otpe
 COMMON /packx / in,iout,ii,nn,incr
 COMMON /zzzzzz/ z(1)
 EQUIVALENCE     (z(1),iz(1)), (coord(1),icoord(1))
 DATA    bgpdt , cstm,dircos   / 101,102,201 /
 DATA    nam   / 4HGENC,4HOS   /
 
!     OPEN CORE AND BUFFERS
 
 lcore = korsz(z)
 buf1  = lcore - ibuf + 1
 lcore = buf1  - 1
 IF (lcore <= 0) GO TO 1008
 
!     CHECK FOR SCALAR POINTS AND SET NSCALE
 
 mcb(1) = bgpdt
 CALL rdtrl (mcb)
 npts = mcb(2)
 CALL gopen (bgpdt,z(buf1),0)
 DO  i = 1,npts
   CALL fread (bgpdt,coord,4,0)
   IF (icoord(1) == -1) GO TO 2
 END DO
 nscale = 0
 GO TO 3
 2 nscale = 1
 3 CALL CLOSE (bgpdt,1)
 
 IF (DIRECT >= 1  .AND. DIRECT <= 3) GO TO 5
 IF (DIRECT /= 12 .AND. DIRECT /= 13 .AND. DIRECT /= 23 .AND.  &
     DIRECT /= 123) GO TO 500
 5 IF (shock <  0) GO TO 500
 ncstm  = 0
 ncount = 0
 all    = .false.
 REC    = .false.
 ndir   = 2
 IF (DIRECT <= 3) ndir = 1
 IF (DIRECT == 123)  ndir = 3
 IF (luset*ndir > lcore) GO TO 1008
 
 SELECT CASE ( ndir )
   CASE (    1)
     GO TO 6
   CASE (    2)
     GO TO 7
   CASE (    3)
     GO TO 8
 END SELECT
 
 6 idir(1) = DIRECT
 GO TO 9
 
 7 IF (DIRECT == 23) GO TO 175
 idir(1) = 1
 idir(2) = 2
 IF (DIRECT == 13) idir(2) = 3
 GO TO 9
 175 idir(1) = 2
 idir(2) = 3
 GO TO 9
 
 8 idir(1) = 1
 idir(2) = 2
 idir(3) = 3
 9 CONTINUE
 
 
!     READ CSTM FOR FETCHING TRANSFORMATION MATRICES
 
 CALL OPEN (*10,cstm,z(buf1),0)
 GO TO 30
 
!     CSTM IS PURGED.  SO, GLOBAL SYSTEM IS BASIC AND SHOCK SYSTEM MUST
!     BE ALSO.  IF SHOCK SYSTEM IS NOT 0, FATAL MESSAGE.  IF IT IS 0,
!     THEN NEED ONLY IDENTITIES.
 
 10 IF (shock == 0) GO TO 25
 WRITE  (otpe,20) ufm
 20 FORMAT (a23,', IN GENCOS, CSTM IS PURGED AND SHOCK COORDINATE ',  &
     'SYSTEM IS NOT BASIC')
 CALL mesage (-61,0,0)
 
!     EVERYTHING IS BASIC - CHECK FOR SCALAR POINTS - IF THEY EXIST,
!     WE MUST READ BGPDT
 
 25 IF (nscale == 1) GO TO 55
 all  = .true.
 isys = 0
 GO TO 130
 
 30 FILE = cstm
 CALL fwdrec (*1002,cstm)
 CALL READ (*1002,*40,cstm,z,lcore,0,ncstm)
 GO TO 1008
 40 CALL CLOSE (cstm,1)
 
!     CHECK FOR ENOUGH OPEN CORE
 
 IF (ncstm+luset*ndir > lcore) GO TO 1008
 CALL pretrs (z(1),ncstm)
 
!     IF SHOCK COORDINATE SYSTEM IS RECTANGULAR, LET'S GET THE TRANS-
!     FORMATION MATRIX ONCE SINCE IT WILL NOT BE POINT-DEPENDENT.
 
 IF (shock == 0) GO TO 55
 DO  i = 1,ncstm,14
   IF (shock /= iz(i)) CYCLE
   IF (iz(i+1) /= 1) GO TO 60
   
!     RECTANGULAR
   
   REC = .true.
   DO  j = 1,9
     tshock(j) = z(i+j+4)
   END DO
   GO TO 60
 END DO
 
!     CAN'T FIND SHOCK COORDINATE SYSTEM
 
 CALL mesage (-30,25,shock)
 
!     SHOCK IS BASIC
 
 55 REC = .true.
 DO  i = 1,9
   tshock(i) = 0.
 END DO
 tshock(1) = 1.
 tshock(5) = 1.
 tshock(9) = 1.
 
!     OPEN BGPDT TO GET GRID POINT OUTPUT COORDINATE SYSTEMS AND
!     BASIC COORDINATES
 
 60 CALL gopen (bgpdt,z(buf1),0)
 FILE = bgpdt
 70 CALL READ (*1002,*210,bgpdt,coord,4,0,iwords)
 isys = icoord(1)
 IF (icoord(1) == -1) GO TO 150
 IF (icoord(1) /=  0) GO TO 80
 
!     IDENTITY - BASIC SYSTEM
 
 DO  i = 1,9
   tpoint(i) = 0.
 END DO
 tpoint(1) = 1.
 tpoint(5) = 1.
 tpoint(9) = 1.
 GO TO 85
 
!     FETCH GLOBAL-TO-BASIC MATRIX FOR THIS POINT
 
 80 CALL transs (coord,tpoint)
 
!     IF SHOCK IS NOT RECTANGULAR, FETCH SHOCK-TO-BASIC FOR THIS POINT
 
 85 IF (REC) GO TO 90
 icoord(1) = shock
 CALL transs (coord,tshock)
 
!     THE MATRIX WE NEED IS (TRANSPOSE(TPOINT))*(TSHOCK)
 
 90 IF (shock == 0) GO TO 100
 IF (isys  == 0) GO TO 110
 
!     NEITHER MATRIX IS NECESSARILY IDENTITY
 
 CALL gmmats (tpoint,3,3,1,tshock,3,3,0,tfinal)
 GO TO 150
 
!     TSHOCK IS IDENTITY
 
 100 IF (isys == 0) GO TO 130
 
!     BUT TPOINT IS NOT
 
 tfinal(1) = tpoint(1)
 tfinal(2) = tpoint(4)
 tfinal(3) = tpoint(7)
 tfinal(4) = tpoint(2)
 tfinal(5) = tpoint(5)
 tfinal(6) = tpoint(8)
 tfinal(7) = tpoint(3)
 tfinal(8) = tpoint(6)
 tfinal(9) = tpoint(9)
 GO TO 150
 
!     TPOINT IS IDENTITY, BUT TSHOCK IS NOT
 
 110 DO  i = 1,9
   tfinal(i) = tshock(i)
 END DO
 GO TO 150
 
!     BOTH ARE IDENTITY
 
 130 DO  i = 1,9
   tfinal(i) = 0.
 END DO
 tfinal(1) = 1.
 tfinal(5) = 1.
 tfinal(9) = 1.
 
!     STORE TFINAL BY INTERNAL ORDERING AND DIRECTIONS REQUESTED START-
!     ING AT Z(NCSTM+1) - MAKE UP TO 3 COLUMNS OF LUSET EACH
 
 150 isub(1) = ncstm   + ncount
 isub(2) = isub(1) + luset
 isub(3) = isub(2) + luset
 
 DO  i = 1,ndir
   ip   = idir(i)
   jsub = isub(i)
   IF (isys == -1) GO TO 195
   z(jsub+1) = tfinal(ip  )
   z(jsub+2) = tfinal(ip+3)
   z(jsub+3) = tfinal(ip+6)
   z(jsub+4) = 0.
   z(jsub+5) = 0.
   z(jsub+6) = 0.
   CYCLE
   
!     SCALAR
   
   195 z(jsub+1) = 1.
 END DO
 
!     GO BACK FOR ANOTHER POINT
 
 ncount = ncount + 6
 IF (isys == -1)  ncount = ncount - 5
 IF (.NOT.all) GO TO 70
 IF (ncount == luset) GO TO 210
 GO TO 150
 
!     DONE WITH ALL POINTS - PACK RESULTS
 
 210 IF (.NOT.all) CALL CLOSE (bgpdt,1)
 CALL gopen (dircos,z(buf1),1)
 in   = 1
 iout = 1
 ii   = 1
 nn   = luset
 incr = 1
 mcb(1) = dircos
 mcb(2) = 0
 mcb(3) = luset
 mcb(4) = 2
 mcb(5) = 1
 mcb(6) = 0
 mcb(7) = 0
 DO  i = 1,ndir
   jsub = ncstm + luset*(i-1)
   CALL pack (z(jsub+1),dircos,mcb)
 END DO
 
 CALL CLOSE (dircos,1)
 CALL wrttrl (mcb)
 RETURN
 
 500 WRITE  (otpe,510) ufm,shock,DIRECT
 510 FORMAT (a23,', SHOCK AND DIRECT ARE',2I10, /10X,'RESPECTIVELY. ',  &
     'SHOCK MUST BE NONNEGATIVE AND DIRECT MUST BE EITHER 1,2',  &
     ',3,12,13,23, OR 123')
 CALL mesage (-61,0,0)
 
 1002 n = -2
 GO TO 1010
 1008 n = -8
 FILE = 0
 1010 CALL mesage (n,FILE,nam)
 RETURN
END SUBROUTINE gencos
