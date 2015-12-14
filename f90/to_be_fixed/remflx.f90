SUBROUTINE remflx (ngrids)
     
!     CHECK FOR REMFLUX IN MAGNETIC FIELD PROBLEMS WHEN COMPUTING
!     PROLATE SPHEROIDAL COEFFICIENTS
 
 
 INTEGER, INTENT(IN)                      :: ngrids
 LOGICAL :: remfl,hitone
 INTEGER :: remfld,scr1,buf1,buf3,buf2,hest,mcb(7),FILE,  &
     pointr(6,19),ipoint(32),dit,ditfil,eltype,estwds
 DIMENSION       nam(2),rem(3),ecpt(200),necpt(200),iz(1),g(3,3), iwork(3,3)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /gpta1 / nelems,last,incr,NE(1)
 COMMON /matin / matid,inflag,eltemp,stress,sinth,costh
 COMMON /hmtout/ xmat(6)
 COMMON /hmatdd/ iihmat,nnhmat,mptfil,ditfil
 COMMON /unpakx/ jout,ii,nn,jncr
 COMMON /system/ ibuf,iout
 COMMON /zzzzzz/ z(1)
 COMMON /biot  / dum(10),buf1,remfl,lcore
 EQUIVALENCE     (z(1),iz(1)),(ecpt(1),necpt(1))
 DATA    remfld, hest,mpt,dit, scr1/ 107   , 108 ,109,110, 301 /
 DATA    nam   / 4HREMF  ,4HLX     /
 
!                   TYPE  ISIL   MID   ITH NGRIDS ITEMP
 
 DATA    pointr/ 1,    2,    4,    0,    2,    17,  &
     3,    2,    4,    0,    2,    16, 6,    2,    6,    5,    3,    27,  &
     9,    2,    6,    5,    3,    21, 10,    2,    4,    0,    2,    17,  &
     16,    2,    7,    6,    4,    26, 17,    2,    6,    5,    3,    21,  &
     18,    2,    7,    6,    4,    26, 19,    2,    7,    6,    4,    32,  &
     34,    2,   16,    0,    2,    42, 36,    2,    6,    5,    3,    19,  &
     37,    2,    7,    6,    4,    24, 39,    3,    2,    0,    4,    23,  &
     40,    3,    2,    0,    6,    33, 41,    3,    2,    0,    8,    43,  &
     42,    3,    2,    0,    8,    43, 65,    2,   10,    0,    8,    48,  &
     66,    2,   22,    0,   20,   108, 67,    2,   34,    0,   32,   168  /
 
 remfl  = .false.
 mcb(1) = remfld
 CALL rdtrl (mcb)
 
!     CHECK FOR ANY REMFLUX
 
 IF (mcb(6) == 0) RETURN
 
!     TOO BAD
 
 remfl  = .true.
 ncol   = mcb(2)
 nrow   = mcb(3)
 ncount = nrow/3
 
!     BRING IN MATERIALS SINCE H=B/MU
 
 iihmat = ngrids
 nnhmat = lcore
 mptfil = mpt
 ditfil = dit
 CALL prehma (z)
 nextz  = nnhmat + 1
 
 buf2   = buf1 - ibuf
 buf3   = buf2 - ibuf
 
!     SET UP POINTERS
!     IHC = START OF RESULTS HC = B/MU
!     IREM= REMFL COLUMN
!     ICT = COUNTER FOR NUMBER OF ELEMENTS AT EACH PROLATE GRID (FOR
!     AVERAGING
 
 ihc  = nextz
 irem = ihc + 3*ngrids
 ict  = irem + nrow
 IF (buf3 < ict+ngrids) GO TO 1008
 
 CALL gopen (scr1,z(buf1),1)
 CALL gopen (remfld,z(buf2),0)
 CALL gopen (hest,z(buf3),0)
 
 ii = 1
 nn = nrow
 jncr = 1
 jout = 1
 jcount = 0
 
 3 DO  i = 1,ngrids
   iz(ict+i) = 0
 END DO
 n3 = 3*ngrids
 DO  i = 1,n3
   z(ihc+i) = 0.
 END DO
 
!     UNPACK A COULMN OF REMFLD
 
 jcount = jcount + 1
 CALL unpack (*20,remfld,z(irem+1))
 GO TO 40
 
!     ZERO COLUMN
 
 20 DO  i = 1,n3
   z(ihc+i) = 0.
 END DO
 GO TO 130
 
!     SINCE THE ELEMENT DATA DO NOT CHANGE WITH REMFLD COLIMN, THIS IS
!     NOT NECESSARILY THE BEST KIND OF LOOPING. BUT OTHER WAYS WOULD
!     NEED MORE CORE AND IF THERE IS MORE THAN ONE REMFLUX CASE, IT
!     WOULD BE A SURPRISE
 
 40 FILE  = hest
 kount = 0
 45 CALL READ (*100,*1003,hest,eltype,1,0,iflag)
 idx   = (eltype-1)*incr
 estwds= NE(idx+12)
 
!     PICK UP ELEMENT TYPE INFO
 
 DO  i = 1,19
   jel = i
   IF (eltype-pointr(1,i) < 0.0) THEN
     GO TO   500
   ELSE IF (eltype-pointr(1,i) == 0.0) THEN
     GO TO    60
   ELSE
     GO TO    50
   END IF
 END DO
 GO TO 500
 
 60 isil  = pointr(2,jel)
 imid  = pointr(3,jel)
 ith   = pointr(4,jel)
 igrids= pointr(5,jel)
 itemp = pointr(6,jel)
 
 65 CALL READ (*1002,*45,hest,ecpt,estwds,0,iflag)
 
!     PICK UP REMFLUX FOR THIS ELEMENT
 
 kount = kount + 1
 nhit  = 0
 hitone= .false.
 DO  i = 1,igrids
   ipoint(i) = 0
 END DO
 DO  i = 1,ngrids
   DO  j = 1,igrids
     ipt = necpt(isil+j-1)
     IF (ipt == iz(i)) GO TO 67
   END DO
   CYCLE
   
!     MATCH
   
   67 hitone = .true.
   nhit   = nhit + 1
   iz(ict+i) = iz(ict+i) + 1
   ipoint(j) = i
   IF (nhit == igrids) GO TO 69
 END DO
 IF (.NOT.hitone) GO TO 65
 69 CONTINUE
 
 isub   = irem + 3*(kount-1)
 rem(1) = z(isub+1)
 rem(2) = z(isub+2)
 rem(3) = z(isub+3)
 
!     PICK UP MATERIALS
 
 matid  = necpt(imid)
 eltemp = ecpt(itemp)
 inflag = 3
 sinth  = 0.
 costh  = 0.
 CALL hmat (necpt(1))
 g(1,1) = xmat(1)
 g(1,2) = xmat(2)
 g(1,3) = xmat(3)
 g(2,1) = xmat(2)
 g(2,2) = xmat(4)
 g(2,3) = xmat(5)
 g(3,1) = xmat(3)
 g(3,2) = xmat(5)
 g(3,3) = xmat(6)
 
!     FOR COMMENTS ON MATERIALS SEE EM2D
 
 IF (ith == 0) GO TO 80
 angle = ecpt(ith)*0.017453293
 IF (xmat(3) == 0. .AND. xmat(5) == 0.) GO TO 70
 GO TO 80
 70 IF (ABS(angle) <= .0001) GO TO 80
 s   = SIN(angle)
 c   = COS(angle)
 csq = c*c
 ssq = s*s
 cs  = c*s
 x2  = 2.*cs*xmat(2)
 g(1,1) = csq*xmat(1) - x2 + ssq*xmat(4)
 g(1,2) = cs*(xmat(1) - xmat(4)) + (csq-ssq)*xmat(2)
 g(2,2) = ssq*xmat(1) + x2 + csq*xmat(4)
 g(2,1) = g(1,2)
 g(3,3) = xmat(6)
 g(1,3) = 0.
 g(2,3) = 0.
 g(3,1) = 0.
 g(3,2) = 0.
 
!     SINCE MAT5 INFO FOR TRAPRG,TRIARG IS GIVEN IN X-Y ORDER,
!     INETRCHANGE YA AND Z
 
 temp   = g(2,2)
 g(2,2) = g(3,3)
 g(3,3) = temp
 temp   = g(1,2)
 g(1,2) = g(1,3)
 g(1,3) = temp
 g(2,1) = g(1,2)
 g(3,1) = g(1,3)
 
!     SOLVE MU*H = B
 
 80 CALL invers (3,g,3,rem,1,det,ising,iwork)
 IF (ising == 2) GO TO 510
 
!     REM NOW HAS HC- CHECK POINTER LIST TO SEE WHICH GRIDS ARE ON THE
!     SPHEROID AND ADD REMFLUX T THOSE ALREADY ACCUMULATED
 
 DO  i = 1,igrids
   IF (ipoint(i) == 0) CYCLE
   isub = ihc + 3*(ipoint(i)-1)
   DO  j = 1,3
     z(isub+j) = z(isub+j) + rem(j)
   END DO
 END DO
 
!     GO BACK FOR ANOTHER ELEMEENT
 
 GO TO 65
 
!     DONE WITH ALL TYPES-AVERAGE THE RESULTS BY NUMBER OF ELEMENTS AT
!     EACH
 
 100 DO  i = 1,ngrids
   den = FLOAT(iz(ict+i))
   IF (den == 0.) CYCLE
   isub = 3*(i-1) + ihc
   DO  j = 1,3
     z(isub+j) = z(isub+j)/den
   END DO
 END DO
 
!     WRITE RESULTS TO SCR1
 
 130 CALL WRITE (scr1,z(ihc+1),3*ngrids,1)
 
!     GO BACK FOR ANOTHER REMFLD RECORD
 
 IF (jcount == ncol) GO TO 140
 CALL REWIND (hest)
 CALL fwdrec (*1002,hest)
 GO TO 3
 
!     DONE
 
 140 CALL CLOSE (scr1,1)
 mcb(1) = scr1
 mcb(2) = ncol
 mcb(3) = 3*ngrids
 DO  i = 4,7
   mcb(i) = 0
 END DO
 CALL wrttrl (mcb)
 CALL CLOSE (hest,1)
 CALL CLOSE (remfld,1)
 RETURN
 
 500 WRITE  (iout,501) ufm
 501 FORMAT (a23,', ILLEGAL ELEMENT TYPE IN REMFLX')
 GO TO  1061
 510 WRITE  (iout,511) ufm,matid
 511 FORMAT (a23,', MATERIAL',i9,' IS SINGULAR IN REMFLX')
 GO TO  1061
 
 1002 n = -2
 GO TO 1010
 1003 n = -3
 GO TO 1010
 1008 n = -8
 FILE = 0
 1010 CALL mesage (n,FILE,nam)
 1061 CALL mesage (-61,0,0)
 RETURN
END SUBROUTINE remflx
