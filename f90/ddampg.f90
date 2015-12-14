SUBROUTINE ddampg
     
!     DDAMPG  MP,PVW/PG/V,N,NMODES/V,N,NDIR $
 
!     MP IS MGG*PHIG, PVW IS (PF)*SSDV*OMEGA, PARTICIPATION FACTORS X
!     SHOCK SPECTRUM DESIGN VALUES X RADIAN FREQUENCIES.
!     MP IS (NXM).  IF PVW IS A VECTOR (MX1), WE WANT TO MULTIPLY THE
!     ITH. TERM INTO THE ITH. COLUMN OF MP.  PG IS THEN NXM.
!     IF PVW IS A MATRIX (MXL), WE REPEAT THE PREVIOUS COMPUTATION FOR
!     EACH OF THE L VECTORS, MAKING PG (NX(MXL)).
!     NMODES IS NUMBER OF MODES. NDIR IS NUMBER OF SHOCK DIRECTIONS
 
 INTEGER :: mp,pvw,pg,buf1,buf2,buf3,FILE
 DIMENSION       nam(2),mcb(7)
 COMMON /unpakx/ jout,iii,nnn,jncr
 COMMON /packx / iin,iout,ii,nn,incr
 COMMON /system/ ibuf(80)
 COMMON /BLANK / nmodes,ndir
 COMMON /zzzzzz/ z(1)
 DATA    mp,pvw, pg /101,102,201/
 DATA    nam   / 4HDDAM,4HPG    /
 
!     SET UP OPEN CORE AND BUFFERS
 
 lcore = korsz(z)
 buf1  = lcore - ibuf(1) + 1
 buf2  = buf1 - ibuf(1)
 buf3  = buf2 - ibuf(1)
 lcore = buf3 - 1
 IF (lcore <= 0) GO TO 1008
 
!     PICK UP ROW AND COLUMN STATISTICS AND SET PACK/UNPACK PARAMETERS
 
 mcb(1) = mp
 CALL rdtrl (mcb)
 ncolmp = mcb(2)
 nmodes = ncolmp
 nrowmp = mcb(3)
 mcb(1) = pvw
 CALL rdtrl (mcb)
 ncolpv = mcb(2)
 ndir   = ncolpv
 nrowpv = mcb(3)
 mcb4   = mcb(4)
 mcb5   = mcb(5)
 
 
 IF (lcore < nrowpv+nrowmp) GO TO 1008
 IF (ncolmp /= nrowpv) GO TO 1007
 mcb(1) = pg
 mcb(2) = 0
 mcb(3) = nrowmp
 mcb(4) = mcb4
 mcb(5) = mcb5
 mcb(6) = 0
 mcb(7) = 0
 
 jout = 1
 iin  = 1
 iout = 1
 ii   = 1
 iii  = 1
 nn   = nrowmp
 incr = 1
 jncr = 1
 
 CALL gopen (mp,z(buf1),0)
 CALL gopen (pvw,z(buf2),0)
 CALL gopen (pg,z(buf3),1)
 
 DO  ijk = 1,ncolpv
   nnn = nrowpv
   CALL unpack (*20,pvw,z(1))
   GO TO 60
   
!     NULL COLUMN FOR PVW-WRITE OUT NCOLMP ZERO COLUMNS OF LENGTH NROWMP
   
   20 DO  k = 1,nrowmp
     z(k) = 0.
   END DO
   DO  k = 1,ncolmp
     CALL pack (z,pg,mcb)
   END DO
   GO TO 125
   
   60 DO  j = 1,ncolmp
     nnn = nrowmp
     CALL unpack (*80,mp,z(nrowpv+1))
     GO TO 100
     
     80 DO  k = 1,nrowmp
       z(nrowpv+k) = 0.
     END DO
     GO TO 115
     
     100 DO  k = 1,nrowmp
       isub = nrowpv + k
       z(isub) = z(isub)*z(j)
     END DO
     115 CALL pack (z(nrowpv+1),pg,mcb)
   END DO
   125 CALL REWIND (mp)
   FILE = mp
   CALL fwdrec (*1002,mp)
 END DO
 
 CALL wrttrl (mcb)
 CALL CLOSE (mp,1)
 CALL CLOSE (pvw,1)
 CALL CLOSE (pg,1)
 
 RETURN
 
!     FATAL ERRORS
 
 1002 n = -2
 GO TO 1010
 1007 n = -7
 GO TO 1010
 1008 n = -8
 FILE = 0
 1010 CALL mesage (n,FILE,nam)
 RETURN
END SUBROUTINE ddampg
