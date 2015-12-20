SUBROUTINE mce1d
     
!     MCE1D SOLVES FOR GM IN THE MATRIX EQUATION RM*GM = -RN
!     WHERE RM IS A DIAGONAL MATRIX.
 
 INTEGER :: sysbuf,eol    ,eor   ,TYPE  ,rdp   ,bcd  ,rm   , rn    ,gm
 REAL :: z(1)  ,a(1)   ,b(1)
 DOUBLE PRECISION :: zd   ,ad     ,bd
 DIMENSION       bcd(2),mcb1(7),mcb2(7)
 COMMON /BLANK / uset  ,rg     ,gm    ,scr1  ,scr2  ,scr3  ,rm  ,  &
     rn    ,l      ,u     ,mcb(7)
 COMMON /zntpkx/ ad (2),i      ,eol   ,eor
 COMMON /zblpkx/ bd (2),j
 COMMON /zzzzzz/ zd (1)
 COMMON /system/ sysbuf,skp(53),ipr
 EQUIVALENCE     (mcb(2),ncol) ,(ad(1),a(1)) ,(mcb(5),TYPE) ,  &
     (zd(1) ,z(1)) ,(bd(1),b(1)) ,(mcb1(2),ncol1)
 DATA    bcd   , rdp   /4HMCE1 ,4HD   ,  2   /
 
!     OPEN RM MATRIX,SKIP HEADER RECORD AND READ MATRIX CONTROL BLOCK
 
 nz = korsz(z)
 n  = nz - sysbuf
 CALL gopen (rm,z(n+1),0)
 mcb(1) = rm
 CALL rdtrl (mcb)
 
!     FORM -RM
 
 ncol = mcb(2)
 DO  k = 1,ncol
   CALL intpk (*83,rm,0,rdp,0)
   CALL zntpki
   IF (i /= k) GO TO 84
   zd(k) = -ad(1)
 END DO
 CALL CLOSE (rm,1)
 
!     OPEN RN MATRIX,SKIP HEADER RECORD AND READ MATRIX CONTROL BLOCK
 
 CALL gopen (rn,z(n+1),0)
 mcb1(1) = rn
 CALL rdtrl (mcb1)
 
!     SET UP MATRIX CONTROL BLOCK BLOCK FOR GM
 
 CALL makmcb (mcb2,gm,mcb1(3),mcb1(4),ipr)
 
!     OPEN OUTPUT FILE FOR GM AND WRITE HEADER RECORD
 
 n1 = n - sysbuf
 CALL gopen (gm,z(n1+1),1)
 
!     FORM GM = -RM(-1)*RN
 
 ncol1 = mcb1(2)
 DO  k = 1,ncol1
   CALL bldpk (rdp,ipr,gm,0,0)
   CALL intpk (*62,rn,0,rdp,0)
   61 CALL zntpki
   j = i
   bd(1) = ad(1)/zd(j)
   CALL zblpki
   IF (eol == 0.0) THEN
     GO TO    61
   ELSE
     GO TO    62
   END IF
   62 CALL bldpkn (gm,0,mcb2)
 END DO
 
!     CLOSE GM AND RM FILES AND WRITE TRAILER FOR GM
 
 CALL CLOSE (gm,1)
 CALL CLOSE (rn,1)
 CALL wrttrl (mcb2)
 RETURN
 
!     CALL MESSAGE WRITER IF FATAL ERROR DETECTED
 
 83 l = -5
 GO TO 86
 84 l = -16
 86 CALL mesage (l,rm,bcd)
 RETURN
END SUBROUTINE mce1d
