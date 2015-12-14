SUBROUTINE dmpy (z,zd)
     
!     DMPY WILL PRE OR POST MULTIPLY AN ARBITRARY MATRIX BY A DIAGONAL
!     MATRIX.
 
!     FILEA = MATRIX CONTROL BLOCK FOR DIAGONAL MATRIX.
!     FILEB = MATRIX CONTROL BLOCK FOR ARBITRARY MATRIX.
!     FILEC = MATRIX CONTROL BLOCK FOR PRODUCT MATRIX.
!     Z     = ADDRESS OF A BLOCK OF CORE FOR WORKING SPACE. ZD IS SAME
!             BLOCK.
!     NZ    = LENGTH OF THIS BLOCK.
!     FLAG .EQ. 0 FOR PRE-MULTIPLICATION BY DIAGONAL.
!     FLAG .NE. 0 FOR POST-MULTIPLICATION BY DIAGONAL.
!     SIGN .EQ. +1 FOR POSITIVE PRODUCT.
!     SIGN .EQ. -1 FOR NEGATIVE PRODUCT.
 
 
 
 INTEGER, INTENT(IN OUT)                  :: z(1)
 DOUBLE PRECISION, INTENT(IN)             :: zd(1)
 INTEGER :: filea ,fileb ,filec ,flag  ,SIGN  ,sysbuf,eol   ,  &
     eor   ,TYPE  ,one   , rd    ,rdrew ,wrt   ,  &
     buf1  ,buf2  ,clsrew,rcc   ,ptype ,qtype ,wrtrew
 DOUBLE PRECISION :: ad    ,xd
 DIMENSION        filea(7)     ,fileb(7)     ,filec(7)
 COMMON /dmpyx /  filea ,fileb ,filec ,nz    ,flag  ,SIGN
 COMMON /names /  rd    ,rdrew ,wrt   ,wrtrew,clsrew
 COMMON /zntpkx/  ad (2),i     ,eol   ,eor
 COMMON /zblpkx/  xd (2),ix
 COMMON /unpakx/  TYPE  ,one   ,nx    ,incr
 COMMON /system/  sysbuf
 
 
!     PERFORM GENERAL INITIALIZATION
 
 buf1 = nz - sysbuf + 1
 buf2 = buf1 - sysbuf
 one  = 1
 incr = 1
 filec(2) = 0
 filec(6) = 0
 filec(7) = 0
 nx = filea(3)
 
!     COMPUTE TYPE OF C MATRIX.
!     RCC = 1 FOR REAL, = 2 FOR COMPLEX
!     QTYPE = 2 FOR RDP, = 4 FOR CDP
 
 rcc = 0
 IF (filea(5) > 2 .OR. fileb(5) > 2) rcc = 2
 qtype = rcc + 2
 IF (rcc == 0) rcc = 1
 TYPE  = qtype*SIGN
 ptype = filec(5)
 
!     OPEN PRODUCT MATRIX AND WRITE HEADER RECORD.
 
 CALL gopen (filec(1),z(buf1),wrtrew)
 
!     UNPACK DIAGONAL MATRIX IN CORE AND OPEN ARBITRARY MATRIX.
 
 CALL gopen (filea(1),z(buf2),rdrew)
 CALL unpack (*130,filea,z)
 CALL CLOSE (filea(1),clsrew)
 CALL gopen (fileb(1),z(buf2),rdrew)
 
!     PERFORM MATRIX MULTIPLICATION.
 
 j  = 1
 60 kr = (j-1)*rcc + 1
 CALL bldpk (qtype,ptype,filec(1),0,0)
 CALL intpk (*90,fileb(1),0,qtype,0)
 70 CALL zntpki
 kl = (i-1)*rcc + 1
 k  = kl
 IF (flag /= 0) k = kr
 xd(1) = zd(k)*ad(1)
 IF (rcc == 1) GO TO 80
 xd(1) = xd(1) - zd(k+1)*ad(2)
 xd(2) = zd(k)*ad(2) + zd(k+1)*ad(1)
 80 ix = i
 CALL zblpki
 IF (eol == 0) GO TO 70
 90 CALL bldpkn (filec(1),0,filec)
 j = j + 1
 IF (j <= fileb(2)) GO TO 60
 GO TO 140
 
!     CODE FOR NULL DIAGONAL MATRIX.
 
 130 CALL bldpkn (filec(1),0,filec)
 IF (filec(2) < fileb(2)) GO TO 130
 
!     CLOSE FILES AND RETURN.
 
 140 CALL CLOSE (filea(1),clsrew)
 CALL CLOSE (fileb(1),clsrew)
 CALL CLOSE (filec(1),clsrew)
 RETURN
END SUBROUTINE dmpy
