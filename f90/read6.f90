SUBROUTINE read6 (irig,gphia,nr,phia)
     
!     ADDS GIVENS EIGENVECTORS TO RIGID BODY MODES ON PHIA
 
 
 INTEGER, INTENT(IN)                      :: irig
 INTEGER, INTENT(IN)                      :: gphia
 INTEGER, INTENT(IN)                      :: nr
 INTEGER, INTENT(IN OUT)                  :: phia
 INTEGER :: sysbuf, mcb(7),FILE
 REAL :: z(3)
 COMMON /system/ sysbuf
 COMMON /zzzzzz/ iz(1)
 COMMON /unpakx/ it2u,iiu,jju,incr1u
 COMMON /packx / it1,it2,ii,jj,incr1
 EQUIVALENCE     (iz(1),z(1))
 
 
 ibuf1 =  korsz(z) - sysbuf + 1
 ibuf2 =  ibuf1 - sysbuf
 mcb(1)= gphia
 CALL rdtrl (mcb)
 ncol = mcb(2) - nr
 ii   = 1
 jj   = mcb(3)
 it1  = mcb(5)
 it2  = mcb(5)
 it2u = mcb(5)
 CALL makmcb (mcb,phia,jj,mcb(4),it1)
 incr1 = 1
 CALL gopen (phia,z(ibuf1),1)
 IF (nr == 0) GO TO 21
 FILE = irig
 CALL gopen (irig,z(ibuf2),0)
 z(1) = 0.0
 z(2) = 0.0
 DO  i = 1,nr
   iiu = 0
   CALL unpack (*11,irig,z(3))
   ii = iiu
   jj = jju
   CALL pack (z(3),phia,mcb)
   CYCLE
   11 ii = 1
   jj = 1
   CALL pack (z,phia,mcb)
 END DO
 CALL CLOSE (irig,1)
 21 CONTINUE
 IF (ncol <= 0) GO TO 31
 CALL gopen (gphia,z(ibuf2),0)
 FILE = gphia
 incr1u = 1
 z(1) = 0.0
 z(2) = 0.0
 CALL skprec (gphia,nr)
 DO  i = 1,ncol
   iiu = 0
   CALL unpack (*35,gphia,z(3))
   ii = iiu
   jj = jju
   CALL pack (z(3),phia,mcb)
   CYCLE
   35 ii = 1
   jj = 1
   CALL pack (z,phia,mcb)
 END DO
 CALL CLOSE (gphia,1)
 31 CALL CLOSE (phia,1)
 CALL wrttrl (mcb)
 RETURN
END SUBROUTINE read6
