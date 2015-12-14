SUBROUTINE frd2e (in,io,nload,nfreq)
     
 
 INTEGER, INTENT(IN)                      :: in
 INTEGER, INTENT(IN)                      :: io
 INTEGER, INTENT(IN)                      :: nload
 INTEGER, INTENT(IN)                      :: nfreq
 INTEGER :: ma(7),mb(7)
 COMMON /system/ isys
 COMMON /unpakx/ iout,inn,nnn,incr1
 COMMON /zzzzzz/ z(1)
 
!     MAKE UHDF FROM IN
 
 incr1 = 1
 ma(1) = in
 mb(1) = io
 ib1   = korsz(z) - isys
 ib2   = ib1 - isys
 CALL rdtrl (ma)
 iout  = ma(5)
 CALL gopen (in,z(ib1),0)
 CALL gopen (io,z(ib2),1)
 CALL makmcb (mb,io,ma(3),ma(4),iout)
 DO  j = 1,nload
   CALL skprec (in,j-1)
   DO  i = 1,nfreq
     CALL cyct2b (in,io,1,z,mb)
     IF (i /= nfreq) CALL skprec (in,nload-1)
   END DO
   CALL REWIND (in)
   CALL skprec (in,1)
 END DO
 CALL CLOSE (in,1)
 CALL CLOSE (io,1)
 CALL wrttrl (mb)
 RETURN
END SUBROUTINE frd2e
