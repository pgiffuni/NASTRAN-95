SUBROUTINE frd2d (in,io,ip)
     
 
 INTEGER, INTENT(IN)                      :: in
 INTEGER, INTENT(IN)                      :: io
 INTEGER, INTENT(IN OUT)                  :: ip
 INTEGER :: sysbuf,ma(7),mb(7)
 COMMON /system/ sysbuf
 COMMON /unpakx/ iout,inn,mnn,incr1
 COMMON /zzzzzz/ z(1)
 
!     ADD IN TO END OF IO
 
 incr1 = 1
 ma(1) = in
 mb(1) = io
 CALL rdtrl (ma)
 iout  = ma(5)
 nc    = korsz(z)
 ib1   = nc  - sysbuf
 ib2   = ib1 - sysbuf
 CALL gopen (in,z(ib1),0)
 IF (ip /= 0) GO TO 10
 CALL gopen (io,z(ib2),1)
 CALL makmcb (mb,io,ma(3),2,iout)
 GO TO 20
 10 CALL gopen (io,z(ib2),3)
 CALL rdtrl (mb)
 20 n = ma(2)
 CALL cyct2b (in,io,n,z,mb)
 CALL CLOSE  (in,1)
 CALL CLOSE  (io,3)
 CALL wrttrl (mb)
 CALL dmpfil (-in,z,nc)
 RETURN
END SUBROUTINE frd2d
