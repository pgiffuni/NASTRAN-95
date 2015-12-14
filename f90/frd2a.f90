SUBROUTINE frd2a (nqhl,qhr,qhi,ih,nfreq)
     
 
 INTEGER, INTENT(IN)                      :: nqhl
 INTEGER, INTENT(IN OUT)                  :: qhr
 INTEGER, INTENT(IN)                      :: qhi
 INTEGER, INTENT(IN)                      :: ih
 INTEGER, INTENT(IN OUT)                  :: nfreq
 INTEGER :: sysbuf,mcb(7),thr(7),thi(7)
 DIMENSION       z(1)
 COMMON /zzzzzz/ z
 COMMON /system/ sysbuf
 COMMON /unpakx/ iout,inn,nnn,incr1
 COMMON /packx / iti,ito,ii,nn,incr
 
!     FIND COLUMN OF NQHL AND COPY REAL TO QHR AND IMAG TO QHI
 
 nz = korsz(z) - sysbuf
 mcb(1) = nqhl
 CALL rdtrl (mcb)
 IF (mcb(2) == 0) GO TO 999
 iout = mcb(5)
 iti  = 1
 IF (iout == 4) iti = 2
 ito  = iti
 nnn  = mcb(3)
 inn  = 1
 incr1= 1
 ii   = 1
 nn   = ih
 incr = 2
 nwc  = 2
 IF (iout == 4) nwc = 4
 ibuf1 = nz
 ibuf2 = ibuf1 - sysbuf
 CALL OPEN (*999,nqhl,z(ibuf1),0)
 CALL READ (*999,*999,nqhl,z(1),-2,1,flag)
 CALL makmcb (thr,qhr,ih,mcb(4),ito)
 CALL makmcb (thi,qhi,ih,mcb(4),ito)
 CALL skprec (nqhl,nfreq-1)
 CALL unpack (*25,nqhl,z(1))
 GO TO 30
 25 CALL zeroc  (z,nnn*nwc)
 30 j = 1
 CALL CLOSE (nqhl,1)
 CALL gopen (qhr,z(ibuf2),1)
 CALL gopen (qhi,z(ibuf1),1)
 DO  i = 1,ih
   CALL pack (z(j),qhr,thr)
   CALL pack (z(j+1),qhi,thi)
   j = j + ih*nwc
 END DO
 CALL CLOSE  (qhr,1)
 CALL CLOSE  (qhi,1)
 CALL wrttrl (thr)
 CALL wrttrl (thi)
 CALL dmpfil (-qhr,z,nz)
 CALL dmpfil (-qhi,z,nz)
 GO TO 1000
 999 CALL makmcb (thr,qhr,0,0,0)
 CALL wrttrl (thr)
 thr(1) = qhi
 CALL wrttrl (thr)
 1000 RETURN
END SUBROUTINE frd2a
