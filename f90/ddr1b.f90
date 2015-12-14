SUBROUTINE ddr1b (in1,in2,iout)
     
!     THIS ROUTINE REPLACES DISPLACEMNTS ON IN1 WITH DISPLACEMENTS ON
!     IN2  AND WRITES ON  IOUT
 
 
 INTEGER, INTENT(IN)                      :: in1
 INTEGER, INTENT(IN OUT)                  :: in2
 INTEGER, INTENT(IN)                      :: iout
 INTEGER :: sysbuf, mcb(7)
 COMMON /system/ sysbuf
 COMMON /unpakx/ itc,ii,jj,incr
 COMMON /zzzzzz/ z(1)
 
 
 nz = korsz(z) - sysbuf
 CALL gopen (in1,z(nz+1),0)
 nz = nz - sysbuf
 CALL gopen (in2,z(nz+1),0)
 nz = nz - sysbuf
 CALL gopen (iout,z(nz+1),1)
 mcb(1) = in1
 CALL rdtrl (mcb)
 mcb(1) = iout
 nd  = mcb(2)/3
 itc = mcb(5)
 incr = 1
 mcb(2) = 0
 mcb(6) = 0
 mcb(7) = 0
 DO  i = 1,nd
   CALL skprec (in1,1)
   CALL cyct2b (in2,iout,1,z,mcb)
   CALL cyct2b (in1,iout,2,z,mcb)
 END DO
 CALL wrttrl (mcb)
 CALL CLOSE (in1,1)
 CALL CLOSE (in2,1)
 CALL CLOSE (iout,1)
 RETURN
END SUBROUTINE ddr1b
