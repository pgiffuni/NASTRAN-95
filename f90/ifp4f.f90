SUBROUTINE ifp4f (ibit,FILE,bit)
     
!     TEST BIT -IBIT- IN TRAILER OF DATA BLOCK -FILE-
 
 
 INTEGER, INTENT(IN)                      :: ibit
 INTEGER, INTENT(IN)                      :: FILE
 LOGICAL, INTENT(OUT)                     :: bit
 EXTERNAL     andf
 
 INTEGER :: two, trail(7), andf
 COMMON /two/ two(32)
 
 trail(1) = FILE
 CALL rdtrl (trail)
 i1 = (ibit-1)/16 + 2
 i2 = ibit - (i1-2)*16 + 16
 IF (andf(trail(i1),two(i2)) == 0.0) THEN
   GO TO    20
 END IF
 10 bit = .true.
 RETURN
 20 bit = .false.
 RETURN
END SUBROUTINE ifp4f
