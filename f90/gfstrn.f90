SUBROUTINE gfstrn(a,at,i,scr1)
     
!     MATRIX TRANSPOSE ROUTINE
 
!     MORE EFFICIENT THEN TRANSPOSE FOR SPARSE MATRICES
 
 
!     TRANSPOSE IS SOLVED BY THE FOLLOWING EQUATION
 
!                                    T
!                  --    --   --   --  --   --
!                  I      I   I     I  I     I
!                  I  AT  I = I  A  I  I  I  I
!                  I      I   I     I  I     I
!                  --    --   --   --  --   --
 
!     WHERE I IS AN IDENITY MATRIX
 
 
 
 INTEGER, INTENT(IN)                      :: a
 INTEGER, INTENT(IN OUT)                  :: at
 INTEGER, INTENT(IN OUT)                  :: i
 INTEGER, INTENT(IN OUT)                  :: scr1
 INTEGER :: mcb(7) ,z        ,sysbuf   ,NAME(2)
 
!     SYSTEM PARAMETERS
 
 COMMON /system /        sysbuf
 
!     PACK COMMON
 
 COMMON / zblpkx /       val(4)   ,irow
 
!     OPEN CORE
 
 COMMON / zzzzzz /        z(1)
 
 DATA NAME / 4HGFST , 4HRN   /
 
!***********************************************************************
 
 nz = korsz(z)
 ibuf = nz - sysbuf
 IF(ibuf < 0) CALL mesage(-8,0,NAME)
 
!     GET MATRIX TRAILER
 
 mcb(1) = a
 CALL rdtrl(mcb)
 IF(mcb(1) < 0) RETURN
 ir = mcb(3)
 
!     GENERATE A SQUARE IDENITY MATRIX   IR BY IR
 
 val(1) = 1.0
 CALL makmcb(mcb,i,ir,2,2)
 CALL gopen(i,z(ibuf),1)
 
 DO  irow=1,ir
   CALL bldpk(1,2,i,0,0)
   CALL zblpki
   CALL bldpkn(i,0,mcb)
 END DO
 CALL CLOSE(i,1)
 CALL wrttrl(mcb)
 
!     PERFORM MULTIPLY
 
 CALL ssg2b(a,i,0,at,1,2,1,scr1)
 
 RETURN
END SUBROUTINE gfstrn
