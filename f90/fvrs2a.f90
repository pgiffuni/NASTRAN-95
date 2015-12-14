SUBROUTINE fvrs2a (FILE,kk1,kk2,noro,buffer)
     
!     GENERATE COLUMN REORDERING MATRIX. THIS MATRIX WILL REORDER
!     COLUMNS OF A MATRIX BY POST-MULTIPLYING THE MATRIX WHOSE
!     COLUMNS ARE TO BE REORDERED BY THE REORDERING MATRIX.
 
!     THE MATRIX WILL BE A REAL SINGLE-PRECISION SQUARE MATRIX.
 
 
 INTEGER, INTENT(IN)                      :: FILE
 INTEGER, INTENT(IN)                      :: kk1
 INTEGER, INTENT(IN)                      :: kk2
 INTEGER, INTENT(OUT)                     :: noro
 INTEGER, INTENT(IN OUT)                  :: buffer(1)
 INTEGER :: row,trl(7), typin,typout
 
 COMMON /packx/ typin,typout,ii,nn,incr
 
 noro = -1
 IF(kk1 == 1 .OR. kk2 == 1) RETURN
 
 noro = 1
 
 typin  = 1
 typout = 1
 incr   = 1
 
 trl(1) = FILE
 trl(2) = 0
 trl(3) = kk1*kk2
 trl(4) = 1
 trl(5) = typout
 trl(6) = 0
 trl(7) = 0
 
 CALL gopen(FILE,buffer,1)
 
 value = 1.0
 
 DO  k1 = 1,kk1
   row = k1
   DO  k2 = 1,kk2
     
     ii = row
     nn = row
     CALL pack(value,FILE,trl)
     
     row = row + kk1
     
   END DO
 END DO
 
 CALL CLOSE(FILE,1)
 CALL wrttrl(trl)
 
 RETURN
END SUBROUTINE fvrs2a
