SUBROUTINE ema1d( j, nsca, scalas, pivot, dict, cgv, kgg, cp, f )
!     SUBROUTINE EMA1S( J, NSCA, SCALAS, PIVOT, DICT, CGV, KGG, CP, F )
!*******
! EMA1S ADDS A COLUMN VECTOR IN REAL SINGLE PRECISION
! EMA1D ADDS A COLUMN VECTOR IN REAL DOUBLE PRECISION
 
!     J        INDEX IN SCALAS TO CURRENT RELATIVE COLUMN NBR
!     NSCA     NBR OF ROWS( TERMS ) PER GRID POINT IN COLUMN VECTOR
!     SCALAS   ARRAY OF RELATIVE ROW/COLUMN NUMBERS
!     PIVOT    6-WORD ARRAY AS FOLLOWS...
!          (1) INTERNAL INDEX OF PIVOT
!          (2) DOF OF PIVOT
!          (3) DOF OF EACH POINT CONNECTED TO PIVOT
!          (4) NBR OF CONNECTED POINTS
!          (5) INTERNAL INDEX OF  1ST CONNECTED POINT
!          (6) INTERNAL INDEX OF LAST CONNECTED POINT
!     DICT     DICTIONARY ENTRY FOR ELEMENT AS FOLLOWS...
!          (1) ELEMENT ID
!          (2) FORM( 1=RECT, 2=DIAG )
!          (3) NBR OF TERMS PER COLUMN
!          (4) COMPONENT CODE( DECODED IN SCALAS ARRAY )
!          (5) GE
!          (6) INTERNAL INDEX OF 1ST POINT
!          (7) GINO ADDR OF 1ST COLUMN PARTITION
!         ....
!     CGV      CONNECTED GRID POINT VECTOR
!     KGG      ADDR OF KGG COLUMNS FOR PIVOT
!     CP       ADDR OF COLUMN PARTITION
!     F        FACTOR( RSP ) TO BE APPLIED TO EACH TERM IN CP
 
!******
 
 INTEGER, INTENT(IN OUT)                  :: j
 INTEGER, INTENT(IN)                      :: nsca
 INTEGER, INTENT(IN)                      :: scalas(1)
 INTEGER, INTENT(IN)                      :: pivot(6)
 INTEGER, INTENT(IN)                      :: dict(7)
 INTEGER, INTENT(IN)                      :: cgv(1)
 DOUBLE PRECISION, INTENT(OUT)            :: kgg(1)
 DOUBLE PRECISION, INTENT(IN)             :: cp(1)
 REAL, INTENT(IN)                         :: f
 
 
!     REAL             KGG(1), CP(1)
 
 
! INITIALIZE
 
 icol0 = scalas(j)*pivot(3)*pivot(4)
 ii0   = pivot(5) - 1
 l     = 1
 IF( dict(2) /= 2 ) GO TO 20
 
! PROCESS DIAGONAL PARTITION
 
 ii = pivot(1)
 imat = icol0 + cgv(ii-ii0)+ scalas(j)
 kgg(imat) = kgg(imat) + f*cp(1)
 RETURN
 
! PROCESS RECTANGULAR PARTITION
 
 20 CONTINUE
 ngrid = 4 + 2*dict(3)/nsca
 DO  i=6,ngrid,2
   k = dict(i)
   IF( k == 0 ) RETURN
   imat = icol0 + cgv(k-ii0)
   DO  k=1,nsca
     m = scalas(k)
     kgg(imat+m) = kgg(imat+m) + f*cp(l)
     l = l + 1
   END DO
 END DO
 RETURN
END SUBROUTINE ema1d
