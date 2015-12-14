SUBROUTINE shxtrs (nrow,ncol,array)
     
!     TO EXTRAPOLATE VALUES IN ARRAY FROM A SET OF EVALUATION POINTS TO
!     THE GRID POINTS OF SPECIFIC SHELL ELEMENTS.
!     THE EXTRAPOLATION IS IN TWO DIMENSIONS.
 
!     INPUT :
!           NROW   - SIZE OF THE SET OF VALUES
!           NCOL   - NUMBER OF EVALUATION POINTS
!           ARRAY  - ARRAY OF DATA TO BE EXTRAPOLATED
 
!     OUTPUT:
!           ARRAY  - ARRAY OF EXTRAPOLATED DATA
 
 
 
 INTEGER, INTENT(IN)                      :: nrow
 INTEGER, INTENT(IN)                      :: ncol
 REAL, INTENT(IN OUT)                     :: array(nrow,1)
 LOGICAL :: tria
 REAL :: temp(4,4),shp(4),tpoint(2,3),qpoint(2,4), xsi,eta
 
 
 
 DATA    tpoint  / 0.0,  0.0,  1.0,  0.0,  0.0,  1.0             /
 DATA    qpoint  /-1.0, -1.0,  1.0, -1.0,  1.0,  1.0, -1.0,  1.0 /
 
!     INITIALIZE
 
 tria   = ncol == 3
 
 DO  i = 1,nrow
   DO  j = 1,ncol
     temp(i,j) = 0.0
   END DO
 END DO
 
!     BEGIN LOOP ON DESTINATION POINTS
 
 DO  i = 1,ncol
   
!     EVALUATE PSEUDO-SHAPE FUNCTIONS
   
   IF (.NOT.tria) GO TO 30
   
   xsi    = tpoint(1,i)
   eta    = tpoint(2,i)
   shp(1) = 1.66666667 - 2.0*xsi - 2.0*eta
   shp(2) = 2.0*xsi - 0.33333333
   shp(3) = 2.0*eta - 0.33333333
   GO TO 40
   
   30 xsi    = qpoint(1,i)
   eta    = qpoint(2,i)
   const  = 0.577350269
   shp(1) = 0.75*(const-xsi)*(const-eta)
   shp(2) = 0.75*(const-xsi)*(const+eta)
   shp(3) = 0.75*(const+xsi)*(const-eta)
   shp(4) = 0.75*(const+xsi)*(const+eta)
   
!     EXTRAPOLATE
   
   40 DO  j = 1,nrow
     DO  k = 1,ncol
       temp(j,i) = temp(j,i) + shp(k)*array(j,k)
     END DO
   END DO
   
 END DO
 
!     COPY THE EXTRAPOLATED DATA BACK INTO ARRAY
 
 DO  j = 1,ncol
   DO  i = 1,nrow
     array(i,j) = temp(i,j)
   END DO
 END DO
 
 RETURN
END SUBROUTINE shxtrs
