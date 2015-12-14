SUBROUTINE ploadx
     
!     PLOADX BUILDS THE PRESSURE LOADS FROM A PLOADX CARD FOR THE
!     TRIAX6 ELEMENT
 
 INTEGER :: slt,igd(1),islc(1),nam(2)
 DIMENSION       slc(5),gd(12),p(3),pn(3)
 COMMON /loadx / lc,slt
 COMMON /ssg1ax/ z(1)
 COMMON /condsa/ pi
 EQUIVALENCE     (slc(1),islc(1),p1),(slc(2),p3),(gd(1),igd(1))
 DATA    nam   / 4HPLOA,4HDX    /
 
 CALL READ (*30,*40,slt,slc,5,0,flag)
 j  = 1
 DO  i = 1,3
   CALL fndpnt (gd(j),islc(i+2))
   j  = j + 4
 END DO
 rl = gd(10) - gd(2)
 zl = gd(12) - gd(4)
 
!     LOADS IN NORMAL DIRECTION
 
 pn(1) = pi/30.*(9.0*gd( 2)*p1 + gd( 2)*p3 + gd(10)*p1 - gd(10)*p3)
 pn(2) = pi/7.5*(3.*(gd( 2)*p1 + gd(10)*p3) + 2.*(gd( 2)*p3 + gd(10)*p1))
 pn(3) = pi/30.*(9.0*gd(10)*p3 + gd( 2)*p3 + gd(10)*p1 - gd( 2)*p1)
 
 j = 1
 DO  i = 1,3
   p(1) =-zl*pn(i)
   p(2) = 0.0
   p(3) = rl*pn(i)
   
!     CONVERT TO GLOBAL IF NEEDED, AND INSERT INTO THE LOAD VECTOR
   
   IF (igd(j) /= 0) CALL basglb (p,p,gd(j+1),igd(j))
   CALL fndsil (islc(i+2))
   k = islc(i+2)
   z(k  ) = z(k  ) + p(1)
   z(k+1) = z(k+1) + p(2)
   z(k+2) = z(k+2) + p(3)
   j = j + 4
 END DO
 GO TO 60
 
!     ERROR MESSAGE
 
 30 j = -1
 GO TO 50
 40 j = -2
 50 CALL mesage (j,slt,nam)
 
 60 RETURN
END SUBROUTINE ploadx
