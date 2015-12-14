SUBROUTINE dmatrs (d,v,c,ca, ca2, va, dm, db, yi)
     
 
! THIS ROUTINE COMPUTES THE STIFFNESS MATRIX IN FIELD COORDINATES FOR
! THE TOROIDAL RING
 
 
! NOTE THE DOUBLE SUBSCRIPTING USED IN DMATRIX SUBROUTINE IS
! COMPATIBLE WITH THE CALLING PROGRAM. THE DELINT ARRAY OF INTEGRALS
! IS A (11X6) SINGLY SUBSCRIPTED ARRAY (STORED ROWWISE) IN THE CALLING
! PROGRAM AND IT IS A (6X11) DOUBLY SUBSCRIPTED ARRAY (STORED
! COLUMNWISE) IN DMATRX ROUTINE.
 
 
 
 
 REAL, INTENT(OUT)                        :: d(10,10)
 REAL, INTENT(IN)                         :: v
 REAL, INTENT(IN)                         :: c
 REAL, INTENT(IN)                         :: ca
 REAL, INTENT(IN OUT)                     :: ca2
 REAL, INTENT(IN)                         :: va
 REAL, INTENT(IN)                         :: dm
 REAL, INTENT(IN)                         :: db
 REAL, INTENT(IN)                         :: yi(6,11)
 
 
!     ------------------------------------------------------------------
 
 d(1,1)   =  dm * (ca2*yi(1,1) + 2.*va*yi(2,1) + yi(3,1))
 d(2,1)   =  dm * (ca2*yi(1,2) + 2.*va*yi(2,2) + yi(3,2))
 d(3,1)   =  dm * (ca2*yi(1,3) + 2.*va*yi(2,3) + yi(3,3))
 d(4,1)   =  dm * (ca2*yi(1,4) + 2.*va*yi(2,4) + yi(3,4))
 d(5,1)   =  dm * (ca2*yi(1,5) + 2.*va*yi(2,5) + yi(3,5))
 d(6,1)   =  dm * (ca2*yi(1,6) + 2.*va*yi(2,6) + yi(3,6))
 d(7,1)   =  dm * (va*yi(4,1)  +  yi(5,1))
 d(8,1)   =  dm * (ca*yi(1,1)  +  va*yi(4,2) + v*yi(2,1) + yi(5,2))
 d(9,1)   =  dm * (2.*ca*yi(1,2) + va*yi(4,3) + 2.*v*yi(2,2) + yi(5 ,3))
 d(10,1)  =  dm * (3.*ca*yi(1,3) + va*yi(4,4) + 3.*v*yi(2,3) + yi(5 ,4))
 d(2,2)   =  db * yi(6,1)  +  d(3,1)
 d(3,2)   =  db * (2.*v*yi(4,1) + 2.*yi(6,2))  +  d(4,1)
 d(4,2)   =  db * (6.*v*yi(4,2) + 3.*yi(6,3))  +  d(5,1)
 d(5,2)   =  db * (12.*v*yi(4,3) + 4.*yi(6,4)) +  d(6,1)
 d(6,2)   =  dm * (ca2*yi(1,7) + 2.*va*yi(2,7) + yi(3,7))  +  &
     db * (20.*v*yi(4,4) + 5.*yi(6,5))
 d (7,2)  =  dm * (va*yi(4,2) + yi(5,2))
 d(8,2)   =  dm * (ca*yi(1,2) + va*yi(4,3) + v*yi(2,2) + yi(5,3))
 d(9,2)   =  dm * (2.*ca*yi(1,3) + va*yi(4,4) + 2.*v*yi(2,3) + yi(5 ,4))
 d(10,2)  =  dm * (3.*ca*yi(1,4) + va*yi(4,5) + 3.*v*yi(2,4) + yi(5 ,5))
 d(3,3)   =  db * 4.*(c*yi(1,1) + 2.*v*yi(4,2) + yi(6,3)) + d(5,1)
 d(4,3)   =  db * 6.*(2.*c*yi(1,2) + 3.*v*yi(4,3) + yi(6,4))+d(6,1)
 d(5,3)   =  dm * (ca2*yi(1,7) + 2.*va*yi(2,7) + yi(3,7))  +  &
     db * 2.*(12.*c*yi(1,3)+ 16.*v*yi(4,4)+ 4.*yi(6,5))
 d(6,3)   =  dm * (ca2*yi(1,8) + 2.*va*yi(2,8) + yi(3,8))  +  &
     db * 10.*(4.*c*yi(1,4) + 5.*v*yi(4,5) + yi(6,6))
 d (7,3)  =  dm * (va*yi(4,3) + yi(5,3))
 d(8,3)   =  dm * (ca*yi(1,3) + va*yi(4,4) + v*yi(2,3) + yi(5,4))
 d(9,3)   =  dm * (2.*ca*yi(1,4) + va*yi(4,5) + 2.*v*yi(2,4) + yi(5 ,5))
 d(10,3)  =  dm * (3.*ca*yi(1,5) + va*yi(4,6) + 3.*v*yi(2,5) + yi(5 ,6))
 d(4,4)   =  dm * (ca2*yi(1,7) + 2.*va*yi(2,7) + yi(3,7))  +  &
     db * 9.*(4.*c*yi(1,3) + 4.*v*yi(4,4) + yi(6,5))
 d(5,4)   =  dm * (ca2*yi(1,8) + 2.*va*yi(2,8) + yi(3,8))  +  &
     db * 12.*(6.*c*yi(1,4) + 5.*v*yi(4,5) + yi(6,6))
 d(6,4)   =  dm * (ca2*yi(1,9) + 2.*va*yi(2,9) + yi(3,9))  +  &
     db * 15.*(8.*c*yi(1,5) + 6.*v*yi(4,6) + yi(6,7))
 d (7,4)  =  dm * (va*yi(4,4) + yi(5,4))
 d(8,4)   =  dm * (ca*yi(1,4) + va*yi(4,5) + v*yi(2,4) + yi(5,5))
 d(9,4)   =  dm * (2.*ca*yi(1,5) + va*yi(4,6) + 2.*v*yi(2,5) + yi(5 ,6))
 d(10,4)  =  dm * (3.*ca*yi(1,6) + va*yi(4,7) + 3.*v*yi(2,6) + yi(5,7))
 d(5,5)   =  dm * (ca2*yi(1,9) + 2.*va*yi(2,9) + yi(3,9))  +  &
     db * 16.*(9.*c*yi(1,5) + 6.*v*yi(4,6) + yi(6,7))
 d(6,5)   =  dm * (ca2*yi(1,10) + 2.*va*yi(2,10) + yi(3,10)) +  &
     db * 20.*(12.*c*yi(1,6) + 7.*v*yi(4,7)  + yi(6,8))
 d (7,5)  =  dm * (va*yi(4,5) + yi(5,5))
 d(8,5) = dm * (ca*yi(1,5) + va*yi(4,6) + v*yi(2,5) + yi(5,6))
 d(9,5)   =  dm * (2.*ca*yi(1,6) + va*yi(4,7) + 2.*v*yi(2,6) + yi(5 ,7))
 d(10,5)  =  dm * (3.*ca*yi(1,7) + va*yi(4,8) + 3.*v*yi(2,7) + yi(5 ,8))
 d(6,6)   =  dm * (ca2*yi(1,11) + 2.*va*yi(2,11) + yi(3,11))  +  &
     db * 25.*(16.*c*yi(1,7) + 8.*v*yi(4,8) + yi(6,9))
 d (7,6)  =  dm * (va*yi(4,6) + yi(5,6))
 d(8,6)   =  dm * (ca*yi(1,6) + va*yi(4,7) + v*yi(2,6) + yi(5,7))
 d(9,6)   =  dm * (2.*ca*yi(1,7) + va*yi(4,8) + 2.*v*yi(2,7) + yi(5 ,8))
 d(10,6)  =  dm * (3.*ca*yi(1,8) + va*yi(4,9) + 3.*v*yi(2,8) + yi(5 ,9))
 d (7,7)  =  dm * yi(6,1)
 d (8,7)  =  dm * (v*yi(4,1) + yi(6,2))
 d (9,7)  =  dm * (2.*v*yi(4,2) + yi(6,3))
 d (10,7) =  dm * (3.*v*yi(4,3) + yi(6,4))
 d(8,8)   =  dm * (c*yi(1,1) + 2.*v*yi(4,2) + yi(6,3))
 d(9,8)   =  dm * (2.*c*yi(1,2) + 3.*v*yi(4,3) + yi(6,4))
 d(10,8)  =  dm * (3.*c*yi(1,3) + 4.*v *yi(4,4)+ yi(6,5))
 d(9,9)   =  dm * (4.*c*yi(1,3) + 4.*v*yi(4,4) + yi(6,5))
 d(10,9)  =  dm * (6.*c*yi(1,4) + 5.*v*yi(4,5) + yi(6,6))
 d(10,10) =  dm * (9.*c*yi(1,5) + 6.*v*yi(4,6) + yi(6,7))
 DO  i=1,10
   DO  j=1,i
     d(j,i)  = d(i,j)
   END DO
 END DO
 RETURN
END SUBROUTINE dmatrs
