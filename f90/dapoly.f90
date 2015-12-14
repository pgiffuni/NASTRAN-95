DOUBLE PRECISION FUNCTION dapoly(n,p)
     
!     CALCULATES AREA OF A POLYGON DESCRIBED BY N POINTS (P)
!        ( N .LE. 10 )
 
!        AREA= -1* LINE INTEGRAL OF Y*DX
 
!     AREA CONTRIBUTION FROM SIDE WHOSE ENDS ARE P(I), P(J):
!        A(I,J)= 0.5 * (Y(I)+Y(J)) * (X(I)-X(J))
 
 
 INTEGER, INTENT(IN)                      :: n
 DOUBLE PRECISION, INTENT(IN)             :: p(2,1)
 
 INTEGER :: kedge(2,10), k(2,10)
 
 DATA kedge/ 1,2,  2,3,  3,4,  4,5,  5,6,  6,7,  7,8,  8,9,  9,10, 10,1/
 
 DO  i=1,2
   DO  j=1,n
     k(i,j)= kedge(i,j)
   END DO
 END DO
 k(2,n)= 1
 dapoly= 0.0
 
 DO   nn= 1,n
   k1= k(1,nn)
   k2= k(2,nn)
   dapoly= dapoly +5.d-1 * (p(2,k1)+p(2,k2)) * (p(1,k1)-p(1,k2))
 END DO
 RETURN
END FUNCTION dapoly
