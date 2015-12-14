SUBROUTINE gauss2 (a,n,n2)
     
 
 COMPLEX, INTENT(IN OUT)                  :: a(20,1)
 INTEGER, INTENT(IN)                      :: n
 INTEGER, INTENT(IN)                      :: n2
 
 DOUBLE COMPLEX da(20,30)
 
 DO  i = 1, n
   DO  j = 1, n2
     da(i,j) = a(i,j)
   END DO
 END DO
 DO  i=1,n
   k=i+1
   DO  j=k,n2
     da(i,j)=da(i,j)/da(i,i)
   END DO
   DO  m=1,n
     IF(m == i) CYCLE
     DO  l=k,n2
       da(m,l)=da(m,l)-da(m,i)*da(i,l)
     END DO
   END DO
 END DO
 DO  i = 1, n
   DO  j = 1, n2
     a(i,j) = da(i,j)
   END DO
 END DO
 RETURN
END SUBROUTINE gauss2
