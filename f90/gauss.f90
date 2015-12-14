SUBROUTINE gauss (a,n,n2)
     
 
 
 COMPLEX, INTENT(OUT)                     :: a(20,1)
 INTEGER, INTENT(IN)                      :: n
 INTEGER, INTENT(IN)                      :: n2
 
 
 DO  i=1,n
   k=i+1
   DO  j=k,n2
     a(i,j)=a(i,j)/a(i,i)
   END DO
   DO  m=1,n
     IF(m == i) CYCLE
     DO  l=k,n2
       a(m,l)=a(m,l)-a(m,i)*a(i,l)
     END DO
   END DO
 END DO
 RETURN
END SUBROUTINE gauss
