FUNCTION area (g,i,j,k)
     
!     THIS ROUTINE IS CALLED BY SFAREA WHICH IS CALLED BY EMGFIN TO
!     COMPUTE THE SURFACE AREAS OF THE SOLID ELEMENTS
 
 
 
 REAL, INTENT(IN OUT)                     :: g(1)
 INTEGER, INTENT(IN OUT)                  :: i
 INTEGER, INTENT(IN OUT)                  :: j
 INTEGER, INTENT(IN OUT)                  :: k
 
 area = 0.5*SQRT(  &
     ((g(j+2)-g(i+2))*(g(k+3)-g(i+3))-(g(j+3)-g(i+3))*(g(k+2)-g(i+2))) **2  &
     +((g(j+3)-g(i+3))*(g(k+1)-g(i+1))-(g(j+1)-g(i+1))*(g(k+3)-g(i+3))) **2  &
     +((g(j+1)-g(i+1))*(g(k+2)-g(i+2))-(g(j+2)-g(i+2))*(g(k+1)-g(i+1))) **2)
 RETURN
END FUNCTION area
