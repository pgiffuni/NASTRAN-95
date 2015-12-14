FUNCTION apdf (f,in,ns)
     
 
 
 REAL, INTENT(IN)                         :: f(1)
 INTEGER, INTENT(IN OUT)                  :: in
 INTEGER, INTENT(IN OUT)                  :: ns
 
 
!     IF (NS .EQ. 0) GO TO 10
!     APDF = FLOAT(IN-1)/FLOAT(NS)
!     RETURN
!  10 APDF = F(IN)
!     RETURN
 
 apdf = f(in)
 IF (ns /= 0) apdf = FLOAT(in-1)/FLOAT(ns)
 RETURN
END FUNCTION apdf
