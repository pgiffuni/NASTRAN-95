COMPLEX FUNCTION traile (x,j,n,p,m,boxl)
     
!     ROUTINE TO FIND PHI FOR TRAILING EDGE
 
 
 REAL, INTENT(IN)                         :: x
 INTEGER, INTENT(IN OUT)                  :: j
 INTEGER, INTENT(IN OUT)                  :: n(1)
 COMPLEX, INTENT(IN)                      :: p(3,m)
 INTEGER, INTENT(IN OUT)                  :: m
 REAL, INTENT(IN)                         :: boxl
 
 
 
!     CHECK TO SEE IF TRAILING EDGE HAS BEEN COMPUTED
 
 IF (n(j) >= 0) GO TO 300
 200  traile = p(1,j)
 RETURN
 
 300  xa = x/boxl + 0.5 - FLOAT(n(j))
 IF (REAL(p(2,j)) == 0.0) GO TO 200
 traile = p(1,j) + xa*(p(1,j) - p(2,j))
 RETURN
END FUNCTION traile
