FUNCTION alg4 (h,s)
     
 
 
 REAL, INTENT(IN OUT)                     :: h
 REAL, INTENT(IN OUT)                     :: s
 COMMON /gas/ g,ej,r,cp,gamma,rojcp
 
 alg4=EXP(ALOG(h/cp)/rojcp-ej/r*s)
 RETURN
END FUNCTION alg4
