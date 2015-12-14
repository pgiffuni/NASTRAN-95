FUNCTION alg2 (s,p)
     
 
 
 REAL, INTENT(IN OUT)                     :: s
 REAL, INTENT(IN OUT)                     :: p
 COMMON /gas/ g,ej,r,cp,gamma,rojcp
 
 alg2=cp*EXP(s/cp+rojcp*ALOG(p))
 RETURN
END FUNCTION alg2
