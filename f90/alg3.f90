FUNCTION alg3 (p,h)
     
 
 
 REAL, INTENT(IN OUT)                     :: p
 REAL, INTENT(IN OUT)                     :: h
 COMMON /gas/ g,ej,r,cp,gamma,rojcp
 
 alg3=cp*ALOG(h/cp)-r/ej*ALOG(p)
 RETURN
END FUNCTION alg3
