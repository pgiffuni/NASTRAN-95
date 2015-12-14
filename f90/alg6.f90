FUNCTION alg6 (p,t)
     
 
 
 REAL, INTENT(IN OUT)                     :: p
 REAL, INTENT(IN)                         :: t
 COMMON /gas/ g,ej,r,cp,gamma,rojcp
 
 alg6=cp*t
 RETURN
END FUNCTION alg6
