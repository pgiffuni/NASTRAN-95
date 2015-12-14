FUNCTION alg5 (h,s)
     
 
 
 REAL, INTENT(IN OUT)                     :: h
 REAL, INTENT(IN OUT)                     :: s
 COMMON /gas/ g,ej,r,cp,gamma,rojcp
 
 alg5=alg4(h,s)/(r*h)*cp
 RETURN
END FUNCTION alg5
