FUNCTION alg9 (h,s,v2)
 
    REAL, INTENT(IN OUT)                     :: h
    REAL, INTENT(IN OUT)                     :: s
    REAL, INTENT(IN)                         :: v2
    COMMON /gas/ g,ej,r,cp,gamma,rojcp

    alg9=cp*v2/(gamma*g*r*h)
 
    RETURN
END FUNCTION alg9
