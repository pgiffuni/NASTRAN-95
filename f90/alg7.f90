FUNCTION alg7 (h,s)
 
    REAL, INTENT(IN)                         :: h
    REAL, INTENT(IN OUT)                     :: s
    COMMON /gas/ g,ej,r,cp,gamma,rojcp
 
    alg7=h/cp

    RETURN
END FUNCTION alg7
