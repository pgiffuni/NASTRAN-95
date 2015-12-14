FUNCTION apdf (f,in,ns)
 
    REAL, INTENT(IN)                         :: f(1)
    INTEGER, INTENT(IN OUT)                  :: in
    INTEGER, INTENT(IN OUT)                  :: ns
 
    apdf = f(in)
    IF (ns /= 0) apdf = FLOAT(in-1)/FLOAT(ns)

    RETURN
END FUNCTION apdf
