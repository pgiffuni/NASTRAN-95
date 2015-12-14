SUBROUTINE alg1 (lnct)
 
    INTEGER, INTENT(OUT)                     :: lnct
    DIMENSION       rdata(4)
    COMMON /gas   / g,ej,r,cp,gamma,rojcp
    COMMON /system/ sysbuf,nout
    COMMON /algino/ dum,nalgdb
    COMMON /ud3prt/ iprtc
 
    log1 = nalgdb
    log2 = nout
 
    CALL fread (log1,rdata,4,1)
 
    cp   = rdata(1)
    r    = rdata(2)
    g    = rdata(3)
    ej   = rdata(4)
 
    IF (cp == 0.0) cp = 0.24
    IF (r  == 0.0) r  = 53.32
    IF (g  == 0.0) g  = 32.174
    IF (ej == 0.0) ej = 778.16
    IF (iprtc == 1) WRITE(log2,10) cp,r,g,ej
 
10  FORMAT (/10X,'SPECIFIC HEAT AT CONSTANT PRESSURE',5X,1H=,f8.5,  &
        /10X,'GAS CONSTANT',27X,1H=,f8.4,  &
        /10X,'GRAVITATIONAL CONSTANT',17X,1H=,f8.4,  &
        /10X,'JOULES EQUIVALENT',22X,1H=,f8.3)
             
    lnct  = lnct + 5
    rojcp = r/(ej*cp)
    gamma = 1.0/(1.0-rojcp)
 
    RETURN
END SUBROUTINE alg1
