SUBROUTINE alg03 (lnct,l)
     
 
    INTEGER, INTENT(OUT)                     :: lnct
    INTEGER, INTENT(IN)                      :: l
    COMMON /upage / limit,lq
    COMMON /ud3prt/ iprtc
 
    lnct=lnct+l
    IF(lnct <= limit)RETURN
    lnct=1+l
    IF (iprtc /= 0) WRITE(lq,100)
100 FORMAT(1H1)

    RETURN
END SUBROUTINE alg03
