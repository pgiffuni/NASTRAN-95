FUNCTION corwds (i,j)
     
 
 
 INTEGER, INTENT(IN OUT)                  :: i
 INTEGER, INTENT(IN OUT)                  :: j
 INTEGER :: corwds
 
 corwds = IABS(locfx(i) - locfx(j)) + 1
 RETURN
END FUNCTION corwds
