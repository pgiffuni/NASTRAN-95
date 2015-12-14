SUBROUTINE eqscod (loc,n,z)
     
 
 INTEGER, INTENT(IN)                      :: loc
 INTEGER, INTENT(IN)                      :: n
 INTEGER, INTENT(IN OUT)                  :: z(1)
 EXTERNAL  lshift,orf
 INTEGER :: orf
 
 i   = loc
 mend= loc+n-1
 1 ist = i
 ng  = 1
 2 CONTINUE
 IF (i >= mend-2) GO TO 3
 IF (z(i+3) /= z(ist)) GO TO 3
 ng = ng+1
 i  = i+3
 GO TO 2
 3 CONTINUE
 IF (ng /= 1) GO TO 4
 i = i+3
 IF (i >= mend-2) GO TO 6
 GO TO 1
 4 DO  j=1,ng
   iloc  = ist+3*(j-1)
   icode = 8*j+ng
   inew  = lshift(icode,26)
   z(iloc+2) = orf(z(iloc+2),inew)
 END DO
 i = i+3
 IF (i >= mend-2) GO TO 6
 GO TO 1
 6 CONTINUE
 RETURN
END SUBROUTINE eqscod
