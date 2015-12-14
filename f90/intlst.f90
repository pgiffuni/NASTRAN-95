SUBROUTINE intlst (list,n,SIGN,n1,n2)
     
 
 INTEGER, INTENT(IN)                      :: list(1)
 INTEGER, INTENT(OUT)                     :: n
 INTEGER, INTENT(OUT)                     :: SIGN
 INTEGER, INTENT(OUT)                     :: n1
 INTEGER, INTENT(OUT)                     :: n2
 INTEGER :: TO,thru
 DATA    TO,thru/ 2HTO,4HTHRU /
 
 SIGN = ISIGN(1,list(n))
 n1 = IABS(list(n))
 IF (list(n+1) == TO .OR. list(n+1) == thru) GO TO 110
 n2 = n1
 n  = n + 1
 GO TO 150
 
 110 n2 = IABS(list(n+2))
 n  = n + 3
 IF (n1 <= n2) GO TO 150
 i  = n1
 n1 = n2
 n2 = i
 
 150 RETURN
END SUBROUTINE intlst
