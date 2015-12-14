INTEGER FUNCTION eject (lines)
     
 
 INTEGER, INTENT(IN)                      :: lines
 COMMON /system/ skp1(8),maxlin,skp2(2),lincnt
 
!     LINES = NUNBER OF LINES TO BE PRINTED.
!     RESULT = 1 IF NEW PAGE IS STARTED.
 
 eject = 0
 IF (lincnt+lines+2 <= maxlin) GO TO 105
 CALL page1
 eject = 1
 105 RETURN
END FUNCTION eject
