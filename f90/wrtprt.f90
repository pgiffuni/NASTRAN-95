SUBROUTINE wrtprt (FILE,list,FORMAT,n)
     
 
 
 INTEGER, INTENT(IN OUT)                  :: FILE
 INTEGER, INTENT(IN OUT)                  :: list(1)
 INTEGER, INTENT(IN OUT)                  :: FORMAT(n)
 INTEGER, INTENT(IN OUT)                  :: n
 
 
 CALL WRITE (FILE,list,list(1)+1,0)
 CALL WRITE (FILE,n,1,0)
 CALL WRITE (FILE,FORMAT,n,0)
 RETURN
END SUBROUTINE wrtprt
