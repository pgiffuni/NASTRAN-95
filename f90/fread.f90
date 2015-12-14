SUBROUTINE fread (FILE,BLOCK,n,eor)
     
 
 INTEGER, INTENT(IN OUT)                  :: FILE
 REAL, INTENT(IN OUT)                     :: BLOCK(1)
 INTEGER, INTENT(IN OUT)                  :: n
 INTEGER, INTENT(IN OUT)                  :: eor
 
 REAL :: subnam(2)
 DATA subnam / 4H fre,4HAD  /
 
 CALL READ (*100,*101,FILE,BLOCK,n,eor,k)
 RETURN
 100 k = -2
 GO TO 110
 101 k = -3
 110 CALL mesage (k,FILE,subnam)
 GO TO 110
END SUBROUTINE fread
