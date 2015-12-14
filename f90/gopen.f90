SUBROUTINE gopen (FILE,buffer,option)
     
 
 INTEGER, INTENT(IN OUT)                  :: FILE
 REAL, INTENT(IN OUT)                     :: buffer(1)
 INTEGER, INTENT(IN)                      :: option
 INTEGER :: ERR,outrew,outnor
 REAL :: subnam(2),header(2)
 DATA    subnam / 4H gop,4HEN  /
 DATA    outrew,inpnor,outnor  / 1,2,3 /
 
 CALL OPEN (*200,FILE,buffer,option)
 IF (option == inpnor .OR. option == outnor) GO TO 150
 IF (option == outrew)  GO TO 110
 CALL READ (*201,*202,FILE,header,2,1,ERR)
 GO TO 150
 110 CALL fname (FILE,header)
 CALL WRITE (FILE,header,2,1)
 150 RETURN
 
 200 ERR = -1
 GO TO 210
 201 ERR = -2
 GO TO 210
 202 ERR = -3
 210 CALL mesage (ERR,FILE,subnam)
 
 RETURN
END SUBROUTINE gopen
