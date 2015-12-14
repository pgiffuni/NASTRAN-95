SUBROUTINE skprec(ifile,k)
     
 
 INTEGER, INTENT(IN OUT)                  :: ifile
 INTEGER, INTENT(IN)                      :: k
 INTEGER :: NAME(2)
 
 DATA NAME /4HSKPR,2HEC /
 
! ----------------------------------------------------------------------
 
 IF( k  < 0) THEN
   GO TO    10
 ELSE IF ( k  == 0) THEN
   GO TO    30
 ELSE
   GO TO    20
 END IF
 
 10 m=IABS(k)
 DO  i=1,m
   CALL bckrec(ifile)
 END DO
 GO TO 30
 
 20 DO  i=1,k
   CALL fwdrec(*40,ifile)
 END DO
 
 30 RETURN
 
 40 CALL mesage(-2,ifile,NAME)
 GO TO 30
 
END SUBROUTINE skprec
