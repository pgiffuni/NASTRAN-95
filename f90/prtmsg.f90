SUBROUTINE prtmsg
     
 COMMON /zzzzzz/ buf(1)
 COMMON /output/ title(32,6)
 
 DATA inprew,msg,BLANK / 0,101,4H    /
 
 CALL OPEN (*110,msg,buf,inprew)
 CALL READ (*110,*110,msg,0,0,1,j)
 DO  j = 4,6
   DO  i = 1,32
     title(i,j) = BLANK
   END DO
 END DO
 CALL wrtmsg (msg)
 110 RETURN
END SUBROUTINE prtmsg
