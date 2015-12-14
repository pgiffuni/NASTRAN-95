SUBROUTINE tabpt
     
!     MODULE DRIVER TO PRINT TABLES
 
 DIMENSION      in(5),itrl(7)
 COMMON /BLANK/ op(2),irc,iwd
 DATA    in   / 101,102,103,104,105 /, BLANK / 4H     /
 
 DO  i = 1,5
   itrl(1) = in(i)
   CALL rdtrl (itrl(1))
   IF (itrl(1) > 0) CALL tabprt (in(i))
 END DO
 op(1) = BLANK
 op(2) = BLANK
 irc = 0
 iwd = 0
 RETURN
END SUBROUTINE tabpt
