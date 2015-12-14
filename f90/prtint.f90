SUBROUTINE prtint
     
 INTEGER :: trlr
 INTEGER :: opt, prt
 
!     OPT = 0 IF MATRIX BY COLUMNS...1 IF BY ROWS.
 
 REAL :: NAME(2)
 
 COMMON /BLANK/ opt, prt
 COMMON /xxmprt/ trlr(7)
 COMMON /zzzzzz/ x(1)
 
 IF (prt < 0)  GO TO 100
 trlr(1) = 101
 CALL rdtrl (trlr)
 IF (trlr(1) <= 0)  GO TO 100
 CALL fname (trlr,NAME)
 CALL intprt (x,opt,1,NAME)
 100 RETURN
END SUBROUTINE prtint
