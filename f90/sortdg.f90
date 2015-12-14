SUBROUTINE sortdg (stk1,stk2,x1,x2,ndeg)
     
!     SORTDG SORTS STK2 BY DEGREE OF THE NODE AND ADDS IT TO THE END
!     OF STK1 IN ORDER OF LOWEST TO HIGHEST DEGREE.  X1 AND X2 ARE THE
!     NUMBER OF NODES IN STK1 AND STK2 RESPECTIVELY.
 
 
 INTEGER, INTENT(OUT)                     :: stk1(1)
 INTEGER, INTENT(IN OUT)                  :: stk2(1)
 INTEGER, INTENT(OUT)                     :: x1
 INTEGER, INTENT(IN)                      :: x2
 INTEGER, INTENT(IN)                      :: ndeg(1)
 INTEGER :: temp
 
 COMMON /bandg /  n,        idpth
 
 ind=x2
 10 itest=0
 ind=ind-1
 IF (ind < 1) GO TO 40
 DO  i=1,ind
   j=i+1
   istk2=stk2(i)
   jstk2=stk2(j)
   IF (ndeg(istk2) <= ndeg(jstk2)) CYCLE
   itest=1
   temp=stk2(i)
   stk2(i)=stk2(j)
   stk2(j)=temp
 END DO
 IF (itest == 1) GO TO 10
 40 DO  i=1,x2
   x1=x1+1
   stk1(x1)=stk2(i)
 END DO
 RETURN
END SUBROUTINE sortdg
