SUBROUTINE step (u2,u1,u0,p,ibuf)
     
!     STEP WILL INTEGRATE FORWARD ONE TIME STEP
 
!     THIS ROUTINE IS SUITABLE FOR SINGLE PRECISION OPERATION
 
 
 REAL, INTENT(IN OUT)                     :: u2(1)
 REAL, INTENT(IN OUT)                     :: u1(1)
 REAL, INTENT(IN OUT)                     :: u0(1)
 REAL, INTENT(IN OUT)                     :: p(1)
 INTEGER, INTENT(IN OUT)                  :: ibuf(1)
 INTEGER :: FILE(7)  ,sqr     ,rsp     ,dum
 
 COMMON /trdxx / dum(21)  ,iscr1   ,iscr2   ,iscr3   ,iscr4   ,  &
     iscr5    ,iscr6   ,iopen   ,isym
 COMMON /names / dumm(7)  ,rsp     ,dumn(3) ,sqr
 COMMON /infbsx/ ifil(7)  ,ifilu(7)
 
 FILE(1) = iscr1
 FILE(2) = dum(2)
 FILE(4) = sqr
 
!     TELL MATVEC/INTFBS FILES ARE OPEN
 
 iopen = 1
 
!     FORM R.H.S. OF THE INTEGRATION EQUATION
 
 CALL matvec (u1(1),p(1),FILE,ibuf)
 FILE(1) = iscr4
 CALL matvec (u0(1),p(1),FILE,ibuf)
 
!     CALL INTFBS/FBSINT TO DO THE FORWARD/BACKWARD PASS
 
 ifil(1)  = iscr2
 ifilu(1) = iscr3
 CALL rdtrl (ifil)
 CALL rdtrl (ifilu)
 ifil(5)  = rsp
 ifil(3)  = dum(2)
 IF (isym == 1) CALL intfbs (p(1),u2(1),ibuf)
 IF (isym == 0) CALL fbsint (p(1),u2(1))
 iopen = 0
 RETURN
END SUBROUTINE step
