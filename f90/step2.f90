SUBROUTINE step2 (u2,u1,u0,p,ibuf)
     
!     STEP2 WILL INTEGRATE FORWARD ONE TIME STEP
 
!     THIS ROUTINE IS SUITABLE FOR DOUBLE PRECISION OPERATION
 
 
 DOUBLE PRECISION, INTENT(IN OUT)         :: u2(1)
 DOUBLE PRECISION, INTENT(IN OUT)         :: u1(1)
 DOUBLE PRECISION, INTENT(IN OUT)         :: u0(1)
 DOUBLE PRECISION, INTENT(IN OUT)         :: p(1)
 INTEGER, INTENT(IN OUT)                  :: ibuf(1)
 INTEGER :: FILE(7) ,sqr     ,rdp     ,dum
 
 COMMON /trdxx /  dum(21) ,iscr1   ,iscr2   ,iscr3   ,iscr4   ,  &
     iscr5   ,iscr6   ,iopen   ,isym
 COMMON /names /  dumm(8) ,rdp     ,dumn(2) ,sqr
 COMMON /infbsx/  ifil(7) ,ifilu(7)
 
 FILE(1) = iscr1
 FILE(2) = dum(2)
 FILE(4) = sqr
 
!     TELL MATVC2/INVFBS FILES ARE OPEN
 
 iopen = 1
 
!     FORM R.H.S. OF THE INTEGRATION EQUATION
 
 CALL matvc2 (u1(1),p(1),FILE,ibuf)
 FILE(1) = iscr4
 CALL matvc2 (u0(1),p(1),FILE,ibuf)
 
!     CALL INVFBS/FBSINT TO DO THE FORWARD/BACKWARD PASS
 
 ifil(1)  = iscr2
 ifilu(1) = iscr3
 CALL rdtrl (ifil)
 CALL rdtrl (ifilu)
 ifil(5) = rdp
 ifil(3) = dum(3)
 IF (isym == 1) iopen = - 20
 IF (isym == 1) CALL invfbs (p(1),u2(1),ibuf)
 IF (isym == 0) CALL fbsint (p(1),u2(1))
 iopen = 0
 RETURN
END SUBROUTINE step2
