SUBROUTINE form22 (uddip1, udiprm, uiprm, piprm, ibuf)
!*******
!     FORM22 GENERATES THE VECTORS NECESSARY TO CHANGE THE TIME STEP
 
!     THIS ROUTINE IS SUITABLE FOR DOUBLE PRECISION OPERATION
!*******
 
 DOUBLE PRECISION, INTENT(IN OUT)         :: uddip1(1)
 DOUBLE PRECISION, INTENT(IN OUT)         :: udiprm(1)
 DOUBLE PRECISION, INTENT(IN OUT)         :: uiprm(1)
 DOUBLE PRECISION, INTENT(IN OUT)         :: piprm(1)
 INTEGER, INTENT(IN OUT)                  :: ibuf(1)
 
 
 
 
 COMMON   /trdxx /  ifilk(7)  ,ifilm(7) ,ifilb(7)
!*******
!     FORM UDOT(I+1), UDDOT(I+1), UDOT-(I), AND U-(I)
!*******
 CALL matvc2(uddip1(1),piprm(1),ifilm(1),ibuf)
 CALL matvc2(udiprm(1),piprm(1),ifilb(1),ibuf)
 CALL matvc2(uiprm(1),piprm(1),ifilk(1),ibuf)
 RETURN
END SUBROUTINE form22
