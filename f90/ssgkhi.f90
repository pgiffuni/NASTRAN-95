SUBROUTINE ssgkhi (treal,tint,fn)
     
!     THIS SUBROUTINE COMPUTES THE (5X1) KHI  VECTOR FOR USE BY TRBSC,
!                                           E
!     TRPLT, AND QDPLT.
 
!     WHEN PROCESSING THE TRPLT OR QDPLT THIS ROUTINE SHOULD BE CALLED
!     AFTER THE FIRST SUBTRIANGLE ONLY DUE TO THE D MATRIX ORIENTATION.
 
 
 REAL, INTENT(IN)                         :: treal(6)
 INTEGER, INTENT(IN OUT)                  :: tint(6)
 REAL, INTENT(IN)                         :: fn
 INTEGER :: INDEX(9)
 REAL :: khi
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /ssgtri/ d(9),khi(5),ks(30),p(6)
 COMMON /matout/ dum(7),alpha1,alpha2,alph12
 COMMON /trimex/ eid
 COMMON /system/ sysbuf,iout
 
!     DETERMINE TYPE OF TEMPERATURE DATA
 
 IF (tint(6) /= 1) GO TO 100
 
!     TEMPERATURE DATA IS TEMPP1 OR TEMPP3 TYPE.
 
 khi(1) = -alpha1*treal(2)*fn
 khi(2) = -alpha2*treal(2)*fn
 khi(3) = -alph12*treal(2)*fn
 GO TO 120
 
!     TEMPERATURE DATA IS TEMPP2 TYPE.
 
!     NO NEED TO COMPUTE DETERMINANT SINCE IT IS NOT USED SUBSEQUENTLY.
 
 100 ising = -1
 CALL invers (3,d(1),3,0,0,determ,ising,INDEX)
 IF (ising /= 2) GO TO 110
 WRITE  (iout,105) ufm,eid
 105 FORMAT (a23,' 4018, A SINGULAR MATERIAL MATRIX -D- FOR ELEMENT',  &
     i9,' HAS BEEN DETECTED BY ROUTINE SSGKHI', /26X,'WHILE ',  &
     'TRYING TO COMPUTE THERMAL LOADS WITH TEMPP2 CARD DATA.')
 CALL mesage (-61,0,0)
 110 CALL gmmats (d(1),3,3,0, treal(2),3,1,0, khi(1))
 khi(1) = khi(1)*fn
 khi(2) = khi(2)*fn
 khi(3) = khi(3)*fn
 120 khi(4) = 0.0
 khi(5) = 0.0
 RETURN
END SUBROUTINE ssgkhi
