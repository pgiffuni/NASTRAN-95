SUBROUTINE timtst
     
!     TIMETEST   /,/ C,N,N / C,N,M / C,N,T / C,N,O1 / C,N,O2 $
 
 INTEGER :: t,o1,o2
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm
 COMMON /BLANK / n,m,t,o1,o2
 COMMON /system/ isysbf, nout
 
 IF (o1 < 1 .OR. o1 > 2) GO TO 9901
 SELECT CASE ( o1 )
   CASE (    1)
     GO TO 100
   CASE (    2)
     GO TO  200
 END SELECT
 
 100 CONTINUE
 CALL timts1
 GO TO 900
 
 200 CONTINUE
 CALL timts2
 
 900 CONTINUE
 RETURN
 
!     ERROR MESSAGES
 
 9901 WRITE  (nout,9951) uwm
 9951 FORMAT (a25,' 2195, ILLEGAL VALUE FOR P4 =',i7)
 
 WRITE  (nout,9996)
 9996 FORMAT ('0*** MODULE TIMETEST TERMINAL ERROR.')
 
 RETURN
 
END SUBROUTINE timtst
