SUBROUTINE dtranp
     
!     DRIVER OF MATRIX TRANSPOSE MODULE
 
!     TRNSP    IA/IAT/C,N,IXX  $
 
!     THE DIAGONALS OF THE LOWER OR UPPER TRIANGULAR MATRICES ARE
!     REPLACED BY UNITY (1.0) IF IXX IS ONE. (DEFAULT IS ZERO)
 
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm
 COMMON /BLANK / ixx
 COMMON /system/ ibuf,nout
 COMMON /zzzzzz/ core(1)
 COMMON /trnspx/ ia(7),iat(7),lcore,nscr,iscr(8)
 DATA    in1   , in2  /101, 201 /
 
 ia(1) = in1
 CALL rdtrl (ia(1))
 IF (ia(1) > 0) GO TO 20
 WRITE  (nout,10) uwm
 10 FORMAT (a25,' FROM TRNSP, MISSING INPUT DATA BLOCK FOR MATRIX ',  &
     'TRANSPOSE')
 GO TO 60
 20 iat(1) = in2
 iat(2) = ia(3)
 iat(3) = ia(2)
 iat(4) = ia(4)
 iat(5) = ia(5)
 iat(6) = 0
 iat(7) = 0
 lcore  = korsz(core)
 nscr = 8
 DO  i = 1,nscr
   iscr(i) = 300 + i
 END DO
 IF (ixx == 1) ixx = -123457890
 CALL trnsp (core(1))
 CALL wrttrl (iat(1))
 
 60 RETURN
END SUBROUTINE dtranp
