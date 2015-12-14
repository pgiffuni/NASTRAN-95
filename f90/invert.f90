SUBROUTINE invert( ia,ib,scr1)
     
!     DRIVER  FOR  INVTR
 
!     INVERTS  LOWER OR UPPER TRIANGLE  IA  ONTO  IB
 
!     SCR1 WILL BE USED  ONLY IF  IA  IS AN  UPPER  TRIANGLE
 
 
 INTEGER, INTENT(IN)                      :: ia
 INTEGER, INTENT(IN)                      :: ib
 INTEGER, INTENT(IN)                      :: scr1
 INTEGER :: fa,fb,scrfil,prec, NAME(2)
 
 COMMON /invtrx/ fa(7),fb(7),scrfil,nx,prec
 COMMON / zzzzzz/ z(1)
 
 DATA NAME /4HINVE,4HRT  /
 
!     FILL  MATRIX  CONTROL  BLOCKS  FOR  A  AND  B
 
 fa(1) = ia
 CALL rdtrl(fa)
 fb(1) = ia
 CALL rdtrl(fb)
 fb(1) = ib
 scrfil = scr1
 prec = fa(5)
 nx  =  korsz(z)
 CALL invtr(*50,z,z)
 CALL wrttrl(fb)
 RETURN
 
!     SINGULAR  MATRIX
 
 50 CALL mesage(-5,fa,NAME)
 GO TO 50
END SUBROUTINE invert
