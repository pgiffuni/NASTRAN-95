SUBROUTINE factor (INPUT,lower,scr1,scr2,scr3,scr4)
     
 
 INTEGER, INTENT(IN)                      :: INPUT
 INTEGER, INTENT(IN OUT)                  :: lower
 INTEGER, INTENT(IN OUT)                  :: scr1
 INTEGER, INTENT(IN)                      :: scr2
 INTEGER, INTENT(IN)                      :: scr3
 INTEGER, INTENT(IN)                      :: scr4
 INTEGER :: bcd(2)
 DOUBLE PRECISION :: det
 COMMON /system/  ksystm(65)
 COMMON /sfact /  filea(7),filel(7),fileu(7),scr1fl,scr2fl,nz    ,  &
     det(2)  ,p       ,scr3fl  ,xx3   ,xx4   ,chl
 COMMON /zzzzzz/  z(1)
 DATA    lowtri/  4 /
 DATA    bcd   /  4HFACT,4HOR   /
 
!     INITIALIZE MATRIX CONTROL BLOCKS AND SFACT COMMON
 
 nz = korsz(z)
 filea(1) = INPUT
 CALL rdtrl (filea)
 CALL makmcb (filel,lower,filea(3),lowtri,filea(5))
 fileu(1) = IABS(scr1)
 scr1fl = scr2
 scr2fl = scr3
 scr3fl = scr4
 chl = 0
 IF (scr1 < 0) chl = 1
 
!     DECOMPOSE INPUT MATRIX INTO LOWER TRIANGULAR FACTOR.
 
 CALL sdcomp (*40,z,z,z)
 
!     WRITE TRAILER FOR LOWER TRIANGULAR FACTOR.
 
 CALL wrttrl (filel)
 RETURN
 
!     FATAL ERROR MESSAGE FOR SINGULAR INPUT MATRIX
 
 40 CALL mesage (-5,INPUT,bcd)
 RETURN
END SUBROUTINE factor
