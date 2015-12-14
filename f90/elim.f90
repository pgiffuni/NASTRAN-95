SUBROUTINE elim (in1,in2,in3,in4,out1,scr1,scr2,scr3)
     
!     ELIM EVALUATES THE MATRIX EQUATION -
 
!     OUT1 = IN1 + IN4(T)*IN2 + IN2(T)*IN4 + IN4(T)*IN3*IN4
 
 
 INTEGER, INTENT(IN)                      :: in1
 INTEGER, INTENT(IN)                      :: in2
 INTEGER, INTENT(IN)                      :: in3
 INTEGER, INTENT(IN)                      :: in4
 INTEGER, INTENT(IN)                      :: out1
 INTEGER, INTENT(IN)                      :: scr1
 INTEGER, INTENT(IN)                      :: scr2
 INTEGER, INTENT(IN)                      :: scr3
 INTEGER :: filea ,fileb ,filec  ,  &
     filed ,t     ,signab,signc ,prec  ,rdp   ,plus   , scrtch
 DIMENSION       filea(7)     ,fileb(7)     ,filec(7)     ,filed(7)
!    1,               MCB(7)
 COMMON /mpyadx/ filea ,fileb ,filec ,filed ,nz    ,t     ,signab ,  &
     signc ,prec  ,scrtch
 COMMON /system/ idum(54)     ,iprec
 COMMON /zzzzzz/ z(1)
 DATA    plus  / +1 /
 
 rdp = iprec
 
!     PERFORM GENERAL INITIALIZATION
 
 nz     = korsz(z)
 signab = plus
 signc  = plus
 prec   = rdp
 scrtch = scr3
 
!     INITIALIZE MATRIX CONTROL BLOCKS FOR IN3,IN4,IN2 AND SCR1
 
 filea(1) = in3
 CALL rdtrl (filea)
 fileb(1) = in4
 CALL rdtrl (fileb)
 filec(1) = in2
 CALL rdtrl (filec)
 filed(1) = scr1
 filed(3) = filec(3)
 filed(4) = filec(4)
 filed(5) = rdp
 
!     COMPUTE SCR1 = IN3*IN4 + IN2
 
 t = 0
 CALL mpyad (z,z,z)
 
!     SAVE MATRIX CONTROL BLOCK FOR SCR1
 
!     DO 41 I = 1,7
 CALL wrttrl (filed)
 
!     INITIALIZE MATRIX CONTROL BLOCKS FOR IN2, IN4, IN1 AND SCR2
 
 DO  i = 1,7
   filea(i) = filec(i)
 END DO
 filec(1) = in1
 CALL rdtrl (filec)
 filed(1) = scr2
 filed(3) = filec(3)
 filed(4) = filec(4)
 
!     COMPUTE SCR2 = IN2(T)*IN4 + IN1
 
 t = 1
 CALL mpyad  (z,z,z)
 CALL wrttrl (filed)
 
!     INITIALIZE MATRIX CONTROL BLOCKS FOR IN4,SCR1,SCR2 AND OUT1
 
 filea(1) = fileb(1)
 fileb(1) = scr1
 filec(1) = filed(1)
 CALL rdtrl (filea)
 CALL rdtrl (fileb)
 CALL rdtrl (filec)
 filed(1) = out1
 filed(3) = filec(3)
 filed(4) = filec(4)
 
!     COMPUTE  OUT1= IN4(T)*SCR1 + SCR2
 
 t = 1
 CALL mpyad (z,z,z)
 
!     WRITE TRAILER FOR OUT1 AND RETURN
 
 CALL wrttrl (filed)
 RETURN
END SUBROUTINE elim
