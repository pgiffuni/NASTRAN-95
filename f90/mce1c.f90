SUBROUTINE mce1c
     
!     MCE1C PERFORMS A FORWARD-BACKWARD SUBSTITUTION WITH THE
!     TRIANGULAR FACTORS OF RM TO SOLVE FOR GM IN THE EQUATION
!     RM*GM = -RN.
 
 
 INTEGER :: uset  , rg    ,gm    ,scr1  ,scr2  ,scr3  ,rm    ,rn     ,  &
     u     , ux    ,rnx   ,gmx   ,prec  ,SIGN
 COMMON /BLANK / uset  ,rg    ,gm    ,scr1  ,scr2  ,scr3  ,rm    ,  &
     rn    ,l     ,u     ,mcb(7)
 COMMON /gfbsx / lx (7),ux (7), rnx(7),gmx(7),nz   ,prec  ,SIGN
 COMMON /zzzzzz/ z(1)
 COMMON /system/ ksystm(65)
 EQUIVALENCE     (ksystm(55),iprec)
 
!     INITIALIZE MATRIX CONTROL BLOCKS
 
 nz = korsz(z)
 lx(1)  = l
 CALL rdtrl (lx)
 ux(1)  = u
 CALL rdtrl (ux)
 rnx(1) = rn
 CALL rdtrl (rnx)
 gmx(1) = gm
 gmx(3) = rnx(3)
 gmx(4) = rnx(4)
 gmx(5) = iprec
 prec   = iprec
 SIGN   =-1
 
!     PERFORM SOLUTION
 
 CALL gfbs (z,z)
 
!     WRITE TRAILER
 
 CALL wrttrl (gmx)
 RETURN
END SUBROUTINE mce1c
