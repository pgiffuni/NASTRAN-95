SUBROUTINE mce1b
     
!     MCE1B DECOMPOSES RM INTO LOWER AND UPPER TRIANGULAR FACTORS
 
 INTEGER :: uset  , rg    ,gm    ,scr1  ,scr2  ,scr3  ,rm    ,rn   ,  &
     a     , ux    ,scrx1 ,scrx2 ,scrx3 ,nam(2),u
 DOUBLE PRECISION :: det  ,mindia
 COMMON /BLANK / uset  ,rg    ,gm    ,scr1  ,scr2  ,scr3  ,rm    ,  &
     rn    ,l     ,u     ,mcb(7)
 COMMON /dcompx/ a(7)  ,lx(7) ,ux(7) ,scrx1 ,scrx2 ,scrx3 ,det   ,  &
     power ,nz    ,mindia
 COMMON /zzzzzz/ z(1)
 DATA    nam   / 4HMCE1,4HB   /
 
!     INITIALIZE MATRIX CONTROL BLOCKS
 
 nz    = korsz(z)
 a(1)  = rm
 CALL rdtrl (a)
 lx(1) = l
 lx(3) = a(3)
 lx(4) = 4
 lx(5) = a(5)
 ux(1) = u
 ux(3) = a(3)
 ux(4) = 5
 ux(5) = a(5)
 scrx1 = scr1
 scrx2 = scr2
 scrx3 = scr3
 
!     PERFORM DECOMPOSITION
 
 CALL decomp (*40,z,z,z)
 
!     WRITE TRAILERS
 
 CALL wrttrl (lx)
 CALL wrttrl (ux)
 RETURN
 
!     FATAL ERROR MESSAGE FOR SINGULAR MATRIX
 
 40 CALL mesage (-5,rm,nam)
 RETURN
END SUBROUTINE mce1b
