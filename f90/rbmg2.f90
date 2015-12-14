SUBROUTINE rbmg2
     
 INTEGER :: scr1,    scr2,    scr3,    scr4
 DOUBLE PRECISION :: det(2)
 COMMON /BLANK /  jpowr,   detrm
 COMMON /sfact /  qq(29)
 EQUIVALENCE      (qq(25), det(1)), (qq(29), ipwr)
 DATA    kll   ,  lll,     scr1,    scr2,    scr3,    scr4  /  &
     101   ,  201,     301,     302,     303,     304   /
 
!     DECOMPOSE KLL INTO LLL
 
 CALL factor (kll,lll,scr1,scr2,scr3,scr4)
 jpowr = ipwr
 detrm = det(1)
 RETURN
END SUBROUTINE rbmg2
