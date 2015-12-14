SUBROUTINE ssg2a (pg,pnbar,pm,pvact)
     
 
 INTEGER, INTENT(IN)                      :: pg
 INTEGER, INTENT(IN)                      :: pnbar
 INTEGER, INTENT(IN)                      :: pm
 INTEGER, INTENT(IN)                      :: pvact
 INTEGER :: rule,pvect(7),core(6)
 COMMON /parmeg/ ia1(7),ia11(7),ia12(7),ia21(7),ia22(7),lcr,rule
 COMMON /patx  / lcore,n,no(4)
 COMMON /zzzzzz/ icore(1)
 EQUIVALENCE     (icore(1),core(1))
 
 
 pvect(1)= pvact
 CALL rdtrl (pvect)
 ia1(1)  = pg
 CALL rdtrl (ia1)
 ia11(1) = pnbar
 ia12(1) = pm
 DO  i = 2,5
   ia11(i) = ia1(i)
   ia12(i) = ia1(i)
 END DO
 ia11(3) = n
 ia12(3) = no(1)
 ia21(1) = 0
 ia22(1) = 0
 rule    = 0
 lcr     = korsz(core)
 core(1) = 0
 core(2) = 1
 core(3) = ia1(2)
 core(4) = 2
 core(5) = 1
 core(6) = 0
 CALL partn (core,pvect,core)
 IF (ia11(1) /= 0) CALL wrttrl (ia11)
 IF (ia12(1) /= 0) CALL wrttrl (ia12)
 RETURN
END SUBROUTINE ssg2a
