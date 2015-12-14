SUBROUTINE eadd (p,prec)
     
 
 DOUBLE PRECISION, INTENT(IN)             :: p(1)
 INTEGER, INTENT(IN)                      :: prec
 INTEGER :: ia(7)    ,ib(7)    ,ic(7)
 DOUBLE PRECISION :: beta(2)  , alpha(2)
 COMMON /regean/  im(7)    ,ik(7)    ,iev(7)    ,ka(5)    ,lc    ,  &
     nn(13)   ,ibuck
 COMMON /BLANK /  xx
 COMMON /saddx /  nomat    ,nz       ,mcbs(67)
 COMMON /zzzzzz/  core(1)
 EQUIVALENCE      (mcbs( 1),ia(1))   ,(mcbs( 8),ialp)     ,  &
     (mcbs( 9),alpha(1)),(mcbs(13),ib(1))    ,  &
     (mcbs(20),ibeta)   ,(mcbs(21),beta(1))  , (mcbs(61),ic(1))
 
 nz = (korsz(core)/2)*2 - lc
 DO  i = 1,7
   ia(i) = im(i)
   ib(i) = ik(i)
   ic(i) = ik(i)
 END DO
 ic(1) = ka(1)
 kprec = ik(5)
 IF (prec >= 1 .AND. prec <= 4) kprec = prec
 ialp = kprec
 alpha(1) = p(1)
 ibeta  = kprec
 beta(1)= 1.0D0
 nomat  = 2
 CALL sadd (core,core)
 CALL wrttrl (ic)
 RETURN
END SUBROUTINE eadd
