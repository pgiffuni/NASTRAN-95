SUBROUTINE frd2b (a,alp,b,bet,c,gam,d,del,e,eps,out)
     
!     ADD UP MATRICIES
 
 
 INTEGER, INTENT(IN)                      :: a
 REAL, INTENT(IN)                         :: alp(2)
 INTEGER, INTENT(IN)                      :: b
 REAL, INTENT(IN)                         :: bet(2)
 INTEGER, INTENT(IN)                      :: c
 REAL, INTENT(IN)                         :: gam(2)
 INTEGER, INTENT(IN)                      :: d
 REAL, INTENT(IN)                         :: del(2)
 INTEGER, INTENT(IN)                      :: e
 REAL, INTENT(IN)                         :: eps(2)
 INTEGER, INTENT(IN OUT)                  :: out
 INTEGER :: typa,typb,typc,typd,TYPE
 REAL :: z(1)
 COMMON /system/ ksystm(54), iprec
 COMMON /zzzzzz/ z
 COMMON /saddx / nomat,lcore,mcba(7),typa,alpha(4),mcbb(7),typb,  &
     beta(4),mcbc(7),typc,gama(4),mcbd(7),typd,  &
     delta(4),mcbe(7),TYPE,epsln(4),mc(7)
 COMMON /frd2bc/ ih
 
 nc    = korsz(z)
 nomat = 5
 lcore = nc
 typa  = 3
 typb  = 3
 typc  = 3
 typd  = 3
 TYPE  = 3
 alpha(1) = alp(1)
 alpha(2) = alp(2)
 beta(1)  = bet(1)
 beta(2)  = bet(2)
 gama(1)  = gam(1)
 gama(2)  = gam(2)
 delta(1) = del(1)
 delta(2) = del(2)
 epsln(1) = eps(1)
 epsln(2) = eps(2)
 mcba(1)  = a
 mcbb(1)  = b
 mcbc(1)  = c
 mcbd(1)  = d
 mcbe(1)  = e
 CALL rdtrl (mcba)
 CALL rdtrl (mcbb)
 CALL rdtrl (mcbc)
 CALL rdtrl (mcbd)
 CALL rdtrl (mcbe)
 ifo = 6
 ity = 3
 IF (ih == 0 .AND. iprec == 2) ity = 4
 
!     IH IN /FRD2BC/ IS INITIALIZED BY ROUTINE FRRD2.
!     (COMPLEX D.P. ARITHMETIC IS USED IF IH = 0)
 
 n = 0
 DO  i = 1,49,12
   IF (mcba(i  ) < 0) mcba(i) = 0
   IF (mcba(i+1) == 0) mcba(i) = 0
   IF (mcba(i  ) == 0) CYCLE
   IF (n == 0) n = mcba(i+1)
   irow = mcba(i+2)
   IF (mcba(i+3) /= 6) ifo  = 1
 END DO
 CALL makmcb (mc,out,irow,ifo,ity)
 mc(2) = n
 CALL sadd (z,z)
 CALL wrttrl (mc)
 CALL dmpfil (-out,z,nc)
 RETURN
END SUBROUTINE frd2b
