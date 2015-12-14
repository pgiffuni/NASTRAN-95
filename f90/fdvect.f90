SUBROUTINE fdvect (delta,pk)
     
 
 REAL, INTENT(IN OUT)                     :: delta
 DOUBLE PRECISION, INTENT(IN)             :: pk
 INTEGER :: sysbuf,icore(1),mcb(7)
!    1,                NAME(2)
 DOUBLE PRECISION :: dcore(1),dmax
 COMMON /system/  ksystm(65)
 COMMON /zzzzzz/  core(1)
 COMMON /regean/  ia(14),ivect(7),ib(5),lc1,ib1(9),loads,lx,icount,  &
     lama,ibuck,nsym
 COMMON /packx /  it1p,it2p,iip,jjp,incrp
 EQUIVALENCE      (icore(1),core(1),dcore(1)),  &
     (ksystm(1),sysbuf),(ksystm(55),iprec)
!     DATA    NAME  /  4HFDVE,4HCT  /
 
 nprob = ia(3)
 kprec = ia(5)
 IF (kprec /= 1 .AND. kprec /= 2) kprec = iprec
 npro2 = nprob
 icnt  = ivect(2)
 im1   = 1
 lcore = (korsz(core)/2)*2 - lc1 - sysbuf
 x     = nprob
 y     = icount
 mcb(1)= ib(5)
 CALL rdtrl (mcb(1))
 
 CALL detfbs (npro2+1,icore(lcore+1),mcb,nprob,icount)
 
!     COPY FX ONTO IVECT + NORMALIZE
 
 ipm1 = ivect(1)
 IF (icnt == 0) GO TO 40
 CALL gopen (ivect(1),icore(lcore+1),0)
 CALL skprec (ivect(1),icnt)
 CALL CLOSE (ivect(1),2)
 im1 = 3
 40 CALL gopen (ivect,icore(lcore+1),im1)
 lcore = lcore - sysbuf
 IF (kprec == 2) GO TO 71
 xmax = 0.0
 DO   i = 1, nprob
   xmax = AMAX1(xmax,ABS(core(i)))
 END DO
 DO  i = 1,nprob
   core(i) = core(i)/xmax
 END DO
 GO TO 73
 71 dmax = 0.0D0
 DO  i = 1,nprob
   IF (DABS(dcore(i)) > dmax) dmax = DABS(dcore(i))
 END DO
 DO  i = 1,nprob
   dcore(i) = dcore(i)/dmax
 END DO
 73 CONTINUE
 it1p = kprec
 it2p = iprec
 iip  = 1
 jjp  = nprob
 incrp = 1
 CALL pack (core,ivect,ivect)
 CALL CLOSE (ivect(1),1)
 ipm1 = lama
 CALL gopen (lama,icore(lcore+1),3)
 dcore(1) = pk
 CALL WRITE (lama,core,iprec,1)
 CALL CLOSE (lama,2)
 RETURN
END SUBROUTINE fdvect
