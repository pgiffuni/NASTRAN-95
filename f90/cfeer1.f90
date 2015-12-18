SUBROUTINE cfeer1
     
!     CFEER1 INITIALIZES AND CALLS SUBROUTINE SADD FOR CFCNTL
 
 LOGICAL :: no_b     ,qpr
 INTEGER :: scr1     ,scr2     ,scr11    ,sqr       ,  &
     typout   ,ifila(7) ,ifilb(7) ,ifilc(7)
 DOUBLE PRECISION :: alpha(2) ,beta(2)  ,lambda   ,dz(1)
 DIMENSION         salpha(4),sbeta(4)
 COMMON  /feeraa/  ik(7)    ,im(7)    ,ib(7)    ,dum(15)   ,  &
     scr1     ,scr2     ,scr(8)   ,scr11     , dumaa(91),mcblmb(7)
 COMMON  /feerxc/  lambda(2),dumxc(12),no_b     ,dumxc2(4) , qpr
 COMMON  /saddx /  nomat    ,nz       ,mcbs(67)
 COMMON  /names /  dumm(10) ,cdp      ,sqr
 COMMON  /zzzzzz/  z(1)
 COMMON  /unpakx/  typout   ,irow     ,nlast    ,incr
 COMMON  /system/  ksystm(65)
 
 EQUIVALENCE      (mcbs( 1), ifila(1)), (mcbs( 8), itypal  ),  &
     (mcbs(61), ifilc(1)), (mcbs(13), ifilb(1)),  &
     (mcbs(20), itypbt  ), (mcbs(21), beta(1) ),  &
     (mcbs( 9), alpha(1)), (iprec,  ksystm(55)),  &
     (alpha(1),salpha(1)), (beta(1),  sbeta(1)),  &
     (z(1)    , dz(1)   ), (nout,   ksystm(2) )
 
!     FORM   -(B + LAMBDA*M)  ON SCR2
 
 itype    = iprec + 2
 nomat    = 2
 DO  i  = 1,7
   ifila(i) = im(i)
   ifilb(i) = ib(i)
 END DO
 IF (iprec == 2) GO TO 15
 salpha(1)=-SNGL(lambda(1))
 salpha(2)=-SNGL(lambda(2))
 salpha(3)= 0.
 salpha(4)= 0.
 sbeta(1) =-1.
 sbeta(2) = 0.
 sbeta(3) = 0.
 sbeta(4) = 0.
 GO TO 16
 15 alpha(1) =-lambda(1)
 alpha(2) =-lambda(2)
 beta(1)  =-1.d0
 beta(2)  =  0.d0
 16 itypal   = itype
 itypbt   = itype
 nz       = korsz(z)
 ifilc(1) = scr2
 ifilc(2) = ik(2)
 ifilc(3) = ik(3)
 ifilc(4) = 1
 ifilc(5) = itype
 IF (no_b) GO TO 100
 CALL sadd (z,z)
 
!---------- SPECIAL PRINT ------------------------------
 
 IF (.NOT.qpr) GO TO 25
 WRITE  (nout,2)
 2 FORMAT (1H0,//7H cfeer1,//)
 typout= itype
 irow  = 1
 nlast = ik(2)
 limit = 2*nlast
 incr  = 1
 ibuf  = nz - ksystm(1) - 2
 CALL gopen (ifilc(1),z(ibuf),0)
 DO  i = 1,nlast
   WRITE  (nout,1) i
   1 FORMAT (7H column,i4)
   CALL unpack (*20,ifilc(1),z)
   IF (iprec == 2) WRITE (nout,3) (dz(j),j=1,limit)
   IF (iprec /= 2) WRITE (nout,5) ( z(j),j=1,limit)
   20 CONTINUE
 END DO
 CALL CLOSE (ifilc(1),1)
 3 FORMAT (1H ,13(10H----------)/(1H ,4D25.16))
 5 FORMAT (1H ,13(10H----------)/(1H ,4E25.16))
 25 CONTINUE
 
 
!     FORM  (LAMBDA**2*M + LAMBDA*B + K)  ON SCR1
 
 DO  i  = 1,7
   ifila(i) = ik(i)
 END DO
 ifilb(1) = ifilc(1)
 ifilb(2) = ik(2)
 ifilb(3) = ik(3)
 ifilb(4) = sqr
 ifilb(5) = itype
 IF (iprec == 2) GO TO 35
 salpha(1) = 1.
 salpha(2) = 0.
 salpha(3) = 0.
 salpha(4) = 0.
 sbeta(1)  =-SNGL(lambda(1))
 sbeta(2)  =-SNGL(lambda(2))
 sbeta(3)  = 0.
 sbeta(4)  = 0.
 GO TO 50
 35 alpha(1) = 1.d0
 alpha(2) = 0.d0
 beta(1)  =-lambda(1)
 beta(2)  =-lambda(2)
 50 ifilc(1) = scr1
 CALL sadd (z,z)
 
!---------- SPECIAL PRINT ------------------------------
 
 IF (.NOT.qpr) GO TO 75
 WRITE  (nout,4)
 4 FORMAT (1H ,13(10H----------),//,19H the dynamic matrix,//)
 CALL gopen (ifilc(1),z(ibuf),0)
 DO  i = 1,nlast
   WRITE (nout,1) i
   CALL unpack (*70,ifilc(1),z)
   IF (iprec == 2) WRITE (nout,3) (dz(j),j=1,limit)
   IF (iprec /= 2) WRITE (nout,5) ( z(j),j=1,limit)
   70 CONTINUE
 END DO
 CALL CLOSE (ifilc(1),1)
 75 CONTINUE
 
!-------------------------------------------------------
!     MCBLMB NOT USED WHEN DAMPING MATRIX ABSENT
 
 DO  i = 1,7
   mcblmb(i) = ifilb(i)
 END DO
 GO TO 200
 
!     DAMPING MATRIX ABSENT
 
 100 DO  i = 1,7
   ifilb(i) = ik(i)
 END DO
 IF (iprec == 2) GO TO 120
 salpha(1) = SNGL(lambda(1)**2 - lambda(2)**2)
 salpha(2) = 2.*SNGL(lambda(1)*lambda(2))
 sbeta(1)  = 1.
 GO TO 130
 120 alpha(1) = lambda(1)**2 - lambda(2)**2
 alpha(2) = 2.d0*lambda(1)*lambda(2)
 beta(1)  = 1.d0
 
!----------- LOGIC FOR SPECIAL PRINT -------------------------
 
 130 IF (.NOT.qpr) GO TO 50
 typout= itype
 irow  = 1
 nlast = ik(2)
 limit = 2*nlast
 incr  = 1
 ibuf  = nz - ksystm(1) - 2
!-------------------------------------------------------------
 
 GO TO 50
 
 200 RETURN
END SUBROUTINE cfeer1
