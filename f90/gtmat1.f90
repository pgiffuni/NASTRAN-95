SUBROUTINE gtmat1 (sym,tt)
     
!     THIS SUBROUTINE PROCESSES TRANSFORMATION MATRICES
!     IT IS CALLED ONLY BY CMSFIL
 
 
 INTEGER, INTENT(IN)                      :: sym
 REAL, INTENT(OUT)                        :: tt(3,3)
 EXTERNAL        rshift   ,andf     ,orf
 INTEGER :: trn      , orf     ,tran    ,ecpt1    ,  &
     andf     ,chk1     ,chk2    ,NAME(2) ,rshift
 DIMENSION       ecpt(4)  ,tid(3,3) , list(32),symm(6,6),  &
     smat(6,3),prod(6)  ,tc(3,3) ,tg6(6,6),tg(3,3)  , t(6,6)
 DIMENSION       acpt(1)
 COMMON /gtmatx/ loc1     ,len1     ,trn     ,tt6(6,6),tc6(6,6)
 COMMON /zzzzzz/ z(1)
 EQUIVALENCE     (ecpt1,ecpt(1))
 EQUIVALENCE     (iflag,rflag)
 DATA    tid   / 1., 0., 0., 0., 1., 0., 0., 0., 1. /
 DATA    smat  /-1., 1., 1., 1.,-1.,-1., 1.,-1., 1.,-1., 1.,-1.,  &
     1., 1.,-1.,-1.,-1., 1.  /
 DATA    NAME  / 4HGTMT, 4H1Z            /
 
 ikind = 0
 DO  i = 1,6
   DO  j = 1,6
     tt6(i,j) = 0.0
   END DO
 END DO
 IF (trn == 0 .AND. sym == 0) GO TO 170
 IF (loc1 == 0 .OR. trn == 0) GO TO 30
 CALL pretrs (z(loc1),len1)
 ikind = orf(ikind,1)
 DO  i = 2,4
   ecpt(i) = 0.0
 END DO
 ecpt1 = trn
 CALL transs (ecpt,tt)
 GO TO 50
 30 DO  i = 1,3
   DO  j = 1,3
     tt(i,j) = tid(i,j)
   END DO
 END DO
 50 DO  i = 1,3
   DO  j = 1,3
     tt6(i  ,j  ) = tt(i,j)
     tt6(i+3,j+3) = tt(i,j)
   END DO
 END DO
 DO  i = 1,6
   DO  j = 1,6
     symm(i,j) = 0.0
   END DO
 END DO
 IF (sym == 0) GO TO 120
 ikind = orf(ikind,1)
 CALL decode (sym,list,ndir)
 DO  i = 1,6
   prod(i) = 1.0
 END DO
 DO  i = 1,ndir
   idir = list(i) + 1
   idir = 4 - idir
   DO  j = 1,6
     prod(j) = prod(j)*smat(j,idir)
   END DO
 END DO
 DO  i = 1,6
   symm(i,i) = prod(i)
 END DO
 GO TO 140
 120 DO  i = 1,6
   symm(i,i) = 1.0
 END DO
 140 CALL gmmats (tt6,6,6,0, symm,6,6,0, t)
 DO  i = 1,6
   DO  j = 1,6
     tt6(i,j) = t(i,j)
   END DO
 END DO
 DO  i = 1,3
   DO  j = 1,3
     tt(i,j) = tt6(i,j)
   END DO
 END DO
 isav = ikind
 RETURN
 
 170 DO  i = 1,6
   tt6(i,i) = 1.0
 END DO
 DO  i = 1,3
   DO  j = 1,3
     tt(i,j) = tid(i,j)
   END DO
 END DO
 isav = ikind
 chk1 = 13579
 RETURN
 
 
 ENTRY gtmat2 (loc2,len2,acpt,tc)
!     ================================
 
 ikind = isav
 DO  i = 1,6
   DO  j = 1,6
     tc6(i,j) = 0.0
   END DO
 END DO
 rflag = acpt(1)
 IF (loc2 == 0 .OR. iflag == 0) GO TO 210
 CALL pretrs (z(loc2),len2)
 CALL transs (acpt,tc)
 ikind = orf(ikind,2)
 GO TO 230
 210 DO  i = 1,3
   DO  j = 1,3
     tc(i,j) = tid(i,j)
   END DO
 END DO
 230 DO  i = 1,3
   DO  j = 1,3
     tc6(i  ,j  ) = tc(i,j)
     tc6(i+3,j+3) = tc(i,j)
   END DO
 END DO
 chk2 = 24680
 RETURN
 
 
 ENTRY gtmat3 (tran,tg,tg6,ihelp)
!     ================================
 
 IF (chk1 /= 13579 .AND. chk2 /= 24680) CALL mesage (-37,0,NAME)
 DO  i = 1,6
   DO  j = 1,6
     tg6(i,j) = 0.0
   END DO
 END DO
 IF (tran < 0.0) THEN
   GO TO   340
 ELSE IF (tran == 0.0) THEN
   GO TO   330
 END IF
 310 CALL pretrs (z(loc1),len1)
 DO  i = 2,4
   ecpt(i) = 0.0
 END DO
 ecpt1 = tran
 ikind = orf(ikind,8)
 IF (tran /= trn) ikind = orf(ikind,16)
 CALL transs (ecpt,tg)
 ikind = orf(ikind,4)
 GO TO 370
 330 ikind = orf(ikind,4)
 340 DO  i = 1,3
   DO  j = 1,3
     tg(i,j) = tid(i,j)
   END DO
 END DO
 IF (andf(rshift(ikind,1),1) /= 1 .OR. tran /= -1) GO TO 360
 CALL gmmats (tt6,6,6,0, tc6,6,6,0, tg6)
 ihelp = ikind
 RETURN
 
 360 CONTINUE
 370 DO  i = 1,3
   DO  j = 1,3
     tg6(i  ,j  ) = tg(i,j)
     tg6(i+3,j+3) = tg(i,j)
   END DO
 END DO
 ihelp = ikind
 RETURN
END SUBROUTINE gtmat1
