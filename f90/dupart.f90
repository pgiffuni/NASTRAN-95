SUBROUTINE dupart
     
!     DRIVER FOR DMAP MODULE UPARTN
 
!     DMAP CALLING SEQUENCE IS
!     UPARTN    USET,KNN/KFF,KSF,KFS,KSS/C,N,N/C,N,F/C,N,S $
 
 INTEGER :: uset,scr1,sub0,sub1,iabit(16),NAME(2),ib(3)
 COMMON /BLANK / major(2),sub0(2),sub1(2)
 COMMON /bitpos/ ibit(32),iabit
 DATA    NAME  / 4HUPAR  ,4HTN   /
 DATA    uset  , knn, kff, ksf, kfs, kss, scr1 /  &
     101   , 102, 201, 202, 203, 204, 301  /
 
 
 nogo = 0
 
!     DECIDE IF CHARACTERS ARE LEGAL BIT NUMBERS
 
 ib(1) = major(1)
 ib(2) = sub0(1)
 ib(3) = sub1(1)
 
 loop30:  DO  j = 1,3
   DO  i = 1,32
     IF (ib(j) /= iabit(i)) CYCLE
     ib(j) = ibit(i)
     CYCLE loop30
   END DO
   
!     INVALID
   CALL mesage (59,ib(j),NAME)
   nogo = 1
 END DO loop30
 
 IF (nogo == 1) CALL mesage (-7,0,NAME)
 
 CALL upart (uset,scr1,ib(1),ib(2),ib(3))
 CALL mpart (knn,kff,ksf,kfs,kss)
 RETURN
END SUBROUTINE dupart
