SUBROUTINE initl2 (offset,deltt)
     
!     INITL2 WILL COMPUTE THE STARTING VALUES FOR THE INTEGRATION
!     ROUTINE
 
!     THIS ROUTINE IS SUITABLE FOR DOUBLE PRECISION OPERATION
 
 
 INTEGER, INTENT(IN)                      :: offset
 REAL, INTENT(IN)                         :: deltt
 INTEGER :: rsp      ,filem    ,fileb    ,  &
     filek    ,sqr      ,FILE     ,ifila(7) , ifilb(7) ,ifilc(7) ,NAME(2)  ,rdp
 DOUBLE PRECISION :: det      ,mindia   ,alpha(2) ,beta(2)
 COMMON /saddx /  nomat    ,nz       ,mcbs(67)
 COMMON /names /  rd       ,rdrew    ,wrt      ,wrtrew   ,  &
     rew      ,norew    ,eofnrw   ,rsp      , rdp      ,csp      ,cdp      ,sqr
 COMMON /sfact /  ifa(7)   ,ifl(7)   ,ifu(7)   ,isc1     ,  &
     isc2     ,nxx      ,id(5)    ,isc3     , id1(2)   ,ichl
 COMMON /dcompx/  ia(7)    ,il(7)    ,iu(7)    ,iscr10   ,  &
     iscr20   ,iscr30   ,det      ,power    , nx       ,mindia
 COMMON /trdxx /  filek(7) ,filem(7) ,fileb(7) ,  &
     iscr1    ,iscr2    ,iscr3    ,iscr4    , iscr5    ,iscr6    ,iopen    ,isym
 COMMON /zzzzzz/  z(1)
 EQUIVALENCE      (mcbs(1),ifila(1)) ,(mcbs(8),itypal)   ,  &
     (mcbs(9),alpha(1)) ,(mcbs(13),ifilb(1)),  &
     (mcbs(20),itypbt)  ,(mcbs(21),beta(1)) , (mcbs(61),ifilc(1))
 DATA    NAME  /  4HINIT,4HL2  /
 
 nomat   = 2
 iprec   = rdp
 alpha(2)= 0.0D0
 beta(2) = 0.0D0
 nx      = korsz(z) - offset
 nz      = nx
 
!     FORM AND DECOMPOSE THE LEFT HAND MATRIX
 
 itypal   = rdp
 itypbt   = rdp
 alpha(1) = 1.0D0/deltt**2
 beta(1)  = 0.5D0/deltt
 ifilc(4) = 6
 DO  i = 1,7
   ifila(i) = filem(i)
   ifilb(i) = fileb(i)
 END DO
 ifilc(2) = filek(2)
 ifilc(1) = iscr2
 IF (filek(1) <= 0) ifilc(1) = iscr1
 ifilc(3) = filek(2)
 IF (ifila(1) /= 0 .AND. ifila(4) /= 6) ifilc(4) = sqr
 IF (ifilb(1) /= 0 .AND. ifilb(4) /= 6) ifilc(4) = sqr
 ifilc(5) = iprec
 IF (filem(1) <= 0 .AND. fileb(1) <= 0) GO TO 60
 CALL sadd (z,z)
 IF (filek(1) <= 0) GO TO 21
 11 DO  i = 1,7
   ifila(i) = ifilc(i)
   ifilb(i) = filek(i)
 END DO
 IF (ifilb(4) /= 6) ifilc(4) = sqr
 ifilc(1) = iscr1
 alpha(1) = 1.0D0
 beta(1)  = 1.0D0/3.0D0
 CALL sadd (z,z)
 21 CONTINUE
 CALL wrttrl (ifilc)
 IF (ifilc(4) /= 6) GO TO 31
 
!     SET UP FOR SYMMETRIC DECOMPOSITION
 
 DO   i = 1,7
   ifa(i) = ifilc(i)
 END DO
 ifl(1) = iscr2
 ifu(1) = iscr3
 isc1   = iscr4
 isc2   = iscr5
 isc3   = iscr6
 ifl(5) = iprec
 ichl   = 0
 nxx    = nx
 FILE   = ifa(1)
 CALL sdcomp (*1030,z,z,z)
 CALL wrttrl (ifl)
 isym   = 0
 GO TO 33
 
!     SET UP FOR UNSYMMETRIC DECOMPOSITION
 
 31 CONTINUE
 isym = 1
 DO  i = 1,7
   ia(i)  = ifilc(i)
 END DO
 il(1)  = iscr2
 iu(1)  = iscr3
 iscr10 = iscr4
 iscr20 = iscr5
 iscr30 = iscr6
 il(5)  = iprec
 FILE   = ia(1)
 CALL decomp (*1030,z(1),z(1),z(1))
 CALL wrttrl (il)
 CALL wrttrl (iu)
 
!     FORM FIRST RIGHT HAND MATRIX
 
 33 CONTINUE
 DO  i = 1,7
   ifila(i) = filem(i)
 END DO
 alpha(1) = 2.0D0/deltt**2
 beta(1)  = -1.0D0/3.0D0
 ifilc(1) = iscr1
 CALL sadd (z,z)
 
!     FORM SECOND RIGHT HAND MATRIX
 
 alpha(1) = -1.0D0/deltt**2
 ifilc(1) = iscr5
 CALL sadd (z,z)
 DO  i = 1,7
   ifila(i) = ifilc(i)
   ifilb(i) = fileb(i)
 END DO
 alpha(1) = 1.0D0
 beta(1)  = 0.5D0/deltt
 ifilc(1) = iscr4
 CALL sadd (z,z)
 RETURN
 
!     ERRORS
 
 1030 ip1 = -5
 1031 CALL mesage (ip1,FILE,NAME(1))
 
!     NO BDD OR MDD
 
 60 IF (filek(1) <= 0) GO TO 70
 ifilc(1) =0
 GO TO 11
 
!     ILLEGAL INPUT.   NO MATRICES
 
 70 ip1 = -7
 GO TO 1031
END SUBROUTINE initl2
