SUBROUTINE cfactr (a,ll,ul,scr1,scr2,scr3,iopt)
     
 
 INTEGER, INTENT(IN)                      :: a
 INTEGER, INTENT(IN)                      :: ll
 INTEGER, INTENT(IN)                      :: ul
 INTEGER, INTENT(IN)                      :: scr1
 INTEGER, INTENT(IN)                      :: scr2
 INTEGER, INTENT(IN)                      :: scr3
 INTEGER, INTENT(OUT)                     :: iopt
 INTEGER :: fa        ,fl       ,fu       ,sr1      , sr2       ,sr3      ,  &
      mcb(7)    ,NAME(2)
 DOUBLE PRECISION :: det       ,mind
 COMMON   /cdcmpx/ fa(7)     ,fl(7)    ,fu(7)    ,sr1      ,  &
     sr2       ,sr3      ,det(2)   ,powr     ,  &
     nx        ,mind     ,ib       ,ibbar
 COMMON   /sfact / mfa(7)    ,mfl(7)   ,mfc(7)   ,m1fil    ,  &
     m2fil     ,mxx      ,d(5)     ,m3fil    , d1(2)     ,ichol
 COMMON   /sdccsp/ jfa(7)    ,jfl(7)   ,jfc(7)   ,j1fil    , j2fil     ,jx
 COMMON   /zzzzzz/ iz(1)
 DATA      NAME  / 4HCFAC,4HTR   /
 
 
 nz = korsz(iz)
 mcb(1) = a
 CALL rdtrl (mcb)
 IF (mcb(4) /= 6) GO TO 200
 
!     SYMMETRIC  COMPLEX
 
 DO  i = 1,7
   mfa(i) = mcb(i)
   mfl(i) = mcb(i)
   mfc(i) = mcb(i)
 END DO
 mfl(1) = ll
 mfc(1) = ul
 mfl(4) = 4
 mfc(4) = 5
 m1fil  = scr1
 m2fil  = scr2
 mxx    = nz
 m3fil  = scr3
 ichol  = 0
 CALL sdcomp (*900,iz,iz,iz)
 CALL wrttrl (mfl)
 iopt  = 2
 GO TO  60
 
!     UNSYMMETRIC  COMPLEX
 
 200 DO  i = 1,7
   fa(i) = mcb(i)
   fl(i) = mcb(i)
   fu(i) = mcb(i)
 END DO
 fl(1) = ll
 fu(1) = ul
 fl(4) = 4
 fu(4) = 5
 sr1   = scr1
 sr2   = scr2
 sr3   = scr3
 nx    = nz
!     IB    = 0
 
!     IF IB IS SET TO ZERO HERE, T08021 PRINTS 27 MORE MESSAGES 3027
!     AND 3028 FROM GENVEC WHICH IS CALLED BY CFACTR, WHCIH IS CALLED BY
!     FRD2C, IN FRRD2 MODULE
 
!IBMI 6/93
 ibbar = 0
 CALL cdcomp (*900,iz,iz,iz)
 CALL wrttrl (fu)
 CALL wrttrl (fl)
 iopt  = 1
 60 RETURN
 
!     ERRORS
 
 900 CALL mesage (-5,a,NAME)
 
 RETURN
END SUBROUTINE cfactr
