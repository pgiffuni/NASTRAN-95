SUBROUTINE rowdyz (nfb,nlb,row,ntzys,d,dx,dy,dz,beta,idzdy,ntape,  &
        sgr,cgr,iprnt,yb,zb,ar,nsbe,xis1,xis2,a0)
     
!     CALCULATE A ROW OF DZ OR DY
 
!     SLENDER BODY
 
!     NFB       FIRST BODY OF THE DESIRED ORIENTATION - Z OR Y -
!     NLB       LAST  BODY OF THE DESIRED ORIENTATION
!     ROW       ROW OF DZ OR DY BEING CALCULATED
!     NTZYS     NO. COLUMNS TO BE CALCULATED
!     D         CALCULATED ROW
!     DX        X - COORD. OF RECEIVING POINT
!     DY        Y - COORD. OF RECEIVING POINT
!     DZ        Z - COORD. OF RECEIVING POINT
!     BETA      EQUALS SQRT(1-M**2)
!     MACH      MACH NO., M
!     IDZDY     FLAG REQUIRED FOR FLLD
 
 
 INTEGER, INTENT(IN)                      :: nfb
 INTEGER, INTENT(IN)                      :: nlb
 INTEGER, INTENT(IN OUT)                  :: row
 INTEGER, INTENT(IN OUT)                  :: ntzys
 REAL, INTENT(OUT)                        :: d(2,ntzys)
 REAL, INTENT(IN OUT)                     :: dx
 REAL, INTENT(IN OUT)                     :: dy
 REAL, INTENT(IN)                         :: dz
 REAL, INTENT(IN OUT)                     :: beta
 INTEGER, INTENT(IN)                      :: idzdy
 INTEGER, INTENT(IN OUT)                  :: ntape
 REAL, INTENT(IN OUT)                     :: sgr
 REAL, INTENT(IN OUT)                     :: cgr
 INTEGER, INTENT(IN OUT)                  :: iprnt
 REAL, INTENT(IN)                         :: yb(1)
 REAL, INTENT(IN)                         :: zb(1)
 REAL, INTENT(IN)                         :: ar(1)
 INTEGER, INTENT(IN)                      :: nsbe(1)
 REAL, INTENT(IN)                         :: xis1(1)
 REAL, INTENT(IN)                         :: xis2(1)
 REAL, INTENT(IN)                         :: a0(1)
 INTEGER :: b1,t1,b,t
 REAL :: kr
 
 
 COMMON /amgmn / mcb(7),nrow,nd,NE,refc,fmach,kr
 COMMON /system/ n,npot
 
 delta  = FLOAT(nd)
 epslon = FLOAT(NE)
 b1     =  0
 t1     =  0
 it1    = 0
 IF (nfb == 1 .OR. idzdy == 0) GO TO 10
 b   = nfb - 1
 DO  t = 1,b
   it1 = it1 + nsbe(t)
 END DO
 10 CONTINUE
 DO  b = nfb,nlb
   b1    = b1 + 1
   dar   = ar(b)
   nsbeb = nsbe(b)
   IF (iprnt /= 0) WRITE (npot,15) b,dy,yb(b),dz,zb(b)
   15 FORMAT (12H rowdyz  b =,i10,4E20.8)
   
!     LOOP FOR EACH ELEMENT IN BODY -B-
   
   DO  t = 1,nsbeb
     t1  = t1  + 1
     it1 = it1 + 1
     d(1,t1) = 0.0
     d(2,t1) = 0.0
     xi1  = xis1(it1)
     xi2  = xis2(it1)
     azro = a0(it1)
     eta  = yb(b)
     zeta = zb(b)
     
!     CHECK TO SEE IF CALCULATIONS ARE TO BE MADE
     
     IF (dy == eta .AND. dz == zeta) GO TO 30
     ASSIGN 20 TO jdzdy
     lhs = 0
     GO TO 100
     20 d(1,t1) = dzyr
     d(2,t1) = dzyi
     
!     SKIP IF NO SYMMETRY
     
     30 CONTINUE
     IF (delta == 0.0) GO TO 70
     eta = -yb(b)
     
!     CHECK TO SEE IF CALCULATIONS ARE TO BE MADE
     
     IF (dy == eta .AND. dz == zeta) GO TO 50
     lhs = 1
     ASSIGN 40 TO jdzdy
     GO TO 100
     40 d(1,t1) = d(1,t1) + delta*dzyr
     d(2,t1) = d(2,t1) + delta*dzyi
     50 CONTINUE
     IF (epslon == 0.0) CYCLE
     
!     CALC. ONLY IF DELTA AND EPSLON  NOT EQUAL ZERO
     
     eta  = -yb(b)
     zeta = -zb(b)
     
!     CHECK TO SEE IF CALCULATIONS ARE TO BE MADE
     
     IF (dy == eta .AND. dz == zeta) GO TO 70
     ASSIGN 60 TO jdzdy
     GO TO 100
     60 d(1,t1) = d(1,t1) + epslon*delta*dzyr
     d(2,t1) = d(2,t1) + epslon*delta*dzyi
     
!     SKIP IF NO GROUND EFFECTS
     
     70 IF (epslon == 0.0) CYCLE
     eta  =  yb(b)
     zeta = -zb(b)
     
!     CHECK TO SEE IF CALCULATIONS ARE TO BE MADE
     
     IF (dy == eta .AND. dz == zeta) CYCLE
     lhs = 1
     ASSIGN 80 TO jdzdy
     GO TO 100
     80 d(1,t1) = d(1,t1) + epslon*dzyr
     d(2,t1) = d(2,t1) + epslon*dzyi
     CYCLE
     
!     CALL SEQUENCE TO DZY
     
     100 CALL dzy (dx,dy,dz,sgr,cgr,xi1,xi2,eta,zeta,dar,azro,kr,refc,  &
         beta,fmach,lhs,idzdy,dzyr,dzyi)
     lhs = 0
     GO TO jdzdy, (20,40,60,80)
     
   END DO
   
!     END OF LOOP FOR ELEMENT
   
!     200 IS END OF LOOP ON SLENDER BODY
   
 END DO
 
!     WRITE ROW ON TAPE, ROW NUMBER, NO. ELEMENTS, DATA
 
 CALL WRITE (ntape,d,2*t1,0)
 IF (iprnt /= 0) WRITE (npot,210) row,t1,d
 210 FORMAT (' ROWDYZ - ROW NO.',i5,1H,,i10,' ELEMENTS',/(1X,6E12.4))
 RETURN
END SUBROUTINE rowdyz
