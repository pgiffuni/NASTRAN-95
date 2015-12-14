SUBROUTINE shctss (ierr,elid,pid,mid,tlam,tmean,tgrad,thetae,  &
        ftherm,epslnt,icore,core)
     
!     SINGLE PRECISION ROUTINE TO EVALUATE THERMAL STRAINS FOR COMPOSITE
!     SHELL ELEMENTS.
 
!     INPUT :
!           ELID   - ELEMENT ID
!           PID    - PROPERTY ID
!           MID    - ARRAY OF LAMINATE MATERIAL ID'S
!           TLAM   - LAMINATE THICKNESS
!           TMEAN  - ELEMENT MEAN TEMPERATURE
!           TGRAD  - THERMAL GRADIENT
!           THETAE - ANGLE FROM ELEMENT X-AXIS TO MATERIAL X-AXIS
!           FTHERM - ARRAY OF THERMAL FORCES CONTAINING THE USER-
!                    DEFINED THERMAL MOMENTS, IF SUPPLIED
!           IPCMPI AND NPCMPI ARE THE STARTING POINT AND THE NUMBER
!           OF WORDS OF PCOMPI DATA IN CORE, AS INPUT BY /SDR2C1/.
!     OUTPUT:
!           EPSLNT - ARRAY OF THERMAL STRAINS FOR THE LAMINATE
 
 
 
 INTEGER, INTENT(OUT)                     :: ierr
 INTEGER, INTENT(IN OUT)                  :: elid
 INTEGER, INTENT(IN OUT)                  :: pid
 INTEGER, INTENT(IN)                      :: mid(4)
 REAL, INTENT(IN)                         :: tlam
 REAL, INTENT(IN)                         :: tmean
 REAL, INTENT(IN)                         :: tgrad
 REAL, INTENT(IN)                         :: thetae
 REAL, INTENT(OUT)                        :: ftherm(6)
 REAL, INTENT(IN OUT)                     :: epslnt(6)
 INTEGER, INTENT(IN)                      :: icore(1)
 REAL, INTENT(IN)                         :: core(1)
 LOGICAL :: nonmem,pcmp,pcmp1,pcmp2
 INTEGER :: pcomp,pcomp1,pcomp2, pidloc,sym,symmem,indx(6,3)
 
 REAL :: minrt,abbd(6,6),stiff(36),glay(9),glayt(9),  &
     gbar(9),gprop(25),alphal(3),alphae(3),galpha(3),  &
     theta,transl(9),tsubo,delta,deltat,zk,zk1,zref,  &
     zsubi,c,c2,s,s2,pi,twopi,raddeg,degrad,determ, dum(6)
 COMMON /condas/  pi,twopi,raddeg,degrad
 COMMON /matin /  matid,inflag,eltemp
 COMMON /sdr2c1/  ipcmp,npcmp,ipcmp1,npcmp1,ipcmp2,npcmp2
 
 
 DATA    pcomp ,  pcomp1,pcomp2 / 0,1,2 /
 DATA    sym   ,  mem   ,symmem / 1,2,3 /
 
!     INITIALIZE
 
 ierr = 0
 DO  ll = 1,6
   DO  mm = 1,6
     abbd(ll,mm) = 0.0
   END DO
 END DO
 
 minrt = tlam*tlam*tlam/12.0
 zref  =-tlam/2.0
 
 inflag = 12
 eltemp = tmean
 
 itype = -1
 lpcomp= ipcmp + npcmp + npcmp1 + npcmp2
 pcmp  = npcmp  > 0
 pcmp1 = npcmp1 > 0
 pcmp2 = npcmp2 > 0
 
!     ISSUE ERROR IF PCOMPI DATA HAS NOT BEEN READ INTO CORE
 
 IF (lpcomp == ipcmp) GO TO 600
 
!     LOCATE PID BY PERFORMING A SEQUENTIAL SEARCH OF THE PCOMPI DATA
!     BLOCK WHICH IS IN CORE.
 
!     SEARCH FOR PID IN PCOMP DATA
 
 IF (.NOT.pcmp) GO TO 110
 ip = ipcmp
 IF (icore(ip) == pid) GO TO 210
 ipc11 = ipcmp1 - 1
 DO  ip = ipcmp,ipc11
   IF (icore(ip) /= -1 .OR. ip >= ipc11) CYCLE
   IF (icore(ip+1) == pid) GO TO 200
 END DO
 
!     SEARCH FOR PID IN PCOMP1 DATA
 
 110 IF (.NOT.pcmp1) GO TO 130
 ip = ipcmp1
 IF (icore(ip) == pid) GO TO 230
 ipc21 = ipcmp2 - 1
 DO  ip = ipcmp1,ipc21
   IF (icore(ip) /= -1 .OR. ip >= ipc21) CYCLE
   IF (icore(ip+1) == pid) GO TO 220
 END DO
 
!     SEARCH FOR PID IN PCOMP2 DATA
 
 130 IF (.NOT.pcmp2) GO TO 150
 ip = ipcmp2
 IF (icore(ip) == pid) GO TO 250
 lpc11 = lpcomp - 1
 DO  ip = ipcmp2,lpc11
   IF (icore(ip) /= -1 .OR. ip >= lpc11) CYCLE
   IF (icore(ip+1) == pid) GO TO 240
 END DO
 
!     PID WAS NOT LOCATED; ISSUE ERROR
 
 150 GO TO 600
 
!     PID WAS LOCATED; DETERMINE TYPE
 
 200 ip     = ip + 1
 210 itype  = pcomp
 pidloc = ip
 nlay   = icore(pidloc+1)
 ipoint = pidloc + 8 + 4*nlay
 GO TO 300
 
 220 ip     = ip + 1
 230 itype  = pcomp1
 pidloc = ip
 nlay   = icore(pidloc+1)
 ipoint = pidloc + 8 + nlay
 GO TO 300
 
 240 ip     = ip + 1
 250 itype  = pcomp2
 pidloc = ip
 nlay   = icore(pidloc+1)
 ipoint = pidloc + 8 + 2*nlay
 
 300 tsubo  = core(ipoint+24)
 delta  = tmean - tsubo
 lamopt = icore(pidloc+8)
 nonmem = lamopt /= mem .AND. lamopt /= symmem
 
!     LAMOPT - LAMINATION GENERATION OPTION
!            = ALL     (ALL PLYS, DEFAULT)
!            = SYM     (SYMMETRIC)
!            = MEM     (MEMBRANE ONLY)
!            = SYMMEM  (SYMMETRIC-MEMBRANE)
 
!     CONSTRUCT THE LAMINATE FORCE-STRAIN MATRIX
 
!     EXTENSIONAL
 
 matid = mid(1)
 CALL mat (elid)
 CALL lprops (gprop)
 
 DO  ll = 1,3
   ii = 3*(ll-1)
   DO  mm = 1,3
     abbd(ll,mm) = gprop(mm+ii)*tlam
   END DO
 END DO
 
!     BENDING
 
 IF (.NOT.nonmem) GO TO 400
 
 matid = mid(2)
 CALL mat (elid)
 CALL lprops (gprop)
 
 DO  ll = 1,3
   ii = 3*(ll-1)
   DO  mm = 1,3
     abbd(ll+3,mm+3) = gprop(mm+ii)*minrt
   END DO
 END DO
 
!     MEMBRANE-BENDING
 
 IF (lamopt == sym) GO TO 400
 
 matid = mid(4)
 CALL mat (elid)
 CALL lprops (gprop)
 
 DO  ll = 1,3
   ii = 3*(ll-1)
   DO  mm = 1,3
     abbd(ll,mm+3) = gprop(mm+ii)*tlam*tlam
     abbd(ll+3,mm) = gprop(mm+ii)*tlam*tlam
   END DO
 END DO
 
 
!     BEGIN THE LOOP OVER LAYERS
 
 400 zk = zref
 DO  k = 1,nlay
   
!     SET THE LAYER-DEPENDENT VARIABLES
   
   zk1 = zk
   IF (itype /= pcomp) GO TO 410
   zk = zk1 + core(pidloc+6+4*k)
   theta = core(pidloc   +7+4*k)
   GO TO 430
   
   410 IF (itype /= pcomp1) GO TO 420
   zk = zk1 + core(pidloc+8  )
   theta = core(pidloc   +8+k)
   GO TO 430
   
   420 IF (itype /= pcomp2) GO TO 430
   zk = zk1 + core(pidloc+7+2*k)
   theta = core(pidloc   +8+2*k)
   
!     LAYER MATERIAL PROPERTIES
   
   430 DO  ir = 1,9
     glay(ir) = core(ipoint+ir)
   END DO
   
   DO  ir = 1,3
     alphal(ir) = core(ipoint+13+ir)
   END DO
   
   ti = zk - zk1
   zsubi = (zk+zk1)/2.0
   deltat = delta + zsubi*tgrad
   
!     TRANSFORM THE LAYER MATERIAL PROPERTIES FROM THE FIBER SYSTEM TO
!     THE ELEMENT SYSTEM
   
   theta = theta*degrad + thetae
   c   = COS(theta)
   c2  = c*c
   s   = SIN(theta)
   s2  = s*s
   
   transl(1)  = c2
   transl(2)  = s2
   transl(3)  = c*s
   transl(4)  = s2
   transl(5)  = c2
   transl(6)  =-c*s
   transl(7)  =-2.0*c*s
   transl(8)  = 2.0*c*s
   transl(9)  = c2 - s2
   
!                _            T
!     CALCULATE [G] = [TRANSL] [GLAY][TRANSL]
   
   CALL gmmats (glay(1),3,3,0,  transl(1),3,3,0, glayt(1))
   CALL gmmats (transl(1),3,3,1, glayt(1),3,3,0, gbar(1))
   
!     CALCULATE [ALPHAE] = [TRANSL]X[ALPHA]
!     MODIFY [TRANSL] FOR TRANSFORMATIONS OF ALPHAS
   
   transl(3) = -transl(3)
   transl(6) = -transl(6)
   transl(7) = -transl(7)
   transl(8) = -transl(8)
   
   CALL gmmats (transl(1),3,3,0, alphal(1),3,1,0, alphae(1))
   
   
!     CALCULATE THERMAL FORCES AND MOMENTS
   
   CALL gmmats (gbar(1),3,3,0, alphae(1),3,1,0, galpha(1))
   
   DO  ir = 1,3
     ftherm(ir) = ftherm(ir  ) + galpha(ir)*deltat*(zk-zk1)
     IF (nonmem)  ftherm(ir+3) = ftherm(ir+3) - galpha(ir)*  &
         deltat*(zk*zk-zk1*zk1)/2.0
   END DO
   
!     CALCULATE CONTRIBUTION FROM SYMMETRIC LAYERS
   
   IF (lamopt /= sym .AND. lamopt /= symmem) GO TO 480
   deltat = delta - zsubi*tgrad
   
   DO  ir = 1,3
     ftherm(ir) = ftherm(ir  ) + galpha(ir)*deltat*(zk-zk1)
     IF (nonmem)  ftherm(ir+3) = ftherm(ir+3) - galpha(ir)*  &
         deltat*(zk1*zk1-zk*zk)/2.0
   END DO
   480 IF (itype == pcomp) ipoint = ipoint + 27
   
 END DO
 
 
!     END OF LOOP OVER THE LAYERS
 
!     COMPUTE THERMAL STRAIN VECTOR
 
!                      -1
!     {EPSLNT} = [ABBD]  {FTHERM}
 
 ising = -1
 CALL invers (6,abbd,6,dum,0,determ,ising,indx)
 
 DO  ll = 1,6
   nn = 6*(ll-1)
   DO  mm = 1,6
     stiff(nn+mm) = abbd(ll,mm)
   END DO
 END DO
 
 CALL gmmats (stiff(1),6,6,0, ftherm(1),6,1,0, epslnt(1))
 GO TO 700
 
 600 ierr = 1
 700 RETURN
END SUBROUTINE shctss
