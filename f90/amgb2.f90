SUBROUTINE amgb2 (INPUT,d1jk,d2jk)
     
!     DRIVER FOR COMPRESSOR BLADE THEORY
!     COMPUTATIONS ARE FOR D1JK AND D2JK MATRICES
!     FOR COMPRESSOR BLADES K-SET = J-SET = NLINES*NSTNS
!     D1JK = F(INVERSE)TRANSPOSE
!     NOTE - AMP MODULE TAKES D1JK(TRANSPOSE) SO OUTPUT
!            F(INVERSE)TRANSPOSE TO GET EFFECT OF F(INVERSE) IN AMP.
 
!     D2JK = NULL
 
 
 INTEGER, INTENT(IN OUT)                  :: INPUT
 INTEGER, INTENT(IN OUT)                  :: d1jk
 INTEGER, INTENT(IN OUT)                  :: d2jk
 LOGICAL :: tsonic,debug
 INTEGER :: td1jk,td2jk,ecore,sysbuf,NAME(2),sln
 REAL :: minmac,maxmac,mach
 DIMENSION       iz(1)
 COMMON /amgp2 / td1jk(7),td2jk(7)
 COMMON /amgmn / mcb(7),nrow,dum(2),refc,sigma,rfreq
 COMMON /condas/ pi,twopi,radeg,degra,s4pisq
 COMMON /packx / iti,ito,ii,nn,incr
 COMMON /system/ sysbuf,iout
 COMMON /zzzzzz/ work(1)
 COMMON /BLANK / nk,nj
 COMMON /bamg2l/ iref,minmac,maxmac,nlines,nstns,refstg,refcrd,  &
     refmac,refden,refvel,refflo,sln,nstnsx,stager,  &
     chord,radius,bspace,mach,den,vel,flowa,amach, redf,blspc,amachr,tsonic,xsign
 COMMON /amgbug/ debug
 EQUIVALENCE     (work(1),iz(1))
 DATA    NAME  / 4HAMGB,4H2   /
 
!     READ PARAMETERS IREF,MINMAC,MAXMAC,NLINES AND NSTNS
 
 CALL fread (INPUT,iref,5,0)
 IF (debug) CALL bug1 ('ACPT-REF  ',5,iref,5)
 
!     READ REST OF ACPT RECORD INTO OPEN CORE AND LOCATE REFERENCE
!     PARAMETERS REFSTG,REFCRD,REFMAC,REFDEN,REFVEL AND REFFLO
 
 ecore = korsz(iz) - 3*sysbuf
 CALL READ (*10,*10,INPUT,iz,ecore,1,nwar)
 GO TO 120
 10 ndata = 3*nstns + 10
 IF (debug) CALL bug1 ('ACPT-REST ',10,iz,nwar)
 irsln = 0
 nline = 0
 DO  i = 1,nwar,ndata
   IF (iref == iz(i)) irsln = i
   nline = nline + 1
 END DO
 
!     DETERMINE DIRECTION OF BLADE ROTATION VIA Y-COORDINATES AT TIP
!     STREAMLINE. USE COORDINATES OF FIRST 2 NODES ON STREAMLINE.
 
 iptr  = ndata*(nlines-1)
 xsign = 1.0
 IF (work(iptr+15) < work(iptr+12)) xsign = -1.0
 
 IF (debug) CALL bug1 ('RADII     ',20,radii,nlines)
 
!     DID IREF MATCH AN SLN OR IS THE DEFAULT TO BE TAKEN (BLADE TIP)
 
 IF (irsln == 0) irsln = (nlines-1)*ndata + 1
 refstg = work(irsln+2)
 refcrd = work(irsln+3)
 refmac = work(irsln+6)
 refden = work(irsln+7)
 refvel = work(irsln+8)
 refflo = work(irsln+9)
 
!     REPOSITION ACPT TO BEGINNING OF COMPRESSOR BLADE DATA
 
 CALL bckrec (INPUT)
 CALL fread (INPUT,0,-6,0)
 
 IF (debug) CALL bug1 ('BAMG2L    ',22,iref,27)
 
!     COMPUTE POINTERS AND SEE IF THERE IS ENOUGH CORE
 
 nsns = nstns*nstns
 ip1  = 1
 ip2  = ip1 + nsns
 next = ip2 + 3*nstns
 IF (next > ecore) GO TO 120
 
!     COMPUTE F(INVERSE) FOR EACH STREAMLINE
 
 nn = ii + nstns - 1
 DO  nline = 1,nlines
   CALL amgb2a (INPUT,work(ip1),work(ip2),work(ip2))
   
!     OUTPUT D1JK (=F(INVERSE)TRANSPOSE) FOR THIS STREAMLINE.
!     NOTE - AMP MODULE TAKES D1JK(TRANSPOSE) SO OUTPUT
!            F(INVERSE)TRANSPOSE TO GET EFFECT OF F(INVERSE) IN AMP.
   
   ip3 = ip2 + nstns - 1
   DO  i = 1,nstns
     k   = i
     DO  j = ip2,ip3
       work(j) = work(k)
       k   = k + nstns
     END DO
     CALL pack (work(ip2),d1jk,td1jk)
     IF (debug) CALL bug1 ('D1JK      ',31,work(ip2),nstns)
   END DO
   ii  = ii +nstns
   IF (nline == nlines) CYCLE
   nn  = nn + nstns
 END DO
 
!     OUTPUT D2JK = NULL
 
 DO  icol = 1,nk
   CALL bldpk (iti,ito,d2jk,0,0)
   CALL bldpkn (d2jk,0,td2jk)
 END DO
 RETURN
 
!     ERROR MESSAGES
 
!     NOT ENOUGH CORE
 
 120 CALL mesage (-8,0,NAME)
 RETURN
END SUBROUTINE amgb2
