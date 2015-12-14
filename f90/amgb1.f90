SUBROUTINE amgb1 (INPUT,matout,skj)
     
    !     DRIVER FOR COMPRESSOR BLADE THEORY.
    !     COMPUTATIONS ARE FOR THE AJJL AND SKJ MATRICES.
    !     FOR COMPRESSOR BLADES K-SET = J-SET = NLINES*NSTNS.
    !     SKJ = W*F(INVERS)TRANSPOSE.
 
    INTEGER, INTENT(IN OUT)                  :: INPUT
    INTEGER, INTENT(IN OUT)                  :: matout
    INTEGER, INTENT(IN OUT)                  :: skj
    LOGICAL :: tsonic,debug
    INTEGER :: ecore,sysbuf,iz(1),NAME(2),sln, tskj
    REAL :: minmac,maxmac,mach,radii(50)
    CHARACTER (LEN=23) :: ufm
    COMMON /xmssg / ufm
    COMMON /amgmn / mcb(7),nrow,dum(2),refc,sigma,rfreq,tskj(7),isk, nsk
    COMMON /condas/ pi,twopi,radeg,degra,s4pisq
    COMMON /bamg1l/ iref,minmac,maxmac,nlines,nstns,refstg,refcrd,  &
        refmac,refden,refvel,refflo,sln,nstnsx,stager,  &
        chord,radius,bspace,mach,den,vel,flowa,amach, redf,blspc,amachr,tsonic,xsign
    COMMON /zzzzzz/ work(1)
    COMMON /system/ sysbuf,iout
    COMMON /packx / iti,ito,ii,nn,incr
    COMMON /BLANK / nk,nj
    COMMON /amgbug/ debug
    EQUIVALENCE     (work(1),iz(1))
    DATA    NAME  / 4HAMGB,4H1   /
 
    !     READ PARAMETERS IREF,MINMAC,MAXMAC,NLINES AND NSTNS
 
    CALL READ (*999,*999,INPUT,iref,5,0,n)
    IF (debug) CALL bug1 ('ACPT-REF  ',5,iref,5)
 
    !     READ REST OF ACPT RECORD INTO OPEN CORE AND LOCATE REFERENCE
    !     PARAMETERS REFSTG,REFCRD,REFMAC,REFDEN,REFVEL AND REFFLO
    !     STORE STREAMLINE RADIUS FOR ALL STREAMLINES
 
    ecore = korsz(iz) - 4*sysbuf
    CALL READ (*10,*10,INPUT,iz,ecore,1,nwar)
    GO TO 998
10  irsln = 0
    IF (debug) CALL bug1 ('ACPT-REST ',10,iz,nwar)
    ntsonx= 0
    ndata = 3*nstns + 10
    nline = 0
    DO  i = 1,nwar,ndata
   
        !     LOCATE REFERENCE STREAMLINE NUMBER (IREF = SLN)
   
        IF (iref == iz(i)) irsln = i
   
        !     STORE AMACH FOR LATER DATA CHECK. COUNT TRANSONIC STREAMLINES
   
        amachl = work(i+6)*COS(degra*(work(i+9)-work(i+2)))
        IF (amachl > maxmac .AND. amachl < minmac) ntsonx = ntsonx + 1
        nline = nline + 1
        work(nwar+nline) = amachl
        radii(nline) = work(i+4)
    END DO
 
    !     DETERMINE DIRECTION OF BLADE ROTATION VIA Y-COORDINATES AT TIP
    !     STREAMLINE. USE COORDINATES OF FIRST 2 NODES ON STREAMLINE.
 
    iptr  = ndata*(nlines-1)
    xsign = 1.0
    IF (work(iptr+15) < work(iptr+12)) xsign = -1.0
 
    IF (debug) CALL bug1 ('RADII     ',25,radii,nlines)
 
    !     INPUT CHECKS -
    !     (1) AMACH MUST INCREASE FROM BLADE ROOT TO BLADE TIP
    !     (2) ALL TRANSONIC AMACH-S ARE NOT ALLOWED AT PRESENT
 
    ibad = 0
    IF (ntsonx < nlines) GO TO 30
    ibad = 1
    WRITE (iout,1001) ufm
30 CONTINUE
   nw1 = nwar + 1
   nw2 = nwar + nlines - 1
   DO  i = nw1,nw2
       IF (work(i) > work(i+1)) GO TO 40
   END DO
   GO TO 45
40 ibad = 1
   isln = (i-nwar-1)*ndata + 1
   WRITE (iout,1002) ufm,iz(isln)
45 IF (ibad /= 0) GO TO 997
 
   !     SET TSONIC IF THERE ARE ANY TRANSONIC STREAMLINES
 
   tsonic = .false.
   IF (ntsonx > 0) tsonic = .true.
 
   !     STORE REFERENCE PARAMETERS
   !     DID IREF MATCH AN SLN OR IS THE DEFAULT TO BE TAKEN  (BLADE TIP)
 
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
   IF (debug) CALL bug1 ('BAMG1L    ',47,iref,26)
 
   !     COMPUTE POINTERS AND SEE IF THERE IS ENOUGH CORE
   !     IP1 AND IP2 ARE COMPLEX POINTERS
 
   najjc  = nstns
   ntsonx = 1
   IF (tsonic) najjc  = nlines*nstns
   IF (tsonic) ntsonx = nlines
   ip1  = 1
   ip2  = ip1 + 2*(nstns*najjc)
   ip3  = ip2 + 2*nstns
   ip4  = ip3 + ntsonx
   ip5  = ip4 + ntsonx
   next = ip5 + ntsonx
   IF (next > ecore) GO TO 998
 
   !     CALL ROUTINE TO COMPUTE AND OUTPUT AJJL.
 
   iti = 3
   ito = 3
 
   CALL amgb1a (INPUT,matout,work(ip1),work(ip2),work(ip3), work(ip4),work(ip5))
   IF (debug) CALL bug1 ('AJJL      ',48,work(ip1),ip2-1)
 
   !     COMPUTE F(INVERSE) AND W(FACTOR) FOR EACH STREAMLINE
 
   !     COMPUTE POINTERS AND SEE IF THERE IS ENOUGH CORE
 
   nsns = nstns*nstns
   ip1  = 1
   ip2  = ip1 + nsns
   next = ip2 + 3*nstns
   IF (next > ecore) GO TO 998
 
   !     REPOSITION ACPT TO BEGINNING OF COMPRESSOR BLADE DATA
 
   CALL bckrec (INPUT)
   CALL fread  (INPUT,0,-6,0)
 
   iti = 1
   ito = 3
 
   ii  = isk
   nsk = nsk + nstns
   nn  = nsk
   DO  nline = 1,nlines
       CALL amgb1s (INPUT,work(ip1),work(ip2),work(ip2),radii,wfact, nline)
   
       !     OUTPUT SKJ (= WFACT*F(INVERS)TRANSPOSE) FOR THIS STREAMLINE
   
       ip3 = ip2 + nstns - 1
       DO  i = 1,nstns
           k  = i
           DO  j = ip2,ip3
               work(j) = work(k)*wfact
               k  = k + nstns
           END DO
           CALL pack (work(ip2),skj,tskj)
           IF (debug) CALL bug1 ('SKJ       ',55,work(ip2),nstns)
       END DO
       ii = ii + nstns
       IF (nline == nlines) CYCLE
       nn = nn + nstns
   END DO
 
   !     UPDATE NROW AND PACK POINTERS
 
   nrow = nrow + nlines*nstns
   IF (debug) CALL bug1 ('NEW-NROW  ',110,nrow,1)
   isk = ii
   nsk = nn
   RETURN
 
   !     ERROR MESSAGES
 
   !     BAD STREAMLINE DATA
 
997 CALL mesage (-61,0,0)
 
   !     NOT ENOUGH CORE
 
998 CALL mesage (-8,0,NAME)
 
   !     INPUT NOT POSITIONED PROPERLY OR INCORRECTLY WRITTEN
 
999 CALL mesage (-7,0,NAME)
   RETURN
 
1001 FORMAT (a23,' -AMG MODULE- ALL TRANSONIC STREAMLINES NOT ALLOWED',  &
       /39X,'CHECK MACH ON STREAML2 BULK DATA CARDS OR', /39X,  &
       'CHANGE PARAMETERS MINMACH AND MAXMACH.')
1002 FORMAT (a23,' -AMG MODULE- MACH NUMBERS MUST INCREASE FROM BLADE',  &
       ' ROOT TO BLADE TIP.', /39X, 'CHECK STREAML2 BULK DATA CARD WITH SLN =',i3)

END SUBROUTINE amgb1
