SUBROUTINE amgt1 (INPUT,matout,skj)
     
!     DRIVER FOR SWEPT TURBOPROP BLADES (AEROELASTIC THEORY 7).
 
!     COMPUTATIONS ARE FOR THE AJJL AND SKJ MATRICES.
!     FOR SWEPT TURBOPROPS K-SET = J-SET = 2*NSTNS*NLINES.
!     SKJ = F(INVERS)TRANSPOSE.
 
 
 INTEGER, INTENT(IN OUT)                  :: INPUT
 INTEGER, INTENT(IN OUT)                  :: matout
 INTEGER, INTENT(IN OUT)                  :: skj
 LOGICAL :: tsonic,debug
 INTEGER :: ecore,sysbuf,iz(1),NAME(2),sln, tskj
 REAL :: minmac,maxmac,mach
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /amgmn / mcb(7),nrow,dum(2),refc,sigma,rfreq,tskj(7),isk, nsk
 COMMON /condas/ pi,twopi,radeg,degra,s4pisq
 COMMON /tamg1l/ iref,minmac,maxmac,nlines,nstns,refstg,refcrd,  &
     refmac,refden,refvel,refswp,sln,nstnsx,stager,  &
     chord,dcbdzb,bspace,mach,den,vel,sweep,amach, redf,blspc,amachr,tsonic,xsign
 COMMON /amgbug/ debug
 COMMON /zzzzzz/ work(1)
 COMMON /system/ sysbuf,iout
 COMMON /packx / iti,ito,ii,nn,incr
 COMMON /BLANK / nk,nj
 EQUIVALENCE     (work(1),iz(1))
 DATA    NAME  / 4HAMGT,4H1   /
 
!     READ PARAMETERS IREF,MINMAC,MAXMAC,NLINES AND NSTNS
 
 CALL READ (*999,*999,INPUT,iref,5,0,n)
 IF (debug) CALL bug1 ('ACPT-REF  ',5,iref,5)
 
!     READ REST OF ACPT RECORD INTO OPEN CORE AND LOCATE REFERENCE
!     PARAMETERS REFSTG,REFCRD,REFMAC,REFDEN,REFVEL AND REFSWP
 
 ecore = korsz(iz) - 4*sysbuf
 CALL READ (*10,*10,INPUT,iz,ecore,1,nwar)
 GO TO 998
 10 irsln = 0
 IF (debug) CALL bug1 ('ACPT-REST ',10,iz,nwar)
 ntsonx= 0
 ndata = 3*nstns + 10
 nline = 0
 DO  i = 1,nwar,ndata
   
!     LOCATE REFERENCE STREAMLINE NUMBER (IREF = SLN)
   
   IF (iref == iz(i)) irsln = i
   
!     STORE MACH NUMBERS FOR LATER DATA CHECK.
   
   mach = work(i+6)
   IF (mach > maxmac .AND. mach < minmac) ntsonx = ntsonx + 1
   nline = nline + 1
   work(nwar+nline) = mach
 END DO
 
!     DETERMINE DIRECTION OF BLADE ROTATION VIA Y-COORDINATES AT TIP
!     STREAMLINE. USE COORDINATES OF FIRST 2 NODES ON STREAMLINE.
 
 iptr  = ndata*(nlines-1)
 xsign = 1.0
 IF (work(iptr+15) < work(iptr+12)) xsign = -1.0
 
!     INPUT CHECKS -
 
!     (1) MACH NUMBERS MUST INCREASE FROM BLADE ROOT TO BLADE TIP.
!         NOTE - THIS CHECK WILL NOT BE MADE FOR SWEPT TURBOPROPS.
!     (2) SUPERSONIC CASCADE CODE HAS BEEN INSTALLED IN SUB.AMGT1C
!     (3) LINEAR INTERPOLATION EXISTS FOR TRANSONIC STREAMLINES
!     (4) ALL TRANSONIC STREAMLINES ARE NEVER ALLOWED.
 
 
!     CHECK FOR ALL TRANSONIC STREAMLINES.
 
 ibad = 0
 IF (ntsonx < nlines ) GO TO 30
 ibad = 1
 WRITE (iout,1001) ufm
 30 CONTINUE
 
!     MACH NUMBERS MUST INCREASE FROM BLADE ROOT TO BLADE TIP.
 
!     NOTE - THIS CHECK WILL NOT BE MADE FOR SWEPT TURBOPROPS.
 
 IF (ibad /= 0) GO TO 997
 
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
 refswp = work(irsln+9)
 
!     REPOSITION ACPT TO BEGINNING OF BLADE DATA.
 
 CALL bckrec (INPUT)
 CALL fread (INPUT,0,-6,0)
 IF (debug) CALL bug1 ('TAMG1L    ',46,iref,26)
 
!     COMPUTE POINTERS AND SEE IF THERE IS ENOUGH CORE.
!     IP1 AND IP2 ARE COMPLEX POINTERS.
 
 nstns2 = 2*nstns
 najjc  = nstns2
 ntsonx = 1
 IF (tsonic) najjc  = nlines*nstns2
 IF (tsonic) ntsonx = nlines
 ip1  = 1
 ip2  = ip1 + 2*(nstns2*najjc)
 ip3  = ip2 + 1
 ip4  = ip3 + ntsonx
 ip5  = ip4 + ntsonx
 next = ip5 + ntsonx
 IF (next > ecore) GO TO 998
 
!     CALL ROUTINE TO COMPUTE AND OUTPUT AJJL.
 
 iti = 3
 ito = 3
 
 CALL amgt1a (INPUT,matout,work(ip1),work(ip3),work(ip4),work(ip5), nstns2)
 IF (debug) CALL bug1 ('AJJL      ',48,work(ip1),ip2-1)
 
!     COMPUTE F(INVERSE) FOR EACH STREAMLINE
 
!     COMPUTE POINTERS AND SEE IF THERE IS ENOUGH CORE
 
 nsns = nstns*nstns
 ip1  = 1
 ip2  = ip1 + nsns
 next = ip2 + 3*nstns
 IF (next > ecore) GO TO 998
 
!     REPOSITION ACPT TO BEGINNING OF BLADE DATA.
 
 CALL bckrec (INPUT)
 CALL fread (INPUT,0,-6,0)
 
 iti = 1
 ito = 3
 
 ii  = isk
 nsk = nsk + nstns
 nn  = nsk
 DO  nline = 1,nlines
   CALL amgt1s (INPUT,work(ip1),work(ip2),work(ip2))
   
!     OUTPUT SKJ (= F(INVERS)TRANSPOSE) FOR THIS STREAMLINE
   
   ip3 = ip2 + nstns - 1
   DO  i = 1,nstns
     k = i
     DO  j = ip2,ip3
       work(j) = work(k)
       k = k + nstns
     END DO
     CALL pack (work(ip2),skj,tskj)
     IF (debug) CALL bug1 ('SKJ       ',55,work(ip2),nstns)
   END DO
   ii = ii + nstns
   nn = nn + nstns
   DO  i = 1,nstns
     k  = i
     DO  j = ip2,ip3
       work(j) = work(k)
       k  = k + nstns
     END DO
     CALL pack (work(ip2),skj,tskj)
     IF (debug) CALL bug1 ('SKJ       ',75,work(ip2),nstns)
   END DO
   ii = ii + nstns
   IF (nline == nlines) CYCLE
   nn = nn + nstns
 END DO
 
!     UPDATE NROW AND PACK POINTERS
 
 nrow = nrow + nlines*nstns2
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
END SUBROUTINE amgt1
