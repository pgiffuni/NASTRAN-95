SUBROUTINE amgb1s (INPUT,fmat,xyzb,INDEX,radii,wfact,nline)
     
!     COMPUTE F(INVERSE) AND WFACT FOR THIS STREAMLINE
 
 
 INTEGER, INTENT(IN OUT)                  :: INPUT
 REAL, INTENT(OUT)                        :: fmat(nstns,nstns)
 REAL, INTENT(IN)                         :: xyzb(3,nstns)
 INTEGER, INTENT(IN OUT)                  :: INDEX(1)
 REAL, INTENT(IN OUT)                     :: radii(1)
 REAL, INTENT(OUT)                        :: wfact
 INTEGER, INTENT(IN)                      :: nline
 LOGICAL :: tsonic,debug
 INTEGER :: sln
 REAL :: minmac,maxmac,mach
 DIMENSION  tbl(3,3)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /system/ sysbuf,iout
 COMMON /amgmn / mcb(7),nrow,dum(2),refc,sigma,rfreq
 COMMON /condas/ pi,twopi,radeg,degra,s4piso
 COMMON /bamg1l/ iref,minmac,maxmac,nlines,nstns,refstg,refcrd,  &
     refmac,refden,refvel,refflo,sln,nstnsx,stager,  &
     chord,radius,bspace,mach,den,vel,flowa,amach, redf,blspc,amachr,tsonic,xsign
 COMMON /amgbug/ debug
 
!     READ STREAMLINE DATA
 
 nstns3 = 3*nstns
 CALL fread (INPUT,sln,10,0)
 CALL fread (INPUT,xyzb,nstns3,0)
 IF (debug) CALL bug1 ('ACPT-SLN  ',11,sln,10)
 IF (debug) CALL bug1 ('XYZB      ',12,xyzb,nstns3)
 
!     COMPUTE PARAMETERS
 
 amach  = mach  *COS(degra*(flowa-stager))
 amachr = refmac*COS(degra*(refflo-refstg))
 IF (debug) CALL bug1 ('AMACH     ',13,amach,1)
 IF (debug) CALL bug1 ('AMACHR    ',14,amachr,1)
 
!     (1) COMPUTE BASIC TO LOCAL TRANSFORMATION
!         XYZB ARRAY CONTAINS X,Y,Z COORDINATES IN BASIC SYSTEM
!         FOR ALL NODES ON THE STREAMLINE LEADING EDGE TO TRAILING EDGE
!     (2) TRANSFORM BASIC X,Y,Z ON STREAMLINE TO LOCAL X,Y,Z-S
!     (3) COMPUTE FMAT(NSTNS X NSTNS)
!     (4) COMPUTE FMAT(INVERS) - USE -
!         CALL INVERS(NSTNS,FMAT,NSTNS,DUM1,0,DETERM,ISING,INDEX)
!     (5) COMPUTE WFACT FOR THIS STREAMLINE - NOTE - ALL RADIUS HAVE
!         ALREADY BEEN STORED IN ARRAY RADII FOR ALL STREAMLINES
 
 xa =  xyzb(1,1)
 ya =  xyzb(2,1)
 za =  xyzb(3,1)
 xb =  xyzb(1,nstns)
 yb =  xyzb(2,nstns)
 zb =  xyzb(3,nstns)
 
!     EVALUATE  TBL  ROW 2
 
 xba = xb - xa
 yba = yb - ya
 zba = zb - za
 al2sq = xba**2 + yba**2
 al2 = SQRT(al2sq)
 al1sq = al2sq  + zba**2
 al1=  SQRT(al1sq)
 tbl(2,1) =-xsign*(yba/al2)
 tbl(2,2) = xsign*(xba/al2)
 tbl(2,3) = 0.0
 
!     EVAL  TBL  ROW 1
 
 tbl(1,1) = xba/al1
 tbl(1,2) = yba/al1
 tbl(1,3) = zba/al1
 
!     EVALUATE TBL  ROW 3
 
 tbl(3,1) =-tbl(1,3)*(xba/al2)
 tbl(3,2) =-tbl(1,3)*(yba/al2)
 tbl(3,3) = al2/al1
 fmat(1,1)= 1.0
 pic =  pi/chord
 ch2 = 2.0/chord
 DO  i = 2,nstns
   x = tbl(1,1)*(xyzb(1,i)-xyzb(1,1)) + tbl(1,2)*(xyzb(2,i)-xyzb(2,1))  &
       + tbl(1,3)*(xyzb(3,i)-xyzb(3,1))
   fmat(1,i) = 0.0
   fmat(i,1) = 1.0
   fmat(i,2) = ch2*x
   DO  j = 3,nstns
     an  = j - 2
     arg = pic*an*x
     fmat(i,j) = SIN(arg)
   END DO
 END DO
 IF (debug) CALL bug1 ('FMAT      ',50,fmat,nstns*nstns)
 ising = -1
 CALL invers (nstns,fmat,nstns,dum1,0,determ,ising,INDEX)
 IF (debug) CALL bug1 ('FMAT-INV  ',60,fmat,nstns*nstns)
 IF (ising == 2) GO TO 80
 k = nline + 1
 l = k - 2
 IF (nline ==      1) l = 1
 IF (nline == nlines) k = nlines
 
!     COMPUT WFACT FOR THIS STREAMLINE
 
 wfact = (den/refden)*(vel/refvel)**2 * ((amach*refmac)/(mach*amachr))**2 *  &
     (radii(k) - radii(l))*0.5
 IF (debug) CALL bug1 ('WFACT     ',70,wfact,1)
 RETURN
 
!     ERROR MESSAGE, SINGULAR MATRIX
 
 80 WRITE  (iout,90) ufm,sln
 90 FORMAT (a23,' -AMG MODULE- SINGULAR MATRIX IN ROUTINE AMGB1S FOR',  &
     ' STREAML2, SLN =',i3, /39X,'CHECK STREAML2 BULK DATA CARD.')
 CALL mesage (-61,0,0)
 RETURN
END SUBROUTINE amgb1s
