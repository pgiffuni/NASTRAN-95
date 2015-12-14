SUBROUTINE amgt1d (ajj,tsonx,tamach,tredf,nstns2)
     
!     TRANSONIC INTERPOLATION CODE FOR SWEPT TURBOPROPS.
 
 
 COMPLEX, INTENT(IN OUT)                  :: ajj(nstns2,1)
 INTEGER, INTENT(IN OUT)                  :: tsonx(1)
 REAL, INTENT(IN OUT)                     :: tamach(1)
 REAL, INTENT(IN OUT)                     :: tredf(1)
 INTEGER, INTENT(IN)                      :: nstns2
 INTEGER :: sln
 
 
 
 
 
 
 COMMON /tamg1l/iref,minmac,maxmac,nlines,nstns,refstg,refcrd,  &
     refmac,refden,refvel,refswp,sln,nstnsx,stag,  &
     chord,dcbdzb,bspace,mach,den,vel,sweep,amach, redf,blspc,amachr,tsonic
 COMMON /amgmn /mcb(7),nrow,dum(2),refc,sigma,rfreq
 
 numm = 2 * nstns2 * nstns2
 DO  nline = 1,nlines
   IF(tsonx(nline) == 0) CYCLE
   ns = 0
   IF(nline == 1) GO TO 90
   IF(tamach(nline) >= 1.0) GO TO 20
!       SUBSONIC
   IF(nline . EQ.2) nline1=1
   IF(nline . EQ.2) GO TO 93
   17 nline1 = nline -2
   nline2 = nline -1
   GO TO 70
!        SUPERSONIC
   20 IF( nline == nlines) GO TO 17
   ns =1
   GO TO 90
   30 IF(nline1 == 0) GO TO 17
   IF(nline2 /= 0) GO TO 70
   nline2 = nline1
   nline1 = nline-1
   70 CALL intert(nline,nline1,nline2,numm,ajj,tamach)
   CYCLE
!       SEARCH FOR 1ST--2--KNOWN STREAMLINES
   90 nline1 = 0
   93 nline2 = 0
   nnline = nline + 1
   DO   i=nnline,nlines
     IF(nline2 /= 0)  EXIT
     IF(tsonx(i) /= 0) CYCLE
     IF(nline1 == 0)  nline1 = i
     IF(nline1 /= i)  nline2 = i
   END DO
   97 IF(ns == 0) GO TO 70
   GO TO 30
 END DO
 RETURN
END SUBROUTINE amgt1d
