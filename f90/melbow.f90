SUBROUTINE melbow
     
!     THIS ROUTINE COMPUTES THE MASS MATRIX M(NPVT,NPVT) FOR AN ELBOW.
 
!     ECPT FOR THE ELBOW
 
!     ECPT( 1)  -  IELID        ELEMENT ID. NUMBER
!     ECPT( 2)  -  ISILNO(2)    * SCALAR INDEX NOS. OF THE GRID POINTS
!     ECPT( 3)  -    ...        *
!     ECPT( 9)  -  A            CROSS-SECTIONAL AREA
!     ECPT(13)  -  NSM          NON-STRUCTURAL MASS
!     ECPT(29)  -  R            RADIUS OF CURVATURE
!     ECPT(30)  -  BETAR        ANGLE FROM GA TO GB
 
 LOGICAL :: heat
 INTEGER :: iz(1),eor,clsrw,clsnrw,frowic,tnrows,outrw,bggind
 REAL :: nsm
 DOUBLE PRECISION :: ta(9),tb(9),dp(6),veci(3),dela(6),delb(6),fl,  &
     m(36),dumdp,fm
 DIMENSION        ecpt(9),iecpt(1)
 COMMON /sma2ht/  heat
 COMMON /sma2io/  ifcstm,ifmpt,ifdit,idum1,ifecpt,igecpt,ifgpct,  &
     iggpct,idum2,idum3,ifmgg,igmgg,ifbgg,igbgg,  &
     idum4,idum5,inrw,outrw,clsnrw,clsrw,neor,eor, mcbmgg(7),mcbbgg(7)
 COMMON /zzzzzz/  z(1)
 
!     SMA2 VARIABLE CORE BOOKKEEPING PARAMETERS
 
 COMMON /sma2bk/  icstm,ncstm,igpct,ngpct,ipoint,npoint,i6x6m,  &
     n6x6m,i6x6b,n6x6b
 
!     SMA2 PROGRAM CONTROL PARAMETERS
 
 COMMON /sma2cl/  iopt4,bggind,npvt,left,frowic,lrowic,nrowsc,  &
     tnrows,jmax,nlinks,link(10),nogo
 
!     ECPT COMMON BLOCK
 
 COMMON /sma2et/  ielid,isilno(2),smallv(3),icssv,imatid,a,i1,i2,  &
     fj,nsm,fe,dum(14),r,betar,dumm(8),tempel
 
!     SMA2 LOCAL VARIABLES
 
 COMMON /sma2dp/  ta,tb,dp,veci,dela,delb,fl,m,dumdp
 
!     INPUT AND OUTPUT BLOCKS FOR SUBROUTINE MAT
 
 COMMON /matin /  matidc,matflg,eltemp,stress,sinth,costh
 COMMON /hmtout/  cp
 COMMON /matout/  rho,prop(8)
 EQUIVALENCE      (z(1),iz(1),dz),(ecpt(1),iecpt(1),ielid)
 DATA    dcr   /  0.01745329/
 
!     COMPUTE LENGTH OF ELBOW, FL
 
 dp(1) = r
 dp(2) = betar
 dp(3) = dcr
 fl    = dp(1)*dp(2)*dp(3)
 IF (fl == 0.0D0) GO TO 200
 IF (heat) GO TO 300
 
!     GET RHO FROM MPT BY CALLING MAT
 
 matidc = imatid
 matflg = 4
 eltemp = tempel
 CALL mat (ecpt(1))
 DO  i = 1,36
   m(i) = 0.0D0
 END DO
 fm = 0.5*fl*(rho*a + nsm)
 
!     PUT MASS IN M-ARRAY
 
 m(1)  = fm
 m(8)  = m(1)
 m(15) = m(1)
 
!     INSERT THE 6 X 6
 
 CALL sma2b (m,npvt,-1,ifmgg,dumdp)
 RETURN
 
 200 CALL mesage (30,26,iecpt(1))
 
!     SET FLAG FOR FATAL ERROR WHILE ALLOWING ERROR MESSAGES TO
!     ACCUMULATE
 
 nogo = 1
 RETURN
 
!     HEAT FORMULATION
 
!     GET CP USING -HMAT- ROUTINE.
 
 300 matidc = imatid
 matflg = 4
 CALL hmat (ielid)
 m(1) = fl*DBLE(ecpt(9))*DBLE(cp)/2.0D0
 
!     OUTPUT THE MASS FOR HEAT PROBLEM.
 
 CALL sma2b (m(1),npvt,npvt,ifbgg,dumdp)
 RETURN
END SUBROUTINE melbow
