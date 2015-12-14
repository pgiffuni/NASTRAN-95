SUBROUTINE timeeq (b,bbar,c,cbar,r,ientry,ncol,tim)
     
!     TIMEEQ SOLVES THE TIME AND CORE FUNCTIONS FOR DECOMP AND CDCOMP
 
 
 REAL, INTENT(IN)                         :: b
 REAL, INTENT(IN)                         :: bbar
 REAL, INTENT(IN)                         :: c
 REAL, INTENT(IN)                         :: cbar
 REAL, INTENT(IN)                         :: r
 INTEGER, INTENT(IN)                      :: ientry
 INTEGER, INTENT(IN OUT)                  :: ncol
 REAL, INTENT(OUT)                        :: tim
 INTEGER :: sysbuf
 REAL :: mb(1),mc(1),k1,k2,k3,k4,k5
 COMMON /ntime / lntime, tcons(15) /system/ ksystm(65)
 
 EQUIVALENCE     (ksystm( 1),sysbuf),(ksystm(40),nbpw),  &
     (ksystm(55),iprec ),(tcons (1) ,aaio),  &
     (tcons ( 2),aapak ),(tcons (8),mb(1)), (tcons (12),mc(1) )
 
 
 iret  = 0
 ientr = ientry
 1 amb   = mb(iprec)
 amc   = mc(iprec)
 IF (nbpw < 60) GO TO 2
 amb   = 3.0*amb
 amc   = 3.0*amc
 2 aio   = aaio
 apak  = aapak
 IF (ientr == 1) GO TO 10
 amb   = 5.*amb
 amc   = 5.*amc
 aio   = aio+aio
 apak  = 1.1*apak
 10 IF (iret == 1) GO TO 20
 tim = FLOAT(ncol)*(amb*bbar*r+amc*(bbar*c+bbar*cbar+b*cbar+  &
     2.0*c*cbar)+aio*bbar*(b+bbar-r-1.0))*1.e-06
 RETURN
 
 
 ENTRY tfin (ab,abbar,ac,acbar,ar,jentry,ancol,timex)
!     ====================================================
 
 iret  = 1
 ientr = jentry
 GO TO 1
 20 timex = 0.
 k1    = ancol - ab - abbar - abbar
 IF (k1 <= 0.) GO TO 30
 timex = k1*(amb*abbar*ar+aio*abbar*(ab+abbar-ar)+apak*(ab+abbar* 2.))
 30 k2  = ab + abbar
 k3  = k2
 IF (ancol >= ab+abbar+abbar) GO TO 35
 k2  = ancol - abbar
 k3  = ab + abbar
 IF (ancol < ab+abbar) k3 = ancol
 35 timex = timex+.5*k2*(abbar*k2*amb+(k3-ar)*(aio-amb)*abbar+  &
     2.*apak*abbar+apak*k2)
 IF (ancol < ab+abbar+abbar) GO TO 40
 k4 = ab + abbar - ar
 k5 = ab + 1.5*abbar
 IF (ab > ar) k4 = abbar
 GO TO 50
 40 k4 = ancol - ar
 k5 = ancol
 IF (ancol-ar > abbar) k4 = abbar
 50 timex = timex+abbar**3/3.*amb+k4**3*.5*aio+apak*abbar*k5
 timex = (timex+(ancol-abbar)*(amc*(abbar*ac+ab*acbar+abbar*acbar+  &
     ac*acbar)+apak*(ac+acbar)))*1.e-06
 RETURN
 
 
 ENTRY rcore (ib,ibbar,ic,icbar,incol,kentry,nx,ir)
!     ==================================================
!     ENTRY FOR THE CORE FUNCTION
 
 ir = (nx-((ib+ibbar+1) +2*kentry*MIN0(incol,ib+ibbar+ibbar)+  &
     2*kentry*ic*(ibbar+2)+2*icbar*kentry*(MIN0(ib+ibbar,incol)+1)  &
     +2*kentry*ic*icbar +ic+icbar*kentry+icbar)-6*sysbuf)/ (2*kentry*ibbar)
 RETURN
END SUBROUTINE timeeq
