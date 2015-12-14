SUBROUTINE ssold2 (itype,ftemp)
     
!     PHASE TWO STRESS DATA RECOVERY FOR THE SOLID ELEMENTS
 
!     ITYPE = 1,2,3,OR4 CORRESPONDING TO THE TETRA,WEDGE,HEXA1,ORHEXA2
!             ELEMENTS
 
!     PHIOUT CONTAINS THE FOLLOWING WHERE N IS THE NUMBER OF CORNERS
 
!             ELEMENT ID
!             N SILS
!             T SUB 0
!             6 THERMAL STRESS COEFFICIENTS
!             N VOLUME RATIO COEFFICIENTS
!             N 6 BY 3 MATRICES RELATING STRESS TO DISPLACEMENTS
 
!  $MIXED_FORMATS
 
 
 INTEGER, INTENT(IN)                      :: itype
 REAL, INTENT(IN OUT)                     :: ftemp(8)
 INTEGER :: nphi(1),eject,ishd(7),typ(8),istyp(2)
 REAL :: frlast(2)
 COMMON /system/ ibfsz,nout,idm(9),line
 COMMON /zzzzzz/ z(1)
 COMMON /sdr2x4/ dummy(35),ivec,ivecn,ldtemp,deform
 COMMON /sdr2x7/ phiout(100),stres(100),forvec(25)
 COMMON /sdr2x8/ temp(6),factor,npts,npoint,ks,kbeta,sigma(9), ctmp(6),csig(7)
 COMMON /sdr2x9/ nchk,isub,ild,frtmei(2),twotop,fnchk
 EQUIVALENCE     (nphi(1),phiout(1)),(ishd(1),lsub),  &
     (ishd(2),lld),(ishd(6),frlast(1))
 DATA    typ   / 4HTETR,1HA, 4HWEDG,1HE, 4HHEXA,1H1, 4HHEXA,1H2  /
 DATA    lld   , lsub,frlast / 2*-1, -1.0E30, -1.0E30  /
 
 SELECT CASE ( itype )
   CASE (    1)
     GO TO 100
   CASE (    2)
     GO TO 110
   CASE (    3)
     GO TO 120
   CASE (    4)
     GO TO 120
 END SELECT
 
 100 npts = 4
 GO TO 130
 110 npts = 6
 GO TO 130
 120 npts = 8
 
 
 130 DO  i = 1,9
   csig(i)  = 0.0
   sigma(i) = 0.0
 END DO
 
!     LOOP ON GRID POINTS, DISPLACEMENT EFFECTS
 
 DO  n = 1,npts
   npoint = ivec + nphi(n+1) - 1
   ks =  18*n + 2*npts - 9
   CALL smmats (phiout(ks),6,3,0, z(npoint),3,1,0, temp,ctmp)
   DO  i = 1,6
     csig (i+1) =  csig(i+1) + ctmp(i)
     sigma(i+1) = sigma(i+1) + temp(i)
   END DO
   
!     TEMPERATURE EFFECTS
   
   IF (ldtemp == -1) CYCLE
   kbeta  = npts + n + 8
   factor = (ftemp(n) - phiout(npts+2))*phiout(kbeta)
   
   DO  i = 1,6
     kbeta = npts + i + 2
     sigma(i+1) = sigma(i+1) - phiout(kbeta)*factor
   END DO
 END DO
 sigma(1) = phiout(1)
 DO  i = 1,7
   stres(i) = sigma(i)
 END DO
 
!     OCTAHEDRAL STRESS AND HYDROSTATIC PRESSURE
 
 stres(8) = SQRT(sigma(2)*(sigma(2) - sigma(3) - sigma(4))*2.0 +  &
     2.0*sigma(3)*(sigma(3) - sigma(4)) + 2.0* sigma(4)**2 +  &
     6.0*(sigma(5)**2 + sigma(6)**2 + sigma(7)**2))/3.0
 stres(9) = -(sigma(2) + sigma(3) + sigma(4))/3.0
 IF (nchk <= 0) GO TO 450
 
!   . CHECK PRECISION
 
 csig(1) = phiout(1)
 k = 0
 
!   . STRESSES
 
 CALL sdrchk (sigma(2),csig(2),6,k)
 IF (k == 0) GO TO 450
 
!   . LIMITS EXCEEDED
 
 j = 2*itype
 istyp(1) = typ(j-1)
 istyp(2) = typ(j  )
 j = 0
 
 IF (lsub == isub .AND. frlast(1) == frtmei(1) .AND.  &
     lld == ild  .AND. frlast(2) == frtmei(2)) GO TO 420
 lsub = isub
 lld  = ild
 frlast(1) = frtmei(1)
 frlast(2) = frtmei(2)
 j = 2
 CALL page1
 400 CALL sd2rhd (ishd,j)
 line = line + 1
 WRITE  (nout,410)
 410 FORMAT (7X,4HTYPE,5X,3HEID,5X,2HSX,5X,2HSY,5X,2HSZ,4X,3HTYZ,4X,  &
     3HTXZ,4X,3HTXY)
 GO TO 430
 420 IF (eject(2) /= 0) GO TO 400
 430 WRITE  (nout,440,ERR=450) istyp,csig
 440 FORMAT (1H0,6X,a4,a1,i7,6F7.1)
 
 450 CONTINUE
 RETURN
END SUBROUTINE ssold2
