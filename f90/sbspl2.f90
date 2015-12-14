SUBROUTINE sbspl2( ntype, ti )
     
!     PHASE TWO STRESS DATA RECOVERY BASIC BENDING TRIANGLE.
 
!     NTYPE = 0  IMPLIES BASIC BENDING TRIANGLE
!     NTYPE = 3 IMPLIES TRI-PLATE IS CALLING
!     NTYPE = 4 IMPLIES QUAD-PLATE IS CALLING
 
 
 INTEGER, INTENT(IN)                      :: ntype
 REAL, INTENT(IN)                         :: ti(6)
 DIMENSION nsil(1), si(1)
 REAL :: sdelta(3),frlast(2)
 INTEGER :: eject    ,ished(7) ,tloads   ,bsc   ,plt   ,qd   ,tr  &
     ,       istyp(2)
 LOGICAL :: flag
 
 COMMON   /system/  ibfsz    ,nout     ,idm(9)   ,line
 COMMON /zzzzzz/ z(1)
 COMMON /sdr2x4/ dummy(35),ivec,dum11(11),tloads
 COMMON /sdr2x7/ est(100),stres(100),forvec(25)
 COMMON /sdr2x8/ eye,i,j,npoint,vec(5),zoveri,temp,delta  &
     ,            cvec(5),cfrvec(12)
 COMMON /sdr2x9/ nchk,isub,ild,frtmei(2),twotop,fnchk
 COMMON /sdr2de/ skp(8),ieltyp
 
 EQUIVALENCE (si(1),est(9)),(nsil(1),est(2))
 EQUIVALENCE (nelid,est(1))
 EQUIVALENCE (f1,n1)  ,  (ished(6),frlast(1) )  &
     ,           (ished(1),lsub)  ,  (ished(2),lld)
 
 DATA tr,qd,bsc,plt / 4H  tr, 4H  qd, 4HBSC , 4HPLT  /
 DATA lld,lsub,frlast / 2*-100, -1.0E30, -1.0E30 /
 
!     ******************************************************************
!  . ZERO OUT FORCE AND PRECISION CHECK VECTOR...
 DO  i = 1,6
   forvec(i) = 0.0E0
   cfrvec(i) = 0.0E0
   cfrvec(i+6) = 0.0E0
 END DO
 forvec(1) = est(1)
 
 npts = 3
 IF( ntype == 4 ) npts = 4
 
!                          NPTS
!         FORCE VECTOR = SUMMATION (S )(U )
!                          I=1       I   I
 
 DO  i = 1,npts
   
!     POINTER TO DISPLACEMENT VECTOR IN VARIABLE CORE
   
   npoint = ivec + nsil(i) - 1
   
   CALL smmats (si(30*i-29),5,6,0,  z(npoint),6,1,0,  vec,cvec )
   
   DO  j = 1,5
     cfrvec(j+1) = cfrvec(j+1) + cvec(j)
     forvec(j + 1) = forvec(j + 1) + vec(j)
   END DO
 END DO
 IF( tloads == 0 ) GO TO 55
 flag = .false.
 jst = 98
 IF( ntype == 4 ) jst = 128
 f1 = ti(6)
 IF( n1 == 1 ) GO TO 50
 forvec(2) = forvec(2) - ti(2)
 forvec(3) = forvec(3) - ti(3)
 forvec(4) = forvec(4) - ti(4)
 IF( ti(5) == 0.0 .AND. ti(6) == 0.0 ) flag = .true.
 GO TO 55
 50 forvec(2) = forvec(2) + ti(2)*est(jst+1)
 forvec(3) = forvec(3) + ti(2)*est(jst+2)
 forvec(4) = forvec(4) + ti(2)*est(jst+3)
 IF( ti(3) == 0.0 .AND. ti(4) == 0.0 ) flag = .true.
 55 CONTINUE
 
!     FORCE VECTOR COMPLETE AND CONTAINS M , M , M  , V , V
!                                         X   Y   XY   X   Y
 
!     AND ALSO INCLUDES THE ELEMENT ID AS THE FIRST ENTRY
!     ******************************************************************
 
!     STRESSES AT FIBER DISTANCES Z1 AND Z2 = - M * Z / I
 
 stres(2) = est(7)
 stres(11) = est(8)
 eye = est(6)
 i = 2
 k = 7
 k1 = 0
 200 zoveri = - stres(i) / eye
 zovi = ABS (zoveri)
 IF( tloads == 0  .OR.  flag ) GO TO 207
 j = 98
 IF( ntype == 4 ) j = 128
 IF( n1 == 1 ) GO TO 205
 
 ff = ti(k1+5) - ti(1)
 sdelta(1) = (est(jst+1)*ff + ti(2)*stres(i))/eye
 sdelta(2) = (est(jst+2)*ff + ti(3)*stres(i))/eye
 sdelta(3) = (est(jst+3)*ff + ti(4)*stres(i))/eye
 GO TO 210
 
 205 ff = (ti(k1+3) - stres(i)*ti(2) - ti(1)) / eye
 sdelta(1) = est(jst+1)*ff
 sdelta(2) = est(jst+2)*ff
 sdelta(3) = est(jst+3)*ff
 GO TO 210
 
 207 sdelta(1) = 0.0
 sdelta(2) = 0.0
 sdelta(3) = 0.0
 210 CONTINUE
 stres(i+1) = forvec(2) * zoveri - sdelta(1)
 stres(i+2) = forvec(3) * zoveri - sdelta(2)
 stres(i+3) = forvec(4) * zoveri - sdelta(3)
 cfrvec(k  ) = cfrvec(2) * zovi
 cfrvec(k+1) = cfrvec(3) * zovi
 cfrvec(k+2) = cfrvec(4) * zovi
 
!     PRINCIPAL STRESSES AND ANGLE OF ACTION PHI
 temp = stres(i+1) - stres(i+2)
 stres(i+7) = SQRT(  (temp/2.0E0)**2  +  stres(i+3)**2  )
 delta = (stres(i + 1) + stres(i + 2) ) / 2.0E0
 stres(i+5) = delta + stres(i+7)
 stres(i+6) = delta - stres(i+7)
 delta = 2.0E0 * stres(i+3)
 IF( ABS(delta) < 1.0E-15 .AND. ABS(temp) < 1.0E-15)GO TO 101
 stres(i+4) = ATAN2( delta,temp ) * 28.6478898E0
 GO TO 100
 101 stres(i+4) = 0.0E0
 100 IF( i == 11 ) GO TO 111
 i = 11
 k1 = 1
 k = 10
 GO TO 200
 111 stres( 1) = est(1)
 
!     ABOVE COMPLETES 2 VECTORS EACH WITH...
 
!     ELEM ID, Z, SIGMA X, SIGMA Y, SIGMA XY, PHI, SIG 1, SIG 2, TAU-MAX
 
!     STRESSES AND FORCES COMPLETE
 
 
!     ADDITON TO ELIMINATE 2ND ELEMENT ID IN OUTPUT
 
 DO  i = 10,17
   stres(i) = stres(i+1)
 END DO
 
!  . STRESS CHECK...
 
 IF (nchk <= 0) GO TO 350
 cfrvec(1) = est(1)
 k = 0
 
!  . FORCES...
 CALL sdrchk (forvec(2),cfrvec(2),5,k)
 
!  . STRESSES...
 CALL sdrchk (stres(3),cfrvec(7),3,k)
 CALL sdrchk (stres(11),cfrvec(10),3,k)
 IF (k == 0) GO TO 350
 
!  . LIMITS EXCEEDED...
 j = 0
 istyp(1) = tr
 IF (ieltyp == 15) istyp(1) = qd
 istyp(2) = plt
 IF (ieltyp == 7) istyp(2) = bsc
 
 IF (lsub == isub .AND. frlast(1) == frtmei(1) .AND.  &
     lld == ild  .AND. frlast(2) == frtmei(2) ) GO TO 320
 
 lsub = isub
 lld = ild
 frlast(1) = frtmei(1)
 frlast(2) = frtmei(2)
 j = 2
 CALL page1
 300 CALL sd2rhd (ished,j)
 line = line + 1
 WRITE(nout,310)
 310 FORMAT (7X,4HTYPE,5X,3HEID,5X,2HMX,5X,2HMY,4X,3HMXY,5X,2HVX,5X,  &
     2HVY,4X,3HSX1,4X,3HSY1,3X,4HSXY1,4X,3HSX2,4X,3HSY2,3X,4HSXY2)
 GO TO 330
 
 320 IF (eject(2) /= 0) GO TO 300
 330 WRITE(nout,340) istyp,nelid,(cfrvec(i),i=2,12)
 340 FORMAT (1H0,4X,a4,a3,i7,11F7.1)
 
 350 CONTINUE
 RETURN
END SUBROUTINE sbspl2
