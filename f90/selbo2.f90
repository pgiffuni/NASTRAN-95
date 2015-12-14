SUBROUTINE selbo2 (ti)
     
!     THIS ROUTINE IS THE PHASE II SUBROUTINE OF STRESS DATA RECOVERY
!     FOR THE BEAM ELEMENT.
 
 
 REAL, INTENT(IN)                         :: ti(14)
 INTEGER :: tloads
 REAL :: i1,i2,l,m1a,m2a,m1b,m2b,i12,k1a,k2a,k1b,k2b,  m2bt
 EQUIVALENCE     (ldtemp,templd),(msten,smten),(mscom,smcom)
 COMMON /zzzzzz/ zz(1)
 COMMON /sdr2x4/ xxxxxx(33),icstm,ncstm,ivec,ivecn,ldtemp,eldefm,  &
     dum8(8),tloads
 
!     THE FIRST 100 LOCATIONS OF THE SDR2X7 BLOCK ARE RESERVED FOR INPUT
!     PARAMETERS, THE SECOND 100 FOR STRESS OUTPUT PARAMETERS, AND FORCE
!     OUTPUT PARAMETERS BEGIN AT LOCATION 201.
 
 COMMON /sdr2x7/ jelid,jsilno(2),sa(36),sb(36),st,sdelta,a,fj,i1,  &
     i2,c,r1,t1,r2,t2,r3,t3,r4,t4,t_sub_0,sigmat, sigmac,l,r,betar,therm(4)
 COMMON /sdr2x7/ iselid,sig1a,sig2a,sig3a,sig4a,sigax,sigamx,  &
     sigamn,msten,sig1b,sig2b,sig3b,sig4b,sigbx, sigbmx,sigbmn,mscom,yyyyyy(83)
 COMMON /sdr2x7/ ifelid,m1a,m2a,v1,v2,fx,t,m1b,m2bt,v1bt,fxbt,tbt
 
 COMMON /sdr2x8/ fa(6),fb(6),idisp,iua,iub,p1,k1a,k2a,k1b,k2b,q,w
 DATA    dcr   / .017453292 /
 
 sid(x) = SIN(x*dcr)
 cod(x) = COS(x*dcr)
 
 x   = 1.0
 yl  = r*(1.-cod(betar))
 xl  = r*sid(betar)
 i12 = 0.
 idisp = ivec - 1
 iua = idisp + jsilno(1)
 CALL gmmats (sa(1),6,6,0, zz(iua),6,1,0, fa(1))
 iub = idisp + jsilno(2)
 CALL gmmats (sb(1),6,6,0, zz(iub),6,1,0, fb(1))
 fx  = -fa(1) - fb(1)
 v1  = -fa(2) - fb(2)
 v2  = -fa(3) - fb(3)
 t   = -fa(4) - fb(4)
 m2a =  fa(5) + fb(5)
 m1a = -fa(6) - fb(6)
 
!     IF LDTEMP = -1, THE LOADING TEMPERATURE IS UNDEFINED
 
 IF (tloads == 0) GO TO 10
 tbar = ti(1)
 dt   = tbar - tsub0
 DO  i = 1,6
   fa(i) = dt*therm(i)
 END DO
 fx  = fx  + fa(1)
 v1  = v1  + fa(2)
 m1a = m1a + fa(6)
 10    m1b = m1a - v1*xl + fx*yl
 m2b = m2a - v2*xl
 tb  = t   - v2*yl
 
!     TRANSFORM FORCES AT B-END TO A COORD. SYS TANGENT TO B-END
 
 fxbt = v1*sid(betar) + fx*ABS(cod(betar))
 v1bt = v1*ABS(cod(betar))  - fx*sid(betar)
 m2bt = m2b*ABS(cod(betar)) + tb*sid(betar)
 tbt  =-m2b*sid(betar) + tb*ABS(cod(betar))
 
!     COMPUTE ELEMENT STRESSES AT 4 POINTS
 
 
!     COMPUTE K1A AND K2A
 
 IF (i12 /= 0.0) GO TO 30
 IF (i1  /= 0.0) GO TO 20
 k1a = 0.0
 GO TO 40
 20 k1a = -m1a/i1
 GO TO 40
 30 k1a = (m2a*i12 - m1a*i2)/(i1*i2 - i12**2)
 k2a = (m1a*i12 - m2a*i1)/(i1*i2 - i12**2)
 GO TO 60
 40 IF (i2 /= 0.0) GO TO 50
 k2a = 0.0
 GO TO 60
 50 k2a = -m2a/i2
 
!     CHANGE STRESS RECOVERY CONSTANTS FROM CYL. TO RECT. COORD.
 
 c1 = r1*sid(t1)
 c2 = r1*cod(t1)
 d1 = r2*sid(t2)
 d2 = r2*cod(t2)
 f1 = r3*sid(t3)
 f2 = r3*cod(t3)
 g1 = r4*sid(t4)
 g2 = r4*cod(t4)
 
!     COMPUTE SIG1A, SIG2A, SIG3A AND SIG4A
 
 60 sig1a = k1a*c1*c + k2a*c2
 sig2a = k1a*d1*c + k2a*d2
 sig3a = k1a*f1*c + k2a*f2
 sig4a = k1a*g1*c + k2a*g2
 
!     COMPUTE K1B AND K2B
 
 IF (i12 /= 0.0) GO TO 80
 IF (i1  /= 0.0) GO TO 70
 k1b = 0.0
 GO TO 90
 70 k1b = -m1b/i1
 GO TO 90
 80 k1b = (m2bt*i12 - m1b *i2)/(i1*i2 - i12**2)
 k2b = (m1b *i12 - m2bt*i1)/(i1*i2 - i12**2)
 GO TO 110
 90 IF (i2 /= 0.0) GO TO 100
 k2b = 0.0
 GO TO 110
 100 k2b = -m2bt/i2
 
!     COMPUTE SIG1B, SIG2B, SIG3B AND SIG4B
 
 110 sig1b = k1b*c1*c + k2b*c2
 sig2b = k1b*d1*c + k2b*d2
 sig3b = k1b*f1*c + k2b*f2
 sig4b = k1b*g1*c + k2b*g2
 IF (tloads == 0) GO TO 115
 
!     TEST IF AT LEAST ONE POINT TEMPERATURE IS GIVEN
 
 DO  i = 7,14
   IF (ti(i) /= 0.0) GO TO 112
 END DO
 GO TO 115
 112 IF (a == 0.0) GO TO 115
 ealf  =-st/a
 sig1a = sig1a + ealf*(ti( 7) - ti(3)*c1*c - ti(5)*c2 - ti(1))
 sig2a = sig2a + ealf*(ti( 8) - ti(3)*d1*c - ti(5)*d2 - ti(1))
 sig3a = sig3a + ealf*(ti( 9) - ti(3)*f1*c - ti(5)*f2 - ti(1))
 sig4a = sig4a + ealf*(ti(10) - ti(3)*g1*c - ti(5)*g2 - ti(1))
 sig1b = sig1b + ealf*(ti(11) - ti(4)*c1*c - ti(6)*c2 - ti(2))
 sig2b = sig2b + ealf*(ti(12) - ti(4)*d1*c - ti(6)*d2 - ti(2))
 sig3b = sig3b + ealf*(ti(13) - ti(4)*f1*c - ti(6)*f2 - ti(2))
 sig4b = sig4b + ealf*(ti(14) - ti(4)*g1*c - ti(6)*g2 - ti(2))
 115 CONTINUE
 
!     COMPUTE AXIAL STRESS
 
 sigax = 0.0
 sigbx = 0.0
 IF (a /= 0.0) sigax = fx/a
 IF (a /= 0.0) sigbx = fxbt/a
 
!     COMPUTE MAXIMA AND MINIMA
 
 sigamx = sigax + AMAX1(sig1a,sig2a,sig3a,sig4a)
 sigbmx = sigbx + AMAX1(sig1b,sig2b,sig3b,sig4b)
 sigamn = sigax + AMIN1(sig1a,sig2a,sig3a,sig4a)
 sigbmn = sigbx + AMIN1(sig1b,sig2b,sig3b,sig4b)
 
!     COMPUTE MARGIN OF SAFETY IN TENSION
 
 IF (sigmat <= 0.0) GO TO 620
 IF (AMAX1(sigamx,sigbmx) <= 0.0) GO TO 620
 q = sigmat/AMAX1(sigamx,sigbmx)
 smten = q - 1.0
 GO TO 630
 620 msten = 1
 
!     COMPUTE MARGIN OF SAFETY IN COMPRESSION
 
 630 IF (sigmac <= 0.0) GO TO 640
 IF (AMIN1(sigamn,sigbmn) >= 0.0) GO TO 640
 w = -sigmac/AMIN1(sigamn,sigbmn)
 smcom  = w - 1.0
 GO TO 150
 640 mscom  = 1
 150 iselid = jelid
 ifelid = jelid
 RETURN
END SUBROUTINE selbo2
