SUBROUTINE trd1e(mhh,bhh,khh,ph,uhv,ngroup)
     
!     THIS ROUTINE SOLVES TRANSIENT PROBLEM ANALYTICALLY IN CASE
!         OF UNCOUPLED MODAL WITH NO NONLINEAR LOADS
 
 
 INTEGER, INTENT(IN)                      :: mhh
 INTEGER, INTENT(IN)                      :: bhh
 INTEGER, INTENT(IN)                      :: khh
 INTEGER, INTENT(IN OUT)                  :: ph
 INTEGER, INTENT(IN OUT)                  :: uhv
 INTEGER, INTENT(IN)                      :: ngroup
 REAL :: mi,ki
 INTEGER :: iz(1),sysbuf,iuhv(7), FILE
 INTEGER :: NAME(2)
 
!RLBNB SPR94003 9/94
 COMMON /BLANK / dummy(4), ncol
!RLBNE
 COMMON /packx/ it1,it2,ii,jj,incur
 COMMON /unpakx/it3,iii,jjj,incur1
 COMMON /zzzzzz/ z(1)
 COMMON /system/ sysbuf
 
 EQUIVALENCE (iz(1),z(1))
 
 DATA NAME/4HTRD1,4HE   /
 DATA epsi/1.0E-8/
!*********
!     DEFINITION OF VARIABLES
!*********
!     IGROUP   POINTER TO TIME STEP DATA  N1,DELTAT,NO
!     NGROUP   NUMBER OF TIME STEP CHANGES
!     MHH      MODAL MASS FILE
!     KHH      MODAL STIFFNESS FILE
!     BHH      MODAL DAMPING FILE
!     PH       LOAD FILE
!     UHV      DISPLACEMENT,VELOCITY, AND ACCELERATION FILE
!     NMODES   ORDER OF MODAL FORMULATION
!     IMII     POINTER TO MASSES
!     IBII     POINTER TO DAMPING
!     IKII     POINTER TO STIFFNESS
!     IF       POINTER TO F-S
!     IFPR     POINTER TO F PRIMES
!     IG       POINTER TO G-S
!     IGPR
!     IA       POINTER TO A-S
!     IAPR
!     IB       POINTER TO B-S
!     IBPR
!     IUJ      POINTER TO OLD  DISP
!     IUJ1             TO NEW  DISP
!     IUDJ     POINTER TO  OLD VELOCITY VECTOR
!     IUDJ1                NEW VELOCITY VECTOR
!     IPHJ     POINTER TO  OLD LOAD VECTOR
!     IPHJ1                NEW LOAD VECTOR
!     NSTEP    NUMBER OF STEPS AT CURRENT INCREMENT
!     H        CURRENT DELTA T
!     NOUT     OUTPUT INCURMENT
!     EPSI     CASE SELTION TOLERANCE
 
!********    HERE WE GO --GET LOTS OF PAPER
 
 lc = korsz(z)
 lc =lc -ngroup*3
 igroup = lc+1
 ist =-1
 ibuf1 =lc -sysbuf
 ibuf2 =ibuf1 -sysbuf
 lc = lc - 2*sysbuf
 iuhv(1)= mhh
 CALL rdtrl(iuhv)
 nmodes = iuhv(2)
 it1=1
 it2=1
 it3=1
 incur=1
 incur1=1
 ii=1
 jj=nmodes
 icrq = 17*nmodes - lc
 IF(icrq > 0) GO TO 340
 
!     BRING IN H MATRICES
 
 
!     BRING IN  MHH
 FILE =mhh
 imii =0
 kk=imii
 ASSIGN 10 TO iretn
 GO TO 280
 
!     BRING IN BHH
 10 DO  j=1,nmodes
   IF(z(j) == 0.0) GO TO 350
 END DO
 FILE = bhh
 ibii= imii+ nmodes
 kk = ibii
 ASSIGN 20 TO iretn
 GO TO 280
 
!     BRING IN KHH
 20 FILE =khh
 ikii = ibii +nmodes
 kk= ikii
 ASSIGN 30 TO iretn
 GO TO 280
 
!     ASSIGN ADDITIONAL POINTERS
 
 30 iii=1
 jjj=nmodes
 IF = ikii + nmodes
 ig = IF   + nmodes
 ia = ig   + nmodes
 ib = ia   + nmodes
 ifpr=ib   + nmodes
 igpr=ifpr + nmodes
 iapr=igpr + nmodes
 ibpr=iapr + nmodes
 iuj =ibpr + nmodes
 iuj1=iuj  + nmodes
 iudj=iuj1 + nmodes
 iudj1=iudj+ nmodes
 iphj =iudj1+nmodes
 iphj1=iphj +nmodes
!RLBNB SPR94003 9/94
 IF (ncol <= 2) GO TO 37
 
!     RETRIEVE OLD DISPLACEMENT AND VELOCITY
!     FROM A PREVIOUSLY CHECKPOINTED RUN
 
 CALL gopen (uhv, iz(ibuf1), 0)
 i = 3*(ncol - 1)
 CALL skprec (uhv, i)
 
!     RETRIEVE OLD DISPLACEMENT
 
 CALL unpack (*31, uhv, z(iuj1+1))
 GO TO 33
 31 DO  i = 1, nmodes
   k = iuj1 + i
   z(k) = 0.0
 END DO
 
!     RETRIEVE OLD VELOCITY
 
 33 CALL unpack (*34, uhv, z(iudj1+1))
 GO TO 36
 34 DO  i = 1, nmodes
   k = iudj1 + i
   z(k) = 0.0
 END DO
 36 CALL CLOSE (uhv, 1)
!RLBNE
 
!     READY UHV
 
!RLBR SPR94003 9/94      CALL GOPEN(UHV,IZ(IBUF1),1)
 37 CALL gopen(uhv,iz(ibuf1),1)
 CALL makmcb(iuhv,uhv,nmodes,2,1)
 
!     READY LOADS
 
 CALL gopen(ph,iz(ibuf2),0)
 CALL unpack(*40,ph,z(iphj1+1))
 GO TO 60
 
!     ZERO LOAD
 
 40 DO  i=1,nmodes
   k = iphj1+i
   z(k) = 0.0
 END DO
!RLBNB SPR94003 9/94
 60 IF (ncol > 2) GO TO 75
!RLBNE
 
!     ZERO INITIAL DISPLACEMENT AND VELOCITY
 
!RLBR SPR 94003 9/94   60 DO 70 I=1,NMODES
 DO  i=1,nmodes
   k = iuj1+i
   z(k) = 0.0
   k = iudj1+i
   z(k) = 0.0
 END DO
 
!     BEGIN LOOP ON EACH DIFFERENT TIME STEP
 
!RLBR SPR 94003 9/94      I = 1
 75 i = 1
 80 nstep = iz(igroup)
 IF(i == 1) nstep = nstep+1
 h     =  z(igroup+1)
 nout = iz(igroup+2)
 igroup = igroup +3
 jk = 1
 IF(i == 1) GO TO 170
 
!     COMPUTE F-S ,G-S,A-S,B-S
 
 90 DO  j=1,nmodes
   k= imii+j
   mi= z(k)
   IF(mi == 0.0) GO TO 350
   k= ibii+j
   bi= z(k)
   k = ikii+j
   ki= z(k)
   wosq =ki/mi
   beta = bi/(2.0*mi)
   betasq =beta*beta
   wsq  = ABS(wosq - betasq)
   w = SQRT(wsq)
   IF(SQRT(wsq + betasq)*h < 1.e-6) GO TO 100
   t1 = ( wosq-betasq ) / wosq
   IF( t1 > epsi ) GO TO 110
   IF( t1 < -epsi) GO TO 130
   
!     CASE  3  CRITICALLY DAMPED
   
   bh = beta*h
   expbh = EXP(-bh)
   t1 = h*ki
   k = IF+j
   
!     COMPUTE F
   
   z(k) = expbh*(1.0 +bh)
   
!     COMPUTE  G
   
   k = ig +j
   z(k)= h*expbh
   
!     COMPUTE A
   
   k = ia +j
   z(k) = (2.0/beta - expbh/beta*(2.0 +2.0*bh + bh*bh))/ t1
   
!     COMPUTE B
   
   k=ib +j
   z(k) = (-2.0 +bh+expbh*(2.0+bh))/(bh*ki)
   
!     COMPUTE  F PRIME
   
   k= ifpr+j
   z(k) = -betasq*h*expbh
   
!     COMPUTE  G PRIME
   
   k = igpr+j
   z(k)= expbh*(1.0- bh)
   
!     COMPUTE A PRIME
   
   k = iapr +j
   z(k) = (expbh*(1.0 + bh + bh*bh)- 1.0)/t1
   
!     COMPUTE  B PRIME
   
   k = ibpr +j
   z(k) = (1.0 -expbh*(bh +1.0))/t1
   CYCLE
   
!     CASE  4   W0 = BETA =0.0
   
   100 k=IF+j
   z(k)=1.0
   k= ig+j
   z(k)=h
   k= ia+j
   z(k)= h*h/(3.0*mi)
   k= ib+j
   z(k)= h*h/(6.0*mi)
   k= ifpr+j
   z(k)=0.0
   k= igpr+j
   z(k)=1.0
   t1 = h/(2.0*mi)
   k = iapr+j
   z(k)= t1
   k=  ibpr+j
   z(k)= t1
   CYCLE
   
!     CASE 1 --UNDERDAMPED
   
   110 wh = w*h
   expbh = EXP(-beta*h)
   sinwh = SIN(wh)
   coswh = COS(wh)
   
!     COMPUTE F
   
   120 k= IF +j
   z(k)= expbh*(coswh +beta/w *sinwh)
   
!     COMPUTE G
   
   k = ig +j
   z(k) = expbh/w*sinwh
   
!     COMPUTE A
   
   k= ia+j
   t1 =(wsq -betasq)/wosq
   t2 = 2.0*w*beta/wosq
   t3 = wh*ki
   z(k)= (expbh*((t1-beta*h)*sinwh-(t2+wh)*coswh)+t2)/t3
   
!     COMPUTE  B
   
   k =ib +j
   z(k) = (expbh*(-t1*sinwh + t2*coswh)+wh- t2)/t3
   
!     COMPUTE  FPRIME
   
   k = ifpr+j
   z(k) = -wosq/w*expbh*sinwh
   
!     COMPUTE G PRIME
   
   k =igpr +j
   z(k) = expbh*(coswh -beta/w *sinwh)
   
!     COMPUTE A PRIME
   
   k = iapr +j
   z(k) =(expbh*((beta +wosq*h)*sinwh +w*coswh)- w)/t3
   
!     COMPUTE B PRIME
   
   k =ibpr +j
   z(k) = (-expbh*(beta*sinwh +w*coswh) + w)/t3
   CYCLE
   
!     CASE  3    W0 - BETASQ L -E
   
   130 wh =w*h
   expbh= EXP(-beta*h)
   sinwh =   SINH(wh)
   coswh =   COSH(wh)
   betasq = -betasq
   GO TO 120
 END DO
 
!     BEGIN LOOP ON INCREMENTS
 
 
!     COMPUTE  NEW DISPLACEMENTS
 
 150 k = iuj1
 kk=iudj1
 DO  l=1,nmodes
   k=k+1
   kk =kk+1
   z(k)=0.0
   z(kk)=0.0
   kkk = IF+l
   kd =  iuj +l
   z(k) =z(kkk)*z(kd) +z(k)
   kkk = ifpr +l
   z(kk) = z(kkk)*z(kd) +z(kk)
   kd= iudj+l
   kkk = ig +l
   z(k) = z(kkk)*z(kd) +z(k)
   kkk = igpr +l
   z(kk) = z(kkk)*z(kd) +z(kk)
   kd = iphj +l
   kkk = ia +l
   z(k) = z(kkk)*z(kd) +z(k)
   kkk  = iapr +l
   z(kk)= z(kkk)*z(kd) + z(kk)
   kd = iphj1+l
   kkk=  ib +l
   z(k) = z(kkk)*z(kd) +z(k)
   kkk  = ibpr +l
   z(kk) = z(kkk)*z(kd) + z(kk)
 END DO
 IF(jk == nstep) GO TO 200
 IF( jk /= 1 .AND. MOD(jk+ist,nout) /= 0) GO TO 180
 
!     TIME TO OUTPUT--YOU LUCKY FELLOW
 
 170 ASSIGN 190 TO iretn
 GO TO 220
 180 ASSIGN 190 TO iretn
 GO TO 240
 190 jk = jk+1
 IF(jk == 2 .AND. i == 1) GO TO 90
 IF(jk <= nstep) GO TO 150
 200 ASSIGN 210 TO iretn
 GO TO 220
 210 i =i+1
 ist = 0
 IF( i <= ngroup) GO TO 80
 CALL CLOSE(ph,1)
 CALL CLOSE(uhv,1)
 CALL wrttrl(iuhv)
 RETURN
 
!     INTERNAL SUBROUTINE FOR OUTPUT AND VELOCITY COMPUTE
 
 220 CALL pack(z(iuj1+1),uhv,iuhv)
 CALL pack(z(iudj1+1),uhv,iuhv)
 
!     COMPUTE  ACCELERATIONS
 
 DO  l=1,nmodes
   k= iudj+l
   kk=iphj1+l
   kkk = imii+l
   kd = ibii+l
   kd1= iudj1+l
   kd2= iuj1 +l
   kd3 = ikii+l
   z(k) = z(kk)/z(kkk)-z(kd)*z(kd1)/z(kkk)-z(kd3)*z(kd2)/z(kkk)
 END DO
 CALL pack(z(iudj+1),uhv,iuhv)
 
!     SWITCH POINTS TO STUFF
 
 240 kd= iuj
 iuj = iuj1
 iuj1=kd
 kd= iudj
 iudj =iudj1
 iudj1=kd
 kd = iphj
 iphj =iphj1
 iphj1= kd
 
!     BRING IN NEXT LOAD VECTOR
 
 CALL unpack(*260,ph,z(iphj1+1))
 250 GO TO iretn,(190,210)
 260 DO  kd=1,nmodes
   k = iphj1 +kd
   z(k) =0.0
 END DO
 GO TO 250
 
!     INTERNAL SUBROUTINE TO BRING  IN H MATRICES
 
 280 CALL OPEN(*302,FILE,iz(ibuf1),0)
 CALL skprec(FILE,1)
 DO  kd=1,nmodes
   iii= kd
   jjj= kd
   kd1= kk+kd
   CALL unpack(*290,FILE,z(kd1))
   CYCLE
   290 z(kd1)= 0.0
 END DO
 CALL CLOSE(FILE,1)
 301 GO TO iretn,(10,20,30)
 
!      ZERO CORE FOR PURGED FILES
 
 302 DO  kd = 1,nmodes
   kd1 = kk + kd
   z(kd1) = 0.0
 END DO
 GO TO 301
 
!     ERROR MESAGES
 
 320 CALL mesage(ip1,FILE,NAME)
 RETURN
 340 ip1 = -8
 FILE = icrq
 GO TO 320
 350 ip1 = -43
 FILE = j
 GO TO 320
END SUBROUTINE trd1e
