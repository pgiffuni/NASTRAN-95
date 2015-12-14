SUBROUTINE dtrmem( iopt )
     
!     DIFFERENTIAL STIFFNESS CALCULATIONS FOR THE TRIANGULAR MEMBRANE
!     ELEMENT.  THREE 6X6 MATRICES FOR THE PIVOT POINT ARE INSERTED.
!     IF THIS ROUTINE IS CALLED FROM DTRIA OR DQUAD ONLY THE IN PLANE
!     EFFECTS ARE GENERATED AND THE STRESS VALUES ARE RETURNED.
 
!     THE VALUE OF IOPT TELLS US WHICH ROUTINE IS CALLING DTRMEM.
!      THE OPTIONS ARE
!                IOPT        ROUTINE
!               ******       *******
!                 0            DSIA
!                 1            DQDMEM
!                 2            DTRIA
!                 3            DQUAD
 
 
!     THIS ROUTINE COMPUTES AN E-MATRIX UNIQUE TO THIS ROUTINE.
 
!                       IX  IY  IZ
!                  E =  JX  JY  JZ
!                       KX  KY  KZ
 
 
 INTEGER, INTENT(IN)                      :: iopt
 DOUBLE PRECISION :: e             ,c ,kd            ,sigx  &
     ,sigy          ,sigxy ,temp1         ,temp2  &
     ,kij ,g             ,xsubb  &
     ,xsubc         ,ysubc ,sum           ,mu  &
     ,lamda         ,delta ,temp          ,gamma1  &
     ,gamma2        ,gamma3 ,disp          ,t
 DOUBLE PRECISION :: areat         ,dumdp
 
 DIMENSION          sum(3)        ,necpt(6)      ,kij(36)
 
 
!     INTERFACE DATA BLOCKS
 
 COMMON /condas/ consts(5)
 COMMON /matin / matid, inflag, eltemp, stress, sinth, costh
 COMMON /matout/ g11,g12,g13,g22,g23,g33,rho,alpha1,alpha2,alph12,  &
     tsub0,gsube,sigten,sigcom,sigshe,g2x211,g2x212,g2x222
 COMMON /ds1aaa/ npvt, icstm, ncstm
 COMMON /ds1aet/ ecpt(21),eldef,ldtemp,sdisp(9)
 COMMON /ds1adp/    e(9)          ,c(54) ,kd(36)        ,temp1(18)  &
     ,temp2(18) ,g(9)          ,t(9)  &
     ,disp(9)       ,mu ,lamda         ,delta  &
     ,temp          ,gamma1 ,gamma2        ,gamma3  &
     ,areat         ,xsubb ,xsubc         ,ysubc  &
     ,dumdp(12)     ,theta ,icstm1        ,npivot  &
     ,idum ,sigx          ,sigy           ,sigxy
 
 EQUIVALENCE ( consts(4) , degra  )
 EQUIVALENCE (ldtemp,ftemp),(necpt(1),ecpt(1)),(sum(1),sigx)
 EQUIVALENCE        (kij(1),kd(1))
 
 
!     ******************************************************************
!     ECPT( 1) = ELEMENT ID
!     ECPT( 2) = GRID POINT A OR 1
!     ECPT( 3) = GRID POINT B OR 2
!     ECPT( 4) = GRID POINT C OR 3
!     ECPT( 5) = THETA = ANGLE OF MATERIAL CUT IF ANISOTROPIC
!     ECPT( 6) = MATERIAL ID
!     ECPT( 7) = THICKNESS
!     ECPT( 8) = NON-STRUCTURAL MASS
!     ECPT( 9) = COORD. SYSTEM ID 1
!     ECPT(10) = X1
!     ECPT(11) = Y1
!     ECPT(12) = Z1
!     ECPT(13) = COORD. SYSTEM ID 2
!     ECPT(14) = X2
!     ECPT(15) = Y2
!     ECPT(16) = Z2
!     ECPT(17) = COORD. SYSTEM ID 3
!     ECPT(18) = X3
!     ECPT(19) = Y3
!     ECPT(20) = Z3
!     ECPT(21) = ELEMENT TEMPERATURE
!     ECPT(22) = ELEMENT DEFORMATION DELTA
!     ECPT(23) = AVG. LOADING TEMPERATURE =(-1) IF NO LOADING TEMP.
!     ECPT(24) = X-TRANS POINT 1
!     ECPT(25) = Y-TRANS POINT 1
!     ECPT(26) = Z-TRANS POINT 1
!     ECPT(27) = X-TRANS POINT 2
!     ECPT(28) = Y-TRANS POINT 2
!     ECPT(29) = Z-TRANS POINT 2
!     ECPT(30) = X-TRANS POINT 3
!     ECPT(31) = Y-TRANS POINT 3
!     ECPT(32) = Z-TRANS POINT 3
!     ******************************************************************
!//////
!     CALL BUG(4HTMET,0,ECPT,32)
!//////
 
 sigx=0.0D0
 sigy=0.0D0
 sigxy=0.0D0
 IF(ecpt(7) == 0.0 .OR. necpt(6) == 0 ) RETURN
!     FILL ELEMENT TO GLOBAL E-TRANSFORMATION MATRIX
 
!     IVEC = E(1). . .E(3)
!     JVEC = E(4). . .E(6)
!     KVEC = E(7). . .E(9)
 
 DO  i=1,3
   e(i) = DBLE( ecpt(i+13) ) - DBLE( ecpt(i+9) )
 END DO
 
!     LENGTH THEN = XSUBB
 
 xsubb = DSQRT( e(1)**2 + e(2)**2 + e(3)**2 )
 
!     R  - R   (INTERMEDIATE STEP) AND NOMALIZE IVECTOR = E(1). . .E(3)
!      C    A
 
 DO  i=1,3
   e(i+3) = DBLE( ecpt(i+17) ) - DBLE( ecpt(i+9) )
   e(i) = e(i) / xsubb
 END DO
 
!     XSUBC = I DOT (R  - R )
!                     C    A
 
 xsubc = e(1) * e(4)  +  e(2) * e(5)  +  e(3) * e(6)
 
!     KVEC = IVEC CROSS (R  - R )
!                         C    A
 
 e(7) = e(2) * e(6)  -  e(3) * e(5)
 e(8) = e(3) * e(4)  -  e(1) * e(6)
 e(9) = e(1) * e(5)  -  e(2) * e(4)
 
!     LENGTH = YSUBC
 
 ysubc = DSQRT(e(7)**2 + e(8)**2 + e(9)**2 )
 
!     NORMALIZE KVECTOR
 e(7) = e(7) / ysubc
 e(8) = e(8) / ysubc
 e(9) = e(9) / ysubc
 
!     JVECTOR = I CROSS K
 
 e(4) = e(3) * e(8)  -  e(2) * e(9)
 e(5) = e(1) * e(9)  -  e(3) * e(7)
 e(6) = e(2) * e(7)  -  e(1) * e(8)
 
!     NORMALIZE JVECTOR TO MAKE SURE
 temp = DSQRT( e(4)**2 + e(5)**2 + e(6)**2 )
 e(4) = e(4) / temp
 e(5) = e(5) / temp
 e(6) = e(6) / temp
 
!     MU, LAMDA, AND DELTA
 
 mu    = 1.0D0 / xsubb
 lamda = 1.0D0 / ysubc
 delta =(xsubc/xsubb) - 1.0D0
 areat = xsubb * ysubc * 0.50D0 * DBLE( ecpt(7) )
 
!     C MATRIX    C  =(3X2) STORED C( 1). . .C( 6)
!                  A
!                 C  =(3X2) STORED C( 7). . .C(12)
!                  B
!                 C  =(3X2) STORED C(13). . .C(18)
!                  C
 
 c( 1) = -mu
 c( 2) =  0.0D0
 c( 3) =  0.0D0
 c( 4) =  lamda * delta
 c( 5) =  c(4)
 c( 6) = -mu
 c( 7) =  mu
 c( 8) =  0.0D0
 c( 9) =  0.0D0
 c(10) = -lamda * mu * xsubc
 c(11) =  c(10)
 c(12) =  mu
 c(13) =  0.0D0
 c(14) =  0.0D0
 c(15) =  0.0D0
 c(16) =  lamda
 c(17) =  lamda
 c(18) =  0.0D0
 
 IF( iopt >= 1 ) GO TO 30
!     THE REASON FOR THIS IS THAT IF THE DQDMEM ROUTINE IS CALLING,
!     EACH INDIVIDUAL SUBTRIANGLE WILL ALREADY HAVE A SINTH AND COSTH.
 
 theta = ecpt(5) * degra
 sinth = SIN( theta )
 costh = COS( theta )
 30 IF( ABS(sinth) < 1.0E-06 ) sinth = 0.0E0
 
 eltemp = ecpt(21)
 matid = necpt(6)
 inflag = 2
 CALL mat( ecpt(1) )
 
!     FILL G-MATRIX WITH OUTPUT FROM MAT ROUTINE.
 
 g(1) = g11
 g(2) = g12
 g(3) = g13
 g(4) = g12
 g(5) = g22
 g(6) = g23
 g(7) = g13
 g(8) = g23
 g(9) = g33
 
!     G, E, C MATRICES ARE COMPLETE
 
!     FOLLOWING COMPUTES SIG , SIG , SIG      (3X1) VECTOR
!                           X     Y     XY
 
!         I=3
!      = (SUM (G)(C )(E)(T )(DISP )) - (S )(LDTEMP - T )
!         I=1      I      I      I       T            0
 
!        WHERE  S  =(G)(ALPHAS)   (3X1)
!                T
 
 sum(1) = 0.0E0
 sum(2) = 0.0E0
 sum(3) = 0.0E0
 
!     MAKE DISPLACEMENT VECTOR DOUBLE PRECISION
 
 DO  i=1,9
   disp(i) = sdisp(i)
 END DO
 
 DO  i=1,3
!     DO WE NEED TRANSFORMATIONS
   
   IF(necpt(4*i+5) == 0) THEN
     GO TO    60
   END IF
   50 CALL transd( necpt(4*i+5),t(1))
   CALL gmmatd( t(1),3,3,0, disp(3*i-2),3,1,0, temp1(1))
   GO TO 80
   
   60 DO  j=1,3
     idum=  3*(i-1)+j
     temp1(j) = disp(idum)
   END DO
   
   80 CALL gmmatd( e(1),2,3,0,temp1(1),3,1,0, temp2(1)  )
   CALL gmmatd( c(6*i-5),3,2,0,  temp2(1),2,1,0,  temp1(1) )
   CALL gmmatd( g(1),3,3,0,    temp1(1),3,1,0,    temp2(1) )
   
   sum(1) = sum(1) + temp2(1)
   sum(2) = sum(2) + temp2(2)
   sum(3) = sum(3) + temp2(3)
   
 END DO
 
 IF( ldtemp == (-1) ) GO TO 110
!     COMPUTE S MATRIX
!               T
 
 temp2(1) = alpha1
 temp2(2) = alpha2
 temp2(3) = alph12
!     ABOVE IS FOR SINGLE TO DOUBLE PRECISION.
 
 CALL gmmatd( g(1),3,3,0,  temp2(1),3,1,0,  temp1(1) )
 temp = ftemp - tsub0
 DO  i=1,3
   sum(i) = sum(i) - temp1(i) * temp
 END DO
 
!//////
!     CALL BUG(4HSUMS,90,SUM,6)
!//////
!  90 AT 90 SIG = SUM(1),  SIG = SUM(2),  SIG   = SUM(3)
!              X              Y              XY
 
!     ABOVE SIMULATES SMA,SDR2-PHASE I+II
!     FROM ABOVE THE E MATRIX (3X3), AND THE SUM (3X1) MATRIX ALONG WITH
!     XSUBB, XSUBC, AND YSUBC ARE NOW USED...
 110 DO   i =1,36
   kd(i) =0.0D0
 END DO
 
 IF( iopt == 3 ) areat=areat/2.0D0
 
 mu = sigx*areat
 lamda = sigy*areat
 delta = sigxy *areat
 
 IF ( iopt >= 2)  GO TO 130
 kd(1) = lamda
 kd(2) =-delta
 kd(7) = kd(2)
 kd(8) = mu
 130 kd(15) = mu+lamda
 kd(16) =-delta
 kd(17) = delta
 kd(18) = mu -lamda
 kd(21) = kd(16)
 kd(27) = kd(17)
 kd(33) = kd(18)
 
!     GENERATE C MATRICES
 
 DO  i=1,54
   c(i) =0.0D0
 END DO
 
!     FILL NON ZERO TERMS
 
 gamma1 = 1.0D0 /xsubb
 gamma2 = 1.0D0 /ysubc
 gamma3 = xsubc /( xsubb*ysubc)
 c(3) = gamma3 -gamma2
 c(6) = gamma1
 c(7) =-c(3)/2.0D0
 c(8) =-gamma1/2.0D0
 c(10)=-gamma1
 c(14)= c(3)
 c(16)=-c(7)
 c(17)= c(8)
 
 c(21)=-gamma3
 c(24)=-gamma1
 c(25)= gamma3/2.0D0
 c(26)=-c(8)
 c(28)= gamma1
 c(32)=-gamma3
 c(34)=-c(25)
 c(35)= c(26)
 
 c(39)=gamma2
 c(43)=-gamma2/2.0D0
 c(50)= gamma2
 c(52)=-c(43)
 
!     REPLACE C MATRICES BY  (C)(E )(T) FOR EACH POINT
 DO  i =1,3
   IF( necpt(4*i+5) == 0) THEN
     GO TO   160
   END IF
   
!     GLOBAL TO BASIC MATRIX T IS GENERATED AGAIN HERE
   
   150 CALL transd( necpt(4*i+5),t(1))
   CALL gmmatd( e(1),3,3,0,  t(1),3,3,0, temp1(1) )
   GO TO 180
   160 DO  j =1,9
     temp1(j) = e(j)
   END DO
   
   180 CALL gmmatd( c(18*i-17),6,3,0,  temp1(1),3,3,0,  temp2(1))
   DO   j=1,18
     idum =  18*(i-1) +j
     c(idum) = temp2(j)
   END DO
   
 END DO
 
 DO  i =1,3
   IF(necpt(i+1) /= npvt) CYCLE
   npivot= i
   GO TO 220
 END DO
 RETURN
 220 CALL gmmatd( c(18*npivot-17),6,3,1, kd(1),6,6,0, temp1(1))
 
!     TEMP1 NOW CONTAINS                   T
!                           ( (C )(E)(T ) ) ( KD)
!                               J      J
!     WHERE J IS THE PIVOT POINT
 
!     GENERATE THE THREE BY THREE PARTITIONS IN GLOBAL COORDINATES HERE
 
 DO  i=1,3
   CALL gmmatd( temp1,3,6,0, c(18*i-17),6,3,0,temp2(1)  )
!//////
!     CALL BUG(4HTRMK,260,TEMP2,18)
!//////
   DO  j=1,36
     kij(j) = 0.0D0
   END DO
   kij( 1) = temp2(1)
   kij( 2) = temp2(2)
   kij( 3) = temp2(3)
   kij( 7) = temp2(4)
   kij( 8) = temp2(5)
   kij( 9) = temp2(6)
   kij(13) = temp2(7)
   kij(14) = temp2(8)
   kij(15) = temp2(9)
   
   CALL ds1b( kij(1), necpt(i+1) )
   
 END DO
 RETURN
END SUBROUTINE dtrmem
