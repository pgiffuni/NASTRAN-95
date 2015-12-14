SUBROUTINE crdrd2 (*,*,mu,indcom,n23)
     
!     WRITE THE RIGID ROD ELEMENT ON THE RG FILE
 
!     EXTERNAL          ORF    ,LSHIFT
!     INTEGER           ORF
 
 , INTENT(IN OUT)                         :: *
 , INTENT(IN OUT)                         :: *
 INTEGER, INTENT(OUT)                     :: mu
 INTEGER, INTENT(IN OUT)                  :: indcom
 INTEGER, INTENT(OUT)                     :: n23
 INTEGER :: geomp  ,bgpdt  ,cstm   ,rgt    ,scr1   ,  &
     buf(20),mask16 ,gpoint ,z(1)   ,mcode(2)
 REAL :: rz(1)
 DOUBLE PRECISION :: indtfm(9),deptfm(9),rodcos(3),idrcos(3), ddrcos(3),  &
     dz(1)  ,xd     ,yd     ,zd     ,rlngth ,cdep
 COMMON /zzzzzz/   z
 COMMON /gp4fil/   geomp  ,bgpdt  ,cstm   ,rgt    ,scr1
 COMMON /gp4prm/   buf    ,buf1   ,buf2   ,buf3   ,buf4   ,knkl1  ,  &
     mask16 ,nogo   ,gpoint ,kn
 EQUIVALENCE       (z(1)  ,dz(1)) ,(z(1)  ,rz(1))
 DATA              mask15 /32767/
 
!     INDTFM = INDEPENDENT GRID POINT TRANSFORMATION MATRIX
!     DEPTFM = DEPENDENT GRID POINT TRANSFORMATION MATRIX
!     RODCOS = BASIC COSINES OF ROD ELEMENT
!     IDRCOS = DIRECTION COSINES OF INDEPENDENT GRID POINT
!     DDRCOS = DIRECTION COSINES OF DEPENDENT GRID POINT
 
!     OBTAIN TRANSFORMATION MATRIX
 
 IF (z(knkl1+3) == 0) GO TO 50
 DO  i = 1,4
   buf(i) = z(knkl1+2+i)
 END DO
 CALL transd (buf,indtfm)
 50 IF (z(knkl1+10) == 0) GO TO 70
 DO  i = 1,4
   buf(i) = z(knkl1+9+i)
 END DO
 CALL transd (buf,deptfm)
 
!     COMPUTE THE LENGTH OF THE RIGID ROD ELEMENT
 
 70 xd = rz(knkl1+11) - rz(knkl1+4)
 yd = rz(knkl1+12) - rz(knkl1+5)
 zd = rz(knkl1+13) - rz(knkl1+6)
 
!     CHECK TO SEE IF LENGTH OF ROD IS ZERO
 
 IF (xd == 0.0D0 .AND. yd == 0.0D0 .AND. zd == 0.0D0) RETURN 1
 rlngth = DSQRT(xd*xd + yd*yd + zd*zd)
 
!     COMPUTE THE BASIC DIRECTION COSINES OF THE RIGID ROD ELEMENT
 
 rodcos (1) = xd/rlngth
 rodcos (2) = yd/rlngth
 rodcos (3) = zd/rlngth
 
!     OBTAIN THE DIRECTION COSINES ASSOCIATED WITH
!     THE INDEPENDENT GRID POINT
 
 IF (z(knkl1+3) /= 0) GO TO 100
 DO  i = 1,3
   idrcos(i) = rodcos(i)
 END DO
 GO TO 200
 100 CALL gmmatd (rodcos,1,3,0,indtfm,3,3,0,idrcos)
 
!     OBTAIN THE DIRECTION COSINES ASSOCIATED WITH
!     THE DEPENDENT GRID POINT
 
 200 IF (z(knkl1+10) /= 0) GO TO 300
 DO  i = 1,3
   ddrcos(i) = rodcos(i)
 END DO
 GO TO 400
 300 CALL gmmatd (rodcos,1,3,0,deptfm,3,3,0,ddrcos)
 
!     DETERMINE THE DEPENDENT SIL AND THE CORRESPONDING COEFFICIENT
 
 400 DO  i = 1,3
   IF (indcom /= i) CYCLE
   idep = z(knkl1+6+i)
   cdep = rodcos(i)
   GO TO 600
 END DO
 
!     CHECK TO SEE IF RIGID ROD IS PROPERLY DEFINED
 
 600 IF (DABS(cdep) < 0.001D0) RETURN 2
 mcode(2) = idep
 IF (idep > mask15) n23 = 3
 DO  i = 1, 3
   mcode(1) = z(knkl1+i-1)
   IF (mcode(1) > mask15) n23 = 3
   coeff = -idrcos(i)/cdep
   CALL WRITE (rgt,mcode,2,0)
   CALL WRITE (rgt,coeff,1,0)
   mcode(1) = z(knkl1+6+i)
   IF (mcode(1) > mask15) n23 = 3
   coeff = ddrcos(i)/cdep
   CALL WRITE (rgt,mcode,2,0)
   CALL WRITE (rgt,coeff,1,0)
 END DO
 z(mu) = idep
 mu = mu - 1
 RETURN
END SUBROUTINE crdrd2
