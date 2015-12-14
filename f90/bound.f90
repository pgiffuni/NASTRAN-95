SUBROUTINE bound (fbrec,afe,nafe,kge,nkge)
     
!     COMPUTES AREA FACTOR AND GRAVITIONAL STIFFNESS MATRICES FOR A FACE
!     OF A INDIVIDUAL FLUID ELEMENT
 
 
 INTEGER, INTENT(IN)                      :: fbrec(12)
 DOUBLE PRECISION, INTENT(OUT)            :: afe(48)
 INTEGER, INTENT(OUT)                     :: nafe
 DOUBLE PRECISION, INTENT(OUT)            :: kge(144)
 INTEGER, INTENT(OUT)                     :: nkge
 LOGICAL :: error    ,grav
 INTEGER :: gf1      ,gf2      ,gf3      ,gf4      ,  &
     gs1      ,gs2      ,gs3      ,gs4      ,grid(3,4),  &
     gsi      ,gsj      ,iz(1)    ,locsof(4),loctof(3),  &
     locfos(4),fledge(2,4)        ,fedge(2,4)         , stedge(2,3)
 REAL :: z
 DOUBLE PRECISION :: in(3)    ,jn(3)    ,  &
     kn(3)    ,r12(3)   ,r13(3)   ,r14(3)   ,r24(3)   ,  &
     h        ,nn       ,ks(3)    ,mag      ,x3       ,  &
     y3       ,y4       ,s(48)    ,x4       ,akj(3,4) ,  &
     aa       ,bb       ,cc       ,a        ,zz       ,  &
     dvmag    ,rhoxg    ,y(3)     ,e(3,2)   ,kii(144) ,  &
     ktwo(2,2),kik(9)   ,t(3,3)   ,ktemp(2,3)         ,  &
     tfst(3,3),z1       ,x1       ,dhalf    ,c1       ,  &
     st(3,4)  ,z2       ,y1       ,eps(2)   ,c2       ,  &
     fl(3,4)  ,z3       ,x2       ,dlb      ,c3       ,  &
     tr(3,3)  ,z4       ,y2       ,dub      ,d1       ,  &
     p(2,7)   ,ss(9)    ,aa2      ,nn1      ,d2       ,  &
     c(4,7)   ,epslon   ,fdet     ,zz1      ,dz       ,  &
     f(3,7)   ,nx       ,akjcon   ,dd       ,aeps     ,  &
     pt(3,4)  ,tria     ,epso10   ,aflel    ,leps     ,  &
     vtemp(3) ,ksb(3)   ,nz       ,kident(3),factii   ,  &
     conii    ,fii      ,dpoly    ,astria   ,astrel   , aflstr   ,dadotb   ,dapoly
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm      ,uwm
 
!     OPEN CORE
 
 COMMON /zzzzzz/ z(1)
 
!     CORE POINTERS
 
 COMMON /flbptr/ error    ,icore    ,lcore    ,ibgpdt   ,nbgpdt   ,  &
     isil     ,nsil     ,igrav    ,ngrav
 
!     MATERIAL PROPERTIES
 
 COMMON /matin / matid    ,inflag
 COMMON /matout/ dum(3)   ,rho
 
!     MODULE PARAMETERS
 
 COMMON /BLANK / nograv
 
!     NASTRAN PARAMETERS
 
 COMMON /system/ sysbuf   ,nout
 EQUIVALENCE     (tfst(1,1),in(1))  ,(x1,fl(1,1))  ,(x2,fl(1,2))  ,  &
     (tfst(1,2),jn(1))  ,(y1,fl(2,1))  ,(y2,fl(2,2))  ,  &
     (tfst(1,3),kn(1))  ,(z1,fl(3,1))  ,(z2,fl(3,2))  ,  &
     (ss(1),bb)         ,(x3,fl(1,3))  ,(x4,fl(1,4))  ,  &
     (ss(2),cc)         ,(y3,fl(2,3))  ,(y4,fl(2,4))  ,  &
     (ss(3),zz)         ,(z3,fl(3,3))  ,(z4,fl(3,4))  ,  &
     (ss(4),nn)         ,(ss(5),nn1)   ,(ss(6),zz1 )  ,  &
     (eps(1),aeps)      ,(eps(2),leps) ,(fii,bb    )  ,  &
     (factii,cc)        ,(conii,akjcon),(z(1),iz(1))
 
!     GRID POINTS TO BE USED IN SUBDIVIDING QUADS INTO TRIANGLES
 
 DATA    grid  /  1     ,2     ,3   , 2     ,3     ,4   ,  &
     3     ,4     ,1   , 4     ,1     ,2   /
 
 DATA    dz, d1,  d2, dhalf / 0.d0,  1.d0, 2.d0, .5D0 /
 DATA    epslon,  epso10    / 1.d-3, 1.d-4            /
 DATA    dlb   ,  dub       /-1.d-3, 1.001D0          /
 DATA    x1, x2,  y1, y2, z1, z2, z3, z4 / 8*0.d0     /
 
 DATA    fedge /  1,2, 2,3, 3,4, 4,1 /
 DATA    stedge/  1,2, 2,3, 3,1      /
 DATA    kident/  0.d0, 0.d0, 1.d0   /
 
 
!     DETERMINE SIZES OF MATRIX PARTITIONS
 
 ngrids = 4
 IF (fbrec( 6) < 0) ngrids = 3
 ngridf = 4
 IF (fbrec(12) < 0) ngridf = 3
 
 nrow = 3*ngrids
 nafe = nrow*ngridf*2
 nkge = 0
 
!     OBTAIN MATERIAL PROPERTY AND GRAVITY DATA IF GRAV ID IS
!     PRESENT
 
 grav   = .false.
 IF (fbrec(7) == 0) GO TO 600
 inflag = 11
 matid  = fbrec(8)
 CALL mat (fbrec(1))
 
 IF (ngrav == 0) GO TO 8013
 lgrav = igrav + ngrav - 1
 DO  i = igrav,lgrav,6
   IF (iz(i) == fbrec(7)) GO TO 400
 END DO
 
 GO TO 8013
 
 400 g     = SQRT(z(i+3)**2 + z(i+4)**2 + z(i+5)**2)
 g     = g*z(i+2)
 rhoxg = DBLE(rho)*DBLE(g)
 nkge  = nrow*nrow*2
 nograv= 1
 grav  = .true.
 
!     NORMILIZE THE GRAVITY VECTOR
 
 e(1,2) = DBLE(z(i+3))
 e(2,2) = DBLE(z(i+4))
 e(3,2) = DBLE(z(i+5))
 CALL dnorm (e(1,2),mag)
 IF (iz(i+1) == 0) GO TO 600
 
!     TRANSFORM GRAVITY VECTOR TO BASIC
 
 j = iz(ibgpdt)
 iz(ibgpdt) = iz(i+1)
 CALL transd (iz(ibgpdt),tr)
 iz(ibgpdt) = j
 CALL gmmatd (tr,3,3,0,e(1,2),3,1,0,vtemp)
 DO  j = 1,3
   e(j,2) = vtemp(j)
 END DO
 
 
!     COMPUTE NEW COORDINATES FOR FLUID FACE BASED ON FLUID COORDINATE
!     SYSTEM - PERFORM THIS ONLY IF THE FLUID FACE HAS CHANGED
!     THESE COMPUTATIONS INCLUDE --
 
!        IN,JN,KN  - NORMAL VECTORS TO DEFINE FLUID COORDINATE SYSTEM
!        X2,X3,X4  - X COORDINATES OF GRID POINTS IN NEW SYSTEM
!                    ( X1 = 0 )
!        Y3,Y4       Y COORDINATES OF GRID POINTS IN NEW SYSTEM
!                    ( Y1,Y2 = 0 )
 
!     NORMAL (UNIT) VECTORS STORED *COLUMN-WISE* IN U --
!           I IN U(L,1), J IN U(L,2), K IN U(L,3), L= 1,3
!        TRANSFORMED FLUID COORDINATES STORED IN FL
 
 
!     LOCATE GRID POINTS COORDINATES FOR THE FLUID GRID POINTS IN THE
!     BGPDT TABLE
 
 600 gf1 = ibgpdt + (fbrec( 9)-1)*4
 gf2 = ibgpdt + (fbrec(10)-1)*4
 gf3 = ibgpdt + (fbrec(11)-1)*4
 gf4 = -1
 IF (ngridf == 4) gf4 = ibgpdt + (fbrec(12)-1)*4
 
 IF (ngridf == 4) GO TO 700
 
!     TRIANGULAR FLUID FACE
 
 DO  i = 1,3
   r12(i) = z(gf2+i) - z(gf1+i)
   in(i)  = r12(i)
   r13(i) = z(gf3+i) - z(gf1+i)
 END DO
 
 CALL dnorm (in,mag)
 x2 = mag
 
 CALL daxb  (r12,r13,kn)
 CALL dnorm (kn,mag)
 
 CALL daxb (kn,in,jn)
 
 x3 = r13(1)*in(1) + r13(2)*in(2) + r13(3)*in(3)
 y3 = r13(1)*jn(1) + r13(2)*jn(2) + r13(3)*jn(3)
 GO TO 1000
 
!     QUADRATIC FLUID FACE
 
 700 DO  i = 1,3
   r12(i) = z(gf2+i) - z(gf1+i)
   r13(i) = z(gf3+i) - z(gf1+i)
   r14(i) = z(gf4+i) - z(gf1+i)
   r24(i) = z(gf4+i) - z(gf2+i)
 END DO
 
 CALL daxb  (r13,r24,kn)
 CALL dnorm (kn,mag)
 
 h = r12(1)*kn(1) + r12(2)*kn(2) + r12(3)*kn(3)
 
 DO  i = 1,3
   in(i) = r12(i) - h*kn(i)
 END DO
 CALL dnorm (in,mag)
 
 x2 = mag
 
 CALL daxb (kn,in,jn)
 
 x3 = r13(1)*in(1) + r13(2)*in(2) + r13(3)*in(3)
 x4 = r14(1)*in(1) + r14(2)*in(2) + r14(3)*in(3)
 y3 = r13(1)*jn(1) + r13(2)*jn(2) + r13(3)*jn(3)
 y4 = r14(1)*jn(1) + r14(2)*jn(2) + r14(3)*jn(3)
 
!     VARIOUS CALCULATIONS DEPENDENT ON FLUID FACE
 
!     INDICES FOR CORNERS OF FLUID ELEMENT
 
 1000 DO  n = 1,2
   DO  j = 1,ngridf
     fledge(n,j) = fedge(n,j)
   END DO
 END DO
 fledge(2,ngridf) = 1
 
!     SET UP FOR FLUID TRIANGLE
 
 c1 = (d1 - fl(1,3)/fl(1,2))/fl(2,3)
 c2 = fl(1,3)/(fl(1,2)*fl(2,3))
 DO  n = 1,3
   r12(n) = fl(n,2) - fl(n,1)
   r13(n) = fl(n,3) - fl(n,1)
 END DO
 CALL daxb (r12,r13,vtemp)
 
 IF (ngridf == 3) GO TO 1040
 
!     SET UP FOR FLUID QUADRANGLE
 
 c1  = fl(2,3) - fl(2,4)
 c2  = fl(1,2)*fl(2,4)
 c3  = fl(1,2) - fl(1,3) + fl(1,4)
 aa  =-fl(1,2)*c1
 aa2 = d2*aa
 
 DO  n = 1,3
   r13(n) = fl(n,3) - fl(n,1)
   r24(n) = fl(n,4) - fl(n,2)
 END DO
 CALL daxb (r13, r24,vtemp)
 1040 aflel = dvmag(vtemp,dz)
 
!     ZERO OUT AREA FACTOR MATRIX
!     AND AREA COMMON TO FLUID AND STRUCTURE ELEMENTS (AFLSTR)
 
 DO  i = 1,48
   afe(i) = dz
   s(i) = 0.0D0
 END DO
 DO  i = 1,144
   kge(i) = 0.0D0
 END DO
 aflstr = 0.0
 
!     DETERMINE NUMBER OF STRUCTURAL TRIANGLES TO BE USED, ITRIA
!     AND CUMULATIVE AREA CONSTANT, TRIA
!        ITRIA= 4, TRIA= .5 WHEN STRUCTURE ELEMENT IS QUADRANGLE
!        ITRIA= 1, TRIA= 1. WHEN STRUCTURE ELEMENT IS TRIANGLE
 
 itria = 1
 tria  = d1
 IF (ngrids == 3) GO TO 1050
 itria = 4
 tria  = dhalf
 
!     TRANSFORM STRUCTURE COORDINATES TO FLUID COORDINATE SYSTEM
 
 1050 gs1 = ibgpdt + (fbrec(3)-1)*4
 gs2 = ibgpdt + (fbrec(4)-1)*4
 gs3 = ibgpdt + (fbrec(5)-1)*4
 gs4 = -1
 IF (ngrids == 4) gs4 = ibgpdt + (fbrec(6)-1)*4
 
 DO  n = 1,3
   pt(n,1) = z(gs1+n) - z(gf1+n)
   pt(n,2) = z(gs2+n) - z(gf1+n)
   pt(n,3) = z(gs3+n) - z(gf1+n)
   pt(n,4) = dz
   IF (ngrids == 4) pt(n,4) = z(gs4+n) - z(gf1+n)
   DO  k = 1,4
     st(n,k) = dz
   END DO
 END DO
 
 DO  k = 1,ngrids
   DO  n = 1,3
     DO  m = 1,3
       st(n,k) =  st(n,k) + pt(m,k)*tfst(m,n)
     END DO
   END DO
 END DO
 DO  n = 1,2
   r12(n) = st(n,2) - st(n,1)
   r13(n) = st(n,3) - st(n,1)
   IF (ngrids == 4) r24(n) = st(n,4) - st(n,2)
 END DO
 CALL daxb (r12,r13,vtemp)
 IF (ngrids == 4) CALL daxb (r12,r24,vtemp)
 astrel = dvmag(vtemp,dz)
 aeps   = dhalf*DMIN1(aflel,astrel)
 leps   = dz
 IF (aeps > dz) leps = epslon*DSQRT(aeps)
 aeps   = epslon*aeps
 
!     LOCATE STRUCTURE ELEMENT GRIDS RELATIVE TO FLUID SURFACE
!     LOCSOF FLAGS STRUCTURE ON FLUID:
!            1= INSIDE, -1= OUTSIDE, 0= ON FLUID EDGE
 
 CALL locpt (ngrids,st,ngridf,fl,fledge,kident,eps,locsof)
 
 
!     LOOP THRU (INCREMENTAL) STRUCTURAL TRIANGLES (ITRIA IS 1 OR 4)
 
 DO  it = 1,itria
   
!     LOCATE COORDINATES OF CURRENT TRIANGLE
   
   gs1 = grid(1,it)
   gs2 = grid(2,it)
   gs3 = grid(3,it)
   
   loctof(1) = locsof(gs1)
   loctof(2) = locsof(gs2)
   loctof(3) = locsof(gs3)
   
!     TRANSFER COORDINATES OF CURRENT STRUCTURE TRIANGLE TO CONTIGUOUS
!     ARRAY, AND DO VARIOUS CALCULATIONS DEPENDENT ON THEM
   
   DO  n = 1,3
     tr(n,1) = st(n,gs1)
     tr(n,2) = st(n,gs2)
     tr(n,3) = st(n,gs3)
     r12(n)  = tr(n,2) - tr(n,1)
     r13(n)  = tr(n,3) - tr(n,1)
   END DO
   
!     OBTAIN KS, UNIT VECTOR NORMAL TO (XY) PLANE OF CURRENT STRUCTURAL
!     TRIANGLE (IN SYSTEM LOCAL TO FLUID ELEMENT)
   
   CALL daxb (r12,r13,ks)
   astria = dvmag(ks,dz)
   CALL dnorm (ks,mag)
   
!     OBTAIN KSB, UNIT VECTOR NORMAL TO (XY) PLANE OF CURRENT STRUCTURE
!     TRIANGLE (IN NASTRAN BASIC COORD SYSTEM)
   
   DO  n = 1,3
     r12(n) = pt(n,gs2) - pt(n,gs1)
     r13(n) = pt(n,gs3) - pt(n,gs1)
   END DO
   
   CALL daxb  (r12,r13,ksb)
   CALL dnorm (ksb,mag)
   
!     CALCULATE EPSLON FUNCTIONS FOR SIGNIFICANCE TESTING
   
   leps = dz
   aeps = dhalf*DMIN1(aflel,astria)*epslon
   IF (aeps > dz)  leps = DSQRT(aeps)
   
!     DETERMINE POINTS DESCRIBING AREA POLYGON COMMON TO BOTH FLUID
!     ELEMENT AND (INCREMENTAL) STRUCTURAL TRIANGLE
   
!        POLYGON POINTS IN   P(2,I)    I .LE. 7
!        FLUID POINTS IN     FL(3,J)   J .LE. 4
!        TRIANGLE POINTS IN  TR(3,K)   K=1,3
   
!     DETERMINE POINTS DESCRIBING POLYGON OF COMMON AREA
   
   
!     LOCATE FLUID ELEMENT POINTS RELATIVE TO BOUNDRY OF THIS STRUCTURAL
!     TRIANGLE
   
   CALL locpt (ngridf,fl,3,tr,stedge,ks,eps,locfos)
   DO  j = 1,ngridf
     IF (locfos(j) < 0) GO TO 1300
   END DO
   
!     FLUID ELEMENT IS COMMON AREA POLYGON WHEN NO FLUID POINTS ARE
!     OUTSIDE BOUNDRY OF THIS STRUCTURAL TRIANGLE
   
   npoly = ngridf
   DO  n = 1,2
     DO  j = 1,ngridf
       p(n,j) = fl(n,j)
     END DO
   END DO
   GO TO 2000
   
!     CALL POLYPT TO DETERMINE POINTS DESCRIBING THE COMMON AREA POLYGON
   
   1300 CALL polypt (loctof,stedge,tr,ngridf,fledge,fl,locfos,eps,npoly,p)
   
!     SKIP TO NEXT (INCREMENTAL) STRUCTURAL TRIANGLE WHEN THIS TRIANGLE
!     IS DISJOINT FROM FLUID ELEMENT
   
   IF (npoly < 3) CYCLE
   
!     AREA OF COMMON POLYGON AND HALVED WHEN OVERLAPPING (INCREMENTAL)
!     STRUCTURE TRIANGLES USED CUMULATIVE AREA OF FLUID/STRUCTURAL
!     ELEMENT OVERLAP
   
   2000 a = tria*dapoly(npoly,p)
   aflstr = aflstr + a
   
!     TERMS FOR LOAD FACTORS
   
   ss(1) =  tr(1,1)*tr(2,2)
   ss(2) = -tr(1,1)*tr(2,3)
   ss(3) =  tr(1,2)*tr(2,3)
   ss(4) = -tr(1,2)*tr(2,1)
   ss(5) =  tr(1,3)*tr(2,1)
   ss(6) = -tr(1,3)*tr(2,2)
   fdet  =  dz
   DO  m = 1,6
     fdet  =  fdet  + ss(m)
   END DO
   ss(1) =  ss(1) + ss(4)
   ss(2) =  ss(2) + ss(5)
   ss(3) =  ss(3) + ss(6)
   ss(4) =  tr(2,2) - tr(2,3)
   ss(5) =  tr(2,3) - tr(2,1)
   ss(6) =  tr(2,1) - tr(2,2)
   ss(7) =  tr(1,3) - tr(1,2)
   ss(8) =  tr(1,1) - tr(1,3)
   ss(9) =  tr(1,2) - tr(1,1)
   
!     GET LOAD DISTRIBUTION FACTORS, F(K,I)
!     - FROM -
!         I -- AREA POLYGON POINT -- P(N,I)
!         K -- STRUCTURE TRIANGLE POINT -- TR(N,K)
   
   DO  i = 1,npoly
     f(1,i) = p(1,i)*ss(4) + p(2,i)*ss(7) + ss(3)
     f(2,i) = p(1,i)*ss(5) + p(2,i)*ss(8) + ss(2)
     f(3,i) = p(1,i)*ss(6) + p(2,i)*ss(9) + ss(1)
   END DO
   
!     GET PRESSURE DISTRIBUTION FACTORS, C(J,I)
!     - FROM -
!         I -- AREA POLYGON POINT  -- P(N,I)
!         J -- FLUID ELEMENT POINT -- FL(N,J)
   
   IF (ngridf == 4) GO TO 2030
   
!     FLUID ELEMENT IS TRIANGLE
   
   DO  i = 1,npoly
     bb     = p(1,i)/fl(1,2)
     c(1,i) = d1 - bb - p(2,i)*c1
     c(2,i) = bb - p(2,i)*c2
     c(3,i) = p(2,i)/fl(2,3)
   END DO
   GO TO 2100
   
!     FLUID ELEMENT IS QUADRANGLE
   
   2030 DO  i = 1,npoly
     bb = p(1,i)*c1 - c2 + p(2,i)*c3
     cc = p(1,i)*fl(2,4) - p(2,i)*fl(1,4)
     IF (bb == dz .OR. DABS(aa) > DABS(bb*epslon)) GO TO 2040
     zz = -cc/bb
     GO TO 2045
     
     2040 dd = DSQRT(bb*bb - d2*aa2*cc)
     zz = (dd-bb)/aa2
     IF (zz > dlb .AND. zz < dub) GO TO 2045
     zz = (-dd-bb)/aa2
     
     2045 nn = p(2,i)/(fl(2,4) + zz*c1)
     IF (nn <= dlb .OR. nn >= dub) GO TO 8005
     
     zz1 = d1 - zz
     nn1 = d1 - nn
     c(1,i) = zz1*nn1
     c(2,i) = zz *nn1
     c(3,i) = zz *nn
     c(4,i) = zz1*nn
   END DO
   
!     CALCULATE AREA TERMS FOR THIS STRUCTURAL TRIANGLE AND INSERT IN
!     MATRIX
   
   2100 dpoly  = npoly
   akjcon = a/(fdet*dpoly)
   dpoly  = npoly - 1
   factii = d1/dpoly
   
   DO  j = 1,ngridf
     jloc = 3*ngrids*(j-1)
     
     DO  k = 1,3
       loc = jloc + 3*(grid(k,it)-1)
       
       akj(k,j) = dz
       DO  i = 1,npoly
         akj(k,j) = akj(k,j) + f(k,i)*c(j,i)
       END DO
       akj(k,j) = akjcon*akj(k,j)
       
       DO  n = 1,3
         s(loc+n) = s(loc+n) + akj(k,j)*ksb(n)
       END DO
     END DO
   END DO
   
   IF (.NOT. grav) CYCLE
   
!     CALCULATE GRAVITATIONAL STIFFNESS TERMS FOR THIS TRIANGLE
!     AND INSERT INTO MATRIX
   
   DO  n = 1,3
     e(n,1) = dz
   END DO
   CALL daxb (e(1,2),ksb,y)
   mag = dadotb(y,y)
   IF (mag > dz) mag = DSQRT(mag)
   IF (mag < epso10) GO TO 2220
   
   CALL daxb  (e(1,2),y,e)
   CALL dnorm (e,mag)
   
   2220 nx = 0.d0
   nz = 0.d0
   DO  n = 1,3
     nx = nx + e(n,1)*ksb(n)
     nz = nz + e(n,2)*ksb(n)
   END DO
   conii = rhoxg*akjcon/(d2*fdet)
   ktwo(1,1) = dz
   
   ktwo(2,1) = nx
   ktwo(1,2) = ktwo(2,1)
   ktwo(2,2) = nz
   CALL gmmatd (e,2,3,1, ktwo,2,2,0, ktemp)
   CALL gmmatd (ktemp,3,2,0, e,2,3,0, kik )
   
   DO  kk1 = 1,3
     k1loc = 9*ngrids*(grid(kk1,it)-1)
     
     DO  kk2 = 1,3
       loc = k1loc + 9*(grid(kk2,it)-1)
       
       h = 0.d0
       DO  i1 = 1,npoly
         DO  i2 = 1,npoly
           fii = f(kk1,i1)*f(kk2,i2)
           IF (i1 /= i2) fii= factii*fii
           h = h + fii
         END DO
       END DO
       
       DO  n = 1,9
         kge(loc+n) = kge(loc+n) - kik(n)*h*conii
       END DO
       
     END DO
   END DO
   
!     END OF (INCREMENTAL) STRUCTURAL TRIANGLE LOOP
   
 END DO
 
!     WARNING MESSAGE WHEN FLUID AND STRUCTURE ELEMENTS ARE DISJOINT
 
 IF (aflstr <= dz) GO TO 8014
 
!     TRANSFORM THE AREA AND STIFFNESS MATRICES TO GLOBAL COORDINATES IF
!     REQUIRED
 
 DO  irow = 1,ngrids
   gsi = ibgpdt + (fbrec(irow+2)-1)*4
   CALL transd (z(gsi),t)
   
!     AREA FACTOR MATRIX
   
   jloc = 3*(irow-1)
   
   DO  icol = 1,ngridf
     iloc = 3*ngrids*(icol-1) + jloc
     
     IF (iz(gsi) == 0) GO TO 2510
     CALL gmmatd (t,3,3,1,s(iloc+1),3,1,0,afe(iloc+1))
     CYCLE
     
     2510 DO  i = 1,3
       afe(iloc+i) = s(iloc+i)
     END DO
     
   END DO
   IF (.NOT.grav) CYCLE
   
!     GRAVITATIONAL STIFFNESS MATRIX
   
   jloc = 9*(irow-1)
   
   DO  icol = 1,ngrids
     iloc = 9*ngrids*(icol-1) + jloc
     
     IF (iz(gsi) == 0) GO TO 2540
     CALL gmmatd (t,3,3,1, kge(iloc+1),3,3,0, kik)
     GO TO 2570
     
     2540 kloc = iloc
     DO  i = 1,9
       kik(i) = kge(kloc+i)
     END DO
     
     2570 gsj = ibgpdt + (fbrec(icol+2)-1)*4
     IF (iz(gsj) == 0) GO TO 2580
     CALL transd (z(gsj),t)
     CALL gmmatd (kik,3,3,0, t,3,3,0, kii(iloc+1))
     CYCLE
     
     2580 kloc = iloc
     DO  i = 1,9
       kii(kloc+i) = kik(i)
     END DO
   END DO
 END DO
 
!     REARANGE THE STORAGE OF THE GRAVITATIONAL STIFFNESS MATRIX
!     TO COLUMNWISE FOR THE USE WITH THE ASSEMBLER
 
 IF (.NOT.grav) RETURN
 
 DO  icol = 1,ngrids
   jloc = 9*ngrids*(icol-1)
   
   DO  irow = 1,ngrids
     iloc = jloc + 9*(irow-1)
     kloc = jloc + 3*(irow-1)
     
     DO  i = 1,3
       kge(kloc+1) = kii(iloc+1)
       kge(kloc+2) = kii(iloc+4)
       kge(kloc+3) = kii(iloc+7)
       kloc = kloc + 3*ngrids
       iloc = iloc + 1
     END DO
   END DO
 END DO
 RETURN
 
!     ERROR CONDITIONS
 
 8005 WRITE (nout,9005) ufm,fbrec(2)
 error = .true.
 GO TO 9000
 8013 WRITE (nout,9013) ufm,fbrec(1),fbrec(7)
 error = .true.
 GO TO 9000
 8014 WRITE (nout,9014) uwm,fbrec(1),fbrec(2)
 9000 RETURN
 
 9005 FORMAT (a23,' 8005. BAD GEOMETRY DEFINED FOR STRUCTURAL ELEMENT ', i8)
 
 9013 FORMAT (a23,' 8013, FLUID ELEMENT',i9,' ON A CFLSTR CARD ',  &
     'REFERENCES UNDEFINED GRAVITY ID',i9)
 
 9014 FORMAT (a25,' 8014, FLUID ELEMENT',i9,' AND STRUCTURE ELEMENT',i9,  &
     ' ARE DISJOINT. CHECK CFLSTR CARDS.')
END SUBROUTINE bound
