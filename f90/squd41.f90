SUBROUTINE squd41
     
!     PHASE 1  STRESS DATA RECOVERY FOR CQUAD4 ELEMENT
 
!                EST  LISTING
 
!     WORD       TYPE         DESCRIPTION
!     --------------------------------------------------------------
!       1          I    ELEMENT ID, EID
!       2 THRU 5   I    SILS, GRIDS 1 THRU 4
!       6 THRU 9   R    MEMBRANE THICKNESSES T AT GRIDS 1 THRU 4
!      10          R    MATERIAL PROPERTY ORIENTATION ANGLE, THETA
!               OR I    COORD. SYSTEM ID (SEE TM ON CQUAD4 CARD)
!      11          I    TYPE FLAG FOR WORD 10
!      12          R    GRID ZOFF  (OFFSET)
!      13          I    MATERIAL ID FOR MEMBRANE, MID1
!      14          R    ELEMENT THICKNESS, T (MEMBRANE, UNIFORMED)
!      15          I    MATERIAL ID FOR BENDING, MID2
!      16          R    BENDING INERTIA FACTOR, I
!      17          I    MATERIAL ID FOR TRANSVERSE SHEAR, MID3
!      18          R    TRANSV. SHEAR CORRECTION FACTOR TS/T
!      19          R    NON-STRUCTURAL MASS, NSM
!      20 THRU 21  R    Z1, Z2  (STRESS FIBRE DISTANCES)
!      22          I    MATERIAL ID FOR MEMBRANE-BENDING COUPLING, MID4
!      23          R    MATERIAL ANGLE OF ROTATION, THETA
!               OR I    COORD. SYSTEM ID (SEE MCSID ON PSHELL CARD)
!      24          I    TYPE FLAG FOR WORD 23
!      25          I    INTEGRATION ORDER
!      26          R    STRESS ANGLE OF ROTATION, THETA
!               OR I    COORD. SYSTEM ID (SEE SCSID ON PSHELL CARD)
!      27          I    TYPE FLAG FOR WORD 26
!      28          R    ZOFF1 (OFFSET)  OVERRIDDEN BY EST(12)
!      29 THRU 44  I/R  CID,X,Y,Z - GRIDS 1 THRU 4
!      45          R    ELEMENT TEMPERATURE
 
 
 LOGICAL :: badjac,membrn,bendng,shrflx,mbcoup,norpth,nocsub
 INTEGER :: nest(45),nphi(2395),sil(4),ksil(4),kcid(8),  &
     igpdt(4,4),elid,scsid,flags,flagm,necpt(4),  &
     INDEX(3,3),mid(4),q4strs,ipn(4),hunmeg,rowflg, TYPE,NAME(2)
 REAL :: bgpdm(3,4),cent(3),gpth(4),gpnorm(4,4),bgpdt(4,4),  &
     matset,mominr,tmpthk(4),tgrid(4,4),epnorm(4,4),  &
     egpdt(4,4),g(6,6),gi(36),shp(4),dshp(8),gge(9),  &
     ggu(9),ptint(2),ptintp(3),tbs(9),teu(9),tse(9),  &
     teb(9),tbg(9),tub(9),tum(9),tsu(9),u(9),GT(9),  &
     tbm(9),tem(9),tmi(9),ecpt(4),gpc(3),xa(4),yb(4),  &
     alfa(3),gpth2(4),relout(300),nunorx,nunory,  &
     ugpdm(3,4),cente(3),bmatrx(192),xybmat(96),  &
     jacob(3,3),phi(9),psitrn(9),tmpshp(4),dshptp(8),  &
     kheat,tms(9),dq(24),jacobu(9),jacbs(9),jacobe(9), zc(4),vnt(3,4)
!WKBNB 11/93 SPR 93020
 REAL :: vd1(3), vd2(3), vkn(3), vks(3)  &
     ,               v12(3), v41(3), vp12(3),vis(3), vjs(3)
!WKBNE 11/93 SPR 93020
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /sdr2x5/ est(100),phiout(2395)
 COMMON /sdr2x6/ ielout(300)
 COMMON /condas/ pi,twopi,raddeg,degrad
 COMMON /system/ systm(100)
 COMMON /matin / matid,inflag,eltemp
 COMMON /matout/ rmtout(25)
 COMMON /q4dt  / detj,hzta,psitrn,nnode,badjac,node
 COMMON /terms / membrn,bendng,shrflx,mbcoup,norpth
 COMMON /hmtout/ kheat(7),TYPE
 COMMON /q4coms/ anglei(4),edgshr(3,4),edgel(4),unv(3,4),  &
     uev(3,4),rowflg,iorder(4)
 EQUIVALENCE     (igpdt(1,1),bgpdt(1,1)),(est(1)  ,nest(1)   ),  &
     (bgpdt(1,1),est(29)   ),(gpth(1) ,est(6)    ),  &
     (elth      ,est(14)   ),(sil(1)  ,nest(2)   ),  &
     (nphi(1)   ,phiout(1) ),(INT     ,nest(25)  ),  &
     (zoff      ,nest(12)  ),(zoff1   ,est(28)   ),  &
     (ielout(1) ,relout(1) ),(matset  ,rmtout(25)),  &
     (necpt(1)  ,ecpt(1)   ),(systm(2),nout      ),  &
     (phiout(65),gpth2(1)  ),(systm(3),nogo      ),  &
     (htcp      ,kheat(4)  ),(itherm  ,systm(56) )
 DATA    eps1  / 1.0E-16/   ,ipn / 1,4,2,3 /
 DATA    NAME  / 4HQUAD,4H4      /
 DATA    hunmeg/ 100000000       /
 DATA    const / 0.57735026918962/
 
!     PHIOUT DATA BLOCK
!     --------------------------------------------------------------
!     PHIOUT(1)                 = ELID (ELEMENT ID)
!     PHIOUT(2-9)               = SIL NUMBERS
!     PHIOUT(10-17)             = ARRAY IORDER
!     PHIOUT(18)                = TSUB0 (REFERENCE TEMP.)
!     PHIOUT(19-20)             = Z1 & Z2 (FIBER DISTANCES)
!     PHIOUT(21)                = AVGTHK  (AVERAGE THICKNESS)
!     PHIOUT(22)                = MOMINR  (MOMENT OF INER. FACTOR)
!     PHIOUT(23-58)             = GBAR (BASIC MAT. PROP. MATRIX)
!                                 (W/O SHEAR)
!     PHIOUT(59-61)             = THERMAL EXPANSION COEFFICIENTS
!                                 FOR MEMBRANE MATERIAL
!     PHIOUT(62-64)             = THERMAL EXPANSION COEFFICIENTS
!                                 FOR BENDING MATERIAL
!     PHIOUT(65-68)             = CORNER NODE THICKNESSES
!     PHIOUT(69-77)             = 3X3 TRANSFORMATION FROM USER TO
!                                 MATERIAL COORD. SYSTEM
!     PHIOUT(78)                = OFFSET OF ELEMENT FROM GP PLANE
!     PHIOUT(79)                = ID OF THE ORIGINAL PCOMP(I)
!                                 PROPERTY ENTRY FOR COMPOSITES
!     PHIOUT(80-(79+9*NNODE))   = 3X3 TRANSFORMATIONS FROM GLOBAL
!                                 TO ELEMENT COORDINATE SYSTEM
!                                 FOR EACH EXISTING NODE
 
!     THE FOLLOWING IS REPEATED FOR EACH EVALUATION POINT AND THE
!     CENTER POINT (10 TIMES). THE EVALUATION POINTS ARE AT THE
!     STANDARD 2X2X2 GAUSSIAN POINTS. THE CHOICE OF THE
!     FINAL STRESS AND FORCE OUTPUT POINTS IS MADE AT THE SUBCASE
!     LEVEL (PHASE 2.)
 
!              1                  THICKNESS OF THE ELEMENT AT THIS
!                                 EVALUATION POINT
!            2 - 10               3X3 TRANSFORMATION FROM TANGENT
!                                 TO STRESS C.S. AT THIS EVAL. PT.
!           11 - 19               CORRECTION TO GBAR-MATRIX FOR
!                                 MEMBRANE-BENDING COUPLING AT THIS
!                                 EVALUATION POINT
!           20 - 28               3X3 TRANSFORMATION FROM MATERIAL
!                                 TO INTEGRATION PT. COORDINATE
!                                 SYSTEM
!           29 - 32               2X2 PROPERTY MATRIX FOR OUT-OF-
!                                 PLANE SHEAR (G3)
!         32+1 - 32+NNODE         ELEMENT SHAPE FUNCTIONS
!   32+NNODE+1 - 32+NNODE+8*NDOF  STRAIN RECOVERY MATRIX
 
 
!              IELOUT DATA BLOCK      (TOTAL OF NWORDS = 102)
!     --------------------------------------------------------------
!              1                  ELEMENT ID
!              2                  AVERAGE THICKNESS
 
!     THE FOLLOWING IS REPEATED FOR EACH CORNER POINT.
 
!         WORD  1                 SIL NUMBER
!         WORD  2-10              TBS TRANSFORMATION FOR Z1
!         WORD 11-19              TBS TRANSFORMATION FOR Z2
!         WORD 20-22              NORMAL VECTOR IN BASIC C.S.
!         WORD 23-25              GRID COORDS IN BASIC C.S.
 
 
 q4strs = 0
 elid   = nest(1)
 nphi(1)= elid
 norpth =.false.
 node   = 4
 nnode  = 4
 ndof   = nnode*6
 nd2    = ndof*2
 nd3    = ndof*3
 nd4    = ndof*4
 nd5    = ndof*5
 nd6    = ndof*6
 nd7    = ndof*7
 nd8    = ndof*8
 
!     FILL IN ARRAY GGU WITH THE COORDINATES OF GRID POINTS 1, 2 AND 4.
!     THIS ARRAY WILL BE USED LATER TO DEFINE THE USER COORD. SYSTEM
!     WHILE CALCULATING  TRANSFORMATIONS INVOLVING THIS COORD. SYSTEM.
 
 DO  i = 1,3
   ii = (i-1)*3
   ij = i
   IF (ij == 3) ij = 4
   DO  j = 1,3
     jj = j + 1
     ggu(ii+j) = bgpdt(jj,ij)
   END DO
 END DO
!WKBD 11/93 SPR93020      CALL BETRNS (TUB,GGU,0,ELID)
!WKBNB 11/93 SPR93020
!    ADD FROM SHEAR ELEMENT
 
!    COMPUTE DIAGONAL VECTORS
 
 DO  i = 1,3
   ii=i+1
   vd1(i) = bgpdt(ii,3) - bgpdt(ii,1)
   vd2(i) = bgpdt(ii,4) - bgpdt(ii,2)
 END DO
 
!    COMPUTE THE NORMAL VECTOR VKN, NORMALIZE, AND COMPUTE THE PROJECTED
!    AREA, PA
 
 vkn(1) = vd1(2)*vd2(3) - vd1(3)*vd2(2)
 vkn(2) = vd1(3)*vd2(1) - vd1(1)*vd2(3)
 vkn(3) = vd1(1)*vd2(2) - vd1(2)*vd2(1)
 vkl = SQRT( vkn(1)**2 + vkn(2)**2 + vkn(3)**2 )
 IF ( vkl == 0. ) WRITE( nout, 2070 ) est(1)
 2070  FORMAT(//,' ILLEGAL GEOMETRY FOR QUAD4 ELEMENT, ID=',i10 )
 vks(1) = vkn(1)/vkl
 vks(2) = vkn(2)/vkl
 vks(3) = vkn(3)/vkl
 pa = vkl/2.
 
!  COMPUTE SIDES -12- AND -41-
 DO  i = 1,3
   ii = i + 1
   v12(i) = bgpdt(ii,2) - bgpdt(ii,1)
   v41(i) = bgpdt(ii,1) - bgpdt(ii,4)
 END DO
 
!  COMPUTE DOT PRODUCT, V12DK, OR V12 AND VK, THE VECTORS VP12, VI, VJ
 
 v12dk   = v12(1)*vks(1) + v12(2)*vks(2) + v12(3)*vks(3)
 vp12(1) = v12(1) - v12dk*vks(1)
 vp12(2) = v12(2) - v12dk*vks(2)
 vp12(3) = v12(3) - v12dk*vks(3)
 vp12l   = SQRT( vp12(1)**2 + vp12(2)**2 + vp12(3)**2 )
 IF ( vp12l == 0. ) WRITE( nout, 2070 ) est(1)
 vis(1) = vp12(1) / vp12l
 vis(2) = vp12(2) / vp12l
 vis(3) = vp12(3) / vp12l
 vjs(1) = vks(2)*vis(3) - vks(3)*vis(2)
 vjs(2) = vks(3)*vis(1) - vks(1)*vis(3)
 vjs(3) = vks(1)*vis(2) - vks(2)*vis(1)
 
!   NORMALIZE J FOR GOOD MEASURE
 
 vjl = SQRT( vjs(1)**2 + vjs(2)**2 + vjs(3)**2 )
 IF ( vjl == 0. ) WRITE ( nout, 2070 ) est(1)
 vjs(1) = vjs(1) / vjl
 vjs(2) = vjs(2) / vjl
 vjs(3) = vjs(3) / vjl
 DO  i = 1,3
   tub(i)   = vis(i)
   tub(i+3) = vjs(i)
   tub(i+6) = vks(i)
 END DO
!WKBNE 11/93 SPR93020
 
 
!     STORE INCOMING BGPDT FOR ELEMENT C.S.
 
 DO  i = 1,3
   i1 = i + 1
   DO  j = 1,4
     bgpdm(i,j) = bgpdt(i1,j)
   END DO
 END DO
 
!     TRANSFORM BGPDM FROM BASIC TO USER C.S.
 
 DO  i = 1,3
   ip = (i-1)*3
   DO  j = 1,4
     ugpdm(i,j) = 0.0
     DO  k = 1,3
       kk = ip + k
       ugpdm(i,j) = ugpdm(i,j) + tub(kk)*((bgpdm(k,j))-ggu(k))
     END DO
   END DO
 END DO
 
!     THE ORIGIN OF THE ELEMENT C.S. IS IN THE MIDDLE OF THE ELEMENT
 
 DO  j = 1,3
   cent(j) = 0.0
   DO  i = 1,4
     cent(j) = cent(j) + ugpdm(j,i)/nnode
   END DO
 END DO
 
!     STORE THE CORNER NODE DIFF. IN THE USER C.S.
 
 x31 = ugpdm(1,3) - ugpdm(1,1)
 y31 = ugpdm(2,3) - ugpdm(2,1)
 x42 = ugpdm(1,4) - ugpdm(1,2)
 y42 = ugpdm(2,4) - ugpdm(2,2)
 aa  = SQRT(x31*x31 + y31*y31)
 bb  = SQRT(x42*x42 + y42*y42)
 
!     NORMALIZE XIJ'S
 
 x31 = x31/aa
 y31 = y31/aa
 x42 = x42/bb
 y42 = y42/bb
 exi = x31 - x42
 exj = y31 - y42
 
!     STORE GGE ARRAY, THE OFFSET BETWEEN ELEMENT C.S. AND USER C.S.
 
 gge(1) = cent(1)
 gge(2) = cent(2)
 gge(3) = cent(3)
 
 gge(4) = gge(1) + exi
 gge(5) = gge(2) + exj
 gge(6) = gge(3)
 
 gge(7) = gge(1) - exj
 gge(8) = gge(2) + exi
 gge(9) = gge(3)
 
!     START FILLING IN IELOUT ARRAY WITH DATA TO BE STORED IN GPSRN
 
 ielout(1) = elid
 DO  i = 1,4
   ielout(3+(i-1)*25) = sil(i)
   DO  j = 1,3
     relout(25*i+j-1) = bgpdt(j+1,i)
   END DO
 END DO
 
!     THE ARRAY IORDER STORES THE ELEMENT NODE ID IN
!     INCREASING SIL ORDER.
 
!     IORDER(1) = NODE WITH LOWEST  SIL NUMBER
!     IORDER(4) = NODE WITH HIGHEST SIL NUMBER
 
!     ELEMENT NODE NUMBER IS THE INTEGER FROM THE NODE LIST G1,G2,G3,G4.
!     THAT IS, THE 'I' PART OF THE 'GI' AS THEY ARE LISTED ON THE
!     CONNECTIVITY BULK DATA CARD DESCRIPTION.
 
 
 DO  i = 1,4
   iorder(i) = 0
   ksil(i) = sil(i)
 END DO
 
 DO  i = 1,4
   itemp = 1
   isil  = ksil(1)
   DO  j = 2,4
     IF (isil <= ksil(j)) CYCLE
     itemp = j
     isil  = ksil(j)
   END DO
   iorder(i)   = itemp
   ksil(itemp) = 99999999
 END DO
 
!     ADJUST EST DATA
 
!     USE THE POINTERS IN IORDER TO COMPLETELY REORDER THE
!     GEOMETRY DATA INTO INCREASING SIL ORDER.
!     DON'T WORRY!! IORDER ALSO KEEPS TRACK OF WHICH SHAPE
!     FUNCTIONS GO WITH WHICH GEOMETRIC PARAMETERS!
 
 DO  i = 1,4
   ksil(i)   = sil(i)
   tmpthk(i) = gpth(i)
   kcid(i)   = igpdt(1,i)
   DO  j = 2,4
     tgrid(j,i) = bgpdt(j,i)
   END DO
 END DO
 DO  i = 1,4
   ipoint  = iorder(i)
   gpth(i) = tmpthk(ipoint)
   igpdt(1,i) = kcid(ipoint)
   sil(i)     = ksil(ipoint)
   nphi(i+1 ) = ksil(ipoint)
   nphi(i+5 ) = 0
   nphi(i+9 ) = ipoint
   nphi(i+13) = 0
   DO  j = 2,4
     bgpdt(j,i) = tgrid(j,ipoint)
   END DO
 END DO
 
 nphi(19) = nest(20)
 nphi(20) = nest(21)
 phiout(18) = 0.0
 offset   = zoff
 IF (zoff == 0.0) offset = zoff1
 phiout(78) = offset
 
!     COMPUTE NODE NORMALS
 
 CALL q4nrms (bgpdt,gpnorm,iorder,iflag)
 IF (iflag == 0) GO TO 130
 WRITE (nout,1710) ufm,elid
 GO TO 1430
 130 CONTINUE
 
!     PUT NORMALS IN IELOUT
 
 DO  i = 1,nnode
   io  = iorder(i)
   iop = (io-1)*25 + 21
   relout(iop+1) = gpnorm(2,i)
   relout(iop+2) = gpnorm(3,i)
   relout(iop+3) = gpnorm(4,i)
 END DO
 
!     COMPUTE NODE NORMALS
 
 avgthk = 0.0
 DO  i = 1,nnode
   io = iorder(i)
   IF (gpth(i) == 0.0) gpth(i) = elth
   IF (gpth(i) > 0.0) GO TO 150
   WRITE (nout,1700) ufm,elid,sil(i)
   GO TO 1430
   150 avgthk = avgthk + gpth(i)/nnode
   gpth2(io) = gpth(i)
 END DO
 
 mominr = 0.0
 tsfact = 5.0/6.0
 nocsub = .false.
 IF (nest(15) /=  0) mominr = est(16)
 IF (nest(17) /=  0) ts = est(18)
 IF ( est(18) == .0) ts = 5.0/6.0
 phiout(21) = avgthk
 phiout(22) = mominr
 
!     SET LOGICAL NOCSUB IF EITHER MOMINR OR TS ARE NOT DEFAULT
!     VALUES. THIS WILL BE USED TO OVERRIDE ALL CSUBB COMPUTATIONS.
!     I.E. DEFAULT VALUES OF UNITY ARE USED.
 
 epsi = ABS(mominr - 1.0)
 epst = ABS(ts  - tsfact)
 eps  = .05
!     NOCSUB = EPSI.GT.EPS .OR. EPST.GT.EPS
 
!     PUT THE AVERAGE THICKNESS IN RELOUT
 
 relout(2) = avgthk
 
!     THE COORDINATES OF THE ELEMENT GRID POINTS HAVE TO BE
!     TRANSFORMED FROM THE BASIC C.S. TO THE ELEMENT C.S.
 
 CALL betrns (teu,gge,0,elid)
 CALL gmmats (teu,3,3,0, tub,3,3,0, teb)
 CALL gmmats (tub,3,3,1, cent,3,1,0, cente)
 
 DO  i = 1,3
   ii = i + 1
   ip = (i-1)*3
   DO  j = 1,nnode
     epnorm(ii,j) = 0.0
     egpdt (ii,j) = 0.0
     DO  k = 1,3
       kk = ip + k
       k1 = k + 1
       cc = bgpdt(k1,j) - ggu(k) - cente(k)
       epnorm(ii,j) = epnorm(ii,j) + teb(kk)*gpnorm(k1,j)
       egpdt (ii,j) = egpdt (ii,j) + teb(kk)*cc
     END DO
   END DO
 END DO
 
!     INITIALIZE MATERIAL VARIABLES
 
!     SET INFLAG = 12 SO THAT SUBROUTINE MAT WILL SEARCH FOR-
!     ISOTROPIC MATERIAL PROPERTIES AMONG THE MAT1 CARDS,
!     ORTHOTROPIC MATERIAL PROPERTIES AMONG THE MAT8 CARDS, AND
!     ANISOTROPIC MATERIAL PROPERTIES AMONG THE MAT2 CARDS.
 
 inflag = 12
 rho    = 0.0
 eltemp =  est(45)
 mid(1) = nest(13)
 mid(2) = nest(15)
 mid(3) = nest(17)
 mid(4) = nest(22)
 membrn = mid(1) > 0
 bendng = mid(2) > 0 .AND. mominr > 0.0
 shrflx = mid(3) > 0
 mbcoup = mid(4) > 0
 
!     CHECK FOR COMPOSITE MATERIAL
 
 nphi(79) = 0
 DO  img = 1,4
   IF (mid(img) > hunmeg) GO TO 190
 END DO
 GO TO 200
 190 nphi(79) = mid(img) - img*hunmeg
 200 CONTINUE
 
!     DETERMINE FACTORS TO BE USED IN CSUBB CALCULATIONS
 
 IF (.NOT.bendng) GO TO 250
 DO  i = 1,4
   DO  j = 1,nnode
     jo = iorder(j)
     IF (i /= jo) CYCLE
     xa(i) = egpdt(2,j)
     yb(i) = egpdt(3,j)
     zc(i) = egpdt(4,j)
     vnt(1,i) = epnorm(2,j)
     vnt(2,i) = epnorm(3,j)
     vnt(3,i) = epnorm(4,j)
   END DO
 END DO
 
 a = 0.5*(xa(2) + xa(3) - xa(1) - xa(4))
 b = 0.5*(yb(4) + yb(3) - yb(1) - yb(2))
 IF (a > b) aspect = b/a
 IF (a <= b) aspect = a/b
 
!     IRREGULAR 4-NODE CODE-  GEOMETRIC VARIABLES
 
!     CALCULATE AND NORMALIZE- UNIT EDGE VECTORS,UNIT NORMAL VECTORS
 
 DO  i = 1,4
   j = i + 1
   IF (j == 5) j = 1
   uev(1,i) = xa(j) - xa(i)
   uev(2,i) = yb(j) - yb(i)
   uev(3,i) = zc(j) - zc(i)
   unv(1,i) = (vnt(1,j)+vnt(1,i))*0.50
   unv(2,i) = (vnt(2,j)+vnt(2,i))*0.50
   unv(3,i) = (vnt(3,j)+vnt(3,i))*0.50
   cc       = uev(1,i)**2 + uev(2,i)**2 + uev(3,i)**2
   IF (cc >= 1.0E-8) cc = SQRT(cc)
   edgel(i) = cc
   uev(1,i) = uev(1,i)/cc
   uev(2,i) = uev(2,i)/cc
   uev(3,i) = uev(3,i)/cc
   cc       = SQRT(unv(1,i)**2 + unv(2,i)**2 + unv(3,i)**2)
   unv(1,i) = unv(1,i)/cc
   unv(2,i) = unv(2,i)/cc
   unv(3,i) = unv(3,i)/cc
 END DO
 
!     CALCULATE INTERNAL NODAL ANGLES
 
 DO  i = 1,4
   j = i - 1
   IF (j == 0) j = 4
   anglei(i) =-uev(1,i)*uev(1,j)-uev(2,i)*uev(2,j)-uev(3,i)*uev(3,j)
   IF (ABS(anglei(i)) < 1.0E-8) anglei(i) = 0.0
 END DO
 250 CONTINUE
 
!     SET THE INTEGRATION POINTS
 
 ptint(1) = -const
 ptint(2) =  const
 
 IF (itherm /= 0) GO TO 1500
 
!     IN PLANE SHEAR REDUCTION
 
 xi  = 0.0
 eta = 0.0
 kpt = 1
 
 CALL q4shps (xi,eta,shp,dshp)
 
!     SORT THE SHAPE FUNCTIONS AND THEIR DERIVATIVES INTO SIL ORDER.
 
 DO  i = 1,4
   tmpshp(i  ) = shp (i  )
   dshptp(i  ) = dshp(i  )
   dshptp(i+4) = dshp(i+4)
 END DO
 DO  i = 1,4
   kk = iorder(i)
   shp( i  ) = tmpshp(kk  )
   dshp(i  ) = dshptp(kk  )
   dshp(i+4) = dshptp(kk+4)
 END DO
 
 DO  izta = 1,2
   zeta = ptint(izta)
   
!     COMPUTE THE JACOBIAN AT THIS GAUSS POINT,
!     ITS INVERSE AND ITS DETERMINANT.
   
   hzta = zeta/2.0
   
   CALL jacobs (elid,shp,dshp,gpth,egpdt,epnorm,jacob)
   IF (badjac) GO TO 1430
   
!     COMPUTE PSI TRANSPOSE X JACOBIAN INVERSE.
!     HERE IS THE PLACE WHERE THE INVERSE JACOBIAN IS FLAGED TO BE
!     TRANSPOSED BECAUSE OF OPPOSITE MATRIX LOADING CONVENTION BETWEEN
!     INVER AND GMMAT.
   
   CALL gmmats (psitrn,3,3,0, jacob,3,3,1, phi)
   
!     CALL Q4BMGS TO GET B MATRIX
!     SET THE ROW FLAG TO 2. IT WILL SAVE THE 3RD ROW OF B-MATRIX AT
!     THE TWO INTEGRATION POINTS.
   
   rowflg = 2
   CALL q4bmgs (dshp,gpth,egpdt,epnorm,phi,xybmat(kpt))
   kpt = kpt + nd2
 END DO
 
!     FETCH MATERIAL PROPERTIES
 
!     SET THE ARRAY OF LENGTH 4 TO BE USED IN CALLING TRANSS.
!     NOTE THAT THE FIRST WORD IS THE COORDINATE SYSTEM ID WHICH
!     WILL BE SET IN POSITION LATER.
 
 290 DO  iec = 2,4
   ecpt(iec) = 0.0
 END DO
 
 
!     EACH MATERIAL PROPERTY MATRIX G HAS TO BE TRANSFORMED FROM
!     THE MATERIAL COORDINATE SYSTEM TO THE ELEMENT COORDINATE
!     SYSTEM. THESE STEPS ARE TO BE FOLLOWED-
 
!     1- IF MCSID HAS BEEN SPECIFIED, SUBROUTINE TRANSS IS CALLED
!        TO CALCULATE TBM-MATRIX (MATERIAL TO BASIC TRANSFORMATION).
!        THIS WILL BE FOLLOWED BY A CALL TO SUBROUTINE BETRNS
!        TO CALCULATE TEB-MATRIX (BASIC TO ELEMENT TRANSFORMATION).
!        TBM-MATRIX IS THEN PREMULTIPLIED BY TEB-MATRIX TO OBTAIN
!        TEM-MATRIX. THEN STEP 3 WILL BE TAKEN.
 
!     2- IF THETAM HAS BEEN SPECIFIED, SUBROUTINE ANGTRS IS CALLED
!        TO CALCULATE TEM-MATRIX (MATERIAL TO ELEMENT TRANSFORMATION).
 
!                          T
!     3-           G   =  U   G   U
!                   E          M
 
 
 flagm = nest(11)
 IF (flagm == 0) GO TO 360
 mcsid = nest(10)
 
!     CALCULATE TUM-MATRIX USING MCSID
 
 310 IF (mcsid > 0) GO TO 330
 DO  i = 1,9
   tem(i) = teb(i)
 END DO
 GO TO 340
 330 necpt(1) = mcsid
 CALL transs (ecpt,tbm)
 
!     MULTIPLY TEB AND TBM MATRICES
 
 CALL gmmats (teb,3,3,0, tbm,3,3,0, tem)
 
!     CALCULATE THETAM FROM THE PROJECTION OF THE X-AXIS OF THE
!     MATERIAL C.S. ON TO THE XY PLANE OF THE ELEMENT C.S.
 
 340 CONTINUE
 xm = tem(1)
 ym = tem(4)
 IF (ABS(xm) > eps1 .OR. ABS(ym) > eps1) GO TO 350
 nest(2) = mcsid
 j = 231
 GO TO 1440
 350 thetam = ATAN2(ym,xm)
 GO TO 370
 
!     CALCULATE TEM-MATRIX USING THETAM
 
 360 thetam = est(10)*degrad
 IF (thetam == 0.0) GO TO 380
 370 CALL angtrs (thetam,1,tum)
 CALL gmmats (teu,3,3,0, tum,3,3,0, tem)
 GO TO 400
 
!     DEFAULT IS CHOSEN, LOOK FOR VALUES OF MCSID AND/OR THETAM
!     ON THE PSHELL CARD.
 
 380 flagm = nest(24)
 IF (flagm == 0) GO TO 390
 mcsid = nest(23)
 GO TO 310
 
 390 thetam = est(23)*degrad
 GO TO 370
 
 400 CONTINUE
 
!     STORE TUM IN PHIOUT
 
 DO  iem = 1,9
   phiout(68+iem) = tum(iem)
 END DO
 
 IF (itherm /= 0) GO TO 1600
 
!     BEGIN THE LOOP TO FETCH PROPERTIES FOR EACH MATERIAL ID
 
 DO  ll = 1,36
   gi(ll) = 0.0
 END DO
 
 m    = 0
 it0  = 0
 igobk= 0
 430 m    = m + 1
 IF (m > 4) GO TO 680
 IF (m == 4 .AND. igobk == 1) GO TO 690
 matid = mid(m)
 IF (matid == 0 .AND. m /= 3) GO TO 430
 IF (matid == 0 .AND. m == 3 .AND. .NOT.bendng) GO TO 430
 IF (matid == 0 .AND. m == 3 .AND. bendng) matid = mid(2)
 
 IF (m-1 < 0) THEN
   GO TO   460
 ELSE IF (m-1 == 0) THEN
   GO TO   450
 END IF
 440 IF (matid == mid(m-1) .AND. igobk == 0) GO TO 460
 450 CALL mat (elid)
 460 CONTINUE
 
 IF (it0 > 0) GO TO 470
 tsub0 = rmtout(11)
 IF (matset == 8.0) tsub0 = rmtout(10)
 phiout(18) = tsub0
 it0 = 1
 470 CONTINUE
 
 coeff = 1.0
!     IF (M .EQ. 2) COEFF = MOMINR
 IF (m == 3) coeff = ts
 lpoint = (m-1)*9 + 1
 
 CALL q4gmgs (m,coeff,gi(lpoint))
 
!WKBDB 11/93 SPR93020
!      IF (M .GT. 0) GO TO 490
!      IF (.NOT.SHRFLX .AND. BENDNG) GO TO 480
!      NEST(2) = MATID
!      J = 231
!      GO TO 1440
 
!  480 M = -M
!  490 CONTINUE
!      MTYPE = IFIX(MATSET+.05) - 2
!      IF (NOCSUB) GO TO 580
!      GO TO (580,500,540,580), M
!C
!  500 IF (MTYPE) 510,520,530
!  510 ENORX = RMTOUT(16)
!      ENORY = RMTOUT(16)
!      GO TO 580
!  520 ENORX = RMTOUT(1)
!      ENORY = RMTOUT(4)
!      GO TO 580
!  530 ENORX = RMTOUT(1)
!      ENORY = RMTOUT(3)
!      GO TO 580
 
!  540 IF (MTYPE) 550,560,570
!  550 GNORX = RMTOUT(6)
!      GNORY = RMTOUT(6)
!      GO TO 580
!  560 GNORX = RMTOUT(1)
!      GNORY = RMTOUT(4)
!      GO TO 580
!  570 GNORX = RMTOUT(6)
!      GNORY = RMTOUT(5)
!      IF (GNORX .EQ. 0.0) GNORX = RMTOUT(4)
!      IF (GNORY .EQ. 0.0) GNORY = RMTOUT(4)
!  580 CONTINUE
!WKBDE 11/93 SPR93020
!WKBNB 11/93 SPR93020
 IF (m > 0) GO TO 490
 IF (.NOT.shrflx .AND. bendng) GO TO 480
 nest(2) = matid
 j = 231
 GO TO 1440
 480 m = -m
 490 CONTINUE
 mtype = IFIX(matset+.05) - 2
 IF (nocsub) GO TO 580
 SELECT CASE ( m )
   CASE (    1)
     GO TO 580
   CASE (    2)
     GO TO 500
   CASE (    3)
     GO TO 540
   CASE (    4)
     GO TO 580
 END SELECT
!WKBNE 11/93 SPR93020
!WKBNB 2/94 SPR93020
 500 IF ( mtype  < 0) THEN
   GO TO   510
 ELSE IF ( mtype  == 0) THEN
   GO TO   520
 ELSE
   GO TO   530
 END IF
 510 enorx = rmtout(16)
 enory = rmtout(16)
 dnux  = gi( lpoint+1 ) / gi( lpoint )
 dnuy  = gi( lpoint+3 ) / gi( lpoint+4 )
 GO TO 580
 520 enorx = rmtout(1)
 enory = rmtout(4)
 dnux  = gi( lpoint+1 ) / gi( lpoint )
 dnuy  = gi( lpoint+3 ) / gi( lpoint+4 )
 GO TO 580
 530 enorx = rmtout(1)
 enory = rmtout(3)
 dnux  = gi( lpoint+1 ) / gi( lpoint )
 dnuy  = gi( lpoint+3 ) / gi( lpoint+4 )
 GO TO 580
 540 IF ( mtype  < 0) THEN
   GO TO   550
 ELSE IF ( mtype  == 0) THEN
   GO TO   560
 ELSE
   GO TO   570
 END IF
 550 gnorx = rmtout(6)
 gnory = rmtout(6)
 GO TO 580
 560 gnorx = rmtout(1)
 gnory = rmtout(4)
 GO TO 580
 570 gnorx = rmtout(6)
 gnory = rmtout(5)
 IF ( gnorx == 0.0D0 ) gnorx = rmtout(4)
 IF ( gnory == 0.0D0 ) gnory = rmtout(4)
 580 CONTINUE
!WKBNE 2/94 SPR93020
 IF (matset == 1.0) GO TO 610
 IF (m      ==   3) GO TO 590
 u(1) = tem(1)*tem(1)
 u(2) = tem(2)*tem(2)
 u(3) = tem(1)*tem(2)
 u(4) = tem(4)*tem(4)
 u(5) = tem(5)*tem(5)
 u(6) = tem(4)*tem(5)
 u(7) = tem(1)*tem(4)*2.0
 u(8) = tem(2)*tem(5)*2.0
 u(9) = tem(1)*tem(5) + tem(2)*tem(4)
 l    = 3
 GO TO 600
 
 590 u(1) = tem(5)*tem(9) + tem(6)*tem(8)
 u(2) = tem(4)*tem(9) + tem(6)*tem(7)
 u(3) = tem(2)*tem(9) + tem(3)*tem(8)
 u(4) = tem(1)*tem(9) + tem(3)*tem(7)
 l    = 2
 
 600 CALL gmmats (u(1),l,l,1, gi(lpoint),l,l,0, GT(1))
 CALL gmmats (GT(1),l,l,0, u(1),l,l,0, gi(lpoint))
 
!     TRANSFORM THERMAL EXPANSION COEFF'S AND STORE THEM IN PHIOUT
 
 610 CONTINUE
 IF (m      > 2 ) GO TO 430
 IF (matset == 2.) GO TO 620
 IF (matset == 8.) GO TO 640
 
!     MAT1
 
 alfa(1) = rmtout(8)
 alfa(2) = rmtout(8)
 alfa(3) = 0.0
 GO TO 650
 
!     MAT2
 
 620 DO  imat = 1,3
   alfa(imat) = rmtout(7+imat)
 END DO
 GO TO 650
 
!     MAT8
 
 640 alfa(1) = rmtout(8)
 alfa(2) = rmtout(9)
 alfa(3) = 0.0
 
 650 mpoint = (m-1)*3 + 59
 IF (matset == 1.0) GO TO 660
 CALL invers (3,u,3,bdum,0,detu,isngu,INDEX)
 CALL gmmats (u,3,3,0, alfa,3,1,0, phiout(mpoint))
 GO TO 430
 660 DO  ialf = 1,3
   mp = mpoint - 1 + ialf
   phiout(mp) = alfa(ialf)
 END DO
 GO TO 430
 680 CONTINUE
 IF (mid(3) < hunmeg) GO TO 690
 IF (gi(19) /= 0. .OR. gi(20) /= 0. .OR. gi(21) /= 0. .OR.  &
     gi(22) /= 0.) GO TO 690
 igobk = 1
 m = 2
 mid(3) = mid(2)
 GO TO 430
 690 CONTINUE
 
 nocsub = enorx == 0.0 .OR. enory == 0.0 .OR.  &
     gnorx == 0.0 .OR. gnory == 0.0 .OR. mominr == 0.0
 
 
!     FILL IN THE BASIC 6X6 MATERIAL PROPERTY MATRIX G
 
 DO  ig = 1,6
   DO  jg = 1,6
     g(ig,jg) = 0.0
   END DO
 END DO
 
 IF (.NOT.membrn) GO TO 720
 DO  ig = 1,3
   ig1 = (ig-1)*3
   DO  jg = 1,3
     jg1 = jg + ig1
     g(ig,jg) = gi(jg1)
   END DO
 END DO
 
 720 IF (.NOT.bendng) GO TO 750
 DO  ig = 4,6
   ig2 = (ig-2)*3
   DO  jg = 4,6
     jg2 = jg + ig2
     g(ig,jg) = gi(jg2)
   END DO
 END DO
 
 IF (.NOT.membrn) GO TO 750
 DO  ig = 1,3
   kg  = ig + 3
   ig1 = (ig-1)*3
   DO  jg = 1,3
     lg  = jg + 3
     jg1 = jg + ig1
     g(ig,lg) = gi(jg1)
     g(kg,jg) = gi(jg1)
   END DO
 END DO
 
!     STORE 6X6 GBAR-MATRIX IN PHIOUT
 
 750 ig1 = 22
 DO  ig = 1,6
   DO  jg = 1,6
     ig1 = ig1 + 1
     phiout(ig1) = g(ig,jg)
   END DO
 END DO
 
 
!     STRESS TRANSFORMATIONS
!     ----------------------
 
!     THE NECESSARY TRANSFORMATIONS ARE PERFORMED IN THE FOLLOWING
!     MANNER-
 
!     1- ALL THE TRANSFORMATIONS ARE CALCULATED IN PHASE I AND THEN
!        TRANSFERED THRU DATA BLOCK 'PHIOUT' TO PHASE II WHERE THE
!        ACTUAL MULTIPLICATIONS ARE PERFORMED.
 
!     2- THE STRAIN RECOVERY MATRIX B
!        IS EVALUATED IN THE ELEMENT COORDINATE SYSTEM IN PHASE I
!        AND TRANSFERED TO PHASE II. THE DISPLACEMENTS, HOWEVER,
!        ENTER PHASE II IN GLOBAL COORDINATES. THEREFORE,
!        2A) 3X3 TRANSFORMATIONS FROM GLOBAL TO ELEMENT COORDINATE
!            SYSTEM (TEG) FOR EACH GRID POINT ARE CALCULATED AND
!            STORED IN  PHIOUT (80 - (79+9*NNODE)).
!            USING THESE TRANSFORMATIONS THE DISPLACEMENTS AT
!            EACH GRID POINT WILL BE EVALUATED IN THE ELEMENT
!            COORDINATE SYSTEM AFTER ENTERING PHASE II.
 
!        2B) A 3X3 TRANSFORMATION FROM THE TANGENT TO THE USER-
!            DEFINED STRESS COORDINATE SYSTEM (TSI) IS CALCULATED
!            FOR EACH INTEGRATION POINT AND STORED ALONG WITH OTHER
!            DATA FOR THAT INTEGRATION POINT AT POSITIONS 2-10 OF
!            THE REPEATED DATA FOR EACH EVALUATION POINT.
!            IT WILL BE USED TO TRANSFORM THE STRESS OUTPUT TO
!            ANY DESIRED COORDINATE SYSTEM.
!            NOTE THAT THESE CALCULATIONS WILL BE PERFORMED INSIDE
!            THE DOUBLE LOOP.
 
!     CALCULATIONS FOR TEG-MATRIX
 
!     CALCULATE  TBG-MATRIX (GLOBAL TO BASIC), THEN
!     MULTIPLY  TEB AND TBG MATRICES  TO GET  TEG-MATRIX
!     FOR THIS GRID POINT AND STORE IT IN PHIOUT.
 
 DO  i = 1,nnode
   ip = 80 + (i-1)*9
   IF (igpdt(1,i) <= 0) GO TO 800
   CALL transs (igpdt(1,i),tbg)
   CALL gmmats (teb,3,3,0, tbg,3,3,0, phiout(ip))
   CYCLE
   
   800 DO  j = 1,9
     phiout(ip+j-1) = teb(j)
   END DO
 END DO
 
!     INITIALIZE THE ARRAYS USED IN THE DOUBLE LOOP CALCULATION.
!     EVALUATION OF STRESSES IS DONE AT 2X2 POINTS AND AT THE
!     CENTER OF THE ELEMENT, AT THE MID-SURFACE.
 
 IF (bendng) GO TO 840
 j = nd3 + 1
 DO  ibmx = j,nd8
   bmatrx(ibmx) = 0.0
 END DO
 840 CONTINUE
 
 icount = -(8*ndof+nnode+32) + 79 + 9*nnode
 
 ptintp(1) =-const
 ptintp(2) = const
 ptintp(3) = 0.0
 
 
!     HERE BEGINS THE TRIPLE LOOP ON STATEMENTS 835 AND 840
!     -----------------------------------------------------
 
 DO  ixsi = 1,3
   xi = ptintp(ixsi)
   
   DO  ieta = 1,3
     eta = ptintp(ieta)
     IF (ixsi == 3 .AND. ieta /= 3) CYCLE
     IF (ixsi /= 3 .AND. ieta == 3) CYCLE
     
     CALL q4shps (xi,eta,shp,dshp)
     
!     SORT THE SHAPE FUNCTIONS AND THEIR DERIVATIVES INTO SIL ORDER.
     
     DO  i = 1,4
       tmpshp(i  ) = shp (i  )
       dshptp(i  ) = dshp(i  )
       dshptp(i+4) = dshp(i+4)
     END DO
     DO  i = 1,4
       kk = iorder(i)
       shp (i  ) = tmpshp(kk  )
       dshp(i  ) = dshptp(kk  )
       dshp(i+4) = dshptp(kk+4)
     END DO
     
     th = 0.0
     DO  ith = 1,nnode
       th = th + shp(ith)*gpth(ith)
     END DO
     reali = mominr*th*th*th/12.0
     tsi = ts*th
     
     IF (nocsub) GO TO 970
     IF (.NOT.bendng) GO TO 970
!      NUNORX = MOMINR*ENORX/(2.0*GNORX) - 1.0
!      NUNORY = MOMINR*ENORY/(2.0*GNORY) - 1.0
!WKBNB 2/94 SPR93020
     nunorx = mominr*enorx/(2.0*gnorx) - 1.0
     nunory = mominr*enory/(2.0*gnory) - 1.0
     IF ( nunorx < 0. ) nunorx = dnux
     IF ( nunory < 0. ) nunory = dnuy
!WKBNE 2/94 SPR93020
!WKBDB 2/94 SPR93020
!      EIX = MOMINR*ENORX
!      EIY = MOMINR*ENORY
!      TGX = 2.0*GNORX
!      TGY = 2.0*GNORY
!      NUNORX = EIX/TGX - 1.0
!      IF (EIX .GT. TGX) NUNORX = 1.0 - TGX/EIX
!      NUNORY = EIY/TGY - 1.0
!      IF (EIY .GT. TGY) NUNORY = 1.0 - TGY/EIY
!      IF (NUNORX .GT. 0.999999) NUNORX = 0.999999
!      IF (NUNORY .GT. 0.999999) NUNORY = 0.999999
!WKBDE 2/94 SPR93020
!     IF (NUNORX .GT. .49) NUNORX = 0.49
!     IF (NUNORY .GT. .49) NUNORY = 0.49
     cc = aspect
     ax = a
     IF (eta < 0.0) ax = a + const*(xa(2)-xa(1)-a)
     IF (eta > 0.0) ax = a + const*(xa(3)-xa(4)-a)
     psiinx = 32.0*reali/((1.0-nunorx)*tsi*ax*ax)
     by = b
     IF (xi < 0.0) by = b + const*(yb(4)-yb(1)-b)
     IF (xi > 0.0) by = b + const*(yb(3)-yb(2)-b)
     psiiny = 32.0*reali/((1.0-nunory)*tsi*by*by)
     IF (.NOT.shrflx) GO TO 930
     tsmfx = psiinx
     tsmfy = psiiny
     IF (tsmfx > 1.0) tsmfx = 1.0
     IF (tsmfy > 1.0) tsmfy = 1.0
     GO TO 980
     930 IF (psiinx >= 1.0) GO TO 940
     tsmfx = psiinx/(1.0-psiinx)
     IF (tsmfx <= 1.0) GO TO 950
     940 tsmfx = 1.0
     950 IF (psiiny >= 1.0) GO TO 960
     tsmfy = psiiny/(1.0-psiiny)
     IF (tsmfy <= 1.0) GO TO 980
     960 tsmfy = 1.0
     GO TO 980
     
     970 tsmfx = 1.0
     tsmfy = 1.0
     980 CONTINUE
     
!     IRREGULAR 4-NODE CODE-  CALCULATION OF NODAL EDGE SHEARS
!                             AT THIS INTEGRATION POINT
     
     
     DO  ij = 1,4
       ii = ij - 1
       IF (ii == 0) ii = 4
       ik = ij + 1
       IF (ik == 5) ik = 1
       
       DO  ir = 1,4
         IF (ij /= iorder(ir)) CYCLE
         ioj = ir
         EXIT
       END DO
       1010 DO  ir = 1,4
         IF (ik /= iorder(ir)) CYCLE
         iok = ir
         EXIT
       END DO
       1030 aa = shp(ioj)
       bb = shp(iok)
       
       DO  is = 1,3
         edgshr(is,ij) = (uev(is,ij)+anglei(ij)*uev(is,ii))*aa/  &
             (1.0-anglei(ij)*anglei(ij))  &
             + (uev(is,ij)+anglei(ik)*uev(is,ik))*bb/ (1.0-anglei(ik)*anglei(ik))
       END DO
     END DO
     
     DO  izta = 1,2
       zeta = ptint(izta)
       hzta = zeta/2.0
       ibot = (izta-1)*nd2
       
!     SET THE PHIOUT POINTER
       
       icount = icount + 32 + nnode + 8*ndof
       
       phiout(icount+1) = th
       
!     STORE SHAPE FUNCTION VALUES IN PHIOUT
       
       DO  i = 1,nnode
         phiout(icount+32+i) = shp(i)
       END DO
       
!     STORE THE CORRECTION TO GBAR-MATRIX IN PHIOUT
       
       ig1 = icount + 10
       ig4 = 28
       DO  ig = 1,9
         ig1 = ig1 + 1
         phiout(ig1) = -gi(ig4)*zeta*6.0
         ig4 = ig4 + 1
       END DO
       
!     STORE G3-MATRIX IN PHIOUT
       
       iph = icount + 28
       phiout(iph+1) = tsmfy*gi(19)
       phiout(iph+2) = SQRT(tsmfx*tsmfy)*gi(20)
       phiout(iph+3) = SQRT(tsmfx*tsmfy)*gi(21)
       phiout(iph+4) = tsmfx*gi(22)
       
!     COMPUTE THE JACOBIAN AT THIS GAUSS POINT,
!     ITS INVERSE AND ITS DETERMINANT.
       
       CALL jacobs (elid,shp,dshp,gpth,egpdt,epnorm,jacob)
       IF (badjac) GO TO 1430
       
!     COMPUTE PSI TRANSPOSE X JACOBIAN INVERSE.
!     HERE IS THE PLACE WHERE THE INVERSE JACOBIAN IS FLAGED TO BE
!     TRANSPOSED BECAUSE OF OPPOSITE MATRIX LOADING CONVENTION BETWEEN
!     INVER AND GMMAT.
       
       CALL gmmats (psitrn,3,3,0, jacob,3,3,1, phi)
       
       CALL gmmats (tem,3,3,1, psitrn,3,3,1, tmi)
       
!     STORE TMI-MATRIX IN PHIOUT
       
       iph = icount + 20
       DO  i = 1,9
         phiout(iph) = tmi(i)
         iph = iph + 1
       END DO
       
!     ARRAY ECPT(4) WHICH IS USED IN TRANSS CONSISTS OF THE C.S. ID
!     AND THE COORDINATES (IN BASIC C.S.) OF THE POINT FROM (OR TO)
!     WHICH THE TRANSFORMATION IS BEING PERFORMED. THE COORDINATES
!     ARE NOT USED IF THE DESIGNATED COORDINATE SYSTEM IS RECTANGULAR.
       
       DO  i = 1,3
         gpc(i) = 0.0
         ii = i + 1
         DO  j = 1,nnode
           gpc(i) = gpc(i) + shp(j)*(bgpdt(ii,j) + hzta*gpth(j)*gpnorm(ii,j))
         END DO
         ecpt(ii) = gpc(i)
       END DO
       
!     CALCULATIONS FOR TSE-MATRIX
       
       flags = nest(27)
       IF (flags == 0) GO TO 1300
       
!     FLAGS IS 1, I.E. SCSID HAS BEEN SPECIFIED.
!     CALCULATE TBS-MATRIX (STRESS TO BASIC)
       
       scsid = nest(26)
       IF (scsid <= 0) GO TO 1200
       necpt(1) = scsid
       CALL transs (ecpt,tbs)
       GO TO 1220
       1200 DO  i = 1,3
         ii = (i-1)*3
         DO  j = 1,3
           jj = (j-1)*3
           tsu(ii+j) = tub(i+jj)
         END DO
       END DO
       GO TO 1230
       
!     MULTIPLY
!               T         T
!            TBS  AND  TUB  TO GET TSU-MATRIX (USER TO STRESS)
       
       1220 CALL gmmats (tbs,3,3,1, tub,3,3,1, tsu)
       
!     CALCULATE THETAS FROM THE PROJECTION OF THE X-AXIS OF THE
!     STRESS C.S. ON TO THE XY PLANE OF THE ELEMENT C.S.
       
       1230 CONTINUE
       xs = tsu(1)
       ys = tsu(2)
       IF (ABS(xs) > eps1 .OR. ABS(ys) > eps1) GO TO 1240
       nest(2) = scsid
       j = 233
       GO TO 1440
       1240 thetas = ATAN2(ys,xs)
       GO TO 1310
       
!     FLAGS IS 0, I.E. THETAS HAS BEEN SPECIFIED.
!     SUBROUTINE ANGTRS RETURNS THE 3X3 TRANSFORMATION USING THETAS.
!     NOTE THAT IF THETAS IS LEFT BLANK (DEFAULT), THE TRANSFORMATION
!     WILL BE IDENTITY,  I.E. THE STRESSES WILL BE OUTPUT IN THE
!     ELEMENT COORDINATE SYSTEM.
!     IF Q4STRS IS SET EQUAL TO 1, STRESSES WILL BE OUTPUT IN THE E C.S.
!     WHICH COOINCIDES WITH MSC'S  VERSION OF ELEMENT COORDINATE SYSTEM.
       
       1300 thetas = est(26)*degrad
       1310 IF (q4strs == 1) GO TO 1320
       CALL angtrs (thetas,0,tsu)
       CALL gmmats (tsu,3,3,0, teu,3,3,1, tse)
       GO TO 1330
       1320 CALL angtrs (thetas,0,tse)
!                                   T
!     CALCULATE  TSI  = TSE X PSITRN  AND STORE IT IN PHIOUT
       
       1330 CALL gmmats (tse,3,3,0, psitrn,3,3,1, phiout(icount+2))
       
!     FOR CORNER POINTS (THE STRESS EVALUATION POINTS EXCEPT FOR THE
!     ONES AT THE CENTER), CALCULATE TSB-MATRIX AND STORE IT IN IELOUT.
       
       IF (ixsi+ieta > 4) GO TO 1340
       ip  = (ixsi-1)*2 + ieta
       ip1 = ipn(ip)
       ip2 = (ip1-1)*25 + 4 + (izta-1)*9
       CALL gmmats (tse,3,3,0, teb,3,3,0, relout(ip2))
       1340 CONTINUE
       
!     CALL Q4BMGS TO GET B MATRIX
!     SET THE ROW FLAG TO 3 TO CREATE THE FIRST 6 ROWS. THEN SET IT
!     TO 1 FOR THE LAST 2 ROWS.
       
       rowflg = 3
       CALL q4bmgs (dshp,gpth,egpdt,epnorm,phi,bmatrx(1))
       DO  ix = 1,ndof
         bmatrx(ix+nd2) = xybmat(ibot+ix)
       END DO
       
       IF (.NOT.bendng) GO TO 1370
       rowflg = 1
       CALL q4bmgs (dshp,gpth,egpdt,epnorm,phi,bmatrx(1+nd6))
       DO  ix = 1,ndof
         bmatrx(ix+nd5) = xybmat(ibot+ix+ndof)
       END DO
       1370 CONTINUE
       
       
!     HERE WE SHIP OUT THE STRAIN RECOVERY MATRIX.
!     --------------------------------------------
       
       kcount = icount + 32 + nnode
       DO  iph = 1,nd8
         phiout(kcount+iph) = bmatrx(iph)
       END DO
     END DO
   END DO
 END DO
 RETURN
 
 1430 nogo = 1
 RETURN
 
 1440 CALL mesage (30,j,NAME)
 GO TO 1430
 
!     BEGINNING OF HEAT RECOVERY.
 
 1500 CONTINUE
 matid    = nest(13)
 inflag   = 2
 nphi(22) = 2
 nphi(23) = nnode
 nphi(24) = NAME(1)
 nphi(25) = NAME(2)
 xi  = 0.0
 eta = 0.0
 CALL q4shps (xi,eta,shp,dshp)
 
!     SORT THE SHAPE FUNCTIONS AND THEIR DERIVATIVES INTO SIL ORDER.
 
 DO  i = 1,4
   tmpshp(i  ) = shp (i  )
   dshptp(i  ) = dshp(i  )
   dshptp(i+4) = dshp(i+4)
 END DO
 DO  i = 1,4
   kk = iorder(i)
   shp (i  ) = tmpshp(kk  )
   dshp(i  ) = dshptp(kk  )
   dshp(i+4) = dshptp(kk+4)
 END DO
 
 hzta = 0.0
 CALL jacobs (elid,shp,dshp,gpth,egpdt,epnorm,jacobe)
 IF (badjac) GO TO 1430
 
 DO  i = 2,4
   ecpt(i) = 0.0
   DO  j = 1,nnode
     ecpt(i) = ecpt(i) + shp(j)*bgpdt(i,j)
   END DO
 END DO
 
 flags = nest(27)
 IF (flags == 0) GO TO 1580
 scsid = nest(26)
 IF (scsid <= 0) GO TO 1540
 necpt(1) = scsid
 CALL transs (ecpt,tbs)
 CALL gmmats (tbs,3,3,1, tub,3,3,1, tsu)
 GO TO 1560
 1540 DO  i = 1,3
   ii = (i-1)*3
   DO  j = 1,3
     jj = (j-1)*3
     tsu(ii+j) = tub(i+jj)
   END DO
 END DO
 1560 CONTINUE
 xs = tsu(1)
 ys = tsu(2)
 IF (ABS(xs) > eps1 .OR. ABS(ys) > eps1) GO TO 1570
 nest(2) = scsid
 j = 233
 GO TO 1440
 1570 thetas = ATAN2(ys,xs)
 GO TO 1590
 1580 thetas = est(26)*degrad
 1590 CALL angtrs (thetas,0,tsu)
 sinmat = 0.0
 cosmat = 1.0
 CALL hmat (elid)
 phiout(26) = kheat(1)
 phiout(27) = kheat(2)
 phiout(28) = kheat(2)
 phiout(29) = kheat(3)
 
!     BRANCH IF THERMAL CONDUCTIVITY KHEAT IS ISOTROPIC.
!     OTHERWISE, FIND TBM, TBS AND TMS AND COMPUTE THE KHEAT
!     TENSOR IN 2-DIMENSIONAL STRESS COORDINATE SYSTEM.
 
!     COMMENTS FROM G.CHAN/UNISYS     10/88
!     HMAT ROUTINE DOES NOT RETURN 'TYPE' IN COSMIC NASTRAN
!     SO WE CAN ONLY ASSUME THERMAL CONDUCTIVITY IS ISOTROPIC AND
!     BRANCH TO 1610 UNCONDITIOANLLY BY SETTING TYPE =-1
 
 TYPE =-1
 
 IF (TYPE == 4 .OR. TYPE == -1) GO TO 1610
 GO TO 290
 1600 CONTINUE
 CALL gmmats (tum,3,3,1, tsu,3,3,1, tms)
 tms(3) = tms(4)
 tms(4) = tms(5)
 CALL gmmats (tms,2,2,1, phiout(26),2,2,0, tum)
 CALL gmmats (tum,2,2,0, tms,2,2,0, phiout(26))
 1610 CONTINUE
 CALL gmmats (teu,3,3,1, jacobe,3,3,0, jacobu)
 CALL gmmats (tsu,3,3,0, jacobu,3,3,0, jacbs)
 DO  j = 1,nnode
   dq(j) = dshp(j)
   jn = j + nnode
   dq(jn) = dshp(j+4)
   jn = jn + nnode
   dq(jn) = 0.0
 END DO
 CALL gmmats (jacbs,3,3,0, dq,3,nnode,0, phiout(35))
 RETURN
 
 1700 FORMAT (a23,', QUAD4 ELEMENT HAS UNDEFINED THICKNESS.  ELEMENT',  &
     ' ID =',i8,', SIL ID =',i8)
 1710 FORMAT (a23,', MODULE SDR2 DETECTS BAD OR REVERSE GEOMETRY FOR ',  &
     'ELEMENT ID =',i8)
END SUBROUTINE squd41
