SUBROUTINE quad4d
     
!     FORMS STIFFNESS AND MASS MATRICES FOR THE QUAD4 PLATE ELEMENT
 
!     DOUBLE PRECISION VERSION
 
!     EST  LISTING
 
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
 
 LOGICAL :: heat,membrn,bendng,shrflx,mbcoup,norpth,badjac, anis,nocsub,nogo
 
 INTEGER :: nest(45),iegpdt(4,4),cpmass,flags,nout,eltype,    &
            elid,estid,sil(4),ksil(4),kcid(4),dict(9),        &
            igpdt(4,4),igpth(4),nam(2),mid(4),TYPE,necpt(4),  &
            rowflg,notran(4),hsil(8),horder(8)
            
 REAL :: tsfact,epsi,epst,eps,gpth(4),matout,egpdt(4,4),  &
         gsube,bgpdm(3,4),gpnorm(4,4),bgpdt(4,4),adamp,   &
         matset,nsm,epnorm(4,4),kheat,htcp,sinmat,cosmat, ecpt(4),SAVE(20)
         
 DOUBLE PRECISION :: amgg(1),akgg,dgpth(4),bmat1(384),xybmat(96),  &
                     zeta,mominr,vol,voli,th,area,area2,detj,      &
                     ptint(2),eps1,xi,eta,zta,hzta,thk,            &
                     xmasso,v(3,3),coeff,xmtmp(16),xmass(16),      &
                     tmpmas(9),jacob(3,3),tmpshp(4),tmpthk(4),     &
                     dshptp(8),psitrn(9),phi(9),shp(4),dshp(8),    &
                     tgrid(4,4),colstf(144),trans(36),trans1(36),  &
                     coltmp(144),avgthk,temp
 
 DOUBLE PRECISION :: vkl, v12dk, vp12l, vjl
 DOUBLE PRECISION :: dnux, dnuy
 
!     DATA FOR ADDING ELEMENT, USER AND MATERIAL COORDINATE SYSTEMS
 
 DOUBLE PRECISION :: aa,bb,cc,x31,y31,x42,y42,exi,exj,ugpdm(3,4),  &
                     cent(3),cente(3),tbm(9),teb(9),tem(9),tub(9), &
                     tum(9),teu(9),tbg(9),gge(9),ggu(9)
 
!     DATA FOR ADDING CSUBB, MIDI, MATERIAL TRANS., AND HEAT
 
 DOUBLE PRECISION :: rho,ts,tsi,reali,rhox,thetam,xm,ym,u(9),a,b,      &
                     aspect,thlen,xa(4),yb(4),GT(9),gi(36),            &
                     enorx,enory,gnorx,gnory,nunorx,nunory,dsub,dsub4, &
                     psiinx,psiiny,tsmfx,tsmfy,curvtr(3,4),curve(3),   &
                     sineax,sineay,w1,pi,twopi,raddeg,degrad,          &
                     htflx(12),htcap(16),htcon(16),dvol,dheat,weitc,   &
                     bterms(32),determ
 
 DOUBLE PRECISION :: vd1(3), vd2(3), vkn(3), vks(3),  &
                     v12(3), v41(3), vp12(3),vis(3), vjs(3)
 
!     DATA FOR IRREGULAR 4-NODE
 
 DOUBLE PRECISION :: zc(4),uev,anglei,edgel,edgshr,unv,vnt(3,4),const, &
                     aspctx,aspcty,gfour(10,10),dfour(7,7),bfour(240), &
                     csubb4,csubx,csuby,csubt,csubtx,csubty,offset,    &
                     sfctr1,sfctr2,sfctx1,sfctx2,sfcty1,sfcty2
     
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg /  ufm,uwm,uim,sfm
 
!     ICORE = FIRST WORD OF OPEN CORE
!     JCORE = NEXT AVAILABLE LOCATION IN OPEN CORE.
!     NCORE = CURRENT LAST AVAILABLE LOCATION IN OPEN CORE
 
 COMMON /emgprm/  icore,jcore,ncore,icstm,ncstm,imat,nmat,ihmat, &
                  nhmat,idit,ndit,icong,ncong,lcong,anycon,      &
                  flags(3),precis,error,heat,cpmass,lcstm,lmat,  &
                  lhmat,kflags(3),l38
 COMMON /emgest/  est(45)
 COMMON /emgdic/  eltype,ldict,nlocs,elid,estid
 COMMON /system/  sys(100)
 COMMON /matin /  matid,inflag,eltemp,dummy,sinmat,cosmat
 COMMON /matout/  matout(25)
 COMMON /hmtout/  kheat(7),TYPE
 COMMON /zzzzzz/  akgg(20000)
 COMMON /q4dt  /  detj,hzta,psitrn,nnode,badjac,n1
 COMMON /terms /  membrn,bendng,shrflx,mbcoup,norpth
 COMMON /q4comd/  anglei(4),edgshr(3,4),edgel(4),unv(3,4),  &
                  uev(3,4),rowflg,iorder(4)
 COMMON /condad/  pi,twopi,raddeg,degrad
 COMMON /comjac/  xi,eta,zeta,determ,dum2,ltypfl
 COMMON /cjacob/  th,vi(3),vj(3),vn(3)
 COMMON /trplm /  ndof,ibot,iptx1,iptx2,ipty1,ipty2
 
 EQUIVALENCE      (sys(01) ,sysbuf  ), (sys(02) ,nout      ),  &
                  (sys(03) ,nogo    ), (sys(55) ,iprec     )
 EQUIVALENCE      (flags(1),kgg1    ), (flags(2),mgg1      ),  &
                  (adamp   ,dict(5) ), (igpth(1),gpth(1)   ),  &
                  (est(1)  ,nest(1) ), (INT     ,nest(25)  ),  &
                  (elth    ,est(14) ), (gpth(1) ,est(6)    ),  &
                  (zoff    ,est(12) ), (zoff1   ,est(28)   ),  &
                  (sil(1)  ,nest(2) ), (matset  ,matout(25)),  &
                  (nsm     ,est(19) ), (amgg(1) ,akgg(1)   ),  &
                  (htcp    ,kheat(4)), (htflx(1),tmpmas(1) ),  &
                  (htcap(1),xmass(1)), (htcon(1),xmtmp(1)  ),  &
                  (necpt(1),ecpt(1) ), (bgpdt(1,1) ,est(29)),  &
                  (iegpdt(1,1),egpdt(1,1)),  &
                  (igpdt(1,1) ,bgpdt(1,1))

 DATA    eps1  /  1.0D-7 /
 DATA    const /  0.57735026918962D0 /
 DATA    nam   /  4HQUAD, 4H4D       /
 
 elid   = nest(1)
 ltypfl = 1
 offset = zoff
 IF (zoff == 0.0) offset = zoff1
 
!     CHECK FOR SUFFICIENT OPEN CORE FOR ELEMENT STIFFNESS
 
 jcored = jcore/iprec + 1
 ncored = ncore/iprec - 1
 IF ((jcored+576) <= ncored .OR. heat .OR. kgg1 == 0) GO TO 10
 GO TO 1730
 
!     COPY THE SILS AND BGPDT DATA INTO SAVE ARRAY SINCE THE DATA
!     WILL BE REORDERED BASED ON INCREASING SILS.
 
 10 j = 1
 DO  i = 1,20
   SAVE(i) = est(i+j)
   IF (i == 4) j = 24
 END DO
 
 nnode = 4
 n1    = 4
 nodesq= nnode*nnode
 ndof  = nnode*6
 ndof3 = nnode*3
 nd2   = ndof*2
 nd3   = ndof*3
 nd4   = ndof*4
 nd5   = ndof*5
 nd6   = ndof*6
 nd7   = ndof*7
 
!     FILL IN ARRAY GGU WITH THE COORDINATES OF GRID POINTS 1, 2 AND 4.
!     THIS ARRAY WILL BE USED LATER TO DEFINE THE USER COORD. SYSTEM
!     WHILE CALCULATING TRANSFORMATIONS INVOLVING THIS COORD. SYSTEM.
 
 DO  i = 1,3
   ii = (i-1)*3
   ij = i
   IF (ij == 3) ij = 4
   DO  j = 1,3
     jj = j + 1
     ggu(ii+j) = bgpdt(jj,ij)
   END DO
 END DO
 
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
 vkl = DSQRT( vkn(1)**2 + vkn(2)**2 + vkn(3)**2 )
 IF ( vkl == 0. ) WRITE( nout, 2070 ) nest(1)
 2070  FORMAT(//,' ILLEGAL GEOMETRY FOR QUAD4 ELEMENT, ID=',i10 )
 vks(1) = vkn(1)/vkl
 vks(2) = vkn(2)/vkl
 vks(3) = vkn(3)/vkl
 pa = vkl/2.d0
 
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
 vp12l   = DSQRT( vp12(1)**2 + vp12(2)**2 + vp12(3)**2 )
 IF ( vp12l == 0. ) WRITE( nout, 2070 ) nest(1)
 vis(1) = vp12(1) / vp12l
 vis(2) = vp12(2) / vp12l
 vis(3) = vp12(3) / vp12l
 vjs(1) = vks(2)*vis(3) - vks(3)*vis(2)
 vjs(2) = vks(3)*vis(1) - vks(1)*vis(3)
 vjs(3) = vks(1)*vis(2) - vks(2)*vis(1)
 
!   NORMALIZE J FOR GOOD MEASURE
 
 vjl = DSQRT( vjs(1)**2 + vjs(2)**2 + vjs(3)**2 )
 IF ( vjl == 0. ) WRITE ( nout, 2070 ) nest(1)
 vjs(1) = vjs(1) / vjl
 vjs(2) = vjs(2) / vjl
 vjs(3) = vjs(3) / vjl
 DO  i = 1,3
   tub(i)   = vis(i)
   tub(i+3) = vjs(i)
   tub(i+6) = vks(i)
 END DO
 
!     STORE INCOMING BGPDT FOR LUMPED MASS AND ELEMENT C.S.
 
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
     ugpdm(i,j) = 0.0D0
     DO  k = 1,3
       kk = ip + k
       ugpdm(i,j) = ugpdm(i,j) + tub(kk)*(DBLE(bgpdm(k,j))-ggu(k))
     END DO
   END DO
 END DO
 
 
!     THE ORIGIN OF THE ELEMENT C.S. IS IN THE MIDDLE OF THE ELEMENT
 
 DO  j = 1,3
   cent(j) = 0.0D0
   DO  i = 1,4
     cent(j) = cent(j)+ugpdm(j,i)/nnode
   END DO
 END DO
 
!     STORE THE CORNER NODE DIFF. IN THE USER C.S.
 
 x31 = ugpdm(1,3) - ugpdm(1,1)
 y31 = ugpdm(2,3) - ugpdm(2,1)
 x42 = ugpdm(1,4) - ugpdm(1,2)
 y42 = ugpdm(2,4) - ugpdm(2,2)
 aa  = DSQRT(x31*x31 + y31*y31)
 bb  = DSQRT(x42*x42 + y42*y42)
 IF (aa == 0.d0 .OR. bb == 0.d0) GO TO 1700
 
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
 
!     THE ARRAY IORDER STORES THE ELEMENT NODE ID IN
!     INCREASING SIL ORDER.
 
!     IORDER(1) = NODE WITH LOWEST  SIL NUMBER
!     IORDER(4) = NODE WITH HIGHEST SIL NUMBER
 
!     ELEMENT NODE NUMBER IS THE INTEGER FROM THE NODE
!     LIST  G1,G2,G3,G4 .  THAT IS, THE 'I' PART
!     OF THE 'GI' AS THEY ARE LISTED ON THE CONNECTIVITY
!     BULK DATA CARD DESCRIPTION.
 
 DO  i = 1,4
   iorder(i) = 0
   horder(i) = 0
   ksil(i) = sil(i)
   hsil(i) = sil(i)
 END DO
 
 DO  i = 1,4
   itemp = 1
   isil  = ksil(1)
   DO  j = 2,4
     IF (isil <= ksil(j)) CYCLE
     itemp = j
     isil  = ksil(j)
   END DO
   iorder(i) = itemp
   horder(i) = itemp
   ksil(itemp) = 99999999
 END DO
 
!     ADJUST EST DATA
 
!     USE THE POINTERS IN IORDER TO COMPLETELY REORDER THE
!     GEOMETRY DATA INTO INCREASING SIL ORDER.
!     DON'T WORRY!! IORDER ALSO KEEPS TRACK OF WHICH SHAPE
!     FUNCTIONS GO WITH WHICH GEOMETRIC PARAMETERS!
 
 DO  i = 1,4
   ksil(i) = sil(i)
   tmpthk(i) = gpth(i)
   kcid(i) = igpdt(1,i)
   DO  j = 2,4
     tgrid(j,i) = bgpdt(j,i)
   END DO
 END DO
 DO  i = 1,4
   ipoint  = iorder(i)
   sil(i)  = ksil(ipoint)
   gpth(i) = tmpthk(ipoint)
   igpdt(1,i) = kcid(ipoint)
   DO  j = 2,4
     bgpdt(j,i) = tgrid(j,ipoint)
   END DO
 END DO
 
!     COMPUTE NODE NORMALS
 
 CALL q4nrmd (bgpdt,gpnorm,iorder,iflag)
 IF (iflag == 0) GO TO 130
 GO TO 1700
 
!     DETERMINE NODAL THICKNESSES
 
 130 avgthk = 0.0D0
 DO  i = 1,nnode
   iord = iorder(i)
   DO  ic = 1,3
     curvtr(ic,iord) = gpnorm(ic+1,i)
   END DO
   
   IF (gpth(i) == 0.0) gpth(i) = elth
   IF (nest(13) == 0 .AND. elth == 0.0) gpth(i) = 1.0E-14
   IF (gpth(i) > 0.0) GO TO 150
   WRITE (nout,2010) ufm,elid
   nogo = .true.
   GO TO 1710
   150 dgpth(i) = gpth(i)
   avgthk = avgthk + dgpth(i)/nnode
 END DO
 
!     NEST(13) = MID1 ID FOR MEMBRANE
!     NEST(15) = MID2 ID FOR BENDING
!     NEST(17) = MID3 ID FOR TRANSVERSE SHEAR
!     NEST(22) = MID4 ID FOR MEMBRANE-BENDING COUPLING
!                MID4 MUST BE BLANK UNLESS MID1 AND MID2 ARE NON-ZERO
!                MID4 ID MUST NOT EQUAL MID1 OR MID2 ID
!     (WHEN LAYER COMPOSITE IS USED, MID ID IS RAISED TO ID*100000000)
!      EST(14) = MEMBRANE THICKNESS, T
!      EST(16) = BENDING STIFFNESS PARAMETER, 12I/T**3
!      EST(18) = TRANSVERSE SHEAR  PARAMETER, TS/T
 
!     0.8333333 = 5.0/6.0
 
 mominr = 0.0D0
 tsfact = .8333333
 nocsub = .false.
 IF (nest(15) /=  0) mominr = est(16)
 IF (nest(17) /=  0) ts = est(18)
 IF ( est(18) == .0) ts = .833333D0
 
!     FIX FOR LAMINATED COMPOSITE WITH MEMBRANE BEHAVIOUR ONLY.
!     REQUIRED TO PREVENT ZERO DIVIDE ERRORS.
 
 IF (nest(15) == 0 .AND. nest(13) > 100000000) ts = .833333D0
 
!     SET LOGICAL NOCSUB IF EITHER MOMINR OR TS ARE NOT DEFAULT
!     VALUES. THIS WILL BE USED TO OVERRIDE ALL CSUBB COMPUTATIONS.
!     I.E. DEFAULT VALUES OF UNITY ARE USED.
 
 epsi = ABS(mominr - 1.0)
 epst = ABS(ts  - tsfact)
 eps  = .05
!     NOCSUB = EPSI.GT.EPS .OR. EPST.GT.EPS
 IF (nest(13) > 100000000) nocsub = .false.
 
!     THE COORDINATES OF THE ELEMENT GRID POINTS HAVE TO BE
!     TRANSFORMED FROM THE BASIC C.S. TO THE ELEMENT C.S.
 
 CALL betrnd (teu,gge,0,elid)
 CALL gmmatd (teu,3,3,0,tub ,3,3,0,teb  )
 CALL gmmatd (tub,3,3,1,cent,3,1,0,cente)
 identt = 0
 IF (teb(1) == 1.d0 .AND. teb(5) == 1.d0 .AND. teb(9) == 1.d0 .AND.  &
     teb(2) == 0.d0 .AND. teb(3) == 0.d0 .AND. teb(4) == 0.d0 .AND.  &
     teb(6) == 0.d0 .AND. teb(7) == 0.d0 .AND. teb(8) == 0.d0 ) identt = 1
 ip = -3
 DO  ii = 2,4
   ip = ip + 3
   DO  j = 1,nnode
     epnorm(ii,j) = 0.0
     egpdt(ii,j)  = 0.0
     DO  k = 1,3
       kk = ip + k
       k1 = k  + 1
       cc = DBLE(bgpdt(k1,j)) - ggu(k)-cente(k)
       epnorm(ii,j) = epnorm(ii,j) + teb(kk)*gpnorm(k1,j)
       egpdt(ii,j)  = egpdt(ii,j)  + SNGL(teb(kk)*cc)
     END DO
   END DO
 END DO
 
!     BEGIN INITIALIZING MATERIAL VARIABLES
 
!     SET INFLAG = 12 SO THAT SUBROUTINE MAT WILL SEARCH FOR-
!     ISOTROPIC MATERIAL PROPERTIES AMONG THE MAT1 CARDS,
!     ORTHOTROPIC MATERIAL PROPERTIES AMONG THE MAT8 CARDS, AND
!     ANISOTROPIC MATERIAL PROPERTIES AMONG THE MAT2 CARDS.
 
 inflag = 12
 rho    = 0.0D0
 eltemp = est(45)
 mid(1) = nest(13)
 mid(2) = nest(15)
 mid(3) = nest(17)
 mid(4) = nest(22)
 membrn = mid(1) > 0
 bendng = mid(2) > 0 .AND. mominr > 0.0D0
 shrflx = mid(3) > 0
 mbcoup = mid(4) > 0
 
!     FIGURE OUT PATH OF THE TRIPLE MULTIPLY AND THE NO. OF ROWS IN
!     THE B-MATRIX (I.E. STRAIN-NODAL DISPLACEMENT MATRIX)
 
!     NORPTH = MID(1).EQ.MID(2) .AND. MID(1).EQ.MID(3) .AND. MID(4).EQ.0
!    1        .AND. DABS(MOMINR-1.0D0).LE.EPS1
 
 norpth = .false.
 
!     DETERMINE FACTORS TO BE USED IN CSUBB CALCULATIONS
 
!     IF (.NOT.BENDNG) GO TO 290
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
 
 a = 0.5D0*DABS(xa(2)+xa(3)-xa(1)-xa(4))
 b = 0.5D0*DABS(yb(4)+yb(3)-yb(1)-yb(2))
 IF (a > b) aspect = b/a
 IF (a <= b) aspect = a/b
 thlen = avgthk/a
 IF (a < b) thlen = avgthk/b
 
!     TORSION-RELATED SHEAR CORRECTION FOR 4-NODE-
!     PRELIMINARY FACTORS
 
 aspctx = a/b
 aspcty = b/a
 csubb4 = 1.6D0
 csubt  = 71.d0*aspect*(1.6D0/csubb4)*(1.d0+415.d0*aspect*thlen**2)
 csubtx = csubt*aspctx**2
 csubty = csubt*aspcty**2
 
 i  = 2
 j  = 2
 jj = 3
 sineax = 0.0D0
 sineay = 0.0D0
 220 CALL daxb (curvtr(1,i-1),curvtr(1,i),curve)
 cc = curve(1)*curve(1) + curve(2)*curve(2) + curve(3)*curve(3)
 IF (cc < eps1) GO TO 230
 cc = 0.5D0*DSQRT(cc)
 230 sineax = sineax + cc
 IF (i /= 2) GO TO 240
 i = 4
 GO TO 220
 
 240 CALL daxb (curvtr(1,j),curvtr(1,jj),curve)
 cc = curve(1)*curve(1) + curve(2)*curve(2) + curve(3)*curve(3)
 IF (cc < eps1) GO TO 250
 cc = 0.5D0*DSQRT(cc)
 250 sineay = sineay+cc
 IF (j /= 2) GO TO 260
 j  = 1
 jj = 4
 GO TO 240
 260 cc = 28.0D0
 sineax = cc*sineax + 1.0D0
 sineay = cc*sineay + 1.0D0
 IF (sineax > sineay) sineay = sineax
 IF (sineay > sineax) sineax = sineay
 
!     IRREGULAR 4-NODE CODE-  GEOMETRIC VARIABLES
 
!     CALCULATE AND NORMALIZE- UNIT EDGE VECTORS, UNIT NORMAL VECTORS
 
 DO  i = 1,4
   j = i + 1
   IF (j == 5) j = 1
   uev(1,i) = xa(j) - xa(i)
   uev(2,i) = yb(j) - yb(i)
   uev(3,i) = zc(j) - zc(i)
   unv(1,i) = (vnt(1,j) + vnt(1,i))*0.50D0
   unv(2,i) = (vnt(2,j) + vnt(2,i))*0.50D0
   unv(3,i) = (vnt(3,j) + vnt(3,i))*0.50D0
   cc = uev(1,i)**2 + uev(2,i)**2 + uev(3,i)**2
   IF (cc == 0.d0) GO TO 1700
   IF (cc >= eps1) cc = DSQRT(cc)
   edgel(i) = cc
   uev(1,i) = uev(1,i)/cc
   uev(2,i) = uev(2,i)/cc
   uev(3,i) = uev(3,i)/cc
   cc = unv(1,i)**2 + unv(2,i)**2 + unv(3,i)**2
   IF (cc == 0.d0) GO TO 1700
   IF (cc >= eps1) cc = DSQRT(cc)
   unv(1,i) = unv(1,i)/cc
   unv(2,i) = unv(2,i)/cc
   unv(3,i) = unv(3,i)/cc
 END DO
 
!     CALCULATE INTERNAL NODAL ANGLES
 
 DO  i = 1,4
   j = i - 1
   IF (j == 0) j = 4
   anglei(i)=-uev(1,i)*uev(1,j) -uev(2,i)*uev(2,j) -uev(3,i)*uev(3,j)
   IF (DABS(anglei(i))  < eps1) anglei(i) = 0.0D0
 END DO
 
!     SET THE INTEGRATION POINTS
 
 ptint(1)  = -const
 ptint(2)  =  const
 IF (heat) GO TO 1790
 
!     TRIPLE LOOP TO SAVE THE LAST 2 ROWS OF B-MATRIX AT 2X2X2
!     INTEGRATION POINTS FOR LATER MANIPULATION.
 
 IF (kgg1 == 0) GO TO 400
 i  = 1
 kpt= 1
 
 DO  ixsi = 1,2
   xi = ptint(ixsi)
   
   DO  ieta = 1,2
     eta = ptint(ieta)
     
     CALL q4shpd (xi,eta,shp,dshp)
     
!     IRREGULAR 4-NODE CODE-  CALCULATION OF NODAL EDGE SHEARS
!                             AT THIS INTEGRATION POINT
     
     DO  ij = 1,4
       ii = ij - 1
       IF (ii == 0) ii = 4
       ik = ij + 1
       IF (ik == 5) ik = 1
       aa = shp(ij)
       bb = shp(ik)
       
       DO  is = 1,3
         edgshr(is,ij) = (uev(is,ij)+anglei(ij)*uev(is,ii))*aa/  &
             (1.0D0-anglei(ij)*anglei(ij))  &
             + (uev(is,ij)+anglei(ik)*uev(is,ik))*bb/  &
             (1.0D0-anglei(ik)*anglei(ik))
       END DO
     END DO
     
!     SORT THE SHAPE FUNCTIONS AND THEIR DERIVATIVES INTO SIL ORDER.
     
     DO  is = 1,4
       tmpshp(is  ) =  shp(is  )
       dshptp(is  ) = dshp(is  )
       dshptp(is+4) = dshp(is+4)
     END DO
     DO  is = 1,4
       kk = iorder(is)
       shp (is  ) = tmpshp(kk  )
       dshp(is  ) = dshptp(kk  )
       dshp(is+4) = dshptp(kk+4)
     END DO
     
     DO  izta = 1,2
       zta = ptint(izta)
       
!     COMPUTE THE JACOBIAN AT THIS GAUSS POINT,
!     ITS INVERSE AND ITS DETERMINANT.
       
       hzta = zta/2.0D0
       CALL jacob2 (elid,shp,dshp,dgpth,egpdt,epnorm,jacob)
       IF (badjac) GO TO 1710
       
!     COMPUTE PSI TRANSPOSE X JACOBIAN INVERSE.
!     HERE IS THE PLACE WHERE THE INVERSE JACOBIAN IS FLAGED TO BE
!     TRANSPOSED BECAUSE OF OPPOSITE MATRIX LOADING CONVENTION
!     BETWEEN INVER AND GMMAT.
       
       CALL gmmatd (psitrn,3,3,0,jacob,3,3,1,phi)
       
!     CALL Q4BMGD TO GET B MATRIX
!     SET THE ROW FLAG TO 1. IT SIGNALS SAVING THE LAST 2 ROWS.
       
       rowflg = 1
       CALL q4bmgd (dshp,dgpth,egpdt,epnorm,phi,bmat1(kpt))
       kpt = kpt + nd2
     END DO
   END DO
 END DO
 
!     IN PLANE SHEAR REDUCTION
 
 xi  = 0.0D0
 eta = 0.0D0
 kpt = 1
 kpnt= nd2
 
 CALL q4shpd (xi,eta,shp,dshp)
 
!     SORT THE SHAPE FUNCTIONS AND THEIR DERIVATIVES INTO SIL ORDER.
 
 DO  i = 1,4
   tmpshp(i  ) =  shp(i  )
   dshptp(i  ) = dshp(i  )
   dshptp(i+4) = dshp(i+4)
 END DO
 DO  i = 1,4
   kk = iorder(i)
   shp(i   ) = tmpshp(kk  )
   dshp(i  ) = dshptp(kk  )
   dshp(i+4) = dshptp(kk+4)
 END DO
 
 DO  izta = 1,2
   zta  = ptint(izta)
   hzta = zta/2.0D0
   CALL jacob2 (elid,shp,dshp,dgpth,egpdt,epnorm,jacob)
   IF (badjac) GO TO 1710
   
   CALL gmmatd (psitrn,3,3,0,jacob,3,3,1,phi)
   
!     CALL Q4BMGD TO GET B-MATRIX
!     SET THE ROW FLAG TO 2. IT WILL SAVE THE 3RD ROW OF B-MATRIX AT
!     THE TWO INTEGRATION POINTS.
   
   rowflg = 2
   CALL q4bmgd (dshp,dgpth,egpdt,epnorm,phi,xybmat(kpt))
   kpt = kpt + kpnt
 END DO
 
!     SET THE ARRAY OF LENGTH 4 TO BE USED IN CALLING TRANSD.
!     NOTE THAT THE FIRST WORD IS THE COORDINATE SYSTEM ID WHICH
!     WILL BE SET IN POSITION LATER.
 
 400 DO  iec = 2,4
   ecpt(iec) = 0.0
 END DO
 
!     FETCH MATERIAL PROPERTIES
 
 
!     EACH MATERIAL PROPERTY MATRIX G HAS TO BE TRANSFORMED FROM
!     THE MATERIAL COORDINATE SYSTEM TO THE ELEMENT COORDINATE
!     SYSTEM. THESE STEPS ARE TO BE FOLLOWED-
 
!     1- IF MCSID HAS BEEN SPECIFIED, SUBROUTINE TRANSD IS CALLED
!        TO CALCULATE TBM-MATRIX (MATERIAL TO BASIC TRANSFORMATION).
!        TBM-MATRIX IS THEN PREMULTIPLIED BY TEB-MATRIX TO OBTAIN
!        TEM-MATRIX.
!        THEN USING THE PROJECTION OF X-AXIS, AN ANGLE IS CALCULATED
!        UPON WHICH STEP 2 IS TAKEN.
 
!     2- IF THETAM HAS BEEN SPECIFIED, SUBROUTINE ANGTRD IS CALLED
!        TO CALCULATE TEM-MATRIX (MATERIAL TO ELEMENT TRANSFORMATION).
 
!                          T
!     3-           G  =   U   G   U
!                   E          M
 
 
 IF (nest(11) == 0) GO TO 470
 mcsid = nest(10)
 
!     CALCULATE TEM-MATRIX USING MCSID
 
 420 IF (mcsid > 0) GO TO 440
 DO  i = 1,9
   tem(i) = teb(i)
 END DO
 GO TO 450
 440 necpt(1) = mcsid
 CALL transd (ecpt,tbm)
 
!     MULTIPLY TEB AND TBM MATRICES
 
 CALL gmmatd (teb,3,3,0,tbm,3,3,0,tem)
 
!     CALCULATE THETAM FROM THE PROJECTION OF THE X-AXIS OF THE
!     MATERIAL C.S. ON TO THE XY PLANE OF THE ELEMENT C.S.
 
 450 CONTINUE
 xm = tem(1)
 ym = tem(4)
 IF (DABS(xm) > eps1 .OR. DABS(ym) > eps1) GO TO 460
 nest(2) = mcsid
 j = 231
 GO TO 1720
 460 thetam = DATAN2(ym,xm)
 GO TO 480
 
!     CALCULATE TEM-MATRIX USING THETAM
 
 470 thetam = DBLE(est(10))*degrad
 IF (thetam == 0.0D0) GO TO 490
 480 CALL angtrd (thetam,1,tum)
 CALL gmmatd (teu,3,3,0,tum,3,3,0,tem)
 GO TO 510
 
!     DEFAULT IS CHOSEN, LOOK FOR VALUES OF MCSID AND/OR THETAM
!     ON THE PSHELL CARD.
 
 490 IF (nest(24) == 0) GO TO 500
 mcsid = nest(23)
 GO TO 420
 
 500 thetam = DBLE(est(23))*degrad
 GO TO 480
 
 510 CONTINUE
 IF (heat) GO TO 1810
 
 DO  m = 1,36
   gi(m)  = 0.0D0
 END DO
 sinmat = 0.
 cosmat = 0.
 igobk  = 0
 
!     BEGIN M-LOOP TO FETCH PROPERTIES FOR EACH MATERIAL ID
 
 m = 0
 610 m = m + 1
 IF (m > 4) GO TO 790
 IF (m == 4 .AND. igobk == 1) GO TO 800
 matid = mid(m)
 IF (matid == 0 .AND. m /= 3) GO TO 610
 IF (matid == 0 .AND. m == 3 .AND. .NOT.bendng) GO TO 610
 IF (matid == 0 .AND. m == 3 .AND. bendng) matid = mid(2)
 
 IF (m-1 < 0) THEN
   GO TO   640
 ELSE IF (m-1 == 0) THEN
   GO TO   630
 END IF
 620 IF (matid == mid(m-1) .AND. igobk == 0) GO TO 640
 630 CALL mat (elid)
 640 CONTINUE
 
 IF (membrn .AND. m == 1) rho = matout(7)
 rhox = rho
 IF (rho == 0.0D0) rhox = 1.0D0
 IF (kgg1 ==    0) GO TO 610
 
 IF (membrn .AND. m /= 1 .OR. .NOT.membrn .AND. m /= 2) GO TO 650
 gsube = matout(12)
 IF (matset == 8.) gsube = matout(16)
 650 CONTINUE
 
 IF (m == 2 .AND. norpth) GO TO 670
 coeff  = 1.0D0
 lpoint = (m-1)*9 + 1
 
 CALL q4gmgd (m,coeff,gi(lpoint))
 
 IF (m == 3) GO TO 770
 u(1) = tem(1)*tem(1)
 u(2) = tem(4)*tem(4)
 u(3) = tem(1)*tem(4)
 u(4) = tem(2)*tem(2)
 u(5) = tem(5)*tem(5)
 u(6) = tem(2)*tem(5)
 u(7) = tem(1)*tem(2)*2.0D0
 u(8) = tem(4)*tem(5)*2.0D0
 u(9) = tem(1)*tem(5) + tem(2)*tem(4)
 l=3
 GO TO 780
 
 770 u(1) = tem(5)*tem(9) + tem(6)*tem(8)
 u(2) = tem(2)*tem(9) + tem(8)*tem(3)
 u(3) = tem(4)*tem(9) + tem(7)*tem(6)
 u(4) = tem(1)*tem(9) + tem(3)*tem(7)
 l    = 2
 
 780 CALL gmmatd (u(1),l,l,1,gi(lpoint),l,l,0,GT(1))
 CALL gmmatd (GT(1),l,l,0,u(1),l,l,0,gi(lpoint))
 
 IF (m > 0) GO TO 670
 IF (.NOT.shrflx .AND. bendng) GO TO 660
 nest(2) = matid
 j = 232
 GO TO 1720
 660 m = -m
 670 CONTINUE
 mtype = IFIX(matset+.05) - 2
 IF (nocsub) GO TO 760
 SELECT CASE ( m )
   CASE (    1)
     GO TO 760
   CASE (    2)
     GO TO 680
   CASE (    3)
     GO TO 720
   CASE (    4)
     GO TO 760
 END SELECT
 680 IF ( mtype  < 0) THEN
   GO TO   690
 ELSE IF ( mtype  == 0) THEN
   GO TO   700
 ELSE
   GO TO   710
 END IF
 690 enorx = matout(16)
 enory = matout(16)
 dnux  = gi( lpoint+1 ) / gi( lpoint )
 dnuy  = gi( lpoint+3 ) / gi( lpoint+4 )
 GO TO 760
 700 enorx = matout(1)
 enory = matout(4)
 dnux  = gi( lpoint+1 ) / gi( lpoint )
 dnuy  = gi( lpoint+3 ) / gi( lpoint+4 )
 GO TO 760
 710 enorx = matout(1)
 enory = matout(3)
 dnux  = gi(lpoint+1)/gi(lpoint)
 dnuy  = gi(lpoint+3)/gi(lpoint+4)
 GO TO 760
 720 IF ( mtype  < 0) THEN
   GO TO   730
 ELSE IF ( mtype  == 0) THEN
   GO TO   740
 ELSE
   GO TO   750
 END IF
 730 gnorx = matout(6)
 gnory = matout(6)
 GO TO 760
 740 gnorx = matout(1)
 gnory = matout(4)
 GO TO 760
 750 gnorx = matout(6)
 gnory = matout(5)
 IF ( gnorx == 0.0D0 ) gnorx = matout(4)
 IF ( gnory == 0.0D0 ) gnory = matout(4)
 760 CONTINUE
 GO TO 610
 
!     END OF M-LOOP
 
 790 CONTINUE
 IF (mid(3) < 100000000) GO TO 800
 IF (gi(19) /= 0.d0 .OR. gi(20) /= 0.d0 .OR. gi(21) /= 0.d0 .OR.  &
     gi(22) /= 0.d0) GO TO 800
 igobk = 1
 m = 2
 mid(3) = mid(2)
 GO TO 610
 800 CONTINUE
 
 nocsub = enorx == 0.0D0 .OR. enory == 0.0D0 .OR.  &
     gnorx == 0.0D0 .OR. gnory == 0.0D0 .OR. mominr == 0.0D0
 
 mattyp = IFIX(matset+.05)
 
!     IF MGG1 IS NON-ZERO AND RHO IS GREATER THAN 0.0,
!     THEN COMPUTE THE MASS MATRIX.
 
 IF (mgg1 == 0) GO TO 810
 IF (jcored+144 <= ncored) GO TO 810
 810 CONTINUE
 
 limit = jcored + ndof*ndof
 DO  i = jcored,limit
   akgg(i)  = 0.0D0
 END DO
 DO  i = 1,nodesq
   xmass(i) = 0.0D0
   xmtmp(i) = 0.0D0
 END DO
 area     = 0.0D0
 vol      = 0.0D0
 
 
!     HERE BEGINS THE TRIPLE LOOP ON STATEMENTS 1310 AND 1300 TO
!     GAUSS INTEGRATE FOR THE ELEMENT MASS AND STIFFNESS MATRICES.
!     -----------------------------------------------------------
 
 DO  ixsi = 1,2
   xi = ptint(ixsi)
   DO  ieta = 1,2
     eta = ptint(ieta)
     CALL q4shpd (xi,eta,shp,dshp)
     
!     SORT THE SHAPE FUNCTIONS AND THEIR DERIVATIVES INTO SIL ORDER.
     
     DO  i = 1,4
       tmpshp(i  ) =  shp(i  )
       dshptp(i  ) = dshp(i  )
       dshptp(i+4) = dshp(i+4)
     END DO
     DO  i = 1,4
       kk = iorder(i)
       shp (i  ) = tmpshp(kk  )
       dshp(i  ) = dshptp(kk  )
       dshp(i+4) = dshptp(kk+4)
     END DO
     CALL gmmatd (shp,1,nnode,0,dgpth,1,nnode,1,thk)
     reali = mominr*thk*thk*thk/12.0D0
     tsi   = ts*thk
     
!     SKIP MASS CALCULATIONS IF NOT REQUESTED
     
     IF (nsm /=  0.0) GO TO 920
     IF (mgg1 ==   0) GO TO 1020
     IF (rho == 0.d0) GO TO 1020
     IF (rho > 0.d0) GO TO 920
     WRITE (nout,2030) uwm,rho,mid(1),nest(1)
     
!     COMPUTE S AND T VECTORS AT THE MID-SURFACE
!     FOR MASS CALCULATIONS ONLY.
     
     920 CONTINUE
     DO  i = 1,2
       ipoint = 4*(i-1)
       DO  j = 1,3
         v(i,j) = 0.0D0
         DO  k = 1,nnode
           ktemp = k + ipoint
           jtemp = j + 1
           v(i,j)= v(i,j) + dshp(ktemp)*bgpdt(jtemp,k)
         END DO
       END DO
     END DO
     
!     COMPUTE S CROSS T AT THE MID-SURFACE FOR MASS CALCULATIONS.
     
     v(3,1) = v(1,2)*v(2,3) - v(2,2)*v(1,3)
     v(3,2) = v(1,3)*v(2,1) - v(2,3)*v(1,1)
     v(3,3) = v(1,1)*v(2,2) - v(2,1)*v(1,2)
     area2  = v(3,1)*v(3,1) + v(3,2)*v(3,2) + v(3,3)*v(3,3)
     
!     AREA2 = NORM OF S CROSS T IS THE AREA OF THE ELEMENT
!     AS COMPUTED AT THIS GAUSS POINT.
     
     IF ( area2 <= 0.0 ) GO TO 1700
     
     area2 = DSQRT(area2)
     area  = area + area2
     voli  = area2*thk
     vol   = vol + voli
     
     IF (mgg1 ==   0) GO TO 1020
     IF (cpmass > 0) GO TO 1000
     i4 = 1
     DO  j4 = 1,nnode
       xmass(i4) = xmass(i4) + voli*rhox*shp(j4)
       i4 = i4 + nnode + 1
     END DO
     GO TO 1020
     
!     COMPUTE CONSISTENT MASS MATRIX
     
!     COMPUTE THE CONTRIBUTION TO THE MASS MATRIX
!     FROM THIS INTEGRATION POINT.
     
     1000 CALL gmmatd (shp,1,nnode,1,shp,1,nnode,0,xmtmp)
     
!     ADD MASS CONTRIBUTION FROM THIS INTEGRATION POINT
!     TO THE ELEMENT MASS MATRIX.
     
     DO  i = 1,nodesq
       xmass(i) = xmass(i) + voli*rhox*xmtmp(i)
     END DO
     1020 IF (kgg1 == 0) GO TO 1330
     
!     BEGIN STIFFNESS COMPUTATIONS
     
!     SET DEFAULT VALUES OF CSUBB FACTORS
     
     sfcty1 = 1.0D0
     sfcty2 = 1.0D0
     sfctx1 = 1.0D0
     sfctx2 = 1.0D0
     tsmfx  = 1.0D0
     tsmfy  = 1.0D0
     IF (nocsub) GO TO 1090
     IF (.NOT.bendng) GO TO 1090
     nunorx = mominr*enorx/(2.0D0*gnorx) - 1.0D0
     nunory = mominr*enory/(2.0D0*gnory) - 1.0D0
     IF (nunorx > 0.999999D0) nunorx = 0.999999D0
     IF (nunory > 0.999999D0) nunory = 0.999999D0
     IF ( nunorx <= 0. ) nunorx = dnux
     IF ( nunory <= 0. ) nunory = dnuy
     cc = aspect
     w1 = 1.0D0 + 4400.0D0*thlen*thlen*thlen*thlen
     IF (cc < 0.2D0) GO TO 1030
     dsub4 = (18.375D0-11.875D0*cc)*w1
     GO TO 1040
     1030 dsub4 = (159.85D0*cc-15.97D0)*w1
     1040 IF (dsub4 < .01D0) dsub4 = 0.01D0
     IF (dsub4 > 2.0D3) dsub4 = 2000.0D0
     dsub  = dsub4
     coeft = const
     ax    = a
     IF (eta  < 0.0D0) ax = a + coeft*(xa(2)-xa(1)-a)
     IF (eta > 0.0D0) ax = a + coeft*(xa(3)-xa(4)-a)
     psiinx = 20.0D0*dsub*reali*sineax*(1.0D0+aspect*aspect)/  &
         (tsi*(1.0D0-nunorx)*ax*ax)
     dsub  = dsub4
     coeft = const
     by    = b
     IF (xi < 0.0D0) by = b + coeft*(yb(4)-yb(1)-b)
     IF (xi > 0.0D0) by = b + coeft*(yb(3)-yb(2)-b)
     psiiny = 20.0D0*dsub*reali*sineay*(1.0D0+aspect*aspect)/  &
         (tsi*(1.0D0-nunory)*by*by)
     IF (.NOT.shrflx) GO TO 1050
     tsmfx = psiinx/(1.0D0+psiinx)
     tsmfy = psiiny/(1.0D0+psiiny)
     GO TO 1060
     1050 tsmfx = psiinx
     tsmfy = psiiny
     
     1060 CONTINUE
     IF (tsmfx <= 0.0D0) tsmfx = eps1
     IF (tsmfy <= 0.0D0) tsmfy = eps1
     
!     FILL IN THE 7X7 MATERIAL PROPERTY MATRIX D FOR NORPTH
     
     IF (.NOT.norpth) GO TO 1090
     DO  ig = 1,7
       DO  jg = 1,7
         dfour(ig,jg) = 0.0D0
       END DO
     END DO
     
     DO  ig = 1,3
       ig1 = (ig-1)*3
       DO  jg = 1,3
         jg1 = jg + ig1
         dfour(ig,jg) = gi(jg1)
       END DO
     END DO
     GO TO 1150
     
!     FILL IN THE 10X10 G-MATRIX WHEN MID4 IS NOT PRESENT
     
     1090 DO  ig = 1,10
       DO  jg = 1,10
         gfour(ig,jg) = 0.0D0
       END DO
     END DO
     IF (mbcoup) GO TO 1150
     
     IF (.NOT.membrn) GO TO 1120
     DO  ig = 1,3
       ig1 = (ig-1)*3
       DO  jg = 1,3
         jg1 = jg + ig1
         gfour(ig,jg) = gi(jg1)
       END DO
     END DO
     
     1120 IF (.NOT.bendng) GO TO 1250
     DO  ig = 4,6
       ig2 = (ig-2)*3
       DO  jg = 4,6
         jg2 = jg + ig2
         gfour(ig,jg) = gi(jg2)*mominr
       END DO
     END DO
     
     IF (.NOT.membrn) GO TO 1150
     DO  ig = 1,3
       ig1 = (ig-1)*3
       kg  = ig + 3
       DO  jg = 1,3
         jg1 = jg + ig1
         lg  = jg + 3
         gfour(ig,lg) = gi(jg1)
         gfour(kg,jg) = gi(jg1)
       END DO
     END DO
     1150 CONTINUE
     
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
       1170 DO  ir = 1,4
         IF (ik /= iorder(ir)) CYCLE
         iok = ir
         EXIT
       END DO
       1190 aa = shp(ioj)
       bb = shp(iok)
       
       DO  is = 1,3
         edgshr(is,ij) = (uev(is,ij)+anglei(ij)*uev(is,ii))*aa/  &
             (1.0D0-anglei(ij)*anglei(ij))  &
             + (uev(is,ij)+anglei(ik)*uev(is,ik))*bb/  &
             (1.0D0-anglei(ik)*anglei(ik))
       END DO
     END DO
     
!     TORSION-RELATED SHEAR CORRECTION FOR 4-NODE-
!     SET-UP OF EXPANDED SHEAR MATERIAL PROPERTY MATRICES (G OR D)
     
     csubx  = 20.0D0*reali/(tsi*(1.0D0-nunorx)*a*a)
     csuby  = 20.0D0*reali/(tsi*(1.0D0-nunory)*b*b)
     sfctr1 = csubb4*csubx
     sfctr2 = csubtx*csubx
     IF (.NOT.shrflx) GO TO 1220
     sfctr1 = sfctr1/(1.0D0+sfctr1)
     sfctr2 = sfctr2/(1.0D0+sfctr2)
     1220 CONTINUE
     sfctx1 = sfctr1 + sfctr2
     sfctx2 = sfctr1 - sfctr2
     sfctr1 = csubb4*csuby
     sfctr2 = csubty*csuby
     IF (.NOT.shrflx) GO TO 1230
     sfctr1 = sfctr1/(1.0D0+sfctr1)
     sfctr2 = sfctr2/(1.0D0+sfctr2)
     1230 CONTINUE
     sfcty1 = sfctr1 + sfctr2
     sfcty2 = sfctr1 - sfctr2
     
!     FILL IN THE EXPANDED MATERIAL PROPERTY MATRIX
     
     IF (norpth) GO TO 1240
     gfour( 7, 7) = 0.25D0*sfcty1*ts*gi(19)
     gfour( 8, 8) = 0.25D0*sfcty1*ts*gi(19)
     gfour( 8, 7) = 0.25D0*sfcty2*ts*gi(19)
     gfour( 7, 8) = gfour(8,7)
     gfour( 9, 9) = 0.25D0*sfctx1*ts*gi(22)
     gfour(10,10) = 0.25D0*sfctx1*ts*gi(22)
     gfour(10, 9) = 0.25D0*sfctx2*ts*gi(22)
     gfour( 9,10) = gfour(10,9)
     gfour( 7, 9) = DSQRT(tsmfx*tsmfy)*ts*gi(20)
     gfour( 9, 7) = gfour(7,9)
     GO TO 1250
     
     1240 dfour(4,4) = 0.25D0*sfcty1*ts*gi(19)
     dfour(5,5) = 0.25D0*sfcty1*ts*gi(19)
     dfour(5,4) = 0.25D0*sfcty2*ts*gi(19)
     dfour(4,5) = dfour(5,4)
     dfour(6,6) = 0.25D0*sfctx1*ts*gi(22)
     dfour(7,7) = 0.25D0*sfctx1*ts*gi(22)
     dfour(7,6) = 0.25D0*sfctx2*ts*gi(22)
     dfour(6,7) = dfour(7,6)
     dfour(4,6) = DSQRT(tsmfx*tsmfy)*ts*gi(20)
     dfour(6,4) = dfour(4,6)
     1250 CONTINUE

     DO  izta = 1,2
       zta  = ptint(izta)
       ibot = (izta-1)*nd2
       
       hzta = zta/2.0D0
       
!     TORSION-RELATED SHEAR CORRECTION FOR 4-NODE-
!     SET-UP OF POINTERS TO THE SAVED B-MATRIX
       
       iptx1 = ((ixsi-1)*2+ieta-1)*2*nd2 + ibot
       iptx2 = ((ixsi-1)*2+2-ieta)*2*nd2 + ibot
       ipty1 = ((ixsi-1)*2+ieta-1)*2*nd2 + ibot
       ipty2 = ((2-ixsi)*2+ieta-1)*2*nd2 + ibot
       
!     FILL IN THE 10X10 G-MATRIX IF MID4 IS PRESENT
       
       IF (.NOT.mbcoup) GO TO 1290
       DO  ig = 1,3
         ig1 = (ig-1)*3
         DO  jg = 1,3
           jg1 = jg  + ig1
           jg4 = jg1 + 27
           gfour(ig,jg) = gi(jg1)
         END DO
       END DO
       
       DO  ig = 4,6
         ig2 = (ig-2)*3
         DO  jg = 4,6
           jg2 = jg  + ig2
           jg4 = jg2 + 18
           gfour(ig,jg) = gi(jg2)*mominr
         END DO
       END DO
       
       DO  ig = 1,3
         ig4 = (ig+8)*3
         kg  = ig  + 3
         DO  jg = 1,3
           jg4 = jg  + ig4
           jg1 = jg4 - 27
           lg  = jg  + 3
           gfour(ig,lg) = -gi(jg4)*zta*6.0D0+gi(jg1)
           gfour(kg,jg) = -gi(jg4)*zta*6.0D0+gi(jg1)
         END DO
       END DO
       1290 CONTINUE
       
!     COMPUTE THE JACOBIAN AT THIS GAUSS POINT,
!     ITS INVERSE AND ITS DETERMINANT.
       
       CALL jacob2 (elid,shp,dshp,dgpth,egpdt,epnorm,jacob)
       IF (badjac) GO TO 1710
       
!     COMPUTE PSI TRANSPOSE X JACOBIAN INVERSE.
!     HERE IS THE PLACE WHERE THE INVERSE JACOBIAN IS FLAGED TO BE
!     TRANSPOSED BECAUSE OF OPPOSITE MATRIX LOADING CONVENTION
!     BETWEEN INVER AND GMMAT.
       
       CALL gmmatd (psitrn,3,3,0,jacob,3,3,1,phi)
       
!     CALL Q4BMGD TO GET B-MATRIX. SET THE ROW FLAG TO 3.
!     IT WILL RETURN THE FIRST 6 ROWS OF B-MATRIX.
       
       rowflg = 3
       CALL q4bmgd (dshp,dgpth,egpdt,epnorm,phi,bfour(1))
       
!     TORSION-RELATED SHEAR CORRECTION FOR 4-NODE -
!     SET-UP OF B-MATRIX AND TRIPLE MULTIPLY
       
       CALL trplmd (gfour,dfour,bfour,bmat1,xybmat,mattyp,jcored,detj)
     END DO
   END DO
 END DO
 
!     EQUALIZE THE OFF- DIAGNOAL TERMS TO GUARANTEE PERFECT SYMMETRIC
!     MATRIX IF NO DAMPING INVLOVED
 
 IF (gsube /= 0.0) GO TO 1330
 ij = jcored - 1
 ndofm1 = ndof - 1
 DO  ii = 1,ndofm1
   ip1 = ii + 1
   im1 = (ii-1)*ndof + ij
   DO  jj = ip1,ndof
     i = im1 + jj
     j = (jj-1)*ndof + ii + ij
     temp = (akgg(i) + akgg(j))*.5D0
     IF (DABS(temp) < 1.0D-17) temp = 0.0D0
     akgg(i) = temp
     akgg(j) = temp
   END DO
 END DO
 
!     END OF STIFFNESS LOOP
 
!     ADD NON-STRUCTURAL MASS
 
 1330 CONTINUE
 IF (mgg1 == 0) GO TO 1410
 IF (rho == 0.d0 .AND. nsm == 0.0) GO TO 1410
!     IF (CPMASS .GT. 0) GO TO 1410
 IF (nsm ==  0.0) GO TO 1410
 IF (vol == 0.d0 .OR. rhox == 0.d0) WRITE (nout,2060) sfm,elid,  &
     area,vol,rhox,mgg1,kgg1
 factor = (vol*rho+nsm*area)/(vol*rhox)
 DO  i = 1,nodesq
   xmass(i) = xmass(i)*factor
 END DO
 1410 CONTINUE
 
!     PICK UP THE GLOBAL TO BASIC TRANSFORMATIONS FROM THE CSTM.
 
 DO  i = 1,36
   trans(i)  = 0.0D0
 END DO
 DO  i = 1,nnode
   notran(i) = 0
   ipoint = 9*(i-1) + 1
   IF (igpdt(1,i) <= 0) GO TO 1420
   igpth(1) = igpdt(1,i)
   gpth (2) = bgpdt(2,i)
   gpth (3) = bgpdt(3,i)
   gpth (4) = bgpdt(4,i)
   
!     NOTE THAT THE 6X6 TRANSFORMATION WHICH WILL BE USED LATER
!     IN THE TRIPLE MULTIPLICATION TO TRANSFORM THE ELEMENT
!     STIFFNESS MATRIX FROM BASIC TO GLOBAL COORDINATES, IS BUILT
!     UPON THE 3X3 TRANSFORMATION FROM GLOBAL TO BASIC TBG-MATRIX.
!     THIS IS DUE TO THE DIFFERENCE IN TRANSFORMATION OF ARRAYS
!     AND MATRICES.
   
   CALL transd (gpth,tbg)
   CALL gmmatd (teb,3,3,0,tbg,3,3,0,trans(ipoint))
   CYCLE
   
   1420 IF (identt /= 1 .OR. offset /= 0.0D0) GO TO 1430
   notran(i) = 1
   CYCLE
   
   1430 DO  j = 1,9
     trans(ipoint+j-1) = teb(j)
   END DO
 END DO
 
 
!     HERE WE SHIP OUT THE STIFFNESS AND DAMPING MATRICES.
!     ----------------------------------------------------
 
 IF (kgg1 == 0) GO TO 1600
 
!     SET UP I-LOOP TO DUMP OUT BASIC TO GLOBAL TRANSFORMED, NODAL
!     PARTITIONED (6 D.O.F. PER NODE) COLUMNS OF THE ELEM. STIFFNESS.
 
!     THIS MEANS WE ARE SENDING TO EMGOUT 6 COLUMNS OF THE ELEMENT
!     STIFFNESS MATRIX AT A TIME.  EACH BUNCH OF 6 COLUMNS CORRESPOND
!     TO ONE PARTICULAR NODE OF THE ELEMENT.  FOR THE MASS MATRIX, WE
!     ONLY SEND 3 COLUMNS PER NODE TO EMGOUT SINCE THE OTHER 3 D.O.F.
!     ARE ZERO ANYWAY.  THE CODE WORD (DICT(4)) TELLS EMGOUT WHICH
!     COLUMNS ARE THE NON ZERO ONES THAT WE ARE SENDING. (SEE SECTION
!     6.8.3.5.1 OF THE PROGRAMMER MANUAL)
 
 
 dict(1) = estid
 dict(2) = 1
 dict(3) = ndof
 dict(4) = 63
 npart   = ndof*6
 DO  i = 1,nnode
   ibegin = 6*(i-1) + jcored - 1
   
!     DUMP AN UNTRANSFORMED NODAL COLUMN PARTITION.
   
   DO  j = 1,ndof
     kpoint = ndof*(j-1) + ibegin
     lpoint = 6*(j-1)
     DO  k = 1,6
       colstf(lpoint+k) = akgg(kpoint+k)
     END DO
   END DO
   IF (notran(i) == 1) GO TO 1515
   
!     THIS COLUMN PARTITION NEEDS TO BE TRANSFORMED TO GLOBAL
!     COORDINATES. (SEE PAGE 2.3-43 OF THE PROGRAMMER MANUAL)
   
!     LOAD THE 6X6 TRANSFORMATION
   
   CALL tldrd (offset,i,trans,trans1)
   
!     TRANSFORM THE NODAL COLUMN PARTITION.
   
   CALL gmmatd (colstf,ndof,6,0,trans1,6,6,0,coltmp)
   DO  ii = 1,npart
     colstf(ii) = coltmp(ii)
   END DO
   
!     NOW TRANSFORM THE ROWS OF THIS PARTITION.
   
   1515 DO  m = 1,nnode
     IF (notran(m) == 1) CYCLE
     mpoint = 36*(m-1) + 1
     
!     LOAD THE 6X6 TRANSFORMATION
     
     CALL tldrd (offset,m,trans,trans1)
     
!     TRANSFORM THE 6 ROWS FOR THIS SUBPARTITION
     
     CALL gmmatd (trans1,6,6,1,colstf(mpoint),6,6,0,coltmp)
     iipnt = mpoint - 1
     DO  ii = 1,36
       colstf(iipnt+ii) = coltmp(ii)
     END DO
   END DO
   
!     HERE WE MUST CHANGE FROM THE ROW LOADING CONVENTION
!     FOR GMMATD TO THE COLUMN LOADING CONVENTION FOR EMGOUT.
   
   DO  ii = 1,6
     ipoint = ndof*(ii-1)
     DO  jj = 1,ndof
       jpoint = 6*(jj-1)
       coltmp(ipoint+jj) = colstf(jpoint+ii)
     END DO
   END DO
   
!     DUMP THE TRANSFORMED NODAL COLUMN PARTITION
   
   ieoe = 0
   IF (i == nnode) ieoe = 1
   adamp = gsube
   
!     INTEGER 1 IN THE NEXT TO LAST FORMAL PARAMETER OF
!     EMGOUT MEANS WE ARE SENDING STIFFNESS DATA.
   
   CALL emgout (coltmp,coltmp,npart,ieoe,dict,1,iprec)
 END DO
 
 
!     HERE WE SHIP OUT THE MASS MATRIX.
!     ---------------------------------
 
 1600 IF (mgg1 == 0) GO TO 1710
 
 ndof  = nnode*3
 npart = ndof*3
 dict(3) = ndof
 dict(4) = 7
 adamp = 0.0D0
 
!     SET UP I-LOOP TO PROCESS AND DUMP THE NODAL COLUMN PARTITIONS.
 
 DO  i = 1,nnode
   DO  ijk = 1,npart
     amgg(jcored-1+ijk) = 0.0D0
   END DO
   
!     SET UP J-LOOP TO LOAD THE UNTRANSFORMED NODAL COLUMN PARTITION.
   
   DO  j = 1,nnode
     ipoint = 9*(j-1) + jcored
     jpoint = ipoint  + 4
     kpoint = ipoint  + 8
     ifrom  = nnode*(j-1) + i
     xmasso = xmass(ifrom)
     amgg(ipoint) = xmasso
     amgg(jpoint) = xmasso
     amgg(kpoint) = xmasso
   END DO
   IF (notran(i) == 1) GO TO 1670
   
!     THIS COLUMN PARTITION NEEDS TO BE TRANSFORMED
!     TO GLOBAL COORDINATES.
   
   DO  m = 1,nnode
     mpoint = 9*(m-1) + jcored
     CALL gmmatd (amgg(mpoint),3,3,0,trans(9*i-8),3,3,0,tmpmas)
     iicore = mpoint - 1
     DO  k = 1,9
       amgg(iicore+k) = tmpmas(k)
     END DO
   END DO
   
!     SET UP M-LOOP TO TRANSFORM THE NODAL ROW PARTITIONS
!     OF THIS NODAL COLUMN PARTITION.
   
   DO  m = 1,nnode
     mpoint = 9*(m-1) + jcored
     
!     TRANSFORM THE 3 ROWS FOR THIS SUBPARTITION.  THIS IS CORRECT
!     (3 ROWS).  REMEMBER THAT FOR THE MASS MATIIX FOR THIS ELEMENT
!     THERE ARE NO MASS MOMENT OF INERTIA TERMS.  THIS GIVES THREE
!     ROWS OF ZERO TERMS INTERSPERSED BETWEEN 3 ROWS OF NONZERO
!     TRANSLATIONAL MASS TERMS FOR EACH NODE.
     
     CALL gmmatd (trans(9*m-8),3,3,1,amgg(mpoint),3,3,0,tmpmas)
     iicore = mpoint - 1
     DO  k = 1,9
       amgg(iicore+k) = tmpmas(k)
     END DO
   END DO
   
!     HERE WE MUST CHANGE FROM THE ROW LOADING CONVENTION
!     FOR GMMATD TO THE COLUMN LOADING CONVENTION FOR EMGOUT.
   
   1670 DO  ii = 1,3
     ipoint = ndof*(ii-1)
     DO  jj = 1,ndof
       jpoint = 3*(jj-1) + jcored - 1
       coltmp(ipoint+jj) = amgg(jpoint+ii)
     END DO
   END DO
   
!     DUMP THIS TRANSFORMED MASS NODAL COLUMN PARTITION.
   
   ieoe = 0
   IF (i == nnode) ieoe = 1
   
!     INTEGER 2 IN THE NEXT TO LAST FORMAL PARAMETER OF
!     EMGOUT MEANS WE ARE SENDING MASS DATA.
   
   CALL emgout (coltmp,coltmp,npart,ieoe,dict,2,iprec)
 END DO
 GO TO 1710
 
 1700 j = 230
 GO TO 1720
 
 1710 CONTINUE
 RETURN
 
 1720 CALL mesage (30,j,nest)
 IF (l38 == 1) CALL mesage (-61,0,0)
 nogo = .true.
 GO TO 1710
 1730 CALL mesage (-30,234,nam)
 
 
!     HEAT FLOW OPTION STARTS HERE.
 
!     WE NEED TO RESTORE THE ORIGINAL ORDER OF SILS AND BGPDT DATA
 
 1790 j = 1
 DO  i = 1,20
   est(i+j) = SAVE(i)
   IF (i == 4) j = 24
 END DO
 
 inflag = 2
 cosmat = 1.0
 sinmat = 0.0
 matid  = nest(13)
 CALL hmat (elid)
 gi(1) = DBLE(kheat(1))
 gi(2) = DBLE(kheat(2))
 gi(3) = gi(2)
 gi(4) = DBLE(kheat(3))
 anis  = TYPE /= 4 .AND. TYPE /= -1
!     COMMENT-  ANIS = .FALSE. MEANS ISOTROPIC THERMAL CONDUCTIVITY.
 
 IF (anis) GO TO 400
 GO TO 1820
 1810 CONTINUE
 tem(3) = tem(4)
 tem(4) = tem(5)
 CALL gmmatd (tem,2,2,0,gi,2,2,0,GT)
 CALL gmmatd (GT,2,2,0,tem,2,2,1,gi)
 1820 CONTINUE
 DO  i = 1,16
   htcon(i) = 0.0D0
   htcap(i) = 0.0D0
 END DO
 DO  i = 5,8
   hsil(i)   = 0
   horder(i) = 0
 END DO
 
 DO  ixsi = 1,2
   xi = ptint(ixsi)
   DO  ieta = 1,2
     eta = ptint(ieta)
     DO  izta = 1,2
       zeta = ptint(izta)
       
       CALL termsd (nnode,dgpth,epnorm,egpdt,horder,hsil,bterms)
       dvol = determ
       
       DO  i = 1,4
         ecpt(i) = gi(i)*dvol
       END DO
       weitc = dvol*htcp
       
       ip = 1
       DO  i = 1,nnode
         idn = i + nnode
         htflx(ip+1) = ecpt(3)*bterms(i) + ecpt(4)*bterms(idn)
         htflx(ip  ) = ecpt(1)*bterms(i) + ecpt(2)*bterms(idn)
         ip = ip + 2
       END DO
       CALL gmmatd (bterms,2,nnode,-1,htflx,nnode,2,1,htcon)
       
     END DO
     IF (htcp == 0.0) CYCLE
     ip = 0
     DO  i = 1,nnode
       dheat = weitc*shp(i)
       DO  j = 1,nnode
         ip = ip + 1
         htcap(ip) = htcap(ip) + dheat*shp(j)
       END DO
     END DO
   END DO
 END DO
 dict(1) = estid
 dict(2) = 1
 dict(3) = nnode
 dict(4) = 1
 IF (htcp == 0.0) GO TO 1900
 adamp = 1.0
 CALL emgout (htcap,htcap,nodesq,1,dict,3,iprec)
 1900 CONTINUE
 adamp = 0.0
 CALL emgout (htcon,htcon,nodesq,1,dict,1,iprec)
 GO TO 1710
 
 2010 FORMAT (a23,', THE ELEMENT THICKNESS FOR QUAD4 EID =',i9,  &
                   ' IS NOT COMPLETELY DEFINED.')
 2030 FORMAT (a25,', RHO = ',1P,d12.4,' IS ILLEGAL FROM MATERIAL ID =',  &
               i9,' FOR QUAD4 EID =',i9)
 2060 FORMAT (a25,', ZERO VOLUME OR DENSITY FOR QUAD4 ELEMENT ID =',i9,  &
                  ', AREA,VOL,RHO =',3D12.3, /70X,'MGG1,KGG1 =',2I8)
                  
END SUBROUTINE quad4d
