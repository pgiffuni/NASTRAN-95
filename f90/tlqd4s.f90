SUBROUTINE tlqd4s
     
!     ELEMENT THERMAL LOAD GENERATOR FOR 4-NODE ISOPARAMETRIC
!     QUADRILATERAL SHELL ELEMENT (QUAD4)
!     (SINGLE PRECISION VERSION)
 
!     COMPLETELY RESTRUCTURED FOR COMPOSITES WITH THE FOLLOWING
!     LIMITATION -
!     1. FOR DIFFERENT GRID POINT TEMPERATURES AN AVERAGE
!        VALUE IS TAKEN.                       HEMANT  2/24/86
 
 
!                 EST LISTING
!     ---------------------------------------------------------
!      1          EID
!      2 THRU 5   SILS, GRIDS 1 THRU 4
!      6 THRU 9   T (MEMBRANE), GRIDS 1 THRU 4
!     10          THETA (MATERIAL)
!     11          TYPE FLAG FOR WORD 10
!     12          ZOFF  (OFFSET)
!     13          MATERIAL ID FOR MEMBRANE
!     14          T (MEMBRANE)
!     15          MATERIAL ID FOR BENDING
!     16          I FACTOR (BENDING)
!     17          MATERIAL ID FOR TRANSVERSE SHEAR
!     18          FACTOR FOR T(S)
!     19          NSM (NON-STRUCTURAL MASS)
!     20 THRU 21  Z1, Z2  (STRESS FIBRE DISTANCES)
!     22          MATERIAL ID FOR MEMBRANE-BENDING COUPLING
!     23          THETA (MATERIAL) FROM PSHELL CARD
!     24          TYPE FLAG FOR WORD 23
!     25          INTEGRATION ORDER
!     26          THETA (STRESS)
!     27          TYPE FLAG FOR WORD 26
!     28          ZOFF1 (OFFSET)  OVERRIDDEN BY EST(12)
!     29 THRU 44  CID,X,Y,Z - GRIDS 1 THRU 4
!     45          ELEMENT TEMPERATURE
 
 
 LOGICAL :: badjac,membrn,bendng,shrflx,mbcoup,norpth,  &
     tempp1,tempp2,pcmp,pcmp1,pcmp2,compos
 INTEGER :: intz(1),nout,nest(45),elid,sil(4),ksil(4),  &
     kcid(4),igpdt(4,4),flag,rowflg,necpt(4),  &
     mid(4),INDEX(3,3),indx(6,3),comps,nam(2),  &
     pcomp,pcomp1,pcomp2,sym,symmem,pid,pidloc
 REAL :: gpth(4),tgrid(4,4),gpnorm(4,4),bgpdt(4,4),  &
     matset,tmpthk(4),ecpt(4),egpdt(4,4),  &
     epnorm(4,4),bgpdm(3,4),alpham(6),tsub0,stemp,z
 REAL :: dgpth(4),thk,eps1,xi,eta,detj,hzta,psitrn(9),  &
     jacob(9),phi(9),mominr,coeff,reali,pi,twopi,  &
     raddeg,degrad,shp(4),dshp(8),tmpshp(4),  &
     dshptp(8),GT(9),g(6,6),gi(36),u(9),trans(36),  &
     ptint(2),gge(9),ggu(9),tbm(9),teb(9),tem(9),  &
     tub(9),tum(9),teu(9),tbg(9),ugpdm(3,4),cente(3),  &
     cent(3),x31,y31,x42,y42,aa,bb,cc,exi,exj,xm,ym,  &
     thetam,bmatrx(144),xybmat(96),alpha(6),alfam(3),  &
     alfab(3),talfam(3),talfab(3),alphad(6),pt(24),  &
     ptg(24),tbar,ttbar,tgrad,thrmom(3),g2i(9),g2(9),  &
     gtemps(6),epsubt(6),gepsbt(6),detu,detg2,determ
 REAL :: abbd(6,6),stiff(36),gprop(25),glay(9),glayt(9),  &
     gbar(9),galpha(3),alphal(3),alphae(3),mintr,  &
     tlam,theta,thetae,transl(9),tsubo,tmean,tempel,  &
     delta,deltat,zk,zk1,zref,zsubi,c,c2,s,s2,  &
     ftherm(6),epslnt(6),offset,const,uev,anglei, edgel,edgshr,unv
!WKBNB 11/93 SPR 93020
 REAL :: vd1(3), vd2(3), vkn(3), vks(3)  &
     ,                v12(3), v41(3), vp12(3),vis(3), vjs(3)
!WKBNE 11/93 SPR 93020
 COMMON /condas/  pi,twopi,raddeg,degrad
 COMMON /trimex/  est(45)
 COMMON /system/  buffer(100)
 COMMON /matin /  matid,inflag,eltemp
 COMMON /matout/  rmtout(25)
 COMMON /sgtmpd/  stemp(8)
!ZZ   COMMON /ZZSSB1/  Z(1)
 COMMON /zzzzzz/  z(20000)
 COMMON /BLANK /  nrowsp,iparam,comps
 COMMON /compst/  ipcmp,npcmp,ipcmp1,npcmp1,ipcmp2,npcmp2
 COMMON /q4dt  /  detj,hzta,psitrn,nnode,badjac,n1
 COMMON /terms /  membrn,bendng,shrflx,mbcoup,norpth
 COMMON /q4coms/  anglei(4),edgshr(3,4),edgel(4),unv(3,4),  &
     uev(3,4),rowflg,iorder(4)
 
 EQUIVALENCE     (z(1)     ,intz(1)), (igpdt(1,1),bgpdt(1,1))
 EQUIVALENCE     (est(1)   ,nest(1)), (bgpdt(1,1),est(29)   )
 EQUIVALENCE     (elth     ,est(14)), (gpth(1)   ,est(6)    )
 EQUIVALENCE     (sil(1)   ,nest(2)), (matset    ,rmtout(25))
 EQUIVALENCE     (zoff     ,est(12)), (zoff1     ,est(28)   )
 EQUIVALENCE     (necpt(1) ,ecpt(1)), (buffer(1) ,sysbuf    )
 EQUIVALENCE     (buffer(2),nout   ), (buffer(3) ,nogo      )
 EQUIVALENCE     (stemp(7) ,flag   ), (alfab(1)  ,alpha(4)  )
 EQUIVALENCE     (alfam(1) ,alpha(1))
 
 DATA eps1     / 1.0E-7 /
 DATA pcomp    / 0 /
 DATA pcomp1   / 1 /
 DATA pcomp2   / 2 /
 DATA sym      / 1 /
 DATA mem      / 2 /
 DATA symmem   / 3 /
 DATA const    / 0.57735026918962 /
 DATA nam      / 4HTLQD,4H4S      /
 
!-----ZERO THE VARIOUS ALPHA ARRAYS
 
 DO  i =1,6
   alpham(i) = 0.0
   alpha(i)  = 0.0
   alphad(i) = 0.0
 END DO
 DO  i =1,3
   talfam(i) = 0.0
   talfab(i) = 0.0
 END DO
 
 elid=nest(1)
 ltypfl=1
 offset=zoff
 IF (zoff == 0.0) offset=zoff1
 
!     TEST FOR COMPOSITE ELEMENT
 
 compos = .false.
 
 pid = nest(13) - 100000000
 compos = comps == -1 .AND. pid > 0
 
!-----CHECK FOR THE TYPE OF TEMPERATURE DATA
!     NOTES-  1- TYPE TEMPP1 ALSO INCLUDES TYPE TEMPP3
!             2- IF NO TEMPPI CARDS, GRID POINT TEMPERATURES
!                ONLY ARE PRESENT
 
 tempp1 = flag == 13
 tempp2 = flag == 2
 
 n1=4
 nnode=4
 ndof=nnode*6
 nd2=ndof*2
 nd3=ndof*3
 nd4=ndof*4
 nd5=ndof*5
 
!     FILL IN ARRAY GGU WITH THE COORDINATES OF GRID POINTS
!     1, 2 AND 4. THIS ARRAY WILL BE USED LATER TO DEFINE
!     THE USER COORDINATE SYSTEM WHILE CALCULATING
!     TRANSFORMATIONS INVOLVING THIS COORDINATE SYSTEM.
 
 DO  i=1,3
   ii=(i-1)*3
   ij=i
   IF (ij == 3) ij=4
   DO  j=1,3
     jj=j+1
     ggu(ii+j)=bgpdt(jj,ij)
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
 
!     STORE INCOMING BGPDT FOR ELEMENT C.S. CALCULATION
 
 DO  i=1,3
   i1=i+1
   DO  j=1,4
     bgpdm(i,j)=bgpdt(i1,j)
   END DO
 END DO
 
!     TRANSFORM BGPDM FROM BASIC TO USER C.S.
 
 DO  i=1,3
   ip=(i-1)*3
   DO  j=1,4
     ugpdm(i,j)=0.0
     DO  k=1,3
       kk=ip+k
       ugpdm(i,j)=ugpdm(i,j)+tub(kk)*((bgpdm(k,j))-ggu(k))
     END DO
   END DO
 END DO
 
 
!     THE ORIGIN OF THE ELEMENT C.S. IS IN THE MIDDLE OF THE ELEMENT
 
 DO  j=1,3
   cent(j)=0.0
   DO  i=1,4
     cent(j)=cent(j)+ugpdm(j,i)/nnode
   END DO
 END DO
 
!     STORE THE CORNER NODE DIFF. IN THE USER C.S.
 
 x31=ugpdm(1,3)-ugpdm(1,1)
 y31=ugpdm(2,3)-ugpdm(2,1)
 x42=ugpdm(1,4)-ugpdm(1,2)
 y42=ugpdm(2,4)-ugpdm(2,2)
 aa=SQRT(x31*x31+y31*y31)
 bb=SQRT(x42*x42+y42*y42)
 
!     NORMALIZE XIJ'S
 
 x31=x31/aa
 y31=y31/aa
 x42=x42/bb
 y42=y42/bb
 exi=x31-x42
 exj=y31-y42
 
!     STORE GGE ARRAY, THE OFFSET BETWEEN ELEMENT C.S. AND USER C.S.
 
 gge(1)=cent(1)
 gge(2)=cent(2)
 gge(3)=cent(3)
 
 gge(4)=gge(1)+exi
 gge(5)=gge(2)+exj
 gge(6)=gge(3)
 
 gge(7)=gge(1)-exj
 gge(8)=gge(2)+exi
 gge(9)=gge(3)
 
 
!     THE ARRAY IORDER STORES THE ELEMENT NODE ID IN
!     INCREASING SIL ORDER.
 
!     IORDER(1) = NODE WITH LOWEST  SIL NUMBER
!     IORDER(4) = NODE WITH HIGHEST SIL NUMBER
 
!     ELEMENT NODE NUMBER IS THE INTEGER FROM THE NODE
!     LIST  G1,G2,G3,G4 .  THAT IS, THE 'I' PART
!     OF THE 'GI' AS THEY ARE LISTED ON THE CONNECTIVITY
!     BULK DATA CARD DESCRIPTION.
 
 
 DO  i=1,4
   iorder(i)=0
   ksil(i)=sil(i)
 END DO
 
 DO  i=1,4
   itemp=1
   isil=ksil(1)
   DO  j=2,4
     IF (isil <= ksil(j)) CYCLE
     itemp=j
     isil=ksil(j)
   END DO
   iorder(i)=itemp
   ksil(itemp)=99999999
 END DO
 
!     ADJUST EST DATA
 
 
!     USE THE POINTERS IN IORDER TO COMPLETELY REORDER THE
!     GEOMETRY DATA INTO INCREASING SIL ORDER.
!     DON'T WORRY!! IORDER ALSO KEEPS TRACK OF WHICH SHAPE
!     FUNCTIONS GO WITH WHICH GEOMETRIC PARAMETERS!
 
 
 DO  i=1,4
   ksil(i)=sil(i)
   tmpthk(i)=gpth(i)
   kcid(i)=igpdt(1,i)
   DO  j=2,4
     tgrid(j,i)=bgpdt(j,i)
   END DO
 END DO
 DO  i=1,4
   ipoint=iorder(i)
   sil(i)=ksil(ipoint)
   gpth(i)=tmpthk(ipoint)
   igpdt(1,i)=kcid(ipoint)
   DO  j=2,4
     bgpdt(j,i)=tgrid(j,ipoint)
   END DO
 END DO
 
!-----SORT THE GRID POINT TEMPERATURES (IN STEMP(1-4))
!     IF PRESENT AND MAKE REAL          THE OTHER
!     KINDS OF TEMPERATURE DATA IF TEMPPI CARDS PRESENT
 
 IF (tempp1 .OR. tempp2) GO TO 150
 
 tempel = 0.0
 DO  i =1,4
   ipnt = iorder(i)
   gtemps(i) = stemp(ipnt)
   tempel = tempel + 0.25 * gtemps(i)
 END DO
 GO TO 170
 
 150 IF (tempp2) GO TO 160
 
 tbar  = stemp(1)
 tgrad = stemp(2)
 GO TO 170
 
 160 tbar = stemp(1)
 thrmom(1) = stemp(2)
 thrmom(2) = stemp(3)
 thrmom(3) = stemp(4)
 170 CONTINUE
 
!     COMPUTE NODE NORMALS
 
 CALL q4nrms (bgpdt,gpnorm,iorder,iflag)
 IF (iflag == 0) GO TO 180
 j = -230
 GO TO 1580
 
!     DETERMINE NODAL THICKNESSES
 
 180 DO  i=1,nnode
   IF (gpth(i) == 0.0) gpth(i)=elth
   IF (gpth(i) > 0.0) GO TO 190
   WRITE (nout,1700) elid
   nogo=1
   GO TO 1600
   190 dgpth(i)=gpth(i)
 END DO
 
 mominr=0.0
 IF (nest(15) /= 0) mominr=est(16)
 
 
!     THE COORDINATES OF THE ELEMENT GRID POINTS HAVE TO BE
!     TRANSFORMED FROM THE BASIC C.S. TO THE ELEMENT C.S.
 
 
 CALL betrns (teu,gge,0,elid)
 CALL gmmats (teu,3,3,0,tub,3,3,0,teb)
 CALL gmmats (tub,3,3,1,cent,3,1,0,cente)
 
 ip = -3
 DO  ii=2,4
   ip=ip+3
   DO  j=1,nnode
     epnorm(ii,j)=0.0
     egpdt(ii,j)=0.0
     DO  k=1,3
       kk=ip+k
       k1=k+1
       cc=(bgpdt(k1,j))-ggu(k)-cente(k)
       epnorm(ii,j)=epnorm(ii,j)+teb(kk)*gpnorm(k1,j)
       egpdt( ii,j)=egpdt(ii,j)+(teb(kk)*cc)
     END DO
   END DO
 END DO
!WKBNB 11/93 SPR93020
 DO  j = 1, 4
   egpdt(4,j) = cent(3)
 END DO
!WKBNE 11/93 SPR93020
 
!     BEGIN INITIALIZING MATERIAL VARIABLES
 
!     SET INFLAG = 12 SO THAT SUBROUTINE MAT WILL SEARCH FOR -
!     ISOTROPIC MATERIAL PROPERTIES AMONG THE MAT1 CARDS,
!     ORTHOTROPIC MATERIAL PROPERTIES AMONG THE MAT8 CARDS, AND
!     ANISOTROPIC MATERIAL PROPERTIES AMONG THE MAT2 CARDS.
 
 inflag=12
 eltemp= est(45)
 mid(1)=nest(13)
 mid(2)=nest(15)
 mid(3)=0
 mid(4)=nest(22)
 membrn=mid(1) > 0
 bendng=mid(2) > 0 .AND. mominr > 0.0
 shrflx=mid(3) > 0
 mbcoup=mid(4) > 0
 norpth=.false.
 
!     SET THE INTEGRATION POINTS
 
 ptint(1) = -const
 ptint(2) =  const
 
!     IN PLANE SHEAR REDUCTION
 
 xi =0.0
 eta=0.0
 kpt=1
 kpt1=nd2
 
 CALL q4shps (xi,eta,shp,dshp)
 
!     SORT THE SHAPE FUNCTIONS AND THEIR DERIVATIVES INTO SIL ORDER.
 
 DO  i=1,4
   tmpshp(i  )=shp(i)
   dshptp(i  )=dshp(i)
   dshptp(i+4)=dshp(i+4)
 END DO
 DO  i=1,4
   kk=iorder(i)
   shp (i  )=tmpshp(kk)
   dshp(i  )=dshptp(kk)
   dshp(i+4)=dshptp(kk+4)
 END DO
 
 DO  izta=1,2
   zta =ptint(izta)
   hzta=zta/2.0
   CALL jacobs (elid,shp,dshp,dgpth,egpdt,epnorm,jacob)
   IF (badjac) GO TO 1600
   
   CALL gmmats (psitrn,3,3,0,jacob,3,3,1,phi)
   
!     CALL Q4BMGS TO GET B MATRIX
!     SET THE ROW FLAG TO 2. IT WILL SAVE THE 3RD ROW OF B AT
!     THE TWO INTEGRATION POINTS.
   
   rowflg = 2
   CALL q4bmgs (dshp,dgpth,egpdt,epnorm,phi,xybmat(kpt))
   kpt=kpt+kpt1
 END DO
 
!     SET THE ARRAY OF LENGTH 4 TO BE USED IN CALLING TRANSS.
!     NOTE THAT THE FIRST WORD IS THE COORDINATE SYSTEM ID WHICH
!     WILL BE SET IN POSITION LATER.
 
 DO  iec=2,4
   ecpt(iec)=0.0
 END DO
 
!     FETCH MATERIAL PROPERTIES
 
!     EACH MATERIAL PROPERTY MATRIX G HAS TO BE TRANSFORMED FROM
!     THE MATERIAL COORDINATE SYSTEM TO THE ELEMENT COORDINATE
!     SYSTEM. THESE STEPS ARE TO BE FOLLOWED -
 
!     1- IF MCSID HAS BEEN SPECIFIED, SUBROUTINE TRANSS IS CALLED
!        TO CALCULATE TBM MATRIX (MATERIAL TO BASIC TRANSFORMATION).
!        THIS WILL BE FOLLOWED BY A CALL TO SUBROUTINE BETRNS
!        TO CALCULATE TEB MATRIX (BASIC TO ELEMENT TRANSFORMATION).
!        TBM IS THEN PREMULTIPLIED BY TEB TO OBTAIN TEM MATRIX.
!        THEN USING THE PROJECTION OF X-AXIS, AN ANGLE IS CALCULATED
!        UPON WHICH STEP 2 IS TAKEN.
 
!     2- IF THETAM HAS BEEN SPECIFIED, SUBROUTINE ANGTRS IS CALLED
!        TO CALCULATE TEM MATRIX (MATERIAL TO ELEMENT TRANSFORMATION).
 
!                        T
!     3-          G  =  U   G   U
!                  E         M
 
 
 IF (nest(11) == 0) GO TO 390
 mcsid=nest(10)
 
!     CALCULATE TEM USING MCSID
 
 340 IF (mcsid > 0) GO TO 360
 DO  i=1,9
   tem(i)=teb(i)
 END DO
 GO TO 370
 360 necpt(1)=mcsid
 CALL transs (ecpt,tbm)
 
!     MULTIPLY TEB AND TBM
 
 CALL gmmats (teb,3,3,0,tbm,3,3,0,tem)
 
!     CALCULATE THETAM FROM THE PROJECTION OF THE X-AXIS OF THE
!     MATERIAL C.S. ON TO THE XY PLANE OF THE ELEMENT C.S.
 
 370 imt=-1
 xm=tem(1)
 ym=tem(4)
 IF (ABS(xm) <= eps1) imt=imt+1
 IF (ABS(ym) <= eps1) imt=imt+2
 IF (imt < 2) GO TO 380
 nest(2) = mcsid
 j = -231
 GO TO 1580
 380 thetam= ATAN2(ym,xm)
 GO TO 400
 
!     CALCULATE TEM USING THETAM
 
 390 thetam = (est(10))*degrad
 IF (thetam == 0.0) GO TO 410
 400 CALL angtrs (thetam,1,tum)
 CALL gmmats (teu,3,3,0,tum,3,3,0,tem)
 GO TO 430
 
!     DEFAULT IS CHOSEN, LOOK FOR VALUES OF MCSID AND/OR THETAM
!     ON THE PSHELL CARD.
 
 410 IF (nest(24) == 0) GO TO 420
 mcsid=nest(23)
 GO TO 340
 
 420 thetam = (est(23))*degrad
 GO TO 400
 
 430 CONTINUE
 
!     BEGIN THE LOOP TO FETCH PROPERTIES FOR EACH MATERIAL ID
 
 m=0
 500 m=m+1
 IF (m > 4) GO TO 690
 matid=mid(m)
 IF (matid == 0) GO TO 500
 
 IF (m-1 < 0) THEN
   GO TO   530
 ELSE IF (m-1 == 0) THEN
   GO TO   520
 END IF
 510 IF (matid == mid(m-1)) GO TO 530
 520 CALL mat (elid)
 530 CONTINUE
 
 tsub0 = rmtout(11)
 IF (matset == 8.0) tsub0 = rmtout(10)
 
 coeff=1.0
 lpoint=(m-1)*9+1
 
 CALL q4gmgs (m,coeff,gi(lpoint))
 
 IF (thetam == 0.0) GO TO 550
 
 u(1)=tem(1)*tem(1)
 u(2)=tem(4)*tem(4)
 u(3)=tem(1)*tem(4)
 u(4)=tem(2)*tem(2)
 u(5)=tem(5)*tem(5)
 u(6)=tem(2)*tem(5)
 u(7)=tem(1)*tem(2)*2.0
 u(8)=tem(4)*tem(5)*2.0
 u(9)=tem(1)*tem(5)+tem(2)*tem(4)
 l=3
 
 CALL gmmats (u(1),l,l,1,gi(lpoint),l,l,0,GT(1))
 CALL gmmats (GT(1),l,l,0,u(1),l,l,0,gi(lpoint))
 
 550 CONTINUE
 
 IF (compos) GO TO 500
 
!-----TRANSFORM THERMAL EXPANSION COEFFICIENTS AND STORE THEM IN ALPHA
 
 IF (m > 2) GO TO 500
 morb = (m-1)*3
 IF (matset == 2.0) GO TO 610
 IF (matset == 8.0) GO TO 630
 
!     MAT1
 
 DO  imat=1,2
   alpham(imat+morb)=rmtout(8)
 END DO
 alpham(3+morb) = 0.0
 GO TO 640
 
!     MAT2
 
 610 DO  imat=1,3
   alpham(imat+morb)=rmtout(7+imat)
 END DO
 GO TO 640
 
!     MAT8
 
 630 alpham(morb+1)=rmtout(8)
 alpham(morb+2)=rmtout(9)
 alpham(morb+3)=0.0
 
!-----SKIP THE TRANSFORMATION OF ALPHAM IF MATSET = 1.0
!     OR THETAM = 0.0
 
 640 CONTINUE
 
 IF (matset == 1.0) GO TO 650
 IF (thetam /= 0.0) GO TO 670
 
 650 DO  ig = 1,3
   alpha(ig+morb) = alpham(ig+morb)
 END DO
 GO TO 500
 
!-----THE ALPHAS NEED TO BE PREMULTIPLIED BY U INVERSE.
!     INCREMENT MORB BY 1 TO INDICATE WHERE TO FILL THE
!     ARRAYS, AND PUT THE SINGLE PREC. ARRAY OF ALPHAM
!     INTO THE DOUBLE PREC. ARRAY OF ALPHAD FOR THE CALL
!     TO GMMATS.
 
 670 morb = morb + 1
 DO  i =1,6
   alphad(i) = alpham(i)
 END DO
 CALL invers (3,u,3,bdum,0,detu,isngu,INDEX)
 CALL gmmats (u,3,3,0,alphad(morb),3,1,0,alpha(morb))
 GO TO 500
 
 690 IF (.NOT.compos) GO TO  1070
 
!****
!      IF LAMINATED COMPOSITE ELEMENT, DETERMINE THE THERMAL
!      STRAIN VECTOR DUE TO THE APPLIED THERMAL LOADING.
!      NOTE THE FOLLOWING -
!         1. DIFFERENT GRID POINT TEMPERATURES ARE NOT SUPPORTED
 
!****
!     LOCATE PID BY CARRYING OUT A SEQUENTIAL SEARCH
!     OF THE PCOMPS DATA BLOCK, AND ALSO DETERMINE
!     THE TYPE OF 'PCOMP' BULK DATA ENTRY.
!****
 
!****
!     POINTER DESCRIPITION
!     --------------------
!     IPCMP  - LOCATION OF START OF PCOMP DATA IN CORE
!     NPCMP  - NUMBER OF WORDS OF PCOMP DATA
!     IPCMP1 - LOCATION OF START OF PCOMP1 DATA IN CORE
!     NPCMP1 - NUMBER OF WORDS OF PCOMP1 DATA
!     IPCMP2 - LOCATION OF START OF PCOMP2 DATA IN CORE
!     NPCMP2 - NUMBER OF WORDS OF PCOMP2 DATA
 
!     ITYPE  - TYPE OF PCOMP BULK DATA ENTRY
 
 
!     LAMOPT - LAMINATION GENERATION OPTION
!            = SYM  (SYMMETRIC)
!            = MEM  (MEMBRANE)
!            = SYMMEM  (SYMMETRIC-MEMBRANE)
 
 
!****
 
 
!**** SET POINTER LPCOMP
 lpcomp = ipcmp + npcmp + npcmp1 + npcmp2
 
!**** SET POINTERS
 itype = -1
 
 pcmp  = .false.
 pcmp1 = .false.
 pcmp2 = .false.
 
 pcmp   = npcmp  > 0
 pcmp1  = npcmp1 > 0
 pcmp2  = npcmp2 > 0
 
!**** CHECK IF NO 'PCOMP' DATA HAS BEEN READ INTO CORE
 
 IF (pcmp .OR. pcmp1 .OR. pcmp2) GO TO 700
 j = -229
 GO TO 1580
 
!**** SEARCH FOR PID IN PCOMP DATA
 
 700 IF (.NOT.pcmp) GO TO 750
 
 ip = ipcmp
 IF (intz(ip) == pid) GO TO 740
 ipc11 = ipcmp1 - 1
 DO  ip = ipcmp,ipc11
   IF (intz(ip) == -1 .AND. ip < (ipcmp1-1)) GO TO 710
   CYCLE
   710 IF (intz(ip+1) == pid) GO TO 730
 END DO
 GO TO 750
 
 730 ip = ip+1
 740 itype = pcomp
 GO TO 860
 
!**** SEARCH FOR PID IN PCOMP1 DATA
 
 750 IF (.NOT.pcmp1) GO TO 800
 ip = ipcmp1
 IF (intz(ip) == pid) GO TO 790
 ipc21 = ipcmp2 - 1
 DO  ip = ipcmp1,ipc21
   IF (intz(ip) == -1 .AND. ip < (ipcmp2-1)) GO TO 760
   CYCLE
   760 IF (intz(ip+1) == pid) GO TO 780
 END DO
 GO TO 800
 
 780 ip = ip+1
 790 itype = pcomp1
 GO TO 860
 
!**** SEARCH FOR PID IN PCOMP2 DATA
 
 800 ip = ipcmp2
 IF (intz(ip) == pid) GO TO 840
 lpc11 = lpcomp - 1
 DO  ip = ipcmp2,lpc11
   IF (intz(ip) == -1 .AND. ip < (lpcomp-1)) GO TO 810
   CYCLE
   810 IF (intz(ip+1) == pid) GO TO 830
 END DO
 GO TO 850
 
 830 ip = ip+1
 840 itype = pcomp2
 GO TO 860
 
 
!**** CHECK IF PID HAS NOT BEEN LOCATED
 
 850 IF (itype /= -1) GO TO 860
 j = -229
 GO TO 1580
 
!**** LOCATION OF PID
 
 860 pidloc = ip
 lamopt = intz(pidloc+8)
 
 
!**** DETERMINE INTRINSIC LAMINATE PROPERTIES
 
!     LAMINATE THICKNESS
 
 tlam = elth
 
!**** LAMINATE EXTENSIONAL, BENDING AND MEMBRANE-BENDING MATRICES
 
 DO  ll = 1,6
   DO  mm = 1,6
     abbd(ll,mm) = 0.0
   END DO
 END DO
 
!     EXTENSIONAL
 
 matid = mid(1)
 CALL mat (elid)
 
 CALL lprops (gprop)
 
 DO  ll = 1,3
   DO  mm = 1,3
     ii = mm + 3*(ll-1)
     abbd(ll,mm) = gprop(ii)*tlam
   END DO
 END DO
 
 IF (lamopt == mem .OR. lamopt == symmem) GO TO 910
 
!     BENDING
 
 matid = mid(2)
 CALL mat (elid)
 
 CALL lprops (gprop)
 
!**** MOMENT OF INERTIA OF LAMINATE
 mintr = (tlam**3)/12.0
 
 DO  ll = 1,3
   DO  mm = 1,3
     ii = mm + 3*(ll-1)
     abbd(ll+3,mm+3) = gprop(ii)*mintr
   END DO
 END DO
 
 IF (lamopt == sym) GO  TO 910
 
!     MEMBRANE-BENDING
 
 matid = mid(4)
 CALL mat (elid)
 
 CALL lprops (gprop)
 
 DO  ll = 1,3
   DO  mm = 1,3
     ii = mm + 3*(ll-1)
     abbd(ll,mm+3) = gprop(ii)*tlam*tlam
     abbd(ll+3,mm) = gprop(ii)*tlam*tlam
   END DO
 END DO
 
 910 CONTINUE
 
!**** REFERENCE SURFACE
 zref = -tlam/2.0
 
!**** NUMBER OF LAYERS
 nlay = intz(pidloc + 1)
 
!**** SET POINTER
 IF (itype == pcomp ) ipoint = (pidloc + 8 + 4*nlay)
 IF (itype == pcomp1) ipoint = (pidloc + 8 + nlay)
 IF (itype == pcomp2) ipoint = (pidloc + 8 + 2*nlay)
 
!****
!     ALLOW FOR THE ORIENTATION OF THE MATERIAL AXIS FROM
!     THE ELEMENT AXIS
!****
 
 thetae = ATAN(tem(2)/tem(1))
 thetae = thetae*degrad
 
 
!**** LAMINATE REFERENCE (OR LAMINATION) TEMPERATURE
 tsubo = z(ipoint+24)
 
 IF (tempp1 .OR. tempp2) GO TO 920
 tmean = tempel
 GO TO 930
 
 920 tmean = stemp(1)
 
 930 delta = tmean - tsubo
 
 DO  ll = 1,6
   ftherm(ll) = 0.0
 END DO
 
!**** ALLOW FOR APPLIED THERMAL MOMENTS
 
 IF (.NOT.tempp2) GO TO 960
 
 DO  ll = 1,3
   ftherm(ll+3) = thrmom(ll)
 END DO
 
 960 CONTINUE
 
 
!     L O O P  O V E R  N L A Y
 
 DO  k = 1,nlay
   
   zk1 = zk
   IF (k == 1) zk1 = zref
   IF (itype == pcomp ) zk = zk1 + z(pidloc + 6 + 4*k)
   IF (itype == pcomp1) zk = zk1 + z(pidloc + 7)
   IF (itype == pcomp2) zk = zk1 + z(pidloc + 7 + 2*k)
   
   zsubi = (zk+zk1)/2.0
   
!**** LAYER THICKNESS
   ti = zk - zk1
   
!**** LAYER ORIENTATION
   IF (itype == pcomp ) theta = z(pidloc + 7 + 4*k)
   IF (itype == pcomp1) theta = z(pidloc + 8 + k)
   IF (itype == pcomp2) theta = z(pidloc + 8 + 2*k)
   
   
   theta = theta * degrad
   
   IF (thetae > 0.0) theta = theta + thetae
   
   c   = COS(theta)
   c2  = c*c
   s   = SIN(theta)
   s2  = s*s
   
   transl(1)  = c2
   transl(2)  = s2
   transl(3)  = c*s
   transl(4)  = s2
   transl(5)  = c2
   transl(6)  =-c*s
   transl(7)  =-2.0*c*s
   transl(8)  = 2.0*c*s
   transl(9)  = c2-s2
   
!**** CALCULATE GBAR = TRANST X GLAY X TRANS
   
   DO  ir = 1,9
     glay(ir) = z(ipoint+ir)
   END DO
   
   CALL gmmats (glay(1),3,3,0,transl(1),3,3,0,glayt(1))
   CALL gmmats (transl(1),3,3,1,glayt(1),3,3,0,gbar(1))
   
!**** CALCULATE ALPHAE = TRANSL X ALPHA
   
   
!     MODIFY TRANSL FOR TRANSFORMATIONS OF ALPHAS
   
   transl(3) = -transl(3)
   transl(6) = -transl(6)
   transl(7) = -transl(7)
   transl(8) = -transl(8)
   
   DO  ir = 1,3
     alphal(ir) = z(ipoint+13+ir)
   END DO
   
   CALL gmmats (transl(1),3,3,0,alphal(1),3,1,0,alphae(1))
   
   
!**** CALCULATE LAMINATE OPERATING TEMPERATURE (ALLOWING FOR
!     TEMPERATURE GRADIENT IF APPLIED)
   
   deltat = delta
   IF (tempp1) deltat = delta + zsubi*tgrad
   
!**** CALCULATE THERMAL FORCES AND MOMENTS
   
   CALL gmmats (gbar(1),3,3,0,alphae(1),3,1,0,galpha(1))
   
   DO  ir = 1,3
     ftherm(ir) = ftherm(ir) + galpha(ir)*deltat*(zk - zk1)
     IF (lamopt == mem .OR. lamopt == symmem) CYCLE
     ftherm(ir+3) = ftherm(ir+3) - galpha(ir)*deltat*((zk**2)-(zk1**2))/2.0
   END DO
   
   IF (lamopt /= sym .AND. lamopt /= symmem) GO TO 1040
   
!**** CALCULATE CONTRIBUTION FROM SYMMETRIC LAYERS
   
   deltat = delta
   IF (tempp1) deltat = delta - zsubi*tgrad
   
   DO  ir = 1,3
     ftherm(ir) = ftherm(ir) + galpha(ir)*deltat*(zk-zk1)
     IF (lamopt == symmem) CYCLE
     ftherm(ir+3) = ftherm(ir+3) - galpha(ir)*deltat*((zk1**2)-(zk**2))/2.0
   END DO
   
   1040 IF (itype == pcomp) ipoint = ipoint + 27
   
 END DO
 
 
!****
!     COMPUTE THERMAL STRAIN VECTOR
!****
!                 -1
!     EPSLN = ABBD   X FTHERM
 
 CALL invers (6,abbd,6,dum,0,determ,ising,indx)
 
 DO  ll = 1,6
   DO  mm = 1,6
     nn = mm + 6*(ll-1)
     stiff(nn) = abbd(ll,mm)
   END DO
 END DO
 
 CALL gmmats (stiff(1),6,6,0,ftherm(1),6,1,0,epslnt(1))
 
 1070 CONTINUE
 
 
!-----INITIALIZE NECESSARY ARRAYS BEFORE STARTING THE
!     DOUBLE INTEGRATION LOOP
 
 DO  i =1,9
   g2(i) = 0.0
 END DO
 DO  i =1,6
   epsubt(i) = 0.0
 END DO
 DO  i =1,ndof
   pt(i)  = 0.0
   ptg(i) = 0.0
 END DO
 
!     FILL IN THE 6X6 MATERIAL PROPERTY MATRIX G
 
 DO  ig=1,6
   DO  jg=1,6
     g(ig,jg)=0.0
   END DO
 END DO
 
 IF (.NOT.membrn) GO TO 1150
 DO  ig=1,3
   ig1=(ig-1)*3
   DO  jg=1,3
     jg1=jg+ig1
     g(ig,jg)=gi(jg1)
   END DO
 END DO
 
 1150 IF (.NOT.bendng) GO TO 1180
 i = 0
 DO  ig=4,6
   ig2=(ig-2)*3
   DO  jg=4,6
     jg2=jg+ig2
     g(ig,jg)=gi(jg2)*mominr
     
!     SAVE THE G-MATRIX FOR BENDING IN G2
     
     i = i + 1
     g2(i) = g(ig,jg)
   END DO
 END DO
 
 IF (.NOT.membrn) GO TO 1180
 IF (mbcoup) GO TO 1180
 DO  ig=1,3
   ig1=(ig-1)*3
   kg=ig+3
   DO  jg=1,3
     jg1=jg+ig1
     lg=jg+3
     g(ig,lg)=gi(jg1)
     g(kg,jg)=gi(jg1)
   END DO
 END DO
 1180 CONTINUE
 
!****
!     HERE BEGINS THE DOUBLE LOOP ON STATEMENT 1470 TO
!     GAUSS INTEGRATE FOR THE ELEMENT STIFFNESS MATRIX.
!****
 
 DO  ixsi=1,2
   xi=ptint(ixsi)
   
   DO  ieta=1,2
     eta=ptint(ieta)
     
     CALL q4shps (xi,eta,shp,dshp)
     
!     SORT THE SHAPE FUNCTIONS AND THEIR DERIVATIVES INTO SIL ORDER.
     
     DO  i=1,4
       tmpshp(i  )=shp(i)
       dshptp(i  )=dshp(i)
       dshptp(i+4)=dshp(i+4)
     END DO
     DO  i=1,4
       kk=iorder(i)
       shp (i  )=tmpshp(kk)
       dshp(i  )=dshptp(kk)
       dshp(i+4)=dshptp(kk+4)
     END DO
     
!     CALCULATE THE ELEMENT THICKNESS AT THIS POINT
     
     thk=0.0
     DO  i=1,nnode
       thk=thk+dgpth(i)*shp(i)
     END DO
     reali=thk*thk*thk/12.0
     
!-----CALCULATE T-BAR FOR THIS INTEGRATION POINT
!     SKIP OVER IF TEMPPI CARDS ARE PRESENT
!     THEN CALCULATE ALPHA*T FOR EACH CASE
     
     IF (compos) GO TO 1370
     
     IF (tempp1 .OR. tempp2) GO TO 1310
     tbar = 0.0
     DO  i =1,nnode
       tbar = tbar + shp(i) * gtemps(i)
     END DO
     1310 CONTINUE
     
     ttbar = tbar - tsub0
     
     IF (.NOT.membrn) GO TO 1330
     DO  i =1,3
       talfam(i) = ttbar * alfam(i)
     END DO
     
     1330 IF (.NOT.bendng) GO TO 1370
     IF (.NOT.tempp1 .AND. .NOT.tempp2) GO TO 1370
     IF (tempp2) GO TO 1350
     DO  i =1,3
       talfab(i) = -tgrad * alfab(i)
     END DO
     GO TO 1370
     
     1350 CONTINUE
     DO  ig2=1,9
       g2i(ig2) = g2(ig2)*reali
     END DO
     CALL invers (3,g2i,3,gdum,0,detg2,isngg2,INDEX)
     CALL gmmats (g2i,3,3,0,thrmom,3,1,0,talfab)
     1370 CONTINUE
     
!     START THE THIRD INTEGRATION LOOP (THRU THE THICKNESS)
     
     DO  izta=1,2
       zta =ptint(izta)
       hzta=zta/2.0
       ibot=(izta-1)*nd2
       
       CALL jacobs (elid,shp,dshp,dgpth,egpdt,epnorm,jacob)
       IF (badjac) GO TO 1600
       
       CALL gmmats (psitrn,3,3,0,jacob,3,3,1,phi)
       
!     CALL Q4BMGS TO GET B MATRIX
!     SET THE ROW FLAG TO 3. IT WILL RETURN THE FIRST 6 ROWS.
       
       rowflg = 3
       CALL q4bmgs (dshp,dgpth,egpdt,epnorm,phi,bmatrx(1))
       DO  ix=1,ndof
         bmatrx(ix+nd2)=xybmat(ibot+ix)
       END DO
       
       IF (.NOT.bendng) GO TO 1410
       DO  ix=1,ndof
         bmatrx(ix+nd5)=xybmat(ibot+ix+ndof)
       END DO
       
!     NOW COMPLETE THE G-MATRIX IF COUPLING EXISTS.
       
       IF (.NOT.mbcoup) GO TO 1410
       DO  ig=1,3
         ig4=(ig+8)*3
         kg=ig+3
         DO  jg=1,3
           jg4=jg+ig4
           jg1=jg4-27
           lg=jg+3
           g(ig,lg)=-gi(jg4)*zta*6.0+gi(jg1)
           g(kg,jg)=-gi(jg4)*zta*6.0+gi(jg1)
         END DO
       END DO
       1410 CONTINUE
       
!-----MULTIPLY DETERMINANT, B-TRANSPOSE, G-MATRIX, & THERMAL
!     STRAIN MATRIX.
       
!                         T
!         P  =  DETERM * B  * G * EPSILON
!          T                             T
       
       IF (compos) GO TO 1430
       DO  i =1,3
         epsubt(i) = detj * talfam(i)
         epsubt(i+3) = - detj * talfab(i) * hzta * thk
       END DO
       GO TO 1450
       
       1430 DO  ir = 1,3
         epsubt(ir  ) = detj*epslnt(ir)
         epsubt(ir+3) =-detj*epslnt(ir+3)*thk*hzta
       END DO
       1450 CONTINUE
       
       CALL gmmats (g,6,6,0,epsubt,6,1,0,gepsbt)
       CALL gmmats (bmatrx,6,ndof,-1,gepsbt,6,1,0,pt)
       
     END DO
   END DO
 END DO
 
!----TRIPLE INTEGRATION LOOP IS NOW FINISHED
 
!****
!     PICK UP THE BASIC TO GLOBAL TRANSFORMATION FOR EACH NODE.
!****
 DO  i=1,36
   trans(i)=0.0
 END DO
 
 DO  i=1,nnode
   ipoint=9*(i-1)+1
   IF (igpdt(1,i) <= 0) GO TO 1510
   
   CALL transs (bgpdt(1,i),tbg)
   GO TO 1530
   
   1510 DO  j=1,9
     tbg(j)=0.0
   END DO
   tbg(1)=1.0
   tbg(5)=1.0
   tbg(9)=1.0
   
   1530 CALL gmmats (teb,3,3,0,tbg,3,3,0,trans(ipoint))
 END DO
 
 
!-----TRANSFORM THE THERMAL LOAD VECTOR INTO THE INDIVIDUAL
!     GLOBAL COORDINATE SYSTEMS OF EACH NODE. NOTE THAT THE
!     TRANSFORMATION MATRICES ARE STORED IN  TRANS = TEG,
!     AND THAT THE 6-DOF LOAD VECTOR FOR EACH NODE USES THE
!     SAME 3X3 TRANSFORMATION MATRIX FOR THE TRANSLATIONAL
!     DOF'S (1-3) AND THE ROTATIONAL DOF'S (4-6).
 
!                         T
!               PT  =  TEG   *  PT
!                 G               E
 
 DO  i =1,nnode
   ipt  = (i-1)*9 + 1
   jpt1 = (i-1)*6 + 1
   jpt2 = jpt1 + 3
   
   CALL gmmats (trans(ipt),3,3,1,pt(jpt1),3,1,0,ptg(jpt1))
   CALL gmmats (trans(ipt),3,3,1,pt(jpt2),3,1,0,ptg(jpt2))
   
 END DO
 
 
!-----WE NOW HAVE THE THERMAL LOAD VECTOR IN GLOBAL COORDINATES,
!     IN PTG. THE NEXT AND LAST STEP IS TO COMBINE IT WITH THE
!     SYSTEM LOAD VECTOR CONTAINED IN Z.
 
 l=0
 DO  i =1,nnode
   k = sil(i) - 1
   DO  j =1,6
     k = k + 1
     l = l + 1
     z(k) = z(k) + ptg(l)
   END DO
 END DO
 GO TO 1600
 
 1580 CALL mesage (30,j,nam)
 nogo = 1
 
 
 1600 CONTINUE
 RETURN
 
 1700 FORMAT ('0*** SYSTEM FATAL ERROR. THE ELEMENT THICKNESS FOR',  &
     ' QUAD4 EID = ',i8,' IS NOT COMPLETELY DEFINED.')
END SUBROUTINE tlqd4s
