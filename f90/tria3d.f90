SUBROUTINE tria3d
     
!     DOUBLE PRECISION ROUTINE TO FORM STIFFNESS, MASS, AND DAMPING
!     MATRICES FOR THE CTRIA3 ELEMENT
 
!                 EST  LISTING
 
!        WORD     TYP       DESCRIPTION
!     ----------------------------------------------------------------
!     ECT:
!         1        I   ELEMENT ID, EID
!         2-4      I   SIL LIST, GRIDS 1,2,3
!         5-7      R   MEMBRANE THICKNESSES T, AT GRIDS 1,2,3
!         8        R   MATERIAL PROPERTY ORIENTAION ANGLE, THETA
!               OR I   COORD. SYSTEM ID (SEE TM ON CTRIA3 CARD)
!         9        I   TYPE FLAG FOR WORD 8
!        10        R   GRID OFFSET, ZOFF
!    EPT:
!        11        I   MATERIAL ID FOR MEMBRANE, MID1
!        12        R   ELEMENT THICKNESS,T (MEMBRANE, UNIFORMED)
!        13        I   MATERIAL ID FOR BENDING, MID2
!        14        R   MOMENT OF INERTIA FACTOR, I (BENDING)
!        15        I   MATERIAL ID FOR TRANSVERSE SHEAR, MID3
!        16        R   TRANSV. SHEAR CORRECTION FACTOR, TS/T
!        17        R   NON-STRUCTURAL MASS, NSM
!        18-19     R   STRESS FIBER DISTANCES, Z1,Z2
!        20        I   MATERIAL ID FOR MEMBRANE-BENDING COUPLING, MID4
!        21        R   MATERIAL ANGLE OF ROTATION, THETA
!               OR I   COORD. SYSTEM ID (SEE MCSID ON PSHELL CARD)
!                      (DEFAULT FOR WORD 8)
!        22        I   TYPE FLAG FOR WORD 21 (DEFAULT FOR WORD 9)
!        23        I   INTEGRATION ORDER FLAG
!        24        R   STRESS ANGLE OF RATATION, THETA
!               OR I   COORD. SYSTEM ID (SEE SCSID ON PSHELL CARD)
!        25        I   TYPE FLAG FOR WORD 24
!        26        R   OFFSET, ZOFF1 (DEFAULT FOR WORD 10)
!    BGPDT:
!        27-38   I/R   CID,X,Y,Z  FOR GRIDS 1,2,3
!    ETT:
!        39        I   ELEMENT TEMPERATURE
 
 
 LOGICAL :: heat,noalfa,needk,needm,sheart, membrn,bendng,shrflx,mbcoup,norpth
 INTEGER :: sysbuf,nout,nogo,prec,hunmeg,nest(39),NAME(2),  &
     necpt(4),dict(11),igpdt(4,3),elid,estid,dmat,  &
     sil(3),iorder(3),cpmass,mid(4),TYPE,INDEX(3,3)
 REAL :: bgpdt(4,3),gpth(3),nsm,ecpt(4),kheat,htcp
 DOUBLE PRECISION :: amgg(1),akgg(1),alpha(1),thetam,cente(3),  &
     dgpth(3),egpdt(4,3),epnorm(4,3),gpnorm(4,3),  &
     area,wtstif,wtmass,rho,xmass(9),xmasso,lx,ly,  &
     eps,offset,shpt(3),weight,g(9,9),gi(36),k11,k22,  &
     jok,jog,zz(9),aic(18),egnor(4),edglen(3),  &
     bmtrx(54),bmatrx(162),bterms(6),bmat1(486),  &
     avgthk,mominr,ts,th,reali,tsi,tsm,bdum(3),  &
     determ,detjac,tbg(9),teb(9),tem(9),teu(9),  &
     tub(9),tum(9),tottrn(324),transk(324),trans(27),  &
     tmptrn(36),htflx(18),htcap(36),htcon(36), dheat,weitc,dvol
 COMMON /system/  sysbuf,nout,nogo,idum(51),prec
 COMMON /matin /  matid,inflag,eltemp,dummy,sinmat,cosmat
 COMMON /hmtout/  kheat(7),TYPE
 COMMON /terms /  membrn,bendng,shrflx,mbcoup,norpth
 COMMON /emgprm/  icore,jcore,ncore,icstm,ncstm,imat,nmat,ihmat,  &
     nhmat,idit,ndit,icong,ncong,lcong,anycon,  &
     kgg1,mgg1,ibgg1,precis,error,heat,cpmass, dumm6(6),l38
 COMMON /emgest/  est(39)
 COMMON /emgdic/  eltype,ldict,nlocs,elid,estid
 COMMON /zzzzzz/  z(1)
 EQUIVALENCE      (est( 1),nest(1)), (est( 2),sil(1)),  &
     (est( 5),gpth(1)), (est(10),zoff), (est(12),elth)   , (est(17),nsm),  &
     (est(23),INT)    , (est(26),zoff1), (est(27),bgpdt(1,1),igpdt(1,1)),  &
     (est(39),tempel) , (dict(5),adamp),  &
     (necpt(1),ecpt(1)),(z(1),amgg(1),akgg(1)),  &
     (kheat(4),htcp)  , (htcap(1),xmass(1))
 DATA     hunmeg, eps / 100000000, 1.0D-7 /
 DATA     NAME  , kmat, mmat, dmat / 4HCTRI,4HA3  , 1, 2, 3 /
 
!     INITIALIZE
 
 elid   = nest(1)
 nnode  = 3
 mominr = 0.0D0
 ts     = 0.0D0
 weight = 1.0D0/6.0D0
 eltemp = tempel
 needk  = kgg1 /= 0 .OR. ibgg1 /= 0
 noalfa = .true.
 sheart = .true.
 ieoe   = 1
 offset = zoff
 IF (zoff == 0.0) offset = zoff1
 
!     CHECK FOR SUFFICIENT OPEN CORE FOR ELEMENT STIFFNESS
 
!     OPEN CORE BEGINS AT JCORE
!     OPEN CORE ENDS   AT NCORE
!     LENGTH OF AVAILABLE WORDS = (NCORE-JCORE-1)/PREC
 
 jcored = jcore/prec + 1
 length = (ncore-jcore-1)/prec
 IF (length < 324 .AND. (.NOT.heat .AND. needk)) GO TO 1100
 
!     SET UP THE ELEMENT FORMULATION
 
 CALL t3setd (ierr,sil,igpdt,elth,gpth,dgpth,egpdt,gpnorm,epnorm,  &
     iorder,teb,tub,cente,avgthk,lx,ly,edglen,elid)
 IF (ierr /= 0) GO TO 1110
 CALL gmmatd (teb,3,3,0, tub,3,3,1, teu)
 area = lx*ly/2.0D0
 
!     SET THE NUMBER OF DOF'S
 
 nnod2 = nnode*nnode
 ndof  = nnode*6
 npart = ndof*ndof
 nd2   = ndof*2
 nd6   = ndof*6
 nd7   = ndof*7
 nd8   = ndof*8
 nd9   = ndof*9
 jend  = jcored + npart - 1
 
!     OBTAIN MATERIAL INFORMATION
 
!     PASS THE LOCATION OF THE ELEMENT CENTER FOR MATERIAL
!     TRANSFORMATIONS.
 
 DO  iec = 2,4
   ecpt(iec) = cente(iec-1)
 END DO
 
!     SET MATERIAL FLAGS
!     5.0D0/6.0D0 = 0.833333333D0
 
 IF (nest(13) /=   0) mominr = est(14)
 IF (nest(13) /=   0) ts = est(16)
 IF ( est(16) == 0.0) ts = 0.833333333D0
 IF (nest(13) == 0 .AND. nest(11) > hunmeg) ts = 0.833333333D0
 
 mid(1) = nest(11)
 mid(2) = nest(13)
 mid(3) = nest(15)
 mid(4) = nest(20)
 
 membrn = mid(1) > 0
 bendng = mid(2) > 0 .AND. mominr > 0.0D0
 shrflx = mid(3) > 0
 mbcoup = mid(4) > 0
 norpth = mid(1) == mid(2) .AND. mid(1) == mid(3) .AND. mid(4) == 0  &
     .AND. DABS(mominr-1.0D0) <= eps
 
!     SET UP TRANSFORMATION MATRIX FROM MATERIAL TO ELEMENT COORD.SYSTEM
 
 CALL shcsgd (*1120,nest(9),nest(8),nest(8),nest(21),nest(20),  &
     nest(20),necpt,tub,mcsid,thetam,tum)
 
!     BRANCH ON FORMULATION TYPE.
 
 IF (heat) GO TO 800
 
!     FETCH MATERIAL PROPERTIES
 
 CALL gmmatd (teu,3,3,0,tum,3,3,0,tem)
 CALL shgmgd (*1130,elid,tem,mid,ts,noalfa,gi,rho,gsube,tsub0, egnor,alpha)
 
!     TURN OFF THE COUPLING FLAG WHEN MID4 IS PRESENT WITH ALL
!     CALCULATED ZERO TERMS.
 
 IF (.NOT.mbcoup) GO TO 120
 DO  i = 28,36
   IF (DABS(gi(i)) > eps) GO TO 120
 END DO
 mbcoup = .false.
 
!     GET THE GEOMETRY CORRECTION TERMS
 
 120 IF (.NOT.bendng) GO TO 130
 CALL t3gemd (ierr,egpdt,iorder,gi(10),gi(19),lx,ly,edglen,shrflx,  &
     aic,jog,jok,k11,k22)
 IF (ierr /= 0) GO TO 1110
 
!     REDUCED INTEGRATION LOOP FOR STIFFNESS
 
 130 IF (.NOT.needk .OR. INT /= 0) GO TO 160
 
!     DETERMINE THE AVERAGE B-MATRIX FOR OUT-OF-PLANE SHEAR
 
 DO  ipt = 1,3
   kpt = (ipt-1)*nd9 + 1
   CALL t3bmgd (ierr,sheart,ipt,iorder,egpdt,dgpth,aic,th,detjac,  &
       shpt,bterms,bmat1(kpt))
   IF (ierr /= 0) GO TO 1110
 END DO
 
 DO  i = 1,ndof
   bmtrx(i     ) = bmat1(i+nd6) +bmat1(i+nd6+nd9) +bmat1(i+nd6+2*nd9)
   bmtrx(i+ndof) = bmat1(i+nd7) +bmat1(i+nd7+nd9) +bmat1(i+nd7+2*nd9)
   bmtrx(i+nd2 ) = bmat1(i+nd8) +bmat1(i+nd8+nd9) +bmat1(i+nd8+2*nd9)
 END DO
 
!     INITIALIZE FOR THE MAIN INTEGRATION LOOP
 
 160 needm = mgg1 /= 0 .AND. (nsm > 0.0 .OR. rho > 0.0D0)
 IF (.NOT.needk .AND. .NOT.needm) GO TO 200
 DO  i = jcored,jend
   akgg(i) = 0.0D0
 END DO
 
 DO  i = 1,9
   xmass(i) = 0.0D0
 END DO
 
!     MAIN INTEGRATION LOOP
 
 200 DO  ipt = 1,3
   
   CALL t3bmgd (ierr,sheart,ipt,iorder,egpdt,dgpth,aic,th,detjac,  &
       shpt,bterms,bmatrx)
   IF (ierr /= 0) GO TO 1110
   
!     PERFORM STIFFNESS CALCULATIONS IF REQUIRED
   
   IF (.NOT.needk) GO TO 400
   wtstif = detjac*weight
   reali  = mominr*th*th*th/12.0D0
   tsi = ts*th
   
   IF (INT /= 0) GO TO 220
   DO  ix = 1,ndof
     bmatrx(ix+nd6) = bmtrx(ix     )
     bmatrx(ix+nd7) = bmtrx(ix+ndof)
     bmatrx(ix+nd8) = bmtrx(ix+nd2 )
   END DO
   
!     FILL IN THE 9X9 G-MATRIX
   
   220 DO  ig = 1,81
     g(ig,1) = 0.0D0
   END DO
   
   IF (.NOT.membrn) GO TO 270
   DO  ig = 1,3
     ig1 = (ig-1)*3
     DO  jg = 1,3
       g(ig,jg) = gi(ig1+jg)*th*wtstif
     END DO
   END DO
   
   270 IF (.NOT.bendng) GO TO 340
   DO  ig = 4,6
     ig2 = (ig-2)*3
     DO  jg = 4,6
       g(ig,jg) = gi(ig2+jg)*reali*wtstif
     END DO
   END DO
   
   tsm   = 1.0D0/(2.0D0*12.0D0*reali)
   zz(1) = (jog/tsi)* gi(22) + tsm*jok*k22
   zz(2) =-(jog/tsi)*(gi(20) + gi(21))/2.0D0
   zz(3) = 0.0D0
   zz(4) = zz(2)
   zz(5) = (jog/tsi)* gi(19) + tsm*jok*k11
   zz(6) = 0.0D0
   zz(7) = 0.0D0
   zz(8) = 0.0D0
   zz(9) = (jog/tsi)*(gi(22) + gi(19))/2.0D0  &
       + tsm*12.0D0*area/DSQRT(gi(10)*gi(14))
   CALL inverd (3,zz,3,bdum,0,determ,ising,INDEX)
   IF (ising /= 1) GO TO 1110
   
   DO  ig = 7,9
     ig3 = (ig-7)*3
     DO  jg = 7,9
       g(ig,jg) = zz(ig3+jg-6)*wtstif
     END DO
   END DO
   
   IF (.NOT.mbcoup) GO TO 340
   DO  ig = 1,3
     ig4 = (ig+8)*3
     DO  jg = 1,3
       g(ig,jg+3) = gi(ig4+jg)*th*th*wtstif
       g(ig+3,jg) = g(ig,jg+3)
     END DO
   END DO
   
!     COMPUTE THE CONTRIBUTION TO THE STIFFNESS MATRIX FROM THIS
!     INTEGRATION POINT.
   
   340 CALL t3bgbd (9,ndof,g,bmatrx,akgg(jcored))
   
   
!     END OF STIFFNESS CALCULATIONS.
!     SKIP MASS CALCULATIONS IF NOT REQUIRED
   
   
   400 IF (.NOT.needm) CYCLE
   wtmass = (rho*th+nsm)*detjac*weight
   IF (cpmass <= 0) GO TO 430
   
!     CONSISTENT MASS FORMULATION (OPTION)
   
   DO  i = 1,nnode
     ii = (i-1)*nnode
     DO  j = 1,nnode
       xmass(ii+j) = xmass(ii+j) + shpt(i)*shpt(j)*wtmass
     END DO
   END DO
   CYCLE
   
!     LUMPED MASS FORMULATION (DEFAULT)
   
   430 i3 = 1
   DO  i = 1,nnode
     xmass(i3) = xmass(i3) + shpt(i)*wtmass
     i3 = i3 + 1 + nnode
   END DO
   
!     END OF MAIN INTEGRATION LOOP
   
 END DO
 
!     PICK UP THE ELEMENT TO GLOBAL TRANSFORMATION FOR EACH NODE.
 
 DO  i = 1,nnode
   ipoint = 9*(i-1) + 1
   CALL transd (igpdt(1,i),tbg)
   CALL gmmatd (teb,3,3,0, tbg,3,3,0, trans(ipoint))
 END DO
 
!     SHIP OUT THE STIFFNESS AND DAMPING MATRICES
 
 IF (.NOT.needk) GO TO 600
 
 dict(1) = estid
 dict(2) = 1
 dict(3) = ndof
 dict(4) = 63
 adamp   = gsube
 
!     BUILD THE 18X18 TRANSFORMATION MATRIX FOR ONE-SHOT MULTIPLY
 
 DO  i = 1,npart
   transk(i) = 0.0D0
   tottrn(i) = 0.0D0
 END DO
 
 ndof66 = 6*ndof + 6
 ii = 1
 DO  i = 1,npart,ndof66
   CALL tldrd (offset,ii,trans,tmptrn)
   DO  jj = 1,36,6
     j  = jj - 1
     kk = i - 1 + j*nnode
     DO  k = 1,6
       tottrn(kk+k) = tmptrn(j+k)
     END DO
   END DO
   ii = ii + 1
 END DO
 
!     PERFORM THE TRIPLE MULTIPLY.
 
 CALL mpya3d (tottrn,akgg(jcored),ndof,6,transk)
 
 CALL emgout (transk,transk,npart,ieoe,dict,kmat,prec)
 
!     SHIP OUT THE MASS MATRIX
 
 600 IF (.NOT.needm) GO TO 730
 ndof    = nnode*3
 npart   = ndof*ndof
 dict(2) = 1
 dict(3) = ndof
 dict(4) = 7
 adamp   = 0.0
 jend    = jcored + npart - 1
 
!     ZERO OUT THE POSITIONS, THEN LOOP ON I AND J TO LOAD THE MASS
!     MATRIX.
 
 DO  ijk = jcored,jend
   amgg(ijk) = 0.0D0
 END DO
 
 ndofp1 = ndof + 1
 DO  ii = 1,nnod2  ,nnode
   i = ii - 1
   DO  j = 1,nnode
     xmasso = xmass(i+j)
     ipoint = (j-1)*3 + i*9 + jcored
     jpoint = ipoint + 3*ndof
     DO  k = ipoint,jpoint,ndofp1
       amgg(k) = xmasso
     END DO
   END DO
 END DO
 
!     BYPASS TRANSFORMATIONS IF LUMPED MASS.
 
 IF (cpmass <= 0) GO TO 700
 
!     BUILD THE 9X9 TRANSFORMATION MATRIX FOR ONE-SHOT MULTIPLY
 
 DO  i = 1,npart
   transk(i) = 0.0D0
   tottrn(i) = 0.0D0
 END DO
 
 ndof33 = 3*ndof + 3
 DO  i = 1,npart  ,ndof33
   ii = ((i-1)/(3*ndof))*9
   DO  jj = 1,9,3
     j  = jj - 1
     kk = i - 1 + j*nnode
     DO  k = 1,3
       tottrn(kk+k) = trans(ii+j+k)
     END DO
   END DO
 END DO
 
!     PERFORM THE TRIPLE MULTIPLY.
 
 CALL mpya3d (tottrn,amgg(jcored),ndof,3,transk)
 GO TO 720
 
!     JUST COPY THE LUMPED MASS MATRIX OUT
 
 700 ii = jcored
 DO  i = 1,npart
   transk(i) = amgg(ii)
   ii = ii + 1
 END DO
 
 720 CALL emgout (transk,transk,npart,ieoe,dict,mmat,prec)
 
 730 CONTINUE
 GO TO 1200
 
!     HEAT CALCULATIONS
 
 800 CONTINUE
 inflag = 2
 sinmat = DSIN(thetam)
 cosmat = DCOS(thetam)
 matid  = nest(11)
 
 CALL hmat (elid)
 
 gi(1) = kheat(1)
 gi(2) = kheat(2)
 gi(3) = gi(2)
 gi(4) = kheat(3)
 
 DO  i = 1,18
   htcon(i) = 0.0D0
   htcap(i) = 0.0D0
 END DO
 
!     BEGIN LOOP ON INTEGRATION POINTS
 
 DO  ipt = 1,3
   CALL t3bmgd (ierr,sheart,ipt,iorder,egpdt,dgpth,aic,th,detjac,  &
       shpt,bterms,bmatrx)
   IF (ierr /= 0) GO TO 1110
   
   dvol = weight*detjac*th
   DO  i = 1,4
     g(i,1) = gi(i)*dvol
   END DO
   weitc = dvol*htcp
   
   ip = 1
   DO  i = 1,nnode
     htflx(ip  ) = g(1,1)*bterms(i) + g(2,1)*bterms(i+nnode)
     htflx(ip+1) = g(3,1)*bterms(i) + g(4,1)*bterms(i+nnode)
     ip = ip + 2
   END DO
   CALL gmmatd (bterms,2,nnode,-1, htflx,nnode,2,1, htcon)
   
!     FINISHED WITH HEAT CONDUCTIVITY MATRIX, DO HEAT CAPACITY IF
!     REQUIRED.
   
   IF (htcp == 0.0) CYCLE
   ip = 1
   DO  i = 1,nnode
     dheat = weitc*shpt(i)
     DO  j = 1,nnode
       htcap(ip) = htcap(ip) + dheat*shpt(j)
       ip = ip + 1
     END DO
   END DO
   
 END DO
 
!     END OF INTEGRATION LOOP, SHIP OUT THE RESULTS.
 
 dict(1) = estid
 dict(2) = 1
 dict(3) = nnode
 dict(4) = 1
 IF (weitc == 0.0D0) GO TO 1000
 adamp = 1.0
 CALL emgout (htcap,htcap,nnod2,ieoe,dict,dmat,prec)
 1000 adamp = 0.0
 CALL emgout (htcon,htcon,nnod2,ieoe,dict,kmat,prec)
 
 GO TO 1200
 
 
!     FATAL ERRORS
 
!     INSUFFICIENT MEMORY IS AVAILABLE
 
 1100 CALL mesage (-30,228,NAME)
 GO TO 1140
 
!     CTRIA3 ELEMENT HAS ILLEGAL GEOMETRY OR CONNECTIONS
 
 1110 j = 224
 GO TO 1140
 
!     THE X-AXIS OF THE MATERIAL COORDINATE SYSTEM HAS NO PROJECTION
!     ON TO THE PLANE OF CTRIA3 ELEMENT
 
 1120 j = 225
 nest(2) = mcsid
 GO TO 1140
 
!     ILLEGAL DATA DETECTED ON MATERIAL ID REFERENCED BY CTRIA3 ELEMENT
!     FOR MID3 APPLICATION
 
 1130 j = 226
 nest(2) = mid(3)
 
 1140 CALL mesage (30,j,nest(1))
 IF (l38 == 1) CALL mesage (-61,0,0)
 nogo = 1
 
 1200 CONTINUE
 RETURN
END SUBROUTINE tria3d
