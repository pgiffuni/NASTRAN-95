SUBROUTINE tltr3d
     
!     DOUBLE PRECISION ROUTINE TO GENERATE EQUIVALENT THERMAL LOADS FOR
!     THE CTRIA3 ELEMENT.
 
!     WAS NAMED T3THLD (LOADVC,INTZ,Z) IN UAI
 
 
 
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
 
 
 
 LOGICAL :: membrn,bendng,shrflx,mbcoup,norpth,compos,  &
     tempp1,tempp2,sheart,noalfa
 INTEGER :: hunmeg,nest(39),elid,pid,mid(4),sil(3),  &
     igpdt(4,3),necpt(4),iorder(3),comps,flag, indxg2(3,3),sysbuf,nout,nogo
 REAL :: gpth(3),bgpdt(4,3),ecpt(4),tsub0,stemp,z, loadvc(1)
 DOUBLE PRECISION :: pt(6,3),ptg(6,3),pi,twopi,raddeg,degrad,  &
     egpdt(4,3),epnorm(4,3),gpnorm(4,3),dgpth(6),rho,  &
     thetam,cente(3),shpt(3),weight,wtstif,lx,ly,  &
     bmatrx(162),bterms(6),detjac,g(6,6),gi(36),  &
     aic(1),detg2,g2(3,3),egnor(4),mominr,ts,th,  &
     reali,avgthk,tem(9),tbg(9),teb(9),teu(9),tub(9),  &
     tum(9),alpha(6),alfam(3),alfab(3),talfam(3),  &
     talfab(3),tbar,tgrad,tmean,ftherm(6),thrmom(3),  &
     gtemps(3),epslnt(6),epsubt(6),gepsbt(6),  &
     trans(27),offset,tmptrn(36),eps,thetae,edglen(3)
 COMMON /system/  sysbuf,nout,nogo
 COMMON /matin /  matid,inflag,eltemp,dummy,sinmat,cosmat
 COMMON /BLANK /  nrowsp,iparam,comps
 COMMON /condad/  pi,twopi,raddeg,degrad
 COMMON /terms /  membrn,bendng,shrflx,mbcoup,norpth
 COMMON /zzzzzz/  z(1)
 COMMON /sgtmpd/  stemp(7)
 COMMON /trimex/  est(39)
 EQUIVALENCE      (est( 1),nest(1)),(est( 2),sil(1)),  &
     (est( 5),gpth(1)),(est(10),zoff  ), (est(12),elth   ),(est(26),zoff1 ),  &
     (est(39),tempel ),(est(27),bgpdt(1,1),igpdt(1,1))
 EQUIVALENCE      (necpt(1),ecpt(1)),(stemp(7),flag), (z(1),loadvc(1))
 DATA    hunmeg,  eps / 100000000, 1.0D-7 /
 
 
!     INITIALIZE
 
 nnode  = 3
 elid   = nest(1)
 weight = 1.0D0/6.0D0
 sheart =.false.
 noalfa =.false.
 tgrad  = 0.0D0
 eltemp = tempel
 offset = zoff
 IF (zoff == 0.0) offset = zoff1
 
 DO  ll = 1,3
   talfam(ll) = 0.0D0
   talfab(ll) = 0.0D0
   ftherm(ll) = 0.0D0
   ftherm(ll+3) = 0.0D0
 END DO
 
!     TEST FOR COMPOSITE ELEMENT
 
 pid    = nest(11) - hunmeg
 compos = comps == -1 .AND. pid > 0
 
!     CHECK FOR THE TYPE OF TEMPERATURE DATA
!     - TYPE TEMPP1 ALSO INCLUDES TYPE TEMPP3.
!     - IF TEMPPI ARE NOT SUPPLIED, GRID POINT TEMPERATURES ARE PRESENT.
 
 tempp1 = flag == 13
 tempp2 = flag == 2
 
!     SET UP THE ELEMENT FORMULATION
 
 CALL t3setd (ierr,sil,igpdt,elth,gpth,dgpth,egpdt,gpnorm,epnorm,  &
     iorder,teb,tub,cente,avgthk,lx,ly,edglen,elid)
 IF (ierr /= 0) GO TO 520
 CALL gmmatd (teb,3,3,0, tub,3,3,1, teu)
 
!     SET THE NUMBER OF DOF'S
 
 nnod2 = nnode*nnode
 ndof  = nnode*6
 npart = ndof*ndof
 nd2   = ndof*2
 nd6   = ndof*6
 nd7   = ndof*7
 nd8   = ndof*8
 
!     OBTAIN MATERIAL INFORMATION
 
!     PASS THE LOCATION OF THE ELEMENT CENTER FOR MATERIAL
!     TRANSFORMATIONS.
 
 DO  iec = 2,4
   ecpt(iec) = cente(iec-1)
 END DO
 
!     SET MATERIAL FLAGS
!     0.833333333D0 = 5.0D0/6.0D0
 
 IF (nest(13) /=  0) mominr = est(14)
 IF (nest(13) /=  0) ts = est(16)
 IF ( est(16) == .0) ts = 0.83333333D0
 IF (nest(13) == 0 .AND. nest(11) > hunmeg) ts = 0.83333333D0
 
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
 
 CALL shcsgd (*530,nest(9),nest(8),nest(8),nest(21),nest(20),  &
     nest(20),necpt,tub,mcsid,thetam,tum)
 CALL gmmatd (teu,3,3,0, tum,3,3,0, tem)
 
!     CALCULATE THE ANGLE BETWEEN THE MATERIAL AXIS AND THE ELEMENT AXIS
 
 thetae = DATAN2(tem(4),tem(1))
 
!     FETCH MATERIAL PROPERTIES
 
 CALL shgmgd (*540,elid,tem,mid,ts,noalfa,gi,rho,gsube,tsub0, egnor,alpha)
 
 DO  ial = 1,3
   alfam(ial) = alpha(ial  )
   alfab(ial) = alpha(ial+3)
 END DO
 
!     TURN OFF THE COUPLING FLAG WHEN MID4 IS PRESENT WITH ALL
!     CALCULATED ZERO TERMS.
 
 IF (.NOT.mbcoup) GO TO 50
 DO  i = 28,36
   IF (DABS(gi(i)) > eps) GO TO 50
 END DO
 mbcoup = .false.
 
!     OBTAIN TEMPERATURE INFORMATION
 
!     IF TEMPP1 DATA, GET AVERAGE TEMP AND THERMAL GRADIENT.
 
 50 IF (.NOT.tempp1) GO TO 60
 tmean = stemp(1)
 tgrad = stemp(2)
 GO TO 90
 
!     IF TEMPP2 DATA, GET THERMAL MOMENTS.
 
 60 IF (.NOT.tempp2) GO TO 70
 tmean = stemp(1)
 
 thrmom(1) = stemp(2)
 thrmom(2) = stemp(3)
 thrmom(3) = stemp(4)
 
 ftherm(4) = thrmom(1)
 ftherm(5) = thrmom(2)
 ftherm(6) = thrmom(3)
 GO TO 90
 
!     TEMPPI TEMPERATURE DATA IS NOT AVAILABLE, THEREFORE SORT THE GRID
!     POINT TEMPERATURES (IN STEMP(1-7)).
 
 70 DO  i = 1,nnode
   ipnt = iorder(i)
   gtemps(i) = stemp(ipnt)
 END DO
 tmean = (gtemps(1)+gtemps(2)+gtemps(3))/3.0D0
 90 tbar = tmean - tsub0
 
!     CALCULATE THERMAL STRAINS FOR COMPOSITE ELEMENTS
 
 IF (.NOT.compos) GO TO 100
 CALL shctsd (ierr,elid,pid,mid,avgthk,tmean,tgrad,thetae,ftherm, epslnt,z,z)
 IF (ierr /= 0) GO TO 500
 
!     INITIALIZE FOR THE MAIN INTEGRATION LOOP
 
 100 DO  i = 1,6
   epsubt(i) = 0.0D0
   DO  j = 1,nnode
     pt (i,j)  = 0.0D0
     ptg(i,j)  = 0.0D0
   END DO
 END DO
 
!     MAIN INTEGRATION LOOP
 
 DO  ipt = 1,nnode
   CALL t3bmgd (ierr,sheart,ipt,iorder,egpdt,dgpth,aic,th,detjac,  &
       shpt,bterms,bmatrx)
   IF (ierr /= 0) GO TO 520
   
   wtstif = detjac*weight
   reali  = mominr*th*th*th/12.0D0
   
!     FILL IN THE 6X6 [G]
   
   DO  ig = 1,6
     DO  jg = 1,6
       g(ig,jg) = 0.0D0
     END DO
   END DO
   
   IF (.NOT.membrn) GO TO 220
   DO  ig = 1,3
     ig1 = (ig-1)*3
     DO  jg = 1,3
       g(ig,jg) = gi(ig1+jg)*th
     END DO
   END DO
   
   220 IF (.NOT.bendng) GO TO 250
   DO  ig = 4,6
     ig2 = (ig-2)*3
     DO  jg = 4,6
       g(ig,jg) = gi(ig2+jg)*reali
     END DO
   END DO
   
   IF (.NOT.mbcoup) GO TO 250
   DO  ig = 1,3
     ig4 = (ig+8)*3
     DO  jg = 1,3
       g(ig,jg+3) = gi(ig4+jg)*th*th
       g(ig+3,jg) = g(ig,jg+3)
     END DO
   END DO
   
!     PREPARE THERMAL STRAINS FOR COMPOSITE ELEMENTS
   
   250 IF (.NOT.compos) GO TO 270
   DO  ir = 1,6
     epsubt(ir) = wtstif*epslnt(ir)
   END DO
   GO TO 370
   
!     CALCULATE THERMAL STRAINS FOR NON-COMPOSITE ELEMENTS
   
   270 IF (.NOT.membrn) GO TO 290
   DO  i = 1,3
     talfam(i) = tbar*alfam(i)
   END DO
   
   290 IF (.NOT.bendng) GO TO 350
   IF (.NOT.tempp1) GO TO 310
   DO  i = 1,3
     talfab(i) = -tgrad*alfab(i)
   END DO
   GO TO 350
   
   310 IF (.NOT.tempp2) GO TO 330
   DO  ig = 1,3
     DO  jg = 1,3
       g2(ig,jg) = g(ig+3,jg+3)
     END DO
   END DO
   
   CALL inverd (3,g2,3,gdum,0,detg2,isngg2,indxg2)
   CALL gmmatd (g2,3,3,0, thrmom,3,1,0, talfab)
   GO TO 350
   
   330 DO  i = 1,3
     talfab(i) = 0.0D0
   END DO
   
   350 DO  i = 1,3
     epsubt(i  ) = wtstif*talfam(i)
     epsubt(i+3) = wtstif*talfab(i)
   END DO
   
!                                T
!     [P]  = [P]  + WTSTIF*[B] [G][EPS]
!        T      T                      T
   
   370 CALL gmmatd (g,6,6,0, epsubt,6,1,0, gepsbt)
   CALL gmmatd (bmatrx,6,ndof,-1, gepsbt,6,1,0, pt)
   
 END DO
 
!     END OF MAIN INTEGRATION LOOP
 
!     PICK UP THE ELEMENT TO GLOBAL TRANSFORMATION FOR EACH NODE.
 
 DO  i = 1,nnode
   ipoint = 9*(i-1)+1
   CALL transd (bgpdt(1,i),tbg)
   CALL gmmatd (teb,3,3,0, tbg,3,3,0, trans(ipoint))
 END DO
 
!     TRANSFORM THE THERMAL LOAD VECTOR INTO THE INDIVIDUAL GLOBAL
!     COORDINATE SYSTEMS OF EACH NODE.
 
!                 T
!     [PT] = [TEG] [PT]
!         G            E
 
 DO  i = 1,nnode
   CALL tldrd  (offset,i,trans,tmptrn)
   CALL gmmatd (tmptrn,6,6,1, pt(1,i),6,1,0, ptg(1,i))
 END DO
 
!     ADD THE THERMAL LOAD VECTOR TO THE GLOBAL LOAD VECTOR WHICH
!     RESIDES IN [LOADVC].
 
 DO  i = 1,nnode
   k = sil(i) - 1
   DO  j = 1,6
     loadvc(k+j) = loadvc(k+j) + SNGL(ptg(j,i))
   END DO
 END DO
 GO TO 600
 
!     FATAL ERRORS
 
 500 WRITE  (nout,510)
 510 FORMAT ('0*** SYSTEM FATAL ERROR.  APPROPRIATE COMPOSITE DATA ',  &
     'NOT FOUND IN MODULE SSG1.')
 GO TO 560
 520 j = 224
 GO TO 550
 530 j = 225
 nest(2) = mcsid
 GO TO 550
 540 j = 226
 nest(2) = mid(3)
 550 CALL mesage (30,j,nest(1))
 560 nogo = 1
 
 600 RETURN
END SUBROUTINE tltr3d
