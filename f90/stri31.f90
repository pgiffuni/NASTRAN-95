SUBROUTINE stri31
     
!     ROUTINE TO RECOVER CTRIA3 ELEMENT FORCES, STRAINS, AND STRESSES.
!     PHASE 1.
 
!     WAS NAMED T3ST1D/S IN UAI CODE
 
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
 
 
!     ****************        RESIDES IN COMMON BLOCK SDR2X5 (AFTER EST)
!     PH1OUT DATA BLOCK       TOTAL NO. OF WORDS = 713
!     ****************
 
!     PH1OUT( 1)    = ELID, ELEMENT ID
!     PH1OUT( 2- 4) = SIL NUMBERS
!     PH1OUT( 5- 7) = ARRAY IORDER
!     PH1OUT( 8)    = TSUB0, REFERENCE TEMP.
!     PH1OUT( 9-10) = Z1 & Z2, FIBER DISTANCES
!     PH1OUT(11)    = ID OF THE ORIGINAL PCOMPI PROPERTY ENTRY
!     PH1OUT(12)    = DUMMY WORD (FOR ALLIGNMENT)
 
!     PH1RST( 1)    = AVGTHK, AVERAGE THICKNESS
!     PH1RST( 2)    = MOMINR, MOMENT OF INER. FACTOR
!     PH1RST( 3-38) = 6X6 MATERIAL PROPERTY MATRIX (NO SHEAR)
!     PH1RST(39-41) = THERMAL EXPANSION COEFFICIENTS FOR MEMBRANE
!     PH1RST(42-44) = THERMAL EXPANSION COEFFICIENTS FOR BENDING
!     PH1RST(45-47) = CORNER NODE THICKNESSES
!     PH1RST(48)    = OFFSET OF ELEMENT FROM GP PLANE
!     PH1RST(49-57) = 3X3 USER-TO-MATERIAL COORD. TRNASF. MATRIX UEM
!     PH1RST(58-66) = 3X3 ELEM-TO-STRESS/STRAIN TRANSF. TENSOR TES
!     PH1RST(67-93) = THREE 3X3 GLOBAL-TO-ELEM COORD. TRANSFORMATION
!                     NODAL MATRICES TEG, ONE FOR EACH NODE
 
!     THE FOLLOWING IS REPEATED FOR EACH EVALUATION POINT (4 TIMES, AT
!     THE CENTER OF THE ELEMENT AND AT 3 STANDARD TRIANGULAR POINTS).
!     THE CHOICE OF THE FINAL STRESS/FORCE OUTPUT POINTS IS MADE AT THE
!     SUBCASE LEVEL (PHASE 2).
 
!              1             ELEMENT THICKNESS AT THIS POINT
!            2 - 5           OUT-OF-PLANE-SHEAR-FORCE/STRAIN MATRIX
!            6 - 8           ELEMENT SHAPE FUNCTION VALUES
!          8+1 - 8+8*NDOF    STRAIN RECOVERY MATRIX
 
 
!     *****************      RESIDES IN COMMON BLOCK SDR2X6
!     IELOUT DATA BLOCK      CONTAINS DATA FOR GPSRN
!     *****************      (TOTAL NO OF WORDS =  77)
 
!              1             ELEMENT ID
!              2             AVERAGE THICKNESS
 
!     THE FOLLOWING IS REPEATED FOR EACH NODE.
 
!         WORD  1            SIL NUMBER
!         WORD  2-10         [TSB] FOR Z1
!         WORD 11-19         [TSB] FOR Z2
!         WORD 20-22         NORMAL VECTOR IN BASIC COORD. SYSTEM
!         WORD 23-25         GRID COORDS   IN BASIC COORD. SYSTEM
 
 
 LOGICAL :: membrn,bendng,shrflx,mbcoup,norpth, sheart,noalfa,userst
 INTEGER :: nphi(100),nest(39),sil(3),silo,igpdt(4,3),elid,  &
     necpt(4),mid(4),pido,scsid,flags,hunmeg,itherm,  &
     NAME(2),INDEX(2,3),sysbuf,nout,nogo,prec
 REAL :: relout(300),gpth(3),bgpdt(4,3),ecpt(4),  &
     kheat,htcp,tsub0,gsube,eltemp,z1o,z2o
 REAL :: ph1rst,dgpth(3),th,avgthk,gpnorm(4,3),  &
     epnorm(4,3),egpdt(4,3),cente(3),offset,alpha(6),  &
     gi(36),rho,jog,jok,k11,k22,zz(4),aic(18),  &
     egnor(4),tss,lx,ly,edglen(3),mominr,ts,reali,tsi,  &
     teu(9),tes(9),teb(9),tbg(9),tus(9),tem(9),tsb(9),  &
     tsm(9),tub(9),tum(9),thetam,thetas,uem(9),vem(4),  &
     determ,bdum(3),shpt(3),bterms(6),bmat1(486), bmatrx(162),bmtrx(36),drkce(33)
 COMMON /system/ sysbuf,nout,nogo,idum(51),prec,itherm
 COMMON /sdr2x5/ est(100),  &
     elid,silo(3),iorder(3),tsub0,z1o,z2o,pido,idumal, ph1rst(701)
 COMMON /sdr2x6/ ielout(300)
 COMMON /matin / matid,inflag,eltemp,dummy,sinmat,cosmat
 COMMON /hmtout/ kheat(7),TYPE
 COMMON /terms / membrn,bendng,shrflx,mbcoup,norpth
 EQUIVALENCE     (est( 1),nest(1) ), (est( 2),sil(1)),  &
     (est( 5),gpth(1) ), (est(10),zoff  ),  &
     (est(12),elth    ), (est(23),INT   ),  &
     (est(26),zoff1   ), (est(39),tempel), (est(27),bgpdt(1,1), igpdt(1,1))
 EQUIVALENCE     (nphi( 1),elid   ), (nphi (27),drkce(1) ),  &
     (necpt(1),ecpt(1)), (ielout(1),relout(1)), (htcp,kheat(4))
 DATA    hunmeg, istart / 100000000, 93 /,  eps / 1.0E-17 /
 DATA    NAME  / 4HTRIA , 4H3           /
 
!     INITIALIZE
 
 nnode  = 3
 mominr = 0.0
 ts     = 0.0
 eltemp = tempel
 elid   = nest(1)
 z1o    = est(18)
 z2o    = est(19)
 pido   = nest(11) - hunmeg
 mcsid  = nest(21)
 scsid  = nest(24)
 flags  = nest(25)
 userst = scsid < 0 .AND. flags == 1
 noalfa = .false.
 sheart = .true.
 offset = zoff
 IF (zoff == 0.0) offset = zoff1
 
!     START FILLING IN THE DATA BLOCKS
 
 ielout(1) = elid
 DO  i = 1,3
   ielout(3+(i-1)*25) = sil(i)
   DO  j = 1,3
     relout(25*i+j-1) = bgpdt(j+1,i)
   END DO
 END DO
 
!     SET UP THE ELEMENT FORMULATION
 
 CALL t3sets (ierr,sil,igpdt,elth,gpth,dgpth,egpdt,gpnorm,epnorm,  &
     iorder,teb,tub,cente,avgthk,lx,ly,edglen,elid)
 IF (ierr /= 0) GO TO 600
 CALL gmmats (teb,3,3,0, tub,3,3,1, teu)
 DO  i = 1,3
   silo(i) = sil(i)
 END DO
 
!     SET THE NUMBER OF DOF'S
 
 nnod2 = nnode*nnode
 ndof  = nnode*6
 npart = ndof*ndof
 nd2   = ndof*2
 nd6   = ndof*6
 nd7   = ndof*7
 nd8   = ndof*8
 nd9   = ndof*9
 
!     PASS THE LOCATION OF THE ELEMENT CENTER FOR TRANSFORMATIONS.
 
 DO  iec = 2,4
   ecpt(iec) = cente(iec-1)
 END DO
 
!     STRESS TRANSFORMATIONS
 
 IF (.NOT.userst) GO TO 50
 est (24) = 0.0
 nest(25) = 0
 50 CALL shcsgs (*620,nest(25),nest(24),est(24),nest(25),nest(24),  &
     est(24),necpt,tub,scsid,thetas,tus)
 CALL gmmats (teu,3,3,0, tus,3,3,0, tes)
 
!     OBTAIN MATERIAL INFORMATION
 
!     SET MATERIAL FLAGS
!     0.83333333 = 5.0/6.0
 
 IF (nest(13) /=   0) mominr = est(14)
 IF (nest(13) /=   0) ts = est(16)
 IF ( est(16) == 0.0) ts = 0.83333333
 IF (nest(13) == 0 .AND. nest(11) > hunmeg) ts = 0.83333333
 
 mid(1) = nest(11)
 mid(2) = nest(13)
 mid(3) = nest(15)
 mid(4) = nest(20)
 
 membrn = mid(1) > 0
 bendng = mid(2) > 0 .AND. mominr > 0.0
 shrflx = mid(3) > 0
 mbcoup = mid(4) > 0
 norpth = .false.
 
!     SET UP TRANSFORMATION MATRIX FROM MATERIAL TO ELEMENT COORD. SYSTM
 
 CALL shcsgs (*610,nest(9),nest(8),nest(8),nest(21),nest(20),  &
     nest(20),necpt,tub,mcsid,thetam,tum)
 
!     BRANCH ON FORMULATION TYPE, HEAT
 
 IF (itherm /= 0) GO TO 500
 
!     FETCH MATERIAL PROPERTIES
 
 CALL gmmats (teu,3,3,0, tum,3,3,0, tem)
 CALL gmmats (tes,3,3,1, tem,3,3,0, tsm)
 CALL shgmgs (*630,elid,tsm,mid,ts,noalfa,gi,rho,gsube,tsub0, egnor,alpha)
 
!     TURN OFF THE COUPLING FLAG WHEN MID4 IS PRESENT WITH ALL
!     CALCULATED ZERO TERMS.
 
 IF (.NOT.mbcoup) GO TO 70
 DO  i = 28,36
   IF (ABS(gi(i)) > eps) GO TO 70
 END DO
 mbcoup = .false.
 70 CONTINUE
 
!     CONTINUE FILLING IN THE DATA BLOCKS
 
 ph1rst( 1) = avgthk
 ph1rst( 2) = mominr
 ph1rst(48) = offset
 relout( 2) = avgthk
 
!     PUT NORMALS IN IELOUT, GRID THICKNESS IN PH1OUT
 
 DO  i = 1,nnode
   io  = iorder(i)
   iop = (io-1)*25 + 21
   relout(iop+1) = gpnorm(2,i)
   relout(iop+2) = gpnorm(3,i)
   relout(iop+3) = gpnorm(4,i)
   ph1rst(44+io) = dgpth(i)
 END DO
 
!     CALCULATE [TSB] AND STORE IT IN IELOUT.
 
 CALL gmmats (tes,3,3,1, teb,3,3,0, tsb)
 nd25 = nnode*25
 DO  ip2 = 3,nd25,25
   DO  ix = 1,9
     relout(ip2+ix  ) = tsb(ix)
     relout(ip2+ix+9) = tsb(ix)
   END DO
 END DO
 
!     STORE ALPHA IN PH1RST(39-44)
 
 DO  ialf = 1,6
   ph1rst(38+ialf) = alpha(ialf)
 END DO
 
!     STORE UEM IN PH1RST(49-57)
!     STORE TES IN PH1RST(58-66)
 
 CALL shstts (tem,uem,vem)
 DO  ll = 1,9
   ph1rst(48+ll) = uem(ll)
   ph1rst(57+ll) = tes(ll)
 END DO
 
!     STORE THE 6X6 [G] IN PH1RST
 
 DO  ig = 3,38
   ph1rst(ig) = 0.0
 END DO
 
 IF (.NOT.membrn) GO TO 160
 DO  ig = 1,3
   ig1 = (ig-1)*6 + 2
   ig2 = (ig-1)*3
   DO  jg = 1,3
     ph1rst(ig1+jg) = gi(ig2+jg)
   END DO
 END DO
 
 160 IF (.NOT.bendng) GO TO 210
 DO  ig = 1,3
   ig1 = (ig-1)*6 + 23
   ig2 = (ig-1)*3 +  9
   DO  jg = 1,3
     ph1rst(ig1+jg) = gi(ig2+jg)*mominr
   END DO
 END DO
 
 IF (.NOT.mbcoup) GO TO 210
 DO  ig = 1,3
   ig3 = (ig-1)*3
   ig1 =  ig3 +  5
   ig2 =  ig3 + 27
   ig3 =  ig3 + 20
   DO  jg = 1,3
     ph1rst(ig1+jg) = gi(ig2+jg)
     ph1rst(ig3+jg) = gi(ig2+jg)
   END DO
 END DO
 210 CONTINUE
 
!     CALCULATE [TEG] FOR EACH NODE AND STORE IT IN PH1RST
 
 DO  i = 1,nnode
   ip = 67 + (i-1)*9
   CALL transs (igpdt(1,i),tbg)
   CALL gmmats (teb,3,3,0, tbg,3,3,0, ph1rst(ip))
 END DO
 
!     GET THE GEOMETRY CORRECTION TERMS
 
 IF (.NOT.bendng) GO TO 230
 CALL t3gems (ierr,egpdt,iorder,gi(10),gi(19),lx,ly,edglen,shrflx,  &
     aic,jog,jok,k11,k22)
 IF (ierr /= 0) GO TO 600
 
!     REDUCED INTEGRATION
 
 230 IF (INT /= 0) GO TO 260
 
!     DETERMINE THE AVERAGE [B] FOR OUT-OF-PLANE SHEAR
 
 DO  ipt = 1,3
   kpt = (ipt-1)*nd9 + 1
   CALL t3bmgs (ierr,sheart,ipt,iorder,egpdt,dgpth,aic,th,detjac,  &
       shpt,bterms,bmat1(kpt))
   IF (ierr /= 0) GO TO 600
 END DO
 
 DO  i = 1,ndof
   bmtrx(i     ) = bmat1(i+nd6) +bmat1(i+nd6+nd9) +bmat1(i+nd6+2*nd9)
   bmtrx(i+ndof) = bmat1(i+nd7) +bmat1(i+nd7+nd9) +bmat1(i+nd7+2*nd9)
 END DO
 
!     STRAIN/STRESS EVALUATION LOOP
 
!     PRESET THE PH1RST COUNTER TO THE START OF THE REPEATED SECTION
!     WHICH WILL NOW BE FILLED.
 
 260 icount =  istart
 
 DO  ipt = 4,7
   
   CALL t3bmgs (ierr,sheart,ipt,iorder,egpdt,dgpth,aic,th,detjac,  &
       shpt,bterms,bmatrx)
   IF (ierr /= 0) GO TO 600
   
   IF (INT /= 0) GO TO 310
   DO  ix = 1,ndof
     bmatrx(ix+nd6) = bmtrx(ix     )
     bmatrx(ix+nd7) = bmtrx(ix+ndof)
   END DO
   
!     FINISH FILLING IN THE DATA BLOCKS
   
!     STORE THICKNESS
   
   310 ph1rst(icount+1) = th
   
!     STORE [G3]
   
   IF (.NOT.bendng) GO TO 330
   reali = mominr*th*th*th/12.0
   tsi   = ts*th
   tss   = 1.0/(2.0*12.0*reali)
   
   zz(1) = (jog/tsi)* gi(22) + tss*jok*k22
   zz(2) =-(jog/tsi)*(gi(20) + gi(21))/2.0
   zz(3) = zz(2)
   zz(4) = (jog/tsi)* gi(19) + tss*jok*k11
   
   CALL invers (2,zz,2,bdum,0,determ,ising,INDEX)
   IF (ising /= 1) GO TO 600
   
   DO  ig = 1,4
     ph1rst(icount+1+ig) = zz(ig)
   END DO
   GO TO 350
   
   330 DO  ig = 1,4
     ph1rst(icount+1+ig) = 0.0
   END DO
   
!     STORE SHAPE FUNCTION VALUES
   
   350 DO  i = 1,nnode
     ph1rst(icount+5+i) = shpt(i)
   END DO
   
!     STORE THE STRAIN RECOVERY MATRIX
   
   DO  iph = 1,nd8
     ph1rst(icount+8+iph) = bmatrx(iph)
   END DO
   
!     END OF THE EVALUATION LOOP
   
!     INCREMENT THE PH1RST COUNTER
   
   icount = icount + 8 + 8*ndof
   
 END DO
 GO TO 700
 
 
!     BEGINNING OF HEAT FORCE RECOVERY
 
 500 CONTINUE
 
!     SET UP FOR THE UNIVERSAL PHASE 2 HEAT RECOVERY
 
 nphi(22) = 2
 nphi(23) = nnode
 nphi(24) = NAME(1)
 nphi(25) = NAME(2)
 
 sheart = .false.
 ipt = 4
 CALL t3bmgs (ierr,sheart,ipt,iorder,egpdt,dgpth,aic,th,detjac,  &
     shpt,bterms,bmatrx)
 IF (ierr /= 0) GO TO 600
 
 matid  = nest(11)
 inflag = 2
 thetas = thetas - thetam
 sinmat = SIN(thetas)
 cosmat = COS(thetas)
 CALL hmat (elid)
 
 drkce(1) = kheat(1)
 drkce(2) = kheat(2)
 drkce(3) = kheat(2)
 drkce(4) = kheat(3)
 
 tes(3) = tes(4)
 tes(4) = tes(5)
 CALL gmmats (tes,2,2,1, bterms,2,nnode,0, drkce(10))
 
 GO TO 700
 
 
!     FATAL ERRORS
 
!     CTRIA3 ELEMENT HAS ILLEGAL GEOMETRY OR CONNECTIONS
 
 600 j = 224
 GO TO 640
 
!     THE X-AXIS OF THE MATERIAL COORDINATE SYSTEM HAS NO PROJECTION
!     ON THE PLANE OF THE CTRIA3 ELEMENT
 
 610 j = 225
 nest(2) = mcsid
 GO TO 640
 
!     THE X-AXIS OF THE STRESS COORDINATE SYSTEM ID HAS NO PROJECTION
!     ON THE PLANE OF THE CTRIA3 ELEMENT
 
 620 j = 227
 nest(2) = scsid
 GO TO 640
 
!     ILLEGAL DATA DETECTED ON MATERIAL ID REFERENCED BY CTRIA3 ELEMENT
!     FOR MID3 APPLICATION
 
 630 j = 226
 nest(2) = mid(3)
 
 640 CALL mesage (30,j,nest(1))
 nogo = 1
 
 700 CONTINUE
 RETURN
END SUBROUTINE stri31
