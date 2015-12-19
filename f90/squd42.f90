SUBROUTINE squd42
     
!     PHASE 2 STRESS RECOVERY FOR 4-NODE ISOPARAMETRIC QUADRILATERAL
!     SHELL ELEMENT (QUAD4)
 
!     NOTE - FOR LAMINATED COMPOSITE ELEMENTS THE FOLLOWING ARE
!            NOT SUPPORTED
 
!         1. VARIABLE GRID POINT THICKNESS
!         3. TEMPERATURE AT 'FIBRE' DISTANCE
 
!         ALSO STRESSES ARE ONLY EVALUATED AT THE ELEMENT CENTRE
!         AND SIMILARILY FOR STRESS RESULTANTS
 
 
!     ALGORITHM -
 
!     1- STRAIN RECOVERY DATA IS SENT BY PHASE 1 THRU 'PHIOUT',
!        WHICH INCLUDES ALL THE NECESSARY TRANSFORMATIONS AND
!        STRAIN RECOVERY MATRICES. THE DATA IS REPEATED FOR EACH
!        STRESS EVALUATION POINT.
!     2- GLOBAL DISPLACEMENT VECTOR ENTERS THE ROUTINE IN CORE.
!     3- BASED ON THE DATA IN /SDR2X4/, LOCATION OF THE GLOBAL
!        DISPLACEMENT VECTOR FOR THE CURRENT SUBCASE IS DETERMINED.
!     4- WORD 132 OF /SDR2DE/ CONTAINS THE STRESS OUTPUT REQUEST
!        OPTION FOR THE CURRENT SUBCASE.
!     5- ELEMENT/GRID POINT TEMPERATURE DATA ENTERS THE ROUTINE
!        THRU /SDR2DE/ (POSITIONS 97-103, 104-129 NOT USED.)
!     6- ELEMENT STRAINS ARE CALCULATED, CORRECTED FOR THERMAL
!        STRAINS, AND PREMULTIPLIED BY G-MATRIX.
 
 EXTERNAL        andf
 LOGICAL :: extrm,layer,compos,grids,intgs,maxsh,vonms,bendng,  &
     trnflx,tempp1,tempp2,snrvrx,snrvry,four,pcmp, pcmp1,pcmp2,debug
!WKBNB NCL93012 3/94
 LOGICAL :: ostrai
 REAL :: epsavg(6)
!WKBNE NCL93012 3/94
 INTEGER :: intz(1),igrid(5),nphi(2395),nstres(86),elid,  &
     ksil(8),iorder(8),center,nfors(46),extrnl, &
!WKBR 3/95 SPR94017 2  INDXG2(3,3),INDX(6,3),OPRQST,FLAG,IPN(5),COMPS,  &
 indxg2(3,3),indx(6,3),flag,ipn(5),comps,  &
     oes1l,oef1l,pcomp,pcomp1,pcomp2,pidloc,sym,symmem,  &
     souti,fthr,strain,elemid,plyid,andf,sdest,fdest
!    5,               GPSTRS,INDEXU(3,3),INDEXV(2,3)
 REAL :: mominr,khit,mintr,tdelta(6),delta(48),tstb(5,5),  &
     tstt(5,5),tstn(50),deltat(8),u(36),g(36),g2(9),  &
     alfam(3),alfab(3),z1(5),z2(5),gpth(4),stres(86),  &
     g3(4),tmi(9),trans(9),strnt(3),strnb(3),strntc(3),  &
     strnbc(3),epst(3),epsb(3),epse(3),epstot(3),fb(2),  &
     epslne(3),stresl(3),strese(3),ezerot(6),alpha(3),  &
     v(2),ei(2),zbar(2),trnar(2),trnshr(2),ultstn(6),  &
     abbd(6,6),stiff(36),mther(6),dumc(6),stemp(8)
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm
 COMMON /zzzzzz/ z(1)
 COMMON /system/ ksystm(60)
 COMMON /sdr2c1/ ipcmp,npcmp,ipcmp1,npcmp1,ipcmp2,npcmp2, nstrop
 COMMON /sdr2x2/ dumm(30),oes1l,oef1l
 COMMON /sdr2x4/ dummy(35),ivec,ivecn,ldtemp
 COMMON /sdr2x7/ phiout(2395)
 COMMON /sdr2x8/ sigma(3),icount,nstot,thikns(5),istres,kpoint,  &
     extrnl(8),tstr(50),xpoint(2),shpfnc(4),epsln(8),  &
     khit(3),g2alfb(30),tst(20),tes(9),tesu(9),tesv(4),  &
     reali(5),GT(36),epslnt(6),tsigma(8),signx(4),  &
     signy(4),vxcntr,vycntr,fxcntr,fycntr,fxycnt,  &
     strx(2),stry(2),strs(2),forsul(46)
 COMMON /sdr2de/ ksdrde(141)
!WKBR NCL93012 3/94      COMMON /BLANK / APP(2),SORT2,IDUM(2),COMPS
 COMMON /BLANK / app(2),sort2,idum(2),comps,skp(4),ostrai
 COMMON /condas/ pi,twopi,raddeg,degrad
 EQUIVALENCE     (z(1)   ,intz(1)   ), (nfors(1) ,forsul(1) ),  &
     (nphi(1),phiout(1) ), (nstres(1),stres(1)  ),  &
     (elid   ,nphi(1)   ), (ksil(1)  ,nphi(2)   ),  &
     (tsub0  ,phiout(18)), (iorder(1),nphi(10)  ),  &
     (avgthk ,phiout(21)), (mominr   ,phiout(22)),  &
     (g(1)   ,phiout(23)), (alfam(1) ,phiout(59)),  &
     (gpth(1),phiout(65)), (alfab(1) ,phiout(62)),  &
     (ipid   ,nphi(79)  ), (kstrs    ,ksdrde(42)),  &
     (kforce ,ksdrde(41)), (stemp(1) ,ksdrde(97)),  &
     (sdest  ,ksdrde(26)), (fdest    ,ksdrde(33)),  &
     (nout   ,ksystm(2) ), (stemp(7) ,flag      )
!    1,               (INDEXU(1,1),INDEXV(1,1))
 DATA    debug / .false.  /
 DATA    center/ 4HCNTR   /
 DATA    const / 0.57735026918962/
 DATA    epss  / 1.0E-11  /
 DATA    epsa  / 1.0E-7   /
 DATA    ipn   / 1,4,2,3,5/
 DATA    pcomp / 0 /
 DATA    pcomp1/ 1 /
 DATA    pcomp2/ 2 /
 DATA    sym   / 1 /
 DATA    mem   / 2 /
 DATA    symmem/ 3 /
 DATA    strain/ 5 /
 
!     DEFINE PHIOUT(2395), THE TRANSMITTED DATA BLOCK
 
!     ADDRESS     DESCRIPTIONS
 
!        1        ELID
!      2 - 9      SIL NUMBERS
!     10 - 17     IORDER
!       18        TREF
!     19 - 20     FIBRE DISTANCES Z1, Z2 AS SPECIFIED ON PSHELL CARD
!       21        AVGTHK- AVERAGE THICKNESS OF THE ELEMENT
!       22        MOMINR- MOMENT OF INERTIA FACTOR
!     23 - 58     GBAR-MATRIX, 6X6 MATRIX OF MATERIAL PROPERTY (W/O G3)
!     59 - 61     THERMAL EXPANSION COEFFICIENTS FOR MEMBRANE
!     62 - 64     THERMAL EXPANSION COEFFICIENTS FOR BENDING
!     65 - 68     CORNER NODE THICKNESSES
!     69 - 77     TUM-MATRIX, 3X3 TRANSFORMATION FROM MATERIAL TO USER
!                 DEFINED COORDINATE SYSTEM
!       78        OFFSET OF ELEMENT FROM GP PLANE
!       79        ORIGINAL PROPERTY ID FOR COMPOSITES
!     80 - 79+9*NNODE
!                 TEG-MATRIX, A 3X3 MATRIX FOR THE TRANSFORMATION
!                 MATRIX FROM GLOBAL COORD TO ELMT COORD FOR
!                 EACH NODE.
!                 TEG-MATRIX, 3X3 DATA ARE REPEATED FOR NNODES
!     --------
!     START FROM PHIOUT(79+9*NNODE+1) AS A REFERENCE ADDRESS
!                       79+9*4    +1= 116
 
!     ADDRESS     DESCRIPTIONS
 
!        1        T, MEMBRANE THICKNESS AT THIS EVALUATION POINT
!      2 - 10     TES-MATRIX, A 3X3 TRANSFORMATION MATRIX FROM ELEM.
!                        C.S. TO USER DEFINED STRESS C.S. AT THIS
!                        EVALUATION POINT
!     11 - 19     CORRECTION TO GBAR-MATRIX FOR MEMBRANE-BENDING
!                        COUPLING AT THIS EVALUATION POINT
!     20 - 28     TMI-MATRIX, 3X3 TRANSFORMATION FROM TANGENT TO MATERIA
!     29 - 32     G3-MATRIX
!     33 - 32+NNODE
!                 ELEMENT SHAPE FUNCTION VALUES AT THIS EVAL. POINT
!     32+NNODE+1 -
!     32+NNODE+8*NDOF
!                 B-MATRIX, 8 X NDOF
 
!     --------    ABOVE DATA BATCH REPEATED 10 TIMES
 
!     TOTAL PHIOUT WORDS = (116-1) + (32+4+8*(6*4))*10
!                        =    115  + (32+4+192)*10 = 115 + 2280 = 2395
 
 
!     DEFINE STRES (TOTAL OF 86 WORDS), THE STRESS OUTPUT DATA BLOCK
 
!     ADDRESS     DESCRIPTIONS
 
!        1        ELID
!     -------------------------------------------------------
!        2        INTEGRATION POINT NUMBER
!     3  - 10     STRESSES FOR LOWER POINTS
!     11 - 18     STRESSES FOR UPPER POINTS
!     ---------   ABOVE DATA REPEATED 4 TIMES
!     70 - 86     STRESSES FOR CENTER POINT
 
!     DEFINE FORSUL (TOTAL OF 46 WORDS), THE FORCE RESULTANT OUTPUT
!     DATA BLOCK.
 
!     ADDRESS    DESCRIPTIONS
 
!        1       ELID
!     ------------------------------------------------
!        2       GRID POINT NUMBER
!      3 - 10    FORCES
!     --------   ABOVE DATA REPEATED 4 TIMES
!     38 - 46    FORCES FOR CENTER POINT
 
!     NSTOT  = NUMBER OF DATA OUTPUT THRU 'STRES'
!     NFORCE = NUMBER OF DATA OUTPUT THRU 'FORSUL'
!     NNODE  = TOTAL NUMBER OF NODES
!     NDOF   = TOTAL NUMBER OF DEGREES OF FREEDOM
!     LDTEMP = SWITCH TO DETERMINE IF THERMAL EFFECTS ARE PRESENT
!     ICOUNT = POINTER FOR PHIOUT DATA
 
!     STAGE 1 -  INITIALIZATION
!     =========================
 
!WKBNB 3/95 SPR94017
 DO  i = 1,6
   epsavg( i ) = 0.
 END DO
!WKBNE 3/95 SPR94017
 nstot = 1 + 5 + 5*2*8
 nforce= 1 + 5*9
 nnode = 0
 DO  ichk = 1,8
   IF (ksil(ichk) > 0) nnode = nnode + 1
   extrnl(ichk) = 0
 END DO
 ndof = 6*nnode
 four = nnode == 4
 
!     COMMENTS FROM G.C. 2/1990
!     EXTRNL ARE SET TO ZEROS ABOVE AND NEVER SET TO ANY VALUE LATER.
!     IT IS THEN USED TO SET IGRID. WHAT'S EXTRNL FOR?
!     THE ANSWER IS THAT EXTRNL AND IGRID ARE USED ONLY WHEN GRIDS FLAG
!     IS TRUE. GRIDS IS FALSE IN COSMIC VERSION.
 
!     ALSO, A MISSING ROUTINE, FNDGID, SUPPOSELY RETURNS EXTERNAL GRID
!     NUMBER FROM SIL INDEX. FNDGID IS LOCATED A FEW LINES BELOW 80
 
!     CHECK THE OUTPUT AND STRESS REQUEST
 
 grids = .false.
 intgs = .true.
 maxsh = andf(nstrop,1) == 0
 vonms = andf(nstrop,1) /= 0
 extrm = andf(nstrop,2) == 0
 layer = andf(nstrop,2) /= 0
 bendng= mominr > 0.0
 
!     NOTE - MAXSH AND EXTRM ARE NO LONGER USED
 
!     IF LAYERED STRESS/STARIN OUTPUT IS REQUESTED, AND THERE ARE NO
!     LAYERED COMPOSITE DATA, SET LAYER FLAG TO FALSE
 
 IF (layer .AND. npcmp+npcmp1+npcmp2 <= 0) layer = .false.
 
!     IF LAYERED OUTPUT IS REQUESTED BUT THE CURRENT ELEMENT IS NOT A
!     LAYERED COMPOSITE, SET LAYER FLAG TO FALSE
 
 IF (layer .AND. ipid < 0) layer = .false.
 
!WKBDB 3/95 SPR94017
!      OPRQST = -2
!      IF (KSTRS  .EQ. 1) OPRQST = OPRQST + 1
!      IF (KFORCE .EQ. 1) OPRQST = OPRQST + 2
!WKBI NCL93012 3/94
!      IF ( OSTRAI ) OPRQST = OPRQST + 1
!      IF (OPRQST .EQ.-2) RETURN
!WKBDE 3/95 SPR94017
!WKBI  3/95 SPR94017
 IF ( ( kstrs  /= 1 ) .AND. ( kforce /= 1 ) .AND.  &
     (.NOT.ostrai)          )RETURN
 
!     CHECK FOR FIBRE DISTANCES Z1 AND Z2 BEING BLANK
 
 logz12 = -4
 IF (nphi(19) == 0) logz12 = logz12 + 2
 IF (nphi(20) == 0) logz12 = logz12 + 4
 
!     CHECK FOR THE TYPE OF TEMPERATURE DATA
!     NOTES  1- TYPE TEMPP1 ALSO INCLUDES TYPE TEMPP3
!            2- IF NIETHER TYPE IS TRUE, GRID POINT TEMPERATURES
!               ARE PRESENT.
 
 tempp1 = flag == 13
 tempp2 = flag ==  2
 
!     CHECK FOR OFFSET AND COMPOSITES
 
 offset = phiout(78)
 compos = comps == -1 .AND. ipid > 0
 
!     ZERO OUT STRESS AND FORCE RESULTANT ARRAYS
 
 DO  k = 1,nstot
   stres(k) = 0.0
 END DO
 DO  i = 1,nforce
   forsul(i)= 0.0
 END DO
 nstres(1)= elid
 nfors(1) = elid
 
!     ZERO OUT THE COPY OF GBAR-MATRIX TO BE USED BY THIS ROUTINE
 
 DO  k = 1,36
   GT(k) = 0.0
 END DO
 
!     STAGE 2 - ARRANGEMENT OF INCOMING DATA
!     ======================================
 
!     SORT THE GRID TEMPERATURE CHANGES INTO SIL ORDER (IF PRESENT)
 
 IF (ldtemp ==     -1) GO TO 60
 IF (tempp1 .OR. tempp2) GO TO 60
 
!     DO 50 K = 1,NNODE
!     KPOINT = IORDER(K)
!  50 DELTAT(K) = STEMP(KPOINT)
 
!     COMMENTS FORM G.CHAN/UNISYS  2/93
!     THE ABOVE DO 50 LOOP DOES NOT WORK SINCE STEMP(2 THRU NNODE) = 0.0
 
 DO  k = 1,nnode
   deltat(k) = stemp(1)
 END DO
 
!     PICK UP THE GLOBAL DISPLACEMENT VECTOR AND TRANSFORM IT
!     INTO THE ELEMENT C.S.
 
 60 DO  idelt = 1,nnode
   jdelt = ivec + ksil(idelt) - 2
   kdelt = 6*(idelt-1)
   DO  ldelt = 1,6
     tdelta(ldelt) = z(jdelt+ldelt)
   END DO
   
!     FETCH TEG-MATRIX 3X3 FOR EACH NODE AND LOAD IT IN A 6X6 MATRIX
!     INCLUDE THE EFFECTS OF OFFSET
   
   CALL tldrs  (offset,idelt,phiout(80),u)
   CALL gmmats (u,6,6,0, tdelta,6,1,0, delta(kdelt+1))
 END DO
 
!     GET THE EXTERNAL GRID POINT ID NUMBERS FOR CORRESPONDING SIL
!     NUMBERS.
 
!     CALL FNDGID (ELID,8,KSIL,EXTRNL)
 
!     STAGE 3 - CALCULATION OF STRAINS
!     ================================
 
!     INTEGRATION DATA IN PHIOUT IS ARRANGED IN ETA, XI INCREASING
!     SEQUENCE.
 
 isig  = 1
 icount= -(8*ndof+nnode+32) + 79 + 9*nnode
 
 DO  inplan = 1,5
   inpln1 = ipn(inplan)
   
!     MATCH GRID ID NUMBER WHICH IS IN SIL ORDER
   
   IF (inplan == 5) GO TO 100
   DO  i = 1,nnode
     IF (iorder(i) /= inpln1) CYCLE
     igrid(inplan) = extrnl(i)
     GO TO 110
   END DO
   GO TO 110
   
   100 igrid(inplan) = center
   110 CONTINUE
   
   DO  izta = 1,2
     zeta = (izta*2-3)*const
     
     icount = icount + 8*ndof + nnode + 32
     IF (izta == 2) GO TO 160
     
!     THICKNESS AND MOMENT OF INERTIA AT THIS POINT
     
     thikns(inplan) = phiout(icount+1)
     IF (grids .AND. inplan /= 5) thikns(inplan) = gpth(inpln1)
     reali(inplan) = mominr*thikns(inplan)**3/12.0
     
!     DETERMINE FIBER DISTANCE VALUES
     
     IF (logz12 == -4) GO TO 150
     IF (logz12 < 0) THEN
       GO TO   120
     ELSE IF (logz12 == 0) THEN
       GO TO   130
     ELSE
       GO TO   140
     END IF
     
     120 z1(inplan) = -0.5*thikns(inplan)
     z2(inplan) = phiout(20)
     GO TO 160
     
     130 z1(inplan) = phiout(19)
     z2(inplan) = 0.5*thikns(inplan)
     GO TO 160
     
     140 z1(inplan) = -0.5*thikns(inplan)
     z2(inplan) = -z1(inplan)
     GO TO 160
     
     150 z1(inplan) = phiout(19)
     z2(inplan) = phiout(20)
     160 CONTINUE
     
!     FIRST COMPUTE LOCAL STRAINS UNCORRECTED FOR THERMAL STRAINS
!     AT THIS EVALUATION POINT.
     
!        EPSLN  = PHIOUT(KSIG) * DELTA
!          EPS  =       B      *   U
!          8X1        8XNDOF    NDOFX1
     
     ksig = icount+nnode+33
     CALL gmmats (phiout(ksig),8,ndof,0, delta(1),ndof,1,0, epsln)
     
!     CALCULATE THERMAL STRAINS IF TEMPERATURES ARE PRESENT
     
     IF (ldtemp == -1) GO TO 260
     DO  iet = 1,6
       epslnt(iet) = 0.0
     END DO
     
!     A) MEMBRANE STRAINS
     
     IF (tempp1 .OR. tempp2) GO TO 190
     
!     GRID TEMPERATURES
     
     kshp = icount + 32
     tbar = 0.0
     DO  ish = 1,nnode
       ksh  = kshp + ish
       tbar = tbar + phiout(ksh)*deltat(ish)
     END DO
     tmean= tbar
     GO TO 200
     
!     ELEMENT TEMPERATURES
     
     190 tbar = stemp(1)
     200 tbar = tbar - tsub0
     DO  ieps = 1,3
       epslnt(ieps) = -tbar*alfam(ieps)
     END DO
     
!     B) BENDING STRAINS (ELEMENT TEMPERATURES ONLY)
     
     IF (.NOT.bendng) GO TO 260
     IF (.NOT.(tempp1 .OR. tempp2)) GO TO 260
     
!     EXTRACT G2-MATRIX FROM GBAR-MATRIX AND CORRECT IT FOR COUPLING
     
     ig21 = 0
     DO  ig2 = 1,3
       ig22 = (ig2-1)*6 + 21
       DO  jg2 = 1,3
         ig21 = ig21 + 1
         jg22 = jg2  + ig22
         g2(ig21) = g(jg22) + phiout(icount+10+ig21)
       END DO
     END DO
     
     ig2ab = (isig*3)/5 + 1
     CALL gmmats (g2,3,3,0, alfab,3,1,0, g2alfb(ig2ab))
     
     IF (tempp1) GO TO 240
     CALL invers (3,g2,3,gdum,0,detg2,isngg2,indxg2)
     CALL gmmats (g2,3,3,0, stemp(2),3,1,0, khit)
     DO  ieps = 4,6
       epslnt(ieps) = khit(ieps-3)*zeta*thikns(inplan)/(2.*reali(inplan))
     END DO
     GO TO 260
     
     240 tprime = stemp(2)
     DO  ieps = 4,6
       epslnt(ieps) = -tprime*alfab(ieps-3)*zeta*thikns(inplan)/2.
     END DO
     
!     MODIFY GBAR-MATRIX
     
     260 i1 = -6
     i2 = 12
     i3 = 11 + icount
     DO  i = 1,3
       i1 = i1 + 6
       i2 = i2 + 6
       DO  j = 1,3
         j1 = j  + i1
         j3 = j1 + 3
         j4 = j  + i2
         j2 = j4 + 3
         GT(j1) = g(j1)
         GT(j2) = g(j2)
         GT(j3) = g(j3) + phiout(i3)
         GT(j4) = g(j4) + phiout(i3)
         i3 = i3 + 1
       END DO
     END DO
     
!     DETERMINE G MATRIX FOR THIS EVALUATION POINT
     
     DO  i = 1,4
       g3(i) = phiout(icount+28+i)
     END DO
     
     IF (ldtemp == -1) GO TO 300
     
!     CORRECT STRAINS FOR THERMAL EFFECTS
     
     DO  i = 1,6
       epsln(i) = epsln(i) + epslnt(i)
     END DO
     
!     CALCULATE STRESS VECTOR
     
     300 CALL gmmats (GT(1),6,6,0, epsln(1),6,1,0, tsigma(1))
     CALL gmmats (g3(1),2,2,0, epsln(7),2,1,0, tsigma(7))
!WKBNB NCL93012 3/94
     IF ( izta /= 1 ) GO TO 303
     DO  iav = 1, 3
       epsavg(iav) = epsavg(iav) + epsln(iav)
     END DO
     DO  iav = 4, 6
       epsavg(iav) = epsavg(iav) + epsln(iav) / const
     END DO
     303   CONTINUE
!WKBNE NCL93012 3/94
     IF (.NOT.bendng) GO TO 320
     
!     COMBINE STRESSES ONLY IF 'BENDING'
     
     DO  i = 1,3
       tsigma(i) = tsigma(i+3)
     END DO
     
     320 CONTINUE
     
!     TRANSFORM STRESSES FROM ELEMENT TO STRESS C.S.
     
     DO  i = 1,9
       tes(i) = phiout(icount+1+i)
     END DO
     
     tesu(1) = tes(1)*tes(1)
     tesu(2) = tes(4)*tes(4)
     tesu(3) = tes(1)*tes(4)
     tesu(4) = tes(2)*tes(2)
     tesu(5) = tes(5)*tes(5)
     tesu(6) = tes(2)*tes(5)
     tesu(7) = tes(1)*tes(2)*2.0
     tesu(8) = tes(4)*tes(5)*2.0
     tesu(9) = tes(1)*tes(5) + tes(2)*tes(4)
     
     CALL gmmats (tesu(1),3,3,1, tsigma(1),3,1,0, tstr(isig))
     
     tesv(1) = tes(5)*tes(9) + tes(6)*tes(8)
     tesv(2) = tes(2)*tes(9) + tes(8)*tes(3)
     tesv(3) = tes(4)*tes(9) + tes(7)*tes(6)
     tesv(4) = tes(1)*tes(9) + tes(3)*tes(7)
     
     isig = isig + 3
     CALL gmmats (tesv(1),2,2,1, tsigma(7),2,1,0, tstr(isig))
     
     isig = isig + 2
   END DO
 END DO
 
!     IF REQUIRED, EXTRAPOLATE STRESSES FROM INTEGRATION POINTS
!     TO CORNER POINTS.
 
!     FIRST EXTRAPOLATE ACROSS ZETA, REGARDLESS OF INPLANE REQUEST
 
 DO  ikk = 1,5
   itb = (ikk-1)*10
   DO  ijj = 1,5
     tstb(ikk,ijj) = tstr(itb+  ijj)
     tstt(ikk,ijj) = tstr(itb+5+ijj)
   END DO
 END DO
 
 x1 = -const
 x2 = -x1
 
 DO  k = 1,2
   ik = 0
   xx = -1.0
   IF (k == 2) xx =-xx
   IF (k == 2) ik = 5
   
   xn22 = (xx-x1)/(x2-x1)
   xn11 = 1.0 - xn22
   
   DO  i = 1,5
     ikkn = (i-1)*10 + ik
     DO  j = 1,5
       tstn(ikkn+j) = tstb(i,j)*xn11 + tstt(i,j)*xn22
     END DO
   END DO
 END DO
 
 DO  ii = 1,50
   tstr(ii) = tstn(ii)
 END DO
 
 IF (intgs .OR. compos) GO TO 540
 
 ixtr = 5
 jxtr = ixtr*4
 
 iz1 = 0
 DO  iz = 1,2
   
   DO  i = 1,jxtr
     tst(i) = 0.0
   END DO
   
!     FOR THE SAKE OF COMPATIBILITY BETWEEN THE CONVENTION FOR
!     SHEAR FORCES, AND THE CONVENTION FOR EXTRAPOLATION, WE MAY
!     HAVE TO CHANGE THE SIGNS AROUND FOR SPECIFIC POINTS. THEY
!     WILL BE RETURNED TO THE ORIGINAL SIGNS AFTER EXTRAPOLATION IS
!     COMPLETE.
   
!WKBR 3/95 SPR94017      IF (OPRQST .LT. 0) GO TO 460
   IF ( kforce /= 1 ) GO TO 460
   DO  i = 1,4
     j = (i-1)*2*ixtr + iz1 + 4
     IF (tstr(j) == 0.0) GO TO 410
     signy(i) = tstr(j)/ABS(tstr(j))
     GO TO 420
     410 signy(i) = 0.0
     420 IF (tstr(j+1) == 0.0) GO TO 430
     signx(i) = tstr(j+1)/ABS(tstr(j+1))
     CYCLE
     430 signx(i) = 0.0
   END DO
   
   snrvry = .false.
   IF (signy(1)*signy(2) <= 0.0 .OR. signy(3)*signy(4) <= 0.0 .OR.  &
       signy(3)*signy(1) <= 0.0) snrvry = .true.
   snrvrx = .false.
   IF (signx(1)*signx(2) <= 0.0 .OR. signx(3)*signx(4) <= 0.0 .OR.  &
       signx(3)*signx(1) <= 0.0) snrvrx = .true.
   
   IF (.NOT.snrvry) GO TO 450
   tstr(iz1+4) = -tstr(iz1+4)
   tstr(iz1+4+4*ixtr) = -tstr(iz1+4+4*ixtr)
   450 IF (.NOT.snrvrx) GO TO 460
   tstr(iz1+5) = -tstr(iz1+5)
   tstr(iz1+5+2*ixtr) = -tstr(iz1+5+2*ixtr)
   460 CONTINUE
   
   xpoint(1) = -1.0
   xpoint(2) = +1.0
   ir = 0
   
   DO  ix = 1,2
     xi = xpoint(ix)
     
     DO  ie = 1,2
       eta = xpoint(ie)
       
       shpfnc(1) = 0.75*(const-xi)*(const-eta)
       shpfnc(2) = 0.75*(const-xi)*(const+eta)
       shpfnc(3) = 0.75*(const+xi)*(const-eta)
       shpfnc(4) = 0.75*(const+xi)*(const+eta)
       
       li = ir*ixtr
       ir = ir + 1
       
       DO  is = 1,4
         lk = (is-1)*2*ixtr + iz1
         
         DO  it = 1,ixtr
           tst(li+it) = tst(li+it) + shpfnc(is)*tstr(lk+it)
         END DO
       END DO
     END DO
   END DO
   
   j1 = 0
   DO  is = 1,4
     j2 = (is-1)*2*ixtr + iz1
     DO  js = 1,ixtr
       j1 = j1 + 1
       j2 = j2 + 1
       tstr(j2) = tst(j1)
     END DO
   END DO
   
!     CHANGE THE SIGNS BACK, IF NECESSARY
   
!WKBR 3/95 SPR94017     IF (OPRQST .LT. 0) GO TO 520
   IF ( kforce /= 1 ) GO TO 520
   IF (.NOT.snrvry) GO TO 510
   tstr(iz1+4) = -tstr(iz1+4)
   tstr(iz1+4+4*ixtr) = -tstr(iz1+4+4*ixtr)
   510 IF (.NOT.snrvrx) GO TO 520
   tstr(iz1+5) = -tstr(iz1+5)
   tstr(iz1+5+2*ixtr) = -tstr(iz1+5+2*ixtr)
   520 CONTINUE
   iz1 = iz1 + ixtr
 END DO
 540 CONTINUE
 
!     STAGE 4 - CALCULATION OF OUTPUT STRESSES
!     ========================================
 
!WKBR 3/95 SPR94017     IF (OPRQST .EQ. 0) GO TO 740
 IF ( (kstrs /= 1) .AND. (.NOT. ostrai) ) GO TO 740
 
!WKBNB NCL93012 3/94
 DO  iav = 1, 3
   epsavg(iav) = epsavg(iav) / 5.
 END DO
 DO  iav = 4, 6
   epsavg(iav) = epsavg(iav) / ( 5. * phiout(21)/2. )
 END DO
!WKBNE NCL93012 3/94
 isig    = 0
 ig2a    = 0
 strx(1) = 0.0
 strx(2) = 0.0
 stry(1) = 0.0
 stry(2) = 0.0
 strs(1) = 0.0
 strs(2) = 0.0
 DO  inplan = 1,5
   inpln1 = inplan
   IF (inplan == 2) inpln1 = 4
   IF (inplan == 3) inpln1 = 2
   IF (inplan == 4) inpln1 = 3
   
   istres = (inpln1-1)*17 + 2
   
   idpont = igrid(inplan)
   IF (intgs) idpont = inpln1
   IF (intgs .AND. inplan == 5) idpont = center
   nstres(istres) = idpont
   thick = thikns(inplan)
   
   DO  iz = 1,2
     IF (iz == 2) istres = istres + 8
     fibre = z1(inplan)
     IF (iz == 2) fibre = z2(inplan)
!WKBNB NCL93012 3/94
     IF ( .NOT. ostrai ) GO TO 545
     IF ( iz /= 1 ) GO TO 542
     nstres( istres+1 ) = 0
     sigma( 1 ) = epsavg( 1 )
     sigma( 2 ) = epsavg( 2 )
     sigma( 3 ) = epsavg( 3 )
     GO TO 630
     542   CONTINUE
     nstres( istres+1 ) = -1
     sigma( 1 ) = epsavg( 4 )
     sigma( 2 ) = epsavg( 5 )
     sigma( 3 ) = epsavg( 6 )
     GO TO 630
     545   CONTINUE
!WKBNE NCL93012 3/94
     stres(istres+1) = fibre
     
!     EVALUATE STRESSES AT THIS FIBRE DISTANCE
     
     DO  i = 1,3
       sigma(i) = (0.5-fibre/thick)*tstr(isig+i) + (0.5+fibre/thick)  &
           *tstr(isig+i+5)
     END DO
     
!     IF TEMPERATURES ARE PRESENT, CORRECT STRESSES FOR THERMAL
!     STRESSES ASSOCIATED WITH THE DATA RELATED TO FIBRE DISTANCES.
     
     IF (ldtemp == -1) GO TO 610
     
!     IF NO BENDING, TREAT IT LIKE GRID POINT TEMPERATURES
     
     IF (.NOT.bendng) GO TO 610
     IF (tempp1) GO TO 560
     IF (tempp2) GO TO 570
     GO TO 610
     
     560 tsubi = stemp(2+iz)
     IF (ABS(tsubi) < epss) GO TO 610
     tsubi = tsubi - tprime*fibre
     GO TO 590
     
     570 tsubi = stemp(4+iz)
     IF (ABS(tsubi) < epss) GO TO 610
     DO  ist = 1,3
       sigma(ist) = sigma(ist) - stemp(ist+1)*fibre/reali(inplan)
     END DO
     590 tsubi = tsubi - tbar
     DO  its = 1,3
       sigma(its) = sigma(its) - tsubi*g2alfb(ig2a+its)
     END DO
     
!     AVERAGE THE VALUES FROM OTHER 4 POINTS FOR THE CENTER POINT
     
     610 IF (inplan == 5) GO TO 620
     strx(iz) = strx(iz) + 0.25*sigma(1)
     stry(iz) = stry(iz) + 0.25*sigma(2)
     strs(iz) = strs(iz) + 0.25*sigma(3)
     GO TO 630
     620 sigma(1) = strx(iz)
     sigma(2) = stry(iz)
     sigma(3) = strs(iz)
     630 DO  is = 1,3
       stres(istres+1+is) = sigma(is)
     END DO
     
!     CALCULATE PRINCIPAL STRESSES
     
     sigavg = 0.5*(sigma(1) + sigma(2))
     proj   = 0.5*(sigma(1) - sigma(2))
     taumax = proj*proj + sigma(3)*sigma(3)
!WKBNB 7/94 SPR94004
     IF ( .NOT. ostrai ) GO TO 645
     taumax = proj*proj + sigma(3)*sigma(3)/4.
     GO TO 649
     645   CONTINUE
!WKBNE 7/94 SPR94004
     IF (ABS(taumax) <= epss) GO TO 650
!WKBI  7/94 SPR94004
     649   CONTINUE
     taumax = SQRT(taumax)
     GO TO 660
     650 taumax = 0.0
     
!     PRINCIPAL ANGLE
     
     660 txy2 = sigma(3)*2.0
     proj = proj*2.0
     IF (ABS(txy2) <= epsa .AND. ABS(proj) <= epsa) GO TO 670
     stres(istres+5) = 28.647890*ATAN2(txy2,proj)
     GO TO 680
     670 stres(istres+5) = 0.0
     680 sigma1 = sigavg + taumax
     sigma2 = sigavg - taumax
     stres(istres+6) = sigma1
     stres(istres+7) = sigma2
     
!     OUTPUT VON MISES YIELD STRESS IF ASKED FOR BY THE USER
     
     IF (vonms) GO TO 690
     stres(istres+8) = taumax
!WKBI NCL93012 3/94
     IF ( ostrai ) stres(istres+8) = 2.*taumax
     GO TO 720
     
     690 sigyp = sigma1*sigma1 + sigma2*sigma2 - sigma1*sigma2
     IF (ABS(sigyp) <= epss) GO TO 700
     sigyp = SQRT(sigyp)
     GO TO 710
     700 sigyp = 0.0
     710 stres(istres+8) = sigyp
     
     720 ig2a = ig2a + 3
   END DO
   isig = isig + 10
 END DO
!WKBNB NCL93012 3/94
 DO  iav = 1, 6
   epsavg( iav ) = 0.
 END DO
!WKBNE NCL93012 3/94
 
!     STAGE 5 - ELEMENT FORCE OUTPUT
!     ==============================
 
 740 IF (layer) GO TO 750
!WKBR 3/95 SPR94017     IF (OPRQST .LT. 0) GO TO 790
 IF ( kforce /= 1 ) GO TO 790
 
 750 CONTINUE
 isig   = 0
 vxcntr = 0.0
 vycntr = 0.0
 fxcntr = 0.0
 fycntr = 0.0
 fxycnt = 0.0
 DO  inplan = 1,5
   inpln1 = inplan
   IF (inplan == 2) inpln1 = 4
   IF (inplan == 3) inpln1 = 2
   IF (inplan == 4) inpln1 = 3
   thick = thikns(inplan)
   
   iforce = (inpln1-1)*9 + 2
   
   idpont = igrid(inplan)
   IF (intgs) idpont = inpln1
   IF (intgs .AND. inplan == 5) idpont = center
   nfors(iforce) = idpont
   
!     CALCULATE FORCES AT MID-SURFACE LEVEL
   
   DO  ifor = 1,3
     forsul(iforce+ifor  )=(tstr(isig+ifor)+tstr(isig+ifor+5))*thick/2.
     forsul(iforce+ifor+3)=(tstr(isig+ifor)-tstr(isig+ifor+5))*  &
         reali(inplan)/thick
   END DO
   
!     INTERCHANGE 7 AND 8 POSITIONS TO BE COMPATIBLE WITH THE
!     OUTPUT FORMAT OF VX AND VY (WE HAVE CALCULATED VY AND VX)
   
   IF (inplan == 5) GO TO 770
   forsul(iforce+7) = (tstr(isig+5) + tstr(isig+10))*thick*0.5
   forsul(iforce+8) = (tstr(isig+4) + tstr(isig+ 9))*thick*0.5
   
!     SUBSTITUTE THE AVERAGE OF CORNER (OR INTEGRATION) POINT
!     MEMBRANE AND SHEAR FORCES FOR THE CENTER POINT
   
   fxcntr = fxcntr + forsul(iforce+1)*0.25
   fycntr = fycntr + forsul(iforce+2)*0.25
   fxycnt = fxycnt + forsul(iforce+3)*0.25
   vxcntr = vxcntr + forsul(iforce+7)*0.25
   vycntr = vycntr + forsul(iforce+8)*0.25
   GO TO 780
   770 CONTINUE
   forsul(iforce+1) = fxcntr
   forsul(iforce+2) = fycntr
   forsul(iforce+3) = fxycnt
   forsul(iforce+7) = vxcntr
   forsul(iforce+8) = vycntr
   
   780 isig = isig + 10
 END DO
 
!     DO NOT WRITE TO PHIOUT IF LAYER STRESSES ARE REQUESTED
!     BECAUSE PHIOUT NEEDS TO BE INTACT
 IF (layer) GO TO 900
 
!     STAGE 7 - SHIPPING OF NORMAL STRESSES
!     =====================================
 
!     STORE THE STRESSES WHERE THE HIGHER LEVEL ROUTINES EXPECT
!     TO FIND THEM.
!     BUT FIRST, MOVE THE CENTER POINT STRESSES TO THE TOP.
 
!WKBR 3/95 SPR94017     IF (OPRQST .EQ. 0) GO TO 840
 IF ( (kstrs /= 1) .AND. (.NOT.ostrai) ) GO TO 840
 790 nphi(101) = nstres(1)
 DO  i = 3,18
   i99 = i + 99
   nphi(i99) = nstres(i+68)
 END DO
 
!     DEBUG PRINTOUT
 
 IF (debug) WRITE (nout,810) (stres(i),i=71,86)
 810 FORMAT (' SQUD42 - STRESSES', (/1X,8E13.5))
 
 DO  i = 19,86
   i99 = i + 99
   nphi(i99) = nstres(i-17)
 END DO
 
!     STORE FORCES IN THEIR APPROPRIATE LOCATION
 
!WKBR 3/95 SPR94017     IF (OPRQST .LT. 0) RETURN
 IF ( kforce /= 1 ) RETURN
 840 nphi(201) = nfors(1)
 DO  i = 3,10
   i199 = i + 199
   nphi(i199) = nfors(i+36)
 END DO
 
!     DEBUG PRINTOUT
 
 IF (debug) WRITE (nout,860) (forsul(i),i=39,46)
 860 FORMAT (' SQUD42 - FORCES', (/1X,8E13.5))
 
 DO  i = 11,46
   i199 = i + 199
   nphi(i199) = nfors(i-9)
 END DO
 
!     PROCESSING FOR NORMAL STRESS REQUEST COMPLETED
 
 GO TO 2100
 
!     ELEMENT LAYER STRESS CALCULATION
 
!     CHECK STRESS AND FORCE OUTPUT REQUEST
 
 900 IF ((kforce /= 0 .OR. kstrs /= 0) .AND. .NOT.compos) GO TO 2220
 
!     WRITE FORCE RESULTANTS TO OEF1L IF REQUESTED
!         1.    10*ELEMENT ID + DEVICE CODE (FDEST)
!        2-9.   FORCE RESULTANTS
!               FX, FY, FXY, MX, MY, MXY, VX, VY
 
 IF (kforce ==  0) GO TO 910
 elemid = 10*elid + fdest
 IF (ldtemp /= -1) GO TO 910
 CALL WRITE (oef1l,elemid,1,0)
 CALL WRITE (oef1l,forsul(39),8,0)
 
 910 IF (kstrs == 0 .AND. ldtemp == -1) RETURN
 elemid = 10*elid + sdest
 
!     LOCATE PID BY CARRYING OUT A SEQUENTIAL SEARCH
!     OF THE PCOMPS DATA BLOCK, AND ALSO DETERMINE
!     THE TYPE OF 'PCOMP' BULK DATA ENTRY.
 
!     SET POINTER LPCOMP
 
 lpcomp = ipcmp + npcmp + npcmp1 + npcmp2
 
 
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
!            = MEM  (MEMBRANE )
!            = SYMMEM  (SYMMETRIC-MEMBRANE)
 
!     FTHR   - FAILURE THEORY
!            = 1    HILL
!            = 2    HOFFMAN
!            = 3    TSAI-WU
!            = 4    MAX-STRESS
!            = 5    MAX-STRAIN
 
!     ULTSTN - ULTIMATE STRENGTH VALUES
 
!     SET POINTERS
 
 itype = -1
 
 pcmp  = .false.
 pcmp1 = .false.
 pcmp2 = .false.
 
 pcmp  = npcmp  > 0
 pcmp1 = npcmp1 > 0
 pcmp2 = npcmp2 > 0
 
!     CHECK IF NO 'PCOMP' DATA HAS BEEN READ INTO CORE
 
 IF (.NOT.pcmp .AND. .NOT.pcmp1 .AND. .NOT.pcmp2) GO TO 2200
 
!     SEARCH FOR PID IN PCOMP DATA
 
 IF (.NOT.pcmp) GO TO 960
 
 ip = ipcmp
 IF (intz(ip) == ipid) GO TO 950
 ipc11 = ipcmp1 - 1
 DO  ip = ipcmp,ipc11
   IF (intz(ip) == -1 .AND. ip < (ipcmp1-1)) GO TO 920
   CYCLE
   920 IF (intz(ip+1) == ipid) GO TO 940
 END DO
 GO TO 960
 
 940 ip = ip + 1
 950 itype = pcomp
 GO TO 1070
 
!     SEARCH FOR PID IN PCOMP1 DATA
 
 960 IF (.NOT.pcmp1) GO TO 1010
 ip = ipcmp1
 IF (intz(ip) == ipid) GO TO 1000
 ipc21 = ipcmp2 - 1
 DO  ip = ipcmp1,ipc21
   IF (intz(ip) == -1 .AND. ip < (ipcmp2-1)) GO TO 970
   CYCLE
   970 IF (intz(ip+1) == ipid) GO TO 990
 END DO
 GO TO 1010
 
 990 ip = ip + 1
 1000 itype = pcomp1
 GO TO 1070
 
!     SEARCH FOR PID IN PCOMP2 DATA
 
 1010 IF (.NOT.pcmp2) GO TO 1060
 
 ip = ipcmp2
 IF (intz(ip) == ipid) GO TO 1050
 lpc11 = lpcomp - 1
 DO  ip = ipcmp2,lpc11
   IF (intz(ip) == -1 .AND. ip < (lpcomp-1)) GO TO 1020
   CYCLE
   1020 IF (intz(ip+1) == ipid) GO TO 1040
 END DO
 GO TO 1060
 
 1040 ip = ip + 1
 1050 itype = pcomp2
 GO TO 1070
 
!     CHECK IF PID HAS NOT BEEN LOCATED
 
 1060 IF (itype == -1) GO TO 2200
 
!     LOCATION OF PID
 
 1070 pidloc = ip
 lamopt = intz(pidloc+8)
 
!     INTILIZE
 
 DO  ir = 1,3
   strnt(ir) = 0.0
   strnb(ir) = 0.0
 END DO
 
!     CALCULATION OF STRAINS
 
!     INTEGRATION DATA IN PHIOUT IS ARRANGED IN ETA,XI INCREASING
!     SEQUENCE.
 
 isig   = 1
 icount = -(8*ndof+nnode+32) + 79 + 9*nnode
 
 DO  inplan = 1,5
   inpln1 = ipn(inplan)
   
!     MATCH GRID ID NUMBER WHICH IS IN SIL ORDER
   
   IF (inplan == 5) GO TO 1100
   DO  i = 1,nnode
     IF (iorder(i) /= inpln1) CYCLE
     igrid(inplan) = extrnl(i)
     GO TO 1110
   END DO
   GO TO 1110
   
   1100 igrid(inplan) = center
   1110 CONTINUE
   
   DO  izta = 1,2
     zeta = (izta*2-3)*const
     
     icount = icount + 8*ndof + nnode + 32
     
!     FIRST COMPUTE LOCAL STRAINS AT THIS EVALUATION POINT
     
!        EPSLN = PHIOUT(KSIG) * DELTA
!          EPS =        B     *   U
!          8X1        8XNDOF    NDOFX1
     
     ksig = icount + nnode + 33
     CALL gmmats (phiout(ksig),8,ndof,0, delta(1),ndof,1,0, epsln)
     
!     TRANSFORM THE STRAINS AT THIS EVALUATION POINT TO THE
!     MATERIAL COORDINATE SYSTEM
     
     DO  ir = 1,9
       tmi(ir) = phiout(icount+19+ir)
     END DO
     
!     TOTAL STRAIN AT EVALUATION POINT = MEMBRANE + BENDING
     
     DO  ir = 1,3
       epstot(ir) = epsln(ir) + epsln(ir+3)
     END DO
     
!     GENERATE TRANS-MATRIX TO TRANSFORM STRAINS FROM I TO M SYSTEM
     
     trans(1)  = tmi(1)*tmi(1)
     trans(2)  = tmi(2)*tmi(2)
     trans(3)  = tmi(1)*tmi(2)
     trans(4)  = tmi(4)*tmi(4)
     trans(5)  = tmi(5)*tmi(5)
     trans(6)  = tmi(4)*tmi(5)
     trans(7)  = 2.0*tmi(1)*tmi(4)
     trans(8)  = 2.0*tmi(2)*tmi(5)
     trans(9)  = tmi(1)*tmi(5) + tmi(2)*tmi(4)
     
!     TRANSFORM TOTAL STRAINS
     
     CALL gmmats (trans(1),3,3,0, epstot(1),3,1,0, epse(1))
     
     IF (inplan == 5) GO TO 1160
     
!     AVERAGE THE STRAIN VECTORS OF THE FOUR INTGS POINTS AT EACH
!     LEVEL TO CALCULATE THE ELEMENT CENTRE STRAIN VECTOR FOR THE
!     UPPER AND BOTTOM LEVELS.
     
     DO  ir = 1,3
       IF (izta == 2) GO TO 1140
       strnb(ir) = strnb(ir) + 0.25*epse(ir)
       CYCLE
       1140 strnt(ir) = strnt(ir) + 0.25*epse(ir)
     END DO
     CYCLE
     
!     TOTAL STRAIN VECTORS AT ELEMENT CENTRE
     
     1160 DO  ir = 1,3
       IF (izta == 2) GO TO 1170
       strnbc(ir) = epse(ir)
       CYCLE
       1170 strntc(ir) = epse(ir)
     END DO
     
   END DO
 END DO
 
!     EXTRAPOLATE STRAINS ACROSS ZETA
 
 DO  ir = 1,3
   epst(ir) = (strnt(ir)-strnb(ir))*(+1.0+const)/(2.0*const) +  strnb(ir)
   epsb(ir) = (strnt(ir)-strnb(ir))*(-1.0+const)/(2.0*const) +  strnb(ir)
 END DO
 
!     CALCULATE LAYER STRESSES AND FAILURE INDICES (IF REQUESTED)
!     AND WRITE TO THE OUTPUT FILE OES1L
!         1.    10*ELEMENT ID + DEVICE CODE (SDEST)
!         2.    NLAYER - NUMBER OF LAYERS FOR LAMINATE
!         3.    TYPE OF FAILURE THEORY SELECTED
 
!         4.    PLY ID
!       5,6,7.  LAYER STRESSES
!         8.    PLY FAILURE INDEX (FP)
!         9.    IFLAG (= 1 IF FP.GE.0.999, DEFAULT = 0)
!       10,11.  INTERLAMINAR SHEAR STRESSES
!        12.    SHEAR BONDING INDEX (SB)
!        13.    IFLAG (= 1 IF SB.GE.0.999, DEFAULT = 0)
!         :     4 - 13 REPEATED FOR THE NUMBER OF LAYERS WITH
!         :           LAYER STRESS REQUEST
!      LAST-1.  MAXIMUM FAILURE INDEX OF LAMINATE  (FIMAX)
!       LAST.   IFLAG (= 1 IF FIMAX.GE.0.999, DEFAULT = 0)
 
!      1-LAST.  REPEAT FOR NUMBER OF ELEMENTS
 
!       (NOTE - ONLY THE ELEMENT CENTRE VALUES ARE CALCULATED)
 
!     == 1.
 
 IF (kstrs == 1) CALL WRITE (oes1l,elemid,1,0)
 
!     DETERMINE INTRINSIC LAMINATE PROPERTIES
 
!     LAMINATE THICKNESS
 
 tlam = phiout(21)
 
!     REFERENCE SURFACE
 
 zref = -tlam/2.0
 
!     NUMBER OF LAYERS
 
 nlay = intz(pidloc+1)
 
!     FOR PCOMP BULK DATA DETERMINE HOW MANY LAYERS HAVE THE STRESS
!     OUTPUT REQUEST (SOUTI)
!     NOTE - FOR PCOMP1 OR PCOMP2 BULK DATA ENTRIES LAYER
!            STRESSES ARE OUTPUT FOR ALL LAYERS.
 
 nlayer = nlay
 
 IF (itype /= pcomp) GO TO 1230
 
 nstrqt = 0
 DO  k = 1,nlay
   IF (intz(pidloc+8+4*k) == 1) nstrqt = nstrqt + 1
 END DO
 nlayer = nstrqt
 
!     WRITE TOTAL NUMBER OF LAYERS WITH STRESS REQ TO OES1L
 
 1230 IF (lamopt == sym .OR. lamopt == symmem) nlayer = 2*nlayer
 
!     == 2.
 
 IF (kstrs == 1) CALL WRITE (oes1l,nlayer,1,0)
 
!     SET POINTER
 
 IF (itype == pcomp ) ipoint = pidloc + 8 + 4*nlay
 IF (itype == pcomp1) ipoint = pidloc + 8 +   nlay
 IF (itype == pcomp2) ipoint = pidloc + 8 + 2*nlay
 
!     FAILURE THEORY TO BE USED IN COMPUTING FAILURE INDICES
 
 fthr = intz(pidloc+5)
 
!     WRITE TO OUTPUT FILE TYPE OF FAILURE THEORY SELECTED
 
!     == 3.
 
 IF (kstrs == 1) CALL WRITE (oes1l,fthr,1,0)
 
!     SHEAR BONDING STRENGTH
 
 sb     = z(pidloc+4)
 findex = 0.0
 fbond  = 0.0
 fpmax  = 0.0
 fbmax  = 0.0
 fimax  = 0.0
 
!     SET TRNFLX IF INTERLAMINAR SHEAR STRESS CALCULATIONS
!     IS REQUIRED
 
 trnflx = .false.
 
!     TRANSVERSE SHEAR STRESS RESULTANTS QX AND QY
 
 v(1) = forsul(45)
 v(2) = forsul(46)
 trnflx = v(1) /= 0.0 .AND. v(2) /= 0.0
 IF (.NOT.trnflx) GO TO 1240
 IF (itype == pcomp) icontr = ipoint + 27*nlay
 IF (itype == pcomp1 .OR. itype == pcomp2) icontr = ipoint + 25 + 2*nlay
 
!     LAMINATE BENDING INERTIA
 
 ei(1)   = z(icontr+1)
 ei(2)   = z(icontr+2)
 
!     LOCATION OF NEUTRAL SURFACE
 
 zbar(1) = z(icontr+3)
 zbar(2) = z(icontr+4)
 
!     INTILIZISE
 
 1240 DO  ll = 1,2
   trnar(ll)  = 0.0
   trnshr(ll) = 0.0
 END DO
 
!     ALLOW FOR THE ORIENTATION OF THE MATERIAL AXIS FROM
!     THE USER DEFINED COORDINATE SYSTEM
 
 thetae = ACOS(phiout(69))
 thetae = thetae*degrad
 
!     SWITCH FOR THEMAL EFFECTS
 
 IF (ldtemp == -1) GO TO 1290
 
!     LAMINATE REFERENCE (OR LAMINATION) TEMPERATURE
 
 tsubo = z(ipoint+24)
 
!     MEAN ELEMENT TEMPERATURE
 
 tbar = tmean
 IF (tempp1 .OR. tempp2) tbar = stemp(1)
 IF (lamopt == mem .OR. lamopt == symmem) GO TO 1290
 IF (.NOT.(tempp1 .OR. tempp2)) GO TO 1290
 IF (.NOT.tempp1) GO TO 1260
 
!     TEMPERATURE GRADIENT TPRIME
 
 tprime = stemp(2)
 
 1260 IF (.NOT.tempp2) GO TO 1290
 
!     COMPUTE REFERENCE SURFACE STRAINS AND CURVATURES
!     DUE TO THERMAL MOMENTS
 
!     MOMENT OF INERTIA OF LAMINATE
 
 mintr = (tlam**3)/12.0
 
!     DETERMINE ABBD-MATRIX FROM PHIOUT(23-58)
 
 icount = 89 + 9*nnode
 DO  ll = 1,3
   DO  mm = 1,3
     nn = mm + 6*(ll-1)
     ii = mm + 3*(ll-1)
     abbd(ll  ,mm  ) = phiout(nn+22)*tlam
     abbd(ll  ,mm+3) = phiout(icount+ii)*(tlam*tlam)/(-6.0*const)
     abbd(ll+3,mm  ) = phiout(icount+ii)*(tlam*tlam)/(-6.0*const)
     abbd(ll+3,mm+3) = phiout(nn+43)*mintr
   END DO
 END DO
 
!     COMPUTE THERMAL REF STRAINS AND CURVATURES
!                                   -1
!        EZEROT-VECTOR =  ABBD-MATRIX   X  MTHR-VECTOR
 
 mther( 1) = 0.0
 mther( 2) = 0.0
 mther( 3) = 0.0
 mther( 4) = stemp(2)
 mther( 5) = stemp(3)
 mther( 6) = stemp(4)
 
 CALL invers (6,abbd,6,dumc,0,detrm,ising,indx)
 
 DO  ll = 1,6
   DO  mm = 1,6
     nn = mm + 6*(ll-1)
     stiff(nn) = abbd(ll,mm)
   END DO
 END DO
 
 CALL gmmats (stiff(1),6,6,0, mther(1),6,1,0, ezerot(1))
 
 1290 CONTINUE
 
 DO  ll = 1,6
   forsul(ll) = 0.0
 END DO
 
!     LOOP OVER NLAY
 
 DO  k = 1,nlay
   
!     ZSUBI -DISTANCE FROM REFERENCE SURFACE TO MID OF LAYER K
   
   zk1 = zk
   IF (k == 1) zk1 = zref
   IF (itype == pcomp ) zk = zk1 + z(pidloc+6+4*k)
   IF (itype == pcomp1) zk = zk1 + z(pidloc+7    )
   IF (itype == pcomp2) zk = zk1 + z(pidloc+7+2*k)
   
   zsubi = (zk+zk1)/2.0
   
!     LAYER THICKNESS
   
   ti = zk - zk1
   
!     CALCULATE STRAIN VECTOR AT STN ZSUBI
   
   DO  ir = 1,3
     epslne(ir) = (.5-zsubi/tlam)*epsb(ir) + (.5+zsubi/tlam)*epst(ir)
   END DO
   
!     LAYER ORIENTATION
   
   IF (itype == pcomp ) theta = z(pidloc+7+4*k)
   IF (itype == pcomp1) theta = z(pidloc+8+  k)
   IF (itype == pcomp2) theta = z(pidloc+8+2*k)
   
!     BUILD TRANS-MATRIX TO TRANSFORM LAYER STRAINS FROM MATERIAL
!     TO FIBRE DIRECTION.
   
   theta = theta*degrad
   
   c   = COS(theta)
   c2  = c*c
   s   = SIN(theta)
   s2  = s*s
   
   trans(1)  = c2
   trans(2)  = s2
   trans(3)  = c*s
   trans(4)  = s2
   trans(5)  = c2
   trans(6)  =-c*s
   trans(7)  =-2.0*c*s
   trans(8)  = 2.0*c*s
   trans(9)  = c2-s2
   
!     TRANSFORM STRAINS FROM ELEMENT TO FIBRE COORD SYSTEM
   
   CALL gmmats (trans(1),3,3,0, epslne(1),3,1,0, epsln(1))
   
!     SWITCH FOR TEMPERATURE EFFECTS
   
   IF (ldtemp == -1) GO TO 1470
   
!     CORRECT LAYER STRAIN VECTOR FOR THERMAL EFFECTS
   
!     LAYER THERMAL COEFFICIENTS OF EXPANSION ALPHA-VECTOR
   
   DO  ll = 1,3
     alpha(ll) = z(ipoint+13+ll)
   END DO
   
!     ELEMENT TEMPERATURE
   
   delt = tbar - tsubo
   
   IF (lamopt == mem .OR. lamopt == symmem) GO TO 1420
   IF (.NOT.tempp1) GO TO 1420
   
!     TEMPERATURE GRADIENT TPRIME
   
   delt = delt + zsubi*tprime
   
   1420 DO  ll = 1,3
     epslnt(ll) = -alpha(ll)*delt
   END DO
   
   IF (lamopt == mem .OR. lamopt == symmem) GO TO 1450
   IF (.NOT.tempp2) GO TO 1450
   
!     COMPUTE STRAIN DUE TO THERMAL MOMENTS
   
   DO  ll = 1,3
     epslnt(ll) = epslnt(ll) + (ezerot(ll) + zsubi*ezerot(ll+3))
   END DO
   
!     COMBINE MECHANICAL AND THERMAL STRAINS
   
   1450 DO  ll = 1,3
     epsln(ll) = epsln(ll) + epslnt(ll)
   END DO
   
   1470 CONTINUE
   
!     CALCULATE STRESS VECTOR STRESL IN FIBRE COORD SYS
   
!     STRESL-VECTOR  =  G-MATRIX  X  EPSLN-VECTOR
   
   CALL gmmats (z(ipoint+1),3,3,0, epsln,3,1,0, stresl(1))
   
!     USE FORCE RESTULANTS CALCULATED PREVIOUSLY
!     I.E. AT EXTREME FIBER STATIONS EXCEPT FOR THERMAL LOADING CASES
   
   IF (ldtemp == -1) GO TO 1490
   IF (kforce ==  0) GO TO 1490
   
!     TRANSFORM LAYER STRESSES TO ELEMENT AXIS
   
   IF (thetae > 0.0) theta = theta + thetae
   
!     BUILD STRESS TRANSFORMATION MATRIX
   
   c   = COS(theta)
   c2  = c*c
   s   = SIN(theta)
   s2  = s*s
   
   trans(1)  = c2
   trans(2)  = s2
   trans(3)  =-2.0*c*s
   trans(4)  = s2
   trans(5)  = c2
   trans(6)  = 2.0*c*s
   trans(7)  = c*s
   trans(8)  =-c*s
   trans(9)  = c2-s2
   
   CALL gmmats (trans(1),3,3,0, stresl(1),3,1,0, strese(1))
   
   DO  ir = 1,3
     forsul(ir) = forsul(ir) + strese(ir)*ti
     IF (lamopt == mem .OR. lamopt == symmem) CYCLE
     forsul(ir+3) = forsul(ir+3) - strese(ir)*ti*zsubi
   END DO
   
   1490 IF (fthr <= 0) GO TO 1530
   
!     WRITE ULTIMATE STRENGTH VALUES TO ULTSTN
   
   DO  ir = 1,6
     ultstn(ir) = z(ipoint+16+ir)
   END DO
   
!     CALL FTHR TO COMPUTE FAILURE INDEX FOR PLY
   
   IF (fthr == strain) GO TO 1510
   CALL failur (fthr,ultstn,stresl,findex)
   GO TO 1520
   
   1510 CALL failur (fthr,ultstn,epsln,findex)
   
!     DETERMINE THE MAX FAILURE INDEX
   
   1520 IF (ABS(findex) >= ABS(fpmax)) fpmax = findex
   
   1530 CONTINUE
   
!     SET POINTERS
   
   IF (itype == pcomp) icontr = ipoint + 25
   IF (itype == pcomp1 .OR. itype == pcomp2) icontr = ipoint + 23 + 2*k
   
   IF (lamopt == mem .OR. lamopt == symmem) GO TO 1570
   IF (.NOT.trnflx) GO TO 1570
   
!     CALCULATE INTERLAMINAR SHEAR STRESSES
   
   DO  ir = 1,2
     trnar(ir) = trnar(ir) + (z(icontr+ir))*ti*(zbar(ir)-zsubi)
   END DO
   
!     THE INTERLAMINAR SHEAR STRESSES AT STN ZSUBI
   
   DO  ir = 1,2
     trnshr(ir) = v(ir)*trnar(ir)/ei(ir)
   END DO
   
!     CALCULATE SHEAR BONDING FAILURE INDEX FB
!     NOTE- SB IS ALWAYS POSITIVE
   
   IF (sb == 0.0) GO TO 1570
   
   DO  ir = 1,2
     fb(ir) = ABS(trnshr(ir))/sb
   END DO
   
   fbond = fb(1)
   IF (fb(2) > fb(1)) fbond = fb(2)
   
!     CALCULATE MAX SHEAR BONDING INDEX
   
   IF (fbond >= fbmax) fbmax = fbond
   
   1570 CONTINUE
   
   IF (kstrs == 0) GO TO 1590
   
!     WRITE TO OUTPUT FILE THE FOLLOWING
!       4.    PLY (OR LAYER) ID
!     5,6,7.  LAYER STRESSES
!       8.    LAYER FAILURE INDEX
!       9.    IFLAG (= 1 IF FP.GE.0.999, DEFAULT = 0)
!     10,11.  INTERLAMINAR SHEAR STRESSES
!      12.    SHEAR BONDING FAILURE INDEX
!      13.    IFLAG (= 1 IF SB.GE.0.999, DEFAULT = 0)
   
!     CHECK LAYER STRESS OUTPUT REQUEST (SOUTI) FOR PCOMP BULK DATA
!     (NOT SUPPORTED FOR PCOMP1 OR PCOMP2 BULK DATA)
   
   IF (itype /= pcomp) GO TO 1580
   souti = intz(pidloc+8+4*k)
   IF (souti == 0) GO TO 1590
   1580 plyid = k
   
!     == 4.
   
   CALL WRITE (oes1l,plyid,1,0)
   
!     == 5,6,7.
   
   CALL WRITE (oes1l,stresl(1),3,0)
   
!     == 8.
   
   CALL WRITE (oes1l,findex,1,0)
   
!     SET IFLAG
   
   iflag = 0
   IF (ABS(findex) >= 0.999) iflag = 1
   
!     == 9.
   
   CALL WRITE (oes1l,iflag,1,0)
   
!     == 10,11.
   
   CALL WRITE (oes1l,trnshr(1),2,0)
   
!     == 12.
   
   CALL WRITE (oes1l,fbond,1,0)
   
!     SET IFLAG
   
   iflag = 0
   IF (ABS(fbond) >= 0.999) iflag = 1
   
!     == 13.
   
   CALL WRITE (oes1l,iflag,1,0)
   
   
!     UPDATE IPOINT FOR PCOMP BULK DATA ENTRY
   
   1590 IF (itype == pcomp .AND. k /= nlay) ipoint = ipoint + 27
   
 END DO
 
!     FALL HERE IF SYMMETRIC OPTION HAS BEEN EXERCISED
 
 IF (lamopt /= sym .AND. lamopt /= symmem) GO TO 2000
 
!     LOOP OVER SYMMETRIC LAYERS
 
 DO  kk = 1,nlay
   k = nlay + 1 - kk
   
!     ZSUBI -DISTANCE FROM REFERENCE SURFACE TO MID OF LAYER K
   
   zk1 = zk
   IF (itype == pcomp ) zk = zk1 + z(pidloc+6+4*k)
   IF (itype == pcomp1) zk = zk1 + z(pidloc+7    )
   IF (itype == pcomp2) zk = zk1 + z(pidloc+7+2*k)
   
   zsubi = (zk+zk1)/2.0
   
!     LAYER THICKNESS
   
   ti = zk - zk1
   
!     CALCULATE STRAIN VECTOR AT STN ZSUBI
   
   DO  ir = 1,3
     epslne(ir) = (.5-zsubi/tlam)*epsb(ir) + (.5+zsubi/tlam)*epst(ir)
   END DO
   
!     LAYER ORIENTATION
   
   IF (itype == pcomp ) theta = z(pidloc+7+4*k)
   IF (itype == pcomp1) theta = z(pidloc+8+  k)
   IF (itype == pcomp2) theta = z(pidloc+8+2*k)
   
!     BUILD TRANS-MATRIX TO TRANSFORM LAYER STRAINS FROM MATERIAL
!     TO FIBRE DIRECTION.
   
   theta = theta*degrad
   c   = COS(theta)
   c2  = c*c
   s   = SIN(theta)
   s2  = s*s
   
   trans(1)  = c2
   trans(2)  = s2
   trans(3)  = c*s
   trans(4)  = s2
   trans(5)  = c2
   trans(6)  =-c*s
   trans(7)  =-2.0*c*s
   trans(8)  = 2.0*c*s
   trans(9)  = c2 - s2
   
!     TRANSFORM STRAINS FROM MATERIAL TO FIBRE COORD SYSTEM
   
   CALL gmmats (trans(1),3,3,0, epslne(1),3,1,0, epsln(1))
   
!     SWITCH FOR TEMPERATURE EFFECTS
   
   IF (ldtemp == -1) GO TO 1770
   
!     CORRECT LAYER STRAIN VECTOR FOR THERMAL EFFECTS
   
!     LAYER THERMAL COEFFICIENTS OF EXPANSION ALPHA-VECTOR
   
   DO  ll = 1,3
     alpha(ll) = z(ipoint+13+ll)
   END DO
   
!     ELEMENT TEMPERATURE
   
   delt = tbar - tsubo
   IF (lamopt == symmem) GO TO 1720
   IF (.NOT.tempp1) GO TO 1720
   
!     TEMPERATURE GRADIENT TPRIME
   
   delt = delt + zsubi*tprime
   
   1720 DO  ll = 1,3
     epslnt(ll) = -alpha(ll)*delt
   END DO
   
   IF (lamopt == symmem) GO TO 1750
   IF (.NOT.tempp2) GO TO 1750
   
!     COMPUTE STRAIN DUE TO THERMAL MOMENTS
   
   DO  ll = 1,3
     epslnt(ll) = epslnt(ll) + (ezerot(ll) + zsubi*ezerot(ll+3))
   END DO
   
!     COMBINE MECHANICAL AND THERMAL STRAINS
   
   1750 DO  ll = 1,3
     epsln(ll)  = epsln(ll) + epslnt(ll)
   END DO
   
   1770 CONTINUE
   
!     CALCULATE STRESS VECTOR STRESL IN FIBRE COORD SYS
   
!     STRESL-VECTOR =  G-MATRIX  X  EPSLN-VECTOR
   
   CALL gmmats (z(ipoint+1),3,3,0, epsln,3,1,0, stresl(1))
   
!     COMPUTE FORCE RESULTANTS IF REQUESTED
   
   IF (ldtemp == -1) GO TO 1790
   IF (kforce ==  0) GO TO 1790
   
!     TRANSFORM LAYER STRESSES TO ELEMENT AXIS
   
   IF (thetae > 0.0) theta = theta + thetae
   
!     BUILD STRESS TRANSFORMATION MATRIX
   
   c   = COS(theta)
   c2  = c*c
   s   = SIN(theta)
   s2  = s*s
   
   trans(1)  = c2
   trans(2)  = s2
   trans(3)  =-2.0*c*s
   trans(4)  = s2
   trans(5)  = c2
   trans(6)  = 2.0*c*s
   trans(7)  = c*s
   trans(8)  =-c*s
   trans(9)  = c2 - s2
   
   CALL gmmats (trans(1),3,3,0, stresl(1),3,1,0, strese(1))
   
   DO  ir = 1,3
     forsul(ir) = forsul(ir) + strese(ir)*ti
     IF (lamopt == symmem) CYCLE
     forsul(ir+3) = forsul(ir+3) - strese(ir)*ti*zsubi
   END DO
   
   1790 IF (fthr <= 0) GO TO 1830
   
!     WRITE ULTIMATE STRENGTH VALUES TO ULTSTN
   
   DO  ir = 1,6
     ultstn(ir) = z(ipoint+16+ir)
   END DO
   
!     CALL FTHR TO COMPUTE FAILURE INDEX FOR PLY
   
   IF (fthr == strain) GO TO 1810
   CALL failur (fthr,ultstn,stresl,findex)
   GO TO 1820
   
   1810 CALL failur (fthr,ultstn,epsln,findex)
   
!     DETERMINE THE MAX FAILURE INDEX
   
   1820 IF (ABS(findex) >= ABS(fpmax)) fpmax = findex
   
   1830 CONTINUE
   
!     SET POINTERS
   
   IF (itype == pcomp) icontr = ipoint + 25
   IF (itype == pcomp1 .OR. itype == pcomp2) icontr = ipoint + 23 + 2*k
   
   IF (lamopt == symmem) GO TO 1870
   IF (.NOT.trnflx) GO TO 1870
   
!     CALCULATE INTERLAMINAR SHEAR STRESSES
   
   DO  ir = 1,2
     trnar(ir) = trnar(ir) + (z(icontr+ir))*ti*(zbar(ir)-zsubi)
   END DO
   
!     THE INTERLAMINAR SHEAR STRESSES AT STN ZSUBI
   
   DO  ir = 1,2
     trnshr(ir) = v(ir)*trnar(ir)/ei(ir)
   END DO
   
!     CALCULATE SHEAR BONDING FAILURE INDEX FB
!     NOTE- SB IS ALWAYS POSITIVE
   
   IF (sb == 0.0) GO TO 1870
   
   DO  ir = 1,2
     fb(ir) = ABS(trnshr(ir))/sb
   END DO
   
   fbond = fb(1)
   IF (fb(2) > fb(1)) fbond = fb(2)
   
!     CALCULATE MAX SHEAR BONDING INDEX
   
   IF (fbond >= fbmax) fbmax = fbond
   
   1870 CONTINUE
   
   IF (kstrs == 0) GO TO 1890
   
!     WRITE TO OUTPUT FILE THE FOLLOWING
!       4.     PLY (OR LAYER) ID
!     5,6,7.   LAYER STRESSES
!       8.     LAYER FAILURE INDEX
!       9.     IFLAG (= 1 IF FP.GE.0.999, DEFAULT = 0)
!     10,11.   INTERLAMINAR SHEAR STRESSES
!      12.     SHEAR BONDING FAILURE INDEX
!      13.     IFLAG (= 1 IF SB.GE.0.999, DEFAULT = 0)
   
!     CHECK LAYER STRESS OUTPUT REQUEST (SOUTI) FOR PCOMP BULK DATA
!     (NOT SUPPORTED FOR PCOMP1 OR PCOMP2 BULK DATA)
   
   IF (itype /= pcomp) GO TO 1880
   souti = intz(pidloc+8+4*k)
   IF (souti == 0) GO TO 1890
   1880 plyid = nlay + kk
   
!     == 4.
   
   CALL WRITE (oes1l,plyid,1,0)
   
!     == 5,6,7
   
   CALL WRITE (oes1l,stresl(1),3,0)
   
!     == 8.
   
   CALL WRITE (oes1l,findex,1,0)
   
!     SET IFLAG
   
   iflag = 0
   IF (ABS(findex) >= 0.999) iflag = 1
   
!     == 9.
   
   CALL WRITE (oes1l,iflag,1,0)
   
!     == 10,11.
   
   CALL WRITE (oes1l,trnshr(1),2,0)
   
!     == 12.
   
   CALL WRITE (oes1l,fbond,1,0)
   
!     SET IFLAG
   
   iflag = 0
   IF (ABS(fbond) >= 0.999) iflag = 1
   
!     == 13.
   
   CALL WRITE (oes1l,iflag,1,0)
   
!     UPDATE IPOINT FOR PCOMP BULK DATA ENTRY
   
   1890 IF (itype == pcomp) ipoint = ipoint - 27
 END DO
 
 2000 IF (fthr <= 0) GO TO 2010
 
!     DETERMINE 'FIMAX' THE MAX FAILURE INDEX FOR THE LAMINATE
 
 fimax = fpmax
 IF (fbmax > ABS(fpmax)) fimax = fbmax
 
!     == LAST-1.
 
 2010 IF (kstrs == 1) CALL WRITE (oes1l,fimax,1,0)
 
 iflag = 0
 IF (ABS(fimax) >= 0.999) iflag = 1
 
!     == LAST.
 
 IF (kstrs == 1) CALL WRITE (oes1l,iflag,1,0)
 
 IF (kforce ==  0) GO TO 2100
 IF (ldtemp == -1) GO TO 2100
 CALL WRITE (oef1l,elemid,1,0)
 CALL WRITE (oef1l,forsul(1),6,0)
 CALL WRITE (oef1l,forsul(45),2,0)
 
 2100 RETURN
 
!     ERROR MESSAGES
 
 2200 WRITE  (nout,2210) uwm
 2210 FORMAT (a25,' - NO PCOMP, PCOMP1 OR PCOMP2 DATA AVAILABLE FOR ',  &
     'LAYER STRESS RECOVERY BY SUBROUTINE SQUD42.')
 GO TO 2100
 2220 WRITE  (nout,2230) ufm
 2230 FORMAT (a23,', LAYER STRESS OR FORCE RECOVERY WAS REQUESTED WHILE'  &
     ,      ' PROBLEM WAS NOT SET UP FOR', /5X,'LAYER COMPUTATION')
 CALL mesage (-61,0,0)
END SUBROUTINE squd42
