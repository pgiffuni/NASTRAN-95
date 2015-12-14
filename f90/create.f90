SUBROUTINE create (gplst,x,u,deform,conmin,conmax,elmtid,store,  &
    lcor,b1,b2)
     
 
    INTEGER, INTENT(IN OUT)                  :: gplst(1)
    REAL, INTENT(IN)                         :: x(3,1)
    REAL, INTENT(IN)                         :: u(2,1)
    INTEGER, INTENT(IN OUT)                  :: deform
    REAL, INTENT(OUT)                        :: conmin
    REAL, INTENT(OUT)                        :: conmax
    INTEGER, INTENT(IN)                      :: elmtid(100)
    REAL, INTENT(OUT)                        :: store(202)
    INTEGER, INTENT(IN OUT)                  :: lcor
    INTEGER, INTENT(IN OUT)                  :: b1
    INTEGER, INTENT(IN OUT)                  :: b2
    INTEGER :: layer,idummy(2),est,prnt,scr1,  &
        oes1,stress, where,DIRECT,sub,estsym(7), skipwd(20),elid,esym,gpts(12),  &
        msg1(20),ERR(2),offset
 
    DIMENSION       isym(14),itype(14),pt(2,5),third(2),c(2,4), centrd(2)
    COMMON /BLANK / skip(12),est,skip1(7),prnt,ski(2),oes1,scr1,scr2, NEW
    COMMON /xxparm/ skip4(157),ncntr,cntr(50),icntvl,where,DIRECT,sub,  &
        flag,value,skp20(20),layer
    EQUIVALENCE     (NEW,newoes), (kq4,isym(13)), (kt3,isym(14))
    DATA   ntypes / 14  /
    DATA   isym   / 2HSH,2HT1,2HTB,2HTP,2HTM,2HQP,2HQM,2HT2,2HQ2,2HQ1,  &
        2HM1,2HM2,2HQ4,2HT3/,    kbar/2HBR/
    DATA   itype  / 4, 6, 7, 8, 9, 15, 16, 17, 18, 19, 62, 63, 64, 83/
    DATA   estsym / 7*0 /
    DATA   skipwd /-5,-6,-7,-1,-2,-3,-4,-5,-6,4*0,0,-1,-2,-5,-6,-7,0 /
    DATA   nmsg1  , msg1/20,4H(48H,4H no ,4HSTRE,4HSS c,4HALCU,4HLATI,  &
        4HON f , 4HOUND,4H for,4H ele,4HMENT,4H num,4HBER ,4H,i8,,  &
        4H19H  , 4H- el,4HEMEN,4HT ig,4HNORE,4HD.)                /
 
    twopi  = 8.0*ATAN(1.0)
    irdest = 0
    CALL gopen (scr1,gplst(b1),1)
    stress = oes1
    IF ((icntvl >= 4 .AND. icntvl <= 9) .AND. DIRECT == 2) stress = newoes
    IF (stress == oes1 .AND. (icntvl == 6 .OR. icntvl == 8 .OR.  &
        icntvl == 9)) GO TO 130
    IF (stress == oes1) skipwd(7) = -3
    CALL OPEN (*130,stress,gplst(b2),0)
    conmin = 0.0
    conmax = 0.0
    IEOR   = 0
 
    !     CREATE A LIST OF ELEMENTS TYPES TO BE PLOTTED IN THIS SET
 
    jtj  = 1
    k    = 1
    kest = 0
    CALL READ (*135,*135,est,esym,1,0,m)
    IF (icntvl == 20) GO TO 7
 
    !     ELIMINATE ALL BUT MAXSHEAR FOR CSHEAR ELEMENT
 
3   IF (esym == isym(1) .AND. icntvl /= 3) GO TO 7
 
    !     ELIMINATE MID STRESS FOR TRIA1, QUAD1, TRPLT, OR QDPLT
 
    IF ((esym == isym(2) .OR. esym == isym(10) .OR. esym == isym(4)  &
        .OR. esym == isym(6)) .AND. where == 3) GO TO 7
 
    !     ELIMINATE Z2 AND AVER STRESS FOR CTRMEM, CQDMEM, MEM1, MEM2
 
    IF ((esym == isym(5) .OR. esym == isym(7) .OR. esym == isym(11)  &
        .OR. esym == isym(12)) .AND. (where == -1 .OR. where == 3)) GO TO 7
 
    !     ELIMINATE Z1, Z2 AND MAX FOR TRIA2 OR TRBSC ELEMENTS
 
    IF ((esym == isym(8) .OR. esym == isym(3)) .AND.  &
        (IABS(where) == 1 .OR. where == 2)) GO TO 7
    DO  i = 1,ntypes
        IF (esym == isym(i)) GO TO 6
    END DO
    GO TO 7
6   estsym(k) = i
    k = k + 1
7   CALL fread (est,ngppe,1,0)
 
8   offset = 0
    IF (esym == kbar) offset = 6
    IF (esym == kt3 .OR. esym == kq4) offset = 1
 
    !     FLUSH TO NEXT SYMBOL
 
9   CALL fread (est,elid,1,0)
    IF (elid == 0) THEN
        SELECT CASE ( jtj )
            CASE (    1)
                GO TO 11
            CASE (    2)
                GO TO 25
        END SELECT
    END IF
    j = 1 + ngppe + offset
    CALL fread (est,0,-j,0)
    GO TO 9
 
    !     READ NEXT SYMBOL
 
11  CALL READ (*12,*12,est,esym,1,0,m)
    GO TO 3
 
    !     LOOP BACK UNTIL ALL EST SYMBOLS ARE IN CORE
 
12  k = k - 1
    CALL bckrec (est)
    jtj = 2
 
    !     NOTE THAT THE ASSUMPTION THAT STRESS AND EST FILES ARE ORDERED IN
    !     THE SAME WAY IS NO LONGER NECESSARY
 
15  IF (IEOR == 0) CALL fwdrec (*125,stress)
    IF (icntvl /= 20) GO TO 20
    CALL fwdrec (*125,stress)
    GO TO 25
20  CALL READ  (*125,*120,stress,idummy,2,0,m)
    CALL fread (stress,ieltyp,1,0)
    CALL fread (stress,isub,1,0)
    CALL fread (stress,detail,1,0)
    CALL fread (stress,eigen,1,0)
    eigen = SQRT(ABS(eigen))/twopi
    CALL fread (stress,0,-3,0)
    CALL fread (stress,nwds,1,1)
    IF (sub > 0 .AND. isub /= sub) GO TO 50
    IF (flag == 1.0 .AND. detail /= value) GO TO 50
    IF (flag == 2.0 .AND. ABS(eigen-value) > 1.0E-5) GO TO 50
    IEOR = 0
    DO  i = 1,k
        j = estsym(i)
        IF (j == 0) CYCLE
        IF (ieltyp == itype(j)) GO TO 25
    END DO
    GO TO 50
 
    !     YES, WE DO WANT THIS ELEMENT TYPES STRESS DATA.  FIND THIS TYPES
    !     ELEMENTS IN THE EST
 
25  CALL READ (*905,*905,est,esym,1,0,m)
    irdest = 1
    CALL fread (est,ngppe,1,0)
    IF (icntvl == 20 .OR. esym == isym(j)) GO TO 27
 
    !     FLUSH THE FILE UNTIL FOUND
 
    GO TO 8
 
27  CALL WRITE (scr1,esym,1,0)
    kest = kest + 1
    mem  = 0
    IF (ieltyp == 9  .OR. ieltyp == 16 .OR. ieltyp == 15 .OR.  &
        ieltyp == 62 .OR. ieltyp == 63) mem = 1
    !         TRMEM(9), QDMEM(16), QDPLT(15), QDMEM1(62), QDMEM2(63)
 
    iwds = skipwd(icntvl)
    IF (icntvl > 13) GO TO 29
    IF (mem == 1) iwds = iwds + 1
    IF (where == -1 .AND. mem /= 1) iwds = iwds - 8
    nwds = -nwds - iwds + 2
    IF (IABS(where) /= 1 .AND. mem /= 1) nwds = nwds + 8
    IF (where == -1 .AND. mem == 1) GO TO 50
    IF (ieltyp == 4 .AND. icntvl /= 3) GO TO 50
    !         SHEAR(4)
 
29  is = 0
 
    !     READ STRESS FILE
 
30  is = is + 1
    CALL READ (*58,*58,stress,elmtid(is),1,0,m)
    IF (icntvl <= 9 .OR. icntvl == 20) GO TO 35
    CALL fread (stress,nlayer,1,0)
    laytot = nlayer*11
    layskp = -((layer-1)*10+2)
    CALL fread (stress,0,layskp,0)
35  IF (ieltyp /= 4) GO TO 40
 
    !     MAXIMUM SHEAR FOR CSHEAR ELEMENT
 
    CALL fread (stress,store(is),1,0)
    CALL fread (stress,0,-2,0)
    IF (is == lcor) GO TO 60
    GO TO 30
40  CALL fread (stress,0,iwds,0)
    CALL fread (stress,store(is),1,0)
    IF (icntvl <= 9 .OR. icntvl == 20) GO TO 41
    nlfin = -(laytot-1+layskp+iwds)
    CALL fread (stress,0,nlfin,0)
    GO TO 30
41  IF (icntvl < 20) GO TO 42
    CALL fread (stress,0,-1,0)
    GO TO 30
42  IF (IABS(where) == 1 .OR. mem == 1) GO TO 45
    CALL fread (stress,0,-7,0)
    CALL fread (stress,contur,1,0)
45  CALL fread (stress,0,nwds,0)
    IF (mem == 1 .AND. is >= lcor) GO TO 60
    IF (mem   == 1) GO TO 30
    IF (where == 2) store(is) = AMAX1(store(is),contur)
    IF (where == 3) store(is) = (store(is)+contur)/2.0
    IF (is >= lcor) GO TO 60
    GO TO 30
 
    !     SKIP THIS TYPE
 
50  CALL fwdrec (*125,stress)
    GO TO 20
 
    !     END OF RECORD ON STRESS FILE
 
58  IEOR = 1
    is   = is - 1
 
    !     STORE STRESS VALUES WITH ELEMENT ID.S
 
60  CALL fread (est,elid,1,0)
    IF (elid == 0) GO TO 90
    CALL fread (est,0,-1,0)
    CALL fread (est,gpts,ngppe+offset,0)
 
    !     THE VERY NEXT LINE WAS ACCIDENTALLY DROPPED IN 88 VERSION
 
    IF (elid > elmtid(is)/10) GO TO 100
    DO  ist = 1,is
        IF (elid == elmtid(ist)/10) GO TO 70
    END DO
    ERR(1) = 1
    ERR(2) = elid
    CALL wrtprt (prnt,ERR,msg1,nmsg1)
    GO TO 60
 
    !     FIND ELEMENTS CENTROID
 
    70 DO  i = 1,ngppe
        ig = gpts(i)
        ig = IABS(gplst(ig))
        IF (deform /= 0) GO TO 74
        pt(1,i) = x(2,ig)
        pt(2,i) = x(3,ig)
        CYCLE
74      pt(1,i) = u(1,ig)
        pt(2,i) = u(2,ig)
    END DO
    third(1) = pt(1,3)
    third(2) = pt(2,3)
    INDEX = 1
    pt(1,ngppe+1) = pt(1,1)
    pt(2,ngppe+1) = pt(2,1)
    IF (ngppe < 4) GO TO 80
    INDEX = 4
    CALL centre (*90,pt(1,1),pt(2,1),pt(1,2),pt(2,2),pt(1,3),pt(2,3),  &
        pt(1,4),pt(2,4),centrd)
    third(1) = centrd(1)
    third(2) = centrd(2)
    80 DO  i = 1,INDEX
        CALL centre (*90,pt(1,i),pt(2,i),pt(1,i+1),pt(2,i+1),(third(1)+  &
            pt(1,i+1))*.5,(third(2)+pt(2,i+1))*.5,  &
            (third(1)+pt(1,i))*.5,(third(2)+pt(2,i))*.5,centrd)
        c(1,i) = centrd(1)
        c(2,i) = centrd(2)
    END DO
    IF (ngppe < 4) GO TO 90
    CALL centre (*90,c(1,1),c(2,1),c(1,2),c(2,2),c(1,3),c(2,3),c(1,4),  &
        c(2,4),centrd)
90  CALL WRITE (scr1,elid,1,0)
    IF (elid /= 0) GO TO 91
905 IF (kest == k) GO TO 120
    CALL bckrec (est)
    irdest = 0
    GO TO 15
91 CONTINUE
   CALL WRITE (scr1,store(ist),1,0)
   CALL WRITE (scr1,centrd,2,0)
   IF (conmin /= 0.0 .OR. conmax /= 0.0) GO TO 92
   conmin = store(ist)
   conmax = conmin
   GO TO 60
92 conmin = AMIN1(conmin,store(ist))
   conmax = AMAX1(conmax,store(ist))
   GO TO 60
 
   !     REFILL STRESS STORAGE AREA
 
100 is = 0
   IF (IEOR == 0) GO TO 30
   ERR(1) = 1
   ERR(2) = elid
   CALL wrtprt (prnt,ERR,msg1,nmsg1)
   GO TO 60
120 IF (stress == newoes) GO TO 126
   CALL READ  (*125,*125,oes1,0,-3,0,m)
   CALL fread (oes1,isub,1,0)
   CALL fread (oes1,detail,1,0)
   CALL fread (oes1,eigen,1,1)
   eigen = SQRT(ABS(eigen))/twopi
   IF (isub /= sub) GO TO 125
   IF (flag == 1.0 .AND. detail /= value) GO TO 125
   IF (flag == 2.0 .AND. ABS(eigen-value) > 1.0E-5) GO TO 125
   CALL fwdrec (*125,oes1)
   GO TO 120
125 CALL bckrec (stress)
126 CALL CLOSE (stress,2)
130 CALL WRITE (scr1,0,0,1)
   CALL CLOSE (scr1,1)
   IF (irdest > 0) THEN
       GO TO   135
   ELSE
       GO TO   140
   END IF
135 CALL bckrec (est)
140 CONTINUE

    RETURN
END SUBROUTINE create
