SUBROUTINE strscn (s_or_f)
     
!     THIS ROUTINE PERFORMS STRESS AND FORCE OUTPUT SCAN.
 
!     ACKNOWLEDGEMENT -
 
!     THIS ROUTINE WAS WRITTEN ORIGINALLY BY LOCKHEED/GEORGIA FOR USER
!     DMAP-ALTER APPLICATION. IT WAS GREATLY MODIFIED BY G.CHAN/SPERRY
!     SO THAT USERS CAN SPECIFY THE OUTPUT SCAN PARAMETERS FROM THE
!     CASE CONTROL SECTION VIA THE SCAN INPUT CARD(S).  ONLY A VERY
!     SMALL PORTION OF THE ORIGINAL PROGRAM REMAINS.  THE DMAP-ALTER
!     APPLICATION IS STILL AVAILABLE TO THE USER
 
!     THIS ROUTINE IS CALLED ONLY BY SCAN
!     IT DOES NOT OPEN NOR CLOSE ANY FILE
 
!     SCAN PARAMETER -
 
!          S OR F  - STRESS (1) OR FORCE (2) SCAN FLAG
!          INFILE  - INPUT FILE, EITHER STRESS OR FORCE OUTPUT FILE
!          OUFILE  - OUTPUT FILE FROM SCAN OPERATION, TO BE PRINTED
!                    AGAIN BY OFP
!          IOPT    - OPTION 1, SCAN BY AMAX-AMIN (.GT.AMAX AND .LT.AMIN)
!                    OPTION 2, SCAN BY NTOP-
!                  . TOP N LARGEST (TENSION) AND SMALLEST (COMPRESSION)
!                    IN STRESS SCAN.
!                  . TOP N LARGEST ONLY IF NO COMPRESSION STRESS PRESENT
!                  . TOP N SMALLEST ONLY IF NO TENSION STRESS PRESENT
!                  . TOP N VALUES SCAN FOR FORCES IF TOP N IS POSITIVE
!                  . LEAST N VALUES SCAN FOR FORCES OR MARGIN (STRESS)
!                    IF TOP N IS NEGATIVE
!                  - IOPT IS INITIALIZED IN SCAN
!                  - STRSCN WILL SET IOPT TO A NEGATIVE NUMBER IF INPUT
!                    FILE IS NOT A STRESS OR FORCE FILE
!          ISET    - A LIST OF ELEMENT IDS TO BE SCANNED
!          IEL     - ELEMENT TYPE (CODE) TO BE SCANNED
!          IELT    - ELEMENT NAME IN 2 BCD WORDS
!          ICASE   - USED LOCALLY FOR SUBCASE NUMBER.
!          SUBC    - CURRENT SUBCASE NO. USED IN SCAN AND STRSCN
!          OSUBC   - SUBCASE NO. PROCESSED PREVIOUSLY
!          ISORT   - SET LOCALLY TO 1 IF INPUT FILE DATA IS IN SORT1
!                    TYPE FORMAT, TO 2 IF IN SORT2
!          DEBUG   - LOCAL DEBUG FLAG, SET BY SCAN
!          OEL     - ELEMENT TYPE PROCESSED PREVIOUSLY
 
!     SEE SUBROUTINE SCAN FOR MORE PARAMETER DEFINITIONS
 
! *** IF SCAN IS CALLED BY USER VIA DMAP ALTER, WE HAVE
 
!          ISET      =-2
!          LBEG=LEND = 0, NOT USED
!          LCSE1 AND = BEGINNING AND ENDING POINTERS TO AN ELEM. LIST
!          LCSE2       (SORT1, ALL SUBCASES), OR A SUBCASE LIST (SORT2,
!                      ALL ELEMS) IF THEY ARE GIVEN. OTHERWISE, LCSE1=-1
!                      AND LCSE2=0
 
! *** IF SCAN IS CALLED BY RIGID FORMAT, WE HAVE
 
!          ISET      =-1  IMPLIES THAT ALL ELEMENTS ARE TE BE SCANNED,
!                         AND LBEG .GT. LEND
!          ISET      = N  IMPLIES THAT ELEM. ID SET N IS REQUESTED. THE
!                         SET DATA IS STORED IN IZ(LBEG) THRU IZ(LEND)
!          ISET      = 0  NOT DEFINED
!          LCSE1 AND =    ARE COMPONENT DUPLICATION FLAG (IDUPL) AND
!          LCSE2          INCREMENT FLAG (INC)
!                    = 0  IMPLIES NO DUPLICATION/INCR SET BY COMPONENT
!          LCSE1     =-2  SET AND USE LOCALLY IF SORT2 AND ELEM. SET ARE
!                         INVOLVED.
!          LBEG AND  =    ARE BEGINNING AND ENDING POINTERS TO THE ELEM.
!          LEND           ID SET, ALL ELEMS. (LBEG .GT. LEND IF ISET=-1)
 
 
 INTEGER, INTENT(IN OUT)                  :: s_or_f
 LOGICAL :: iopen,    jopen,    any,      debug
 LOGICAL :: layerd,   stress,   force
 CHARACTER (LEN=100) :: chead
 CHARACTER (LEN=12) :: field,    scnfld(6)
 INTEGER :: oufile,   head,     ihead(25),id(50),   iz(2),  &
            eor,      oel,      nam(2),   sortx(3), iscan(10),  &
            subc,     osubc,    quad4,    tria3,    topn,  iblank
 COMMON /output/ head(96)
 COMMON /xscanx/ infile,   oufile,   lcore,    lbeg,     lend,  &
                 iopen,    jopen,    iel,      iopt,     iset,  &
                 isort,    itrl3,    subc,     osubc,    oel,  &
                 debug,    lloop,    quad4,    tria3,    stress, force,    layerd
 COMMON /BLANK / ielt(2),  icomp,    topn,     amax,     amin,  &
                 lcse1,    lcse2,    icompx
 COMMON /system/ ibuf,     nout,     SPACE(6), nlpp,     dum, npage,    line
 COMMON /names / rd,       rdrew,    wrt,      wrtrew,   rew, norew,    eofnrw
 COMMON /zzzzzz/ z(2)
 EQUIVALENCE     (iz(1),z(1))
 EQUIVALENCE     (chead, ihead(1))
 DATA            nam,                sortx                       /  &
                 4HSTRS,   4HCN  ,   4HSORT,   4H1   ,   4H2     /
 DATA            t24,      eor,      noeor,    iblank            /  &
                 1.e+24,   1,        0    ,    1H                /
 
! *** SET ISCAN ARRAY FROM COMPONENT SPECIFICATION
 
 chead = ' '
 nscan=0
 ntop =IABS(topn)
!      PRINT *,' ENTERRING STRSCN,NTOP,ICOMP=',NTOP,ICOMP
 DO  i=1,30
   j=2**(i-1)
   IF (MOD(icomp,2*j) < j) GO TO 10
   nscan=nscan+1
   IF (i == 1) GO TO 670
   iscan(nscan)=i
   10   IF (icompx == 0) CYCLE
   IF (MOD(icompx,2*j) < j) CYCLE
   nscan=nscan+1
   iscan(nscan)=i+31
 END DO
!      PRINT *,' AFTER 20, ICOMP=',ICOMP
 j=2*j
 IF (icomp < j) GO TO 24
 nscan=nscan+1
 iscan(nscan)=31
 24   IF (icompx < j) GO TO 26
 nscan=nscan+1
 iscan(nscan)=62
 26   IF (icompx /= 0) CALL sort (0,0,1,1,iscan,nscan)
!      DEBUG = .TRUE.
!      PRINT *,' AFTER 26,NSCAN=',NSCAN
!      PRINT *,' AFTER 26,ISCAN=',(ISCAN(KB),KB=1,NSCAN)
 IF (.NOT.debug) GO TO 40
 WRITE (nout,30) iopen,jopen,ielt,iel,iset,icomp,icompx,lcse1,  &
     lcse2,isort,subc,itrl3,lbeg,lend,nscan
 30  FORMAT (//2X,'DEBUG/STRSCN',/,2(2X,L1),2X,2A4,13I8)
 IF (iopt == 2) WRITE (nout,33) ntop,(iscan(j),j=1,nscan)
 IF (iopt == 1) WRITE (nout,35) amax,amin,(iscan(j),j=1,nscan)
 33   FORMAT (5X,i9,31I3)
 35   FORMAT (5X,2E10.3,31I3)
 IF (lend > lbeg) WRITE (nout,38) iset,(iz(j),j=lbeg,lend)
 38   FORMAT (/5X,'SET',i8, (/5X,15I7))
 IF (nscan > 10) GO TO 590
 IF (nscan ==  0) GO TO 670
 
! *** INITIALIZATION
 
 40   idupl=lcse1
 inc  =lcse2
 jnc  =0
 IF (iset+1 < 0) THEN
   GO TO    70
 ELSE IF (iset+1 == 0) THEN
   GO TO    50
 ELSE
   GO TO    60
 END IF
 50   lcse1=0
 lcse2=0
 GO TO 70
 60   lcse1=iz(lbeg)
 lcse2=IABS(iz(lend))
 70   ns   =-1
 IF (lcse1 > lcse2) GO TO 640
 IF (.NOT.iopen .OR. .NOT.jopen) GO TO 600
 
 lbuf1=1
 lbuf3=0
 lbuf0=lbuf1-1
 il2  =0
 nrew =0
 icase=-1
 any  =.false.
 IF (osubc == 0) CALL fwdrec (*610,infile)
 IF (iset == -2 .OR. isort == 2 .OR. subc /= osubc) GO TO 100
 DO  i=1,3
   CALL bckrec (infile)
 END DO
 GO TO 90
 
! *** READ INPUT FILE ID RECORD AND SET ISORT FLAG FOR SORT1 OR SORT2
!     DATA TYPE
!     AT THIS TIME, ISORT MAY BE ALREADY SET BY PREVIOUS SCAN, OR ZERO
 
 85   IF (isort == 2 .OR. nrew >= 2) GO TO 490
 nrew=nrew+1
 CALL REWIND (infile)
 90   CALL fwdrec (*610,infile)
 100  CALL READ (*85,*110,infile,id,50,0,iwds)
 CALL READ (*610,*620,infile,head,96,1,iwds)
 isort=1
 IF (id(2) >= 2000) isort=2
 IF (id(2) >= 3000) GO TO 500
 
! *** SYNCHRONIZE SUBCASE ID (WHICH MAY NOT BE IN ASCENDING ORDER)
 
 IF (iset == -2 .OR. isort == 2 .OR. id(4) == subc) GO TO 130
 110  CALL REWIND (infile)
 nrew=nrew+1
 IF (nrew > 2) GO TO 490
 120  CALL fwdrec (*610,infile)
 CALL READ (*610,*110,infile,id,10,1,iwds)
 IF (id(4) /= subc) GO TO 120
 CALL bckrec (infile)
 GO TO 100
 
! *** SYNCHRONIZE ELEMENT TYPE (WHICH MAY NOT BE IN ASCENDING ORDER)
 
 130  IF (id(3)-iel == 0) THEN
   GO TO   140
 ELSE
   GO TO    90
 END IF
 140  oel =iel
 nrew=0
 
! *** POSITION DATA BLOCK FOR FIRST CASE AND BEGIN SCAN
 
 i=140
 IF (debug) WRITE (nout,145) i,iset,isort,icase,lcse1,lcse2,subc
 145  FORMAT (/9X,'DEBUG/STRSCN',I4,1H-,/2X,I9,11I7,3X,L1)        
 nwds =id(10)
 layerd = .false.

! FOR LAYERED QUAD4 AND TRIA3 IEL WILL BE EITHER 64 OR 83 RESPECTIVELY
! AND ID(10) WILL BE 10 (10 IS THE NUMBER OF WORDS PER LINE TO BE PRINTED).

 IF ( (iel == 64 .OR. iel == 83) .AND. id(10) == 10 ) layerd = .true.
 IF ( .NOT. layerd ) GO TO 144

! TO DETERMINE THE NUMBER OF WORDS PER EACH ELEMENT, WILL NEED TO DETERMINE
! HOW MANY LAYERS ARE PRESENT (Z(LBUF1+1)) AND ALLOW 10 WORDS PER LAYER
! PLUS A THREE WORD HEADER AND TWO EXTRA WORDS ON THE END.

 CALL READ (*610,*340,infile,iz(lbuf1),3,0,iwds)
 nwds = 3 + 10*iz(lbuf1+1) + 2

!      PRINT *,' COMPUTED NWDS=',NWDS
 CALL bckrec ( infile )

 144  nwds1=nwds+1
 lbuf2=lbuf1+nwds1
 lbuf3=lbuf2
 ih1  =lbuf2
 146  ih2  =lbuf2+nwds1*ntop-1
 il1  =ih2+1
 il2  =ih2+nwds1*ntop
 IF (il2 > lcore) GO TO 660
 ii   =0
 jj   =0
 kk   =0
 mm   =ih1
 nn   =il1
 idelm=-1
 icase=0
 any  =.false.
 IF (lcse1 == -2) lcse1=0
 IF (lcse1 > -2) ns=lbeg
 150  ns=ns-1
 ns=MIN0(ns,lbeg-1)
 IF (isort == 2 .AND. idelm /= -1) GO TO 90
 160  CALL READ (*610,*340,infile,iz(lbuf1),nwds,0,iwds)
 
! *** CHECK WHETHER THIS ELEMENT IS NEEDED FOR SCAN
!     WALK THROUGH SET ARRAY IF IT IS NECESSARY TO DO SO (R.F. ONLY)
!     CHECK SUBCASE NO. INSTEAD OF ELEM. ID IF THIS IS A USER DAMP ALTER
!     RUN WITH SORT2 TYPE DATA
 
 IF (iset  /= -2) GO TO 170
 IF (lcse1 <= -1) GO TO 200
 icase=iz(lbuf1)
 IF (isort == 1) icase=icase/10
 IF (icase >= lcse1) IF (icase-lcse2) 200,200,330
 GO TO 160
 170  IF (iset == -1 .OR. lcse1 == -2) GO TO 200
 idelm=iz(lbuf1)/10
 IF (isort == 2) idelm=id(5)/10
 180  ns=ns+1
 IF (ns > lend) GO TO 330
 izn=iz(ns)
 IF (izn >= 0) IF (idelm-izn) 150,190,180
 izn =IABS(izn)
 IF (idelm == izn) GO TO 190
 izn1=iz(ns-1)
 IF (izn1 <= 0 .OR. izn1 > izn) GO TO 640
 IF (idelm > izn ) GO TO 180
 IF (idelm < izn1) GO TO 640
 ns=ns-1
 190  IF (isort == 2) lcse1=-2
 
! *** MAKE SURE DEVICE CODE IS SET TO PRINT (SORT1 ONLY)
!     SET UP COMPONENT DUPLICATION/INC LOOP IF THEY ARE VALID
 
 200  IF (isort == 1) iz(lbuf1)=(iz(lbuf1)/10)*10 + 1
 i=200
 IF (debug) WRITE (nout,145) i,iz(lbuf1),icase,lcse1,lcse2,idupl,  &
     inc,ns,isort,nwds,iset,subc,iopt,any
 jdupl=1
 jnc  =0
 IF (iset == -2 .OR. idupl <= 0) GO TO 210
 jdupl=idupl
 jnc  =inc
 
! *** PICKUP MAX AND MIN OF CURRENT ELEMENT DATA
!     SAVE THESE MAX, MIN AS KEYS FOR SORTING LATER
 
 210  bmax=-t24
 bmin= t24
! QUAD4 (=64) AND TRIA3 (=83) WILL HAVE JDUPL NE 0 FOR LAMINATED
! CASE (I.E., WHEN LAYERD IS TRUE) FOR STRESS CASES
 IF ( ( iel == 64 .OR. iel == 83 ) .AND. jdupl == 49  &
     .AND. .NOT. layerd .AND. stress ) GO TO 492
 IF ( ( iel == 64 .OR. iel == 83 ) .AND. jdupl /= 49  &
     .AND.       layerd .AND. stress ) GO TO 492
!
! SET QUAD4 OR TRIA3 TO FALSE TO INDICATE TO SUBROUTINE SCAN THAT
! DATA FOR THESE ELEMENTS HAS BEEN FOUND IN EITHER OES1 OR OES1L FILES.
 IF ( iel == 64 ) quad4 = 1
 IF ( iel == 83 ) tria3 = 1
! IF JDUPL IS 49 THAN THIS IS A QUAD4 OR TRIA3 LAYERED ELEMENT, GET
! VALUE AFTER ELEMENT ID IN RECORD TO DETERMINE THE NUMBER OF LAYERS IN
! IN THE ELEMENT.
 IF ( jdupl == 49 ) jdupl = iz(lbuf0+2)
!
!      PRINT *,' BEFORE 230,JDUPL,JNC,NSCAN=',JDUPL,JNC,NSCAN
!      PRINT *,' BEFORE 230,ISCAN=',(ISCAN(KB),KB=1,NSCAN)
 DO  j=1,nscan
   i=iscan(j)
   IF (i > nwds) CYCLE
   kk=0
   DO  k=1,jdupl
     zk=z(lbuf0+i+kk)
     IF (zk > bmax) bmax=zk
     IF (zk < bmin) bmin=zk
     kk=kk+jnc
   END DO
 END DO
 
 IF (iopt == 2) GO TO 250
 
! *** OPTION ONE (IOPT=1, BY MAX-MIN)
!     ===============================
 
!     LBUF2 AND LBUF3 ARE BEGINNING AND ENDING POINTERS TO THE SCANNED
!     DATA ARRAY
 
 IF (bmax < amax .AND. bmin > amin) GO TO 160
 IF (lbuf3+nwds1 > lcore) GO TO 630
 any=.true.
 DO  i=1,nwds
   z(lbuf3+i)=z(lbuf0+i)
 END DO
 z(lbuf3)=bmax
 IF (bmin <= amin) z(lbuf3)=bmin
 lbuf3=lbuf3+nwds1
 GO TO 160
 
! *** OPTION TWO (IOPT=2)
!     ===================
 
!     TOP AND BOTTOM N VALUES FOR STRESSES
!     TOP VALUE SCAN FOR FORCES IF TOPN IS POSITIVE
!     BOTTEM VALUE SCAN FOR FORCES AND MARGIN ETC. IF TOPN IS NEGATIVE
 
!     II AND JJ ARE TOP AND BOTTOM ARRAY COUNTERS
!     MM IS POINTER TO THE SMALLEST OF THE TOP VALUSES
!     NN IS POINTER TO THE BIGGEST OF THE BOTTOM VALUSES
 
!     WHEN TOP AND BOTTOM ARRAYS ARE FILLED UP COMPLETELY WITH SCANNED
!     DATA (II=JJ=NTOP), IH1 AND IH2 ARE BEGINNING AND ENDING POINTERS
!     TO THE TOP VALUES, SIMILARY, IL1 AND IL2 ARE FOR THE BOTTOM VALUES
 
!     REMEMBER, SORF=1 FOR STRESS SCAN, SORF=2 FOR FORCE SCAN
!               NTOP=IABS(TOPN)
 
 250  any=.true.
 IF ( sorf == 2 .AND. topn <= 0) GO TO 290
 IF ((sorf == 1 .AND. bmax < 0.0) .OR.  &
     (ii >= ntop .AND. bmax < z(mm))) GO TO 290
 DO  i=1,nwds
   z(mm+i)=z(lbuf0+i)
 END DO
 z(mm)=bmax
 IF (ii >= ntop) GO TO 270
 ii=ii+1
 mm=mm+nwds1
 IF (ii < ntop) GO TO 290
 270  mm  =ih1
 bmax=+t24
 DO  i=ih1,ih2,nwds1
   IF (z(i) > bmax) CYCLE
   bmax=z(i)
   mm  =i
 END DO
 
 290  IF ( sorf == 2 .AND. topn >= 0) GO TO 160
 IF ((sorf == 1 .AND. bmin > 0 .AND. topn > 0) .OR.  &
     (jj >= ntop .AND. bmin > z(nn))) GO TO 160
 DO  i=1,nwds
   z(nn+i)=z(lbuf0+i)
 END DO
 z(nn)=bmin
 IF (jj >= ntop) GO  TO 310
 jj=jj+1
 nn=nn+nwds1
 IF (jj < ntop) GO TO 160
 310  nn  =il1
 bmin=-t24
 DO  i=il1,il2,nwds1
   IF (z(i) < bmin) CYCLE
   bmin=z(i)
   nn  =i
 END DO
 GO TO 160
 
! *** ELEM. ID LIST, OR SUBCASE LIST, HAS BEEN EXHAULSTED
!     (NOTE - SHOULD RETURN WITHOUT FWDREC HERE.  IF STRSCN IS CALLED
!             AGAIN, FWDREC WILL BE DONE AT 90)
 
 330  i=330
 IF (debug) WRITE (nout,145) i,idelm,icase,isort,ns,lbeg,lend, lcse1,lcse2
 IF (any) GO TO 350
 GO TO 510
 
! *** EOR READ (FROM 160)
 
 340  id(11)=0
 IF (any) GO TO 350
 id(11)=1
 id(10)=1
 nwds  =1
 il2   =ih1+1
 iz(il2)=01
 IF (isort == 2) iz(il2)=0
 iz(  2)=1
 
! *** THIS ELEMENT TYPE IS DONE.  BEGIN OUTPUT PROCEDURE
!     MAKE SURE DEVICE CODE IS SET TO PRINT, ALWAYS
!     ADD SCAN HEADER TO LABEL LINE
 
 350  id(1)=(id(1)/10)*10 + 1
 IF (isort == 2) id(5)=(id(5)/10)*10 + 1
 CALL WRITE (oufile,id(1),50,noeor)
 GO TO 530
 360  IF (iopt == 1) GO TO 370
 WRITE(chead(69:100), 8004 ) ntop
 8004  FORMAT('TOP AND BOTTOM  ',i4,' VALUES')
 GO TO 380
 370  WRITE(chead(69:100), 8005 ) amin, amax
 8005  FORMAT('EXCLUDING TO ',2(f8.1))
 380  CONTINUE
 head(73)=iblank
 DO  i = 1, 25
   head(i+64) = ihead(i)
 END DO
 head(95)=sortx(1)
 head(96)=sortx(2)
 IF (isort == 2) head(96)=sortx(3)
 CALL WRITE (oufile,head,96,eor)
 
 kk=1
 j =2
 IF (.NOT.any) GO TO 460
 
! *** (IOPT=2 ONLY) IF TOP AND BOTTOM ARRAYS ARE NOT FULL (I.E. II AND/
!     OR JJ ARE  .LT. NTOP), WE NEED TO SQUEEZE OUT SOME EMPTY CELLS IN
!     THE SPACE FROM Z(IH1) THRU Z(IL2) BEFORE SORTING THE SCANNED DATA
 
 IF (iopt == 2) GO TO 410
 il2=lbuf3
 GO TO 430
 410  IF (ii+jj == 2*ntop) GO TO 430
 kk =(ntop-ii)*nwds1
 il2=ih2+jj*nwds1
 DO  i=il1,il2
   z(i-kk)=z(i)
 END DO
 il2=ih1+(ii+jj)*nwds1-1
 
! *** MOVE MAX-MIN KEYS BEHIND IL2 SPACE AND BEGIN A 2-COLUMN SORT
!     THUS AVOID MASSIVE DATA TRANSFER DURING SORTING IF THE ORIGINAL
!     MULTI-COLUMNS SCANNED DATA WERE USED.
 
 430  kk=(il2-ih1+1)/nwds1
 IF (il2+2*kk > lcore) IF (iopt-1) 630,630,660
 i =ih1
 j =il2-2
 k =0
 440  j =j+2
 k =k+1
 z (j+1)=z(i)
 iz(j+2)=k
 i =i+nwds1
 IF (i < il2) GO TO 440
 k =2*kk
 CALL sortf (0,0,2,1,z(il2+1),k)
 
! *** BEGIN OUTPUT SCANNED DATA
 
 j =j+2
 IF (debug) WRITE (nout,450) j,kk,(z(il2+i),iz(il2+i+1),i=1,k,2)
 450  FORMAT (/9X,'DEBUG/STRSCN 450-',2I7,(/15X,E11.3,I5))        
 460  DO  k=1,kk
   i =ih1+(iz(j)-1)*nwds1
   CALL WRITE (oufile,iz(i+1),nwds,noeor)
   j =j-2
 END DO
 CALL WRITE (oufile,0,0,eor)
 itrl3=itrl3+2
 j    =kk*nwds
 IF (debug) WRITE (nout,750) j,itrl3,ii,jj
 IF (.NOT.any) GO TO 680
 
!*** EOF ON INPUT FILE (FROM 100)
!     NEXT ACTION WILL BE LOOP-BACK FOR MORE OR RETURN TO SCAN
 
!           R.F. (ISET.NE.-2)        I     USER DMAP ALTER (ISET=-2)
!     -------------------------------+----------------------------------
!     SORT1 - RETURN TO SCAN FOR     I  SORT1 - LOOP BACK FOR NEXT SUB-
!             NEXT SUBCASE           I          CASE, DISREGARDING THE
!                                    I          ELEM ID LIST
!     SORT2 - LOOP BACK FOR NEXT     I  SORT2 - LOOP BACK FOR NEXT ELEM,
!             ELEM. IF NO ELEM. LIST I          DISREGARDING THE SUBCASE
!           - IF ELEM. LIST EXISTS,  I          LIST
!             LOOP BACK ONLY IF MORE I
!             ELEM. TO BE PROCESSED  I
!             OTHERWISE, RETURN      I
 
 480  IF (il2 <= 0) GO TO 510
 il2 =-1
 nrew=0
 IF (iset  == -2) GO TO 100
 IF (isort ==  1) GO TO 510
 IF (lend > lbeg) IF (ns-lend) 100,510,510
 GO TO 100
 
! *** COULD NOT FIND ELEMENT OR SUBCASE
 
 490  IF (il2 /= 0) GO TO 510
 492  IF ( iel == 64 .AND. quad4 == 0 ) quad4 = -1
 IF ( iel == 83 .AND. tria3 == 0 ) tria3 = -1
 IF ( iel == 64 .OR. iel == 83 ) GO TO 510

 CALL fname (infile,z(1))
 WRITE (nout,710) ielt,z(1),z(2),nrew
 GO TO 510
 
! *** JOB DONE
 
 500  iopt=-id(2)
 510  DO  i=1,16
   head(i+73)=iblank
 END DO
 head(  95)=iblank
 head(  96)=iblank
 osubc=subc
 RETURN
 
! *** INTERNAL ROUTINE TO SYNTHESIZE THE COMPONENTS FOR HEADING
 
 530  IF (jnc <= 0) GO TO 550
 550   numfld = 0
!      PRINT *,' STRSCN,INC,IDUPL,NSCAN=',INC,IDUPL,NSCAN
!      PRINT *,' ISCAN=',(ISCAN(KB),KB=1,NSCAN)
 loop580:  DO  i = 1, nscan
!      PRINT *,' STRSCN CALLING STRNAM,I,ISCAN=',I,ISCAN(I)
   IF ( stress ) CALL strnam ( iel, iscan(i), field )
   IF ( force  ) CALL fornam ( iel, iscan(i), field )
   IF ( field  == ' ' ) CYCLE loop580
   IF ( numfld == 0   ) GO TO 570
   DO  k = 1, numfld
     IF ( field == scnfld( k ) ) CYCLE loop580
   END DO
   570   IF ( numfld >= 6 ) EXIT loop580
   numfld = numfld + 1
   scnfld( numfld ) = field
 END DO loop580
 585   IF ( numfld == 1 ) chead(1:19) = 'scanned BY FIELD:  '
 IF ( numfld /= 1 ) chead(1:19) = 'scanned BY FIELDS: '
 istr = 20
 DO  i = 1, numfld
   LEN = INDEX( scnfld(i), ' ' )
   iend = istr + LEN - 1
   IF ( iend > 51 ) GO TO 587
   IF ( i == 1 ) chead( istr:iend ) = scnfld(i)(1:LEN)
   IF ( i > 1 ) chead( istr:iend+2 ) = ', '//scnfld(numfld)(1:LEN)
   istr = iend
   IF ( i > 1 ) istr = iend + 2
   CYCLE
   587   chead( istr:51) = ',...'
   EXIT
 END DO
 589   CONTINUE
 IF ( iset <= 0 ) GO TO 360
 WRITE ( chead(52:68), 8001 ) iset
 8001  FORMAT(' SET:',i8 )
 GO TO 360
 
! *** FILE ERRORS
 
 590  WRITE (nout,720) ielt
 GO TO 670
 600  WRITE (nout,700) iopen,jopen
 GO TO 510
 610  IF (.NOT.any) GO TO 680
 j=-2
 GO TO 650
 620  j=-3
 GO TO 650
 630  WRITE (nout,730) ielt
 GO TO 510
 640  WRITE (nout,740) iset,lcse1,lcse2,lbeg,lend,ns,(iz(j),j=lbeg,lend)
 GO TO 510
 650  CALL mesage (j,infile,nam)
 GO TO 510
 660  j=(lcore-lbuf2+1)/(2*nwds1)
 WRITE (nout,760) ielt,ntop,j
 ntop=j
 GO TO 146
 670  WRITE (nout,770) icomp,icompx,ielt
 GO TO 510
 680  IF (debug) WRITE (nout,780) ielt,subc
 CALL mesage (30,220,ielt)
 GO TO 480
 
 700  FORMAT (//5X,'SYSTEM ERROR/STRSCN.  INPUT OR OUTPUT FILE NOT READY', 2(2X,L1))        
 710  FORMAT (//5X,'ELEMENT ',2A4,', OR SUBCASE, NOT IN DATA BLOCK ', &
                2A4,I7,' REWINDS')        
 720  FORMAT (//5X,'TOO MANY COMPONENTS SPECIFIED FOR ',2A4)        
 730  FORMAT (//5X,'INSUFFICIENT CORE TO PROCESS OUTPUT SCAN', &
               /5X,'SMALL VALUES OF AMAX-AMIN REQUIRE LARGE CORE REQUIREMENT') 
 740  FORMAT (//5X,'SYSTEM ERROR/STRSCN 740',7X,6I7, /,(5X,12I10))     
 750  FORMAT (/,I9,' WORDS WRITTEN TO OUTPUT FILE, RECORD',I5,9X,2I5)  
 760  FORMAT (//5X,'INSUFFICIENT CORE TO PROCESS OUTPUT SCAN FOR ',2A4, &
               /5X,'LARGE TOPN VALUE REQUIRES EXCESSIVE CORE REQUIREMENT.', &
                   'TOP2N IS AUTOMATICALLY REDUCED FROM',I5,' TO',I5)        
 770  FORMAT (//5X,'FIELD COMPONENT ERROR, CASE ABORT/STRSCN',5X,2I9,  &
                1X,2A4)        
 780  FORMAT (//5X,'NO APPLICABLE ELEMT OR SUBCASE/STRSCN',3X,2A4,I8)  
 
END SUBROUTINE strscn
