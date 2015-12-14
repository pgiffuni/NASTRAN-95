SUBROUTINE xrgdfm (newsol,oldsol,iapp,iufile,iopen,isize,iscr,  &
        nogo)
     
!     XRGDFM READS AND PROCESSES RIGID FORMATS
 
!     WRITTEN BY  RPK CORPORATION; DECEMBER, 1983
 
!     INPUT
!       ARGUMENTS
!         IAPP        =1, FOR DMAP APPROACH; =2, DISPLACEMENT APPRAOCH
!                     =3, HEAT APPROACH    ; =4, AERO APPROACH
!         IOPEN       ARRAY FROM OPEN CORE TO CONTAIN THE MODULE
!                     EXECUTION DECISION TABLE
!         ISIZE       NUMBER OF WORDS AVAILABLE IN THE IOPEN ARRAY
!         IUFILE      NAME OF USER'S FILE CONTAINING THE RIGID FORMAT
!         NEWSOL      ARRAY CONTAINING THE SOLUTION NUMBER FOLLOWED
!                     BY ALL SUBSET NUMBERS GIVEN BY THE USER
!         OLDSOL      SOLUTION ON PREVIOUS RUN IF THIS IS A RESTART
!       OTHER
!       /XRGDXX/
!         IRESTR      RESTART FLAG - NON-ZERO IF RUN IS A RESTART
!         NSUBST      NUMBER OF SUBSETS GIVEN BY THE USER
!         RECORD      ARRAY CONTAINING THE CARD IMAGE IN 20A4 FORMAT
!       /SYSTEM/
!         IDATE       ARRAY CONTAINING MONTH AND YEAR OF NASTRAN LEVEL
!         OPTAPE      UNIT USED FOR THE OUTPUT PRINT FILE
!       /TWO/
!         TWO         ARRAY CONTAINING THE VALUES OF THE POWERS OF 2.
!       /MEDMSK/
!         N1          NUMBER OF WORDS USED FOR THE CARD NAME RESTART
!                     TABLE
!         N2          NUMBER OF WORDS USED FOR THE FILE NAME RESTART
!                     TABLE
!         N3          NUMBER OF WORDS USED FOR THE RIGID FORMAT
!                     CHANGE RESTART TABLE
 
!     OUTPUT
!       ARGUMENTS
!         IOPEN       ARRAY CONTAINING THE MODULE EXECUTION DECISION
!                     TABLE
!       OTHER
!         /MEDMSK/
!           MEDMSK    MODULE EXECUTION DECISION MASK - SET IF SOLUTION
!                     CHANGE OCCURRED ON A RESTART
!         /SYSTEM/
!           ITHRML    SET TO NON-ZERO FOR A HEAT APPROACH
!         /PHAS11/
!           IPAS11    ARRAY FOR SUBSTRUCTURE CONTROLS-SET TO ZERO
!         /PHAS25/
!           IPAS25    SAME AS IPAS11
!         /PHAS28/
!           IPAS28    SAME AS IPAS11
!         /PHAS31/
!           IPAS31    SAME AS IPAS11
!         /PHAS37/
!           IPAS37    SAME AS IPAS11
!         /XRGDXX/
!           IDMAP     DMAP SEQUENCE NUMBER
!           IGNORE    FLAG SET TO IGNORE ANY CONTROL CARDS FOR THE
!                     CURRENT DMAP STATEMENT - IS SET WHEN THE DMAP
!                     STATEMENT IS TO BE DELETED BY THE SUBSET
!           IPHASE    PHASE NUMBER ASSOCIATED WITH THE ****PHS-
!                     CONTROL CARD
!           ITYPE     SET TO 'FILE' OR 'CARD' FOR TYPE OF CONTROL CARD
!           LIMIT     LOWER/UPPER LIMITS ASSOCIATED WITH THE VALUES
!                     OF A PARTICULAR CARD TYPE
!           MEMBER    NAME OF USER'S FILE CONTAINING A RIGID FORMAT
!                     THIS IS A 2-WORD ARRAY IN 2A4 FORMAT
!           NUMENT    NUMBER OF WORDS PER ENTRY IN THE MODULE EXECUTION
!                     DECISION TABLE
 
 
!     LOCAL VARIABLES
!       ASTERS        VARIABLE CONTAINING THE VALUE OF 4H****
!       CARD          VARIABLE CONTAINING THE VALUE OF 4HCARD
!       COMENT        VARIABLE CONTAINING THE VALUE OF 4H$$$$
!       DOLACR        VARIABLE CONTAINING THE VALUE OF 4H$*CA
!       DOLAFL        VARIABLE CONTAINING THE VALUE OF 4H$*FI
!       FILE          VARIABLE CONTAINING THE VALUE OF 4HFILE
!       FILTYP        ARRAY CONTAINING ACRONYMS FOR APPROACH
!       IBIT          BIT NUMBER TO SET IN THE MEDMSK
!       IFILL         VALUE TO BE USED TO INITIALIZE THE MODULE
!                     EXECUTION DECISION TABLE; =0, IF RESTART;
!                     =1, OTHERWISE
!       LU            FORTRAN LOGICAL UNIT NUMBER AS RETURN FROM RFOPEN
!                     =0, IF OPEN IS NOT SUCCESSFUL
!       INDEX         INDEX INTO CURRENT ENTRY OF MODULE EXEC.
!                     DECISION TABLE
!       ISOL          SOLUTION NUMBER
!       IWORD         WORD IN MEDMSK TO BE SET FOR RESTART FLAG
!       NEXT          FLAG INDICATING THAT A NEW DMAP STATEMENT IS
!                     TO BE PROCESSED; =0, IF NEW DMAP STATEMENT;
!                     =1, IF PROCESSING THE SAME DMAP STATEMENT
!       NUMSOL        ARRAY CONTAINING THE RESTART BITS ASSOCIATED
!                     WITH A RIGID FORMAT SWITCH DURING RESTART
!       MAXSOL        MAX. SOLUTION NUMBER
!       PHASE         ARRAY CONTAINING 'PHS1', PHS2', AND 'PHS3'
!       RFMT          VARIABLE CONTAINING THE VALUE 4HRFMT
!       SOLNUM        ARRAY CONTAINING THE ALPHA REPRESENTATIONS OF
!                     THE SOLUTION NUMBERS
 
!     FUNCTIONS
!       1. INITIALIZES SUBSTRUCTURE CONTROLS TO ZERO
!       2. CHECKS FOR USER SUPPLIED RIGID FORMAT
!       3. IF STANDARD RIGID FORMAT, VALIDATES SOLUTION NUMBER,
!          SETS MEDMSK IF A RESTART OCCURRED ON A DIFFERENT
!          RIGID FORMAT
!       4. SETS NUMENT=1 AND IFILL=1 IS NO RESTART - OTHERWISE
!          NUMENT=N1+N2+N3 AND IFILL=0
!       5. CALLS RFOPEN TO OPEN THE RIGID FORMAT
!       6. READS A CARD IMAGE FROM THE RIGID FORMAT FILE -
!          THE DATE AND YEAR OF THE RIGID FORMAT IS VALIDATED AGAINST
!          THAT THE LEVEL OF NASTRAN
!          RE-DEFINE NO. OF LINES PER OUTPUT PAGE IF 4TH WORD IS
!          PRESENT, .GT.20 .AND. .LE.99,  NO DATE CHECK IF THE ORD WROD
!          IS ****
!       7. READS A CARD FROM THE RIGID FORMAT FILE AND DOES THE
!          FOLLOWING DEPENDING ON THE TYPE OF CARD READ:
!          - FOR '$$$$' COMMENT CARDS, NEXT IS RESET
!          - FOR '****SBST' CARDS SUBROUTINE XRGSUB IS CALLED
!          - FOR '****CARD' CARDS SUBROUTINE XRGDCF IS CALLED
!          - FOR '****FILE' CARDS SUBROUTINE XRGDCF IS CALLED
!          - FOR '****RFMT' CARDS SUBROUTINE XRGDCF IS CALLED
!          - FOR '****PHS-' CARDS SUBROUTINE XRGSST IS CALLED
!          - OTHERWISE, THE CARD IS A DMAP AND WRITEN TO SCRATCH 315
!          (NOTE- FOR NON RESTARTS, THE ****CARD,****FILE,****RFMT
!          CARDS ARE BYPASSED.  FOR DMAP STATEMENTS THAT ARE
!          DELETED BY SUBSET CONTROLS, NO CONTROL CARDS ARE
!          PROCESSED EXCEPT FOR ****PHS- CARDS UNTIL THE NEXT
!          DMAP STATEMENT IS ENCOUNTERED)
!       8. WHEN A '$*CA' OR A '$*FI' CARD IS READ, PROCESSING OF
!          DMAP STATEMENTS TERMINATES - IF THE JOB IS NOT A RESTART
!          XRGDFM RETURNS.  OTHERWISE, A CHECK IS MADE TO ENSURE
!          THAT THE CARD NAME TABLE IS GIVEN FIRST FOLLOWED BY
!          THE FILE NAME TABLE.  SUBROUTINE XRGDTB IS CALLED TO
!          PROCESS BOTH TABLES.  AFTER THESE TABLES ARE PROCESSED,
!          XRGDFM RETURNS.
 
!     SUBROUTINES CALLED - RFOPEN,READ,WRITE,XRGSUB,XRGDCF,XRGSST,
!                          XRGDTB,MESAGE,RFCLOS
 
!     COMMENTS FROM G.C./UNISYS - ALL THE MACHINE DEPENDENT DSX* SUB-
!     ROUTINES ARE NO LONGER USED. SEE RFOPEN.  10/1990
 
!     CALLING SUBROUTINE - XCSA
 
!     ERROR MESSAGES 8023,504,8025,8026,8024,8037 MAY BE ISSUED
 
 
 INTEGER, INTENT(IN)                      :: newsol(12)
 INTEGER, INTENT(IN OUT)                  :: oldsol(12)
 INTEGER, INTENT(IN OUT)                  :: iapp
 INTEGER, INTENT(IN)                      :: iufile(2)
 INTEGER, INTENT(OUT)                     :: iopen(100)
 INTEGER, INTENT(IN OUT)                  :: isize
 INTEGER, INTENT(IN)                      :: iscr
 INTEGER, INTENT(OUT)                     :: nogo
 EXTERNAL        orf
 INTEGER :: record, BLANK, orf, two,  asters, sub(2),  &
            optape, card, FILE, rfmt, coment, subset, dolafl, dolacr, idate(3),  &
            filtyp(4),  solnum(20), numsol(50),oldnum,  phase(3), oldind
 INTEGER :: altfil
 DIMENSION       ioutbf(200)
 COMMON /altrxx/ altfil, newalt
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /xrgdxx/ irestr, nsubst, iphase, icol   , NUMBER, itype  ,  &
                 istate, ierror, num(2), ind    , nument         ,  &
                 record(20)    , ICHAR(80)      , limit(2)       ,  &
                 icount, idmap , iscrx , NAME(2), member(2)      , ignore
 COMMON /system/ ksystm(100)
 COMMON /two   / two(31)
 COMMON /xmdmsk/ n1,n2,n3,medmsk(7)
 COMMON /phas11/ ipas11(8)
 COMMON /phas25/ ipas25(14)
 COMMON /phas28/ ipas28(14)
 COMMON /phas31/ ipas31(2)
 COMMON /phas37/ ipas37(6)
 EQUIVALENCE     (ksystm( 2), optape), (ksystm(56), ithrml) ,  &
     (ksystm(42), idate(1)), (ksystm(69), isubal), (ksystm( 9), nlpp  )
 DATA    filtyp/ 4HDMAP, 4HDISP, 4HHEAT, 4HAERO    /
 DATA    solnum/ 1H1 ,1H2 ,1H3 ,1H4 ,1H5 ,1H6 ,1H7 ,  &
     1H8 ,1H9 ,2H10,2H11,2H12,2H13,2H14, 2H15,2H16,2H17,2H18,2H19,2H20     /
 DATA    card  / 4HCARD /, FILE   / 4HFILE /
 DATA    rfmt  / 4HRFMT /, BLANK  / 4H     /
 DATA    asters/ 4H**** /, coment / 4H$$$$ /
 DATA    subset/ 4HSBST /, dolacr / 4H$*ca /
 DATA    dolafl/ 4H$*fi /
 DATA    phase / 4HPHS1,   4HPHS2, 4HPHS3  /
 DATA    sub   / 4HXRGD,   4HFM            /
 DATA    nas   / 4HNAS  /, maxsol / 19     /
 
!     IN THE FOLLOWING TABLE, VALUES 187-209 ARE FOR STATICS,
!     210-213 ARE FOR HEAT, AND 214-217 ARE FOR AERO -
!     THIS PROVIDES FOR 31 DIFFERENT VALUES IN TOTAL (1 WORD)
 
 DATA    numsol/ 187, 188, 189, 190, 191, 192, 193, 194,  &
     195, 196, 197, 198, 199, 200, 201, 202,  &
     203, 204, 205,  -1,  -1,  -1,  -1, 210,  &
     -1, 211,  -1,  -1,  -1,  -1,  -1, 212,  &
     -1,  -1,  -1,  -1,  -1,  -1, 216, 214, 215, 9*-1 /
! WAS:
!     DATA    NUMSOL/
!    1                187, 188, 189, 190, 191, 192, 193, 194,
!    2                195, 196, 197, 198, 199, 200, 201, 202,
!    3                 -1,  -1,  -1,  -1, 207,  -1, 208,  -1,
!    4                 -1,  -1,  -1,  -1, 209,  -1,  -1,  -1,
!    5                 -1,  -1,  -1,  -1,  -1,  -1, 216, 214,
!    6                215, 9*-1 /
 
 iscrx = iscr
 idmap = 0
 DO  k = 1,8
   ipas11(k) = 0
 END DO
 DO  k = 1,14
   ipas25(k) = 0
   ipas28(k) = 0
 END DO
 DO  k = 1,2
   ipas31(k) = 0
 END DO
 DO  k = 1,6
   ipas37(k) = 0
 END DO
 IF (iufile(1) == 0) GO TO 100
 member(1) = iufile(1)
 member(2) = iufile(2)
 GO TO 210
 100  isol = newsol(1)
 SELECT CASE ( iapp )
   CASE (    1)
     GO TO 700
   CASE (    2)
     GO TO 120
   CASE (    3)
     GO TO 130
   CASE (    4)
     GO TO 140
 END SELECT
 120  IF (isol >= 1 .AND. isol <= maxsol) GO TO 200
 GO TO 700
 130  ithrml = 1
 isol = isol - 23
 IF (isol == 1 .OR. isol == 3 .OR. isol == 9) GO TO 200
 GO TO 700
 140  isol = isol - 30
 IF (isol == 9 .OR. isol == 10 .OR. isol == 11) GO TO 200
 GO TO 700
 200  member(1) = filtyp(iapp)
 member(2) = solnum(isol)
 210  CONTINUE
 
 oldind = oldsol(1)
 IF (oldind == 0 .OR. oldind == newsol(1)) GO TO 270
 
!     MAKE SURE CHECKPOINT TAPE FROM OLDER VERSION IS COMPATIBLE WITH
!     NEW CHANGE MADE IN 1991.
 
 IF (oldind /= 21 .AND. oldind /= 23 .AND. oldind /= 29) GO TO 220
 oldind = oldind + 3
 oldsol(1) = oldind
 IF (oldind == newsol(1)) GO TO 270
 
 220  oldnum = numsol(oldind)
 IF (oldnum <= 0) GO TO 270
 iword = ((oldnum-1)/31) + 1
 ibit  = oldnum - 31*(iword-1) + 1
 medmsk(iword) = orf(medmsk(iword),two(ibit))
 WRITE  (optape,240) oldsol(1),newsol(1),oldnum
 240  FORMAT (51H0*** switched solution for restart - OLD solution =,i4,  &
     16H, NEW solution =,i4,14H, bit NUMBER =,i4)
 270  IF (irestr /= 0) GO TO 280
 nument = 1
 ifill  = 1
 GO TO 290
 280  nument = n1 + n2 + n3
 ifill  = 0
 290  CONTINUE
 idmap = 0
 DO  kb = 1,nument
   iopen(kb) = ifill
 END DO
 INDEX = 1 - nument
 next  = 0
 CALL rfopen (member,lu)
 ignore = 0
 IF (lu == 0) GO TO 790
 READ (lu,305,ERR=720,END=730) record
 305  FORMAT (20A4)
 
!     BLANK OUT THE 19TH AND 20TH WORDS AS THEY
!     MAY CONTAIN SEQUENCE INFORMATION
 
 record(19) = BLANK
 record(20) = BLANK
 
!     ALLOW OPTIONS TO CHANGE NLPP LOCALLY, AND NOT TO CHECK RF DATE.
!     (THE NLPP OPTION HERE IS OBSOLETE. CAN BE EASILY DONE VIA NASINFO
!     FILE - 7/90)
 
 IF (record(3) == asters) GO TO 310
 IF (                           record(2) /= idate(3)) GO TO 770
 310  READ (lu,305,ERR=720,END=730) record
 
!     BLANK OUT THE 19TH AND 20TH WORDS AS THEY
!     MAY CONTAIN SEQUENCE INFORMATION
 
 record(19) = BLANK
 record(20) = BLANK
 IF (record(1) /= coment) GO TO 315
 IF (next      == 0     ) GO TO 310
 next = 0
 IF (INDEX <= isize) GO TO 310
 GO TO 740
 315  IF (record(1) == asters) GO TO 330
 IF (record(1) == dolacr .OR. record(1) == dolafl) GO TO 400
 IF (next == 1) GO TO 325
 IF (newalt == 0) GO TO 317
 CALL xrcard (ioutbf, 200, record)
 CALL WRITE (altfil, ioutbf(2), 2, 0)
 317 CONTINUE
 next   = 1
 idmap  = idmap + 1
 INDEX  = INDEX + nument
 DO  kb = 1,nument
   iopen(kb+INDEX-1) = ifill
 END DO
 325  CONTINUE
 CALL WRITE (iscr,record,18,0)
 ignore = 0
 GO TO 310
 330  IF (record(2) /= subset) GO TO 340
 IF (nsubst == 0) GO TO 310
 CALL xrgsub (iopen(INDEX),newsol(2))
 IF (ierror /= 0) nogo = 3
 GO TO 310
 340  IF (record(2) /= card) GO TO 350
 IF (irestr == 0 .OR. ignore == 1) GO TO 310
 limit(1) = 1
 limit(2) = n1*31
 CALL xrgdcf (iopen(INDEX))
 IF (ierror /= 0) nogo = 3
 GO TO 310
 350  IF (record(2) /= FILE) GO TO 360
 IF (irestr == 0 .OR. ignore == 1) GO TO 310
 limit(1) =  n1*31 + 1
 limit(2) = (n1+n2)*31
 CALL xrgdcf (iopen(INDEX))
 IF (ierror /= 0) nogo = 3
 GO TO 310
 360  IF (record(2) /= rfmt) GO TO 365
 IF (irestr == 0 .OR. ignore == 1) GO TO 310
 limit(1) = (n1+n2)*31 + 1
 limit(2) = (n1+n2+n3)*31
 CALL xrgdcf (iopen(INDEX))
 IF (ierror /= 0) nogo = 3
 GO TO 310
 365  DO  k = 1,3
   IF (record(2) /= phase(k)) CYCLE
   iphase = k
   CALL xrgsst (newsol)
   IF (ierror /= 0) nogo = 3
   GO TO 310
 END DO
 GO TO 750
 400  CALL WRITE (iscr,0,0,1)
 IF (newalt == 0) GO TO 500
 CALL WRITE (altfil, 0, 0, 1)
 CALL CLOSE (altfil, 1)
 500 CONTINUE
 CALL WRITE (iscr,iopen(1),INDEX+nument-1,1)
 IF (irestr == 0) GO TO 800
 itype = card
 IF (record(1) /= dolacr) GO TO 760
 limit(1) = 1
 limit(2) = n1*31
 CALL xrgdtb (lu)
 IF (ierror /= 0) nogo = 3
 itype = FILE
 IF (record(1) /= dolafl) GO TO 760
 limit(1) =  n1*31 + 1
 limit(2) = (n1+n2)*31
 CALL xrgdtb (lu)
 IF (ierror /= 0) nogo = 3
 GO TO 800
 
!     ERRORS
 
 700  WRITE  (optape,710) ufm,isol,filtyp(iapp)
 710  FORMAT (a23,' 8023, SOLUTION NUMBER',i4,' IS ILLEGAL FOR APPROACH'  &
     ,       a4)
 720  WRITE  (optape,725) ufm,member
 725  FORMAT (a23,' 8025, READ ERROR ON FILE ',2A4)
 GO TO 790
 730  WRITE  (optape,735) ufm,member
 735  FORMAT (a23,' 8025, UNEXPECTED EOF ENCOUNTERED ON FILE ',2A4)
 GO TO 790
 740  CALL mesage (-8,0,sub)
 GO TO 800
 750  WRITE  (optape,755) ufm,record
 755  FORMAT (a23,' 8026, THE FOLLOWING CARD HAS AN UNIDENTIFIED ',  &
     'FUNCTION AFTER ',6H'****', //20X,20A4)
 nogo = 3
 GO TO 310
 760  WRITE  (optape,765) ufm,itype,record
 765  FORMAT (a23,' 8024, EXPECTED A ',3H'$*,A4,1H',' CARD.',  &
     ' INSTEAD THE FOLLOWING CARD IS READ', //20X,20A4)
 GO TO 790
 770  WRITE  (optape,775) ufm,idate(1),idate(3),record(1),record(2)
 775  FORMAT (a23,' 8037, NASTRAN IS LEVEL ',2A4,  &
     ' BUT THE RIGID FORMAT DATA BASE IS LEVEL ',2A4)
 790  nogo = 3
 800  CALL rfclse (lu)
 RETURN
END SUBROUTINE xrgdfm
