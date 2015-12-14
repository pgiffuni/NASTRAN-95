BLOCK DATA semdbd
 
!     *****  PRINCIPAL BLOCK DATA PROGRAM FOR NASTRAN  *****
!     (NOTE - MACHINE DEPENDENT CONSTANTS ARE INITIALIZED IN BTSTRP)
 
!     REVISED 7/91 BY G.CHAN/UNISYS
!     MAKE SURE THERE IS NO VARIABLES OR ARRAYS NOT INITIALIZED. GAPS
!     OR MISSING INITIALIZED DATA MAY CAUSE PROBLEMS IN SOME MACHINES.
 
 IMPLICIT INTEGER (a-z)
 LOGICAL :: bitpas, first,  notyet, opnsof, pass,  pct,  star
 INTEGER :: ksystm(100)
 REAL :: oscar,  otapid, tapid , timdta, time,  tolel, x
 CHARACTER (LEN=31) :: sim
 CHARACTER (LEN=27) :: swm
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 EQUIVALENCE (ksystm(1),sysbuf)
  
!     -------------------     /GINOX  /     ----------------------------
 
!     GINOX WORDS USED IN GINO
!     VAX AND UNIX USE 636 WORDS (SEE GINO.MDS)
!DC   CDC USES ONLY 244 WORDS (SEE CDC IO6600). CDC IS CORE THIRSTY AND
!DC   THE 392 WORDS IN OTHERS HERE, AND ON THE DATA LINE BELOW, CAN BE
!DC   COMMENTED OUT TO SAVE CORE SPACE FOR THE CDC MACHINE.
 
 COMMON /ginox /    cdc(244) ,               others(392)
!            --------------------
!                     TOTAL= 636
 
!     -------------------     /XMSSG  /     ----------------------------
 
!     USER FATAL/WARNING/INFO AND SYSTEM FATAL/WARNING/INFO MESSAGES
 
 COMMON /xmssg / ufm, uwm, uim, sfm, swm, sim
 
!     -------------------     /NUMTPX /     ----------------------------
 
!     BCD-LOOK ALIKE FLOATING NUMBERS, USED ONLY BY NUMTYP SUBROUTINE.
!     DATA WILL BE LOADED FROM NASINFO.DOC FILE BY NSINFO
 
 COMMON /numtpx/ nbcd, bcd(19)
 
!     -------------------     / BLANK  /     ---------------------------
 
 COMMON /BLANK / iblnk(100)
 
!     -------------------     / NTIME  /     ---------------------------
 
!     THE NTIME COMMON BLOCK CONTAINS TIMING CONSTANT DATA FOR THE
!     CURRENTLY RUNNING MACHINE CONFIGURATION AS DETERMINED BY THE
!     TMTSIO AND TMTSLP SUBROUTINES
 
 COMMON /ntime / lntime, timdta(23)
 
!     -------------------     / XLINK  /     ---------------------------
 
!     SPECIFIES MODULE LINK RESIDENCE, AND PROVIDES LINK SWITCHING
!     INFORMATION FOR LINK DRIVER SUBROUTINES, XSEMi.
!     LXLINK = NUMBER OF WORDS IN MXLINK.
!     MAXLNK = MAX NUMBER OF LINKS - SEE XGPIBS IF THIS NUMBER IS
!              INCREASED.
!     MXLINK = MODULE LINK SPECIFICATION TABLE - THIS TABLE IS
!              INITIALIZED BY SUBROUTINE XGPIBS.
 
 COMMON /xlink / lxlink, maxlnk, mxlink(220)
 
!     -------------------     / SEM    /     ---------------------------
 
!     SEM DEFINES DATA FOR THE LINK DRIVERS (XSEMI).
!     MASK   = OSCAR MASK
!     MASK2,MASK3 =  OSCAR MASKS (MACHINE DEPENDENT).
!     NAME   = ARRAY OF LINK NAMES
 
 COMMON /sem   / mask, mask2, mask3, NAME(30)
 
!     -------------------     / SYSTEM /     ---------------------------
 
!     SYSTEM DEFINES VARIOUS MACHINE DEPENDENT, OPERATING SYSTEM AND
!     NASTRAN PARAMETERS.
!   1-
!     SYSBUF = (MACHINE DEPENDENT) NO. OF WORDS IN A GINO BUFFER.
!     OUTTAP = (MACHINE DEPENDENT) FORTRAN LOGICAL UNIT NO. FOR SYSTEM
!              PRINT OUTPUT
!     NOGO   = FLAG DEFINING EXECUTION STATUS DURING -FRONT END-.
!     INTP   = (MACHINE DEPENDENT) FORTRAN LOGICAL UNIT NO. FOR SYSTEM
!              INPUT
!     MPC    = MULTI-POINT CONSTRAINT SET ID  FOR CURRENT SUBCASE.
!     SPC    = SINGLE-POINT CONSTRAINT SET ID  FOR CURRENT SUBCASE.
!     LOGFL  = CONSOLE/LOGFILE MESSAGE CONTROL.
!     LOAD   = POINTER TO FIRST RECORD IN CASE CONTROL DATA BLOCK
!              FOR CURRENT SUBCASE.
!     NLPP   = (MACHINE DEPENDENT) NUMBER OF LINES PER PAGE OF PRINTED
!              OUTPUT.
!     MTEMP  = MATERIAL TEMPERATURE SET ID.
!  11-
!     NPAGES = CURRENT  PAGE COUNT.
!     NLINES = CURRENT NUMBER OF LINES ON CURRENT PAGE.
!     TLINES = TOTAL NUMBER OF LINES PRINTED IN JOB.
!     MXLINS = MAXIMUM NO. OF LINES OF PRINTED OUTPUT FOR THE PROBLEM.
!     DATE   = TODAY-S DATE, INTEGERS, 2 DIGITS EACH
!      (3)
!     TIMEZ  = CPU TIME IN SECONDS, WHEN PROBLEM BEGAN. NOT NECESSARY
!              ZERO. TIMEZ IS USED IN TMTOGO
!     ECHOF  = NUMBER INDICATING FORM OF BULK DATA ECHO.
!     PLOTF  = FLAG INDICATING REQUEST FOR STRUCTURAL PLOTS (NON-ZERO=
!              PLOT, SEE PLOTOPT IN SUBROUTINE NASCAR FOR MORE DETAILS)
!  21-
!     APPRCH = APPROACH FLAG (1 = FORCE, 2 = DISPL , 3 = DMAP).
!              APPRCH .LT. 0 MEANS THIS IS A RESTART.
!     LINKNO = CURRENT LINK NO. (IN BCD, E.G. NSXX)  INITIALLY SET TO
!              NS01 IN SUBROUTINE BTSTRP.  SUBSEQUENTLY SET TO THE
!              CORRECT LINK NO. IN SUBROUTINE ENDSYS.
!            = WAS THE MACHINE TYPE, MACH (LEVEL 17 AND EALIDER VERSION)
!     LSYSTM = LENGTH OF SYSTEM COMMON BLOCK.
!     ICFIAT = REPLACING EDTUMF FLAG HERE, WHICH IS NO LONGER USED.
!    (EDTUMF)  ICFIAT IS THE NUMBER OF WORDS PER FIAT ENTRY.
!            . IF ICFIAT=8, DATA BLOCK GINO TRAILER 6 WORDS ARE PACKED
!              INTO 4TH, 5TH, AND 6TH WORDS OF EACH FIAT ENTRY.
!            . IF ICFIAT=11, NO PACKING IN FIAT ENTRY, AND THE TRAILER 6
!              WORDS ARE SAVED IN 4TH THRU 6TH, AND 9TH THRU 11TH WORDS.
!              THE TRAILER WORDS ARE THEREFORE NOT BOUNDED BY SIZE
!              LIMITATION OF 65535 (HALF OF A 32-BIT WORD).
!            . THE FIAT POINTERS IN /XFIST/ MUST BE IN COMPLETE AGREE-
!              MENT WITH THE SELECTION OF ICFIAT=8, OR ICFIAT=11.
!              (SEE THE DATA SETTING OF /XFIST/ BELOW)
!            . THE REST OF NASTRAN .MIS ROUTINES ARE CODED TO HANDLE
!              ICFIAT=8 OR 11 AUTOMATICALLY. THE .MDS ROUTINES ARE NOT
!              AFFECTED SINCE THE 7TH AND 8TH WORDS OF THE FIAT ENTRY
!              REMAIN UNCHANGED.
!            = WAS EDTUMF FLAG, USED IN PRE-1987 NASTRAN VERSION
!     RFFLAG = RIGID FORMAT FLAG
!     CPPGCT = PAGE COUNT USED BY XCHK ROUTINE
!     MN     = NUMBER OF RINGS/NUMBER OF HARMONICS FOR AXISYMMETRIC
!              SHELL.
!     DUMMYI = (UNUSED WORD)
!     MAXFIL = MAXIMUM NUMBER OF UNITS TO BE ALLOCATED TO FIAT.
!     MAXOPN = MAXIMUM NUMBER OF FILES OPEN AT ONE TIME.
!  31-
!     HICORE = HI-CORE LENGTH FOR UNIVAC AND VAX
!     TIMEW  = PROBLEM START TIME (INTEGR SECONDS AFTER MIDNITE)
!     OFPFLG = OFP OPERATE FLAG - SET NON-ZERO WHEN OFP OPERATES
!     NBRCBU = (CDC ONLY) LENGTH OF FET + DUMMY INDEX
!              UNIVAC DRUM FILE ALLOCATION (1 FOR POSITION, 2 FOR TRACK)
!     LPRUS  = (CDC ONLY) NUMBER OF WORDS PER PHYSICAL RECORD UNIT (PRU)
!     NPRUS  = (CDC ONLY) NUMBER OF PRU-S PER GINO RECORD BLOCK
!     KSYS37 = ERROR CONTROL WORD, USED LOCALLY BY QPARMD AND QPARMR.
!              ALSO USED LOCALLY IN LINK1 FOR NASINFO FILE UNIT NO.
!     QQ     = HYDROELASTIC PROBLEM FLAG.
!     NBPC   = (MACHINE DEPENDENT) NO. OF BITS PER CHARACTER.
!     NBPW   = (MACHINE DEPENDENT) NO. OF BITS PER WORD.
!  41-
!     NCPW   = (MACHINE DEPENDENT) NO. OF CHARACTERS PER WORDS.
!     SYSDAT = THREE BCD WORD ARRAY CONTAINING MONTH, ' 19', AND LAST
!      (3)     TWO DIGITS OF YEAR OF SYSTEM GENERATION DATE.
!              THESE CELLS ARE SET BY SUBROUTINE NASCAR.
!     TAPFLG = WORD SET BY NASTRAN CARD TO INDICATE FILES TO BE TAPES
!              WHETHER OR NOT THEY ARE ON DISK.  BITS TURNED ON COUNTING
!              FROM RIGHT REPRESENT THE FILES IN XXFIAT.
!     ADUMEL = NINE WORD ARRAY CONTAINING DATA EXTRACTED FROM THE ADUM-I
!      (9)     CARDS BY IFP.
!  55-
!     IPREC  = PRECISION FLAG, 1=SP, 2=DP.
!     ITHRML = THERMAL ANALYSIS FLAG,  0=STRUCTURAL ANALYSIS,
!                                      1=THERMAL ANALYSIS.
!     MODCOM = NINE WORD ARRAY FOR MODULE COMMUNICATIONS.
!      (9)     SYSTEM(58), PRE-SELECT METHOD FOR MPYAD (1,2,3,DEFAULT=0)
!              SYSTEM(59), PLOT TAPE TRACK SIZE
!  66-
!     HDY    = THREE WORD ARRAY ALA SEW
!      (3)
!     SSCELL = MULTILEVEL SUBSTRUCTURE ANALYSIS COMMUNICATION CELL.
!     TOLEL  = SINGULARITY TOLERANCE FOR SMA1,EMG. RESET BY NASCAR.
!  71-
!     MESDAY = DAYFILE MESSAGE FLAG
!     BITPAS = CDC TAPE PROCESSING BIT - FALSE FOR LINK1 ONLY
!     PASS   = CDC MESSAGE AND TIMING FLAG - FALSE FOR LINK1 ONLY
!     ITIME  = WAS: WALL TIME ELAPSED SINCE PROBLEM START (SECONDS)
!            = IS : PROBLEM START TIME IN SECONDS SINCE JAN-1-1970,
!                   GREENWICH-MEAN-TIME (GMT)
!     CTIME  = WAS: CENTRAL PROCESSOR TIME SINCE PROBLEM START (SECONDS)
!            = IS : PRINT FLAG FOR DMAP SEQUENCE NO. AND NAME, ALL LINKS
!                   (SEE NASTRN OR NAST01.MDS)
!     NOSBE  = (CDC ONLY) FLAG FOR NOS(0) OR NOSBE(1)
!     BANDIT = BANDIT OPTION FLAG (SEE BANDIT FOR MORE DETAILS)
!     PZEL   = PIEZOELECTRIC PROBLEM FLAG (INPUT VIA NASTRAN SYSTEM(78))
!     SWITCH = SENSE SWITCH BITS FOR DIAG CARD AND USED BY SSWTCH
!      (3)
!  82-
!     ICPFLG = CHECKPOINT FLAG (0 = NO CHECKPOINT, 1 = CHECKPOINT)
!     JRUN   = JRUN FOR VARIAN  (HEAT PROBLEM)
!     JMAX   = JMAX FOR VARIAN  (HEAT PROBLEM)
!     LINTC  = MAX. ALLOWABLE LINES OF INTERSECTION USED IN HDPLOT
!     INTRA  = INTERACTIVE REQUEST FLAG FOR PLOT, OUTPUT, AND SCAN
!              (0=NONE, 1=PLOT ONLY, 2=OUTPUT PRINT AND SCAN ONLY,
!               3=BOTH)
!     OSPCNT = BAR OFFSET WARNING MESSAGE IF OFFSET BAR LENGTH EXCEEDS
!              NON-OFFEST LENGTH BY THIS LIMIT (DEFAULT IS 15 PERCENT)
!     K88 90 = 3 WORDS RESERVED FOR USER. WILL NOT BE USED BY COSMIC
!              =========================
!  91-
!     LPCH   = (MACHINE DEPENDENT) FORTRAN LOGICAL UNIT NO. FOR PUNCH
!     LDICT  = FORTRAN LOGICAL UNIT NO. FOR RESTART DICTIONARY PUNCH
!     IAEROT = INTEGER FLAG INDICATING AERODYNAMIC THEORY
!              (SPECIFIED VIA NASTRAN CARD AND USED ONLY IN APDB MODULE
!              0 FOR COMPRESSOR BLADES, THEORY 6, DEFAULT,
!              1 FOR SWEPT TURBOPROP. BLADES, THEORY 7)
!     KSYS94 = FLAG FOR REMOVALS OF MPYDRI(1), MPY4T(10), NEW FBS(100),
!              TRNSPS(1000), AND NEW FBS IN FEER(10000)
!     SPERLK = NASTRAN SUPERLINK FLAG. SET BY SEMDBD OR NASTRN
!              FOR UNIX MACHINE
!     LEFT   = (85 UNUSED WORDS).  KSYS99 USED IN ERRTRC
 
 
 DIMENSION       date(3),sysdat(3),adumel(9),modcom(9),hdy(3),  &
                 switch(3),k88_90(3) ,left(56),left2(28)
 COMMON /system/ sysbuf,outtap,nogo  ,intp  ,mpc   ,spc   ,logfl ,  &
                 load  ,nlpp  ,mtemp ,npages,nlines,tlines,mxlins,  &
                 date  ,timez ,echof ,plotf ,apprch,linkno,lsystm,  &
                 icfiat,rfflag,cppgct,mn    ,dummyi,maxfil,maxopn,  &
                 hicore,timew ,ofpflg,nbrcbu,lprus ,nprus ,ksys37,  &
                 qq    ,nbpc  ,nbpw  ,ncpw  ,sysdat,tapflg,adumel,  &
                 iprec ,ithrml,modcom,hdy   ,sscell,tolel ,mesday,  &
                 bitpas,pass  ,itime ,ctime ,nosbe ,bandit,pzel  ,  &
                 switch,icpflg,jrun  ,jmax  ,lintc ,intra ,ospcnt,  &
                 k88_90,lpch  ,ldict ,iaerot,ksys94,sperlk,left,loglin,left2
 
!     -------------------     / XFIST  /     ---------------------------
 
!     XFIST IS THE FILE STATUS TABLE (FIST).
!     NFIST  = TOTAL NO. OF ENTRIES IN FIST.
!     LFIST  = NO. OF ENTRIES IN THE CURRENT FIST.
!     FIST   = TABLE OF TWO-WORD ENTRIES.
!              FIRST WORD IS GINO FILE NAME.
!              SECOND WORD POINTS TO XFIAT IF .GT. 0 (I.E. NON-PERMANENT
!              ENTRY), OR POINTS TO XXFIAT IF .LE. 0 (I.E. PERMANENT
!              ENTRY). SIGN BIT MUST BE SET FOR ZERO POINTER ON 7094.
 
 COMMON /xfist / nfist, lfist, fist(112)
 
!     -------------------     / XPFIST /     ---------------------------
 
!     XPFIST DEFINES THE NO. OF PERMANENT ENTRIES IN THE FIST.
 
 COMMON /xpfist/ npfist
 
!     -------------------     / XXFIAT /     ---------------------------
 
!     XXFIAT IS EXECUTIVE FILE ALLOCATION TABLE.
 
 COMMON /xxfiat/ xxfiat(24)
 
!     -------------------     / XFIAT  /     ---------------------------
 
!     XFIAT IS THE MODULE FILE ALLOCATION TABLE (FIAT).
!     MFIAT  = NO. OF UNIQUE FILES IN FIAT.
!     NFIAT  = TOTAL NO. OF ENTRIES IN FIAT.
!     LFIAT  = NO. OF ENTRIES IN CURRENT FIAT.
!     FIAT   = TABLE OF 8 OR 11 WORDS PER ENTRY OF GINO FILES
!              (DEFAULT IS SET BY ICFIAT, THE 24TH WORD OF /SYSTEM/)
!            . 1ST WORD DEFINES THE FILE + PURGE,EQVIV,SETUP,ETC INFO.
!            . 2ND AND 3RD WORDS DEFINE THE DATA BLOCK NAME (IN BCD)
!              WHICH IS ATTACHED TO THE FILE.
!            . SEE ICFIAT (24TH WORD OF /SYSTEM/) FOR THE DESCRIPTION
!              OF 4TH THRU 8TH (OR 11TH) WORDS.
!            . SET FIAT(880) IF 11-WORD/ENTRY TABLE IS USED, AND
!              SET FIAT(640) IF  8-WORD/ENTRY TABLE IS USED
 
 COMMON /xfiat / mfiat, nfiat, lfiat, fiat(1100)
 
!     -------------------     / OSCENT /     ---------------------------
 
!     OSCENT DEFINES A STORAGE BLOCK FOR THE CURRENT OSCAR ENTRY.
 
 COMMON /oscent/ oscar(200)
 
!     -------------------     / OUTPUT /     ---------------------------
 
!     OUTPUT DEFINES A STORAGE BLOCK WHERE PROBLEM TITLE, SUBTITLE
!     AND LABEL ARE STORED.
 
 COMMON /output/ output(224)
 
!     -------------------     / XDPL   /     ---------------------------
 
!     XDPL DEFINES THE DATA POOL DICTIONARY.
!     MDPL   = POINTER TO NEXT AVAILABLE FILE.
!     NDPL   = TOTAL NO. OF ENTRIES IN DPL.
!     LDPL   = CURRENT NO. OF ENTRIES IN DPL.
!     DPL    = TABLE OF THREE-WORD ENTRIES
!              1ST + 2ND WORDS ARE DATA BLOCK NAME
!              3RD WORD DEFINES EQUIV FLAG, APPROX SIZE OF DATA BLOCK
!              AND FILE NO. IN THE POOL.
 
 COMMON /xdpl  / mdpl, ndpl, ldpl, dpl(240)
 
!     -------------------     / XVPS   /     ---------------------------
 
!     XVPS IS THE VARIABLE PARAMETER STORAGE TABLE.
!     VPS(1) = TOTAL NO. OF WORDS IN VPS.
!     VPS(2) = POINTER TO LAST WORD USED IN VPS.
!     VPS(3) = TABLE FOR STORAGE OF PARAMETERS
!              (VARIABLE NO OF WORDS/ENTRY).
 
 COMMON /xvps  / vps(600)
 
!     -------------------     / STAPID /     ---------------------------
 
!     STAPID CONTAINS THE I.D. FOR THE NEW AND OLD PROBLEM TAPES.
!     TAPID  = SIX-WORD I.D. FOR NEW PROBLEM TAPE.
!     OTAPID = SIX-WORD I.D. FOR OLD PROBLEM TAPE.
!     IDUMF  = (OBSOLETE) ID FOR USER-S MASTER FILE.
 
 COMMON /stapid/ tapid(6), otapid(6), idumf
 
!     -------------------     / STIME  /     ---------------------------
 
!     STIME DEFINES USER-S ESTIMATED PROBLEM SOLUTION TIME.
 
 COMMON /stime / time(2)
 
!     -------------------     / XCEITB /     ---------------------------
 
!     XCEITB DEFINES LOOP CONTROL PARAMETERS FOR THE CONTROL ENTRY
!     INTERP.
!     CEI(1) = TOTAL NO. OF WORDS IN TABLE.
!     CEI(2) = POINTER TO LAST WORD USED.
!     CEI(3) = TABLE OF FOUR-WORD ENTRIES.
 
 COMMON /xceitb/ cei(42)
 
!     -------------------     / XMDMSK /     ---------------------------
 
!     XMDMSK DEFINES MASK FOR MODIFIED RESTART.
!     NMSKCD = NUMBER OF MASK WORDS FOR CARDS (CURRENTLY SET TO 3)
!     NMSKFL = NUMBER OF MASK WORDS FOR FILES (CURRENTLY SET TO 3)
!     NMSKRF = NUMBER OF MASK WORDS FOR RIGID FORMATS (CURRENTLY 1)
!     MSK    = MASK OF (NMSKCD+NMSKFL+NMSKRF) WORDS (31 BITS/WORD)
 
 COMMON /xmdmsk/ nmskcd, nmskfl, nmskrf, msk(7)
 
!     -------------------     / MSGX   /     ---------------------------
 
!     MSGX DEFINES A TABLE WHERE MESSAGES ARE QUEUED.
!     NMSG   = NUMBER OF MESSAGES CURRENTLY IN QUEUE.
!     MSGLG  = TOTAL NO. OF ENTRIES IN THE MESSAGE QUEUE.
!     MSG    = TABLE OF FOUR-WORD ENTRIES WHERE MESSAGES ARE STORED.
 
 COMMON /msgx  / nmsg, msglg, msg(4,40)
 
!     -------------------     / DESCRP /     ---------------------------
 
!     COMMENTS FROM G.CHAN/UNISYS, 7/1991
!     LABEL COMMON /DESCRP/ APPEARS IN THE FOLLOWING SUBROUTINES. BUT
!     IT IS ACTUALLY NEVER USED.
!        SEMDBD, DECOMP, GENVEC, GFBS,   TRANSP, CDCOMP, CTRNSP, INVTR,
!        MTIMSU, MTMSU1, CDIFBS, INTFBS, MATVC2, MATVEC, ORTCK,  CINFBS,
!        CMTIMU, AND INVFBS
!     IN ADDITION, INTPK IN VAX, IBM, CDC AND UNIVAC, DOES NOT USE THIS
!     /DESCRP/ LABEL COMMON.
!     STARTING IN 1992 VERSION, LABEL COMMON /DESCRP/ IS COMPLETELY
!     REMOVED FROM ALL NASTRAN SUBROUTINES.
 
!     DESCRP IS A STORAGE BLOCK USED BY SUBROUTINE INTPK.
!     LENGTH = TOTAL NO. OF WORDS IN BLOCK.
 
!     COMMON /DESCRP/ LENGTH, BLOCK(50)
 
!     -------------------     / TWO    /     ---------------------------
 
!     TWO DEFINES THE BITS IN A 32-BIT COMPUTER WORD (FROM LEFT TO RT).
!     MZERO = WILL BE SET TO -0.0 (= LSHIFT(1,NBPW-1) = SIGN BIT ON) BY
!             BTSTRP, AND WILL BE USED BY NUMTYP
 
 COMMON /two   / two(32),mzero
 
!     -------------------     / NAMES  /     ---------------------------
 
!     NAMES DEFINES VALUES FOR GINO FILE OPTIONS,ARITHMETIC TYPES
!     AND MATRIX FORMS.
 
 COMMON /names /rd    ,rdrew ,wrt   ,wrtrew,rew   ,norew ,eofnrw,  &
                rsp   ,rdp   ,csp   ,cdp   ,square,rect  ,diag  ,  &
                lower ,upper ,sym   ,row   ,ident
 
!     -------------------     / TYPE   /     ---------------------------
 
!     TYPE DEFINES PROPERTIES AS A FUNCTION OF ARITHMETIC TYPE.
!     PRC    = PRECISION (1=SP, 2=DP).
!     NWDS   = NO. OF WORDS PER ELEMENT.
!     RC     = ARITHMETIC (1=REAL, 2=COMPLEX).
!     X      = PAD TO DEFINE WORK AREA.
 
 COMMON /TYPE  / prc(2), nwds(4), rc(4), x(6)
 
!     -------------------     / BITPOS /     ---------------------------
 
!     BITPOS DEFINES THE BIT POSITIONS FOR THE DEGREES-OF-FREEDOM IN
!     USET, AND HOLLERITH CHARACTERS DESCRIBING DEGREES-OF-FREEDOM.
 
 COMMON /bitpos/ um  ,uo  ,ur  ,usg ,usb ,ul  ,ua  ,uf  ,us  ,un  ,  &
                 ug  ,ue  ,up  ,une ,ufe ,ud  ,ups ,usa ,uk  ,upa ,  &
                 u21 ,u22 ,u23 ,ux  ,uy  ,ufr ,uz  ,uab ,ui  ,u30 , u31 ,u32 , &
                 hm  ,ho  ,hr  ,hsg ,hsb ,hl  ,ha  ,hf  ,hs  ,hn  ,  &
                 hg  ,he  ,hp  ,hne ,hfe ,hd  ,hps ,hsa ,hk  ,hpa ,  &
                 h21 ,h22 ,h23 ,hx  ,hy  ,hfr ,hz  ,hab ,hi  ,h30 , h31 ,h32
 
!     -------------------     / SOFCOM /     ---------------------------
 
!     SOFCOM DEFINES THE NAMES AND SIZES OF THE SOF FILES AND THE STATE
!     OF THE SOF
!     NFILES = NUMBER OF FILES ALLOCATED TO THE SOF (MAX 10)
!     FILNAM = 4 CHAR. BCD NAMES OF THE SOF FILES
!     FILSIZ = SIZES OF THE SOF FILES EXPRESSED IN AN EVEN NUMBER OF
!              BLOCKS
!     STATUS = SOF STATUS.  0 - SOF IS EMPTY.  1 - SOF IS NOT EMPTY.
!     PSSWRD = BCD PASSWORD FOR THE SOF.  EACH RUN USING THE SAME SOF
!              MUST USE THE SAME PASSWORD.
!     FIRST  = .TRUE. IF SOFINT HAS NOT YET BEEN CALLED TO INITIALIZE
!              THE SOF FOR THIS RUN.  OTHERWISE .FALSE.
!     OPNSOF = .TRUE. IF THE SOF IS OPEN.  .FALSE. IF IT IS CLOSED.
!     ASOFCB = ADDRESS OF SOF CONTROL BLOCKS ON IBM 360/370 COMPUTERS
 
 COMMON /sofcom/ nfiles, filnam(10), filsiz(10), STATUS, psswrd(2),  &
                 first , opnsof    , asofcb
 
!     --------------------    / XXREAD /     ---------------------------
 
!     INFLAG AND INSAVE ARE USED IN READFILE COMMAND. IRRX USED IN
!     FFREAD
 
 COMMON /xxread/ inflag, insave, ixxr(3)
 
!     --------------     /XECHOX/ AND /XREADX/     ---------------------
 
!     IECHO      = USED IN FREE-FIELD INPUT FOR INPUT CARD ECHO CONTROL
!     IXSORT,IWASFF,NCARD = USED LOCALLY AMONG XSORT, XREAD, AND FFREAD
!     NOECHO     = USED IN FFREAD AND XCSA ROUTINES
 
!     SCREEN,PROM= LOGICAL UNIT FOR TERMINAL SCREEN AND PROMPT SYMBOL
!     CONTRL NOTYET,STAR,PCT = FREE-FIELD INPUT FLAGS USED IN XREAD
!     LOOP,KOUNT = LOOP COUNT USED IN XREAD
!     ICONT      = 36 CONTROL WORDS USED IN FREE-FILED INPUT NOT TO BE
!                  DESTROYED
 
 COMMON /xechox/ iecho(4), ixsort, iwasff, ncard(3), noecho
 COMMON /xreadx/ screen  , loop  , kount , prom    , notyet, star,  &
                 pct     , icont(36)
 
!     --------------     /MACHIN/ AND /LHPWX/     ---------------------
 
!     6 MACHINE CONSTANTS IN /MACHIN/ AND 7 IN /LHPWX/ WILL BE
!     INITIALZED BY BTSTRP. THESE CONSTANTS NEED TO BE SAVED IN THE ROOT
!     LEVEL OF ALL LINKS
 
 COMMON /machin/ ma(6)
 COMMON /lhpwx / lh(7)
 
 
!     ==================================================================
 
!     -------------------     / GINOX  /     ---------------------------
 
 DATA    cdc   / 244*0   /
 DATA    others/ 392*0   /
 
!     -------------------     / XMSSG  /     ---------------------------
!                               1         2         3
!                      1234567890123456789012345678901
 DATA    ufm   / '0*** USER FATAL MESSAGE'        /
 DATA    uwm   / '0*** USER WARNING MESSAGE'      /
 DATA    uim   / '0*** USER INFORMATION MESSAGE'  /
 DATA    sfm   / '0*** SYSTEM FATAL MESSAGE'      /
 DATA    swm   / '0*** SYSTEM WARNING MESSAGE'    /
 DATA    sim   / '0*** SYSTEM INFORMATION MESSAGE'/
 
!     -------------------     /NUMTPX  /     --------------------------
 
 DATA    nbcd  /  0      /
 DATA    bcd   /  19*0   /
 
!     -------------------     /BLANK  /     --------------------------
 
 DATA    iblnk /  96*0, 4H cea, 4HSE e, 4HMPIR, 4HE >  /
 
!     -------------------     / NTIME  /     ---------------------------
 
 DATA    lntime/ 16      /
 DATA    timdta/ 23*0.   /
 
!     USE A NASTRAN BULKDATA=-3 CARD TO ACTIVATE TIME CONSTANTS COMPUTA-
!     TION AND PRINT OUT FROM SUBROUTINES TMTSIO AND TMTSLP.
 
!     EXAMPLE - THE FOLLOWING CARDS CAN BE USED FOR UNIVAC 1100/82 IN
!     MSFC TO ELIMINATE HARDWARE TIME COMPUTATION IN EVERY NASTRAN RUN.
 
!     DATA    TIMDTA/  0.51, 15.73, 15.00, 11.10, 10.00,  2.20,  2.23,
!    1                 4.00,  5.28, 15.90, 19.40,  4.00,  5.80, 16.45,
!    2                20.16,  0.00/
 
 
!     NOTE - STARTING 1991 VERSION, THESE TIMINGS CONSTANTS CAN BE FED
!     ****   DIRECTLY INTO NASTRAN EXECUTABLE VIA THE NASINFO FILE.
!            THUS, EACH COMPUTER CENTER CAN EASILY SUPPLY THE CORRECT
!            TIMING CONSTANTS FOR ITS MACHINE. (SEE THE WRITE-UP IN THE
!            NASINFO FILE)
!            THE 16TH TIMING IS FOR READING STRING BACKWARD. CURRENTLY
!            NOT USED
 
!     -------------------     / XLINK  /     ---------------------------
 
 DATA    lxlink/ 220 /,  maxlnk / 15/,  mxlink / 220*0 /
 
!     -------------------     / SEM    /     ---------------------------
 
 DATA    mask  / 65535/, mask2 , mask3 / 2*0  /,         NAME  /  &
         4HNS01, 4HNS02, 4HNS03, 4HNS04, 4HNS05, 4HNS06, 4HNS07,  &
         4HNS08, 4HNS09, 4HNS10, 4HNS11, 4HNS12, 4HNS13, 4HNS14,  &
         4HNS15, 4HNS16, 4HNS17, 4HNS18, 4HNS19, 4HNS20, 4HNS21,  &
         4HNS22, 4HNS23, 4HNS24, 4HNS25, 4HNS26, 4HNS27, 4HNS28,  &
         4HNS29, 4HNS30 /
 
!     -------------------     / SYSTEM /     ---------------------------
 
 DATA    sysbuf,outtap,nogo  ,intp  ,mpc   ,spc   ,logfl   /  &
         0     ,0     ,0     ,0     ,0     ,0     ,0       /, &
!            USED ONLY IN MSFC, UNIVAC VERSION - LOGFL = 190  
         load  ,nlpp  ,mtemp ,npages,nlines,tlines,mxlins  /  &
         1     ,0     ,0     ,0     ,0     ,0     ,20000   /, &
         date  ,timez ,echof ,plotf ,apprch,linkno,lsystm  /  &
         3*0   ,0     ,2     ,0     ,0     ,0     ,180     /, &
         icfiat,rfflag,cppgct,mn    ,dummyi,maxfil,maxopn  /  &
         11    ,0     ,0     ,0     ,0     ,35    ,16      /, &
         hicore,timew ,ofpflg,nbrcbu,lprus ,nprus ,ksys37  /  &
         85000 ,0     ,0     ,15    ,64    ,0     ,0       /, &
!        VAX: HICORE IS SET TO 50,000 BY BTSTRP  &
         qq    ,nbpc  ,nbpw  ,ncpw  ,sysdat,tapflg,adumel  /  &
         0     ,0     ,0     ,0     ,3*0   ,0     ,9*0     /, &
         iprec ,ithrml,modcom,hdy   ,sscell,tolel ,mesday  /  &
         0     ,0     ,9*0   ,3*0   ,0     ,0.01  ,0       /
 DATA    bitpas,pass  ,itime ,ctime ,nosbe ,bandit,pzel    /  &
         2*.false.    ,0     ,0     ,0     ,0     ,0       /, &
         switch,icpflg,jrun  ,jmax  ,lintc ,intra ,ospcnt  /  &
         3*0   ,0     ,0     ,0     ,800   ,0     ,15      /, &
         k88_90,lpch  ,ldict ,iaerot,ksys94,sperlk,left  ,loglin,left2 /  &
         3*0   ,0     ,0     ,0     ,0     ,0     ,56*0  ,0     ,28*0  /
 
!     -------------------     / XFIST  /     ---------------------------
 
 DATA    nfist / 56 /, lfist /  56 /,  fist /  &
         4HPOOL,  0,4HOPTP, -1,4HNPTP, -2,4HUMF , -3,4HNUMF, -4,  &
         4HPLT1, -5,4HPLT2, -6,4HINPT, -7,4HINP1, -8,4HINP2, -9,  &
         4HINP3,-10,4HINP4,-11,4HINP5,-12,4HINP6,-13,4HINP7,-14,  &
         4HINP8,-15,4HINP9,-16,4HXPTD,-17,4HSOF ,-18,4HUT1 ,-19,  &
         4HUT2 ,-20,4HUT3 ,-21,4HUT4 ,-22,4HUT5 ,-23,             &
 
!     USE VALUES BELOW WHEN ICFIAT (24TH WORD OF /SYSTEM/) IS 8
!    6           201,  3,   202, 11,   203, 19,   204, 27,   205, 35,
!    7        4HCASE, 43,   207, 51,4HPCDB, 59,   208, 67,   209, 75,
!    8           210, 83,   211, 91,   213, 99,   214,107,   215,115,
!    9           216,123,   301,131,   302,  3,   303, 11,   304, 19,
!    O           305, 27,   306, 35,4HXYCD,139,   307, 67,   308, 75,
!    1           309, 83,   310, 91,   311, 99,   312,107,   313,115,
!    2           314,123,   315,147/
 
!     USE VALUES BELOW WHEN ICFIAT (24TH WORD OF /SYSTEM/) IS 11  &
         201,  3,   202, 14,   203, 25,   204, 36,   205, 47,     &
         4HCASE, 58,   207, 69,4HPCDB, 80,   208, 91,   209,102,  &
         210,113,   211,124,   213,135,   214,146,   215,157,     &
         216,168,   301,179,   302,  3,   303, 14,   304, 25,     &
         305, 36,   306, 47,4HXYCD,190,   307, 91,   308,102,     &
         309,113,   310,124,   311,135,   312,146,   313,157,     &
         314,168,   315,201/
 
!     -------------------     / XPFIST /     ---------------------------
 
 DATA    npfist/ 24   /
 
!     -------------------     / XXFIAT /     ---------------------------
 
 DATA    xxfiat/ 24*0 /
 
!     -------------------     / XFIAT  /     ---------------------------
 
!     USE 8*0 EACH INSTEAD OF 5*0 WHEN ICFIAT = 11
 
 DATA    mfiat / 0 /, nfiat / 100/, lfiat / 0 /, fiat /  &
     0, 4HGEOM, 4H1   , 8*0 , 0, 4HEPT , 4H    , 8*0 ,  &
     0, 4HMPT , 4H    , 8*0 , 0, 4HEDT , 4H    , 8*0 ,  &
     0, 4HDIT , 4H    , 8*0 , 0, 4HCASE, 4HCC  , 3*7 , 2*0 , 3*7 ,  &
     0, 4HDYNA, 4HMICS, 8*0 , 0, 4HPCDB, 4H    , 8*0 ,  &
     0, 4HGEOM, 4H2   , 8*0 , 0, 4HGEOM, 4H3   , 8*0 ,  &
     0, 4HGEOM, 4H4   , 8*0 , 0, 4HGEOM, 4H5   , 8*0 ,  &
     0, 4HFORC, 4HE   , 8*0 , 0, 4HMATP, 4HOOL , 8*0 ,  &
     0, 4HAXIC, 4H    , 8*0 , 0, 4HIFPF, 4HILE , 8*0 ,  &
     0, 4HSCRA, 4HTCH1, 8*0 , 0, 4HXYCD, 4HB   , 8*0 ,  &
     0, 4HSCRA, 4HTC15, 8*0 , 671*0, 220*0 /
 
!     -------------------     / OSCENT /     ---------------------------
 
 DATA    oscar / 200*4H    /
 
!     -------------------     / OUTPUT /     ---------------------------
 
 DATA    output/ 224*4H    /
 
!     -------------------     / XDPL   /     ---------------------------
 
 DATA    mdpl  / 1 /,  ndpl / 80 /,  ldpl  / 0 /,  dpl / 240*0 /
 
!     -------------------     / XVPS   /     ---------------------------
 
 DATA    vps   / 600,  2,   598*0    /
 
!     -------------------     / STAPID /     ---------------------------
 
 DATA    tapid / 6*0.0 /,   otapid  / 6*0.0 /
 DATA    idumf /   0   /
 
!     -------------------     / STIME  /     ---------------------------
 
 DATA    time  / 2*0.0 /
 
!     -------------------     / XCEITB /     ---------------------------
 
 DATA    cei   / 42,   2,   40*0  /
 
!     -------------------     / XMDMSK /     ---------------------------
 
 DATA    nmskcd, nmskfl, nmskrf, msk / 3, 3, 1, 7*0  /
 
!     -------------------     / MSGX   /     ---------------------------
 
 DATA    nmsg  / 0 /,  msglg /  40 /,  msg / 160*0   /
 
!     -------------------     / DESCRP /     ---------------------------
 
!     DATA    LENGTH / 50 /,   BLOCK / 50*0 /
 
!     -------------------     / TWO    /     ---------------------------
 
!     TWO(1)  = LSHIFT(1,31), IS MACHINE DEPENDENT (SET BY BTSTRP)
!     MZERO   = WILL BE SET TO LSHIFT(1,NBPW-1) BY BTSTRP
 
 DATA    two   /  0, 1073741824, 536870912, 268435456, 134217728,  &
                       67108864,  33554432,  16777216,   8388608,  &
                      4194304,   2097152,   1048576,    524288,    &
                      262144,    131072,     65536,     32768,     &
                      16384,      8192,      4096,      2048,      &
                      1024,       512,       256,       128,       &
                      64,        32,        16,         8, 4,      &
                      2,         1/
 DATA    mzero / 0 /
 
!     -------------------     / NAMES  /     ---------------------------
 
 DATA    rd    / 2 /,  rdrew / 0 /,  wrt   / 3 /,  wrtrew/ 1 /,  &
         rew   / 1 /,  norew / 2 /,  eofnrw/ 3 /,  rsp   / 1 /,  &
         rdp   / 2 /,  csp   / 3 /,  cdp   / 4 /,  square/ 1 /,  &
         rect  / 2 /,  diag  / 3 /,  lower / 4 /,  upper / 5 /,  &
         sym   / 6 /,  row   / 7 /,  ident / 8 /
 
!     -------------------     / TYPE   /     ---------------------------
 
 DATA    prc   /     1,     2/      , nwds  /     1,     2,     2,     4/,  &
         rc    /     1,     1,     2,     2/, x     /     6*0.0   /
 
!     -------------------     / BITPOS /     ---------------------------
 
 DATA    um  /32/  ,  hm  /2HM /  ,     ups /16/  ,  hps /2HPS/  ,  &
         uo  /30/  ,  ho  /2HO /  ,     usa /15/  ,  hsa /2HSA/  ,  &
         ur  /29/  ,  hr  /2HR /  ,     uk  /14/  ,  hk  /2HK /  ,  &
         usg /23/  ,  hsg /2HSG/  ,     upa /13/  ,  hpa /2HPA/  ,  &
         usb /22/  ,  hsb /2HSB/  ,     u21 /10/  ,  h21 /4HXXXX/,  &
         ul  /24/  ,  hl  /2HL /  ,     u22 /11/  ,  h22 /4HYYYY/,  &
         ua  /25/  ,  ha  /2HA /  ,     u23 /12/  ,  h23 /4HZZZZ/,  &
         uf  /26/  ,  hf  /2HF /  ,     ux  / 9/  ,  hx  /2HX /  ,  &
         us  /31/  ,  hs  /2HS /  ,     uy  / 8/  ,  hy  /2HY /  ,  &
         un  /27/  ,  hn  /2HN /  ,     ufr / 7/  ,  hfr /2HFR/  ,  &
         ug  /28/  ,  hg  /2HG /  ,     uz  / 6/  ,  hz  /2HZ /  ,  &
         ue  /21/  ,  he  /2HE /  ,     uab / 5/  ,  hab /2HAB/  ,  &
         up  /20/  ,  hp  /2HP /  ,     ui  / 4/  ,  hi  /2HI /  ,  &
         une /19/  ,  hne /2HNE/  ,     u30 / 3/  ,  h30 /2HU3/  ,  &
         ufe /18/  ,  hfe /2HFE/  ,     u31 / 2/  ,  h31 /2HU2/  ,  &
         ud  /17/  ,  hd  /2HD /  ,     u32 / 1/  ,  h32 /2HU1/
 
!     -------------------     / SOFCOM /     ---------------------------
 
 DATA    nfiles/ 1          /
 DATA    filnam/ 4HINPT,9*0 /
 DATA    filsiz/ 100   ,9*0 /
 DATA    STATUS/ 1          /
 DATA    psswrd/ 2*4H       /
 DATA    first / .true.     /
 DATA    opnsof/ .false.    /
 DATA    asofcb/ 0          /
 
!     --------------------    / XXREAD /     ---------------------------
 
 DATA    inflag, insave, ixxr / 5*0 /
 
!     --------------     /XECHOX/ AND /XREADX/     ---------------------
 
 DATA    iecho , ixsort, iwasff, ncard, noecho                /  &
         4*0,    0,      0,      3*0,   0                     /
 DATA    screen, loop, kount, prom, notyet, star, pct, icont  /  &
         6,      -1,   0,     0,    3*.false.,         36*0   /
 
!     --------------     /MACHIN/ AND /LHPWX/     ---------------------
 
 DATA    ma    / 6*0 /,  lh / 7*0 /
 
 END