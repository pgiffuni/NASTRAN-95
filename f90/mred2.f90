SUBROUTINE mred2
     
!     THIS SUBROUTINE IS THE MRED2 MODULE WHICH PERFORMS THE MAJOR
!     COMPUTATIONS FOR THE REDUCE COMMAND.
 
!     DMAP CALLING SEQUENCE
!     MRED2    CASECC,LAMAMR,PHISS,EQST,USETMR,KAA,MAA,BAA,K4AA,PAA,DMR,
!              QSM/KHH,MHH,BHH,K4HH,PHH,POVE/STEP/S,N,DRY/POPT $
 
!     12 INPUT DATA BLOCKS
!     GINO -   CASECC - CASE CONTROL DATA
!              LAMAMR - EIGENVALUE TABLE FOR SUBSTRUCTURE BEING REDUCED
!              PHISS  - EIGENVECTORS FOR SUBSTRUCTURE BEING REDUCED
!              EQST   - EQSS DATA FOR BOUNDARY SET FOR SUBSTRUCTURE
!                       BEINGREDUCED
!              USETMR - USET TABLE FOR REDUCED SUBSTRUCTURE
!              KAA    - SUBSTRUCTURE STIFFNESS MATRIX
!              MAA    - SUBSTRUCTURE MASS MATRIX
!              BAA    - SUBSTRUCTURE VISCOUS DAMPING MATRIX
!              K4AA   - SUBSTRUCTURE STRUCTURE DAMPINF MATRIX
!              PAA    - SUBSTRUCTURE LOAD MATRIX
!              DMR    - FREE BODY MATRIX
!              QSM    - MODEL REACTION MATRIX
!     SOF  -   LAMS   - EIGENVALUE TABLE FOR ORIGINAL SUBSTRUCTURE
!              PHIS   - EIGENVECTOR TABLE FOR ORIGINAL SUBSTRUCTURE
!              LMTX   - STIFFNESS DECOMPOSITION PRODUCT FOR ORIGINAL
!                       SUBSTRUCTURE
!              GIMS   - G TRANSFORMATION MATRIX FOR BOUNDARY POINTS FOR
!                       ORIGINAL SUBSTRUCTURE
!              HORG   - H TRANSFORMATION MATRIX FOR ORIGINAL
!                       SUBSTRUCTURE
 
!     6 OUTPUT DATA BLOCKS
!     GINO -   KHH    - REDUCED STIFFNESS MATRIX
!              MHH    - REDUCED MASS MATRIX
!              BHH    - REDUCED VISCOUS DAMPING MATRIX
!              K4HH   - REDUCED STRUCTURE DAMPING MATRIX
!              PHH    - REDUCED LOAD MATRIX
!              POVE   - INTERIOR POINT LOAD MATRIX
!     SOF  -   LAMS   - EIGENVALUE TABLE FOR ORIGINAL SUBSTRUCTURE
!              PHIS   - EIGENVECTOR TABLE FOR ORIGINAL SUBSTRUCTURE
!              LMTX   - STIFFNESS DECOMPOSITION PRODUCT FOR ORIGINAL
!                       SUBSTRUCTURE
!              GIMS   - G TRANSFORMATION MATRIX FOR BOUNDARY POINTS FOR
!                       ORIGINAL SUBSTRUCTURE
!              HORG   - H TRANSFORMATION MATRIX FOR ORIGINAL
!                       SUBSTRUCTURE
!              UPRT   - PARTITIONING VECTOR FOR MREDUCE FOR ORIGINAL
!                       SUBSTRUCTURE
!              POVE   - INTERNAL POINT LOADS FOR ORIGINAL SUBSTRUCTURE
!              POAP   - INTERNAL POINTS APPENDED LOADS FOR ORIGINAL
!                       SUBSTRUCTURE
!              EQSS   - SUBSTRUCTURE EQUIVALENCE TABLE FOR REDUCED
!                       SUBSTRUCTURE
!              BGSS   - BASIC GRID POINT DEFINITION TABLE FOR REDUCED
!                       SUBSTRUCTURE
!              CSTM   - COORDINATE SYSTEM TRANSFORMATION MATRICES FOR
!                       REDUCED SUBSTRUCTURE
!              LODS   - LOAD SET DATA FOR REDUCED SUBSTRUCTURE
!              LOAP   - APPENDED LOAD SET DATA FOR REDUCED SUBSTRUCTURE
!              PLTS   - PLOT SET DATA FOR REDUCED SUBSTRUCTURE
!              KMTX   - STIFFNESS MATRIX FOR REDUCED SUBSTRUCTURE
!              MMTX   - MASS MATRIX FOR REDUCED SUBSTRUCTURE
!              PVEC   - LOAD MATRIX FOR REDUCED SUBSTRUCTURE
!              PAPD   - APPENDED LOAD MATRIX FOR REDUCED SUBSTRUCTURE
!              BMTX   - VISCOUS DAMPING MATRIX FOR REDUCED SUBSTRUCTURE
!              K4MX   - STRUCTURE DAMPING MATRIX FOR REDUCED
!                       SUBSTRUCTURE
 
!     11 SCRATCH DATA BLOCKS
 
!     PARAMETERS
!     INPUT  - STEP   - CONTROL DATA CASECC RECORD (INTEGER)
!              POPT   - PVEC OR PAPP OPTION FLAG (BCD)
!     OUTPUT - DRY    - MODULE OPERATION FLAG (INTEGER)
!     OTHERS - GBUF   - GINO BUFFERS
!              SBUF   - SOF BUFFERS
!              INFILE - INPUT FILE NUMBERS
!              OTFILE - OUTPUT FILE NUMBERS
!              ISCR   - ARRAY OF SCRATCH FILE NUMBERS
!              ISCR11 - LII PARTITION MATRIX USED IN MRED2B AND MRED2F
!              KORLEN - LENGTH OF OPEN CORE
!              KORBGN - BEGINNING ADDRESS OF OPEN CORE
!              OLDNAM - NAME OF SUBSTRUCTURE BEING REDUCED
!              NEWNAM - NAME OF REDUCED SUBSTRUCTURE
!              FREBDY - FREE BODY MODES CALCULATION FLAG
!              RANGE  - RANGE OF FREQUENCIES TO BE USED
!              NMAX   - MAXIMUM NUMBER OF FREQUENCIES TO BE USED
!              USRMOD - USERMODES CALCULATION FLAG
!              IO     - IO OPTIONS FLAG
!              BOUNDS - OLDBOUNDS OPTION FLAG
!              MODES  - OLDMODES OPTION FLAG
!              RSAVE  - SAVE REDUCTION PRODUCT FLAG
!              LAMSAP - BEGINNING ADDRESS OF MODE USE DESCRIPTION ARRAY
!              MODPTS - NUMBER OF MODAL POINTS
!              MODLEN - LENGTH OF MODE USE ARRAY
 
 EXTERNAL        orf
 LOGICAL :: frebdy,bounds,modes,rsave,ponly
 INTEGER :: step,dry,popt,gbuf1,gbuf2,sbuf1,sbuf2,sbuf3,  &
     otfile,oldnam,usrmod,gbuf3,z,sysbuf,casecc,orf
 DIMENSION       modnam(2),nmonic(10),rz(1),itrlr(7)
 COMMON /BLANK / step,dry,popt,gbuf1,gbuf2,gbuf3,sbuf1,sbuf2,sbuf3,  &
     infile(12),otfile(6),iscr(10),korlen,korbgn,  &
     oldnam(2),newnam(2),frebdy,range(2),nmax,usrmod,  &
     io,bounds,modes,rsave,lamsap,modpts,modlen,ponly, lstzwd,iscr11
 COMMON /zzzzzz/ z(1)
 COMMON /system/ sysbuf,iprntr
 EQUIVALENCE     (casecc,infile(1)), (rz(1),z(1))
 DATA    nmonic/ 4HNAMA,4HNAMB,4HFREE,4HRANG,4HNMAX,4HUSER,4HOUTP,  &
     4HOLDB,4HOLDM,4HRSAV/
 DATA    iblank, nhlods,nhloap/4H    ,4HLODS,4HLOAP/
 DATA    modnam/ 4HMRED,4H2   /
 DATA    itrlr / 106   ,6*0   /
 
!     COMPUTE OPEN CORE AND DEFINE GINO, SOF BUFFERS
 
 IF (dry == -2) RETURN
 nozwds = korsz(z(1))
 lstzwd = nozwds- 1
 gbuf1  = nozwds- sysbuf - 2
 gbuf2  = gbuf1 - sysbuf
 gbuf3  = gbuf2 - sysbuf
 sbuf1  = gbuf3 - sysbuf
 sbuf2  = sbuf1 - sysbuf - 1
 sbuf3  = sbuf2 - sysbuf
 korlen = sbuf3 - 1
 korbgn = 1
 IF (korlen <= korbgn) GO TO 290
 
!     INITIALIZE SOF
 
 CALL sofopn (z(sbuf1),z(sbuf2),z(sbuf3))
 
!     INITIALIZE CASE CONTROL PARAMETERS
 
 DO  i = 1,12
   infile(i) = 100 + i
 END DO
 DO  i = 1,6
   otfile(i) = 200 + i
 END DO
 DO  i = 1, 10
   iscr(i) = 300 + i
 END DO
 iscr11  = 311
 DO  i = 1, 2
   oldnam(i) = iblank
   newnam(i) = iblank
 END DO
 range(1) = 0.0
 range(2) = 1.0E+35
 frebdy = .false.
 nmax   = 2147483647
 usrmod = -1
 io     = 0
 nrange = 0
 bounds = .false.
 modes  = .false.
 rsave  = .false.
 ponly  = .false.
 
!     ** PROCESS CASE CONTROL
 
 ifile  = casecc
 CALL OPEN (*260,casecc,z(gbuf2),0)
 IF (step == 0.0) THEN
   GO TO    40
 END IF
 20 DO  i = 1, step
   CALL fwdrec (*280,casecc)
 END DO
 
!     READ CASECC
 
 40 CALL READ (*270,*280,casecc,z(korbgn),2,0,nwdsrd)
 nwdscc = z(korbgn+1)
 DO  i = 1,nwdscc,3
   CALL READ (*270,*280,casecc,z(korbgn),3,0,nwdsrd)
   
!     TEST CASE CONTROL MNEMONICS
   
   DO  j = 1,10
     IF (z(korbgn) == nmonic(j)) GO TO 60
   END DO
   CYCLE
   
!     SELECT DATA TO EXTRACT
   
   60 SELECT CASE ( j )
     CASE (    1)
       GO TO 70
     CASE (    2)
       GO TO 90
     CASE (    3)
       GO TO 110
     CASE (    4)
       GO TO 120
     CASE (    5)
       GO TO 140
     CASE (    6)
       GO TO 150
     CASE (    7)
       GO TO 160
     CASE (    8)
       GO TO 170
     CASE (    9)
       GO TO 180
     CASE (   10)
       GO TO 190
   END SELECT
   
!     EXTRACT NAME OF SUBSTRUCTURE BEING REDUCED
   
   70 DO  k = 1,2
     oldnam(k) = z(korbgn+k)
   END DO
   CYCLE
   
!     EXTRACT NAME OF REDUCED SUBSTRUCTURE
   
   90 DO  k = 1,2
     newnam(k) = z(korbgn+k)
   END DO
   CYCLE
   
!     EXTRACT FREEBODY MODES FLAG
   
   110 frebdy = .true.
   CYCLE
   
!     EXTRACT FREQUENCY RANGE
   
   120 IF (nrange == 1) GO TO 130
   nrange = 1
   range(1) = rz(korbgn+2)
   CYCLE
   130 range(2) = rz(korbgn+2)
   CYCLE
   
!     EXTRACT MAXIMUM NUMBER OF FREQUENCIES
   
   140 IF (z(korbgn+2) == 0) CYCLE
   nmax = z(korbgn+2)
   CYCLE
   
!     EXTRACT USERMODE FLAG
   
   150 usrmod = z(korbgn+2)
   CYCLE
   
!     EXTRACT OUTPUT FLAGS
   
   160 io = orf(io,z(korbgn+2))
   CYCLE
   
!     EXTRACT OLDBOUND FLAG
   
   170 bounds = .true.
   CYCLE
   
!     EXTRACT OLDMODES FLAG
   
   180 modes = .true.
   CYCLE
   
!     EXTRACT REDUCTION SAVE FLAG
   
   190 rsave = .true.
   
 END DO
 CALL CLOSE (casecc,1)
 
!     TEST FOR RUN = GO
 
 mrd2g = 1
 IF (dry == 0) GO TO 230
 
!     CHECK FOR USERMODE = TYPE 2
 
 IF (usrmod == 2) GO TO 210
 
!     CHECK FOR STIFFNESS PROCESSING
 
 CALL rdtrl (itrlr)
 IF (itrlr(1) > 0) GO TO 208
 
!     CHECK FOR LOADS ONLY
 
 CALL sfetch (newnam,nhlods,3,itest)
 IF (itest == 3) GO TO 204
 CALL sfetch (newnam,nhloap,3,itest)
 IF (itest == 3) GO TO 204
 mrd2g = 4
 GO TO 230
 204 mrd2g = 3
 ponly = .true.
 GO TO 230
 
!     PROCESS STIFFNESS MATRIX
 
 208 mrd2g = 2
 CALL mred2a
 
!     PROCESS OLDBOUND FLAG
 
 CALL mred2b
 
!     PROCESS OLDMODES FLAG
 
 CALL mred2c (1)
 GO TO 220
 
!     PROCESS USERMODES FLAG
 
 210 CALL mred2d
 CALL mred2c (3)
 GO TO 240
 
!     CALCULATE MODAL TRANSFORMATION MATRIX
 
 220 CALL mred2e
 CALL mred2c (2)
 
!     CALCULATE FREE BODY EFFECTS
 
 CALL mred2f
 
!     CALCULATE STRUCTURAL MATRICES
 
!     MRD2G .EQ. 1, M,B,K4,P/PA PROCESSING (RUN = GO)
!     MRD2G .EQ. 2, K,M,B,K4,P/PA PROCESSING
!     MRD2G .EQ. 3, P/PA PROCESSING (ONLY)
!     MRD2G .EQ. 4, M,B,K4,P/PA PROCESSING (RUN = STEP)
 
 230 CALL mred2g (mrd2g)
 IF (mrd2g == 1) GO TO 250
 
!     PROCESS NEW TABLE ITEMS
 
 240 CALL mred2h
 
!     CLOSE ANY OPEN FILES
 
 250 CALL sofcls
 IF (dry == -2) WRITE (iprntr,900)
 RETURN
 
!     PROCESS SYSTEM FATAL ERRORS
 
 260 imsg = -1
 GO TO 300
 270 imsg = -2
 GO TO 300
 280 imsg = -3
 GO TO 300
 290 imsg = -8
 ifile = 0
 300 CALL sofcls
 CALL mesage (imsg,ifile,modnam)
 RETURN
 
 900 FORMAT (//,'  MODULE MREDUCE TERMINATING DUE TO ABOVE ERRORS.')
 
END SUBROUTINE mred2
