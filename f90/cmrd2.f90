SUBROUTINE cmrd2
     
!     THIS SUBROUTINE IS THE CMRED2 MODULE WHICH PERFORMS THE MAJOR
!     COMPUTATIONS FOR THE COMPLEX MODAL REDUCE COMMAND.
 
!     DMAP CALLING SEQUENCE
!     CMRED2   CASECC,LAMAMR,PHISSR,PHISSL,EQST,USETMR,KAA,MAA,BAA,K4AA,
!              PAA/KHH,MHH,BHH,K4HH,PHH,POVE/STEP/S,N,DRY/POPT $
 
!     INPUT  DATA
!     GINO - CASECC - CASE CONTROL DATA
!            LAMAMR - EIGENVALUE TABLE FOR SUBSTRUCTURE BEING REDUCED
!            PHISSR - RIGHT HAND EIGENVECTORS FOR SUBSTRUCTURE BEING
!                     REDUCED
!            PHISSL - LEFT HAND EIGENVECTORS FOR SUBSTRUCTURE BEING
!                     REDUCED
!            EQST   - EQSS DATA FOR BOUNDARY SET FOR SUBSTRUCTURE BEING
!                     REDUCED
!            USETMR - USET TABLE FOR REDUCED SUBSTRUCTURE
!            KAA    - SUBSTRUCTURE STIFFNESS MATRIX
!            MAA    - SUBSTRUCTURE MASS MATRIX
!            BAA    - SUBSTRUCTURE VISCOUS DAMPING MATRIX
!            K4AA   - SUBSTRUCTURE STRUCTURE DAMPINF MATRIX
!            PAA    - SUBSTRUCTURE LOAD MATRIX
!     SOF  - LAMS   - EIGENVALUE TABLE FOR ORIGINAL SUBSTRUCTURE
!            PHIS   - RIGHT HAND EIGENVECTOR TABLE FOR ORIGINAL
!                     SUBSTRUCTURE
!            PHIL   - LEFT HAND EIGENVECTOR TABLE FOR ORIGINAL
!                     SUBSTRUCTURE
!            HORG   - RIGHT HAND H TRANSFORMATION MATRIX FOR ORIGINAL
!                     SUBSTRUCTURE
!            HLFT   - LEFT HAND H TRANSFORMATION MATRIX FOR ORIGINAL
!                     SUBSTRUCTURE
 
!     OUTPUT DATA
!     GINO - KHH    - REDUCED STIFFNESS MATRIX
!            MHH    - REDUCED MASS MATRIX
!            BHH    - REDUCED VISCOUS DAMPING MATRIX
!            K4HH   - REDUCED STRUCTURE DAMPING MATRIX
!            PHH    - REDUCED LOAD MATRIX
!            POVE   - INTERIOR POINT LOAD MATRIX
!     SOF  - LAMS   - EIGENVALUE TABLE FOR ORIGINAL SUBSTRUCTURE
!            PHIS   - RIGHT HAND EIGENVECTOR TABLE FOR ORIG.SUBSTRUCTURE
!            PHIL   - LEFT HAND EIGENVECTOR TABLE FOR ORIG. SUBSTRUCTURE
!            GIMS   - G TRANSFORMATION MATRIX FOR BOUNDARY POINTS FOR
!                     ORIGINAL SUBSTRUCTURE
!            HORG   - RIGHT HAND H TRANSFORMATION MATRIX FOR ORIGINAL
!                     SUBSTRUCTURE
!            HLFT   - LEFT HAND H TRANSFORMATION MATRIX FOR ORIGINAL
!                     SUBSTRUCTURE
!            UPRT   - PARTITIONING VECTOR FOR CREDUCE FOR ORIGINAL
!                     SUBSTRUCTURE
!            POVE   - INTERNAL POINT LOADS FOR ORIGINAL SUBSTRUCTURE
!            POAP   - INTERNAL POINTS APPENDED LOADS FOR ORIGINAL
!                     SUBSTRUCTURE
!            EQSS   - SUBSTRUCTURE EQUIVALENCE TABLE FOR REDUCED
!                     SUBSTRUCTURE
!            BGSS   - BASIC GRID POINT DEFINITION TABLE FOR REDUCED
!                     SUBSTRUCTURE
!            CSTM   - COORDINATE SYSTEM TRANSFORMATION MATRICES FOR
!                     REDUCED SUBSTRUCTURE
!            LODS   - LOAD SET DATA FOR REDUCED SUBSTRUCTURE
!            LOAP   - APPENDED LOAD SET DATA FOR REDUCED SUBSTRUCTURE
!            PLTS   - PLOT SET DATA FOR REDUCED SUBSTRUCTURE
!            KMTX   - STIFFNESS MATRIX FOR REDUCED SUBSTRUCTURE
!            MMTX   - MASS MATRIX FOR REDUCED SUBSTRUCTURE
!            PVEC   - LOAD MATRIX FOR REDUCED SUBSTRUCTURE
!            PAPD   - APPENDED LOAD MATRIX FOR REDUCED SUBSTRUCTURE
!            BMTX   - VISCOUS DAMPING MATRIX FOR REDUCED SUBSTRUCTURE
!            K4MX   - STRUCTURE DAMPING MATRIX FOR REDUCED SUBSTRUCTURE
 
!     PARAMETERS
!     INPUT  - STEP   - CONTROL DATA CASECC RECORD (INTEGER)
!              POPT   - PVEC OR PAPP OPTION FLAG (BCD)
!     OUTPUT - DRY    - MODULE OPERATION FLAG (INTEGER)
!     OTHERS - GBUF   - GINO BUFFERS
!              SBUF   - SOF BUFFERS
!              INFILE - INPUT FILE NUMBERS
!              OTFILE - OUTPUT FILE NUMBERS
!              ISCR   - ARRAY OF SCRATCH FILE NUMBERS
!              KORLEN - LENGTH OF OPEN CORE
!              KORBGN - BEGINNING ADDRESS OF OPEN CORE
!              OLDNAM - NAME OF SUBSTRUCTURE BEING REDUCED
!              NEWNAM - NAME OF REDUCED SUBSTRUCTURE
!              SYMTRY - SYMMETRY FLAG
!              RANGE  - RANGE OF FREQUENCIES TO BE USED
!              NMAX   - MAXIMUM NUMBER OF FREQUENCIES TO BE USED
!              IO     - IO OPTIONS FLAG
!              MODES  - OLDMODES OPTION FLAG
!              RSAVE  - SAVE REDUCTION PRODUCT FLAG
!              LAMSAP - BEGINNING ADDRESS OF MODE USE DESCRIPTION ARRAY
!              MODLEN - LENGTH OF MODE USE ARRAY
!              MODPTS - NUMBER OF MODAL POINTS
 
 EXTERNAL        orf
 LOGICAL :: symtry,modes,rsave,ponly
 INTEGER :: step,dry,popt,gbuf1,gbuf2,gbuf3,sbuf1,sbuf2,sbuf3,  &
     otfile,oldnam,z,sysbuf,casecc,yes,phissl,orf
 DIMENSION       modnam(2),nmonic(8),rz(1),itrlr(7)
 COMMON /BLANK / step,dry,popt,gbuf1,gbuf2,gbuf3,sbuf1,sbuf2,sbuf3,  &
     infile(11),otfile(6),iscr(11),korlen,korbgn,  &
     oldnam(2),newnam(2),symtry,range(2),nmax,io,modes,  &
     rsave,lamsap,modpts,modlen,ponly,lstzwd
 COMMON /zzzzzz/ z(1)
 COMMON /system/ sysbuf,iprntr
 EQUIVALENCE     (casecc,infile(1)),(phissl,infile(4)),(rz(1),z(1))
 DATA    nmonic/ 4HNAMA,4HNAMB,4HSYMF,4HRANG,4HNMAX,4HOUTP,4HOLDM, 4HRSAV/
 DATA    kaa   / 107 /, iblank,yes /4H    , 4HYES /
 DATA    modnam/ 4HCMRD,4H2   /
 DATA    nhlods, nhloap,nhhorg,nhhlft /4HLODS,4HLOAP,4HHORG,4HHLFT/
 
!     COMPUTE OPEN CORE AND DEFINE GINO, SOF BUFFERS
 
 IF (dry == -2) RETURN
 nozwds = korsz(z(1))
 lstzwd = nozwds - 1
 gbuf1  = nozwds - sysbuf - 2
 gbuf2  = gbuf1  - sysbuf
 gbuf3  = gbuf2  - sysbuf
 sbuf1  = gbuf3  - sysbuf
 sbuf2  = sbuf1  - sysbuf - 1
 sbuf3  = sbuf2  - sysbuf
 korlen = sbuf3  - 1
 korbgn = 1
 IF (korlen <= korbgn) GO TO 290
 
!     INITIALIZE SOF
 
 CALL sofopn (z(sbuf1),z(sbuf2),z(sbuf3))
 
!     INITIALIZE CASE CONTROL PARAMETERS
 
 DO  i = 1,11
   IF (i > 6) GO TO 2
   infile(i) = 100 + i
   otfile(i) = 200 + i
   iscr(i) = 300 + i
   CYCLE
   2 infile(i) = 100 + i
   iscr(i) = 300 + i
 END DO
 DO  i = 1,2
   oldnam(i) = iblank
   newnam(i) = iblank
 END DO
 range(1) = -1.0E+35
 range(2) =  1.0E+35
 symtry = .false.
 nmax   = 2147483647
 io     = 0
 modes  = .false.
 rsave  = .false.
 nrange = 0
 ponly  = .false.
 
!     PROCESS CASE CONTROL
 
 ifile = casecc
 CALL OPEN (*260,casecc,z(gbuf2),0)
 IF (step == 0.0) THEN
   GO TO    40
 END IF
 20 DO  i = 1,step
   CALL fwdrec (*280,casecc)
 END DO
 
!     READ CASECC
 
 40 CALL READ (*270,*280,casecc,z(korbgn),2,0,nwdsrd)
 nwdscc = z(korbgn+1)
 DO  i = 1,nwdscc,3
   CALL READ (*270,*280,casecc,z(korbgn),3,0,nwdsrd)
   
!     TEST CASE CONTROL MNEMONICS
   
   DO  j = 1,8
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
       GO TO 160
     CASE (    7)
       GO TO 180
     CASE (    8)
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
   
!     EXTRACT SYMMETRY FLAG
   
   110 IF (z(korbgn+1) /= yes) CYCLE
   symtry = .true.
   CYCLE
   
!     EXTRACT FREQUENCY RANGE
   
   120 IF (nrange == 1) GO TO 125
   nrange = 1
   range(1) = rz(korbgn+2)
   CYCLE
   125 range(2) = rz(korbgn+2)
   CYCLE
   
!     EXTRACT MAXIMUM NUMBER OF FREQUENCIES
   
   140 IF (z(korbgn) == 0) CYCLE
   nmax = z(korbgn+2)
   CYCLE
   
!     EXTRACT OUTPUT FLAGS
   
   160 io = orf(io,z(korbgn+2))
   CYCLE
   
!     EXTRACT OLDMODES FLAG
   
   180 IF (z(korbgn+1) /= yes) CYCLE
   modes = .true.
   CYCLE
   
!     EXTRACT REDUCTION SAVE FLAG
   
   190 IF (z(korbgn+1) /= yes) CYCLE
   rsave = .true.
 END DO
 CALL CLOSE (casecc,1)
 
!     CHECK FOR SYMMETRY
 
 itrlr(1) = phissl
 CALL rdtrl (itrlr)
 npass = 2
 IF (itrlr(1) > 0) GO TO 204
 symtry = .true.
 npass = 1
 
!     CHECK FOR RUN = GO
 
 204 ihorg = 0
 IF (dry == 0) GO TO 240
 
!     CHECK FOR STIFFNESS PROCESSING
 
 itrlr(1) = kaa
 CALL rdtrl (itrlr)
 IF (itrlr(1) > 0) GO TO 208
 
!     CHECK FOR LOADS ONLY PROCESSING
 
 CALL sfetch (newnam,nhlods,3,itest)
 IF (itest == 3) ponly = .true.
 CALL sfetch (newnam,nhloap,3,itest)
 IF (itest == 3) ponly = .true.
 GO TO 240
 
!     PROCESS STIFFNESS MATRIX
 
 208 CALL cmrd2a
 
!     BEGIN COMPLEX MODAL REDUCTION
!     NPASS .EQ. 1, SYMMETRIC REDUCTION
!     NPASS .EQ. 2, UNSYMMETRIC REDUCTION
 
 DO  j = 1,npass
   
!     TEST FOR H TRANSFORMATION MATRICES
   
   SELECT CASE ( j )
     CASE (    1)
       GO TO 212
     CASE (    2)
       GO TO 214
   END SELECT
   212 CALL softrl (oldnam,nhhorg,itrlr)
   IF (itrlr(1) == 1) CYCLE
   ihorg = ihorg + 1
   GO TO 216
   214 CALL softrl (oldnam,nhhlft,itrlr)
   IF (itrlr(1) == 1) CYCLE
   ihorg = ihorg + 2
   
!     PREFORM GUYAN REDUCTION
   
   216 CALL cmrd2c (j)
   
!     PROCESS OLDMODES FLAG
   
   CALL cmrd2b (j)
   
!     CALCULATE MODAL TRANSFORMATION MATRIX
   
   CALL cmrd2d (j)
   IF (j == 1) CALL cmrd2b (3)
   
!     CALCULATE H TRANSFORMATION MATRIX
   
   CALL cmrd2e (j)
 END DO
 
!     CALCULATE STRUCTURAL MATRICES
!     IHORG .EQ. 0, BOTH HORG, HLFT ON SOF
!     IHORG .EQ. 1, HORG CALCULATED, HLFT ON SOF
!     IHORG .EQ. 2, HORG ON SOF, HLFT CALCULATED
!     IHORG .EQ. 3, BOTH HORG, HLFT CALCULATED
 
 240 CALL cmrd2f (ihorg)
 IF (ihorg == 0) GO TO 250
 
!     PROCESS NEW TABLE ITEMS
 
 CALL cmrd2g
 
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
 
 900 FORMAT (50H0  module creduce terminating due TO above errors.)
 
END SUBROUTINE cmrd2
