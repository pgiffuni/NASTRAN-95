SUBROUTINE invpwr
     
!     GIVEN A REAL SYMETRIC MATRIX, INVPWR WILL SOLVE FOR ALL OF THE
!     EIGENVALUES AND EIGENVECTORS WITHIN A SPECIFIED RANGE
 
!     DEFINITION OF INPUT AND OUTPUT PARAMETERS
 
!     FILEK(7) =  MATRIX CONTROL BLOCK FOR THE INPUT STIFFNESS MATRIX K
!     FILEM(7) =  MATRIX CONTROL BLOCK FOR THE INPUT MASS MATRIX M
!     FILELM(7)=  MATRIX CONTROL BLOCK FOR THE OUTPUT EIGENVALUES
!     FILEVC(7)=  MATRIX CONTROL BLOCK FOR THE OUTPUT EIGENVECTORS
!     SR1FIL-
!     SR7FIL   =  SCRATCH FILES REQUIRED INTERNALLY
!     LAMMIN   =  MINIMUM VALUE FOR THE EIGENVALUE
!     LAMMAX   =  MAXIMUM VALUE FOR THE EIGENVALUE
!     NOEST    =  NUMBER OF ESTIMATED EIGENVALUES WITHIN THE SPECIFIED
!                 RANGE
!     NDPLUS   =  NUMBER OF DESIRED EIGENVALUES IN THE POSITIVE RANGE
!     NDMNUS   =  NUMBER OF DESIRED EIGENVALUES IN THE NEGATIVE RANGE
!     EPS      =  CONVERGENCE CRITERIA
 
!     FILELM AND FILEVC WILL BE USED AS SR1FIL AND SR2FIL WHILE THE
!     EIGENVALUES AND EIGENVECTORS WILL BE STORED ON THE ACTUAL SR1FIL
!     AND SR2FIL. THE ORDERING OF THE EIGENVALUES AND EIGENVECTORS WILL
!     PUT THEM ON FILELM AND FILEVC IN THE CORRECT SEQUENCE AT THE END
!     OF THE SUBROUTINE
 
!     SR1FIL-FILELM CONTAINS (K-LAMBDA*M)
!     SR2FIL-FILEVC CONTAINS THE LOWER TRIANGLE L
!     SR3FIL        CONTAINS THE UPPER TRIANGLE U
!     SR4FIL        IS USED AS SCRATCH IN DECOMP
!     SR5FIL        IS USED AS SCRATCH IN DECOMP
!     SR6FIL        IS USED AS SCRATCH IN DECOMP
!     SR7FIL        CONTAINS THE VECTORS WHICH ARE USED TO ORTHOGONALIZE
!                   THE CURRENT ITERATE
 
 EXTERNAL          norm11    ,sub1     ,mtmsu1   ,xtrny1   ,  &
     norm1     ,sub      ,mtimsu   ,xtrnsy
 INTEGER :: sysbuf    ,comflg   ,filek    ,NAME(2)  ,  &
     switch    ,dmpfil   ,iz(12)   ,sturm    , t1        ,t2       ,timed
 INTEGER :: wrtrew    ,rew      ,sr1fil   ,sr2fil
 REAL :: lammin    ,lammax   ,lmin     ,z        , zz(1)
 DOUBLE PRECISION :: lambda    ,lmbda
 COMMON   /dcompx/ dumxx(35) ,isym
 COMMON   /invpwx/ filek(7)  ,filem(7) ,filelm(7),filevc(7),  &
     sr1fil    ,sr2fil   ,sr3fil   ,sr4fil   ,  &
     sr5fil    ,sr6fil   ,sr7fil   ,sr8fil   ,  &
     dmpfil    ,lammin   ,lammax   ,noest    ,  &
     ndplus    ,ndmnus   ,eps      ,northo
 COMMON   /sturmx/ sturm     ,shftpt   ,KEEP(2)
 COMMON   /invpxx/ lambda    ,comflg   ,iter     ,timed    ,  &
     nopos     ,rzero    ,neg      ,nochng   ,  &
     ind       ,lmbda    ,switch   ,nzero    ,  &
     noneg     ,ivect    ,ireg     ,istart
 COMMON   /system/ ksystm(65)
 COMMON   /names / rd        ,rdrew    ,wrt      ,wrtrew   ,  &
     rew       ,norew    ,eofnrw
 COMMON   /zzzzzz/ z(1)
 EQUIVALENCE       (zz(1),z(1))
 EQUIVALENCE       (iz(1),z(1)), (ksystm(1),sysbuf), (ksystm(55),iprec)
 DATA      NAME  / 4HINVP,4HWR  /
 
!     DEFINITION OF INTERNAL PARAMETERS
 
!     NSHIFT =  NUMBER OF SHIFT POINTS
!     ISHIFT =  CURRENT SHIFT REGION
!     NOVECT =  NUMBER OF EIGENVECTORS FOUND IN A GIVEN REGION
!     NOSKIP =  NUMBER OF VECTORS TO SKIP TO REACH THE LAST SHIFT REGION
!     NEG    =  1 = FIND NEGATIVE ROOTS
!               0 = FIND ONLY POSITIVE ROOTS
!              -1 = WE ARE NOW SEARCHING FOR THE NEGATIVE ROOTS
!     LAMBDA =  THE CURRENT SHIFT POINT
!     RZEROP =  THE CURRENT EIGENVALUE MUST BE .LT. LAMBDA + RZEROP
!     RZEROM =  THE CURRENT EIGENVALUE MUST BE .GT. LAMBDA - RZEROM
!     LMBDA  =  THE ORIGINAL VALUE OF LAMBDA IN A GIVEN REGION
!     COMFLG =  0 = INITIAL ENTRY WITH NEW LAMBDA
!               1 = NEW SHIFT POINT WITHIN THE SEARCH REGION
!               2 = NEW SHIFT DUE TO CLOSENESS TO AN EIGENVALUE
!               3 = NUMBER OF DESIRED POSITIVE ROOTS FOUND
!               4 = NUMBER FOUND EXCEEDS 3*NOEST
!     ISING  =  SINGULARITY FLAG  0 = NO SINGULARITY
!                                 1 = SINGULAR MATRIX - CHANGE LAMBDA
!                                     AND TRY ONE MORE TIME
!     ITER   =  TOTAL NUMBER OF ITERATIONS
!     NOCHNG =  NUMBER OF SHIFTS WITHIN ONE REGION
!     TIMED  =  TIME REQUIRED TO FORM AND DECOMPOSE (K-LAMBDA*M)
!     NFIRST =  NUMBER OF VECTORS IN THE FIRST POSITIVE SEARCH REGION
 
 isym   = 1
 nshift = (noest+5)/6
 mxchng = MAX0 (10, nshift)
 ncol   = filek(2)
 ncol2  = 2*ncol
 ishift = 1
 nz     = korsz(zz(1))
 icrq   = ncol*(1+7*iprec) + 4*sysbuf - nz
 IF (icrq > 0) GO TO 220
 nz     = korsz(z(1))
 ibuf1  = nz - sysbuf
 icrq   = ncol2 - ibuf1
 IF (ibuf1 <= ncol2) GO TO 220
 nopos  = northo
 noneg  = 0
 neg    = 0
 ind    = 0
 iter   = 0
 nodcmp = 0
 nostrt = 0
 nomovs = 0
 IF (northo > 0) GO TO 20
 CALL gopen (sr1fil,z(ibuf1),wrtrew)
 CALL CLOSE (sr1fil,norew)
 CALL gopen (sr2fil,z(ibuf1),wrtrew)
 CALL CLOSE (sr2fil,norew)
 20 lmin   = lammin
 IF (lammin >= 0.0) GO TO 30
 lmin   = 0.
 neg    = 1
 IF (lammax > 0.0) GO TO 30
 lmin   = lammax
 neg    =-1
 dellam = lammin - lammax
 GO TO 40
 
!     EVALUATE THE VALUE OF LAMBDA IN THE CENTER OF THE CURRENT SEARCH
!     REGION
 
 30 dellam = lammax - lmin
 40 lambda = lmin + (ishift - 0.5)*dellam/nshift
 rzero  = ABS(0.55*dellam/nshift)
 nostrt = nostrt + 1
 50 comflg = 0
 lmbda  = lambda
 
!     INITIATE CLOCK TIME
 
 CALL klock (istart)
 nochng = 0
 switch = 0
 ivect  = 0
 ireg   = 0
 ind    = ind + 1
 IF (IABS(ind) == 13) ind = 1
 ising  = 0
 GO TO 90
 70 ising  = 0
 switch = 1
 90 IF (nochng >= mxchng) GO TO 160
 nochng = nochng + 1
 CALL klock (t1)
 
!     CALL IN ADD LINK TO FORM  (K-LAMBDA*M)
 
 CALL invp1
 
!     CALL IN DECOMP TO DECOMPOSE THIS MATRIX
 
 nodcmp = nodcmp + 1
 shftpt = lambda
 CALL invp2 (*100)
 CALL klock (t2)
 GO TO 110
 
!     SINGULAR MATRIX. INCREMENT LAMBDA AND TRY ONCE MORE
 
 100 IF (ising == 1) GO TO 150
 ising  = 1
 lambda = lambda + .02*rzero
 GO TO 90
 
!     DETERMINE THE TIME REQUIRED TO FORM AND DECOMPOSE (K-LAMBDA*M)
 
 110 timed  = t2 - t1
 
!     CALL IN THE MAIN LINK TO ITERATE FOR EIGENVALUES
 
 IF (iprec  == 1) CALL invp3 (norm11,sub1,mtmsu1,xtrny1)
 IF (iprec  == 2) CALL invp3 (norm1 ,sub ,mtimsu,xtrnsy)
 IF (comflg == 2) GO TO 200
 IF (comflg == 1) GO TO 70
 IF (comflg == 3) GO TO 130
 IF (comflg == 0) GO TO 120
 GO TO 170
 120 ishift = ishift + 1
 IF (ishift > nshift) GO TO 130
 GO TO 40
 130 IF (neg > 0) THEN
   GO TO   140
 ELSE
   GO TO   180
 END IF
 
!     INITIALIZE PARAMETERS TO SOLVE FOR NEGATIVE EIGENVALUES
 
 140 x      = nshift*(-lammin/lammax)
 ix     = x
 y      = ix
 IF (x /= y) ix = ix + 1
 nshift = ix
 neg    =-1
 dellam = lammin
 ishift = 1
 GO TO 40
 150 iterm  = 1
 GO TO 190
 160 iterm  = 2
 GO TO 190
 170 iterm  = comflg
 GO TO 190
 180 iterm  = 3
 
!     RE-ORDER EIGENVALUES AND EIGENVECTORS
 
 190 CALL gopen (dmpfil,z(ibuf1),wrtrew)
 iz( 1) = 2
 iz( 2) = northo
 iz( 3) = nostrt
 iz( 4) = nomovs
 iz( 5) = nodcmp
 iz( 6) = iter
 iz( 7) = 0
 iz( 8) = iterm
 iz( 9) = 0
 iz(10) = 0
 iz(11) = 0
 iz(12) = 0
 CALL WRITE (dmpfil,iz,12,1)
 CALL CLOSE (dmpfil,rew)
 RETURN
 200 nomovs = nomovs + 1
 GO TO 50
 220 no     =-8
 ifile  = icrq
 CALL mesage (no,ifile,NAME)
 RETURN
END SUBROUTINE invpwr
