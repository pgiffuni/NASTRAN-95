SUBROUTINE cinvpr (eed,method,nfound)
     
!     GIVEN REAL OR COMPLEX MATRICIES, CINVPR WILL SOLVE FOR ALL OF
!     THE EIGENVALUES AND EIGENVECTORS WITHIN A SPECIFIED REGION
 
!     DEFINITION OF INPUT AND OUTPUT PARAMETERS
 
!     FILEK(7) = MATRIX CONTROL BLOCK FOR THE INPUT STIFFNESS MATRIX K
!     FILEM(7) = MATRIX CONTROL BLOCK FOR THE INPUT MASS MATRIX M
!     FILEB(7) = MATRIX CONTROL BLOCK FOR THE INPUT DAMPING MATRIX B
!     FILELM(7)= MATRIX CONTROL BLOCK FOR THE OUTPUT EIGENVALUES
!     FILEVC(7)= MATRIX CONTROL BLOCK FOR THE OUTPUT EIGENVECTORS
!     DMPFIL   = FILE CONTAINING THE EIGENVALUE SUMMARY
!     SR1FIL-  = SCRATCH FILES USED INTERNALLY
!     SR0FIL
!     EPS      = CONVERGENCE CRITERIA
!     NOREG    = NUMBER OF REGIONS INPUT
!     REG(1,I) = X1 FOR REGION I
!     REG(2,I) = Y1 FOR REGION I
!     REG(3,I) = X2 FOR REGION I
!     REG(4,I) = Y2 FOR REGION I
!     REG(5,I) = L1 FOR REGION I
!     REG(6,I) = NO. OF DESIRED  ROOTS FOR REGION I
!     REG(7,I) = NO. OF ESTIMATED ROOTS IN REGION I
 
 
 
 INTEGER, INTENT(IN OUT)                  :: eed
 INTEGER, INTENT(IN)                      :: method
 INTEGER, INTENT(OUT)                     :: nfound
 LOGICAL :: noleft
 INTEGER :: eigc(2)  ,phidli  , switch    ,scrfil   ,ihead(10),ireg(7,1)
 INTEGER :: NAME(2)   ,filelm   ,filevc   ,  &
     REAL :: ,rdp      ,typek    ,typem    ,  &
     typeb     ,comflg   ,iz(1)    ,dmpfil   ,  &
     timed     ,filek    ,t1       ,t2       , fileb     ,filem
 REAL :: l         ,l1       ,maxmod
 DOUBLE PRECISION :: lam1      ,dz(1)    ,mindia
 DOUBLE PRECISION :: lambda    ,lmbda    ,dtemp(2)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg /  ufm
 COMMON /cdcmpx/  idum(30)  ,mindia
 COMMON /cinvpx/  filek(7)  ,filem(7) ,fileb(7) ,filelm(7),  &
     filevc(7) ,dmpfil   ,scrfil(11),noreg   , eps       ,reg(7,10),phidli
 COMMON /names /  rd        ,rdrew    ,wrt       ,wrtrew  ,  &
     rew       ,norew    ,eofnrw    ,rsp     , rdp
 COMMON /system/  ksystm(65)
 COMMON /zzzzzz/  z(1)
 COMMON /output/  head(1)
 COMMON /cinvxx/  lambda(2) ,switch   ,comflg    ,lmbda(2),  &
     iter      ,timed    ,nochng    ,rzero   ,  &
     ind       ,ivect    ,kreg      ,REAL    ,  &
     left      ,northo   ,noroot    ,nzero   ,  &
     lam1(2)   ,maxmod   ,nodes     ,noest   ,  &
     istart    ,ind1     ,iterx     ,isym
 EQUIVALENCE      (ksystm(1),isys )   ,(ireg(1,1),reg(1,1))
 EQUIVALENCE      (filek(5) ,typek)   ,(filem(5),typem)   ,  &
     (fileb(5) ,typeb)   ,(iz(1),z(1))
 EQUIVALENCE      (anodes   ,nodes)   ,(anoest,noest)     ,  &
     (z(1)     ,dz(1))   ,(ksystm(2),nout )
 DATA     ihead/  0,1009,2,7*0 /
 DATA     eigc /  207,2        /
 DATA     NAME /  4HCINV,4HPR  /
 DATA     SIGN /  1.0          /
 
!     DEFINITION OF INTERNAL PARAMETERS
 
!     REAL     = 0 - ALL MATRICIES ARE REAL
!                1 - AT LEAST ONE MATRIX IS COMPLEX
!     NSHIFT   = NO. OF SHIFT POINTS IN A REGION
!     NODES    = NO. OF DESIRED ROOTS IN A REGION
!     NOEST    = NO. OF ESTIMATED ROOTS IN A REGION
!     SHIFT    = INDEX OF THE CURRENT SHIFT POINT
!     ISHIFT   = INDEX OF THE CURRENT SHIFT POINT
!     IMIN     = LOWEST INDEX OF THE COMPLETED SHIFT POINTS
!     IMAX     = HIGHEST INDEX OF COMPLETED SHIFT POINTS
 
!     FILE ALLOCATION
 
!     SR1FIL CONTAINS (LAMBDA**2*M + LAMBDA*B + K)
!     SR2FIL CONTAINS -(B+LAMBDA*M)
!     SR3FIL CONTAINS THE LOWER TRIANGLE OF THE DECOMPOSED DYNAMIC MTRX
!     SR4FIL CONTAINS THE UPPER TRIANGLE OF THE DECOMPOSED DYNAMIC MTRX
!     SR5FIL IS USED AS A SCRATCH FOR CDCOMP
!     SR6FIL IS USED AS A SCRATCH FOR CDCOMP
!     SR7FIL IS USED AS A SCRATCH FOR CDCOMP
!     SR8FIL CONTAINS THE LOWER TRIANGLE L
!     SR9FIL CONTAINS THE UPPER TRIANGLE U
!     SR0FIL CONTAINS THE LEFT EIGENVECTORS
!     S11FIL CONTAINS  -(B + LAMBDA*M)
 
!     DEFINITION OF INTERNAL PARAMETERS
 
!     IND      = AN INDEX FOR GENERATING VARIOUS STARTING VECTORS
!     ITER     = TOTAL NUMBER OF ITERATIONS
!     NODCMP   = TOTAL NUMBER OF DECOMPOSITIONS
!     NOSTRT   = NUMBER OF STARTING POINTS USED
!     NOMOVS   = NUMBER OF TIMES A STARTING POINT HAD TO BE MOVED
!     RZERO    = DISTANCE FROM THE STARTING POINT TO THE CORNER OF THE
!                PARALELOGRAM
!     NOCHNG   = COUNT OF THE NUMBER OF MOVES WHILE LOOKING FOR ONE ROO
!     COMFLG   = 0 -
!              = 1 -
!              = 2 -
!              = 3 -
!              = 4 -
!              = 5 -
!              = 6 -
!     SWITCH   =
!     IVECT    =
!     KREG     = 0-NO VECTORS FOUND IN SEARCH AREA YET
!                1- A VECTOR HAS BEEN FOUND IN THE SEARCH AREA
!     ISING    = SINGULARITY FLAG
!     ITERM    = REASON FOR TERMINATING
!              = 1 - 2SINGULARITIES IN A ROW
!              = 2 - 4 MOVES WHILE TRACKING ONE ROOT
!              = 3 - ALL REGIONS COMPLETED
!              = 4 - 3*NOEST FOUND
!              = 5 - ALL ROOTS FOUND
!              = 8 - 200 ITERATIONS WITH ONE MOVE WITHOUR CONVERGING
!     TIMED    = TIME TOO FORM AND DECOMPOSE THE DYNAMIC MATRIX
!     LEFT     = 1 - DECOMPOSE MATRIX FOR THE COMPUTATION OF THE LEFT
!                EIGENVECTORS
 
 CALL sswtch (12,idiag)
 ind1 = 0
 nz   = korsz(z)
 CALL klock (istart)
 ibuf = nz - isys - 2
 ifile= filelm(1)
 CALL OPEN  (*500,filelm,z(ibuf),wrtrew)
 CALL CLOSE (filelm,rew)
 ifile = filevc (1)
 CALL OPEN  (*500,filevc,z(ibuf),wrtrew)
 CALL CLOSE (filevc,rew)
 CALL gopen (dmpfil,z(ibuf),wrtrew)
 CALL CLOSE (dmpfil,eofnrw)
 ifile = scrfil(10)
 CALL OPEN (*500,ifile,z(ibuf),wrtrew)
 CALL CLOSE (ifile,rew)
 noleft = .false.
 iz(1)  = 204
 CALL rdtrl (iz)
 IF (iz(1) < 0) noleft = .true.
 northo = 0
 nrow   = 2*filek(3)
 nrow2  = 2*nrow
 isym   = 1
 IF (filek(1) /= 0 .AND. filek(4) /= 6) GO TO 2
 IF (filem(1) /= 0 .AND. filem(4) /= 6) GO TO 2
 IF (fileb(1) /= 0 .AND. fileb(4) /= 6) GO TO 2
 isym = 0
 2 CONTINUE
 
!     PICK UP REGION PARAMETERS
 
 CALL preloc (*500,z(ibuf),eed)
 CALL locate (*500,z(ibuf),eigc(1),flag)
 6 CALL fread (eed,ireg,10,0)
 IF (method == ireg(1,1) .OR. method == -1) GO TO 8
 7 CALL fread (eed,ireg,7,0)
 IF (ireg(6,1) /= -1) GO TO 7
 GO TO 6
 8 jreg = 1
 eps  = .0001
 IF (reg(1,2) /= 0.) eps = reg(1,2)
 11 CALL fread (eed,ireg(1,jreg),7,0)
 IF (ireg(6,jreg) == -1) GO TO 9
 jreg = jreg + 1
 IF (jreg > 10) GO TO 9
 GO TO 11
 9 CALL CLOSE (eed,rew)
 noreg = jreg - 1
 jreg  = 0
 
!     PICK UP PARAMETERS FOR REGION I
 
 5 jreg   = jreg + 1
 iter   = 0
 nodcmp = 0
 nostrt = 0
 nomovs = 0
 x1     = reg(1,jreg)
 y1     = reg(2,jreg)
 x2     = reg(3,jreg)
 y2     = reg(4,jreg)
 l      = reg(5,jreg)
 anoest = reg(6,jreg)
 anodes = reg(7,jreg)
 IF (nodes == 0) nodes = 3*noest
 nshift = SQRT((x1-x2)**2+(y1-y2)**2)/l + 1.
 l1     = l*.5
 noroot = 0
 
 
!     FIND SHIFT POINT CLOSEST TO THE ORIGIN
 
 r = SQRT((x1-x2)**2 + (y1-y2)**2)
 IF (r > 0.0) THEN
   GO TO    15
 END IF
 10 WRITE  (nout,12) ufm
 12 FORMAT (a23,' 2366, REGION IMPROPERLY DEFINED ON EIGC CARD.')
 CALL mesage (-61,0,0)
 15 CONTINUE
 d = (FLOAT(nshift)*l-r)/2.0
 xx = x1 + d*(x1-x2)/r
 x2 = x2 + d*(x2-x1)/r
 x1 = xx
 yy = y1 + d*(y1-y2)/r
 y2 = y2 + d*(y2-y1)/r
 y1 = yy
 IF (idiag == 0) GO TO 7000
 WRITE  (nout,1000) x1,y1,x2,y2,l1,nodes,noest,nshift
 1000 FORMAT (1H1,5F10.2,3I5)
 7000 CONTINUE
 deltx = (x1-x2)/FLOAT(nshift)
 delty = (y1-y2)/FLOAT(nshift)
 xx    = x2 + deltx/2.
 yy    = y2 + delty/2.
 range = xx**2  + yy**2
 n     = nshift - 1
 shift = 1.
 IF (deltx /= 0.) GO TO 20
 anum1 = l1
 anum2 = 0.
 GO TO 25
 20 slope = delty/deltx
 arg   = SQRT(1.+slope**2)
 anum1 = slope*l1/arg
 anum2 = l1/arg
 25 CONTINUE
 IF (n == 0) GO TO 40
 DO  i = 1,n
   xx    = xx + deltx
   yy    = yy + delty
   rang  = xx**2 + yy**2
   IF (rang >= range) GO TO 40
   range = rang
   shift = i + 1
 END DO
 
!     COMPUTE COORDINATES OF CORNERS OF THE REGION
 
 40 xl2   = x2 + anum1
 yl2   = y2 - anum2
 imin  = shift
 imax  = shift
 
!     FIND THE MAXIMUM MODULUS OF THE SEARCH REGION
 
 maxmod = xl2**2 + yl2**2
 xx     = x2 - anum1
 yy     = y2 + anum2
 maxmod = AMAX1(maxmod,xx**2+yy**2)
 xx     = x1 + anum1
 yy     = y1 - anum2
 maxmod = AMAX1(maxmod,xx**2+yy**2)
 xx     = x1 - anum1
 yy     = y1 + anum2
 maxmod = AMAX1(maxmod,xx**2+yy**2)
 
!     INITIALIZE
 
 ind    = 0
 left   = 0
 45 ishift = shift
 
!     EVALUATE THE VALUE OF LAMBDA IN THE CENTER OF THE CURRENT SEARCH
!     REGION
 
 lambda(1) = x2 + (shift-.5)*deltx
 lambda(2) = y2 + (shift-.5)*delty
 IF (lambda(2) == 0.0D0) lambda(2) = .01*delty
 
!     COMPUTE DISTANCE TO FARTHEST CORNER OF THE SQUARE SEARCH REGION
 
 xx = xl2 + shift*deltx
 yy = yl2 + shift*delty
 rzero = (lambda(1)-xx)**2 + (lambda(2)-yy)**2
 rzero = SQRT(rzero)*1.05
 IF (idiag == 0) GO TO 7001
 WRITE  (nout,1216)rzero
 1216 FORMAT (//,10H rzero =     ,f10.4)
 7001 CONTINUE
 nostrt = nostrt + 1
 comflg = 0
 61 lmbda(1) = lambda(1)
 lmbda(2) = lambda(2)
 nochng = 0
 switch = 0
 ivect  = 0
 kreg   = 0
 ind    = ind + 1
 IF (IABS (ind) == 13) ind = 1
 ising  = 0
 GO TO 100
 80 ising  = 0
 switch = 1
 100 IF (nochng >= 4) GO TO 220
 nochng = nochng + 1
 CALL klock (t1)
 
!     CALL IN ADD LINK TO FORM (LAMBDA**2*M + LAMBDA*B + K)
 
 CALL cinvp1
 
!     CALL IN CD COMP TO DECOMPOSE THE MATRIX
 
 IF (idiag == 0) GO TO 7002
 WRITE  (nout,1001) lambda
 1001 FORMAT (10H1LAMBDA =   ,2D15.5)
 7002 CONTINUE
 nodcmp = nodcmp + 1
 CALL cinvp2 (*110)
 CALL klock (t2)
 GO TO 120
 110 IF (ising == 1) GO TO 210
 
!     SINGULAR MATRIX. INCREMENT LAMBDA AND TRY ONCE MORE
 
 ising = 1
 lambda(1) = lambda(1) + .02*rzero
 lambda(2) = lambda(2) + .02*rzero
 GO TO 100
 
!     DETERMINE THE TIME REQUIRED TO FORM AND DECOMPOSE THE DYNAMIC
!     MATRIX
 
 120 timed = t2 - t1
 IF (timed == 0) timed = 1
 
!     CALL IN MAIN LINK TO ITERATE FOR EIGENVALUES
 
 121 CALL cinvp3
 IF (left   == 1) GO TO 130
 IF (comflg == 2) GO TO 125
 IF (comflg == 1) GO TO 80
 GO TO 140
 125 nomovs = nomovs + 1
 GO TO 61
 
!     CALL IN LINK TO COMPUTE THE LEFT EIGENVECTOR
 
 130 dtemp(1)  = lambda(1)
 dtemp(2)  = lambda(2)
 lambda(1) = lam1(1)
 lambda(2) = lam1(2)
 131 switch    = -1
 CALL cinvp1
 
!     DECOMPOSE THE DYNAMIC MATRIX AT THE EIGENVALUE TO OBTAIN THE LEFT
!     EIGENVECTOR BY THE DETERMINATE METHOD
 
 IF (idiag == 0) GO TO 132
 WRITE (nout,1001) lambda
 132 CALL cinvp2 (*138)
 
!     BUILD LOAD FOR FBS
 
 d1 = nrow/2
 d2 = northo
 DO  i = 1,nrow,2
   k  = (i+1)/2
   dz(i) = SIGN*mindia/(1.d0+(1.d0-FLOAT(k)/d1)*d2)
   dz(i+1) = 0.0D0
 END DO
 SIGN = -SIGN
 CALL cdifbs (dz(1),z(ibuf))
 lambda(1) = dtemp(1)
 lambda(2) = dtemp(2)
 switch = 0
 
!     NORMALIZE AND STORE THE LEFT EIGENVECTOR
 
 CALL cnorm1 (dz(1),filek(2))
 IF (idiag == 0) GO TO 135
 WRITE  (nout,134) (dz(i),i=1,nrow)
 134 FORMAT (///,15H left vector   ,//,(10D12.4))
 135 CONTINUE
 IF (noleft .OR. isym == 0) GO TO 136
 ifile = phidli
 CALL OPEN  (*500,ifile,z(ibuf),wrt)
 CALL WRITE (ifile,dz(1),nrow2,1)
 CALL CLOSE (ifile,norew)
 136 ifile = scrfil(10)
 CALL gopen (ifile,z(ibuf),rd)
 CALL bckrec (ifile)
 CALL fread (ifile,dz(nrow+2),nrow2,1)
 CALL bckrec (ifile)
 CALL CLOSE (ifile,norew)
 
!     COMPUTE REAL LEFT VECTOR SCALED
 
 CALL cx trn y (dz(1),dz(nrow+2),dtemp)
 CALL cdivid (dz(1),dz(1),dtemp,nrow)
 CALL OPEN (*500,ifile,z(ibuf),wrt)
 CALL WRITE (ifile,dz(1),nrow2,1)
 CALL CLOSE (ifile,rew)
 GO TO 121
 138 lambda(1) = 1.01*lambda(1)
 lambda(2) = 1.01*lambda(2)
 GO TO 131
 140 IF (comflg >= 3) GO TO 200
 IF (comflg == 0) GO TO 160
 IF (idiag  == 0) GO TO 150
 WRITE  (nout,145) noreg,jreg
 145 FORMAT (2I10)
 150 IF (noreg == jreg) RETURN
 GO TO 5
 
!     FIND NEXT SHIFT POINT WHICH IS CLOSEST TO THE ORIGIN
 
 160 IF (imin /= 1) GO TO 170
 IF (imax == nshift) GO TO 250
 165 shift = shift + 1.
 imax  = imax  + 1
 lambda(1) = lmbda(1) + deltx
 lambda(2) = lmbda(2) + delty
 GO TO 45
 170 IF (imax /= nshift) GO TO 180
 175 shift = shift - 1.
 imin  = imin  - 1
 lambda(1) = lmbda(1) - deltx
 lambda(2) = lmbda(2) - delty
 GO TO 45
 180 xx   = lmbda(1) - deltx
 yy   = lmbda(2) - delty
 rang = xx**2 + yy**2
 xx   = lmbda(1) + deltx
 yy   = lmbda(2) + delty
 range= xx**2 + yy**2
 IF (range-rang > 0.0) THEN
   GO TO   165
 ELSE
   GO TO   175
 END IF
 200 iterm = comflg
 GO TO 260
 
!     SINGULARITY ENCOUNTERED TWICE IN A ROW
 
 210 iterm = 1
 GO TO 260
 
!     4 MOVES WHILE TRACKING ONE ROOT
 
 220 iterm = 2
 GO TO 260
 
!     REGIONS COMPLETED
 
 250 iterm = 3
 
!     SET UP THE SUMMARY FILE
 
 260 ifile = dmpfil
 CALL OPEN  (*500,dmpfil,z(ibuf),wrt)
 CALL WRITE (dmpfil,ihead(1),10,0)
 i       = 0
 iz(i+2) = northo
 iz(i+3) = nostrt
 iz(i+4) = nomovs
 iz(i+5) = nodcmp
 iz(i+6) = iter
 iz(i+7) = iterm
 DO  i = 8,12
   iz(i) = 0
 END DO
 i = 2
 CALL WRITE (dmpfil,iz(i),40,0)
 CALL WRITE (dmpfil,head(1),96,1)
 CALL WRITE (dmpfil,iz(1),0,1)
 CALL CLOSE (dmpfil,eofnrw)
 
!     WRITE DUMMY TRAILER
 ixx = filek(1)
 filek(1) = dmpfil
 CALL wrttrl (filek(1))
 filek(1) = ixx
 nfound = northo
 IF (idiag == 0) GO TO 350
 j = 12
 WRITE  (nout,300)(iz(i),i=1,j)
 300 FORMAT (///,12I10)
 350 CONTINUE
 IF (iterm == 5) RETURN
 GO TO 150
 
 500 CALL mesage (-1,ifile,NAME)
 RETURN
END SUBROUTINE cinvpr
