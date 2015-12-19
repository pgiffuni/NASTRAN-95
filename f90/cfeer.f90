SUBROUTINE cfeer (eed,method,nfound)
     
!     PREVIOUS THIS ROUITNE IS CALLED CFCNTL
 
!     GIVEN REAL OR COMPLEX MATRICES, CFEER WILL SOLVE FOR THE
!     REQUESTED NUMBER OF EIGENVALUES AND EIGENVECTORS CLOSEST TO A
!     SPECIFIED POINT IN THE COMPLEX PLANE, FOR UP TO TEN POINTS,
!     VIA THE TRIDIAGONAL REDUCTION (FEER) METHOD.
!     THE SUBROUTINE NAME  CFEER  STANDS FOR COMPLEX FEER CONTROL.
 
!     DEFINITION OF INPUT AND OUTPUT PARAMETERS
 
!     IK(7)    = MATRIX CONTROL BLOCK FOR THE INPUT STIFFNESS MATRIX K
!     IM(7)    = MATRIX CONTROL BLOCK FOR THE INPUT MASS      MATRIX M
!     IB(7)    = MATRIX CONTROL BLOCK FOR THE INPUT DAMPING   MATRIX B
!     ILAM(7)  = MATRIX CONTROL BLOCK FOR THE OUTPUT EIGENVALUES
!     IPHI(7)  = MATRIX CONTROL BLOCK FOR THE OUTPUT EIGENVECTORS
!     IDMPFL   = FILE CONTAINING THE EIGENVALUE SUMMARY
!     ISCR(11) = SCRATCH FILES USED INTERNALLY
!     REG(1,I) = INPUT REAL      PART OF CENTER I (LAMBDA)
!     REG(2,I) = INPUT IMAGINARY PART OF CENTER I (LAMBDA)
!     REG(5,I) = PROBLEM SIZE MAXIMUM FOR SETTING QPR
!     REG(6,I) = SUPPRESSES ANY SPECIAL SYMMETRY LOGIC
!     REG(7,I) = NUMBER OF DESIRED ROOTS AROUND CENTER I
!     REG(8,1) = CONVERGENCE CRITERION (EQUIV. TO REG(1,2) TEMPORARILY)
 
 
 INTEGER, INTENT(IN)                      :: eed
 INTEGER, INTENT(IN OUT)                  :: method
 INTEGER, INTENT(OUT)                     :: nfound
 LOGICAL :: no_b     ,symmet   ,qpr
 INTEGER :: NAME(2)  ,iz(1)     , eigc(2)  ,want(10) ,have(10)
 DOUBLE PRECISION :: lambda   ,eps
 DIMENSION         ireg(7,1),ihead(10)
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON  /xmssg /  ufm      ,uwm      ,uim
 COMMON  /feeraa/  ik(7)    ,im(7)    ,ib(7)    ,ilam(7)   ,  &
     iphi(7)  ,idmpfl   ,iscr(11) ,reg(7,10) ,  &
     mcblt(7) ,mcbut(7) ,mcbvec(7),mcblmb(7)
 COMMON  /feerxc/  lambda(2),symmet   ,mreduc   ,nord      ,  &
     idiag    ,eps      ,northo   ,nord2     ,  &
     nord4    ,nordp1   ,nswp     ,jskip     ,  &
     no_b     ,it       ,ten2mt   ,tenmht    ,  &
     nstart   ,qpr      ,jreg     ,noreg     ,  &
     nzero    ,tenmtt   ,minopn   ,numort    , numran
 COMMON  /zzzzzz/  z(1)
 COMMON  /names /  rd       ,rdrew    ,wrt      ,wrtrew    ,  &
     rew      ,norew    ,eofnrw   ,rsp       , rdp
 COMMON  /system/  ksystm(65)
 COMMON  /output/  head(1)
 EQUIVALENCE       (ireg(1,1),reg(1,1)),(anodes,nodes)     ,  &
     (ksystm(2),nout)    ,(ksystm(40),nbpw)  , (asym,nonsym)       ,(z(1),iz(1))
 DATA    eigc   /  207,2/
 DATA    ihead  /  0,1009,2,7*0/
 DATA    NAME   /  4HCFCN,4HTL  /
 
!     FILE ALLOCATION
 
!     ISCR( 1)  CONTAINS  (LAMBDA**2*M + LAMBDA*B + K) = DYNAMIC MATRIX
!     ISCR( 2)  CONTAINS -(LAMBDA*M + B) = NOT REQUIRED WHEN B = 0
!     ISCR( 3)  CONTAINS LOWER TRIANGLE OF DECOMPOSED DYNAMIC MATRIX
!     ISCR( 4)  CONTAINS UPPER TRIANGLE OF DECOMPOSED DYNAMIC MATRIX
!     ISCR( 5)  CONTAINS REDUCED TRIDIAGONAL MATRIX ELEMENTS
!     ISCR( 6)  CONTAINS SPECIAL UPPER TRIANGLE FOR TRANSPOSED SWEEP
!     ISCR( 7)  CONTAINS THE ORTHOGONAL VECTORS
!     ISCR( 8)  CONTAINS OUTPUT EIGENVALUES , FOR INPUT TO CEAD1A
!     ISCR( 9)  CONTAINS OUTPUT EIGENVECTORS, FOR INPUT TO CEAD1A
!     ISCR(10)  SCRATCH FILE USED IN CFEER4
!     ISCR(11)  NOT USED
 
!     DEFINITION OF INTERNAL PARAMETERS
 
!     NODES  = NUMBER OF DESIRED ROOTS IN CURRENT NEIGHBORHOOD
!     EPS    = ACCURACY CRITERION - USED FOR REJECTING EIGENSOLUTIONS
!     NOREG  = TOTAL NUMBER OF CENTERS (NEIGHBORHOODS) INPUT,
!              EQUIVALENT TO THE NUMBER OF EIGC CONTINUATION CARDS
!     JREG   = COUNTER FOR CURRENT NEIGHBORHOOD
!     MREDUC = SIZE OF THE REDUCED PROBLEM IN CURRENT NEIGHBORHOOD
!     NFOUND = ACCUMULATED NUMBER OF ACCEPTABLE EIGENSOLUTIONS
!     NORD   = 2*N IF B.NE.0 AND = N IF B.EQ.0, WHERE B IS THE
!              DAMPING MATRIX AND N IS THE PROBLEM SIZE
!     NORD2  = VECTOR SIZE OF ORIGINAL PROBLEM (COMPLEX SINGLE
!              PRECISION OR COMPLEX DOUBLE PRECISION)
!     NSWP   = COMPLEX VECTOR SIZE FOR SWEEP ALGORITHM
!     NO B   = LOGICAL INDICATOR FOR ABSENCE OF DAMPING MATRIX B
!     SYMMET = LOGICAL INDICATOR FOR SYMMETRIC DYNAMIC MATRIX
!     NONSYM = PROGRAM INPUT WHICH FORCES THE PROGRAM TO CONSIDER
!              THE DYNAMIC MATRIX AS NON-SYMMETRIC
!     IT     = NUMBER OF DECIMAL DIGITS OF ACCURACY FOR THE COMPUTER
!     TEN2MT = 10**(2-T) CONVERGENCE CRITERION
!     TENMHT = 10**(-HALF*T) CONVERGENCE CRITERION
!     TENMTT = 10**(-THIRD*T) RIGID BODY ROOT CRITERION
!     NORTHO = TOTAL CURRENT NUMBER OF ORTHOGONAL VECTOR PAIRS ON
!              ORTHOGONAL VECTOR FILE. INITIALIZED TO NUMBER OF
!              EIGENVECTOR PAIRS ON THE RESTART FILE.
!     MINOPN = MINIMUM OPEN CORE NOT USED (WORDS)
!     NSTART = NUMBER OF INITIAL REORTHOGONALIZATION ATTEMPTS
!     IDIAG  = DIAG 12 PRINT CONTROL
!     QPR    = LOGICAL INDICATOR FOR VERY DETAILED PRINTOUT
!     WANT   = ARRAY OF DESIRED NUMBER OF ROOTS IN EACH NEIGHBORHOOD
!     HAVE   = ARRAY OF ACTUAL  NUMBER OF ROOTS IN EACH NEIGHBORHOOD
 
 northo = 0
 nfound = northo
 nzero  = northo
 jskip  = 0
 CALL sswtch (12,idiag)
 
!     TEST COMPUTING MACHINE TYPE AND SET PRECISION PARAMETERS
 
 IF (nbpw >= 60) GO TO 20
 it = 8*ksystm(55)
 GO TO 21
 20 it = 14*ksystm(55)
 21 ten2mt = 10.**(2-it)
 tenmht = 10.**(-it/2)
 tenmtt = 10.**(-it/3)
 ik(1)  = 101
 CALL rdtrl (ik)
 im(1)  = 103
 CALL rdtrl (im)
 ib(1)  = 102
 CALL rdtrl (ib)
 IF (ib(1) < 0 .OR. ib(6) == 0) ib(1) = 0
 
!     DETERMINE IF THE DYNAMIC MATRIX IS SYMMETRIC
 
 symmet = .false.
 IF (ik(1) /= 0 .AND. ik(4) /= 6) GO TO 30
 IF (im(1) /= 0 .AND. im(4) /= 6) GO TO 30
 IF (ib(1) /= 0 .AND. ib(4) /= 6) GO TO 30
 symmet = .true.
 30 DO  i = 1,11
   iscr(i)= 300+i
 END DO
 idmpfl = 203
 nz     = korsz(z)
 ibuf   = nz - ksystm(1) - 2
 limsum = 12
 iopn   = ibuf - limsum
 IF (idiag /= 0) WRITE (nout,600) iopn
 IF (iopn  <= 0) CALL mesage (-8,0,NAME)
 minopn = iopn
 ilam(1)= 308
 iphi(1)= 309
 ifile  = ilam(1)
 CALL OPEN (*500,ilam,z(ibuf),wrtrew)
 CALL CLOSE (ilam,rew)
 ifile  = iphi(1)
 CALL OPEN (*500,iphi,z(ibuf),wrtrew)
 CALL CLOSE (iphi,rew)
 CALL gopen (idmpfl,z(ibuf),wrtrew)
 CALL CLOSE (idmpfl,eofnrw)
 
!     PROCURE DATA FROM MAIN EIGC CARD
 
 ifile = eed
 CALL preloc (*500,z(ibuf),eed)
 CALL locate (*500,z(ibuf),eigc(1),flag)
 50 CALL fread (eed,ireg,10,0)
 IF (ireg(1,1) == method) GO TO 70
 60 CALL fread (eed,ireg,7,0)
 IF (ireg(6,1) /= -1) GO TO 60
 GO TO 50
 70 jreg  = 1
 eps   =.1D0/ik(2)/100.d0
 IF (reg(1,2) > 0.) eps = DBLE(reg(1,2))/100.d0
 unidum= SNGL(eps)*100.
 IF (idiag /= 0) WRITE (nout,75) unidum,reg(1,2)
 75 FORMAT (1H0,5HCFEER,6X,18HACCURACY criterion,1P,e16.8,  &
     8X,12H(INPUT value,e16.8,1H))
 
!     PROCURE DATA FROM EIGC CONTINUATION CARDS
 
 80 CALL fread (eed,ireg(1,jreg),7,0)
 IF (ireg(6,jreg) == -1) GO TO 90
 jreg  = jreg + 1
 IF (jreg > 10) GO TO 90
 GO TO 80
 90 CALL CLOSE (eed,rew)
 noreg  = jreg - 1
 nodcmp = 0
 numort = 0
 numran = 0
 jreg   = 0
 
!     PICK UP PARAMETERS FOR NEIGHBORHOOD I
 
 100 jreg = jreg + 1
 IF (jreg <= noreg) GO TO 105
 jreg = noreg
 IF (nzero > 0) jskip = -1
 GO TO 175
 105 x1 = reg(1,jreg)
 y1 = reg(2,jreg)
 anodes = reg(7,jreg)
 asym   = reg(6,jreg)
 IF (nonsym /= 0) symmet = .false.
 nprint = IFIX(reg(5,jreg))
 qpr = .false.
 IF (idiag /= 0 .AND. nprint >= ik(2)) qpr = .true.
 IF (idiag /= 0) WRITE (nout,110) jreg,x1,y1,nodes,nonsym
 110 FORMAT (1H0,5HCFEER,6X,12HNEIGHBORHOOD,i3,8X,8HCENTER =,2F18.8,  &
     8X,15HNO. des. rts. =,i5,8X,8HNONSYM =,i2/1H )
 
!     TEST IF USER PICKED THE ORIGIN
 
 IF (x1 /= 0. .OR. y1 /= 0.) GO TO 120
 x1 = x1 + .001
 WRITE (nout,601) uwm
 120 IF (nodes > 0) GO TO 130
 WRITE (nout,602) uwm,nodes
 nodes = 1
 130 want(jreg) = nodes
 have(jreg) = 0
 nord  = 2*ik(2)
 no_b  = .false.
 IF (ib(1) > 0) GO TO 140
 no_b  = .true.
 nord  = ik(2)
 140 nswp  = ik(2)
 nord2 = 2*nord
 nord4 = 2*nord2
 nordp1= nord + 1
 mreduc= 2*nodes + 10
 nomnf = nord - nfound
 IF (mreduc > nomnf) mreduc = nomnf
 lambda(1) = x1
 lambda(2) = y1
 IF (nodes > nord) WRITE (nout,606) uwm,nodes,jreg,noreg,lambda, nord
 ising = 0
 
!      FORM (LAMBDA**2*M + LAMBDA*B + K) = THE DYNAMIC MATRIX
 
 150 CALL cfeer1
 
!     CALL IN CDCOMP TO DECOMPOSE THE DYNAMIC MATRIX
 
 nodcmp = nodcmp + 1
 CALL cfeer2 (iret)
 IF (iret /= 0) GO TO 160
 GO TO 170
 160 iret = iret + ising
 WRITE (nout,603) uwm,iret,lambda
 IF (ising == 1) GO TO 100
 
!     SINGULAR MATRIX. INCREMENT LAMBDA AND TRY ONCE MORE.
 
 ising = 1
 lambda(1) = lambda(1) + .02D0
 lambda(2) = lambda(2) + .02D0
 GO TO 150
 
!     CALL IN DRIVER TO GENERATE REDUCED TRIDIAGONAL MATRIX
 
 170 CALL cfeer3
 IF (nstart > 2) GO TO 100
 
!     OBTAIN EIGENVALUES AND EIGENVECTORS
 
 CALL cfeer4
 have(jreg) = mreduc
 IF (mreduc <= nodes) GO TO 180
 i = mreduc - nodes
 WRITE (nout,607) uim,i,nodes,jreg,noreg,lambda
 180 nfound = nfound + mreduc
 IF (jreg < noreg .AND. nfound < nord) GO TO 100
 
!     FEER IS FINISHED. PERFORM WRAP-UP OPERATIONS.
 
 175 IF (jskip  < 0) CALL cfeer4
 IF (nfound == 0) GO TO 250
 IF (nfound >= nord) GO TO 220
 200 DO  i = 1,jreg
   IF (have(i) < want(i)) GO TO 240
 END DO
 GO TO 230
 
!     ALL SOLUTIONS FOUND
 
 220 WRITE (nout,604) uim
 IF (jreg < noreg) GO TO 240
 GO TO 200
 
!     EACH REQUESTED NEIGHBORHOOD HAS THE DESIRED NUMBER OF ROOTS
 
 230 iterm = 0
 GO TO 260
 
!     AT LEAST ONE REQUESTED NEIGHBORHOOD FAILS TO HAVE THE DESIRED
!     NUMBER OF ROOTS
 
 240 iterm = 1
 GO TO 260
 
!     ABNORMAL TERMINATION. NO ROOTS FOUND.
 
 250 iterm = 2
 
!     WRITE INFORMATION ON NASTRAN SUMMARY FILE
 
 260 ifile = idmpfl
 CALL OPEN (*500,idmpfl,z(ibuf),wrt)
 DO  i = 1,limsum
   iz(i) = 0
 END DO
 i = 0
 iz(i+2) = northo
 iz(i+3) = numran
 iz(i+5) = nodcmp
 iz(i+6) = numort
 iz(i+7) = iterm
 iz(i+8) = 1
 i = 2
 CALL WRITE (idmpfl,ihead(1),10,0)
 CALL WRITE (idmpfl,iz(i),40,0)
 CALL WRITE (idmpfl,head(1),96,1)
 CALL WRITE (idmpfl,iz(1),0,1)
 CALL CLOSE (idmpfl,eofnrw)
 
!     WRITE DUMMY TRAILER
 
 ixx   = ik(1)
 ik(1) = idmpfl
 CALL wrttrl (ik(1))
 ik(1) = ixx
 
!     INFORM USER IF RUN REGION SIZE CAN BE REDUCED
 
 IF (nbpw-36 < 0) THEN
   GO TO   300
 ELSE IF (nbpw-36 == 0) THEN
   GO TO   310
 ELSE
   GO TO   320
 END IF
 300 i = 4
 GO TO 330
 310 i = 6
 GO TO 330
 320 i = 10
 IF (nbpw == 64) i = 8
 330 i = (i*minopn)/1000
 IF (i < 0) i = 0
 WRITE (nout,605) uim,minopn,i
 RETURN
 
 500 CALL mesage (-1,ifile,NAME)
 RETURN
 
 
 600 FORMAT (1H1,27X,'*****  F E E R  *****  (FAST EIGENVALUE',  &
     ' EXTRACTION ROUTINE)  *****',  ////,1H ,i10,' SINGLE ',  &
     'PRECISION WORDS OF OPEN CORE, NOT USED (SUBROUTINE ', 'CFEER)', //)
 601 FORMAT (a25,' 3149',//5X,'USER SPECIFIED NEIGHBORHOOD CENTERED AT'  &
     ,      ' ORIGIN NOT ALLOWED, CENTER SHIFTED TO THE RIGHT .001',//)
 602 FORMAT (a25,' 3150',//5X,'DESIRED NUMBER OF EIGENVALUES',i8,3X,  &
     'INVALID. SET = 1.',//)
 603 FORMAT (a25,' 3151',//5X,'DYNAMIC MATRIX IS SINGULAR (OCCURRENCE',  &
     i3,') IN NEIGHBORHOOD CENTERED AT ',1P,2D16.8,//)
 604 FORMAT (a29,' 3159',//5X,'ALL SOLUTIONS HAVE BEEN FOUND.',//)
 605 FORMAT (a29,' 3160',//5X,'MINIMUM OPEN CORE NOT USED BY FEER',i9,  &
     ' WORDS (',i9,'K BYTES).',//)
 606 FORMAT (a25,' 3161',//5X,'DESIRED NUMBER OF EIGENSOLUTIONS',i5,  &
     ' FOR NEIGHBORHOOD',i3,' OF',i3,' CENTERED AT ',1P,2D16.8,  &
     //5X,'EXCEEDS THE EXISTING NUMBER',i5,  &
     ', ALL EIGENSOLUTIONS WILL BE SOUGHT.',//)
 607 FORMAT (a29,' 3166',//1X,i5,' MORE ACCURATE EIGENSOLUTIONS THAN ',  &
     'THE',i5,' REQUESTED HAVE BEEN FOUND FOR NEIGHBORHOOD',i3,  &
     ' OF',i3, //5X,'CENTERED AT ',1P,2D16.8,  &
     '. USE DIAG 12 TO DETERMINE ERROR ESTIMATES.',//)
END SUBROUTINE cfeer
