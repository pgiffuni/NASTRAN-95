SUBROUTINE scheme (ig,inv,ii3,INT,icc,ild,norig,ip,un,z)
     
 
 INTEGER, INTENT(IN OUT)                  :: ig(1)
 INTEGER, INTENT(IN)                      :: inv(1)
 INTEGER, INTENT(IN)                      :: ii3
 INTEGER, INTENT(OUT)                     :: INT(1)
 INTEGER, INTENT(IN OUT)                  :: icc(1)
 INTEGER, INTENT(OUT)                     :: ild(1)
 INTEGER, INTENT(IN OUT)                  :: norig(1)
 INTEGER, INTENT(IN OUT)                  :: ip(1)
 REAL, INTENT(IN OUT)                     :: un(1)
 INTEGER, INTENT(OUT)                     :: z(1)
 INTEGER :: scr1,     rd,       rdrew,    wrt, wrtrew,   rew ,     sub(2)
 
 COMMON /banda / ibuf1,    dum4a(4), method,   icrit
 COMMON /bandb / nbitin,   kore,     ifl,      ngrid,    ipass, nw,       kdim
 COMMON /bandd / dum7d(7), neq,      neqr
 COMMON /bands / nn,       mm,       dum2(2),  maxgrd,   maxdeg,  &
     kmod,     mach,     mindeg,   nedge
 COMMON /geomx / gdum(3),  scr1
 COMMON /system/ ibuf,     nout
 COMMON /names / rd,       rdrew,    wrt,      wrtrew,   rew
 DATA            sub /     4HSCHE,4HME    /
 
!     THIS ROUTINE IS USED ONLY IN BANDIT MODULE
 
!     ZERO OUT CORE SPACE AND SET BANDWIDTH IMPROVEMENT FLAG, JUMP
!     JUMP = 1,  NO IMPROVEMENT OF CRITERION SELECTED
!          = 0,  IMPROVEMENT
 
 DO  i = 1,kore
   z(i) = 0
 END DO
 jump = 1
 
!     READ ELEMENT DATA FROM GEOM2 FILE AND SET UP CONNECTION TABLE IG.
!     ALSO, EXAMINE MPC EQUATIONS.
 
 CALL bread (ig,inv,ii3,norig,z)
 IF (ngrid <= 0) RETURN
 
!     NGRID = NO. OF GRID POINTS IN THE PROBLEM
!           =  0, ONE OR MORE SEQGP CARD IS PRESENT IN NASTRAN INPUT
!                 DECK, AND/OR QDSEP ELEMENTS
!           = -1, INSUFFICIENT CORE SPACE (IG TABLE TOO SMALL)
!           = -2, INSUFFICIENT SCRATCH AREA WHILE USING CTHMCK
!           = -3, INSUFFICIENT SCRATCH AREA WHILE USING GIBSTK
 
!     MODIFY IG TO ACCOUNT FOR MPC EQUATIONS AND RIGID ELEMENTS
 
 IF (neq+neqr /= 0) CALL tiger (ig,icc,inv,ii3,norig,z,un)
 
!     SORT ORIGINAL GRID NOS. AND OUTPUT THE LIST IN INT, WHERE INT(I)
!     IS THE I-TH ORIGINAL GRID NUMBER.
!     ALSO OUTPUT ILD, WHERE IDL(I) = SORTED INTERNAL NO. CORRESPONDING
!     TO THE UNSORTED BANDIT INTERNAL LABEL I.
 
!     CALL BRIGIT (INV,II3,INT,ILD)
!     BRIGIT AND INTERN ARE NOW REPLACED BY 17 LINES BELOW /G.CHAN 1988
 
 k = 0
 DO  i = 1,ii3
   IF (inv(i) == 0) CYCLE
   k = k + 1
   INT(k) = inv(i)
 END DO
 CALL sort (0,0,1,1,INT,nn)
 DO  i = 1,nn
   j = INT(i)
   IF (j <= 0) GO TO 120
   loc = j - 1
   16   loc = MOD(loc,kmod) + 1
   IF (inv(loc) == 0) GO TO 120
   IF (inv(loc) /= j) GO TO 16
   j = inv(loc+ii3)
   ild(j) = i
 END DO
 
!     METHOD WAS SET IN BANDIT -
!     METHOD = -1, CM ONLY,    = +1, GPS ONLY,    = 0, BOTH METHODS.
 
 IF (method /= 0) GO TO 20
 
!     SAVE ORIGINAL GRID POINT ORDERING (ILD) IN SCR1 FILE
 
 CALL OPEN (*70,scr1,z(ibuf1),wrtrew)
 CALL WRITE (scr1,ild,nn,1)
 CALL CLOSE (scr1,rew)
 
!     RE-SEQUENCE GRIDS WITH CUTHILL-MCKEE ALGORITHM
 
 20   i = maxgrd + 2
 j = i + maxgrd
 IF (maxdeg > maxgrd) j = j + maxdeg - maxgrd
 k = j + maxgrd
 CALL cthmck (80,1,2,icrit,ig,inv,inv(i),inv(j),inv(k),INT,icc,  &
     ild,ip,jump,un,z)
 ngrid1 = ngrid
 IF (method < 0) THEN
   GO TO    60
 ELSE IF (method == 0) THEN
   GO TO    25
 ELSE
   GO TO    30
 END IF
 
!     READ ORIGINAL SEQUENCE BACK IF CTHMCK MAKES NO IMPROVEMENT
 
 25   IF (jump == 0) GO TO 30
 CALL OPEN (*70,scr1,z(ibuf1),rdrew)
 CALL READ (*80,*80,scr1,ild,nn,1,m)
 CALL CLOSE (scr1,rew)
 30   DO  k1 = 1,nn
   INT(k1) = ild(k1)
 END DO
 
!     RESEQUENCE NODES WITH GPS ALGORITHM.
 
 k1 = 1
 k2 = k1 + kdim
 k3 = k2 + kdim
 k4 = k3 + kdim
 k5 = k4 + kdim/2
 CALL gibstk (ig,INT,ild,inv(i),inv,inv(j),inv(k),icc,jump,icrit,  &
     z(k1),z(k2),z(k3),z(k4),z(k5),un,kdim)
 
!     GENERATE SEQGP CARDS AND OUTPUT THEM TO GEOM1 FILE
 
 60   CALL bseqgp (norig,ild,jump)
 IF (ngrid1 == -2 .OR. ngrid == -3) GO TO 100
 RETURN
 
!     SCRATCH FILE ERROR
 
 70   k = -1
 GO TO 90
 80   k = -2
 90   CALL mesage (k,scr1,sub)
 
 100  WRITE  (nout,110) kdim
 110  FORMAT ('0*** BANDIT SCRATCH ARRAY OF,I5,20H WORDS IS TOO SMALL.' &
      ,/5X,'USER COULD USE ONE OF THE FOLLOWING OPTIONS AND RESUBMIT ', &
      'JOB. (USERS MANUAL P.2.1-1)', /5X, &
      'INCREASE SCRATCH ARRAY BY NASTRAN BANDTDIM OPTION, OR',/5X, &
      'SWITCH TO CUTHILL-MCKEE METHOD ONLY BY  BANDTMTH=1 OR',/5X, &
      'SKIP BANDIT COMPUTATION BY SETTING NASTRAN CARD BANDIT=-1',//)
 GO TO 140
 
 120  WRITE  (nout,130) k,nn,ii3,kmod,maxgrd,maxdeg
 130  FORMAT ('0*** BANDIT FATAL ERROR - TRY TO RERUN JOB WITH ',  &
     22H'NASTRAN BANDTDIM = N',' WHERE N = 3,4,...,OR 9', //5X,  &
     '@17/  K,NN,II3,KMOD,MAXGRD,MAXDEG =',6I8)
 140  CALL mesage (-37,sub,sub)
 RETURN
END SUBROUTINE scheme
