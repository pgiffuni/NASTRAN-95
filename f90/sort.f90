SUBROUTINE sort (idum,jdum,nr,keywd,z,nwds)
     
!     THE ORIGINAL NASTRAN SORT ROUTINE FOR IN-CORE SORTING AND FILE
!     SORT IS NOW RENAMED SORTI
!     (ONLY 5 PERCENT OF NASTRAN ROUTINES ACTUALLY CALL SORTI, WITH NON-
!     ZERO IDUM AND JDUM)
 
!     THIS NEW SORT ROUTINE WITH IDUM=JDUM=0, PERFORMS ONLY IN-CORE SORT
!     FOR INTEGERS, FLOATING POINT NUMBERS, AND BCD WORDS, BY THE
!     MODIFIED SHELL METHOD
!     IT USES MUCH LESS CORE SPACE
 
!     ARRAY Z IS NR-ROWS BY (NWDS/NR)-COLUMNS IN SIZE
!     DATA STORED ROW-WISE IN Z, AND TO BE SORTED BY KEYWD-TH ROW
 
!     USE A NEGATIVE KEYWD  IF THE ORIGINAL ORDER OF THE TABLE ENTRIES
!     ARE TO BE PRESERVED AND THE COLUMN OF KEYWORDS CONTAINS DUPLICATES
!     (INTEGER SORT ONLY)    E.G.
 
!     ORIGINAL TABLE     SORTED(KEYWD=+1)       SORTED(KEYWD=-1)
!     ---------------    ----------------       ----------------
!       1      4             1      4               1      4
!       2      2             1     10               1      3
!       1      3             1      3               1     10
!       1     10             2      2               2      2
 
 
!     THIS ROUTINE WOULD SWITCH BACK TO THE OLD SHUTTLE EXCHANGE METHOD
!     NUMBERS OVERFLOW DUE TO THE REQUIREMENT THAT ORIGINAL ORDER MUST B
!     MAINTAINED
 
!     ENTRY POINTS
 
!     SORT   - TABLE SORT BY INTEGER
!     SORTF  - TABLE SORT BY F.P. NUMBER
!     SORTA  - TABLE SORT BY ALPHABETS, 4-BCD CHARACTERS
!     SORTA8 - TABLE SORT BY ALPHABETS, 8-BCD CHAR. (KEYWD AND KEYWD+1)
!     SORTA7 - SAME AS SORTA8, EXCEPT LEADING CHAR. IS IGNORED
!     SORT2K - 2-KEYWORD SORT, SORT BY KEYWD AND KEYWD+1, INTEGER OR
!              REAL NUMBER KEYS. NEGATIVE KEYWD IS IGNORED
 
!     THE TWO SORT CALLS OF THE FOLLOWING FORM CAN BE REPLACED BY ONE CA
!     TO SORT2K, WHICH IS FASTER, NO DANGER OF NUMBER OVERFLOW, AND THE
!     ORIGINAL SEQUENCE WILL NOT CHANGE WHEN THERE ARE DUPLICATES.
 
!         CALL SORT (0,0,N1,-(N2+1),TABLE,N3)
!         CALL SORT (0,0,N1,-N2,    TABLE,N3)
!              CAN BE REPLACED BY
!         CALL SORT2K (0,0,N1,N2,TABLE,N3)
 
 
!     WRITTEN BY G.CHAN/SPERRY, 3/1987
 
 
 INTEGER, INTENT(IN OUT)                  :: idum
 INTEGER, INTENT(IN)                      :: jdum
 INTEGER, INTENT(IN)                      :: nr
 INTEGER, INTENT(IN)                      :: keywd
 INTEGER, INTENT(IN OUT)                  :: z(nr,1)
 INTEGER, INTENT(IN)                      :: nwds
 LOGICAL :: rvsbcd
 INTEGER :: zi,      zn,      temp, two31,   two, subr(6)
 REAL :: ri,      rn
 COMMON /system/ ibuf,    nout,    dm37(37),nbpw
 COMMON /machin/ mach,    ijhlf(2),lqro
 COMMON /two   / two(16)
 EQUIVALENCE     (zi,ri), (zn,rn)
 DATA    subr  / 2H  ,    2HF ,    2HA ,    2HA8,    2HA7,    2H2K/
 
!     CHECK ERROR, CHECK DATA TYPE, AND PREPARE FOR SORT
 
 isort = 1
 GO TO 10
 
 ENTRY sortf (idum,jdum,nr,keywd,z,nwds)
!     =======================================
 isort = 2
 GO TO 10
 
 ENTRY sorta (idum,jdum,nr,keywd,z,nwds)
!     =======================================
 isort = 3
 GO TO 10
 
 ENTRY sorta8 (idum,jdum,nr,keywd,z,nwds)
!     ========================================
 isort = 4
 GO TO 10
 
 ENTRY sorta7 (idum,jdum,nr,keywd,z,nwds)
!     ========================================
 isort = 5
 GO TO 10
 
 ENTRY sort2k (idum,jdum,nr,keywd,z,nwds)
!     ========================================
 isort = 6
 
 10 IF (nwds == 0) GO TO 330
 IF (idum /= 0 .OR. jdum /= 0) GO TO 300
 rvsbcd = MOD(lqro,10) == 1
 key  = IABS(keywd)
 IF (key > nr) GO TO 280
 nc = nwds/nr
 IF (nc*nr /= nwds) GO TO 280
 m  = nc
 IF (isort /= 1 .OR. keywd >= 0) GO TO 30
 
!                     - INTEGER SORT ONLY -
!     IF ORIGINAL ORDER IS TO BE MAINTAINED WHERE DUPLICATE KEYWORDS MAY
!     OCCUR, ADD INDICES TO THE KEYWORDS (GOOD FOR BOTH POSITIVE AND
!     NEGATIVE RANGES, AND BE SURE THAT KEYWORDS ARE NOT OVERFLOWED),
!     SORT THE DATA, AND REMOVE THE INDICES LATER
 
!     IF KEYWORD OVERFLOWS, SWITCH TO SHUTTLE EXCHANGE METHOD
 
 IF (nc >= two(16) .AND. nbpw <= 32) GO TO 200
 j = 30
 IF (nbpw >= 60) j = 62
 two31 = 2**j
 limit = (two31-nc)/nc
 DO  i = 1,nc
   j = z(key,i)
   IF (IABS(j) > limit) GO TO 200
   j = j*nc + i
   k =-1
   IF (j < 0) k =-nc
   20 z(key,i) = j + k
 END DO
 
!     SORT BY
!     MODIFIED SHELL METHOD, A SUPER FAST SORTER
 
 30 m = m/2
 IF (m == 0) GO TO 180
 j = 1
 k = nc - m
 40 i = j
 50 n = i + m
 zi= z(key,i)
 zn= z(key,n)
 SELECT CASE ( isort )
   CASE (    1)
     GO TO 60
   CASE (    2)
     GO TO 80
   CASE (    3)
     GO TO 90
   CASE (    4)
     GO TO 90
   CASE (    5)
     GO TO 90
   CASE (    6)
     GO TO 60
 END SELECT
!           INT FP A4 A8 A7 2K
 
 60 IF (zi-zn < 0.0) THEN
   GO TO   170
 ELSE IF (zi-zn == 0.0) THEN
   GO TO    70
 ELSE
   GO TO   150
 END IF
 70 IF (isort == 1) GO TO 170
 IF (z(key+1,i)-z(key+1,n) > 0.0) THEN
   GO TO   150
 ELSE
   GO TO   170
 END IF
 80 IF (ri-rn > 0.0) THEN
   GO TO   150
 ELSE
   GO TO   170
 END IF
 90 kk = 1
 IF (isort == 5) GO TO 110
 
!     COMPARE 1ST BYTE, THEN COMPARE 2ND, 3RD, AND 4TH BYTES TOGETHER
!     IF MACHINE DOES NOT USE REVERSED BCD ORDER. THOSE MACHINES WITH
!     REVERSED BCD ORDER (VAX, ULTRIX, S/G) MUST COMPARE EACH BYTE
!     SEPERATELY BECAUSE OF THE SIGN BIT
 
 100 IF (khrfn1(zero,4,zi,1) - khrfn1(zero,4,zn,1) < 0) THEN
   GO TO   170
 ELSE IF (khrfn1(zero,4,zi,1) - khrfn1(zero,4,zn,1) == 0) THEN
   GO TO   110
 ELSE
   GO TO   150
 END IF
 110 IF (.NOT.rvsbcd) IF (khrfn1(zi,1,zero,4)-khrfn1(zn,1,zero,4))  &
     170,140,150
 IF (khrfn1(zero,4,zi,2) - khrfn1(zero,4,zn,2) < 0) THEN
   GO TO   170
 ELSE IF (khrfn1(zero,4,zi,2) - khrfn1(zero,4,zn,2) == 0) THEN
   GO TO   120
 ELSE
   GO TO   150
 END IF
 120 IF (khrfn1(zero,4,zi,3) - khrfn1(zero,4,zn,3) < 0) THEN
   GO TO   170
 ELSE IF (khrfn1(zero,4,zi,3) - khrfn1(zero,4,zn,3) == 0) THEN
   GO TO   130
 ELSE
   GO TO   150
 END IF
 130 IF (khrfn1(zero,4,zi,4) - khrfn1(zero,4,zn,4) < 0) THEN
   GO TO   170
 ELSE IF (khrfn1(zero,4,zi,4) - khrfn1(zero,4,zn,4) == 0) THEN
   GO TO   140
 ELSE
   GO TO   150
 END IF
 140 IF (isort <= 3 .OR. kk == 2) GO TO 170
 zi = z(key+1,i)
 zn = z(key+1,n)
 kk = 2
 GO TO 100
 150 DO  l = 1,nr
   temp = z(l,i)
   z(l,i) = z(l,n)
   z(l,n) = temp
 END DO
 i = i - m
 IF (i >= 1) GO TO 50
 170 j = j + 1
 IF (j-k > 0) THEN
   GO TO    30
 ELSE
   GO TO    40
 END IF
 180 IF (isort /= 1 .OR. keywd >= 0) GO TO 330
 DO  i = 1,nc
   z(key,i) = z(key,i)/nc
 END DO
 GO TO 330
 
!     SORT BY
!     SHUTTLE EXCHANGE THETHOD, A SLOW SORTER
!     (THIS WAS NASTRAN ORIGINAL SORTER, MODIFIED FOR 2-D ARRAY OPERATIO
!     WITH 20-COLUMN LIMITATION REMOVED)
 
 200 IF (i <= 1) GO TO 220
 j = i - 1
 DO  i = 1,j
   z(key,i) = z(key,i)/nc
 END DO
 
 220 DO  ii = 2,nc
   zi = z(key,ii)
   jj = ii - 1
   IF (zi >= z(key,jj)) CYCLE
   230 jj = jj - 1
   IF (jj > 0) IF (zi-z(key,jj)) 230,240,240
   240 jj = jj + 2
   DO  i = 1,nr
     temp = z(i,ii)
     m = ii
     DO  j = jj,ii
       z(i,m) = z(i,m-1)
       m = m - 1
     END DO
     z(i,jj-1) = temp
   END DO
 END DO
 GO TO 330
 
!     ERROR. FORCING A WALK BACK
 
 280 WRITE  (nout,290) subr(isort),nr,key,nwds,nc
 290 FORMAT ('0*** ERROR IN SORT',a2,4I8)
 GO TO 320
 300 WRITE  (nout,310)
 310 FORMAT ('0*** CALLING ROUTINE SHOULD CALL SORTI')
!WKBR  320 CALL ERRTRC ('SORT    ',320)
 320 CONTINUE
 330 RETURN
END SUBROUTINE sort
