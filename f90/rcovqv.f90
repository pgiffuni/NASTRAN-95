SUBROUTINE rcovqv
     
!     THIS SUBROUTINE CALCULATES THE REACTION FORCES FOR THE REQUESTED
!     SUBSTRUCTURE
 
 LOGICAL :: reqf       ,reigen
 INTEGER :: fss        ,rfno       ,ua         ,rss       ,  &
     rc         ,tflag      ,signab     ,signc     ,  &
     prec       ,scrm       ,NAME(2)    ,scr1      ,  &
     scr2       ,scr4       ,scr5       ,scr6      ,  &
     scr7       ,scr8       ,mgg        ,kgg       ,  &
     bgg        ,buf1       ,buf2       ,buf3      ,  &
     buf4       ,sof1       ,sof2       ,sof3      ,  &
     bmtx       ,qvec       ,mmtx       ,kmtx      , sysbuf     ,pa         ,qa
!    9,               UVEC       ,SRD
 CHARACTER (LEN=27) :: swm
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm 
 COMMON /xmssg / ufm        ,uwm        ,uim        ,sfm        , swm
 COMMON /BLANK / dry        ,loop       ,step       ,fss(2)     ,  &
     rfno       ,neigv      ,lui        ,uinms(2,5) ,  &
     nosort     ,uthres     ,pthres     ,qthres
 COMMON /rcovcr/ icore      ,lcore      ,buf1       ,buf2       ,  &
     buf3       ,buf4       ,sof1       ,sof2       , sof3
 COMMON /rcovcm/ mrecvr     ,ua         ,pa         ,qa         ,  &
     iopt       ,rss(2)     ,energy     ,uimpro     ,  &
     range(2)   ,ireq       ,lreq       ,lbasic
 COMMON /system/ sysbuf     ,nout
 COMMON /mpyadx/ mcba(7)    ,mcbb(7)    ,mcbc(7)    ,mcbd(7)    ,  &
     mpyz       ,tflag      ,signab     ,signc      , prec       ,scrm
 COMMON /zzzzzz/ z(1)
 DATA    qvec  , kmtx,   mmtx,   bmtx,   k4mx                   /  &
     4HQVEC, 4HKMTX, 4HMMTX, 4HBMTX, 4HK4MX                 /
!     DATA    UVEC  , SRD   / 4HUVEC, 1                              /
 DATA    kgg   , mgg,    bgg,   k4gg,           NAME            /  &
     103   , 104,    109,    110,         4HRCOV,   4HQV    /
 DATA    scr1  , scr2,   scr4,   scr5,   scr6,   scr7,   scr8   /  &
     301   , 302,    304,    305,    306,    307,    308    /
 
 
!     CHECK TO SEE IF QVEC HAS ALREADY BEEN CALCULATED
 
 CALL mtrxi (scr4,rss,qvec,0,rc)
 IF (rc /= 1) GO TO 10
 qa = scr4
 RETURN
 
!     INITILIZE FOR QVEC CALCULATIONS
 
 10 prec  = 0
 tflag = 0
 scrm  = scr5
 signab= 1
 mpyz  = korsz(z(1)) - lreq
 reqf  = .false.
 IF (fss(1) == rss(1) .AND. fss(2) == rss(2)) reqf = .true.
 reigen = .false.
 malcom = 0
 IF (ua /= scr1) GO TO 30
 scr2 = 301
 scr1 = 302
 
!     CHECK THE DISPLACEMENT MATRIX
 
 30 mcbb(1) = ua
 CALL rdtrl (mcbb)
 IF (mcbb(1) <= 0) GO TO 9200
 
!     BRANCH ON RIGID FORMAT
 
 IF (rfno > 9) GO TO 9007
 SELECT CASE ( rfno )
   CASE (    1)
     GO TO 100
   CASE (    2)
     GO TO 100
   CASE (    3)
     GO TO 200
   CASE (    4)
     GO TO 9007
   CASE (    5)
     GO TO 9007
   CASE (    6)
     GO TO 9007
   CASE (    7)
     GO TO 9007
   CASE (    8)
     GO TO 400
   CASE (    9)
     GO TO 400
 END SELECT
 
!     STATIC SOUTION
 
!     Q = KU - P
 
!     SET UP LOAD VECTOR FOR SUBSTRACTION
 
 100 signc = -1
 mcbc(1) = pa
 IF (pa > 0) CALL rdtrl (mcbc)
 GO TO 500
 
!     NORMAL MODES
 
!     CHECK IF THE EIGEN VECTORS ARE COMPLEX
 
 200 IF (mcbb(5) >= 3) GO TO 300
 
!     REAL NORMAL MODES
 
!      Q = KU + MA   WHERE A = -(2*PI*FREQ)**2 * U
 
 reigen = .true.
 
!     MALCOM TAGG OF MDC, IN MSFC, RECOMMANDED THAT FOR RIGID FORMAT 3
!     THE SPC REACTION FORCE SHOULD NOT CONTAIN THE MASS TERM.  JULY/86
!     I.E.    Q = KU ONLY   (DROP THE MA TERM)
!     THUS,   GO TO 250 THEN TO 500
 
!     MARCH 1989 - MALCOM RECOMMENDATION REMOVED. IT CAUSES IMBALANCED
!     SPC FORCES
 
!     MALCOM = 1
 IF (malcom == 1) GO TO 250
 
!     CALCULATE THE ACCLERATION VECTOR FOR REAL NORMAL MODES
 
 in = ua
 CALL rcovva (in,0,0,0,0,scr8,rss,z(1),z(1),z(1))
 IF (in <= 0) GO TO 9200
 
 
!     INDICATE A POSITIVE SIGN ON THE M * A MULTIPLY
 
 250 signab  = 1
 mcbc(1) = 0
 IF (malcom == 1) GO TO 500
 GO TO 420
 
!     COMPLEX NORMAL MODES
 
!     Q = KU + BV + MA
 
!     CALCULATE THE COMPLEX VELOCITIES AND ACCLERATION VECTORS FOR
!     THE EIGENVECTORS
 
 300 in = ua
 
!     SEE MALCOM TAGG RECOMMENDATION, 25 LINES ABOVE
 
 CALL sofcls
 IF (malcom == 1) GO TO 445
 
 CALL rcovva (in,0,0,scr6,scr7,scr8,rss,z(1),z(1),z(1))
 IF (in <= 0) GO TO 9200
 
!     INDICATE ZERO LOAD VECTOR FOR NORMAL MODES
 
 mcbc(1) = 0
 GO TO 420
 
!     DYNAMIC ANALYSIS
 
!     Q = KU + BV + MA - P
 
 
!     SPLIT DISPLACEMENT, VELOCITIES AND ACCELERATIONS ONTO SEPERATE
!     FILES
 
 400 in = ua
 CALL rcovva (in,1,0,scr6,scr7,scr8,rss,z(1),z(1),z(1))
 IF (in <= 0) GO TO 9200
 
!     SETUP TO SUBTRACT LOAD VECTOR
 
 signc = -1
 mcbc(1) = pa
 IF (pa > 0) CALL rdtrl (mcbc)
 
!     COMMON PROCESSING FOR DYNAMICS AND NORMAL MODES
 
 
!     MULTIPLY AND ADD    SCR1 = MA - P
 
 420 CONTINUE
 IF (.NOT.reqf) GO TO 430
 mcba(1) = mgg
 CALL rdtrl (mcba)
 IF (mcba(1) > 0) GO TO 440
 430 CALL mtrxi (scr4,rss,mmtx,0,rc)
 IF (rc /= 1) GO TO 460
 mcba(1) = scr4
 CALL rdtrl (mcba)
 440 mcbb(1) = scr8
 CALL rdtrl (mcbb)
 CALL makmcb (mcbd,scr1,mcbb(3),mcbb(4),mcbb(5))
 CALL sofcls
 
 CALL mpyad (z(1),z(1),z(1))
 
 445 DO  i = 1,7
   mcbc(i) = mcbd(i)
 END DO
 signc = 1
 CALL sofopn (z(sof1),z(sof2),z(sof3))
 
!     MULTIPLY AND ADD   SCR8 = K4V + MCBC
 
 460 IF (reigen .OR. rfno == 9) GO TO 464
 IF (.NOT.reqf) GO TO 461
 mcba(1) = k4gg
 CALL rdtrl (mcba)
 IF (mcba(1) > 0) GO TO 462
 461 CALL mtrxi (scr4,rss,k4mx,0,rc)
 IF (rc /= 1) GO TO 464
 mcba(1) = scr4
 CALL rdtrl (mcba)
 462 mcbb(1) = scr7
 CALL rdtrl (mcbb)
 CALL makmcb (mcbd,scr8,mcbb(3),mcbb(4),mcbb(5))
 signab = 1
 CALL sofcls
 
 CALL mpyad (z(1),z(1),z(1))
 
 DO  i = 1,7
   mcbc(i) = mcbd(i)
 END DO
 signc = 1
 CALL sofopn (z(sof1),z(sof2),z(sof3))
 464 CONTINUE
 IF (reigen) GO TO 500
 
!     MULTIPLY AND ADD   SCR1 = BV + MCBC
 
 IF (.NOT.reqf) GO TO 470
 mcba(1) = bgg
 CALL rdtrl (mcba)
 IF (mcba(1) > 0) GO TO 480
 470 CALL mtrxi (scr4,rss,bmtx,0,rc)
 IF (rc /= 1) GO TO 500
 mcba(1) = scr4
 CALL rdtrl (mcba)
 480 mcbb(1) = scr7
 CALL rdtrl (mcbb)
 CALL makmcb (mcbd,scr1,mcbb(3),mcbb(4),mcbb(5))
 signab = 1
 CALL sofcls
 
 CALL mpyad (z(1),z(1),z(1))
 
 DO  i = 1,7
   mcbc(i) = mcbd(i)
 END DO
 signc = 1
 CALL sofopn (z(sof1),z(sof2),z(sof3))
 
!     COMMON PROCESSING FOR ALL RIGID FORMATS
 
 
!     MULTIPLY AND ADD  Q = KU + MCBC
 
 500 IF (.NOT.reqf) GO TO 520
 mcba(1) = kgg
 CALL rdtrl (mcba)
 IF (mcba(1) > 0) GO TO 540
 520 item = kmtx
 FILE = scr7
 CALL mtrxi (scr7,rss,kmtx,0,rc)
 IF (rc /= 1) GO TO 6000
 mcba(1) = scr7
 CALL rdtrl (mcba)
 540 mcbb(1) = scr6
 IF (reigen .OR. rfno <= 2) mcbb(1) = ua
 CALL rdtrl (mcbb)
 CALL makmcb (mcbd,scr4,mcbb(3),mcbb(4),mcbb(5))
 signab = 1
 CALL sofcls
 
 CALL mpyad (z(1),z(1),z(1))
 
 CALL wrttrl (mcbd)
 CALL sofopn (z(sof1),z(sof2),z(sof3))
 
!     COPY REACTIONS TO SOF
 
 CALL mtrxo (scr4,rss,qvec,0,rc)
 qa = scr4
 
 RETURN
 
!     ERRORS
 
 6000 IF (rc == 6) GO TO 9100
 CALL smsg (rc-2,item,rss)
 GO TO 9200
 9007 n = 7
 9100 CALL mesage (n,0,NAME)
 9200 qa = 0
 WRITE  (nout,6318) swm
 6318 FORMAT (a27,' 6318, OUTPUT REQUEST FOR REACTIONS FORCES IGNORED.')
 RETURN
 
END SUBROUTINE rcovqv
