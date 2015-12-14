SUBROUTINE rcovr3
     
!     THE RCOVR3 MODULE RECOVERS DATA FOR SUBSTRUCTURE PHASE 3.
 
!     DISPLACEMENTS AND REACTIONS ARE COPIED FROM THE SOF TO GINO FILES.
!     FOR NORMAL MODES, LAMA IS CREATED FROM THE SOLN ITEM.
!     FOR STATICS, THE LOADS AND ENFORCED DISPLACEMENTS ARE FACTORED
!     AND COMBINED TO CORRESPOND WITH THE PHASE 2 SOLUTION SUBCASES.
 
!     JANUARY 1974
 
 LOGICAL :: first
 INTEGER :: rfno     ,trl      ,sysbuf   ,titles   ,otypp    ,  &
     otypun   ,iz(10)   ,pg       ,ps       ,po       ,  &
     ys       ,uas      ,qas      ,pgs      ,pss      ,  &
     pos      ,yss      ,lama     ,ivec(4)  ,ovec(4)  ,  &
     soln     ,uvec     ,qvec     ,initm(3) ,outdb(3) ,  &
     subr(2)  ,srd      ,here     ,buf1     ,buf2     ,  &
     buf3     ,buf4     ,rc       ,fss(2)   ,FILE     , scr1     ,scr2     ,scr3
 DIMENSION       mcbtrl(7)
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm      ,uwm      ,uim      ,sfm
 COMMON /BLANK / rfno     ,NAME(2)  ,noue     ,trl(7)   ,here(6)  , ibuf(3)
 COMMON /system/ sysbuf   ,nout
 COMMON /names / rd       ,rdrew    ,wrt      ,wrtrew   ,rew      , norew
 COMMON /output/ titles(1)
 COMMON /packx / itypp    ,otypp    ,irowp    ,nrowp    ,incp
 COMMON /unpakx/ otypun   ,irowun   ,nrowun   ,incun
 COMMON /zzzzzz/ z(1)
 EQUIVALENCE     (pg  ,ivec(1) )    ,(pgs ,ovec(1) )    ,  &
     (ps  ,ivec(2) )    ,(pss ,ovec(2) )    ,  &
     (po  ,ivec(3) )    ,(pos ,ovec(3) )    ,  &
     (ys  ,ivec(4) )    ,(yss ,ovec(4) )    ,  &
     (soln,initm(1))    ,(lama,outdb(1))    ,  &
     (uvec,initm(2))    ,(uas ,outdb(2))    ,  &
     (qvec,initm(3))    ,(qas ,outdb(3))    , (z(1),iz(1))
 DATA     pg   , ps  ,po  ,ys  ,uas ,qas ,pgs ,pss ,pos ,yss ,lama/  &
     101  , 102 ,103 ,104 ,201 ,202 ,203 ,204 ,205 ,206 ,207 /
 DATA     scr1 , scr2,scr3     / 301  , 302 ,303      /
 DATA     soln , uvec  , qvec   ,iblank  / 4HSOLN, 4HUVEC, 4HQVEC ,4H      /
 DATA     subr         , srd    / 4HRCOV, 4HR3  , 1      /
 
!     INITIALIZATION
 
 lcore = korsz(z)
 buf1  = lcore - sysbuf + 1
 buf2  = buf1  - sysbuf - 1
 buf3  = buf2  - sysbuf
 buf4  = buf3  - sysbuf
 lcore = buf4  - 1
 IF (lcore <= 0) CALL mesage (-8,0,subr)
 nogo  = 0
 itypp = 1
 otypp = 1
 irowp = 1
 incp  = 1
 otypun= 1
 irowun= 1
 incun = 1
 first = .false.
 CALL sofopn (z(buf1),z(buf2),z(buf3))
 DO  i = 1,6
   here(i) = 0
 END DO
 
!     CHECK DATA
 
!     NO EXTRA POINTS
 
 IF (noue /= -1) GO TO 6372
 
!     SUBSTRUCTURE NAME
 
 CALL fdsub (NAME,rc)
 IF (rc == -1) CALL smsg (-2,iblank,NAME)
 
!     PAIRS OF INPUT ITEMS AND OUTPUT BLOCKS
 
 CALL sfetch (NAME,soln,srd,rc)
 IF (rc /= 1) CALL smsg (2-rc,soln,NAME)
 IF (rfno == 1 .OR. rfno == 2) GO TO 15
 trl(1) = outdb(1)
 CALL rdtrl (trl)
 IF (trl(1) > 0) GO TO 15
 CALL mesage (1,outdb(1),subr)
 nogo = 1
 15 DO  i = 2,3
   IF (i == 1 .AND. (rfno == 1 .OR. rfno == 2)) CYCLE
   CALL softrl (NAME,initm(i),mcbtrl)
   rc = mcbtrl(1)
   IF (rc /= 1) CYCLE
   trl(1) = outdb(i)
   CALL rdtrl (trl)
   IF (trl(1) > 0) GO TO 20
   CALL mesage (1,outdb(i),subr)
   nogo = 1
   CYCLE
   20 here(i-1) = 1
 END DO
 
!     PAIRS OF DATA BLOCKS
 
 IF (rfno == 3 .OR. rfno == 8) GO TO 60
 DO  i = 1,4
   trl(1) = ivec(i)
   CALL rdtrl (trl)
   IF (trl(1) < 0)  CYCLE
   IF (i == 4 .AND. trl(6) == 0) CYCLE
   trl(1) = ovec(i)
   CALL rdtrl (trl)
   IF (trl(1) > 0) GO TO 40
   CALL mesage (1,ovec(i),subr)
   nogo = 1
   CYCLE
   40 here(i+2) = 1
 END DO
 
!     TERMINATE IF THERE WERE ERRORS
 
 60 IF (nogo /= 0) GO TO 9037
 
!     COPY DISPLACEMENTS AND REACTIONS FROM SOF TO GINO FILES
 
 IF (here(1) == 1) CALL mtrxi (uas,NAME,uvec,z(buf4),rc)
 IF (here(2) == 1) CALL mtrxi (qas,NAME,qvec,z(buf4),rc)
 
!     BRANCH ON RIGID FORMAT NUMBER
 
 IF (rfno == 3) GO TO 140
 
!     RIGID FORMAT  1 -- STATIC
!     RIGID FORMAT  2 -- INERTIAL RELIEF
!     RIGID FORMAT  8 -- FREQUENCY RESPONSE
!     RIGID FORMAT  9 -- TRANSIENT RESPONSE
!     *************************************
 
!     FETCH SOLN ITEM AND PROCESS GROUP 0 DATA
 
 CALL sfetch (NAME,soln,srd,rc)
 IF (rc /= 1) CALL smsg (2-rc,soln,NAME)
 CALL suread (fss,2,n,rc)
 WRITE (nout,63210) uim,fss,NAME
 CALL suread (ibuf,3,n,rc)
 IF (ibuf(1) /= rfno) GO TO 6322
 IF (ibuf(2) /= 1) GO TO 6324
 nc = ibuf(3)
 
!     WRITE NULL REACTIONS MATRIX TO PREVENT ERROR 3007 IN UMERGE
 
 IF (here(2) == 1) GO TO 80
 nrowp = 1
 CALL makmcb (trl,qas,1,2,1)
 CALL gopen (qas,z(buf4),wrtrew)
 DO  i = 1,nc
   CALL pack (0,qas,trl)
 END DO
 CALL CLOSE (qas,rew)
 CALL wrttrl (trl)
 
!     COPY FREQUENCIES ONTO PPF OR TIME STEPS ONTO TOL
 
 80 IF (rfno < 8) GO TO 120
 j = 1
 CALL sjump (j)
 FILE = lama
 CALL OPEN (*9001,lama,z(buf4),wrtrew)
 CALL fname (lama,ibuf)
 CALL WRITE (lama,ibuf,2,0)
 90 CALL suread (z,lcore,n,rc)
 CALL WRITE (lama,z,n,0)
 IF (rc == 1) GO TO 90
 CALL WRITE (lama,0,0,1)
 
!     WRITE NULL DYNAMIC LOADS MATRIX ONTO PPF
 
 CALL makmcb (trl,lama,1,2,1)
 IF (rfno == 9) GO TO 110
 DO  i = 1,nc
   CALL pack (0,lama,trl)
 END DO
 110 CALL wrttrl (trl)
 CALL CLOSE (lama,rew)
 
!     FOR EACH SUBCASE READ FROM THE SOLN, FORM A COMBINED VECTOR FROM
!     THE VECTORS OF THE APPLIED LOADS OR ENFORCED DISPLACEMENTS DATA
!     BLOCKS
 
 120 lcore = buf3 - 1
 DO  i = 1,4
   IF (here(i+2) == 0) CYCLE
   CALL rcovsl (NAME,0,ivec(i),scr1,scr2,scr3,ovec(i),z,z,lcore, first,rfno)
   IF (ovec(i) /= 0) first= .true.
 END DO
 GO TO 5000
 
!     RIGID FORMAT  3 -- NORMAL MODES
!     *******************************
 
!     WRITE NULL REACTIONS MATRIX TO PREVENT ERROR 3007 IN UMERGE
 
 140 IF (here(2) == 1) GO TO 150
 nrowp = 1
 CALL makmcb (trl,qas,1,2,1)
 CALL gopen (qas,z(buf4),wrtrew)
 CALL pack (0,qas,trl)
 CALL CLOSE (qas,rew)
 CALL wrttrl (trl)
 
!     GENERATE OFP ID RECORD FOR LAMA
 
 150 IF (lcore < 146) GO TO 9008
 CALL gopen (lama,z(buf4),wrtrew)
 DO  i = 3,50
   iz(i) = 0
 END DO
 iz( 1) = 21
 iz( 2) = 6
 iz(10) = 7
 DO  i = 1,96
   iz(i+50) = titles(i)
 END DO
 CALL WRITE (lama,z,146,1)
 
!     GET SOLN ITEM AND CHECK GROUP 0 DATA
 
 CALL sfetch (NAME,soln,srd,rc)
 IF (rc /= 1) CALL smsg (2-rc,soln,NAME)
 CALL suread (fss,2,n,rc)
 WRITE (nout,63210) uim,fss,NAME
 CALL suread (ibuf,-1,n,rc)
 IF (ibuf(1) /= rfno) GO TO 6322
 neigv = ibuf(2)
 IF (neigv > 0) GO TO 180
 
!     NO EIGENVALUES.  WRITE ZERO TRAILER TO INDICATE LAMA IS PURGED
 
 CALL CLOSE (lama,rew)
 CALL makmcb (trl,lama,0,0,0)
 CALL wrttrl (trl)
 GO TO 6323
 
!     COPY SOLN GROUP 1 TO LAMA RECORD 2 AND WRITE NON-ZERO TRAILER
 
 180 CALL suread (z,lcore,n,rc)
 CALL WRITE (lama,z,n,0)
 IF (rc == 1) GO TO 180
 CALL WRITE (lama,0,0,1)
 CALL CLOSE (lama,rew)
 CALL makmcb (trl,lama,0,0,0)
 trl(2) = 1
 CALL wrttrl (trl)
 
!     NORMAL MODULE EXITS
 
 5000 CALL sofcls
 RETURN
 
!     ABNORMAL MODULE EXITS
 
 6372 WRITE (nout,63720) ufm
 GO TO 9061
 9001 n = -1
 GO TO 9200
 6322 WRITE (nout,63220) sfm,ibuf(1),rfno
 GO TO 9061
 6323 WRITE (nout,63230) uwm
 GO TO 9300
 6324 WRITE (nout,63240) ufm,NAME
 GO TO 9061
 9008 n = -8
 GO TO 9200
 9037 n = -37
 GO TO 9200
 9061 n = -61
 9200 CALL sofcls
 CALL mesage (n,FILE,subr)
 9300 CALL sofcls
 RETURN
 
!     FORMAT STATEMENTS FOR DIAGNOSTIC MESSAGES
 
 63210 FORMAT (a29,' 6321, SUBSTRUCTURE PHASE 3 RECOVER FOR FINAL SOLUT',  &
     'ION STRUCTURE ',2A4, /35X,' AND BASIC SUBSTRUCTURE ',2A4)
 63220 FORMAT (a25,' 6322, SOLN HAS INCORRECT RIGID FORMAT NUMBER.',/32X,  &
     'PHASE 2 RIGID FORMAT WAS',i3,' AND PHASE 3 IS',i3)
 63230 FORMAT (a25,' 6323, NO EIGENVALUES FOR THIS SOLUTION')
 63240 FORMAT (a23,' 6324, PHASE 3 RECOVER ATTEMPTED FOR NON-BASIC ',  &
     'SUBSTRUCTURE ',2A4)
 63720 FORMAT (a23,' 6372, NO EXTRA POINTS ALLOWED IN PHASE 3 ',  &
     'SUBSTRUCTURING.')
END SUBROUTINE rcovr3
