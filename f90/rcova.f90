SUBROUTINE rcova
     
!     RCOVA CREATES THE SOLN ITEM FOR A FINAL SOLUTION STRUCTURE (FSS)
!     IN PHASE 2 OF SUBSTRUCTURING
 
 LOGICAL :: mrecvr
 INTEGER :: iz(1)      ,NAME(2)    ,soln       ,dry        ,  &
     step       ,fss        ,rfno       ,buf1       ,  &
     buf2       ,buf3       ,sysbuf     ,rc         ,  &
     sof1       ,sof2       ,sof3       ,km(5)      ,  &
     kmu(5)     ,schk       ,uvec       ,phis       , scr1
 COMMON /BLANK / dry        ,loop       ,step       ,fss(2)     ,  &
     rfno       ,neigv      ,lui        ,uinms(2,5) ,  &
     nosort     ,uthres     ,pthres     ,qthres
 COMMON /rcovcr/ icore      ,lcore      ,buf1       ,buf2       ,  &
     buf3       ,buf4       ,sof1       ,sof2       , sof3
 COMMON /rcovcm/ mrecvr     ,ua         ,pa         ,qa         ,  &
     iopt       ,rss(2)     ,energy     ,uimpro     ,  &
     range(2)   ,ireq       ,lreq       ,lbasic
 COMMON /system/ sysbuf     ,nout
 COMMON /zzzzzz/ z(1)
 EQUIVALENCE     (z(1),iz(1))
 DATA    soln  , uvec,  phis / 4HSOLN,4HUVEC,4HPHIS /
 DATA    km    / 4HKMTX,4HMMTX,4HUVEC,4HBMTX,4HK4MX /
 DATA    kmu   / 103,104,106,109,110 /
 DATA    schk  / 3   /
 DATA    scr1  / 301 /
 DATA    NAME  / 4HRCOV,4HA    /
 
!     INITIALIZE
 
 sof1 = korsz(z) - lreq - sysbuf + 1
 sof2 = sof1 - sysbuf - 1
 sof3 = sof2 - sysbuf
 buf1 = sof3 - sysbuf
 buf2 = buf1 - sysbuf
 buf3 = buf2 - sysbuf
 icore= 1
 lcore= buf3 - 1
 IF (lcore <= 0) CALL mesage (-8,0,NAME)
 CALL sofopn (z(sof1),z(sof2),z(sof3))
 
!     COPY KGG, MGG, UVEC, BGG AND K4GG TO THE SOF IF THEY ARNT THERE
 
 DO  i = 1,5
   IF (km(i) == uvec .AND. mrecvr) CYCLE
   IF (dry < 0) GO TO 5
   3 CALL mtrxo (kmu(i),fss,km(i),z(buf1),rc)
   SELECT CASE ( rc )
     CASE (    1)
       GO TO 20
     CASE (    2)
       GO TO 15
     CASE (    3)
       GO TO 20
     CASE (    4)
       GO TO 10
     CASE (    5)
       GO TO 10
     CASE (    6)
       GO TO 20
   END SELECT
   5 rc = 2
   CALL mtrxo (-1,fss,km(i),0,rc)
   CYCLE
   10 CALL smsg (rc-2,km(i),fss)
   CYCLE
   15 CALL DELETE (fss,km(i),rc)
   GO TO 3
   20 CONTINUE
 END DO
 IF (dry < 0.0) THEN
   GO TO   440
 END IF
 
!     IF MODAL RECOVER, COPY PHIS ITEM TO UVEC
 
 70 IF (.NOT.mrecvr) GO TO 90
 rfno = 3
 CALL mtrxi (scr1,fss,phis,0,rc)
 IF (rc == 1) GO TO 80
 CALL smsg (rc-2,phis,fss)
 GO TO 9100
 80 CALL mtrxo (scr1,fss,uvec,0,rc)
 
!     ATTEMPT TO FETCH SOLN ITEM FOR FSS.  IF IT ALREADY EXISTS, RETURN
 
 90 CALL sfetch (fss,soln,schk,rc)
 IF (rc == 1) GO TO 440
 IF (rc == 3) GO TO 100
 CALL smsg (rc-2,soln,fss)
 GO TO 440
 
!     CREATE SOLN ITEM FOR PROPER RIGID FORMAT
 
 100 IF (rfno < 0 .OR. rfno > 9) GO TO 9007
 SELECT CASE ( rfno )
   CASE (    1)
     GO TO 110
   CASE (    2)
     GO TO 110
   CASE (    3)
     GO TO 130
   CASE (    4)
     GO TO 9007
   CASE (    5)
     GO TO 9007
   CASE (    6)
     GO TO 9007
   CASE (    7)
     GO TO 9007
   CASE (    8)
     GO TO 180
   CASE (    9)
     GO TO 180
 END SELECT
 
!     STATIC SOLUTION - R.F. 1 AND 2
 
 110 CALL rcovss
 GO TO 440
 
!     MODAL SOLUTION - R.F. 3
 
 130 CALL rcovms
 GO TO 440
 
!     DYNAMIC SOLUTION - R.F. 8 AND 9
 
 180 CALL rcovds
 GO TO 440
 
!     FINISHED
 
 440 CALL sofcls
 RETURN
 
!     DIAGNOSTICS
 
 9007 CALL mesage (7,0,NAME)
 9100 iopt = -1
 CALL sofcls
 RETURN
END SUBROUTINE rcova
