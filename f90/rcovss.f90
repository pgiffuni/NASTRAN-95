SUBROUTINE rcovss
     
!     THIS ROUTINE GENERATES THE STATIC SOLUTION ITEM FOR RIGID FORMATS
!     1 AND 2
 
 INTEGER :: dry        ,step       ,fss        ,rfno       ,  &
     sysbuf     ,iz(5)      ,rd         ,rdrew      ,  &
     wrt        ,wrtrew     ,rew        ,eofnrw     ,  &
     lod(4)     ,soln       ,eqss       ,loadc(2)   ,  &
     srd        ,swrt       ,eog        ,eoi        ,  &
     casess     ,geom4      ,scr1       ,rc         ,  &
     buf1       ,buf2       ,buf3       ,cc(2)      ,  &
     FILE       ,NAME(2)    ,casecc(2)
 REAL :: clod(4)
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm        ,uwm
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
 COMMON /names / rd         ,rdrew      ,wrt        ,wrtrew     ,  &
     rew        ,norew      ,eofnrw
 EQUIVALENCE     (z(1),iz(1)), (lod(1),clod(1))
 DATA    NAME  / 4HRCOV,4HSS   /
 DATA    soln  , eqss, lods / 4HSOLN,4HEQSS,4HLODS /
 DATA    casess, geom4,scr1 / 101,102,301 /
 DATA    loadc / 500  ,5    /
 DATA    srd   , swrt, eog,eoi / 1,2,2,3  /
 DATA    casecc/ 4HCASE,4HCC   /
 
!     CREATE SOLN FOR RIGID FORMAT 1 OR 2
 
!     GET NUMBER OF BASIC SUBSTRUCTURES (NS) FROM EQSS AND CREATE
!     GROUP 0 OF SOLN AT TOP OF OPEN CORE
 
 CALL sfetch (fss,eqss,srd,rc)
 IF (rc == 1) GO TO 110
 CALL smsg (rc-2,eqss,fss)
 GO TO 440
 110 CALL suread (z,2,nwds,rc)
 CALL suread (ns,1,nwds,rc)
 IF (lcore < 3*ns+5) GO TO 9008
 CALL suread (z,1,nwds,rc)
 iz(1) = fss(1)
 iz(2) = fss(2)
 iz(3) = rfno
 iz(4) = ns
 
!     GET SUBSTRUCTURE NAMES FROM EQSS
 
 DO  i = 1,ns
   CALL suread (z(3*i+3),2,nwds,rc)
 END DO
 
!     COUNT NUMBER OF SUBCASES (NC) ON CASECC
 
 CALL gopen (casess,z(buf2),rdrew)
 nskip = 1
 130 CALL fread (casess,cc,2,1)
 nskip = nskip + 1
 IF (cc(1) /= casecc(1)) GO TO 130
 IF (cc(2) /= casecc(2)) GO TO 130
 nc = 0
 140 CALL fwdrec (*150,casess)
 nc = nc + 1
 GO TO 140
 150 CALL REWIND (casess)
 iz(5) = nc
 
!     GET NUMBER OF LOAD VECTORS FOR EACH SUBSTRUCTURE FROM LODS
 
 CALL sfetch (fss,lods,srd,rc)
 IF (rc == 1) GO TO 160
 CALL smsg (rc-2,lods,fss)
 GO TO 9200
 160 j = 1
 CALL sjump (j)
 DO  i = 1,ns
   CALL suread (z(3*i+5),1,nwds,rc)
   CALL sjump (j)
 END DO
 
!     SOLN GROUP 0 COMPLETE.  WRITE IT ON SCR1
 
 j = 3
 CALL gopen (scr1,z(buf3),wrtrew)
 CALL WRITE (scr1,z,3*ns+5,1)
 
!     COMPRESS SUBSTRUCTURE NAMES AT TOP OF OPEN CORE
 
 DO  i = 1,ns
   iz(2*i-1) = iz(3*i+3)
   iz(2*i  ) = iz(3*i+4)
 END DO
 
!     PREPARE TO LOOP OVER ALL SUBCASES
 
 icase = 2*ns + 1
 ilods = icase + 166
 IF (ilods > lcore) GO TO 9008
 lodsin= 0
 nlods = ilods - 1
 FILE  = casess
 DO  i = 1,nskip
   CALL fwdrec (*9002,casess)
 END DO
 noldc = 1
 CALL preloc (*195,z(buf1),geom4)
 CALL locate (*195,z(buf1),loadc,i)
 noldc = 0
 
!     BEGIN SUBCASE LOOP.  FOR EACH SUBCASE, BUILD ONE GROUP OF SOLN
 
 195 DO  isc = 1,nc
   CALL fread (casess,z(icase),166,0)
   nlds = 0
   IF (iz(icase+15) /= 0) GO TO 310
   FILE = casess
   CALL fwdrec (*9002,casess)
   FILE = geom4
   
!     PROCESS REGULAR SUBCASE.  IF LODS ITEM NOT IN CORE, GET IT.
   
   IF (iz(icase+3) == 0) GO TO 300
   IF (noldc  == 1) GO TO 300
   IF (lodsin == 1) GO TO 205
   CALL sfetch (fss,lods,srd,rc)
   i = 1
   CALL sjump (i)
   i = ilods
   DO  j = 1,ns
     CALL suread (z(i),1,nwds,rc)
     nlods = i+iz(i)
     IF (nlods > lcore) GO TO 9008
     CALL suread (z(i+1),-1,nwds,rc)
     i = nlods + 1
   END DO
   lodsin = 1
   
!     LODS ITEM IN CORE.  FIND MATCH ON LOADC CARD WITH LOAD SET ID
!     FROM CASECC
   
   205 jsoln = nlods + 2
   210 CALL READ (*9002,*300,geom4,lod,2,0,nwds)
   IF (lod(1) == iz(icase+3)) GO TO 230
   220 CALL fread (geom4,lod,4,0)
   IF (lod(4) == -1) GO TO 210
   GO TO 220
   
!     FOUND MATCH ON LOADC CARD
   
   230 sfac = clod(2)
   
!     LOOP OVER BASIC SUBSTRUCTURES ON THE LOADC CARD
   
   240 CALL fread (geom4,lod,4,0)
   IF (lod(4)  ==    -1) GO TO 290
   IF (jsoln+1 > lcore) GO TO 9008
   
!     FIND BASIC SUBSTRUCTURE NUMBER BY MATCHING ITS NAME WITH THOSE
!     FROM EQSS.  THEN DETERMINE LOAD VECTOR NUMBER BY MATCHING THE
!     BASIC SUBSTRUCTURE LOAD SET ID WITH THOSE IN LODS DATA IN CORE.
   
   DO  i = 1,ns
     IF (lod(1) /= iz(2*i-1)) CYCLE
     k = i
     IF (lod(2) == iz(2*i)) GO TO 250
   END DO
   WRITE (nout,6315) uwm,lod(1),lod(2),lod(3),fss
   250 n = 0
   i = ilods
   j = 1
   260 IF (j == k) GO TO 265
   n = n + iz(i)
   i = i + iz(i) + 1
   j = j + 1
   GO TO 260
   265 j = iz(i)
   DO  k = 1,j
     n = n + 1
     IF (iz(i+k) == lod(3)) GO TO 280
   END DO
   WRITE (nout,6316) uwm,lod(3),lod(1),lod(2),fss
   
!     BUILD SOLN GROUP IN OPEN CORE FOLLOWING LODS DATA
   
   280 iz(jsoln ) = n
   z(jsoln+1) = sfac*clod(4)
   jsoln = jsoln + 2
   nlds  = nlds  + 1
   GO TO 240
   290 iz(nlods+1) = nlds
   jsoln = nlods + 1
   GO TO 385
   
!     NO LOADS FOR THIS SUBCASE
   
   300 nlds = 0
   GO TO 290
   
!     PROCESS SYMCOM OR SUBCOM SUBCASE
   
!     READ SYMSEQ OR SUBSEQ INTO OPEN CORE AT ISEQ
   
   310 lcc   = iz(icase+165)
   lskip = 167 - lcc
   CALL fread (casess,0,lskip,0)
   CALL fread (casess,lsem,1,0)
   320 IF (lsem+nlods < lcore) GO TO 340
   IF (lodsin == 0) GO TO 9008
   
!     SHORT OF CORE.  WIPE OUT LODS DATA AND RE-USE SPACE
   
   330 lodsin = 0
   nlods  = ilods - 1
   GO TO 320
   340 iseq = nlods + 1
   CALL fread (casess,z(iseq),lsem,1)
   
!     READ THE PREVIOUS LSEM GROUPS OF SOLN INTO OPEN CORE FOLLOWING SEQ
   
   jsoln = iseq + lsem
   k = jsoln + 1
   CALL CLOSE (scr1,eofnrw)
   FILE = scr1
   CALL OPEN (*9001,scr1,z(buf3),rd)
   nrec = 1
   nlds = 0
   DO  i = 1,lsem
     342 DO  j = 1,nrec
       CALL bckrec (scr1)
     END DO
     CALL fread (scr1,n,1,0)
     nrec = 2
     IF (n < 0) GO TO 342
     IF (k+2*n-1 < lcore) GO TO 360
     IF (lodsin == 0) GO TO 9008
     
!     SHORT OF CORE.  REPOSITION CASESS, WIPE OUT LODS DATA, AND TRY
!     AGAIN
     
     CALL bckrec (casess)
     CALL fread (casess,0,-166,0)
     GO TO 330
     360 CALL fread (scr1,z(k),2*n,1)
     
!     SCALE LOAD FACTORS BY SYMSEQ OR SUBSEQ FACTORS
     
     DO  j = 1,n
       z(k+2*j-1) = z(iseq+lsem-i)*z(k+2*j-1)
     END DO
     k = k + 2*n
     nlds = nlds + n
   END DO
   iz(jsoln) = -nlds
   
!     COMBINATION GROUP COMPLETE.  REPOSITION SCR1
   
   FILE = scr1
   381 CALL fwdrec (*382,scr1)
   GO TO 381
   382 CALL skpfil (scr1,-1)
   CALL CLOSE (scr1,norew)
   CALL OPEN (*9001,scr1,z(buf3),wrt)
   
!     GROUP COMPLETE IN CORE.  SORT ON LOAD VECTOR NUMBERS
   
   385 CALL sort (0,0,2,1,z(jsoln+1),2*nlds)
   
!     WRITE GROUP ON SCR1 AND POSITION GEOM4 TO BEGINNING OF LOADC CARDS
   
   CALL WRITE (scr1,z(jsoln),2*nlds+1,1)
   IF (noldc == 1) CYCLE
   CALL bckrec (geom4)
   CALL fread (geom4,0,-3,0)
   
!     END OF LOOP OVER SUBCASES
   
 END DO
 CALL CLOSE (casess,rew)
 CALL CLOSE (geom4,rew)
 CALL CLOSE (scr1,rew)
 
!     COPY SOLN FROM SCR1 TO SOF
 
 CALL gopen (scr1,z(buf1),rdrew)
 rc = 3
 CALL sfetch (fss,soln,swrt,rc)
 392 CALL READ (*396,*394,scr1,z,lcore,1,nwds)
 GO TO 9008
 394 CALL suwrt (z,nwds,eog)
 GO TO 392
 396 CALL CLOSE (scr1,rew)
 
!     FINISH
 
 CALL suwrt (0,0,eoi)
 440 CONTINUE
 RETURN
 
!     DIAGNOSTICS
 
 6315 FORMAT (a25,' 6315, RCOVR MODULE IS UNABLE TO FIND SUBSTRUCTURE ',  &
     2A4,' AMONG THOSE ON EQSS.' /32X,'LOAD SET',i9,  &
     ' FOR THAT SUBSTRUCTURE WILL BE IGNORED IN CREATING', /32X,  &
     'THE SOLN ITEM FOR FINAL SOLUTION STRUCTURE ',2A4)
 6316 FORMAT (a25,' 6316, RCOVR MODULE IS UNABLE TO FIND LOAD SET',i9,  &
     ' FOR SUBSTRUCTURE ',2A4, /32X,'AMONG THOSE ON LODS.  ',  &
     'IT WILL BE IGNORED IN CREATING THE SOLN ITEM FOR FINAL',  &
     /32X,'SOLUTION STRUCTURE ',2A4)
 
 9001 n = 1
 GO TO 9100
 9002 n = 2
 GO TO 9100
 9008 n = 8
 9100 CALL mesage (n,FILE,NAME)
 9200 CALL sofcls
 iopt = -1
 CALL CLOSE (casess,rew)
 CALL CLOSE (geom4,rew)
 CALL CLOSE (scr1,rew)
 
 RETURN
END SUBROUTINE rcovss
