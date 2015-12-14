SUBROUTINE gp3a
     
!     GP3A BUILDS THE STATIC LOADS TABLE (SLT).
!     FORCE, FORCE1, FORCE2, MOMENT, MOMNT1, MOMNT2, GRAV, PLOAD, SLOAD
!     AND LOAD CARDS ARE READ. EXTERNAL GRID NOS. ARE CONVERTED TO
!     INTERNAL INDICES. EACH LOAD SET ID (EXCEPT ON LOAD CARD) IS
!     WRITTEN IN THE HEADER RECORD OF THE SLT. THE SLT THEN COMPRISES
!     ONE LOGICAL RECORD PER LOAD SET. THE LAST RECORD OF THE SLT
!     CONTAINS THE LOAD CARDS. RFORCE CARD ADDED IN AUGUST, 1968.
!     PLOAD3 CARD ADDED ON HALLOWEEN 1972
 
 LOGICAL :: piez
 INTEGER :: geom3 ,eqexin,slt   ,gptt  ,buf1  ,buf2  ,buf   ,  &
     z     ,rd    ,rdrew ,wrt   ,wrtrew,clsrew,cardid,  &
     carddt,STATUS,FILE  ,gpoint,scr1  ,scr2  ,first ,  &
     setid ,flag  ,nam(2),ksystm(80)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /system/ isb   ,iptr  ,idm(6),nlpp  ,idum(2),lines
 COMMON /BLANK / nograv,noload,notemp
 COMMON /gp3com/ geom3 ,eqexin,geom2 ,slt   ,gptt  ,scr1  ,scr2  ,  &
     buf1  ,buf2  ,buf(50)      ,cardid(60)   ,idno(30)  &
     ,carddt(60)   ,mask(60)     ,STATUS(60)   ,ntypes,  &
     ipload,igrav ,pload2(2)    ,load(2)      ,nopld2,  &
     temp(2)      ,tempd(2)     ,tempp1(2)    ,  &
     tempp2(2)    ,tempp3(2)    ,temprb(2)    ,buf3  , pload3(2)    ,ipld3
 COMMON /names / rd    ,rdrew ,wrt   ,wrtrew,clsrew
 COMMON /zzzzzz/ z(1)
 EQUIVALENCE     (ksystm(1),isb)
 DATA    nam   / 4HGP3A,4H    / ,irfrc / 9  /
 
!     READ EQEXIN INTO CORE. INITIALIZE BINARY SEARCH ROUTINE.
 
 FILE = eqexin
 CALL OPEN   (*570,eqexin,z(buf1),rdrew)
 CALL fwdrec (*580,eqexin)
 CALL READ   (*580,*20,eqexin,z,buf2,1,neqx)
 CALL mesage (-8,0,nam)
 20 CALL CLOSE  (eqexin,clsrew)
 kn   = neqx/2
 nogo = 0
 
!     INITIALIZE POINTERS AND OPEN SCR1 AND GEOM3.
 
 iset  = buf2 - 2
 kset  = iset
 ilist = neqx + 1
 klist = ilist
 ktabl = 1
 first = 1
 FILE  = scr1
 CALL OPEN (*570,scr1,z(buf2),wrtrew)
 
!     IF PLOAD2 CARDS PRESENT, INITIALIZE TO READ PLOAD DATA FROM SCR2
!     INSTEAD OF GEOM3.
 
 IF (nopld2 == 0) GO TO 40
 FILE = scr2
 CALL OPEN (*570,scr2,z(buf1),rdrew)
 GO TO 60
 40 first = 0
 50 FILE  = geom3
 CALL OPEN   (*570,geom3,z(buf1),rdrew)
 CALL fwdrec (*580,geom3)
 
!     READ 3-WORD RECORD ID. IF ID BELONGS TO LOAD SET, TURN NOLOAD FLAG
!     OFF.
!     SET 1ST WORD IN STATUS ENTRY TO CURRENT POINTER IN LIST TABLE.
!     SET PARAMETERS FOR CONVERSION OF GRID NOS. TO INTERNAL INDICES.
 
 60 CALL READ (*170,*60,FILE,buf,3,0,flag)
 DO  i = 1,ntypes,2
   IF (buf(1) == cardid(i) .AND. buf(2) == cardid(i+1)) GO TO 90
 END DO
 80 CALL fwdrec (*170,FILE)
 GO TO 60
 90 noload = 1
 IF (first == 1) GO TO 100
 
!     IF I POINTS TO PLOAD RECORD AND PLOAD2 CARDS ARE PRESENT, THEN
!     PLOAD DATA IS ALREADY PROCESSED. IN THIS CASE, SKIP PLOAD RECORD.
!     IF I POINTS TO PLOAD3 RECORD ON GEOM3, SKIP RECORD.
 
 IF (i == ipload .AND. nopld2 /= 0 .AND. nopld2 /= 2) GO TO 80
 IF (i == ipld3) GO TO 80
 100 CONTINUE
 STATUS(i) = klist - ilist + 1
 nwds  = carddt(i)
 nwds1 = nwds - 1
 jx  = carddt(i+1)
 jj1 = jx + 1
 jjn = jx + mask(jx)
 id  = 0
 
!     READ A LOAD CARD. IF SET ID IS DIFFERENT FROM LAST READ (OR 1ST
!     ONE) STORE SET ID IN POINTER LIST AND IN SET LIST. STORE POINTER
!     IN POINTER LIST. IF NOT FIRST CARD OF TYPE, STORE WORD COUNT IN
!     POINTER LIST.
 
 110 CALL READ (*580,*160,FILE,buf,nwds,0,flag)
 IF (buf(1) == id) GO TO 120
 z(klist  ) = buf(1)
 z(klist+1) = ktabl
 IF (id /= 0) z(klist-1) = n
 id = buf(1)
 n  = 0
 klist   = klist + 3
 z(kset) = buf(1)
 kset    = kset - 1
 
!     CONVERT EXTERNAL GRID NOS. ON CARD TO INTERNAL NOS. INCREMENT
!     WORD COUNT. WRITE LOAD CARD (WITHOUT SET ID) ON SCR1.
 
 120 IF (jx == 0) GO TO 150
 jj = jj1
 jstop = 0
 130 IF (jstop == 0) GO TO 135
 jx = jx + 1
 GO TO 136
 135 jx = mask(jj)
 IF (jx > 0) GO TO 136
 jx = -jx
 jstop = 1
 136 gpoint = buf(jx)
 piez = .false.
 IF (gpoint < 0 .AND.  ksystm(78) == 1) piez = .true.
 IF (piez) gpoint = -gpoint
 IF (gpoint == -1 .AND. (cardid(i) == 3209 .OR. cardid(i) == 3409)) GO TO 140
 IF (gpoint /= 0) GO TO 450
 140 IF (piez) gpoint = -gpoint
 buf(jx) = gpoint
 jj = jj + 1
 IF (jj <= jjn) GO TO 130
 
!     CHECK FOR PLOAD4 CARD
 
 150 IF (i /= 49) GO TO 152
 
!     CHECK FOR THRU OPTION ON PLOAD4 CARD
 
 IF (buf(7) == 0) GO TO 153
 
 152 CALL WRITE (scr1,buf(2),nwds1,0)
 GO TO 158
 
!     PROCESS PLOAD4 DATA FOR ALL ELEMENT IDS IMPLIED BY THE THRU OPTION
 
 153 iii = buf(2)
 jjj = buf(8)
 buf(7) =-1
 buf(8) = 0
 DO  kkk = iii,jjj
   buf(2) = kkk
   CALL WRITE (scr1,buf(2),nwds1,0)
   n = n + nwds1
   ktabl = ktabl + nwds1
 END DO
 GO TO 110
 
 158 n = n + nwds1
 ktabl = ktabl + nwds1
 GO TO 110
 
!     HERE WHEN ALL CARDS OF CURRENT CARD TYPE HAVE BEEN READ.
!     STORE WORD COUNT FOR LAST SET IN POINTER LIST. STORE POINTER
!     TO LAST ENTRY FOR CARD TYPE IN 2ND WORD OF STATUS ENTRY.
!     LOOP BACK TO READ NEXT CARD TYPE.
 
 160 z(klist-1) = n
 STATUS(i+1) = klist - ilist - 2
 GO TO 60
 170 IF (first == 0) GO TO 175
 first = 0
 CALL CLOSE (scr2,clsrew)
 GO TO 50
 
!     HERE WHEN END-OF-FILE ON GEOM3 ENCOUNTERED. IF ERROR CONDITION
!     NOTED, CALL PEXIT. IF NO LOAD CARDS FOUND, CLOSE FILES AND RETURN.
 
 175 IF (nogo /= 0) CALL mesage (-61,0,0)
 IF (noload /= -1) GO TO 180
 CALL CLOSE (geom3,clsrew)
 CALL CLOSE (scr1,clsrew)
 RETURN
 
!     IF GRAVITY LOADS WERE READ, TURN NOGRAV FLAG OFF.
!     CLOSE FILES AND MOVE POINTER LIST TO BEGINNING OF CORE.
 
 180 IF (STATUS(igrav) > 0 .OR. STATUS(irfrc) > 0) nograv = +1
 CALL WRITE (scr1,0,0,1)
 CALL CLOSE (geom3,clsrew)
 CALL CLOSE (scr1, clsrew)
 n = klist - ilist
 DO  i = 1,n
   k = ilist + i
   z(i)  = z(k-1)
 END DO
 ilist = 1
 nlist = n - 2
 
!     CHECK UNIQUENESS OF LOAD SETS WITTH RESPECT TO GRAVITY LOAD SETS
 
 IF (STATUS(igrav) < 0) GO TO 200
 k1 = STATUS(igrav  )
 k2 = STATUS(igrav+1)
 DO  i = ilist,nlist,3
   IF (i >= k1 .AND. i <= k2) CYCLE
   setid = z(i)
   DO  k = k1,k2,3
     IF (z(k) /= setid) CYCLE
     nogo = 1
     CALL mesage (30,134,setid)
   END DO
 END DO
 
!     SORT THE SET LIST AND DISCARD DUPLICATE SET NOS.
 
 200 n = iset - kset
 kset = kset + 1
 CALL sort (0,0,1,1,z(kset),n)
 z(iset+1) = 0
 k = nlist + 3
 DO  i = kset,iset
   IF (z(i) == z(i+1)) CYCLE
   z(k) = z(i)
   k = k + 1
 END DO
 iset  = nlist + 3
 nset  = k - 1
 itabl = nset
 
!     OPEN SCRATCH FILE AND SLT FILE.
!     WRITE SET LIST IN HEADER RECORD OF THE SLT.
 
 CALL OPEN (*570,scr1,z(buf1),rdrew)
 FILE = slt
 CALL OPEN (*570,slt,z(buf2),wrtrew)
 CALL fname (slt,buf)
 CALL WRITE (slt,buf,2,0)
 n = nset - iset + 1
 CALL WRITE (slt,z(iset),n,1)
 
!     IF ALL LOAD CARDS WILL FIT IN CORE, READ THEM IN.
 
 nwds  = ktabl - 1
 ncore = itabl + ktabl
 IF (ncore >= buf2) GO TO 370
 FILE  = scr1
 CALL READ (*580,*590,scr1,z(itabl+1),nwds,1,flag)
 CALL CLOSE (scr1,clsrew)
 
!     FOR EACH LOAD SET IN THE SET LIST, LOOP THRU THE STATUS TABLE.
!     FOR EACH CARD TYPE PRESENT IN THE STATUS TABLE, PICK UP POINTERS
!     TO THE POINTER LIST. SEARCH THE POINTER LIST FOR A SET ID MATCH.
!     IF FOUND, PICK UP POINTERS TO THE DATA IN CORE. SORT THE DATA ON
!     INTERNAL INDEX (EXCEPT GRAV AND PLOAD CARDS).
!     THEN, WRITE CARD TYPE ID, NO. OF CARDS IN THE SET, AND THE DATA
!     ON THE CARDS. THUS, THE SLT IS COMPRISED OF ONE LOGICAL RECORD PER
!     SET DATA WITHIN EACH RECORD IS GROUPED BY CARD TYPE, AND, WITHIN
!     THE GROUP, IS SORTED BY INTERNAL INDEX (WHERE DEFINED).
 
 DO  k = iset,nset
   setid = z(k)
   ii = 1
   DO  i = 1,ntypes,2
     IF (STATUS(i) < 0) GO TO 270
     jj1 = STATUS(i  )
     jjn = STATUS(i+1)
     DO  jj = jj1,jjn,3
       IF (z(jj) == setid) GO TO 260
     END DO
     GO TO 270
     
     260 CONTINUE
     jx   = itabl + z(jj+1)
     nwds = z(jj+2)
     n    = carddt(i) - 1
     nkey = 1
     IF (idno(ii) == 20) nkey = 5
     IF (idno(ii) == 21) GO TO 265
     IF (idno(ii) >= 22 .AND. idno(ii) <= 24) GO TO 265
     IF (i == ipload .OR. i == ipld3 .OR. i == igrav) GO TO 265
     CALL sort (0,0,n,nkey,z(jx),nwds)
     265 buf(1) = idno(ii)
     buf(2) = nwds/n
     CALL WRITE (slt,buf,2,0)
     CALL WRITE (slt,z(jx),nwds,0)
     270 ii = ii + 1
   END DO
   CALL WRITE (slt,0,0,1)
 END DO
 
!     IF COMBINATION LOADS ARE PRESENT, SET IDS ARE CHECKED TO ASSURE
!     THAT THEY ARE UNIQUE WITH RESPECT TO LOAD CARDS.  THE SET IDS
!     SPECIFIED ON THE LOAD CARD ARE THEN CHECKED AGAINST THOSE IN THE
!     SET LIST TO VERIFY THAT ALL ARE AVAILABLE AND AGAINST EACH OTHER
!     TO ENSURE THAT NO DUPLICATE SPECIFICATIONS EXIST.  THE COMBINATION
!     LOADS ARE WRITTEN AS THE LAST LOGICAL RECORD OF THE SLT.
 
 290 FILE = geom3
 CALL preloc (*570,z(buf1),geom3)
 CALL locate (*360,z(buf1),load,flag)
 300 CALL READ   (*580,*350,geom3,buf,2,0,flag)
 CALL WRITE  (slt,buf,2,0)
 DO  i = iset,nset
   IF (buf(1) == z(i)) GO TO 330
 END DO
 GO TO 340
 330 nogo = 1
 CALL mesage (30,106,buf)
 340 lset = nset + 1
 mset = nset
 idcmld = buf(1)
 341 CALL READ (*580,*350,geom3,buf,2,0,flag)
 CALL WRITE (slt,buf,2,0)
 IF (buf(1) == -1) GO TO 300
 DO  i = iset,nset
   IF (buf(2) == z(i)) GO TO 343
 END DO
 nogo = 1
 WRITE  (iptr,3178) ufm,buf(2),idcmld
 3178 FORMAT (a23,' 3178, LOAD SET',i9,' NOT FOUND.  REQUIRED FOR ',  &
     'DEFINITION OF COMBINATION LOAD',i9)
 lines = lines + 2
 IF (lines >= nlpp) CALL page
 GO TO 341
 343 IF (mset == nset) GO TO 345
 DO  i = lset,mset
   IF (buf(2) == z(i)) GO TO 346
 END DO
 345 mset = mset + 1
 z(mset) = buf(2)
 GO TO 341
 346 nogo = 1
 WRITE  (iptr,3179) ufm,buf(2),idcmld
 3179 FORMAT (a23,' 3179, DUPLICATE LOAD SET',i9,' FOUND IN DEFINITION',  &
     ' OF COMBINATION LOAD',i9)
 lines = lines + 2
 IF (lines >= nlpp) CALL page
 GO TO 341
 350 CALL WRITE (slt,0,0,1)
 360 CALL CLOSE (geom3,clsrew)
 CALL CLOSE (slt,clsrew)
 buf(1) = slt
 buf(2) = nset - iset + 1
 DO  i = 3,7
   buf(i) = 0
 END DO
 CALL wrttrl (buf)
 IF (nogo /= 0) CALL mesage (-61,0,0)
 RETURN
 
!     HERE IF CORE WILL NOT HOLD ALL LOAD CARDS.
!     CODE IS SIMILAR TO THAT ABOVE EXCEPT THAT POINTER LIST NOW POINTS
!     TO THE DATA ON THE SCRATCH FILE INSTEAD OF IN CORE. THEREFORE, THE
!     SCRATCH FILE WILL HAVE TO BE PASSED ONCE FOR EACH SET IN THE SET
!     LIST.
 
 370 FILE  = scr1
 DO  k = iset,nset
   setid = z(k)
   ii    = 1
   nread = 0
   DO  i = 1,ntypes,2
     IF (STATUS(i) < 0) GO TO 420
     jj1 = STATUS(i  )
     jjn = STATUS(i+1)
     DO  jj = jj1,jjn,3
       IF (z(jj) == setid) GO TO 390
     END DO
     GO TO 420
     390 nskip = z(jj+1) - nread - 1
     nwds  = z(jj+2)
     n = carddt(i) - 1
     IF (nskip < 0) THEN
       GO TO   440
     ELSE IF (nskip == 0) THEN
       GO TO   410
     END IF
     400 CALL READ (*580,*590,scr1,0,-nskip,0,flag)
     410 CALL READ (*580,*590,scr1,z(itabl+1),nwds,0,flag)
     nread = z(jj+1) + nwds - 1
     nkey  = 1
     IF (idno(ii) == 20) nkey = 5
     IF (idno(ii) == 21) GO TO 415
     IF (idno(ii) >= 22 .AND. idno(ii) <= 24) GO TO 415
     IF (i == ipload .OR. i == ipld3 .OR. i == igrav) GO TO 415
     CALL sort (0,0,n,nkey,z(itabl+1),nwds)
     415 buf(1) = idno(ii)
     buf(2) = nwds/n
     CALL WRITE (slt,buf,2,0)
     CALL WRITE (slt,z(itabl+1),nwds,0)
     420 ii = ii + 1
   END DO
   CALL WRITE (slt,0,0,1)
   CALL REWIND (scr1)
 END DO
 CALL CLOSE (scr1,clsrew)
 GO TO 290
 440 CALL mesage (-61,0,0)
 
!     BINARY SEARCH ROUTINE
 
 450 klo = 1
 khi = kn
 460 k = (klo+khi+1)/2
 470 IF (gpoint-z(2*k-1) < 0.0) THEN
   GO TO   480
 ELSE IF (gpoint-z(2*k-1) == 0.0) THEN
   GO TO   540
 ELSE
   GO TO   490
 END IF
 480 khi = k
 GO TO 500
 490 klo = k
 500 IF (khi-klo-1 < 0) THEN
   GO TO   550
 ELSE IF (khi-klo-1 == 0) THEN
   GO TO   510
 ELSE
   GO TO   460
 END IF
 510 IF (k == klo) GO TO 520
 k = klo
 GO TO 530
 520 k = khi
 530 klo = khi
 GO TO 470
 540 gpoint = z(2*k)
 GO TO 140
 550 buf(2) = gpoint
 nogo   = 1
 CALL mesage (30,8,buf)
 GO TO 140
 
!     FATAL FILE ERRORS
 
 560 CALL mesage (n,FILE,nam)
 570 n = -1
 GO TO 560
 580 n = -2
 GO TO 560
 590 n = -3
 GO TO 560
END SUBROUTINE gp3a
