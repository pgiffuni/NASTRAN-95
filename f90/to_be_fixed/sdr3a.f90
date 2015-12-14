SUBROUTINE sdr3a (ofpfil)
     
!     SORT-2  MODULE
 
 
 INTEGER, INTENT(OUT)                     :: ofpfil(6)
 INTEGER :: iname(2),trail(7),id(146),idtemp(146),vector(50),  &
     scrtch(8),ofile(6),ifile(6),buff(10), z,  &
     words,eof,core,buff9,buff10,group,oufile,FILE,recs  &
     ,               recpt,eor,outrwd,rwd,full,ovrlap,v in bk,w per bk,  &
     v per bk,total1,total2,ahead,entrys(85)
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm
 COMMON /system/ ibufsz, l
 COMMON /zzzzzz/ z(1)
 EQUIVALENCE     (nwds,id(10))
 
!     IF THE NUMBER OF SCRATCH FILES CHANGE, ONE SHOULD SET NSCRAT EQUAL
!     TO THE NEW NUMBER AND INCREASE THE DATA BELOW
 
!     NFILES BELOW EQUALS THE NUMBER OF INPUT FILES AND ALSO EQUALS
!     THE NUMBER OF OUTPUT FILES.  IF NFILES CHANGES, CHANGE THE DATA
!     BELOW TO CONFORM...
!     ALSO CHANGE DIMENSIONS OF BUFF,IFILE,OFILE,SCRTCH, AS REQUIRED...
 
 DATA    ifile / 101,102,103,104,105,106 /
 DATA    ofile / 201,202,203,204,205,206 /
 DATA    scrtch/ 301,302,303,304,305,306,307,308 /
 DATA    trail / 0,1,2,3,4,5,6 /
 DATA    eor   , noeor,rwd,inprwd,outrwd / 1,0,1,0,1 /
 
 nfiles = 6
 nscrat = 8
 DO  i = 1,6
   ofpfil(i) = 0
 END DO
 DO  i = 1,146
   idtemp(i) = 0
 END DO
 
!     BUFFERS AND OPEN CORE
 
 core = korsz(z)
 
 buff(1) = core - ibufsz + 1
 DO  i = 2,10
   buff(i) = buff(i-1) - ibufsz
 END DO
 buff9  = buff( 9)
 buff10 = buff(10)
 core   = buff(10) - 1
 IF (core < 1) GO TO 700
 
!     OPEN SCRATCH FILES FOR OUTPUT
 
 ierror = 0
 DO  i = 1,nscrat
   ibuff = buff(i)
   CALL OPEN (*20,scrtch(i),z(ibuff),outrwd)
   CYCLE
   20 ierror = 1
   WRITE  (l,21) uwm,i
   21 FORMAT (a25,' 985, SDR3 FINDS SCRATCH',i1,' PURGED.')
 END DO
 
!     EXECUTE FOR NFILES FILES
 
 DO  FILE = 1,nfiles
   eof = 0
   infile = ifile(FILE)
   oufile = ofile(FILE)
   
   CALL OPEN (*460,infile,z(buff9),inprwd)
   CALL OPEN (*480,oufile,z(buff10),outrwd)
   CALL fwdrec (*520,infile)
   
!     HEADER RECORD FOR OUFILE
   
   CALL fname (oufile,iname(1))
   CALL WRITE (oufile,iname(1),2,eor)
   
!     WRITE SOME JUNK IN TRAILER FOR NOW
   
   trail(1) = oufile
   CALL wrttrl (trail(1))
   no frq = 0
   
!     PROCEED WITH TRANSPOSE OF DATA = SORT-2
   
!     GROUP WILL BE THE NUMBER OF THE FIRST REC IN THE PRESENT GROUP OF
!     DATA BLOCKS BEING OPERATED ON, LESS 1
   
   nrecs = 1
   
   50 ASSIGN 120 TO iretrn
   
   60 CALL READ (*80,*620,infile,id(1),146,eor,iamt)
   IF (id(1)/10 == 1) nofrq = 1
   idata = 1
   
   70 icore = core
   recs  = 0
   group = nrecs
   
!     READ FIRST DATA BLOCK INTO CORE
   
   CALL READ (*530,*90,infile,z(1),icore,noeor,iamt)
   
!     INSUFFICIENT CORE,  IF FALL HERE, TO DO SORT II ON THIS FILE..
   
   GO TO 670
   80 CALL CLOSE (infile,rwd)
   CALL CLOSE (oufile,rwd)
   CYCLE
   
   90 IF (iamt == 0) GO TO 441
   entrys(1) = iamt/nwds
   
!     SET UP IN-CORE ENTRY BLOCKS
!     SPOT FOR TRANSPOSE HEADING DATA IS AT Z(ICORE-ENTRYS(1)+1)
   
   ihd2  = icore + 1
   icore = icore - entrys(1)
   ihead = icore
   IF (icore < iamt) GO TO 680
   
!     NOTATION - W PER BK = WORDS PER ENTRY BLOCK
!                V PER BK = VECTORS PER ENTRY BLOCK
!                V IN  BK = VECTORS NOW IN ENTRY BLOCKS
   
   w per bk = icore/entrys(1)
   v per bk = w per bk/nwds
   w per bk = v per bk*nwds
   IF (v per bk < 1) GO TO 690
   
!     DISTRIBUTE FIRST DATA BLOCK TO INCORE ENTRY BLOCKS (BOTTOM TO TOP)
   
   nentry = entrys(1)
   total1 = w per bk*entrys(1) + 1
   total2 = nwds*entrys(1) + 1
   DO  i = 1,nentry
     n1  = total1 - w per bk*i
     n2  = total2 - nwds    *i
     ihd = ihd2 - i
     z(ihd) = z(n2)
     z(n1 ) = id(5)
     
!     SAVE TRANSPOSE HEADING
     
     DO  j = 2,nwds
       n1 = n1 + 1
       n2 = n2 + 1
       z(n1) = z(n2)
     END DO
   END DO
   
   v in bk = 1
   GO TO iretrn, (120,150)
   
   120 ntypes = 1
   130 CALL READ (*159,*630,infile,idtemp(1),146,eor,iamt)
   IF ((id(2) == idtemp(2) .AND. id(3) == idtemp(3) .AND.  &
       id(5) /= idtemp(5)) .OR. (id(5) /= idtemp(5)) .OR.  &
       (id(4) /= idtemp(4))) GO TO 160
   
   ntypes = ntypes + 1
   nwords = idtemp(10)
   
!     WILL READ DATA AND COUNT ENTRYS
   
   IF (ntypes > 30) GO TO 472
   entrys(ntypes) = 0
   140 CALL READ (*540,*130,infile,idtemp(1),nwords,noeor,iamt)
   entrys(ntypes) = entrys(ntypes) + 1
   GO TO 140
   
   150 ahead = 2*ntypes - 2
   IF (ndata == 1) GO TO 260
   GO TO 170
   
!     AT THIS POINT IT IS KNOWN HOW MANY TYPES ARE IN THE PRESENT GROUP
!     OF DATA BLOCKS AND ALSO HOW MANY ENTRYS IN EACH TYPE
   
   159 IF (ntypes == 1) eof = 1
   160 itype = 1
   ndata = 1
   idata = 1
   
!     POSITION TO READ 2-ND ID OF TYPE(ITYPE) IF NOT JUST READ
   
   IF (ntypes == 1) GO TO 200
   CALL REWIND (infile)
   ahead = group + 2*ntypes
   170 DO  i = 1,ahead
     CALL fwdrec (*550,infile)
   END DO
   
   190 CALL READ (*260,*640,infile,idtemp(1),146,eor,iamt)
   
!     CHECK FOR BREAK POINT
   
   200 IF (nofrq == 1) GO TO 201
   IF (id(4) /= idtemp(4)) GO TO 270
   201 CONTINUE
   IF (eof == 1) GO TO 270
   IF (itype == 1) ndata = ndata + 1
   idata  = idata + 1
   nentry = entrys(itype)
   
!     CHECK TO SEE IF THERE IS ENOUGH ROOM IN EACH OF THE INCORE
!     ENTRY BLOCKS FOR ANOTHER VECTOR
!     IF NOT DO SCRATCH FILE  OPERATIONS
   
   IF (v in bk < v per bk) GO TO 220
   
!     NOT ENOUGH ROOM THUS DUMP CORE ENTRY BLOCKS ONTO SCRATCH FILES
   
   IF (ierror == 1) GO TO 451
   npoint = 1
   nfile  = nscrat
   DO  i = 1,nentry
     nfile  = nfile + 1
     IF (nfile > nscrat) nfile = 1
     CALL WRITE (scrtch(nfile),z(npoint),w per bk,eor)
     npoint = npoint + w per bk
   END DO
   recs = recs + nentry
   
!     IN CORE ENTRY BLOCKS ARE NOW EMPTY
   
   v in bk = 0
   
!     DISTRIBUTE DATA TO INCORE ENTRY BLOCKS
   
   220 npoint = v in bk*nwds + 1
   DO  i = 1,nentry
     IEOR = i/nentry
     CALL READ (*560,*650,infile,z(npoint),nwds,IEOR,iamt)
     z(npoint) = idtemp(5)
     npoint = npoint + w per bk
   END DO
   v in bk = v in bk + 1
   
   IF (ntypes == 1) GO TO 190
   IF (itype  == 1) GO TO 240
   IF (idata == ndata) GO TO 270
   
!     NOW POSITION AHEAD TO READ NEXT ID FOR TYPE(ITYPE)
   
   240 ahead = 2*ntypes - 2
   DO  i = 1,ahead
     CALL fwdrec (*570,infile)
   END DO
   GO TO 190
   
!     ONE DATA TYPE IN THIS GROUP IS COMPLETE
   
!     OUTPUT IS IN CORE, AND ON SCRATCH FILES IF RECS IS NOT 0
   
!     NOW DUMP SCRATCHES AND (OR JUST) CORE ONTO FINAL OUTPUT TAPE
   
!     ID WILL BE WRITTEN BEFORE EACH ENTRY INSERTING INTO IT THE NEW
!     HEADER VALUE REPLACING FREQUENCY OR TIME ETC
   
   260 eof = 1
   270 IF (recs == 0) GO TO 290
   
   layers = recs/nentry
   
!     CLOSE SCRATCH FILES AND OPEN AS INPUT FILES
   
   DO  i = 1,nscrat
     CALL CLOSE (scrtch(i),rwd)
     ibuff = buff(i)
     CALL OPEN (*500,scrtch(i),z(ibuff),inprwd)
   END DO
   
!     COMPUTE OVERLAPS PER LAYER
   
   ovrlap = (nentry-1)/nscrat
   
!     COMPUTE HOW MANY TAPES HAVE ALL THE OVERLAPS
   
   full = nentry - ovrlap*nscrat
   
   
!     WRITE FINAL FILE THEN
   
   290 nfile = 0
   id(2) = id(2) + 2000
   DO  i = 1,nentry
     nfile = nfile + 1
     IF (nfile > nscrat) nfile = 1
     
     npoint = ihead + i
     id(5)  = z(npoint)
     CALL WRITE (oufile,id(1),146,eor)
     
!     ANYTHING ON SCRATCH FILES IS NOW WRITTEN
     
     IF (recs == 0) GO TO 390
     
     DO  j = 1,layers
       
!     FORWARD REC IF NECESSARY
       
       IF (j > 1) GO TO 320
       
!     AHEAD TO FIRST PART IF NECESSARY
       
       IF (layers == 1) GO TO 350
       ahead = (i-1)/nscrat
       
       IF (ahead == 0.0) THEN
         GO TO   350
       END IF
       300 DO  k = 1,ahead
         CALL fwdrec (*590,scrtch(nfile))
       END DO
       GO TO 350
       
       320 recpt = ovrlap
       IF (nfile > full) recpt = recpt - 1
       IF (recpt > 0.0) THEN
         GO TO   330
       ELSE
         GO TO   350
       END IF
       330 DO  k = 1,recpt
         CALL fwdrec (*600,scrtch(nfile))
       END DO
       
!     COPY RECORD FROM SCRTCH TO OUTFILE
       
       350 DO  k = 1,vperbk
         IEOR = k/v per bk
         CALL READ (*610,*660,scrtch(nfile),vector(1),nwds,IEOR,iamt)
         CALL WRITE (oufile,vector(1),nwds,noeor)
       END DO
     END DO
     IF (layers > 1) CALL REWIND (scrtch(nfile))
     
!     COPY INCORE VECTORS TO OUTFILE
     
     390 words  = v in bk*nwds
     npoint = w per bk*i - w per bk + 1
     CALL WRITE (oufile,z(npoint),words,eor)
   END DO
   IF (recs == 0) GO TO 420
   
!     CLOSE SCRTCH FILES AND OPEN AS OUTPUT FILES
   
   DO  i = 1,nscrat
     CALL CLOSE (scrtch(i),rwd)
     ibuff = buff(i)
     CALL OPEN (*490,scrtch(i),z(ibuff),outrwd)
   END DO
   
   420 IF (itype == ntypes) GO TO 440
   
   itype = itype + 1
   CALL REWIND (infile)
   ahead = group + itype*2 - 2
   DO  i = 1,ahead
     CALL fwdrec (*580,infile)
   END DO
   ASSIGN 150 TO iretrn
   eof = 0
   GO TO 60
   
!     THIS GROUP IS ABSOLUTELY COMPLETE AND WE ARE AT BREAK POINT
   
   440 IF (eof == 1) GO TO 80
   nrecs = nrecs + 2*ndata*ntypes
   IF (ntypes > 1) GO TO 50
   441 CONTINUE
   DO  i = 1,146
     id(i) = idtemp(i)
   END DO
   GO TO 70
   
   
!     ERROR CONDITIONS FOR THIS DATA BLOCK
   
!     FORMAT OF INPUT DATA BLOCK MAY BE INCORRECT (N=TRACEBACK CODE)
   
   490 n = 23
   GO TO 452
   500 n = 3
   GO TO 452
   520 n = 4
   GO TO 452
   530 n = 5
   GO TO 452
   540 n = 6
   GO TO 452
   550 n = 7
   GO TO 452
   560 n = 8
   GO TO 452
   570 n = 9
   GO TO 452
   580 n = 10
   GO TO 452
   590 n = 11
   GO TO 452
   600 n = 12
   GO TO 452
   610 n = 13
   GO TO 452
   620 n = 14
   GO TO 452
   630 n = 15
   GO TO 452
   640 n = 16
   GO TO 452
   650 n = 17
   GO TO 452
   660 n = 18
   GO TO 452
   452 ofpfil(FILE) = n
   WRITE  (l,453) uwm,FILE
   453 FORMAT (a25,' 982, FORMAT OF SDR3 INPUT DATA BLOCK ',i3,  &
       ' DOES NOT PERMIT SUCCESSFUL SORT-2 PROCESSING.')
   GO TO 80
   472 WRITE  (l,475) ufm,ntypes
   475 FORMAT (a23,' 3129, SDR3 CAN ONLY PROCESS 30 ELEMENT TYPES, ',  &
       'PROBLEM HAS',i5)
   CALL mesage (-61,0,0)
   
!     CORRESPONDING OUTPUT FILE IS PURGED.
   
   480 ofpfil(FILE) = 2
   WRITE  (l,481) uwm,FILE
   481 FORMAT (a25,' 984,  SDR3 FINDS OUTPUT DATA-BLOCK',i4,' PURGED.')
   GO TO 80
   
!     ATTEMPT TO USE SCRATCH FILES 1 OR MORE OF WHICH ARE PURGED.
   
   451 ofpfil(FILE) = 1
   GO TO 80
   
!     INSUFFICIENT CORE
   
   670 n = 19
   GO TO 701
   680 n = 20
   GO TO 701
   690 n = 21
   GO TO 701
   701 WRITE  (l,702) uwm,FILE
   702 FORMAT (a25,' 983, SDR3 HAS INSUFFICIENT CORE TO PERFORM SORT-2',  &
       ' ON INPUT DATA BLOCK',i4, /5X, 'OR DATA-BLOCK IS NOT IN CORRECT FORMAT.')
   ofpfil(FILE) = n
   GO TO 80
   
 END DO
 
!     CLOSE SCRATCH FILES
 
 DO  i = 1,nscrat
   CALL CLOSE (scrtch(i),rwd)
 END DO
 
 GO TO 801
 700 DO  i = 1,5
   ofpfil(i) = 22
 END DO
 WRITE  (l,703) uwm
 703 FORMAT (a25,' 986, INSUFFICIENT CORE FOR SDR3.')
 801 RETURN
END SUBROUTINE sdr3a
