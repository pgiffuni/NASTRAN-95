SUBROUTINE ddrmm
     
!     DYNAMIC-DATA-RECOVERY-MATRIX-METHOD
 
!     DMAP SEQUENCES. ONLY SORT2 IS USED
 
!     (TRANSIENT RESPONCE)
!     ====================
!     DDRMM    CASEXX,UHVT,PPT,IPHIP2,IQP2,IES2,IEF2,XYCDB,EST,MPT,DIT/
!              ZUPV2,ZQP2,ZES2,ZEF2, $
 
!     (FREQUENCY RESPONSE)
!     ====================
!     DDRMM    CASEXX,UHVF,PPF,IPHIP1,IQP1,IES1,IEF1,XYCDB,EST,MPT,DIT/
!              ZUPVC1,ZQPC1,ZESC1,ZEFC1, $
!       OR
!     DDRMM    CASEXX,UHVF,PPF,IPHIP2,IQP2,IES2,IEF2,XYCDB,EST,MPT,DIT/
!              ZUPVC2,ZQPC2,ZESC2,ZEFC2, $
 
 LOGICAL :: trnsnt  ,sort2    ,col1     ,frstid
 INTEGER :: outpt   ,sysbuf   ,z        ,rd       ,wrt      ,cls
 INTEGER :: FILE    ,buf(150) ,uv       ,rdrew    ,wrtrew   ,clsrew
 INTEGER :: casecc  ,phase    ,subr(2)  ,buff     ,dhsize
 INTEGER :: pp      ,dva(3)   ,scrt1    ,scrt2    ,scrt3    ,scrt4
 INTEGER :: scrt5   ,scrt6    ,frqset   ,entrys   ,scrt     ,buf1
 INTEGER :: buf2    ,buf3     ,buf4     ,buf5     ,buf6
 INTEGER :: ifile(4),ofile(4) ,filnam   ,sets     ,passes
 INTEGER :: eor     ,outfil   ,setid    ,uvsol    ,scrt7
 INTEGER :: savdat  ,subcas   ,inblk(15),oublk(15)
 REAL :: rz(1)   ,lambda   ,ridrec(146)
 CHARACTER (LEN=27) :: swm
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON  /xmssg /  ufm      ,uwm      ,uim,sfm  ,swm
 COMMON  /system/  sysbuf   ,outpt    ,xsys(22) ,iswtch
 COMMON  /names /  rd       ,rdrew    ,wrt      ,wrtrew   ,clsrew ,cls
 COMMON  /zzzzzz/  z(1)
 COMMON  /ddrmc1/  idrec(146),buff(6) ,passes   ,outfil   ,jfile  &
     ,mcb(7)   ,entrys   ,sets(5,3),infile   ,lambda  &
     ,FILE     ,sort2    ,col1     ,frstid   ,ncore  &
     ,nsols    ,dhsize   ,filnam(2),rbuf(150),idout  &
     ,icc      ,ncc      ,ilist    ,nlist    ,nwds  &
     ,setid    ,trnsnt   ,i1       ,i2       ,phase  &
     ,itype1   ,itype2   ,nptsf    ,lsf      ,nwdsf  &
     ,scrt(7)  ,ierror   ,itemp    ,device   ,FORM  &
     ,istlst   ,lstlst   ,uvsol    ,nlambs   ,nwords ,omega    ,ipass    ,subcas
 COMMON  /stdata/  lminor   ,nstxtr   ,npos     ,savdat(110)
 EQUIVALENCE       (rz(1),z(1))   , (rbuf(1),buf(1))
 EQUIVALENCE       (scrt1,scrt(1)), (scrt2,scrt(2))
 EQUIVALENCE       (scrt3,scrt(3)), (scrt4,scrt(4))
 EQUIVALENCE       (scrt5,scrt(5)), (scrt6,scrt(6))
 EQUIVALENCE       (scrt7,scrt(7)), (ridrec(1),idrec(1))
 EQUIVALENCE       (buf1,buff(1)) , (buf2,buff(2))
 EQUIVALENCE       (buf3,buff(3)) , (buf4,buff(4))
 EQUIVALENCE       (buf5,buff(5)) , (buf6,buff(6))
 DATA     ifrout/  145 /,  dva   / 20, 32, 29 /
 DATA     istres,  iforce, ispcf / 23, 26, 35 /
 DATA     ilsym /  166 /
 DATA     subr  /  4HDDRM, 4HM        /
 DATA     eor   ,  noeor  / 1,  0        /
 DATA     casecc,  uv, pp /  101, 102, 103 /
 DATA     ifile /  104, 105, 106, 107 /
 DATA     ofile /  201, 202, 203, 204 /
 
!     DETERMINE OPEN CORE AVAILABLE AND ALLOCATE BUFFERS.
 
 DO  i = 1,100
   savdat(i) = 0
 END DO
 DO  i = 6,8
   savdat(i   ) = 102
   savdat(i+11) = 102
 END DO
 savdat(  15) = 102
 savdat(  76) = 2
 savdat(  77) = 10
 DO  i = 1,7
   scrt(i) = i + 300
 END DO
 ncore = korsz(z)
 DO  i = 1,6
   buff(i) = ncore - sysbuf - 2
   ncore   = buff(i) - 1
 END DO
 
!     GET FIRST SUBCASE OF CASE CONTROL INTO CORE
 
 ierror = 0
 subcas = 1
 icc    = 1
 FILE   = casecc
 CALL OPEN (*480,casecc,z(buf1),rdrew)
 CALL fwdrec (*490,casecc)
 CALL READ (*490,*30,casecc,z(icc),ncore-icc,noeor,nwds)
 ierror = 1
 GO TO 510
 
 30 ncc = icc + nwds - 1
 CALL CLOSE (casecc,cls)
 
!     READ TRAILER OF SOLUTION DATA BLOCK. IF SOLUTION IS
!     COMPLEX, THEN FREQUENCY RESPONCE IS ASSUMED. IF REAL, THEN
!     TRANSIENT RESPONSE IS ASSUMED.
 
 mcb(1) = uv
 CALL rdtrl (mcb)
 trnsnt = .true.
 IF (mcb(5) > 2) trnsnt = .false.
 
!     SET NUMBER OF EIGENVALUES = ROWS IN SOLUTION DATA BLOCK
 
 nlambs = mcb(3)
 
!     SET NUMBER OF SOLUTIONS.(TIME STEPS X 3, OR FREQUENCYS)
 
 nsols = mcb(2)
 
!     OPEN UV AND POSITION OVER HEADER RECORD.
 
 FILE = uv
 CALL OPEN (*480,uv,z(buf1),rdrew)
 CALL fwdrec (*490,uv)
 CALL CLOSE (uv,cls)
 
!     READ LIST OF FREQUENCYS OR TIME STEPS FROM INPUT LOAD MATRIX
!     HEADER.
 
 33 ilist = ncc + 1
 FILE  = pp
 CALL OPEN (*480,pp,z(buf1),rdrew)
 ierror = 2
 CALL READ (*490,*500,pp,buf(1),-2,noeor,nwds)
 CALL READ (*490,*35,pp,z(ilist),ncore-ilist,noeor,entrys)
 GO TO 510
 
 35 nlist = ilist + entrys - 1
 CALL CLOSE (pp,clsrew)
 
!     IF FREQUENCY RESPONSE PROBLEM, AND USER HAS SPECIFIED A LIST OF
!     FREQUENCYS TO BE USED AS A GUIDE IN DETERMINING A SUBSET OF
!     SOLUTIONS FOR OUTPUT PURPOSES, AND NOT ALL SOLUTIONS WILL BE
!     OUTPUT, THEN A MODIFIED SOLUTION MATRIX IS NOW FORMED ON
!     SCRATCH-1. THIS WILL ELIMINATE UNNECESSARY MATRIX-MULTIPLIES LATER
 
!     IN ANY EVENT THE NEXT SUBCASE-S SOLUTIONS ARE PLACED ON SCRT1.
 
 uvsol = uv
 IF (trnsnt) GO TO 190
 
!     EXPAND LIST OF FREQS PLACING A FLAG AFTER EACH.
 
 j = nlist
 nlist  = nlist + entrys
 ierror = 4
 IF (nlist > ncore) GO TO 510
 k = nlist - 1
 DO  i = 1,entrys
   z(k  ) = z(j)
   z(k+1) = 0
   k = k - 2
   j = j - 1
 END DO
 
!     SET FLAGS OF FREQUENCYS TO BE OUTPUT.
 
 INDEX  = icc + ifrout - 1
 frqset = z(INDEX)
 IF (frqset <= 0) GO TO 60
 INDEX = icc + ilsym - 1
 INDEX = z(INDEX) + 1
 50 isetx = INDEX + 2
 nsetx = isetx + z(INDEX+1) - 1
 IF (z(INDEX) == frqset) GO TO 80
 INDEX = nsetx + 1
 IF (INDEX < ncc) GO TO 50
 frqset = -1
 
!     ALL FREQUENCYS TO BE OUTPUT.
 
 60 DO  i = ilist,nlist,2
   z(i+1) = 1
 END DO
 GO TO 110
 
!     COMPARE REQUESTED FREQS WITH ACTUAL FREQS.
 
 80 DO  i = isetx,nsetx
   k    = 0
   diff = 1.0E+25
   frq  = rz(i)
   DO  j = ilist,nlist,2
     IF (z(j+1) /= 0) CYCLE
     diff1 = ABS(rz(j)-frq)
     IF (diff1 >= diff) CYCLE
     diff = diff1
     k    = j
   END DO
   IF (k /= 0) z(k+1) = 1
 END DO
 
 110 FILE = uv
 ierror = 5
 CALL OPEN (*480,uv,z(buf1),rd)
 FILE = scrt1
 CALL OPEN (*480,scrt1,z(buf2),wrtrew)
 CALL fname (scrt1,filnam)
 CALL WRITE (scrt1,filnam,2,eor)
 FILE = uv
 
!     COPY SOLUTION COLUMNS TO BE USED BY NOTEING FREQS MARKED FOR USE.
 
 nsols    = 0
 inblk(1) = uv
 oublk(1) = scrt1
 DO  i = ilist,nlist,2
   IF (z(i+1) == 0.0) THEN
     GO TO   120
   ELSE
     GO TO   130
   END IF
   120 CALL fwdrec (*490,uv)
   CYCLE
   
!     BLAST COPY THIS SOLUTION.
   
   130 icol = (i-ilist)/2 + 1
   CALL cpystr (inblk,oublk,0,icol)
   nsols = nsols + 1
 END DO
 
!     RESET -UV- DATA BLOCK DESIGNATOR TO POINT TO SCRT1, AND WRITE
!     A TRAILER.
 
 CALL CLOSE (uv,cls)
 CALL CLOSE (scrt1,clsrew)
 mcb(1) = uv
 CALL rdtrl (mcb)
 mcb(1) = scrt1
 mcb(2) = nsols
 CALL wrttrl (mcb)
 uvsol = scrt1
 
!     SHRINK UP THE FREQUENCY LIST TO MATCH SOLUTION MATRIX
 
 j = ilist - 1
 DO   i = ilist,nlist,2
   IF (z(i+1) == 0.0) THEN
     GO TO   180
   END IF
   170 j = j + 1
   z(j) = z(i)
   180 CONTINUE
 END DO
 nlist = j
 
!     IF THIS IS A TRANSIENT RESPONSE PROBLEM, THE SOLUTION MATRIX IS
!     NOW PARTITIONED INTO 3 SOLUTION MATRICES FOR DISP, VEL, AND ACCEL.
 
 190 IF (.NOT. trnsnt) GO TO 260
 FILE   = uv
 ierror = 6
 CALL OPEN (*480,uv,z(buf1),rd)
 mcb(1)   = uv
 inblk(1) = uv
 CALL rdtrl (mcb)
 DO   i = 1,3
   FILE = scrt(i)
   ibuf = buff(i+1)
   CALL OPEN (*480,FILE,z(ibuf),wrtrew)
   CALL fname (FILE,filnam)
   CALL WRITE (FILE,filnam,2,eor)
   mcb(1) = FILE
   mcb(2) = nsols/3
   CALL wrttrl (mcb)
 END DO
 ierror = 7
 FILE   = uv
 DO  i = 1,nsols,3
   DO  j = 1,3
     oublk(1) = scrt(j)
     CALL cpystr (inblk,oublk,0,i)
   END DO
 END DO
 CALL CLOSE (uv,clsrew)
 nsols = nsols/3
 
 DO  i = 1,3
   CALL CLOSE (scrt(i),clsrew)
 END DO
 
!     SDR2 FORMED MODAL SOLUTIONS FOR DISPLACEMENTS, SINGLE-POINT-
!     CONSTRAINT-FORCES, ELEMENT STRESSES, AND ELEMENT FORCES MAY BE
!     PRESENT. (ALL WILL BE SORT1-REAL, OR SORT2-REAL)
 
!     IF THIS IS A TRANSIENT PROBLEM, THE SOLUTIONS PRESENT HAVE BEEN
!     PARTITIONED INTO THE DISPLACEMENT, VELOCITY, AND ACCELERATION
!     SUBSETS. ONLY WHEN OPERATING ON THE MODAL DISPLACEMENTS WILL THE
!     VELOCITY AND ACCELERATION SOLUTION SUBSET MATRICES BE USED.
 
 260 jfile  = 1
 270 infile = ifile(jfile)
 
!     CHECK FOR EXISTENCE OF MODAL SOLUTION -INFILE-.
 
 CALL OPEN (*470,infile,z(buf1),rdrew)
 CALL fwdrec (*460,infile)
 
!     INFILE DOES EXIST.SET PARAMETERS FOR PROCESSING
 
 
!     OPEN OFP-FORMAT OUTPUT FILE FOR THIS INFILE.
 
 outfil = ofile(jfile)
 iwrt = wrtrew
 IF (subcas > 1) iwrt = wrt
 CALL OPEN (*280,outfil,z(buf4),iwrt)
 IF (subcas > 1) GO TO 305
 GO TO 300
 280 WRITE  (outpt,290) uwm,infile
 290 FORMAT (a25,' 2331. (DDRMM-2) OUTPUT DATA BLOCK CORRESPONDING TO',  &
     ' INPUT MODAL SOLUTION DATA BLOCK',i4, /5X,  &
     'IS NOT PRESENT.  INPUT DATA BLOCK IGNORED.')
 GO TO 460
 
 300 CALL fname (outfil,filnam)
 CALL WRITE (outfil,filnam,2,eor)
 305 CALL CLOSE (outfil,cls)
 
!     READ FIRST OFP-ID RECORD AND DETERMINE WHAT THE HELL IS REALLY
!     PRESENT.
 
 ierror = 14
 CALL READ (*460,*460,infile,idrec,146,eor,nwds)
 
!     MAJOR ID AND SORT1 OR SORT2 DETERMINATION.
 
 itype1 = idrec(2)/1000
 sort2  = .false.
 IF (itype1 > 1) sort2 = .true.
 itype1 = idrec(2) - itype1*1000
 
!     BRANCH ON MAJOR ID
 
 IF (itype1 < 1 .OR. itype1 > 7) GO TO 410
 passes = 1
 SELECT CASE ( itype1 )
   CASE (    1)
     GO TO 410
   CASE (    2)
     GO TO 410
   CASE (    3)
     GO TO 370
   CASE (    4)
     GO TO 380
   CASE (    5)
     GO TO 390
   CASE (    6)
     GO TO 410
   CASE (    7)
     GO TO 310
 END SELECT
 
!     MODAL DISPLACEMENTS = EIGENVECTORS ARE ON INFILE.
 
 310 passes = 3
 nwords = 2
 
!     DETERMINE DISP, VEL, AND ACCEL SET REQUESTS.
 
 312 ibase = icc + ilsym - 1
 ibase = z(ibase) + 1
 DO  i = 1,passes
   IF (passes == 1) GO TO 315
   itemp = icc + dva(i) - 1
   315 sets(1,i) = z(itemp)
   sets(2,i) = z(itemp+1)
   sets(3,i) = IABS(z(itemp+2))
   sets(4,i) = 0
   sets(5,i) = 0
   IF (sets(1,i) > 0.0) THEN
     GO TO   320
   ELSE
     GO TO   350
   END IF
   320 INDEX = ibase
   330 isetx = INDEX + 2
   IF (z(INDEX) == sets(1,i)) GO TO 340
   INDEX = isetx + z(INDEX+1)
   IF (INDEX < ncc) GO TO 330
   sets(1,i) = -1
   CYCLE
   340 sets(4,i) = isetx
   sets(5,i) = z(INDEX+1)
   350 CONTINUE
 END DO
 GO TO 430
 
!     MODAL SPCF-S ARE ON INFILE.
 
 370 itemp  = icc + ispcf - 1
 nwords = 2
 GO TO 312
 
!     MODAL FORCES ARE ON INFILE.
 
 380 itemp  = icc + iforce - 1
 nwords = 1
 GO TO 312
 
!     MODAL STRESSES ARE ON INFILE.
 
 390 itemp  = icc + istres - 1
 nwords = 1
 GO TO 312
 
!     ILLEGAL INFILE DATA.
 
 410 WRITE  (outpt,420) uwm,infile
 420 FORMAT (a25,' 2332.  (DDRMM-4) INVALID INPUT DATA DETECTED IN ',  &
     'DATA BLOCK',i5,'. PROCESSING STOPPED FOR THIS DATA BLOCK')
 GO TO 460
 
!     CALL PROCESSOR TO BUILD DATA-MATRIX ON SCRT5 AND MAPPING-DATA ON
!     SCRT4, AND THEN PERFORM OUTPUT OF RESULTS TO OUTFIL.
 
 430 IF (sort2) GO TO 440
 
!     SORT1 PROCESSOR
 
 CALL ddrmm1 (*480,*490,*500,*510)
 GO TO 450
 
!     SORT2 PROCESSOR
 
 440 CALL ddrmm2 (*480,*490,*500,*510)
 
!     WRAP UP PROCESSING FOR THIS INFILE.
 
 450 mcb(1) = outfil
 mcb(2) = 1
 CALL wrttrl (mcb)
 460 CALL CLOSE (outfil,clsrew)
 CALL CLOSE (infile,clsrew)
 CALL CLOSE (scrt4, clsrew)
 CALL CLOSE (scrt5, clsrew)
 
!     PROCESS NEXT MODAL SOLUTION INPUT.
 
 470 jfile = jfile + 1
 IF (jfile <= 4) GO TO 270
 
!     ALL WORK COMPLETE FOR THIS SUBCASE. IF FREQUENCY RESPONSE PROCESS
!     NEXT SUBCASE.
 
 IF (trnsnt) GO TO 479
 FILE   = casecc
 ierror = 471
 CALL OPEN (*480,casecc,z(buf1),rd)
 CALL READ (*479,*472,casecc,z(icc),ncore-icc,noeor,nwds)
 GO TO 510
 
 472 ncc    = icc + nwds - 1
 subcas = subcas + 1
 CALL CLOSE (casecc,cls)
 GO TO 33
 
!//// SUBCAS  NUMBER NEEDS TO GET INTO OUTPUT BLOCKS
 
 479 ierror = 511
 CALL CLOSE (casecc,clsrew)
 GO TO 560
 
!     ERRORS FORCING TERMINATION OF THIS MODULE.
 
 480 kk = 1
 GO TO 520
 490 kk = 2
 GO TO 520
 500 kk = 3
 GO TO 520
 510 kk = 8
 GO TO 520
 520 CALL mesage (kk,FILE,subr)
 WRITE  (outpt,530) swm,ierror
 530 FORMAT (a27,' 2333.  (DDRMM-1) MODULE DDRMM TERMINATED WITH ',  &
     'VARIABLE IERROR =',i10)
 
!     INSURE ALL FILES CLOSED BEFORE RETURNING.
 
 DO  l = 100,300,100
   DO  m = 1,11
     jfile = m + l
     CALL CLOSE (jfile,clsrew)
   END DO
 END DO
 
!     INSURE ALL OUT-FILES HAVE AN EOF.
 
 560 DO  l = 201,204
   CALL OPEN (*570,l,z(buf1),wrt)
   CALL CLOSE (l,clsrew)
   570 CONTINUE
 END DO

 RETURN
END SUBROUTINE ddrmm
