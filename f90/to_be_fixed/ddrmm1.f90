SUBROUTINE ddrmm1 (*,*,*,*)
     
!     PERFORMS SORT1 TYPE PROCESSING FOR MODULE DDRMM.
 
 
 , INTENT(IN OUT)                         :: *
 , INTENT(IN OUT)                         :: *
 , INTENT(IN OUT)                         :: *
 , INTENT(IN OUT)                         :: *
 LOGICAL :: sort2    ,col1     ,frstid   ,idout    ,trnsnt  , anyxy    ,lminor
 INTEGER :: buf1     ,buf2     ,buf3     ,buf4     ,buf5    ,  &
     buf6     ,buff     ,eor      ,rd       ,rdrew   ,  &
     wrt      ,wrtrew   ,cls      ,clsrew   ,elem    ,  &
     ia(4)    ,sets     ,entrys   ,sysbuf   ,outpt   ,  &
     passes   ,outfil   ,FILE     ,dhsize   ,filnam  ,  &
     setid    ,FORM     ,device   ,phase    ,scrt    ,  &
     scrt1    ,scrt2    ,scrt3    ,scrt4    ,scrt5   ,  &
     scrt6    ,scrt7    ,dvamid(3),buf(150) ,z(1)    ,  &
     uvsol    ,subcas   ,savdat   ,savpos   ,bufsav
 REAL :: ridrec(6),lambda
 CHARACTER (LEN=27) :: swm
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm      ,uwm      ,uim      ,sfm      ,swm
 COMMON /stdata/ lminor   ,nstxtr   ,npos     ,savdat(75)        ,  &
     savpos(25)         ,bufsav(10)
 COMMON /system/ sysbuf   ,outpt
 COMMON /names / rd       ,rdrew    ,wrt      ,wrtrew   ,clsrew  , cls
 COMMON /zblpkx/ a(4)     ,irow
 COMMON /zntpkx/ aout(4)  ,irowo    ,ieol     , IEOR
 COMMON /gpta1 / nelem    ,last     ,incr     ,elem(1)
 COMMON /zzzzzz/ rz(1)
 COMMON /mpyadx/ mcba(7)  ,mcbb(7)  ,mcbc(7)  ,mcbd(7)  ,lz      ,  &
     itflag   ,isinab   ,isinc    ,iprec    ,iscrt
 COMMON /ddrmc1/ idrec(146),buff(6) ,passes   ,outfil   ,jfile   ,  &
     mcb(7)   ,entrys   ,sets(5,3),infile   ,lambda  ,  &
     FILE     ,sort2    ,col1     ,frstid   ,ncore   ,  &
     nsols    ,dhsize   ,filnam(2),rbuf(150),idout   ,  &
     icc      ,ncc      ,ilist    ,nlist    ,nwds    ,  &
     setid    ,trnsnt   ,i1       ,i2       ,phase   ,  &
     itype1   ,itype2   ,nptsf    ,lsf      ,nwdsf   ,  &
     scrt(7)  ,ierror   ,itemp    ,device   ,FORM    ,  &
     istlst   ,lstlst   ,uvsol    ,nlambs   ,nwords  , omega    ,ipass    ,subcas
 COMMON /condas/ pi       ,twopi
 EQUIVALENCE     (scrt1,scrt(1)), (scrt2,scrt(2)), (scrt3,scrt(3)),  &
     (scrt4,scrt(4)), (scrt5,scrt(5)), (scrt6,scrt(6)),  &
     (scrt7,scrt(7)), (buf1 ,buff(1)), (buf2 ,buff(2)),  &
     (buf3 ,buff(3)), (buf4 ,buff(4)), (buf5 ,buff(5)),  &
     (buf6 ,buff(6)), (a(1) ,  ia(1)), (z(1) ,  rz(1)),  &
     (idrec(1),ridrec(1)), (buf(1),rbuf(1))
 DATA    eor   , noeor / 1, 0 /, dvamid / 1, 10, 11 /
 
!     FORMATION OF DATA-MATRIX AND SUBSEQUENT MULTIPLY BY SOLUTION-
!     MATRIX AND ULTIMATE OUTPUT OF TRANSIENT OR FREQUENCY SOLUTIONS.
 
 ipass  = 1
 20 col1   = .true.
 frstid = .true.
 setid  = sets(1,ipass)
 device = sets(2,ipass)
 FORM   = sets(3,ipass)
 istlst = sets(4,ipass)
 lstlst = sets(5,ipass)
 
!     GET LIST OF XYPLOT REQUESTED IDS FOR CURRENT SUBCASE AND
!     OUTFIL TYPE.
 
 SELECT CASE ( jfile )
   CASE (    1)
     GO TO 22
   CASE (    2)
     GO TO 23
   CASE (    3)
     GO TO 24
   CASE (    4)
     GO TO 25
 END SELECT
 
!     DISPLACEMENT, VELOCITY, ACCELERATION
 
 22 ixytyp = ipass
 GO TO 26
 
!     SPCF
 
 23 ixytyp = 4
 GO TO 26
 
!     STRESS
 
 24 ixytyp = 6
 GO TO 26
 
!     FORCE
 
 25 ixytyp = 7
 GO TO 26
 
 26 ixy = nlist + 1
 CALL ddrmmp (*380,z(ixy),buf3-ixy,lxy,ixytyp,subcas,z(buf3),anyxy)
 IF (.NOT.anyxy .AND. setid == 0) GO TO 280
 nxy = ixy + lxy - 1
 
!     INITIALIZE DATA MATRIX FILE(SCRT5), AND MAPPING TABLE FILE(SCRT4).
 
 ierror = 22
 FILE = scrt4
 CALL OPEN (*350,scrt4,z(buf3),wrtrew)
 FILE = scrt5
 CALL OPEN (*350,scrt5,z(buf2),wrtrew)
 CALL fname (scrt5,filnam)
 CALL WRITE (scrt5,filnam,2,eor)
 
!     GENERAL LOGIC TO BUILD SORT1 FORMAT DATA MATRIX.
 
!     EACH COLUMN WRITTEN HERE REPRESENTS ONE EIGENVALUE.
 
!          COMPONENTS FOR FIRST ID   *
!              .                      *
!              .                       *
!              .                        *
!          COMPONENTS FOR NEXT ID        * ONE COLUMN
!              .                        *  OF DATA FOR EACH EIGENVALUE.
!              .                       *
!              .                      *
!             ETC                    *
!     --------------------------------------------- EOR
 
!          IDENTICAL COMPONENTS ARE REPRESENTED IN EACH COLUMN.
 
 
!     READ AN OFP-ID-RECORD AND SET PARAMETERS.
!     (ON ENTRY TO THIS PROCESSOR THE FIRST ID RECORD IS AT HAND)
 
 FILE   = infile
 mcb(1) = scrt5
 mcb(2) = 0
 mcb(3) = 0
 mcb(4) = 2
 mcb(5) = 1
 mcb(6) = 0
 mcb(7) = 0
 IF (ipass == 1 .AND. frstid) GO TO 50
 40 CALL READ (*130,*130,infile,idrec,146,eor,nwds)
 majid = MOD(idrec(2),1000)
 IF (majid /= itype1) GO TO 310
 
!     IF FIRST COLUMN, OFP-ID RECORD IS WRITTEN AS IS TO MAP FILE.
 
 50 IF (.NOT. col1) GO TO 60
 IF (.NOT.frstid .AND. ridrec(6) /= lambda) GO TO 60
 CALL WRITE (scrt4,idrec,146,eor)
 60 lentry = idrec(10)
 i1 = nwords + 1
 i2 = lentry
 minor = idrec(3)
 
!     IF SAME EIGENVALUE AS THAT OF LAST OFP-ID RECORD THEN CONTINUE.
 
 IF (frstid) GO TO 70
 IF (ridrec(6) == lambda) GO TO 80
 
!     NEW EIGENVALUE. COMPLETE CURRENT DATA MATRIX COLUMN AND START
!     NEW COLUMN. PASS ONE IS NOW COMPLETE.
 
 CALL bldpkn (scrt5,0,mcb)
 IF (col1) irow1 = irow
 IF (irow /= irow1) GO TO 290
 col1 = .false.
 
!     START NEW COLUMN.
 
 70 CALL bldpk (1,1,scrt5,0,0)
 irow   = 0
 frstid = .false.
 lambda = ridrec(6)
 
!     READ A POINT OR ELEMENT ENTRY.
 
 80 CALL READ (*360,*120,infile,buf,lentry,noeor,nwds)
 id = buf(1)/10
 
!     CHECK FOR ID IN OUTPUT REQUEST LIST
 
 idvice = device
 IF (setid < 0.0) THEN
   GO TO   100
 ELSE IF (setid == 0.0) THEN
   GO TO    95
 END IF
 
!//// NEXT MAY NOT NEED TO BE INITIALIZED EVERY TIME.
 
 90 next = 1
 CALL setfnd (*95,z(istlst),lstlst,id,next)
 GO TO 100
 95 IF (.NOT. anyxy) GO TO 80
 CALL bisloc (*80,id,z(ixy),1,lxy,jp)
 idvice = 0
 
!     THIS ID IS TO BE OUTPUT.
 
 100 IF (.NOT. col1) GO TO 105
 buf(1) = 10*id + idvice
 CALL WRITE (scrt4,buf(1),nwords,noeor)
 nstxtr = 0
 IF (itype1 /= 5 .OR. savdat(minor) == 0) GO TO 104
 npos   = savdat(minor)/100
 nstxtr = savdat(minor) - npos*100
 DO  i = 1,nstxtr
   j = savpos(npos+i-1)
   bufsav(i) = buf(j)
 END DO
 CALL WRITE (scrt4,bufsav(1),nstxtr,noeor)
 104 CONTINUE
 
!     OUTPUT TO DATA MATRIX THE COMPONENTS OF THIS ENTRY.
 
 105 DO  i = i1,i2
   irow = irow + 1
   a(1) = rbuf(i)
   
!     GET RID OF INTEGERS.
   
!     OLD LOGIC -
!     IF (MACH.NE.5 .AND.  IABS(IA(1)) .LT.   100000000) A(1) = 0.0
!     IF (MACH.EQ.5 .AND. (IA(1).LE.127.AND.IA(1).GE.1)) A(1) = 0.0
!     OLD LOGIC SHOULD INCLUDE ALPHA MACHINE (MACH=21)
   
!     NEW LOGIC BY G.CHAN/UNISYS, 8/91 -
   IF (numtyp(ia(1)) <= 1) a(1) = 0.0
   
   CALL zblpki
 END DO
 GO TO 80
 
!     END OF CURRENT OFP-DATA RECORD ENCOUNTERED.
!     IF NEXT OFP-ID-RECORD INDICATES ANOTHER OFP-DATA RECORD FOR
!     THIS SAME EIGENVALUE (I.E. A CHANGE IN ELEMENT TYPE) THEN
!     FURTHER CONSTRUCTION OF THE DATA MATRIX COLUMN TAKES PLACE.
 
 120 IF (col1) CALL WRITE (scrt4,0,0,eor)
 GO TO 40
 
!     END OF FILE ENCOUNTERED ON INFILE.
!     DATA MATRIX AND MAPING FILE ARE COMPLETE.
 
 130 CALL CLOSE (infile,clsrew)
 CALL CLOSE (scrt4 ,clsrew)
 
!     COMPLETE LAST COLUMN OF DATA MATRIX WRITTEN.
 
 IF (col1) irow1 = irow
 IF (irow /= irow1) GO TO 310
 CALL bldpkn (scrt5,0,mcb)
 mcb(3) = irow
 CALL wrttrl (mcb)
 CALL CLOSE (scrt5,clsrew)
 
!     TO GET SOLUTION MATRIX BASED ON SORT-1 INFILE.
 
!     SOLVE,
!              (DATA MATRIX)     X    (MODAL SOLUTION MATRIX)
!             NCOMPS X NLAMBS           NLAMBS X NSOLUTIONS
!             ===============         =======================
 
!     RESULTANT MATRIX IS NCOMPS BY NSOLUTIONS IN SIZE.
 
 
!     MATRIX MULTIPLY SETUP AND CALL.
 
 mcba(1) = scrt5
 CALL rdtrl (mcba)
 mcbb(1) = uvsol
 IF (trnsnt) mcbb(1) = scrt(ipass)
 CALL rdtrl (mcbb)
 mcbc(1) = 0
 mcbd(1) = scrt6
 mcbd(2) = 0
 mcbd(3) = irow
 mcbd(4) = 2
 mcbd(5) = 1
 mcbd(6) = 0
 mcbd(7) = 0
 IF (.NOT.trnsnt) mcbd(5) = 3
 nxy1    = nxy + 1
 IF (MOD(nxy1,2) == 0) nxy1 = nxy1 + 1
 lz      = korsz(z(nxy1))
 itflag  = 0
 isinab  = 1
 isinc   = 1
 iprec   = 1
 iscrt   = scrt7
 CALL mpyad (z(nxy1),z(nxy1),z(nxy1))
 mcbd(1) = scrt6
 CALL wrttrl (mcbd)
 
!     PRODUCT MATRIX IS NOW OUTPUT, USING THE MAP ON SCRT4 FOR EACH
!     COLUMN.  (SORT-1)  PRODUCT MATRIX IS ON SCRATCH DATA BLOCK 6.
 
 ierror = 10
 FILE = outfil
 CALL OPEN (*350,outfil,z(buf1),wrt)
 FILE = scrt4
 CALL OPEN (*350,scrt4,z(buf2),rdrew)
 FILE = scrt6
 CALL OPEN (*350,scrt6,z(buf3),rdrew)
 CALL fwdrec (*360,scrt6)
 jlist = ilist
 
!     LOOP ON COLUMNS OF SCRT6.
 
 140 CALL ddrmma (.true.)
 
!     READ AN OFP-ID-RECORD FROM THE MAP.
 
 FILE = scrt4
 150 CALL READ (*270,*370,scrt4,idrec,146,eor,nwds)
 
!     SET THE FREQUENCY OR TIME AND CLOBBER THE EIGENVALUE.
 
 ridrec(5) = rz(jlist)
 ridrec(6) = 0.0
 idout = .false.
 minor = idrec(3)
 
!     SET NUMBER OF STRESS OR FORCE WORDS AND COMPLEX POINTERS IF
!     NECESSARY
 
 itype2 = idrec(3)
 IF (itype1 == 3 .OR. itype1 == 7) GO TO 200
 ielem = (itype2-1)*incr
 IF (itype1 == 4) GO TO 180
 IF (itype1 == 5) GO TO 190
 WRITE  (outpt,170) swm,itype1,itype2,infile
 170 FORMAT (a27,' 2334.  (DDRMM-3) ILLEGAL MAJOR OR MINOR OFP-ID ',  &
     'IDENTIFICATIONS =',2I10, /5X,'DETECTED IN DATA BLOCK',i5,  &
     '. PROCESSING OF SAID DATA BLOCK DISCONTINUED.')
 GO TO 340
 
!     FORCES ASSUMED.
 
 180 lsf   = elem(ielem+19)
 nptsf = elem(ielem+21)
 nwdsf = lsf
 GO TO 220
 
!     STRESSES ASSUMED.
 
 190 lsf   = elem(ielem+18)
 nptsf = elem(ielem+20)
 nwdsf = lsf
 GO TO 220
 
!     SPCF OR DISPLACEMENTS ASSUMED
 
 200 IF (.NOT.trnsnt) GO TO 210
 nwdsf = 8
 GO TO 220
 210 nwdsf = 14
 
!     SET OMEGA IF THIS IS THE VELOCITY OR ACCELERATION PASS
 
 SELECT CASE ( ipass )
   CASE (    1)
     GO TO 220
   CASE (    2)
     GO TO 211
   CASE (    3)
     GO TO 212
 END SELECT
 
!     OMEGA FOR VELOCITY PASS
 
 211 omega = twopi*rz(jlist)
 GO TO 220
 
!     OMEGA FOR ACCELERATION PASS
 
 212 omega = -((twopi*rz(jlist))**2)
 
 220 lentry = idrec(10)
 i1 = nwords + 1
 i2 = lentry
 
 
!     SET DISPLACEMENT, VELOCITY, OR ACCELERATION OFP MAJOR ID IF INFILE
!     IS MODAL DISPLACEMENTS.
 
 IF (itype1 /= 7) GO TO 230
 idrec(2) = dvamid(ipass)
 230 IF (.NOT.trnsnt) idrec(2) = idrec(2) + 1000
 
!     RESET APPROACH CODE FROM EIGENVALUE TO TRANSIENT OR FREQUENCY
 
 iapp = 5
 IF (trnsnt) iapp = 6
 idrec(1) = 10*iapp + device
 
!     FILL TITLE, SUBTITLE, AND LABEL FROM CASECC FOR THIS SUBCASE.
 
 DO  i = 1,96
   idrec(i+50) = z(icc+i+37)
 END DO
 idrec(4) = subcas
 
!     READ FIRST WORDS OF OUTPUT ENTRY FROM MAP.
 
 240 CALL READ (*360,*260,scrt4,buf,nwords,noeor,nwds)
 lminor = .true.
 IF (itype1 /= 5 .OR. savdat(minor) == 0) GO TO 241
 npos   = savdat(minor)/100
 nstxtr = savdat(minor) - npos*100
 CALL READ (*360,*370,scrt4,bufsav(1),nstxtr,noeor,nwds)
 lminor = .false.
 241 CONTINUE
 
!     GET BALANCE USING UTILITY WHICH WILL COLLECT AND MAP TOGETHER
!     AS REQUIRED REAL OR COMPLEX, AND GENERATE MAGNITUDE/PHASE IF
!     REQUIRED.  (THIS ROUTINE WILL BUFFER DATA IN FROM SCRT6 AS IT
!     NEEDS IT.)
 
 CALL ddrmma (.false.)
 
!     CALL DDRMMS TO RECOMPUTE SOME ELEMENT STRESS QUANTITIES
!     IN TRANSIENT PROBLEMS ONLY.
 
 IF (trnsnt .AND. itype1 == 5) CALL ddrmms (buf,itype2,buf4,buf5)
 IF (idout) GO TO 250
 idrec( 9) = FORM
 idrec(10) = nwdsf
 CALL WRITE (outfil,idrec,146,eor)
 idout = .true.
 
!     OUTPUT THE COMPLETED ENTRY TO OFP OUTFIL.
 
 250 CALL WRITE (outfil,buf,nwdsf,noeor)
 GO TO 240
 
!     END OF ENTRIES FOR ONE ID-REC HIT.  IF NO EOF ON MAP WITH
!     NEXT READ, THEN CONTINUE OUTPUT OF THIS SOLUTION COLUMN.
 
 260 CALL WRITE (outfil,0,0,eor)
 GO TO 150
 
!     END OF FILE ON MAP.  THUS START NEXT COLUMN IF REQUIRED.
 
 270 jlist = jlist + 1
 IF (jlist > nlist) GO TO 280
 CALL REWIND (scrt4)
 GO TO 140
 
!     ALL DATA OF SOLUTION PRODUCT MATRIX HAS NOW BEEN OUTPUT.
 
 280 CALL CLOSE (outfil,cls)
 CALL CLOSE (infile,clsrew)
 CALL CLOSE (scrt4,clsrew)
 CALL CLOSE (scrt6,clsrew)
 ipass = ipass + 1
 IF (ipass > passes) GO TO 340
 
!     PREPARE FOR ANOTHER PASS
 
 FILE = infile
 CALL OPEN (*350,infile,z(buf1),rdrew)
 CALL fwdrec (*360,infile)
 GO TO 20
 
!     DATA INCONSISTENCY ON -INFILE-.
 
 290 WRITE  (outpt,300) swm,infile
 300 FORMAT (a27,' 2335.  (DDRMM1-1) THE AMOUNT OF DATA IS NOT ',  &
     'CONSISTENT FOR EACH EIGENVALUE IN DATA BLOCK',i5, /5X,  &
     'PROCESSING OF THIS DATA BLOCK TERMINATED.')
 GO TO 330
 
!     CHANGE IN MAJOR OFP-ID DETECTED ON -INFILE-.
 
 310 WRITE  (outpt,320) swm,infile
 320 FORMAT (a27,' 2336.  (DDRMM1-2) A CHANGE IN WORD 2 OF THE OFP-ID',  &
     ' RECORDS OF DATA BLOCK',i5, /5X,'HAS BEEN DETECTED. ',  &
     ' POOCESSING OF THIS DATA BLOCK HAS BEEN TERMINATED.')
 330 ipass = 3
 GO TO 280
 
!     COMPLETION OF PASS FOR INPUT MODAL SOLUTION -FILE-.
 
 340 RETURN
 
!     UNDEFINED FILE.
 
 350 RETURN 1
 
!     END OF FILE HIT.
 
 360 RETURN 2
 
!     END OF RECORD HIT.
 
 370 RETURN 3
 
!     INSUFFICIENT CORE.
 
 380 RETURN 4
END SUBROUTINE ddrmm1
