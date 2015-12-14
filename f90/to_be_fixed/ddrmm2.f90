SUBROUTINE ddrmm2 (*,*,*,*)
     
!     PERFORMS SORT2 TYPE PROCESSING FOR MODULE DDRMM.
 
 
 , INTENT(IN)                             :: *
 , INTENT(IN)                             :: *
 , INTENT(IN)                             :: *
 , INTENT(IN)                             :: *
 LOGICAL :: sort2    ,col1     ,frstid   ,idout    ,trnsnt  , lminor   ,anyxy
 INTEGER :: buf1     ,buf2     ,buf3     ,buf4     ,buf5    ,  &
     buf6     ,buff     ,eor      ,rd       ,rdrew   ,  &
     wrt      ,wrtrew   ,cls      ,clsrew   ,elem    ,  &
     ia(4)    ,sets     ,entrys   ,sysbuf   ,outpt   ,  &
     passes   ,outfil   ,FILE     ,dhsize   ,filnam  ,  &
     setid    ,FORM     ,device   ,phase    ,scrt    ,  &
     scrt1    ,scrt2    ,scrt3    ,scrt4    ,scrt5   ,  &
     scrt6    ,scrt7    ,typout   ,dvamid(3),buf(150),  &
     z(1)     ,uvsol    ,bufa(75) ,bufb(75) ,complx  ,  &
     subcas   ,savdat   ,savpos   ,bufsav  ,elwork(300)
 REAL :: ridrec(1),lambda
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
 COMMON /zntpkx/ aout(4)  ,irowo    ,ieol     ,IEOR
 COMMON /gpta1 / nelem    ,last     ,incr     ,elem(1)
 COMMON /zzzzzz/ rz(1)
 COMMON /mpyadx/ mcba(7)  ,mcbb(7)  ,mcbc(7)  ,mcbd(7)  ,lz      ,  &
     itflag   ,isinab   ,isinc    ,iprec    ,iscrt
 COMMON /clstrs/ complx(1)
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
 EQUIVALENCE    (scrt1,scrt(1)), (scrt2,scrt(2)), (scrt3,scrt(3)),  &
     (scrt4,scrt(4)), (scrt5,scrt(5)), (scrt6,scrt(6)),  &
     (scrt7,scrt(7)), (buf1 ,buff(1)), (buf2 ,buff(2)),  &
     (buf3 ,buff(3)), (buf4 ,buff(4)), (buf5 ,buff(5)),  &
     (buf6 ,buff(6)), (a(1) ,  ia(1)), (z(1) ,  rz(1)),  &
     (buf(1),rbuf(1),bufa(1)), (bufb(1),buf(76)), (idrec(1),ridrec(1))
 
 DATA    eor  , noeor / 1, 0 /, dvamid / 2001, 2010, 2011 /
 
!     FORMATION OF DATA-MATRIX AND SUBSEQUENT MULTIPLICATION BY SAME OF
!     THE SOLUTION MATRIX (TRNNSPOSED), AND ULTIMATE OUTPUT OF TRANSIENT
!     OR FREQUENCY SOLUTIONS.
 
 ipass  = 1
 iomega = nlist  + 1
 nomega = iomega - 1
 minor  = 0
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
 
 26 ixy  = nomega + 1
 CALL ddrmmp (*480,z(ixy),buf3-ixy,lxy,ixytyp,subcas,z(buf3),anyxy)
 IF (.NOT.anyxy .AND. setid == 0) GO TO 400
 nxy  = ixy + lxy - 1
 ierror = 23
 FILE = scrt4
 CALL OPEN (*450,scrt4,z(buf3),wrtrew)
 FILE = scrt5
 CALL OPEN (*450,scrt5,z(buf2),wrtrew)
 CALL fname (scrt5,filnam)
 CALL WRITE (scrt5,filnam,2,eor)
 
!     LOGIC TO BUILD SORT-2 FORMAT DATA MATRIX.
 
!     EACH COLUMN WRITTEN HERE ENCOMPASSES ALL EIGENVALUES FOR
!     ONE COMPONENT OF ONE ID.  THE NUMBER OF COLUMNS THUS EQUALS
!     THE SUM OF ALL COMPONENTS OF ALL REQUESTED ID-S.
 
!     READ AN OFP-ID RECORD AND SET PARAMETERS.
!     (ON ENTRY TO THIS PROCESSOR ONE ID-RECORD IS AT HAND)
 
 FILE   = infile
 ierror = 19
 mcb(1) = scrt5
 mcb(2) = 0
 mcb(3) = nlambs
 mcb(4) = 2
 mcb(5) = 1
 mcb(6) = 0
 mcb(7) = 0
 IF (ipass == 1 .AND. frstid) GO TO 50
 40 CALL READ (*160,*160,infile,idrec,146,eor,nwds)
 
!     OFP-ID RECORD IS WRITTEN TO THE MAP FILE ONLY ON CHANGE OF
!     MINOR ID.
 
 50 major = MOD(idrec(2),1000)
 IF (major /= itype1) GO TO 420
 idvice = device
 id   = idrec(5)/10
 IF (setid < 0.0) THEN
   GO TO    80
 ELSE IF (setid == 0.0) THEN
   GO TO    65
 END IF
 60 next = 1
 CALL setfnd (*65,z(istlst),lstlst,id,next)
 GO TO 80
 65 IF (.NOT.anyxy) GO TO 70
 CALL bisloc (*70,id,z(ixy),1,lxy,jp)
 idvice = 0
 GO TO 80
 
!     ID IS NOT TO BE OUTPUT THUS SKIP UPCOMING OFP-DATA-RECORD.
 
 70 CALL fwdrec (*460,infile)
 GO TO 40
 
!     ID IS TO BE OUTPUT THUS CONTINUE.
 
 80 numwds = nlambs*idrec(10)
 idata  = nxy + 1
 ndata  = idata + numwds - 1
 IF (ndata < buf3) GO TO 100
 
!     INSUFFICIENT CORE
 
 insuf = ndata - buf3
 WRITE  (outpt,90) uwm,infile,insuf
 90 FORMAT (a25,' 2337.  (DDRMM2-2)  DATA BLOCK',i5,' CAN NOT BE ',  &
     'PROCESSED DUE TO', /5X,'A CORE INSUFFICIENCY OF APPROXI',  &
     'MATELY',i11,' DECIMAL WORDS.')
 GO TO 440
 100 IF (.NOT.frstid) GO TO 110
 
!     VERY FIRST ID RECORD,  THUS SET MINOR ID.
 
 frstid = .false.
 GO TO 120
 110 IF (idrec(3) == minor) GO TO 130
 
!     CHANGE IN MINOR ID, I.E. NEW ELEMENT TYPE.  COMPLETE CURRENT
!     RECORD OF MAP AND OUTPUT ANOTHER ID-RECORD.
 
 CALL WRITE (scrt4,0,0,eor)
 120 CALL WRITE (scrt4,idrec,146,eor)
 minor = idrec(3)
 
!     SAME TYPE OF DATA THUS CONTINUE ON.
 
 130 lentry = idrec(10)
 i1 = nwords + 1
 i2 = lentry
 
!     READ AND OUTPUT ONE FULL OFP-DATA RECORD.
 
 CALL READ (*460,*470,infile,z(idata),numwds,eor,nwds)
 DO  i = i1,i2
   
!     START NEW COLUMN
   
   CALL bldpk (1,1,scrt5,0,0)
   irow  = 0
   jdata = idata + i - 1
   kdata = ndata - lentry + i
   DO  j = jdata,kdata,lentry
     irow = irow + 1
     a(1) = rz(j)
     
!     ELIMINATE INTEGERS
     
!     OLD LOGIC -
!     IF (MACH.NE.5 .AND. IABS(IA(1)).LT.100000000) A(1) = 0.0
!     IF (MACH.EQ.5 .AND. (IA(1).LE.127 .AND. IA(1).GE.1)) A(1) = 0.0
!     OLD LOGIC SHOULD INCLUDE ALPHA MACHINE (MACH=21)
     
!     NEW LOGIC, BY G.CHAN/UNISYS  8/91 -
     IF (numtyp(ia(1)) <= 1) a(1) = 0.0
     
     CALL zblpki
   END DO
   
!     COMPLETE COLUMN
   
   CALL bldpkn (scrt5,0,mcb)
 END DO
 
!     OUTPUT TO MAP THE ID PLUS ANY OTHER DATA NECESSARY.
 
 buf(1) = 10*id + idvice
 IF (nwords == 2) buf(2) = z(idata+1)
 nstxtr = 0
 IF (itype1 /= 5 .OR. savdat(minor) == 0) GO TO 155
 npos   = savdat(minor)/100
 nstxtr = savdat(minor) - npos * 100
 DO  i = 1,nstxtr
   j = savpos(npos+i-1)
   buf(i+1) = z(idata+j-1)
 END DO
 155 CALL WRITE (scrt4,buf,nwords+nstxtr,noeor)
 
!     GO FOR NEXT ID.
 
 GO TO 40
 
!     END OF FILE ON INFILE.  MAP AND DATA MATRIX NOW COMPLETE.
 
 160 CALL wrttrl (mcb)
 CALL CLOSE (scrt5,clsrew)
 CALL CLOSE (infile,clsrew)
 CALL WRITE (scrt4,0,0,eor)
 CALL CLOSE (scrt4,clsrew)
 
!     SOLUTION MATRIX MAY BE FOUND BASED ON SORT-2 INFILE.
 
!     SOLVE,
!                               T
!        (MODAL SOLUTION MATRIX)     X      (DATA MATRIX)
!          NLAMBS X NSOLUTIONS             NLAMBS X NCOMPS
!        =======================           ===============
 
!     RESULTANT MATRIX IS NSOLUTIONS BY NCOMPS IN SIZE.
 
 
!     MATRIX MULTIPLY SETUP AND CALL.
 
 mcba(1) = uvsol
 IF (trnsnt) mcba(1) = scrt(ipass)
 CALL rdtrl (mcba)
 mcbb(1) = scrt5
 CALL rdtrl (mcbb)
 mcbc(1) = 0
 mcbd(1) = scrt6
 mcbd(2) = 0
 mcbd(3) = nsols
 mcbd(4) = 2
 mcbd(5) = 1
 mcbd(6) = 0
 mcbd(7) = 0
 IF (.NOT.trnsnt) mcbd(5) = 3
 itflag  = 1
 nxy1    = nxy + 1
 IF (MOD(nxy1,2) == 0) nxy1 = nxy1 + 1
 lz      = korsz(z(nxy1))
 isinab  = 1
 isinc   = 1
 iprec   = 1
 iscrt   = scrt7
 CALL mpyad (z(nxy1),z(nxy1),z(nxy1))
 mcbd(1) = scrt6
 CALL wrttrl (mcbd)
 
!     PRODUCT MATRIX IS NOW OUTPUT USING THE MAP ON SCRT4.
!     EACH COLUMN OF SCRT6 CONTAINS ALL THE TIME OR FREQUENCY STEP
!     VALUES FOR ONE COMPONENT OF ONE ID.
 
!     THUS A NUMBER OF COLUMNS ENCOMPASSING THE COMPONENTS OF ONE ID
!     MUST FIT IN CORE.
 
 ierror = 20
 FILE = outfil
 CALL OPEN (*450,outfil,z(buf1),wrt)
 FILE = scrt4
 CALL OPEN (*450,scrt4,z(buf2),rdrew)
 FILE = scrt6
 CALL OPEN (*450,scrt6,z(buf3),rdrew)
 CALL fwdrec (*460,scrt6)
 
!     READ AN OFP-ID-RECORD FROM THE MAP, AND ALLOCATE SPACE NEEDED
!     FOR SOLUTION DATA.
 
 FILE = scrt4
 170 CALL READ (*400,*470,scrt4,idrec,146,eor,nwds)
 minor = idrec(3)
 
 
!     SET DISPLACEMENT, VELOCITY, OR ACCELERATION OFP MAJOR-ID IF
!     INFILE IS MODAL DISPLACEMETNS = EIGENVECTORS...
 
 IF (itype1 /= 7) GO TO 175
 idrec(2) = dvamid(ipass)
 175 IF (.NOT.trnsnt) idrec(2) = idrec(2) + 1000
 
!     RESET APPROACH CODE FROM EIGENVALUE TO TRANSIENT OR FREQUENCY
 
 iapp = 5
 IF (trnsnt) iapp = 6
 idrec(1) = 10*iapp + device
 lentry = idrec(10) - nwords
 ncols = lentry
 IF (.NOT.trnsnt) lentry = lentry + lentry
 
!     IF FREQUENCY RESPONSE PROBLEM AND THIS IS THE VELOCITY OR
!     ACCELERATION PASS THEN MOVE DOWN ANY XY LIST OF POINTS AND
!     ADD AN OMEGA TABLE.  SOMETIMES THE MOVEDOWN OF THE XY LIST IS
!     REDUNDANT.
 
!     XY LIST IS MOVED FROM BOTTOM UP INCASE XY LIST IS LONGER THAN
!     THE OMEGA LIST WILL BE.
 
 IF (trnsnt .OR. ipass == 1) GO TO 177
 nomega = iomega + nsols - 1
 IF (lxy == 0) GO TO 177
 jxy = nxy
 kxy = nomega + lxy
 DO  i = 1,lxy
   z(kxy) = z(jxy)
   jxy = jxy - 1
   kxy = kxy - 1
 END DO
 
 177 ixy = nomega + 1
 nxy = ixy + lxy - 1
 idata = nxy + 1
 ndata = idata + lentry*nsols - 1
 typout= 3
 IF (trnsnt) typout = 1
 
!     FILL TITLE, SUBTITLE, AND LABEL FROM CASECC FOR THIS SUBCASE.
 
 DO  i = 1,96
   idrec(i+50) = z(icc+i+37)
 END DO
 idrec(4) = subcas
 
!     CHECK FOR SUFFICIENT CORE.
 
 IF (ndata < buf3) GO TO 190
 insuf = ndata - buf3
 WRITE  (outpt,180) uwm,outfil,insuf
 180 FORMAT (a25,' 2338.  (DDRMM2-3)  DATA BLOCK',i5,  &
     ' MAY NOT BE FULLY COMPLETED DUE TO A CORE INSUFFICIENCY',  &
     /5X,'OF APPROXIMATELY',i11,' DECIMAL WORDS.')
 GO TO 440
 
!     LOOP ON ID-S AVAILABLE FROM THE MAP
 
 
!     COMPUTE OMEGAS IF NECESSARY
!     (NOTE, VELOCITY PASS MAY NOT ALWAYS OCCUR)
 
 190 IF (trnsnt .OR. ipass == 1) GO TO 195
 jlist = iomega - 1
 DO  i = ilist,nlist
   jlist = jlist + 1
   rz(jlist) = rz(i)*twopi
 END DO
 IF (ipass == 2) GO TO 195
 DO  i = iomega,nomega
   rz(i) = -rz(i)**2
 END DO
 
 195 CALL READ (*460,*170,scrt4,buf,nwords,noeor,nwds)
 lminor = .true.
 IF (itype1 /= 5 .OR. savdat(minor) == 0)  GO TO 196
 npos   = savdat(minor)/100
 nstxtr = savdat(minor) - npos*100
 CALL READ (*460,*470,scrt4,bufsav(1),nstxtr,noeor,nwds)
 lminor = .false.
 196 CONTINUE
 
!     PREPARE AND OUTPUT THE OFP-ID-RECORD AFTER FIRST ENTRY IS COMBINED
!     AS IN THE CASE OF A FREQUENCY COMPLEX PROBLEM.
 
 idout = .false.
 idrec(5) = buf(1)
 
!     SET STRESS OR FORCE COMPLEX DATA PTRS IF NECESSARY.
 
 IF (trnsnt) GO TO 220
 IF (itype1 == 4) GO TO 200
 IF (itype1 == 5) GO TO 210
 GO TO 220
 
!     FORCES ASSUMED
 
 200 ielem = (idrec(3)-1)*incr
 lsf   = elem(ielem+19)
 nptsf = elem(ielem+21)
 GO TO 220
 
!     STRESSES ASSUMED
 
 210 ielem = (idrec(3)-1)*incr
 lsf   = elem(ielem+18)
 nptsf = elem(ielem+20)
 GO TO 220
 
!     UNPACK DATA FOR ALL COMPONENTS AND ALL SOLUTION STEPS
!     FOR THIS ID.  (NCOLS COLUMNS ARE NEEDED)
 
 
!     ZERO THE DATA SPACE
 
 220 DO  i = idata,ndata
   z(i) = 0
 END DO
 
!     UNPACK NOW-ZERO TERMS.
 
 jdata = idata - lentry
 DO  i = 1,ncols
   CALL intpk (*260,scrt6,0,typout,0)
   
!     COLUMN I HAS ONE OR MORE NON-ZEROES AVAILABLE.
   
   240 CALL zntpki
   itemp = jdata + irowo*lentry
   IF (.NOT.trnsnt) THEN
      SELECT CASE ( ipass )
       CASE (    1)
         GO TO 246
       CASE (    2)
         GO TO 247
       CASE (    3)
         GO TO 248
     END SELECT
   END IF
   
!     TRANSIENT OUTPUTS
   
   rz(itemp) = aout(1)
   IF (ieol > 0) THEN
     GO TO   260
   ELSE
     GO TO   240
   END IF
   
!    DISPLACEMENTS, AND SPCFS (FREQ RESPONSE)
   
   246 rz(itemp) = aout(1)
   itemp     = itemp + ncols
   rz(itemp) = aout(2)
   IF (ieol > 0) THEN
     GO TO   260
   ELSE
     GO TO   240
   END IF
   
!     VELOCITIES  (FREQ RESPONSE)
   
   247 klist     = iomega + irowo - 1
   rz(itemp) =-rz(klist)*aout(2)
   itemp     = itemp + ncols
   rz(itemp) = rz(klist)*aout(1)
   IF (ieol > 0) THEN
     GO TO   260
   ELSE
     GO TO   240
   END IF
   
!     ACCELERATIONS (FREQ RESPONSE)
   
   248 klist     = iomega + irowo - 1
   rz(itemp) = rz(klist)*aout(1)
   itemp     = itemp + ncols
   rz(itemp) = rz(klist)*aout(2)
   IF (ieol > 0) THEN
     GO TO   260
   ELSE
     GO TO   240
   END IF
   260 jdata = jdata + 1
 END DO
 
!     OUTPUT LINES OF DATA COMBINING THEM FOR COMPLEX REAL/IMAGINARY OR
!     MAG/PHASE OFP FORMATS IF NECESSARY.
 
 jlist = ilist - 1
 DO  i = idata,ndata,lentry
   jwords = nwords
   ij = i + ncols - 1
   DO  j = i,ij
     jwords = jwords + 1
     buf(jwords) = z(j)
     IF (trnsnt) CYCLE
     itemp = j + ncols
     buf(jwords+75) = z(itemp)
   END DO
   
!     IF TRANSIENT, ENTRY IS NOW READY FOR OUTPUT.
   
   IF (trnsnt)  GO TO 365
   
!     MAP COMPLEX OUTPUTS TOGETHER PER -COMPLX- ARRAY.
   
   IF (itype1 == 4 .OR. itype1 == 5) GO TO 300
   
!     POINT DATA
   
   DO  k = 3,8
     IF (FORM == 3) CALL magpha (bufa(k),bufb(k))
     bufa(k+6) = bufb(k)
   END DO
   jwords = 14
   GO TO 370
   
!     ELEMENT STRESS OR FORCE DATA.
   
   300 iout = 0
   l = nptsf
   IF (lminor)  GO TO 310
   DO  k = 1,nstxtr
     j = savpos(npos+k-1)
     buf(j) = bufsav(k)
   END DO
   310 npt = complx(l)
   IF (npt < 0) THEN
     GO TO   320
   ELSE IF (npt == 0) THEN
     GO TO   350
   ELSE
     GO TO   340
   END IF
   320 npt = -npt
   IF (FORM /= 3) GO TO 340
   
!     COMPUTE MAGNITUDE/PHASE
   
   CALL magpha (bufa(npt),bufb(npt))
   330 iout = iout + 1
   elwork(iout) = bufa(npt)
   l = l + 1
   GO TO 310
   340 IF (npt <= lsf) GO TO 330
   npt  = npt - lsf
   iout = iout + 1
   elwork(iout) = bufb(npt)
   l = l + 1
   GO TO 310
   
!     MOVE OUTPUT DATA
   
   350 DO  l = 1,iout
     buf(l) = elwork(l)
   END DO
   jwords = iout
   GO TO 370
   365 CONTINUE
   IF (lminor)  GO TO 370
   DO  k = 1,nstxtr
     j = savpos(npos+k-1)
     buf(j) = bufsav(k)
   END DO
   
!     CALL DDRMMS TO RECOMPUTE SOME ELEMENT STRESS QUANTITIES
!     IN TRANSIENT PROBLEMS ONLY.
   
   370 IF (trnsnt .AND. itype1 == 5) CALL ddrmms (buf,idrec(3),buf4,buf5)
   IF (idout) GO TO 380
   idrec( 9) = FORM
   idrec(10) = jwords
   CALL WRITE (outfil,idrec,146,eor)
   idout = .true.
   380 jlist = jlist + 1
   rbuf(1) = rz(jlist)
   CALL WRITE (outfil,buf,jwords,noeor)
 END DO
 CALL WRITE (outfil,0,0,eor)
 
!     GO FOR NEXT OUTPUT ID
 
 GO TO 190
 
!  END OF DATA ON MAP FILE (SCRT4).
 
 400 CALL CLOSE (outfil,cls)
 CALL CLOSE (infile,clsrew)
 CALL CLOSE (scrt4,clsrew)
 CALL CLOSE (scrt6,clsrew)
 ipass = ipass + 1
 IF (ipass > passes) GO TO 410
 
!     PREPARE FOR ANOTHER PASS
 
 FILE = infile
 CALL OPEN (*450,infile,z(buf1),rdrew)
 CALL fwdrec (*460,infile)
 GO TO 20
 410 RETURN
 
!     CHANGE IN MAJOR OFP-ID DETECTED ON -INFILE-.
 
 420 WRITE  (outpt,430) swm,infile
 430 FORMAT (a27,' 2339.  (DDRMM2-1) A CHANGE IN WORD 2 OF THE OFP-ID',  &
     ' RECORDS OF DATA BLOCK',i5, /5X,'HAS BEEN DETECTED. ',  &
     ' POOCESSING OF THIS DATA BLOCK HAS BEEN TERMINATED.')
 440 ipass = 3
 GO TO 400
 
!     UNDEFINED FILE.
 
 450 RETURN 1
 
!     END OF FILE
 
 460 RETURN 2
 
!     END OF RECORD.
 
 470 RETURN 3
 
!     INSUFFICIENT CORE
 
 480 RETURN 4
END SUBROUTINE ddrmm2
