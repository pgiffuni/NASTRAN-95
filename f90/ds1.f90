SUBROUTINE ds1 (iarg)
     
!     THIS ROUTINE CREATES THE SCRATCH FILE ECPTDS BY APPENDING TO EACH
!     ELEMENT IN THE ECPT AN ELEMENT DEFORMATION, AN AVERAGE ELEMENT
!     LOADING TEMPERATURE, AND THE PROPER COMPONENTS OF THE DISPLACEMENT
!     VECTORS. SUBROUTINE DS1A READS THE ECPTDS IN THE SAME WAY AS SMA1A
!     READS THE ECPT IN ORDER TO CREATE A SECOND ORDER APPROXIMATION TO
!     THE KGG, WHICH IS CALLED KGGD.
!     IF DS1 CANNOT FIND ANY ELEMENTS IN THE ECPT WHICH ARE IN THE SET
!     OF ELEMENTS FOR WHICH DIFFERENTIAL STIFFNESS IS DEFINED, IARG IS
!     RETURNED CONTAINING A ZERO TO THE CALLING ROUTINE, DSMG1.
 
 
 INTEGER, INTENT(OUT)                     :: iarg
 EXTERNAL        rshift
 LOGICAL :: dstype,eorflg,endid,record
 INTEGER :: buffr1,buffr2,eor,clsrw,outrw,casecc,gptt,edt,  &
     ugv,ecpt,ecptds,FILE,tsetno,dsetno,tmpset,  &
     recno,edtloc,edtbuf,eltype,elid,bufloc,dfmset,  &
     rshift,jsil(2),oldel,ccbuf,oldeid,buffr3
 DIMENSION       tgrid(33),iz(1),xecpt(328),iecpt(328),ccbuf(2),  &
     gptbf3(3),NAME(2),edtbuf(3),edtloc(2),mcbugv(7)
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim,sfm
 COMMON /machin/ mach,ihalf
 COMMON /gpta1 / nelems,last,incr,NE(1)
 COMMON /system/ isys,sysdum(25),mn,xxx18(18),ndum(9)
 COMMON /zzzzzz/ z(1)
 COMMON /unpakx/ itypeb,iunpk,junpk,incupk
 COMMON /ds1ett/ eltype,oldel,eorflg,endid,bufflg,tsetno,fdfalt,  &
     iback,record,oldeid
 COMMON /BLANK / dscset
 EQUIVALENCE     (z(1),iz(1))       ,(xecpt(1),iecpt(1)),  &
     (gptbf3(1),tmpset) ,(gptbf3(2),idfalt) ,  &
     (gptbf3(3),recno)  ,(edtbuf(1),dfmset) ,  &
     (edtbuf(2),elid)   ,(edtbuf(3),deform) , (ioutpt,sysdum(1))
 DATA            edtloc/104,1 /,     nskip/ 137  /
 DATA            casecc,gptt,edt,ugv,ecpt,ecptds /  &
     101,   102, 104,105,108, 301    /
 DATA            NAME  /4HDS1 ,4H         /
 DATA            inrw,outrw,eor,neor,clsrw/ 0,1,1,0,1 /
 
!     SET IARG TO ZERO
 
 CALL delset
 iarg = 0
 
!     DETERMINE SIZE OF AVAILABLE CORE, DEFINE 2 BUFFERS AND INITIALIZE
!     OPEN CORE POINTERS AND COUNTERS.
 
 izmax  = korsz(z)
 buffr1 = izmax  - isys
 buffr2 = buffr1 - isys
 buffr3 = buffr2 - isys
 bufloc = izmax  - isys - 3
 ileft  = buffr3 - 1
 left   = ileft  - nelems - 2
 isil   = 0
 nsil   = 0
 iedt   = 0
 nedt   = 0
 
!     SET DIFFERENTIAL STIFFNESS FLAGS FOR ALL ELEMENT TYPES TO ZERO
 
 DO  i = 1,nelems
   iz(left+i) = 0
 END DO
 
!     OPEN CASECC, SKIP HEADER, SKIP 5 WORDS AND READ DEFORMATION SET
!     NUMBER AND LOADING TEMPERATURE SET NUMBER.
 
 CALL gopen (casecc,z(buffr1),inrw)
 CALL fread (casecc,0,-5,neor)
 CALL fread (casecc,ccbuf,2,neor)
 dsetno = ccbuf(1)
 tsetno = ccbuf(2)
 
!     STORE THE DIFFERENTIAL STIFFNESS COEFFICIENT (BETA) SET NUMBER
!     IN COMMON.  THIS WORD IS THE 138TH WORD OF THE 2ND RECORD OF CASE
!     CONTROL.
 
 FILE = casecc
 CALL fwdrec (*400,casecc)
 CALL fread  (casecc,0,-nskip,neor)
 CALL fread  (casecc,dscset,1,neor)
 CALL CLOSE  (casecc,clsrw)
 
!     IS THERE A TEMPERATURE LOAD
 
 record =.false.
 iback  = 0
 IF (tsetno <= 0) GO TO 60
 
!     THERE IS. OPEN THE GPTT, SKIP FIRST TWO WORDS OF THE HEADER RECORD
!     AND READ 3 WORD ENTRIES OF THE HEADER RECORD UNTIL A SET NUMBER
!     MATCHES THE SET NUMBER READ IN THE CASE CONTROL RECORD.
 
 FILE = gptt
 CALL OPEN  (*400,gptt,z(buffr3),inrw)
 CALL fread (gptt,0,-2,neor)
 20 CALL fread (gptt,gptbf3,3,neor)
 IF (tmpset == tsetno) GO TO 30
 GO TO 20
 30 fdfalt = gptbf3(2)
 IF (recno  /=  0) GO TO 40
 IF (idfalt == -1) CALL mesage (-30,29,tsetno)
 CALL CLOSE (gptt,clsrw)
 GO TO 60
 
!     POSITION GPTT TO DESIRED TEMPERATURE RECORD
 
 40 CALL REWIND (gptt)
 DO  i = 1,recno
   CALL fwdrec (*410,gptt)
 END DO
 record =.true.
 
!     READ SETID AND VERIFY FOR CORRECTNESS
 
 CALL fread (gptt,idset,1,0)
 IF (tsetno /= idset) CALL mesage (-30,29,tsetno)
 
!     INITIALIZE /DS1ETT/ VARIABLES
 
 oldeid = 0
 oldel  = 0
 eorflg =.false.
 endid  =.true.
 
!     DETERMINE IF AN ENFORCED DEFORMATION SET IS CALLED FOR.
 
 60 iedt = isil
 i    = isil
 IF (dsetno <= 0) GO TO 90
 FILE = edt
 CALL preloc (*90,z(bufloc),edt)
 CALL locate (*450,z(bufloc),edtloc,iflag)
 70 CALL READ (*410,*80,edt,edtbuf,3,neor,iflag)
 IF (dfmset /= dsetno) GO TO 70
 iz(i+1) = elid
 z (i+2) = deform
 nedt = nedt + 2
 i    = i + 2
 left = left - 2
 IF (left <= 0) CALL mesage (-8,0,NAME)
 GO TO 70
 80 CALL CLOSE (edt,clsrw)
 low = iedt + 1
 lim = iedt + nedt
 
!     READ THE UGV INTO CORE.
 
 90 CALL gopen (ugv,z(buffr1),inrw)
 idisp = iedt + nedt
 mcbugv(1) = ugv
 CALL rdtrl (mcbugv(1))
 IF (left < mcbugv(3)) CALL mesage (-8,0,NAME(1))
 itypeb = 1
 iunpk  = 1
 junpk  = mcbugv(3)
 incupk = 1
 CALL unpack (*460,ugv,z(idisp+1))
 CALL CLOSE  (ugv,clsrw)
 
!     OPEN THE ECPTDS AND ECPT FILES.
 
 CALL gopen (ecptds,z(buffr2),outrw)
 CALL gopen (ecpt,z(buffr1),inrw)
 
!     READ THE PIVOT POINT (1ST WORD).
 
 100 FILE   = ecpt
 imhere = 100
 eltype = -1
 j      = -1
 CALL READ (*390,*430,ecpt,npvt,1,neor,iflag)
 ind = 0
 110 dstype =.false.
 
!     READ ELEMENT TYPE (2ND WORD)
 
 CALL READ (*410,*370,ecpt,eltype,1,neor,iflag)
 IF (eltype < 1 .OR. eltype > nelems) GO TO 480
 
!     READ ELEMENT ID (3RD WORD, BEGINNING OF J NO. OF WORDS)
 
 imhere = 115
 CALL READ (*410,*430,ecpt,iecpt,1,neor,iflag)
 IF (iback == 0) GO TO 120
 IF (eltype == oldel .AND. iecpt(1) >= oldeid) GO TO 130
 CALL bckrec (gptt)
 
!     RESET /DS1ETT/ VARIABLES
 
 iback  = 0
 oldeid = 0
 oldel  = 0
 eorflg =.false.
 endid  =.true.
 CALL READ (*410,*420,gptt,idset,1,0,flag)
 IF (tsetno /= idset) CALL mesage (-30,29,tsetno)
 
 120 idx = (eltype-1)*incr
 ntemp = 1
!                IS2D8              IHEX1              IHEX3
 IF (eltype == 80 .OR. (eltype >= 65 .AND. eltype <= 67))  &
     ntemp = NE(idx+15) - 1
 
!     READ ECPT ENTRY FOR THIS ELEMENT (J-1 WORDS)
 
 130 j = NE(idx+12)
 IF (NE(idx+24) /= 0) dstype = .true.
 imhere = 130
 CALL READ (*410,*430,ecpt,xecpt(2),j-1,neor,iflag)
 
!     IS THIS ELEMENT IN THE SET OF DS ELEMENTS.
 
 IF (dstype) GO TO 150
 IF (iz(left+eltype) == 1) GO TO 110
 iz(left+eltype) = 1
 CALL page2 (-2)
 WRITE  (ioutpt,140) uwm,NE(idx+1),NE(idx+2),eltype
 140 FORMAT (a25,' 3117, DIFFERENTIAL STIFFNESS CAPABILITY NOT DEFINED'  &
     ,      ' FOR ',2A4,' ELEMENTS (ELEMENT TYPE ',i3,2H).)
 GO TO 110
 150 iarg = 1
 
!     DETERMINE IF THE ELEMENT IS A CONE.  IF IT IS, IT MUST HAVE A
!     NONZERO MEMBRANE THICKNESS FOR IT TO BE ADMISSIBLE TO THE ECPTDS.
 
 IF (eltype /= 35) GO TO 170
!                 CONEAX
 ntemp = 2
 IF (xecpt(5) == 0.0) GO TO 110
 
!     DETERMINE THE NUMBER OF RINGAX POINTS FROM THE 27TH WORD OF
!     /SYSTEM/.
 
 nrngax = rshift(mn,ihalf)
 
!     DETERMINE THE HARMONIC NUMBER, IHARM, FROM THE ELEMENT IDENT.
!     NUMBER, IECPT(1)
 
 itemp = iecpt(1)/1000
 iharm = iecpt(1) - itemp*1000 - 1
 
!     DETERMINE THE SIL NUMBERS, SIL(1) AND SIL(2), WHICH WILL BE USED
!     TO APPEND TEMPERATURES AND DISPLACEMENT VECTORS.
 
 IF (iharm /= 0) GO TO 160
 jsil(1) = iecpt(2)
 jsil(2) = iecpt(3)
 GO TO 180
 160 itemp   = 6*iharm*nrngax
 jsil(1) = iecpt(2) - itemp
 jsil(2) = iecpt(3) - itemp
 GO TO 180
 
!     IF WE ARE DEALING WITH A TRIA1 OR QUAD1 ELEMENT, IT MUST HAVE A
!     NONZERO MEMBRANE THICKNESS FOR IT TO BE ADMISSIBLE TO THE ECPTDS.
 
 170 IF (eltype /= 6 .AND. eltype /= 19) GO TO 180
!               TRIA1              QUAD1
 kk = 7
 IF (eltype == 19) kk = 8
!                  QUAD1
 IF (xecpt(kk) == 0.0) GO TO 110
 
!     WRITE PIVOT POINT
 
 180 IF (ind == 0) CALL WRITE (ecptds,npvt,1,neor)
 ind = 1
 IF (eltype /= 34) GO TO 200
!                    BAR
 
!     THE ELEMENT IS A BAR.  THE ECPT ENTRY WILL BE REARRANGED SO THAT
!     THE DBAR SUBROUTINE MAY BE CALLED IN SUBROUTINE DS1A.
 
 eltype = 2
!           BEAM
 
!     IF THE COUPLED MOMENT OF INERTIA TERM I12 (=ECPT(33)) IS NON-ZERO
!     SET I12 = 0.0, WRITE WARNING MESSAGE AND PROCEED.
 
 IF (xecpt(33) == 0.0) GO TO 190
 xecpt(33) = 0.0
 CALL mesage (30,111,iecpt(1))
 190 xecpt(47) = xecpt(42)
 xecpt(46) = xecpt(41)
 xecpt(45) = xecpt(40)
 xecpt(44) = xecpt(39)
 xecpt(43) = xecpt(38)
 xecpt(42) = xecpt(37)
 xecpt(41) = xecpt(36)
 xecpt(40) = xecpt(35)
 xecpt(39) = xecpt(34)
 xecpt(29) = xecpt(31)
 xecpt(30) = xecpt(32)
 xecpt(28) = xecpt(21)
 xecpt(27) = xecpt(20)
 xecpt(25) = xecpt(19)
 xecpt(24) = xecpt(18)
 xecpt(21) = xecpt(17)
 xecpt(20) = xecpt(16)
 j = 47
 
!     WRITE ELEMENT TYPE
 
 200 CALL WRITE (ecptds,eltype,1,neor)
 
!     ATTACH THE ELEMENT DEFORMATION TO THE XECPT ARRAY.
 
 j = j + 1
 nogpts = NE(idx+10)
 xecpt(j) = 0.0
 IF (dsetno > 0) GO TO 210
 GO TO 230
 
!     SEARCH THE EDT TO FIND AN ELEMENT NO. IN THE TABLE CORRESPONDING
!     TO THE CURRENT ELEMENT NO., IECPT(1).  IF IT CANNOT BE FOUND NOTE
!     THE ELEMENT DEFORMATION, IECPT(J), HAS BEEN SET TO ZERO.
 
 210 DO  i = low,lim,2
   IF (iz(i) /= iecpt(1)) CYCLE
   xecpt(j) = z(i+1)
   GO TO 230
 END DO
 
!     APPEND THE LOADING TEMPERATURE(S) TO THE XECPT ARRAY
 
 230 IF (eltype == 2) eltype = 34
!                  BEAM          BAR
 CALL ds1etd (iecpt(1),tgrid,ntemp)
 IF (eltype /= 34) GO TO 240
!                    BAR
 eltype = 2
 IF (tsetno <= 0) GO TO 240
 tgrid(1) = (tgrid(1) + tgrid(2))*0.5
 240 iii = 1
 IF (eltype /= 80) GO TO 250
!                  IS2D8
 j = j + 1
 iecpt(j) = tsetno
 iii = 2
 250 CONTINUE
 DO  i = iii,ntemp
   j = j + 1
   xecpt(j) = tgrid(i)
 END DO
 
!     NOW ATTACH THE DISPLACEMENT VECTORS
 
 j = j + 1
 IF (eltype == 35) GO TO 330
!                 CONEAX
 IF (eltype == 2 .OR. eltype == 75) GO TO 290
!                 BEAM             TRSHL
 IF (eltype < 53 .OR. eltype > 61) GO TO 280
!                 DUM1              DUM9
 
 
!     DUMMY ELEMENTS
 
 IF (MOD(ndum(eltype-52),10) == 6) GO TO 290
 280 nwds = 3
 GO TO 300
 290 nwds = 6
 300 DO  i = 1,nogpts
   INDEX = idisp + iecpt(i+1)
   DO  i1 = 1,nwds
     xecpt(j) = z(INDEX)
     INDEX = INDEX + 1
     j = j + 1
   END DO
 END DO
 GO TO 360
 
!     APPEND THE ZERO HARMONIC COMPONENTS OF THE DISPLACEMENT VECTOR.
!     NOTE THAT FOR A CONICAL SHELL ELEMENT DIRECT POINTERS INTO THE
!     DISPLACEMENT VECTOR ARE SIL(1) AND SIL(2).
 
 330 DO  j1 = 1,2
   DO  i  = 1,6
     INDEX = idisp + jsil(j1) + i - 1
     xecpt(j) = z(INDEX)
     j = j + 1
   END DO
 END DO
 
!     THE APPENDED ECPT, ECPTDS, IS NOW COMPLETE.
 
 360 CALL WRITE (ecptds,xecpt,j-1,neor)
 GO TO 110
 
!    IF IND = 0, THEN NO ELEMENTS IN THE CURRENT ECPT RECORD ARE IN THE
!    DS ELEMENT SET.  WRITE A -1 FOR THIS PIVOT POINT.
 
 370 IF (ind /= 0) GO TO 380
 CALL WRITE (ecptds,-1,1,eor)
 GO TO 100
 
!     WRITE AN EOR ON THE ECPTDS FILE
 
 380 CALL WRITE (ecptds,0,0,eor)
 GO TO 100
 
!     CLOSE BOTH FILES
 
 390 CALL CLOSE (ecpt,clsrw)
 CALL CLOSE (gptt,clsrw)
 CALL CLOSE (ecptds,clsrw)
 RETURN
 
!     FATAL ERROR RETURNS
 
 400 j = -1
 GO TO 470
 410 j = -2
 GO TO 470
 420 FILE = gptt
 430 j = -3
 IF (FILE == ecpt) WRITE (ioutpt,440) imhere,eltype,j
 440 FORMAT (/,'0*** DS1/IMHERE,ELTYPE,J = ',3I5)
 GO TO 470
 450 j = -4
 GO TO 470
 460 CALL mesage (-30,83,NAME(1))
 470 CALL mesage (j,FILE,NAME)
 480 WRITE  (ioutpt,490) sfm,eltype
 490 FORMAT (a25,' 2147, ILLEGAL ELEMENT TYPE =',i10,  &
     ' ENCOUNTERED BY DSMG1 MODULE.')
 CALL mesage (-61,0,NAME)
 RETURN
END SUBROUTINE ds1
