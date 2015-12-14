SUBROUTINE ta1a
     
!     TA1A BUILDS THE ELEMENT SUMMARY TABLE (EST).
!     THE EST GROUPS ECT, EPT, BGPDT AND ELEMENT TEMP. DATA FOR EACH
!     SIMPLE ELEMENT OF THE STRUCTURE. THE EST CONTAINS ONE LOGICAL
!     RECORD PER SIMPLE ELEMENT TYPE.
 
 IMPLICIT INTEGER (a-z)
 LOGICAL :: eorflg,endid ,record,frstim,q4t3
 INTEGER :: zeros(4)     ,buf(50)      ,nam(2),gpsav(34)    ,  &
     pcomp(2)     ,pcomp1(2)    ,pcomp2(2)           , ipshel(16)
 REAL :: deftmp,tlam  ,zoffs ,zz(1) ,bufr(50)            ,  &
     tgrid(33)    ,rpshel(16)
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm   ,uwm   ,uim   ,sfm
 COMMON /BLANK / luset ,nosimp,nosup ,nogenl,genl  ,comps
 COMMON /ta1com/ nsil  ,ect   ,ept   ,bgpdt ,sil   ,gptt  ,cstm  ,  &
     mpt   ,est   ,gei   ,gpect ,ecpt  ,gpct  ,mptx  ,  &
     pcomps,eptx  ,scr1  ,scr2  ,scr3  ,scr4
 COMMON /system/ ksystm(65)
 COMMON /gpta1 / nelem ,last  ,incr  ,elem(1)
 COMMON /names / rd    ,rdrew ,wrt   ,wrtrew,clsrew,cls
 COMMON /ta1ett/ eltype,oldel ,eorflg,endid ,bufflg,itemp ,idftmp,  &
     iback ,record,oldeid
 COMMON /ta1acm/ ig(90)
 COMMON /two   / ktwo(32)
 COMMON /zzzzzz/ z(1)
 EQUIVALENCE     (ksystm( 1),sysbuf ) , (ksystm( 2),nout ) ,  &
     (ksystm(10),tempid ) , (ksystm(56),iheat) ,  &
     (idftmp    ,deftmp ) , (bufr(1)   ,buf(1)),  &
     (z(1)      ,zz(1)  ) , (ipshel( 1),rpshel(1))
 DATA    nam   / 4HTA1A,4H    /
 DATA    zeros / 4*0   /
 DATA    bar   / 34    /
 DATA    hbdy  / 52    /
 DATA    qdmem2/ 63    /
 DATA    quad4 / 64    /
 DATA    tria3 / 83    /
 DATA    pcomp / 5502, 55 /
 DATA    pcomp1/ 5602, 56 /
 DATA    pcomp2/ 5702, 57 /
 DATA    sym   / 1     /
 DATA    mem   / 2     /
 DATA    symmem/ 3     /
 
!     PERFORM GENERAL INITIALIZATION
 
 IF (nelem > 90) GO TO 1190
 buf1  = korsz(z) - sysbuf - 2
 buf2  = buf1 - sysbuf - 3
 buf3  = buf2 - sysbuf
 frstim= .true.
 lstprp= 0
 kscalr= 0
 itabl = 0
 nosimp=-1
 nogox = 0
 nogo  = 0
 m8    =-8
 comps = 1
 nopshl=-1
 npshel= 0
 oldid = 0
 CALL sswtch (40,l40)
 
!     READ THE ELEMENT CONNECTION TABLE.
!     IF PROPERTY DATA IS DEFINED FOR THE ELEMENT, READ THE EPT INTO
!     CORE AND SORT IF REQUIRED. THEN FOR EACH ECT ENTRY, LOOK UP AND
!     ATTACH THE PROPERTY DATA. WRITE ECT+EPT ON SCR1.
!     IF PROPERTY DATA NOT DEFINED FOR ELEMENT, COPY ECT DATA ON SCR1.
!     IF NO SIMPLE ELEMENTS IN ECT, RETURN.
 
!     FOR THE PLATE AND SHELL ELEMENTS REFERENCING PCOMP, PCOMP1 OR
!     PCOMP2 BULK DATA ENTRIES, PROPERTY DATA IN THE FORM OF PSHELL
!     BULK DATA ENTRY IS CALLED AND WRITTEN TO SCR1
 
 FILE = ect
 CALL OPEN (*540,ect,z(buf1),rdrew)
 CALL skprec (ect,1)
 FILE = scr1
 CALL OPEN (*1100,scr1,z(buf3),wrtrew)
 buf(1) = ept
 CALL rdtrl(buf)
 noept = buf(1)
 IF (buf(1) < 0) GO TO 10
 CALL preloc (*10,z(buf2),ept)
 
!     LOCATE, ONE AT A TIME, SIMPLE ELEMENT TYPE IN ECT. IF PRESENT,
!     WRITE POINTER ON  SCR1. SET POINTERS AND, IF DEFINED, LOCATE AND
!     READ ALL PROPERTY DATA FROM EPT.
 
 10 CALL ectloc (*200,ect,buf,i)
 id = -1
 eltype = elem(i+2)
 CALL WRITE (scr1,i,1,0)
 q4t3 = .false.
 IF (eltype == quad4 .OR. eltype == tria3) q4t3 = .true.
 IF (elem(i+10) == 0) kscalr = 1
 m  = elem(i+5)
 mm = elem(i+8)
 IF (mm == 0) GO TO 120
 mx = mm
 noprop = 0
 IF (elem(i+6) /= lstprp) GO TO 20
 IF (eltype    == qdmem2) noprop = 1
 GO TO 80
 20 IF (noept < 0) GO TO 1130
 
!     LOCATE PROPERTY CARD
 
 ll = 0
 CALL locate (*50,z(buf2),elem(i+6),flag)
 noprop = 1
 30 lstprp = elem(i+6)
 40 IF (ll+mm >= buf3) CALL mesage (-8,0,nam)
 IF (noprop == 0) GO TO 80
 CALL READ (*1110,*55,ept,z(ll+1),mm,0,flag)
 ll = ll + mm
 GO TO 40
 
!     CHECK FOR QUAD4 AND TRIA3 ELEMENTS WITH ONLY PCOMP CARDS
 
!     SET POINTER FOR NO PSHELL DATA, AND
!     READ PCOMP, PCOMP1 AND PCOMP2 DATA INTO CORE
 
 50 IF (.NOT.q4t3) GO TO 60
 nopshl = 1
 GO TO 700
 
!     CHECK FOR QUAD4 AND TRIA3 ELEMENTS WITH BOTH PCOMP AND PSHELL
!     CARDS
 
!     IF LL.GT.0 HERE, PSHELL DATA IS PRESENT,
!     NEED TO CHECK THE PRESENCE OF PCOMP DATA, AND RESET NOPSHL POINTER
!     IF NECCESSARY
 
!     EVENTUALLY, WE WILL HAVE
 
!     NOPSHL =-1, LOGIC ERROR FOR QUAD4/TRIA3 PROPERTY DATA
!            = 0, ONLY PSHELL DATA PRESENT
!            = 1, ONLY PCOMP TYPE DATA PRESENT
!            = 2, BOTH PSHELL AND PCOMP DATA PRESENT (SEE STA.760)
 
 55 IF (.NOT.q4t3) GO TO 70
 IF (ll <= 0) GO TO 700
 nopshl = 0
 GO TO 70
 
 60 IF (noprop == 0) GO TO 1130
 
!     Z(1) THRU Z(LL) CONTAIN PROPERTY DATA
 
 70 IF (mm <= 4) CALL sort (0,0,mm,1,z(1),ll)
 kn = ll/mm
 IF (nopshl == 0) GO TO 700
 
!     READ ECT DATA FOR ELEMENT. LOOK UP PROPERTY DATA IF CURRENT ELEM.
!     HAS A PROPERTY ID DIFFERNENT FROM THAT OF THE PREVIOUS ELEMENT.
!     WRITE ECT + EPT (OR NEW GENERATED PSHELL) DATA ON SCR1.
 
 80 CALL READ (*1110,*140,ect,buf,m,0,flag)
 nosimp = nosimp + 1
 IF (buf(2) /= id) noprop = 1
 id = buf(2)
 buf(2) = buf(1)
 buf(1) = m + mm - 2
 IF (noprop == 0) GO TO 90
 IF (q4t3 .AND. nopshl == 1) GO TO 800
 npshel = 0
 GO TO 600
 90 CALL WRITE (scr1,buf(1),m,0)
 IF (.NOT.q4t3) GO TO 100
 IF (npshel == 1) GO TO 110
 100 CALL WRITE (scr1,z(kx+2),mm-1,0)
 npshel = 0
 noprop = 0
 GO TO 80
 110 CALL WRITE (scr1,ipshel(1),mm-1,0)
 noprop = 0
 GO TO 80
 
!     EPT DATA NOT DEFINED FOR ELEMENT. COPY ECT DATA ON SCR1.
 
 120 buf(1) = m
 m1 = m + 1
 130 CALL READ  (*1110,*140,ect,buf(2),m,0,flag)
 CALL WRITE (scr1,buf(1),m1,0)
 nosimp = nosimp + 1
 GO TO 130
 140 CALL WRITE (scr1,0,0,1)
 GO TO 10
 
!     HERE WHEN ALL ELEMENTS HAVE BEEN PROCESSED.
!     IF NONE FOUND, EXIT.
 
 200 CONTINUE
 IF (noept >= 0) CALL CLOSE (ept,clsrew)
 CALL CLOSE (scr1,clsrew)
 IF (nosimp == -1) RETURN
 nosimp = nosimp + 1
 
!     READ THE BGPDT INTO CORE (UNLESS ALL SCALAR PROBLEM).
!     READ THE SIL INTO CORE.
 
 nbgp = 0
 IF (kscalr == 0) GO TO 220
 FILE = bgpdt
 CALL OPEN   (*1100,bgpdt,z(buf1),rdrew)
 CALL fwdrec (*1110,bgpdt)
 CALL READ   (*1110,*210,bgpdt,z(1),buf2,1,nbgp)
 CALL mesage (m8,0,nam)
 210 CALL CLOSE  (bgpdt,clsrew)
 220 FILE = sil
 CALL OPEN   (*1100,sil,z(buf1),rdrew)
 CALL fwdrec (*1110,sil)
 CALL READ   (*1110,*230,sil,z(nbgp+1),buf2-nbgp,1,nsil)
 CALL mesage (m8,0,nam)
 230 CALL CLOSE  (sil,clsrew)
 
!     IF TEMP DEPENDENT MATERIALS IN PROBLEM,
!     OPEN GPTT AND POSITION TO PROPER THERMAL RECORD
 
 record = .false.
 itemp  = tempid
 IF (tempid == 0) GO TO 310
 FILE = gptt
 CALL OPEN (*1160,gptt,z(buf3),rdrew)
 itmpid = nbgp+nsil+3
 CALL READ (*1110,*240,gptt,z(itmpid-2),buf2-itmpid,1,nid)
 CALL mesage (-8,0,nam)
 240 ntmpid = itmpid - 5 + nid
 DO  i = itmpid,ntmpid,3
   IF (tempid == z(i)) GO TO 260
 END DO
 GO TO 1160
 260 idftmp = z(i+1)
 IF (idftmp /= -1) deftmp = zz(i+1)
 n = z(i+2)
 IF (n == 0) GO TO 310
 record =.true.
 n = n - 1
 IF (n == 0) GO TO 280
 DO  i = 1,n
   CALL fwdrec (*1110,gptt)
 END DO
 
!     READ SET ID AND VERIFY FOR CORRECTNESS
 
 280 CALL READ (*1110,*1120,gptt,iset,1,0,flag)
 IF (iset == tempid) GO TO 300
 WRITE  (nout, 290) sfm,iset,tempid
 290 FORMAT (a25,' 4020, TA1A HAS PICKED UP TEMPERATURE SET',i9,  &
     ' AND NOT THE REQUESTED SET',i9)
 CALL mesage (-61,0,0)
 
!     INITIALIZE /TA1ETT/ VARIABLES
 
 300 oldeid = 0
 oldel  = 0
 eorflg =.false.
 endid  =.true.
 
!     LOOP THRU THE ECT+EPT DATA
!     CONVERT INTERNAL GRID POINT INDICES TO SIL VALUES FOR EACH NON-
!     SCALER ELEMENT, ATTACH THE BGPDT DATA AND,
!     IF A TEMPERATURE PROBLEM, COMPUTE THE ELEMENT TEMP FROM THE GPTT
!     DATA OR SUBSTITUTE THE DEFAULT TEMP.
!     WRITE THE RESULT ON THE EST, ONE RECORD PER ELEMENT TYPE
 
 310 CALL OPEN  (*1100,scr1,z(buf1),rdrew)
 CALL OPEN  (*1100,est,z(buf2),wrtrew)
 CALL fname (est,buf)
 CALL WRITE (est,buf,2,1)
 locbgp = 1
 
!     RESET SOME OF THE /TA1ACM/ VALUES IF IT IS A -HEAT- FORMULATION
 
 IF (iheat <= 0) GO TO 320
 
!     TRIARG ELEMENT (TYPE 36)
 ig(36) = 14
 
!     TRAPRG ELEMENT (TYPE 37)
 ig(37) = 14
 
!     REPLACE QDMEM1 ELEMENT (TYPE 62) BY QDMEM ELEMENT (TYPE 16)
 ig(62) = 14
 
!     REPLACE QDMEM2 ELEMENT (TYPE 63) BY QDMEM ELEMENT (TYPE 16)
 ig(63) = 14
 
!     READ POINTER FROM SCR1. WRITE ELEMENT TYPE ON EST.
!     SET POINTERS FOR CONVERSION OF GRID NOS TO SIL VALUES.
 
 320 CALL READ (*500,*1120,scr1,i,1,0,flag)
 eltype = elem(i+2)
 CALL WRITE (est,eltype,1,0)
 
!     ELEMENT TYPE  USED TO INDEX INTO /TA1ACM/
!     AND SET USED  /OPEN CORE/  BLOCKS NEGATIVE
 
 ig(eltype) = -ig(eltype)
 NAME  = elem(i   )
 jscalr= elem(i+10)
 mm    = elem(i+ 9)
 lx    = elem(i+12)
 IF (elem(i+8) == 0) lx = lx + 1
 mm    = lx + mm - 1
 jtemp = elem(i+13)
 ntemp = 1
 IF (jtemp == 4) ntemp = elem(i+14) - 1
!         IHEX1/2/3,TRIM6
 
!     READ ECT + EPT DATA FOR ELEMENT FROM SCR1.
 
 330 CALL READ (*1110,*400,scr1,buf,1,0,flag)
 CALL READ (*1110,*1120,scr1,buf(2),buf(1),0,flag)
 
 IF (nogo /= 0 .OR. nogox /= 0) GO TO 350
 IF (eltype /= bar) GO TO 350
 
!     FOR BAR AND BEAM ELEMENTS, STORE COORDINATES AND
!     COORDINATE SYSTEM ID FOR ORIENTATION VECTOR.
 
 kx = 4*(buf(3)-1) + locbgp
 IF (buf(8) == 1) GO TO 340
 buf(8) = buf(5)
 IF (buf(8) == 0) GO TO 340
 k = 4*(buf(8)-1)  + locbgp
 bufr(5) = zz(k+1) - zz(kx+1)
 bufr(6) = zz(k+2) - zz(kx+2)
 bufr(7) = zz(k+3) - zz(kx+3)
 buf(8)  = 0
 GO TO 350
 340 buf(8)  = z(kx)
 
!     SAVE INTERNAL GRID NOS, THEN CONVERT TO SIL NOS
!     AND WRITE ECT + EPT DATA ON EST.
 
 350 DO  l = lx,mm
   gpsav(l) = 0
   IF (buf(l) == 0) CYCLE
   gpsav(l) = buf(l)
   k = gpsav(l) + nbgp
   buf(l) = z(k)
 END DO
 CALL WRITE (est,buf(2),buf(1),0)
 
!     IF NOT SCALAR ELEMENT, PICK UP BGPDT DATA AND WRITE ON EST.
 
 IF (jscalr /= 0) GO TO 330
 DO  l = lx,mm
   IF (gpsav(l) == 0) GO TO 370
   k = (gpsav(l)-1)*4
   CALL WRITE (est,z(k+1),4,0)
   IF (z(k+1) >= 0) CYCLE
   IF (eltype == hbdy .AND. l > lx+3) CYCLE
   nogo = 1
   CALL mesage (30,131,buf(2))
   CYCLE
   370 CALL WRITE (est,zeros,4,0)
 END DO
 
!     ELEMENT TEMP. IS NOT USED IN CONM1 AND CONM2 (ELEM TYPES 29 30)
 
 tgrid(1) = 0.
 IF (eltype == 29 .OR. eltype == 30) GO TO 390
 
!     IF NOT SCALAR ELEMENT, COMPUTE AND WRITE ELEMENT TEMP ON EST.
 
 CALL ta1etd (buf(2),tgrid,ntemp)
 IF (eltype == bar) tgrid(1) = (tgrid(1)+tgrid(2))/2.0
 390 CALL WRITE (est,tgrid,ntemp,0)
 GO TO 330
 
!     CLOSE EST RECORD AND RETURN FOR ANOTHER ELEMENT TYPE.
 
 400 CALL WRITE (est,0,0,1)
 GO TO 320
 
!     ALL ELEMENTS HAVE BEEN PROCESSED-- CLOSE FILES, WRITE TRAILER AND
!     EXIT
 
 500 CALL CLOSE (scr1,clsrew)
 CALL CLOSE (est,clsrew)
 CALL CLOSE (gptt,clsrew)
 buf(1) = est
 buf(2) = nosimp
 IF (nogox /= 0) nogo = 1
 IF (nogo  /= 0) CALL mesage (-61,0,0)
 DO  i = 3,7
   buf(i) = 0
 END DO
 
!     PROCESS /TA1ACM/ LOAD EST TRAILER WITH FLAGS
!     TO THE USED /OPEN CORE/ BLOCKS
 
 DO  i = 1,nelem
   IF (ig(i) >= 0) CYCLE
   k = ig(i)
   DO  j = i,nelem
     IF (ig(j) == k) ig(j) = -ig(j)
   END DO
   j = ig(i)
   IF (j > 48) CALL mesage (-61,i,j)
   k = (j-1)/16
   j = j - k*16
   buf(k+5) = buf(k+5) + ktwo(j+16)
 END DO
 CALL wrttrl (buf)
 540 RETURN
 
!     **************************************************
 
!     INTERNAL BINARY SEARCH ROUTINE
 
 600 klo = 1
 khi = kn
 610 k   = (klo+khi+1)/2
 620 kx  = (k-1)*mx + itabl
 IF (id-z(kx+1) < 0) THEN
   GO TO   630
 ELSE IF (id-z(kx+1) == 0) THEN
   GO TO    90
 ELSE
   GO TO   640
 END IF
 630 khi = k
 GO TO 650
 640 klo = k
 650 IF (khi-klo-1  < 0) THEN
   GO TO   690
 ELSE IF (khi-klo-1  == 0) THEN
   GO TO   660
 ELSE
   GO TO   610
 END IF
 660 IF (k == klo) GO TO 670
 k   = klo
 GO TO 680
 670 k   = khi
 680 klo = khi
 GO TO 620
 690 IF (q4t3 .AND. nopshl >= 1) GO TO 800
 GO TO 1140
 
!     **************************************************
 
!     PROCESSING FOR LAMINATED COMPOSITES
 
!     INTERNAL SUBROUTINE TO READ PCOMP, PCOMP1 AND PCOMP2 DATA INTO
!     CORE
 
 
!     INITIALIZE VARIABLES AND SET POINTERS
 
 700 npc    = 0
 npc1   = 0
 npc2   = 0
 typc   = 0
 typc1  = 0
 typc2  = 0
 n      = buf3 - ll
 
!     LOCATE PCOMP DATA AND READ INTO CORE
 
 ipc  = ll + 1
 CALL locate (*720,z(buf2),pcomp,flag)
 CALL READ   (*1110,*710,ept,z(ipc),n,0,npc)
 CALL mesage (-8,0,nam)
 710 IF (npc > 0) typc = 1
 n = n - npc
 
!     LOCATE PCOMP1 DATA AND READ INTO CORE
 
 720 ipc1 = ipc + npc
 CALL locate (*740,z(buf2),pcomp1,flag)
 CALL READ   (*1110,*730,ept,z(ipc1),n,0,npc1)
 CALL mesage (-8,0,nam)
 730 IF (npc1 > 0) typc1 = 1
 n = n - npc1
 
!     LOCATE PCOMP2 DATA AND READ INTO CORE
 
 740 ipc2 = ipc1 + npc1
 CALL locate (*760,z(buf2),pcomp2,flag)
 CALL READ   (*1110,*750,ept,z(ipc2),n,0,npc2)
 CALL mesage (-8,0,nam)
 750 IF (npc2 > 0) typc2 = 1
 
!     SET SIZE OF LPCOMP. NUMBER OF WORDS READ INTO CORE
 
 760 lpcomp = ipc2 + npc2
 IF (lpcomp-1 > ll) comps = -1
 
!     CHECK FOR NO PCOMP, PCOMP1 OR PCOMP2 DATA
!     SET NOPSHL TO 2 IF BOTH 'PCOMP' AND PSHELL DATA ARE PRESENT
 
 IF (nopshl == 1 .AND. comps == 1) GO TO 1130
 IF (nopshl == 0 .AND. comps == -1) nopshl = 2
 GO TO 80
 
!     ***************************************************************
 
!     INTERNAL SUBROUTINE TO LOCATE A PARTICULAR PROPERTY ID FROM THE
!     'PCOMP' DATA AND TO CONVERT THE DATA TO PSHELL DATA FORMAT
 
 800 CONTINUE
 
!     Z(LL+1) THRU Z(LPCOMP) CONTAIN 'PCOMP' DATA
 
!     SET POINTERS
 
 kpc    = 4
 kpc2   = 2
 LEN    = 0
 nlay   = 0
 eoeloc = 0
 pidloc = 0
 itype  =-1
 
!     SEARCH FOR PID IN PCOMP DATA
 
 IF (typc == 0) GO TO 850
 z(lpcomp+1) = ipc
 npcomp = 0
 n = 2
 
 lpc = ipc1 - 1
 DO  iip = ipc,lpc
   IF (z(iip) /= -1) CYCLE
   z(lpcomp+n  ) = iip
   z(lpcomp+n+1) = iip + 1
   IF (iip == lpc) z(lpcomp+n+1) = 0
   n = n + 2
   npcomp = npcomp + 1
 END DO
 IF (lpcomp+n-2 >= buf3) CALL mesage (-8,0,nam)
 
!     LOCATE PARTICULAR PID
 
 DO  iip = 1,npcomp
   eoeloc = z(lpcomp+2*iip  )
   pidloc = z(lpcomp+2*iip-1)
   IF (z(pidloc) == id) GO TO 840
 END DO
 GO TO 850
 
 840 LEN  = eoeloc - pidloc
 nlay = (LEN-8)/kpc
 itype= 0
 GO TO 940
 
!     SEARCH FOR PID IN PCOMP1 DATA
 
 850 IF (typc1 == 0) GO TO 890
 
 z(lpcomp+1) = ipc1
 npcomp = 0
 n = 2
 
 lpc1 = ipc2 - 1
 DO  iip1 = ipc1,lpc1
   IF (z(iip1) /= -1) CYCLE
   z(lpcomp+n  ) = iip1
   z(lpcomp+n+1) = iip1 + 1
   IF (iip1 == lpc1) z(lpcomp+n+1) = 0
   npcomp = npcomp + 1
   n = n + 2
 END DO
 IF (lpcomp+n-2 >= buf3) CALL mesage (-8,0,nam)
 
!     LOCATE PARTICULAR PID
 
 DO  iip1 = 1,npcomp
   eoeloc = z(lpcomp+2*iip1  )
   pidloc = z(lpcomp+2*iip1-1)
   IF (z(pidloc) == id) GO TO 880
 END DO
 GO TO 890
 
 880 LEN  = eoeloc - pidloc
 nlay = LEN - 8
 itype= 1
 GO TO 940
 
!     SEARCH FOR PID IN PCOMP2 DATA
 
 890 IF (typc2 == 0) GO TO 930
 
 z(lpcomp+1) = ipc2
 npcomp = 0
 n = 2
 
 lpc2 = lpcomp - 1
 DO  iip2 = ipc2,lpc2
   IF (z(iip2) /= -1) CYCLE
   z(lpcomp+n  ) = iip2
   z(lpcomp+n+1) = iip2 + 1
   IF (iip2 == lpc2) z(lpcomp+n+1) = 0
   npcomp = npcomp + 1
   n = n + 2
 END DO
 IF (lpcomp+n-2 >= buf3) CALL mesage (-8,0,nam)
 
!     LOCATE PARTICULAR PID
 
 DO  iip2 = 1,npcomp
   eoeloc = z(lpcomp+2*iip2  )
   pidloc = z(lpcomp+2*iip2-1)
   IF (z(pidloc) == id) GO TO 920
 END DO
 GO TO 930
 
 920 LEN  = eoeloc - pidloc
 nlay = (LEN-8)/kpc2
 itype= 2
 GO TO 940
 
!     CHECK IF PID HAS BEEN FOUND IN 'PCOMP' DATA
 
 930 IF (itype < 0) GO TO 1140
 
!     DETERMINE DATA TO BE WRITTEN IN THE FORM OF PSHELL AND
!     WRITE TO SCR1
 
!     ITYPE  = 0,  PCOMP  ENTRY
!            = 1,  PCOMP1 ENTRY
!            = 2,  PCOMP2 ENTRY
 
!     CALCULATE LAMINATE THICKNESS - TLAM
 
 940 tlam = 0.
 
!     NOTE - IF Z(PIDLOC+7) IS EQUAL TO SYM OR SYMMEM, THE OPTION
!            TO MODEL EITHER A SYMMETRICAL OR SYMMETRICAL-MEMBRANE
!            LAMINATE HAS BEEN EXERCISED.  THEREFORE, THE TOTAL
!            THICKNESS IS TLAM = 2.0*TLAM
 
!     SET LAMOPT
 
 lamopt = z(pidloc+7)
 
!     PCOMP DATA
 
 IF (itype > 0) GO TO 960
 DO  k = 1,nlay
   ii   = (pidloc+5) + 4*k
   tlam = tlam + zz(ii)
 END DO
 IF (lamopt == sym .OR. lamopt == symmem) tlam = 2.0*tlam
 GO TO 1000
 
!     PCOMP1 DATA
 
 960 IF (itype > 1) GO TO 970
 ii   = pidloc + 6
 tlam = zz(ii)*nlay
 IF (lamopt == sym .OR. lamopt == symmem) tlam = 2.0*tlam
 GO TO 1000
 
!     PCOMP2 DATA
 
 970 DO  k = 1,nlay
   ii   = (pidloc+6) + 2*k
   tlam = tlam + zz(ii)
 END DO
 IF (lamopt == sym .OR. lamopt == symmem) tlam = 2.0*tlam
 
 
!     CREATE NEW PSHELL DATA AND WRITE TO ARRAY IPSHEL
!     NOTE - PID IS NOT WRITTEN TO IPSHEL
 
!     IPSHEL DATA TO BE WRITTEN TO SCR1
!     ============================================================
!     IPSHEL( 1)     = MID1     MEMBRANE MATERIAL
!     IPSHEL( 2)     = T        DEFAULT MEMBRANE THICKNESS
!     IPSHEL( 3)     = MID2     BENDING MATERIAL
!     IPSHEL( 4)     = 12I/T**3 BENDING STIFFNESS PARAMETER
!     IPSHEL( 5)     = MID3     TRANVERSE SHEAR MATERIAL
!     IPSHEL( 6)     = TS/T     SHEAR THICKNESS FACTOR
!     IPSHEL( 7)     = NSM      NON-STRUCTURAL MASS
!     IPSHEL(8,9)    = Z1,Z2    FIBRE DISTANCES
!     IPSHEL(10)     = MID4     MEMBRANE-BENDING COUPLING MATERIAL
!     IPSHEL(11)     = MCSID OR THETAM   //DATA FROM PSHELL
!     IPSHEL(12)     = FLAGM               OVERRIDDEN BY EST(18-19)//
!     IPSHEL(13)     = INTEGRATION ORDER (SET TO 0)
!                     (THE INTEGRATION ORDER IS NOT USED IN THE PROGRAM,
!                      BUT THIS WORD IS REQUIRED BECAUSE OF THE DESIGN
!                      OF THE EST DATA FOR THE CQUAD4/TRIA3 ELEMENTS.)
!     IPSHEL(14)     = SCSID OR THETAS   //DATA FROM PSHELL
!     IPSHEL(15)     = FLAGS               OVERRIDDEN BY EST(20-21)//
!     IPSHEL(16)     = ZOFF
 
!     CALCULATE ZOFFS
 
 1000 IF (z(pidloc+1) /= 0) zoffs = zz(pidloc+1) + 0.5*tlam
 IF (z(pidloc+1) == 0) zoffs = 0.0
 IF (ABS(zoffs)  <= 0.001) zoffs = 0.0
 
!     SET POINTER TO INDICATE NEW PSHELL DATA CREATED
 
 npshel = 1
 
!     INITIALIZE IPSHEL ARRAY
 
 DO  kk = 1,16
   ipshel(kk) = 0
 END DO
 
 rpshel( 4) = 1.0
 rpshel( 6) = 1.0
 
 ipshel( 1) = id + 100000000
 rpshel( 2) = tlam
 IF (lamopt == mem .OR. lamopt == symmem) GO TO 1020
 ipshel( 3) = id + 200000000
 ipshel( 5) = id + 300000000
 1020 rpshel( 7) = zz(pidloc+2)
 rpshel( 9) = 0.5*tlam
 rpshel( 8) =-rpshel(9)
 IF (lamopt /= sym .AND. lamopt /= mem .AND. lamopt /= symmem)  &
     ipshel(10) = id + 400000000
 ipshel(13) = 0
 rpshel(16) = zoffs
 
!     DO NOT WRITE TO OUTPUT FILE IF PREVIOUS ID IS SAME AS NEW ID.
!     OTHERWISE, WRITE THE NEWLY CREATED PSHELL BULK DATA ENTRY TO
!     OUTPUT FILE IF DIAG 40 IS TURNED ON
 
 IF (oldid == id) GO TO 1060
 IF (  .NOT.frstim) GO TO 1040
 frstim = .false.
 IF (l40 == 0) GO TO 1060
!WKBR CALL PAGE (3)
 CALL page2 (3)
 WRITE  (nout,1030)
 1030 FORMAT (//9X,'THE INPUT PCOMP, PCOMP1 OR PCOMP2 BULK DATA',  &
     ' ENTRIES HAVE BEEN REPLACED BY THE FOLLOWING PSHELL',  &
     ' AND MAT2 ENTRIES.',//)
 1040 IF (l40 == 0) GO TO 1060
 WRITE (nout,1050) id,ipshel( 1),rpshel( 2),ipshel( 3),  &
     rpshel( 4),ipshel( 5),rpshel( 6), rpshel( 7),rpshel( 8),rpshel( 9),  &
     ipshel(10),rpshel(11),rpshel(14), rpshel(16)
 1050 FORMAT (' PSHELL',i14,i12,1X,1P,e11.4,i12,1X,1P,e11.4,i12,  &
     2(1X,1P,e11.4), /9X,2(1X,1P,e11.4),i12,2(1X,f11.1),1X, 1P,e11.4)
 
!     SET OLDID TO ID
 
 1060 oldid = id
 GO TO 90
 
!     FATAL ERROR MESSAGES
 
 1100 j = -1
 GO TO 1150
 1110 j = -2
 GO TO 1150
 1120 j = -3
 GO TO 1150
 1130 buf(1) = elem(i  )
 buf(2) = elem(i+1)
 nogox  = 1
 CALL mesage (30,11,buf)
 kx = itabl
 GO TO 30
 1140 ksavew = buf(3)
 buf(3) = id
 nogo = 1
 CALL mesage (30,10,buf(2))
 kx = itabl
 buf(3) = ksavew
 GO TO 90
 1150 CALL mesage (j,FILE,nam)
 1160 buf(1) = tempid
 buf(2) = 0
 CALL mesage (-30,44,buf)
 RETURN
 
!     ARRAY IG IS FIRST DIMENSIONED IN TA1ABD
 
 1190 WRITE  (nout,1200) sfm
 1200 FORMAT (a25,', IG ARRAY IN TA1A TOO SMALL')
 CALL mesage (-61,0,0)
 
END SUBROUTINE ta1a
