SUBROUTINE vdra
     
!     VDRA PROCESSES THE CASE CONTROL AND XYCDB DATA BLOCKS. IF XYCDB
!     IS PURGED, NO ACTION IS TAKEN. OTHERWISE, OUTPUT REQUESTS IN
!     CASE CONTROL ARE COMPARED WITH XY REQUESTS IN XYCDB. FOR EACH
!     SUBCASE AND EACH REQUEST TYPE, CASE CONTROL IS MODIFIED TO REFLECT
!     THE UNION OF THE REQUESTS. THE NEW CASE CONTROL IS WRITTEN ON A
!     SCRATCH FILE AND THE POINTER TO CASE CONTROL SWITCHED.
 
 INTEGER :: buf   ,casecc,xycdb ,scr1  ,scr3  ,z     ,app   ,  &
     rd    ,rdrew ,wrt   ,wrtrew,clsrew,sysbuf,xsetno,  &
     buf1  ,buf2  ,buf3  ,subcse,anynew,FILE  ,dbname,  &
     setno ,arg   ,sdr2  ,xset0 ,vdrcom,xycdbf,vdrreq, trn   ,FORMAT,frq   ,sort2
 DIMENSION       nam(2),buf(50)      ,masks(6)     ,cei(2),frq(2),  &
     trn(2),modal(2)     ,DIRECT(2)    ,vdrcom(1)
 COMMON /vdrcom/ vdrcom,idisp ,ivel  ,iacc  ,ispcf ,iloads,istr  ,  &
     ielf  ,iadisp,iavel ,iaacc ,ipnl  ,ittl  ,ilsym ,  &
     ifrout,idload,casecc,eqdyn ,usetd ,infile,oeigs ,  &
     pp    ,xycdb ,pnl   ,outfle,opnl1 ,scr1  ,scr3  ,  &
     buf1  ,buf2  ,buf3  ,nam   ,buf   ,masks ,cei   ,  &
     frq   ,trn   ,DIRECT,xset0 ,vdrreq,modal
 COMMON /zzzzzz/ z(1)
 COMMON /BLANK / app(2),FORM(2),sort2,output,sdr2  ,imode
 COMMON /names / rd    ,rdrew ,wrt   ,wrtrew,clsrew
 COMMON /system/ sysbuf
 
!     SET BUFFER POINTERS AND PERFORM GENERAL INITIALIZATION.
 
 buf1   = korsz(z) - sysbuf
 buf2   = buf1 - sysbuf
 buf3   = buf2 - sysbuf
 ixy    = 1
 lastxy = 0
 anynew = 0
 sdr2   =-1
 vdrreq = 0
 xsetno = xset0
 imstr  = 1
 master = 1
 
!     OPEN XYCDB. IF PURGED, RETURN.
 
 CALL OPEN (*1036,xycdb,z(buf1),rdrew)
 FILE = xycdb
 CALL fwdrec (*1036,xycdb)
 CALL fwdrec (*1036,xycdb)
 
!     READ FIRST LINE OF XYCDB. IF SUBCASE = 0 (MEANING DATA APPLIES
!     TO ALL SUBCASES), READ IN DATA FOR ZERO SUBCASE.
 
 last   = 0
 xycdbf = xycdb
 CALL READ (*1035,*1035,xycdb,buf,6,0,flag)
 subcse = buf(1)
 IF (subcse /= 0) GO TO 1013
 i = imstr
 1011 z(i  ) = buf(2)
 z(i+1) = buf(3)
 i = i + 2
 CALL READ (*2002,*1012,xycdb,buf,6,0,flag)
 IF (buf(1) == 0) GO TO 1011
 nmstr = i - 2
 ixysc = i
 GO TO 1019
 
!     HERE IF MASTER SUBCASE IS THE ONLY SUBCASE IN XYCDB.
 
 1012 nmstr  = i - 2
 nxysc  = nmstr
 master = 0
 lastxy = 1
 
!     REDUCE LIST TO UNIQUE PAIRS
 
 IF (imstr == nmstr) GO TO 1019
 nmstr = nmstr - 2
 j = imstr
 DO  i = imstr,nmstr,2
   IF (z(i+2) == z(j) .AND. z(i+3) == z(j+1)) CYCLE
   z(j+2) = z(i+2)
   z(j+3) = z(i+3)
   j = j + 2
 END DO
 nmstr = j
 nxysc = nmstr
 GO TO 1019
 
!     HERE IF NO MASTER SUBCASE -- CREATE A DUMMY MASTER.
 
 1013 nmstr = imstr
 ixysc = imstr + 2
 z(imstr  ) = 9999
 z(imstr+1) = 0
 
!     OPEN CASE CONTROL AND SCRATCH FILE FOR MODIFIED CASE CONTROL
 
 1019 CALL gopen (casecc,z(buf2),0)
 CALL gopen (scr3,z(buf3),1)
 
!     READ DATA FOR ONE SUBCASE. STORE DATA BLOCK AND ID IN OPEN CORE.
 
 1020 IF (master == 0 .OR. lastxy /= 0) GO TO 1030
 subcse = buf(1)
 i = ixysc
 1021 z(i  ) = buf(2)
 z(i+1) = buf(3)
 i = i + 2
 CALL READ (*1035,*1023,xycdbf,buf,6,0,flag)
 IF (buf(1) == subcse) GO TO 1021
 GO TO 1025
 1023 lastxy = 1
 
!     COPY DATA FROM MASTER SUBCASE AFTER CURRENT SUBCASE.
!     THEN SORT DATA TOGETHER TO FORM SORTED UNION.
 
 1025 DO  j = imstr,nmstr,2
   z(i  ) = z(j  )
   z(i+1) = z(j+1)
   i = i + 2
 END DO
 n = i - ixysc
 CALL sort2k (0,0,2,1,z(ixysc),n)
 
!     REDUCE LIST TO UNIQUE PAIRS.
 
 nxysc = i - 4
 j = ixysc
 DO  i = ixysc,nxysc,2
   IF (z(i+2) == z(j) .AND. z(i+3) == z(j+1)) CYCLE
   z(j+2) = z(i+2)
   z(j+3) = z(i+3)
   j = j + 2
 END DO
 nxysc = j
 
!     READ A RECORD IN CASE CONTROL. SET POINTERS FOR XYCDB DATA TO
!     EITHER MASTER SUBCASE OR CURRENT SUBCASE IN CORE.
 
 1030 icc = nxysc + 1
 CALL READ (*1033,*1031,casecc,z(icc+1),buf3-icc,1,ncc)
 CALL mesage (-8,0,nam)
 1031 IF (master == 0 .OR. z(icc+1) /= subcse) GO TO 1032
 ixy = ixysc
 nxy = nxysc
 GO TO 1040
 1032 ixy = imstr
 nxy = nmstr
 GO TO 1040
 
!     TERMINATE PROCESSING.
 
 1033 CONTINUE
 1035 CALL CLOSE (casecc,clsrew)
 CALL CLOSE (xycdbf,clsrew)
 CALL CLOSE (scr3  ,clsrew)
 IF (anynew /= 0) casecc = scr3
 RETURN
 
 1036 vdrreq = 1
 CALL CLOSE (xycdb,clsrew)
 RETURN
 
!     PICK UP POINTER TO CURRENT OUTPUT REQUEST.
!     DETERMINE IF XYCDB REQUEST EXISTS.
 
 1040 loop   = 1
 1041 dbname = loop
 ireq   = icc + vdrcom(loop+1)
 setno  = z(ireq)
 DO  j = ixy,nxy,2
   IF (z(j) == dbname) GO TO 1043
 END DO
 GO TO 1095
 1043 ixyset = j
 DO  j = ixyset,nxy,2
   IF (z(j) /= dbname) GO TO 1045
 END DO
 nxyset = nxy
 GO TO 1050
 1045 nxyset = j - 2
 
!     BRANCH ON CASECC REQUEST-- NOTE, NO ACTION IF REQUEST = ALL.
 
 1050 IF (loop > 7) GO TO 1051
 sdr2 = +1
 GO TO 1100
 1051 vdrreq = 1
 sort2  =+1
 IF (setno < 0.0) THEN
   GO TO  1098
 ELSE IF (setno == 0.0) THEN
   GO TO  1060
 ELSE
   GO TO  1070
 END IF
 
!     HERE IF NO CASECC REQUEST.
!     BUILD XYCDB SET IN CASECC SET FORMAT. ADD SET TO
!     CASECC RECORD AND TURN ON CASECC REQUEST FOR SET.
 
 1060 xsetno  = xsetno + 1
 z(ireq) = xsetno
 z(ireq+1) = 0
 FORMAT  = -2
 IF (app(1) == trn(1)) FORMAT = -1
 z(ireq+2) = FORMAT
 ix = icc + ncc + 1
 z(ix) = xsetno
 jx = ix + 2
 z(jx) = z(ixyset+1)
 IF (ixyset == nxyset) GO TO 1066
 ixyset = ixyset + 2
 n = 1
 DO  j = ixyset,nxyset,2
   IF (z(j+1)-z(jx) == n) GO TO 1064
   IF (n /= 1) GO TO 1062
   jx = jx + 1
   z(jx)= z(j+1)
   CYCLE
   1062 z(jx+1) = -z(j-1)
   jx = jx + 2
   z(jx) = z(j+1)
   n = 1
   CYCLE
   1064 n = n + 1
 END DO
 IF (n == 1) GO TO 1066
 jx = jx + 1
 z(jx  ) = -z(nxyset+1)
 1066 z(ix+1) = jx - ix - 1
 ncc = ncc + z(ix+1) + 2
 anynew = 1
 GO TO 1100
 
!     HERE IF CASECC SET AND XYCDB SET EXIST.
!     FIRST, LOCATE CASECC SET.
 
 1070 ilist = icc + ncc + 3
 ix    = icc + ilsym
 isetno= ix  + z(ix) + 1
 1071 iset  = isetno + 2
 nset  = z(isetno+1) + iset - 1
 IF (z(isetno) == setno) GO TO 1080
 isetno = nset + 1
 IF (isetno < ilist) GO TO 1071
 GO TO 1100
 
!     COMPARE EACH POINT IN XYCDB REQUEST WITH CASECC SET.
!     ADD ANY POINTS IN XYCDB NOT IN CASECC TO CASECC SET.
 
 1080 i = iset
 j = ixyset
 k = ilist
 l = iset
 1081 arg = z(j+1)
 1082 IF (i-nset < 0) THEN
   GO TO  1083
 ELSE IF (i-nset == 0) THEN
   GO TO  1085
 ELSE
   GO TO  1088
 END IF
 1083 IF (z(i+1) > 0) GO TO 1085
 n = 2
 IF (arg-z(i  ) < 0.0) THEN
   GO TO  1088
 ELSE IF (arg-z(i  ) == 0.0) THEN
   GO TO  1091
 END IF
 1084 IF (arg+z(i+1) < 0.0) THEN
   GO TO  1091
 ELSE IF (arg+z(i+1) == 0.0) THEN
   GO TO  1087
 ELSE
   GO TO  1086
 END IF
 1085 n = 1
 IF (arg-z(i  ) < 0.0) THEN
   GO TO  1088
 ELSE IF (arg-z(i  ) == 0.0) THEN
   GO TO  1087
 END IF
 1086 i = i + n
 GO TO 1082
 1087 i = i + n
 GO TO 1091
 1088 IF (l == i) GO TO 1090
 ln = i - 1
 ll = l
 DO  l = ll,ln
   z(k) = z(l)
   k = k + 1
 END DO
 l = i
 1090 z(k) = arg
 k = k + 1
 1091 j = j + 2
 IF (j <= nxyset) GO TO 1081
 n = k - ilist
 IF (n ==    0) GO TO 1100
 IF (l > nset) GO TO 1094
 DO  ll = l,nset
   z(k) = z(ll)
   k = k + 1
 END DO
 n = k - ilist
 
!     IF NO NEW POINTS IN SET, CURRENT CASECC SET IS UNION.
!     OTHERWISE, NEW SET IS UNION. TURN ON REQUEST FOR IT AND
!     EXTEND END OF CASECC RECORD.
 
 1094 xsetno    = xsetno + 1
 z(ireq  ) = xsetno
 z(ireq+1) = 10*setno + z(ireq+1)
 z(ireq+2) =-IABS(z(ireq+2))
 z(ilist-2)= xsetno
 z(ilist-1)= n
 ncc       = ncc + n + 2
 anynew    = 1
 GO TO 1100
 
!     HERE IF NO XYCDB REQUEST EXISTS.
 
 1095 IF (setno == 0) GO TO 1100
 IF (loop  > 7) GO TO 1096
 sdr2 = 1
 GO TO 1100
 1096 vdrreq = 1
 GO TO 1100
 
!     HERE IF CASECC SET = ALL AND XY REQUEST EXISTS - TURN SORT 2 ON.
 
 1098 z(ireq+2) = -IABS(z(ireq+2))
 
!     TEST FOR COMPLETION OF ALL CASECC REQUESTS FOR CURRENT SUBCASE.
!     WHEN COMPLETE, WRITE CURRENT SUBCASE ON SCRATCH FILE.
 
 1100 loop = loop + 1
 IF (loop <= 11) GO TO 1041
 CALL WRITE (scr3,z(icc+1),ncc,1)
 
!     RETURN TO READ ANOTHER RECORD IN CASE CONTROL OR ANOTHER XYCDB
!     SUBCASE
 
 IF (master == 0) GO TO 1030
 IF (subcse <= z(icc+1)) GO TO 1020
 GO TO 1030
 
!     FATAL FILE ERROR
 
 2000 CALL mesage (n,FILE,nam)
 2002 n = -2
 GO TO 2000
END SUBROUTINE vdra
