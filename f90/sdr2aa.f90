SUBROUTINE sdr2aa
     
!     SDR2AA PROCESSES THE CASE CONTROL AND XYCDB DATA BLOCKS. IF XYCDB
!     IS PURGED, NO ACTION IS TAKEN. OTHERWISE, OUTPUT REQUESTS IN
!     CASE CONTROL ARE COMPARED WITH XY REQUESTS IN XYCDB. FOR EACH
!     SUBCASE AND EACH REQUEST TYPE, CASE CONTROL IS MODIFIED TO
!     REFLECT THE UNION OF THE REQUESTS. THE NEW CASE CONTROL IS
!     WRITTEN ON A SCRATCH FILE AND THE POINTER TO CASE CONTROL SWITCHED
 
 INTEGER :: tab   ,sdr2x1,buf   ,casecc,xycdb ,scr3  ,z     ,  &
     app   ,rd    ,rdrew ,wrt   ,wrtrew,clsrew,sysbuf,  &
     xsetno,buf1  ,buf2  ,buf3  ,subcse,anynew,FILE  ,  &
     dbname,setno ,arg   ,esta  ,xycdbf,trn   ,frq   , cei   ,formt ,sort2
 DIMENSION       sdr2x1(1)    ,tab(14)      ,buf(10)       ,nam(2)
 COMMON /sdr2x1/ sdr2x1,ieldef,itload,isymfl,iloads,idispl,istr  ,  &
     ielf  ,iacc  ,ivel  ,ispcf ,ittl  ,ilsym ,ifrout, isload,idload
 COMMON /sdr2x2/ casecc,cstm  ,mpt   ,dit   ,eqexin,sil   ,gptt  ,  &
     edt   ,bgpdt ,pg    ,qg    ,ugv   ,est   ,phig  ,  &
     eigr  ,opg1  ,oqg1  ,ougv1 ,oes1  ,oef1  ,pugv1 ,  &
     oeigr ,ophig ,pphig ,esta  ,gptta ,harms ,xycdb , scr3
 COMMON /sdr2x4/ x4(72),frq(2),trn(2),bkl(4),cei(2)
 COMMON /BLANK / app(2),sort2
 COMMON /zzzzzz/ z(1)
 COMMON /system/ sysbuf
 COMMON /names / rd    ,rdrew ,wrt   ,wrtrew,clsrew
 DATA    tab   / 1,    6, 2,    10,  &
     3,    9, 4,    11,  &
     5,    5, 6,    7,  &
     7,    8      /, xsetno/ 100000000    /,  &
     nam   / 4HSDR2,4HAA  /
 
!     SET BUFFER POINTERS AND PERFORM GENERAL INITIALIZATION.
 
 buf1  = korsz(z) - sysbuf
 buf2  = buf1 - sysbuf
 buf3  = buf2 - sysbuf
 imstr = 1
 master= 1
 lastxy= 0
 anynew= 0
 sort2 =-1
 
!     OPEN XYCDB. IF PURGED, RETURN.
 
 CALL OPEN (*1034,xycdb,z(buf1),rdrew)
 FILE = xycdb
 CALL fwdrec (*1035,xycdb)
 CALL fwdrec (*1035,xycdb)
 
!     READ FIRST LINE OF XYCDB. IF SUBCASE = 0 (MEANING DATA APPLIES
!     TO ALL SUBCASES), READ IN DATA FOR ZERO SUBCASE.
 
 last   = 0
 xycdbf = xycdb
 CALL READ (*1035,*1035,xycdb,buf,6,0,flag)
 sort2  = 0
 subcse = buf(1)
 IF (subcse /= 0) GO TO 1013
 i      = imstr
 1011 z(i  ) = buf(2)
 z(i+1) = buf(3)
 i      = i + 2
 CALL READ (*2002,*1012,xycdb,buf,6,0,flag)
 IF (buf(1) == 0) GO TO 1011
 nmstr  = i - 2
 ixysc  = i
 GO TO 1019
 
!     HERE IF MASTER SUBCASE IS THE ONLY SUBCASE IN XYCDB.
 
 1012 nmstr  = i - 2
 nxysc  = nmstr
 master = 0
 lastxy = 1
 
!     REDUCE LIST TO UNIQUE PAIRS
 
 IF (imstr == nmstr) GO TO 1019
 nmstr  = nmstr - 2
 j      = imstr
 DO  i = imstr,nmstr,2
   IF (z(i+2) == z(j) .AND. z(i+3) == z(j+1)) CYCLE
   z(j+2) = z(i+2)
   z(j+3) = z(i+3)
   j      = j + 2
 END DO
 nmstr  = j
 nxysc  = nmstr
 GO TO 1019
 
!     HERE IF NO MASTER SUBCASE -- CREATE A DUMMY MASTER.
 
 1013 nmstr  = imstr
 ixysc  = imstr + 2
 z(imstr  ) = 9999
 z(imstr+1) = 0
 master = -1
 GO TO 1019
 
!     OPEN CASE CONTROL AND SCRATCH FILE FOR MODIFIED CASE CONTROL
 
 1019 CALL gopen (casecc,z(buf2),rdrew)
 FILE = scr3
 CALL OPEN (*2001,scr3,z(buf3),wrtrew)
 CALL fname (casecc,buf(9))
 CALL WRITE (scr3,buf(9),2,1)
 
!     READ DATA FOR ONE SUBCASE. STORE DATA BLOCK AND ID IN OPEN CORE.
 
 1020 IF (master == 0 .OR. lastxy /= 0) GO TO 1030
 subcse = buf(1)
 i      = ixysc
 1021 z(i  ) = buf(2)
 z(i+1) = buf(3)
 i      = i + 2
 CALL READ (*1035,*1023,xycdbf,buf,6,0,flag)
 IF (buf(1) == subcse) GO TO 1021
 GO TO 1025
 1023 lastxy = 1
 
!     COPY DATA FROM MASTER SUBCASE AFTER CURRENT SUBCASE.
!     THEN SORT DATA TOGETHER TO FORM SORTED UNION.
 
 1025 DO  j = imstr,nmstr,2
   z(i  ) = z(j  )
   z(i+1) = z(j+1)
   i      = i + 2
 END DO
 n = i - ixysc
 CALL sort (0,0,2,-2,z(ixysc),n)
 CALL sort (0,0,2,-1,z(ixysc),n)
 
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
 CALL READ (*1035,*1031,casecc,z(icc+1),buf3-icc,1,ncc)
 CALL mesage (-8,0,nam)
 1031 IF (subcse == z(icc+1)) GO TO 10311
 IF (master /=       -1) GO TO 1032
 IF (subcse > z(icc+1)) GO TO 1030
 IF (lastxy ==        0) GO TO 1020
 IF (anynew ==        0) GO TO 1035
 CALL WRITE (scr3,z(icc+1),ncc,1)
 GO TO 1030
 10311 ixy = ixysc
 nxy = nxysc
 GO TO 1040
 1032 ixy = imstr
 nxy = nmstr
 GO TO 1040
 
!     TERMINATE PROCESSING.
 
 1035 CALL CLOSE (casecc,clsrew)
 CALL CLOSE (xycdbf,clsrew)
 CALL CLOSE (scr3  ,clsrew)
 IF (anynew /= 0) casecc = scr3
 1034 RETURN
 
!     PICK UP POINTER TO CURRENT OUTPUT REQUEST.
!     DETERMINE IF XYCDB REQUEST EXISTS.
 
 1040 loop   = 1
 1041 dbname = tab(loop  )
 ix     = tab(loop+1)
 ireq   = icc + sdr2x1(ix)
 setno  = z(ireq)
 DO  j = ixy,nxy,2
   IF (z(j) == dbname) GO TO 1043
 END DO
 GO TO 1100
 1043 ixyset = j
 DO  j = ixyset,nxy,2
   IF (z(j) /= dbname) GO TO 1045
 END DO
 nxyset = nxy
 GO TO 1050
 1045 nxyset = j - 2
 
!     BRANCH ON CASECC REQUEST - NOTE, NO ACTION IF REQUEST = ALL.
 
 1050 IF (setno < 0.0) THEN
   GO TO  1098
 ELSE IF (setno == 0.0) THEN
   GO TO  1060
 ELSE
   GO TO  1070
 END IF
 
!     HERE IF NO CASECC REQUEST.
!     BUILD XYCDB SET IN CASECC SET FORMAT. ADD SET TO
!     CASECC RECORD AND TURN ON CASECC REQUEST FOR SET.
 
 1060 xsetno = xsetno + 1
 z(ireq  ) = xsetno
 z(ireq+1) = 0
 formt  = -2
 IF (app(1) == trn(1)) formt = -1
 z(ireq+2) = formt
 sort2  = 0
 ix     = icc + ncc + 1
 z(ix)  = xsetno
 jx     = ix + 2
 z(jx)  = z(ixyset+1)
 IF (ixyset == nxyset) GO TO 1066
 ixyset = ixyset + 2
 n      = 1
 DO  j = ixyset,nxyset,2
   IF (z(j+1)-z(jx) == n) GO TO 1064
   IF (n /= 1) GO TO 1062
   jx     = jx + 1
   z(jx)  = z(j+1)
   CYCLE
   1062 z(jx+1)=-z(j-1)
   jx     = jx + 2
   z(jx)  = z(j+1)
   n      = 1
   CYCLE
   1064 n      = n + 1
 END DO
 IF (n == 1) GO TO 1066
 jx     = jx + 1
 z(jx  )= -z(nxyset+1)
 1066 z(ix+1)= jx - ix - 1
 ncc    = ncc + z(ix+1) + 2
 anynew = 1
 GO TO 1100
 
!     HERE IF CASECC SET AND XYCDB SET EXIST.
!     FIRST, LOCATE CASECC SET.
 
 1070 ilist  = icc + ncc + 3
 ix     = icc + ilsym
 isetno = ix  + z(ix) + 1
 1071 iset   = isetno + 2
 nset   = z(isetno+1) + iset - 1
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
 IF (arg-z(i) < 0.0) THEN
   GO TO  1088
 ELSE IF (arg-z(i) == 0.0) THEN
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
 z(ireq)   = xsetno
 z(ireq+1) = 10*setno + z(ireq+1)
 z(ireq+2) = -IABS(z(ireq+2))
 sort2     = 0
 z(ilist-2)= xsetno
 z(ilist-1)= n
 ncc       = ncc + n + 2
 anynew    = 1
 GO TO 1100
 
!     HERE IF CASECC SET = ALL AND XY REQUEST EXISTS - TURN SORT2 ON.
 
 1098 z(ireq+2) = -IABS(z(ireq+2))
 sort2     = 0
 
!     TEST FOR COMPLETION OF ALL CASECC REQUESTS FOR CURRENT SUBCASE.
!     WHEN COMPLETE, WRITE CURRENT SUBCASE ON SCRATCH FILE.
 
 1100 loop = loop + 2
 IF (loop <= 13) GO TO 1041
 CALL WRITE (scr3,z(icc+1),ncc,1)
 
!     RETURN TO READ ANOTHER RECORD IN CASE CONTROL OR ANOTHER XYCDB
!     SUBCASE
 
 IF (master ==        0) GO TO 1030
 IF (subcse <= z(icc+1)) GO TO 1020
 GO TO 1030
 
!     FATAL FILE ERRORS
 
 2000 CALL mesage (n,FILE,nam)
 2001 n = -1
 GO TO 2000
 2002 n = -2
 GO TO 2000
END SUBROUTINE sdr2aa
