SUBROUTINE ta1h
     
!     FOR LEVEL 16 A MAJOR REVISION HAS BEEN MADE TO TA1B. THE ECPT AND
!     GPCT ARE NO LONGER CONTSTRUCTED BUT, INSTEAD, THE GPECT IS BUILT.
!     THE GPECT IS ESSENTIALLY A TRUNCATED VERSION OF THE OLD ECPT. IT
!     CONTAINS ONE LOGICAL RECORD FOR EACH GRID OR SCALAR POINT IN THE
!     MODEL. EACH LOGICAL RECORD CONTAINS THE CONNECTION DATA FOR EACH
!     ELEMENT CONNECTED TO THE GRID POINT.
 
 EXTERNAL        andf
 INTEGER :: genl  ,ect   ,ept   ,bgpdt ,sil   ,gptt  ,cstm  ,  &
     est   ,gei   ,ecpt  ,gpct  ,scr1  ,scr2  ,scr3  ,  &
     scr4  ,z     ,sysbuf,tempid,elem  ,tempsz,elemid,  &
     outpt ,rd    ,rdrew ,wrt   ,wrtrew,clsrew,cls   ,  &
     buf   ,flag  ,buf1  ,buf2  ,buf3  ,op    ,two24 ,  &
     scri  ,scro  ,blk   ,oufile,andf  ,out(3),gpect , eqexin
 DIMENSION       buf(50)      ,bufr(50)     ,nam(2),blk(2),zz(1)
 COMMON /BLANK / luset ,nosimp,nosup ,nogenl,genl  ,comps
 COMMON /ta1com/ nsil  ,ect   ,ept   ,bgpdt ,sil   ,gptt  ,cstm  ,  &
     mpt   ,est   ,gei   ,gpect ,ecpt  ,gpct  ,mptx  ,  &
     pcomps,eptx  ,scr1  ,scr2  ,scr3  ,scr4  ,eqexin
 COMMON /system/ ksystm(65)
 COMMON /gpta1 / nelem ,jlast ,incr  ,elem(1)
 COMMON /names / rd    ,rdrew ,wrt   ,wrtrew,clsrew,cls
 COMMON /ta1ab / tempsz
!ZZ   COMMON /ZZTAA2/ Z(1)
 COMMON /zzzzzz/ z(20000)
 EQUIVALENCE     (ksystm( 1),sysbuf) ,(ksystm(2),outpt),  &
     (ksystm(10),tempid) ,(buf(1),bufr(1)) , (z(1),zz(1))        ,(blk(2),n)
 DATA    nam   / 4HTA1H, 4H    /     ,two24  / 4194304 /
 
!     PERFORM GENERAL INITIALIZATION
 
 n2   = 2*nelem - 1
 n21  = n2 + 1
 buf1 = korsz(z) - sysbuf - 2
 buf2 = buf1 - sysbuf
 buf3 = buf2 - sysbuf
!WKBR spr 93012      NEQ1 = NSIL + 1
 neq1 = nsil + 2
 neq2 = 0
 
!     THE GRID POINT COUNTER(GPC)HAS ONE ENTRY PER GRID OR SCALAR POINT
!     IN THE STRUCTURE. EACH ENTRY CONTAINS THE NUMBER OF STRUCTURAL
!     ELEMENTS CONNECTED TO THE POINT.
 
 DO  i = 1,nsil
   z(i+1) = 0
 END DO
 
!     OPEN THE ECT. INITIALIZE TO LOOP THRU BY ELEMENT TYPE.
 
 FILE = ect
 CALL gopen (ect,z(buf1),rdrew)
 noect = 1
 
!     IGNORE PLOTEL AND REACT ELEMENTS. OTHERWISE, LOCATE AN ELEMENT
!     TYPE. IF PRESENT, READ ALL ELEMENTS OF THAT TYPE AND INCREMENT
!     THE GPC ENTRY FOR EACH POINT TO WHICH THE ELEMENT IS CONNECTED.
 
 2012 CALL ectloc (*2026,ect,buf,i)
 noect = 0
 lx = elem(i+12)
 mm = lx+elem(i+9) - 1
 m  = elem(i+5)
 2021 CALL READ (*3201,*2012,ect,buf,m,0,flag)
 DO  l = lx,mm
   k  = buf(l)
   IF (k /= 0) z(k+1) = z(k+1) + 1
 END DO
 GO TO 2021
 2026 CONTINUE
 IF (noect /= 0) GO TO 3209
 
!     REPLACE ENTRIES IN THE GPC BY A RUNNING SUM
!     THUS CREATING POINTERS INTO ECPT0
!     QUEUE WARNING MESSAGES FOR GRID PTS. WITH NO ELEMENTS CONNECTED.
 
 z(1)  = 1
 maxel = 0
 DO  i = 1,nsil
   maxel = MAX0(maxel,z(i+1))
   IF (z(i+1) /= 0) GO TO 2037
   
   j = 0
   IF (neq2 < 0) THEN
     GO TO  2035
   ELSE IF (neq2 == 0) THEN
     GO TO  2031
   ELSE
     GO TO  2033
   END IF
   2031 neq2 = -1
   z(neq1) = eqexin
   CALL rdtrl (z(neq1))
   IF (z(neq1) <= 0) GO TO 2035
   FILE = eqexin
   CALL gopen (eqexin,z(buf2),rdrew)
   CALL READ (*3200,*2032,eqexin,z(neq1),buf3,1,neq2)
   2032 CALL CLOSE (eqexin,clsrew)
   CALL sort (0,0,2,2,z(neq1),neq2)
   2033 j = z((i-1)*2+neq1)
   
   2035 buf(1) = i
   buf(2) = j
   CALL mesage (30,15,buf)
   2037 z(i+1) = z(i) + z(i+1)
 END DO
 
!     DETERMINE BAND OF ENTRIES IN ECPT0 WHICH WILL FIT IN CORE
!     NDX1 = POINTER IN GPC TO 1ST  ENTRY FOR CURRENT PASS.
!     NDX2 = POINTER IN GPC TO LAST ENTRY FOR CURRENT PASS.
 
 ndx1 = 1
 ndx2 = nsil
 llx  = 1
 iecpt0 = nsil + 2
 length = buf1 - iecpt0
 op = wrtrew
 2042 IF (z(ndx2+1)-z(ndx1)+2 <= length) GO TO 2050
 ndx2 = ndx2 - 1
 GO TO 2042
 
!     PASS THE ECT. FOR EACH GRID PT IN RANGE ON THIS PASS,
!     STORE ELEMENT POINTER = 2**K * J + WORD POSITION IN ECT RECORD
!     WHERE K=22 FOR LEVEL 16 AND J = ENTRY NBR OF ELEMENT IN /GPTA1/
!     (WHICH IS SAME AS ELEMENT TYPE AS OF LEVEL 15)
 
 2050 FILE = ect
 CALL gopen (ect,z(buf1),rdrew)
 izero = z(ndx1)
 2051 CALL ectloc (*2059,ect,buf,i)
 j  = (i-1)/incr + 1
 idcntr = two24*j
 m  = elem(i+5)
 lx = elem(i+12)
 mm = lx + elem(i+9) - 1
 2052 CALL READ (*3201,*2051,ect,buf,m,0,flag)
 DO  l = lx,mm
   k  = buf(l)
   IF (k < ndx1 .OR. k > ndx2) CYCLE
   ix = z(k) - izero + iecpt0
   z(ix) = idcntr
   z(k)  = z(k) + 1
 END DO
 idcntr = idcntr + m
 GO TO 2052
 2059 CONTINUE
 
!     WRITE ECPT0 AND TEST FOR ADDITIONAL PASSES
!     ECPT0 CONTAINS ONE LOGICAL RECORD FOR EACH GRID OR SCALAR POINT.
!     EACH LOGICAL RECORD CONTAINS N PAIRS OF(-1,ELEMENT POINTER)WHERE
!     N= NUMBER OF ELEMENTS CONNECTED TO THE PIVOT.
!     IF NO ELEMENTS CONNECTED TO POINT, RECORD IS ONE WORD = 0.
 
 FILE = scr1
 CALL OPEN (*3200,scr1,z(buf1),op)
 elemid =  1
 buf(1) = -1
 lj = iecpt0 - 1
 DO  i = ndx1,ndx2
   m  = z(i) - llx
   IF (m /= 0) GO TO 2063
   CALL WRITE (scr1,0,1,1)
   GO TO 2062
   2063 DO  j = 1,m
     lj = lj + 1
     buf(2) = z(lj)
     CALL WRITE (scr1,buf,2,0)
   END DO
   CALL WRITE (scr1,0,0,1)
   2062 llx = z(i)
 END DO
 IF (ndx2 >= nsil) GO TO 2070
 CALL CLOSE (scr1,cls)
 ndx1 = ndx2 + 1
 ndx2 = nsil
 op   = wrt
 GO TO 2042
 
!     READ AS MUCH OF ECT AS CORE CAN HOLD
!     FIRST N21 CELLS OF CORE CONTAIN A POINTER TABLE WHICH HAS TWO
!     ENTRIES PER ELEMENT TYPE. 1ST ENTRY HAS POINTER TO 1ST WORD OF
!     ECT DATA IN CORE FOR AN ELEMENT TYPE  2ND ENTRY HAS WORD POSITION
!     IN ECT RECORD OF THAT TYPE FOR LAST ENTRY READ ON PREVIOUS PASS.
 
 2070 CALL CLOSE (scr1,clsrew)
 scri = scr1
 scro = scr2
 FILE = ect
 CALL gopen (ect,z(buf1),rdrew)
 DO  j = 1,n21
   z(j) = 0
 END DO
 l = n21 + 1
 2072 CALL ectloc (*2080,ect,buf,ielem)
 i = 2*((ielem-1)/incr + 1) - 1
 z(i) = l
 ll   = 0
 m    = elem(ielem+5)
 last = buf3-m
 2073 IF (l > last) GO TO 2080
 CALL READ (*3201,*2074,ect,z(l),m,0,flag)
 z(l)   = elemid
 elemid = elemid +1
 l  = l  + m
 ll = ll + m
 GO TO 2073
 2074 CONTINUE
 GO TO 2072
 
!     PASS ECPT0 ENTRIES LINE BY LINE
!     ATTACH EACH REFERENCED ECT ENTRY WHICH IS NOW IN CORE
 
 2080 FILE = scri
 CALL OPEN (*3200,scri,z(buf2),rdrew)
 CALL OPEN (*3200,scro,z(buf3),wrtrew)
 2082 CALL READ (*2090,*2086,scri,buf,1,0,flag)
 IF (buf(1) < 0.0) THEN
   GO TO  2083
 ELSE IF (buf(1) == 0.0) THEN
   GO TO  2087
 ELSE
   GO TO  2085
 END IF
 2083 CALL READ (*3201,*3202,scri,buf(2),1,0,flag)
 khr = buf(2)/two24
 ktwo24 = khr*two24
 k  = 2*khr - 1
 idptr = buf(2) - ktwo24
 kk = z(k) + idptr - z(k+1)
 IF (z(k) == 0 .OR. kk > last) GO TO 2084
 j  = (khr-1)*incr + 1
 mm = elem(j+5)
 buf(1) = mm
 buf(2) = andf(z(kk),two24-1) + ktwo24
 CALL WRITE (scro,buf,2,0)
 CALL WRITE (scro,z(kk+1),mm-1,0)
 GO TO 2082
 2084 CALL WRITE (scro,buf,2,0)
 GO TO 2082
 2085 CALL READ  (*3201,*3202,scri,buf(2),buf(1),0,flag)
 CALL WRITE (scro,buf,buf(1)+1,0)
 GO TO 2082
 2086 CALL WRITE (scro,0,0,1)
 GO TO 2082
 2087 CALL WRITE (scro,0,1,1)
 CALL fwdrec (*3201,scri)
 GO TO 2082
 
!     TEST FOR COMPLETION OF STEP
!     IF INCOMPLETE, SET FOR NEXT PASS
 
 2090 CALL CLOSE (scri,clsrew)
 CALL CLOSE (scro,clsrew)
 IF (ielem == 0) GO TO 2100
 k = scri
 scri = scro
 scro = k
 l = n21 + 1
 DO  j = 1,n21
   z(j) = 0
 END DO
 z(i) = l
 z(i+1) = ll
 GO TO 2073
 
!     READ THE SIL INTO CORE. OPEN ECPT0 AND GPECT.
!     WRITE HEADER RECORD ON GPECT - 3RD WORD = NO OF ENTRIES IN /GPTA1/
 
 2100 FILE = sil
 CALL gopen (sil,z(buf1),rdrew)
 CALL fread (sil,z,nsil,1)
 z(nsil+1) = luset + 1
 CALL CLOSE (sil,clsrew)
 infile = scro
 oufile = gpect
 maxdof = 0
 FILE   = infile
 CALL OPEN (*3200,infile,z(buf1),rdrew)
 CALL OPEN (*3200,oufile,z(buf2),wrtrew)
 CALL fname (oufile,buf)
 buf(3) = nelem
 CALL WRITE (oufile,buf,3,1)
 
!     PASS ECPT0 LINE BY LINE. FOR EACH LINE -
!     1. CONVERT GRID NBRS TO SIL VALUES
!     2. SORT SIL NBRS AND DISCARD MISSING ONES
!     3. WRITE LINE ON GPECT
 
 DO  ll = 1,nsil
   
!     WRITE SIL AND DOF FOR PIVOT
   
   buf(1) = z(ll)
   buf(2) = z(ll+1) - z(ll)
   CALL WRITE (oufile,buf,2,0)
   
!     READ AN ECT LINE FROM ECPT0. SET POINTERS AS A FUNCTION OF ELEM
!     TYPE.
   
   2140 CALL READ (*3201,*2154,infile,buf,1,0,flag)
   IF (buf(1) < 0.0) THEN
     GO TO  3207
   ELSE IF (buf(1) == 0.0) THEN
     GO TO  2150
   END IF
   2142 CALL READ (*3201,*3202,infile,buf(2),buf(1),0,flag)
   khr   = buf(2)/two24
   ielem = (khr-1)*incr + 1
   ngrids= elem(ielem+9)
   igr1  = elem(ielem+12) + 1
   igr2  = igr1 + ngrids - 1
   maxel = 0
   
!     CONVERT GRID NUMBERS TO SIL VALUES. DISCARD ANY MISSING (ZERO)
!     GRID POINTS THEN SORT LIST ON SIL VALUE
   
   DO  ii = igr1,igr2
     k = buf(ii)
     IF (k /= 0) GO TO 2145
     buf(ii) = 2147483647
     ngrids  = ngrids - 1
     CYCLE
     2145 buf(ii) = z(k)
     maxel   = MAX0(maxel,z(k+1)-z(k))
   END DO
   CALL sort (0,0,1,1,buf(igr1),elem(ielem+9))
   maxdof = MAX0(maxdof,ngrids*maxel)
   
!     WRITE A LINE ON GPECT.
!     - NUMBER OF WORDS IN ENTRY (NOT INCLUDING THIS WORD)
!       ELEMENT ID
!       ELEMENT TYPE
!       THE SORTED SIL LIST FOR THE GRID POINTS
   
   out(1) = -(ngrids+2)
   out(2) = buf(2) - khr*two24
   out(3) = elem(ielem+2)
   CALL WRITE (oufile,out,3,0)
   CALL WRITE (oufile,buf(igr1),ngrids,0)
   GO TO 2140
   
!     HERE IF NO ELEMENTS CONNECTED TO PIVOT.
   
   2150 CALL WRITE (oufile,0,0,1)
   CALL fwdrec (*3202,infile)
   CYCLE
   
!     HERE WHEN ALL ELEMENTS COMPLETE FOR CURRENT PIVOT
   
   2154 CALL WRITE (oufile,0,0,1)
 END DO
 
!     CLOSE FILES, WRITE TRAILER AND RETURN.
 
 CALL CLOSE (infile,clsrew)
 CALL CLOSE (oufile,clsrew)
 buf(1) = oufile
 buf(2) = nelem
 buf(3) = nsil
 buf(4) = maxel
 buf(5) = maxdof
 buf(6) = 0
 buf(7) = 0
 CALL wrttrl (buf)
 RETURN
 
!     FATAL ERROR MESAGES
 
 3200 j = -1
 GO TO 3220
 3201 j = -2
 GO TO 3220
 3202 j = -3
 GO TO 3220
 3207 buf(1) = 0
 buf(2) = 0
 CALL mesage (-30,14,buf)
 3209 buf(1) = 0
 buf(2) = 0
 CALL mesage (-30,13,buf)
 buf(1) = tempid
 buf(2) = 0
 n = 44
 GO TO 3219
 3219 CALL mesage (-30,n,buf)
 3220 CALL mesage (j,FILE,nam)
 RETURN
END SUBROUTINE ta1h
