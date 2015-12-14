SUBROUTINE ta1b
     
!     TA1B BUILDS THE ELEMENT CONNECTION AND PROPERTIES TABLE (ECPT)
!     AND THE GRID POINT CONNECTION TABLE. THE ECPT CONTAINS ONE LOGICAL
!     RECORD FOR EACH GRID OR SCALAR POINT IN THE STRUCTURE.  EACH
!     LOGICAL RECORD CONTAINS EST TYPE DATA FOR ELEMENTS CONNECTED TO
!     THE GRID OR SCALAR POINT THE GPCT IS A SUMMARY OF THE ECPT.  EACH
!     LOGICAL RECORD CONTAINS ALL GRID POINTS CONNECTED TO THE PIVOT (BY
!     MEANS OF STRUCTURAL ELEMENTS).
 
 EXTERNAL        andf
 LOGICAL :: eorflg,endid ,record
 INTEGER :: andf  ,genl  ,ect   ,ept   ,bgpdt ,sil   ,gptt  ,  &
     cstm  ,est   ,gei   ,ecpt  ,gpct  ,scr1  ,scr2  ,  &
     scr3  ,scr4  ,z     ,sysbuf,tempid,elem  ,elemid,  &
     outpt ,cbar  ,plotel,rd    ,rdrew ,wrt   ,wrtrew,  &
     clsrew,cls   ,buf   ,gpsav ,flag  ,buf1  ,buf2  ,  &
     buf3  ,FILE  ,ret   ,ret1  ,op    ,two24 ,scri  ,  &
     scro  ,blk   ,ret2  ,oufile,gpect ,eltype,oldel ,  &
     oldeid,buf4  ,eqexin,zeros(4)     ,quadts,triats,  &
     plot  ,react ,shear ,twist ,bar   ,ppse
 DIMENSION       nam(2),blk(2),zz(1) ,tgrid(33)    ,buf(50),  &
     bufr(50)     ,gpsav(34)
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm   ,uwm   ,uim   ,sfm
 COMMON /BLANK / luset ,nosimp,nosup ,nogenl,genl  ,comps
 COMMON /ta1com/ nsil  ,ect   ,ept   ,bgpdt ,sil   ,gptt  ,cstm  ,  &
     mpt   ,est   ,gei   ,gpect ,ecpt  ,gpct  ,mptx  ,  &
     pcomps,eptx  ,scr1  ,scr2  ,scr3  ,scr4  ,eqexin
 COMMON /system/ ksystm(65)
 COMMON /gpta1 / nelem ,jlast ,incr  ,elem(1)
 COMMON /names / rd    ,rdrew ,wrt   ,wrtrew,clsrew,cls
 COMMON /ta1ett/ eltype,oldel ,eorflg,endid ,bufflg,itemp ,idftmp,  &
     iback ,record,oldeid
 COMMON /zzzzzz/ z(1)
 EQUIVALENCE     (ksystm( 1),sysbuf) ,(ksystm(2),outpt ),  &
     (ksystm(10),tempid) ,(idftmp   ,deftmp)
 EQUIVALENCE     (blk(1),npvt), (buf(1),bufr(1)),(z(1),zz(1)), (blk(2),n)
 DATA    nam   / 4HTA1B,3H   /,  cbar/ 4HBAR  /, plot/ 4HPLOT    /
 DATA    two24 / 8388608     /, zeros/ 0,0,0,0/, ppse/ 4303      /
 DATA    plotel, react,shear,twist,ihex2,ihex3,quadts,triats,bar /  &
     5201  , 5251, 4,    5,    66,   67,   68,    69,    34  /
 
!     PERFORM GENERAL INITIALIZATION
 
 n2   = 2*nelem - 1
 n21  = n2 + 1
 buf1 = korsz(z) - sysbuf - 2
 buf2 = buf1 - sysbuf
 buf3 = buf2 - sysbuf
 buf4 = buf3 - sysbuf
 neq1 = nsil + 1
 neq2 = 0
 kscalr = 0
 
!     THE GRID POINT COUNTER (GPC) HAS ONE ENTRY PER GRID OR SCALAR
!     POINT IN THE STRUCTURE. EACH ENTRY CONTAINS THE NUMBER OF
!     STRUCTURAL ELEMENTS CONNECTED TO THE POINT.
 
 DO  i = 1,nsil
   z(i+1) = 0
 END DO
 
!     OPEN THE ECT. INITIALIZE TO LOOP THRU BY ELEMENT TYPE.
 
 FILE = ect
 CALL preloc (*3200,z(buf1),ect)
 noect = 1
 DO  i = 1,jlast,incr
   
!     IGNORE PLOTEL ELEMENTS.  OTHERWISE, LOCATE AN ELEMENT TYPE.
!     IF PRESENT, READ ALL ELEMENTS OF THAT TYPE AND INCREMENT THE GPC
!     ENTRY FOR EACH POINT TO WHICH THE ELEMENT IS CONNECTED.
   
   IF (elem(i) == plot) CYCLE
   CALL locate (*2026,z(buf1),elem(i+3),flag)
   noect = 0
   lx = elem(i+12)
   mm = lx + elem(i+9) - 1
   m  = elem(i+5)
   IF (elem(i+10) == 0) kscalr = 1
   2021 CALL READ (*3201,*2026,ect,buf,m,0,flag)
   DO  l = lx,mm
     k  = buf(l)
     IF (k /= 0) z(k+1) = z(k+1) + 1
   END DO
   GO TO 2021
 END DO
 CALL CLOSE (ect,clsrew)
 IF (noect /= 0) GO TO 3209
 
!     REPLACE ENTRIES IN THE GPC BY A RUNNING SUM THUS CREATING POINTERS
!     INTO ECPT0.  QUEUE WARNING MESSAGES FOR GRID PTS. WITH NO ELEMENTS
!     CONNECTED.
!     (BRING IN EQEXIN AND ECHO OUT EXTERNAL GRID PT. ID  G.C/UNISYS 91)
 
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
!WKBR CALL GOPEN (EQEXIN,EQEXIN,Z(BUF1),RDREW)
   CALL gopen (eqexin,z(buf1),rdrew)
   CALL READ (*3200,*2032,eqexin,z(neq1),buf4,1,neq2)
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
 op   = wrtrew
 2042 IF (z(ndx2+1)-z(ndx1)+2 <= length) GO TO 2050
 ndx2 = ndx2 - 1
 GO TO 2042
 
!     PASS THE ECT. FOR EACH GRID PT IN RANGE ON THIS PASS,
!     STORE ELEMENT POINTER = 2**24 * J + WORD POSITION IN ECT RECORD.
!     WHERE J= (POINTER IN ELEM TABLE - 1)/INCR * 2 +1
 
 2050 CALL preloc (*3200,z(buf1),ect)
 izero = z(ndx1)
 j = 1
 DO  i = 1,jlast,incr
   IF (elem(i) == plot) GO TO 2055
   idcntr = two24*j
   CALL locate (*2055,z(buf1),elem(i+3),flag)
   m  = elem(i+ 5)
   lx = elem(i+12)
   mm = lx + elem(i+9) - 1
   2052 CALL READ (*3201,*2055,ect,buf,m,0,flag)
   DO  l = lx,mm
     k  = buf(l)
     IF (k < ndx1 .OR. k > ndx2) CYCLE
     ix = z(k) - izero + iecpt0
     z(ix) = idcntr
     z(k ) = z(k) + 1
   END DO
   idcntr = idcntr + m
   GO TO 2052
   2055 j = j + 2
 END DO
 CALL CLOSE (ect,clsrew)
 
!     WRITE ECPT0 AND TEST FOR ADDITIONAL PASSES
!     ECPT0 CONTAINS ONE LOGICAL RECORD FOR EACH GRID OR SCALAR POINT.
!     EACH LOGICAL RECORD CONTAINS N PAIRS OF(-1,ELEMENT POINTER)WHERE
!     N= NUMBER OF ELEMENTS CONNECTED TO THE PIVOT.
!     IF NO ELEMENTS CONNECTED TO POINT, RECORD IS ONE WORD = 0.
 
 FILE = scr1
 CALL OPEN (*3200,scr1,z(buf1),op)
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
 CALL preloc (*3200,z(buf1),ect)
 i = 1
 ielem = 1
 DO  j = 1,n21
   z(j) = 0
 END DO
 l = n21 + 1
 2072 IF (elem(ielem+3) == plotel .OR. elem(ielem+3) == react) GO TO 2074
 CALL locate (*2074,z(buf1),elem(ielem+3),flag)
 z(i) = l
 ll   = 0
 m    = elem(ielem+5)
 last = buf3 - m
 2073 IF (l > last) GO TO 2080
 CALL READ (*3201,*2074,ect,z(l),m,0,flag)
 l  = l  + m
 ll = ll + m
 GO TO 2073
 2074 i  = i + 2
 ielem = ielem + incr
 IF (ielem <= jlast) GO TO 2072
 
!     PASS ECPT0 ENTRIES LINE BY LINE
!     ATTACH EACH REFERENCED ECT ENTRY WHICH IS NOW IN CORE
 
 2080 CALL OPEN (*3200,scri,z(buf2),rdrew)
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
 k = buf(2)/two24
 ktwo24 = k*two24
 idptr  = buf(2) - ktwo24
 kk = z(k) + idptr - z(k+1)
 IF (z(k) == 0 .OR. kk > last) GO TO 2084
 j  = ((k-1)/2)*incr + 1
 mm = elem(j+5)
 buf(1) = mm
 buf(2) = z(kk) + ktwo24
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
 IF (i > n2) GO TO 2100
 k    = scri
 scri = scro
 scro = k
 l = n21 + 1
 DO  j = 1,n21
   z(j) = 0
 END DO
 z(i) = l
 z(i+1) = ll
 GO TO 2073
 
!     READ THE EPT INTO CORE (IF PRESENT)
!     FIRST N21 CELLS OF CORE CONTAINS PROPERTIES POINTER TABLE WHICH
!     HAS TWO WORDS PER ELEMENT TYPE, 1ST WORD HAS POINTER TO 1ST WORD
!     OF PROPERTY DATA FOR THAT ELEMENT TYPE. 2ND WORD HAS NUMBER OF
!     PROPERTY CARDS FOR THAT TYPE.
 
 2100 CALL CLOSE (ect,clsrew)
 DO  i = 1,n21
   z(i) = 0
 END DO
 l = 1
 CALL preloc (*2120,z(buf1),ept)
 ielem  = 1
 lstprp = 0
 l = n21 + 1
 DO  ii = 1,n2,2
   IF (elem(ielem+6) == lstprp .AND. lstprp /= ppse) GO TO 2106
   CALL locate (*2107,z(buf1),elem(ielem+6),flag)
   lstprp = elem(ielem+6)
   m      = elem(ielem+8)
   eltype = elem(ielem+2)
   z(ii)  = l
   2102 IF (l+m >= buf3) CALL mesage (-8,0,nam)
   CALL READ (*3201,*2103,ept,z(l),m,0,flag)
   l = l + m
   GO TO 2102
   2103 n = l - z(ii)
   z(ii+1) = n/m
   IF (eltype == shear .OR. eltype == twist) GO TO 2104
   IF (m > 4) GO TO 2107
   2104 i = z(ii)
   CALL sort (0,0,m,1,z(i),n)
   GO TO 2107
   2106 n = 4
   IF (eltype == ihex2 .OR. eltype == ihex3) n = 2
   z(ii  ) = z(ii-n  )
   z(ii+1) = z(ii-n+1)
   2107 ielem   = ielem + incr
 END DO
 CALL CLOSE (ept,clsrew)
 
!     DETERMINE IF THE BGPDT AND SIL
!     WILL FIT IN CORE ON TOP OF THE EPT.
 
 NUMBER = 4*kscalr + 1
 iback  = 0
 length = buf4 - l - 4*maxel
 IF (NUMBER*nsil > length) GO TO 2150
 
!     IF YES, READ THE BGPDT,SIL AND GPTT INTO CORE
 
 2120 ASSIGN 2130 TO ret
 ipass = 1
 GO TO 3050
 
!     PASS ECPT0 LINE BY LINE
!     FOR EACH ECT ENTRY, 1. ATTACH PROPERTY DATA (IF DEFINED)
!     2. ATTACH BASIC GRID POINT DATA (UNLESS SCALER ELEMENT), AND
!     3. CONVERT GRID PT NOS TO SIL VALUES
!     4. IF TEMPERATURE PROBLEM, ATTACH ELEMENT TEMP(UNLESS SCALAR ELEM)
 
 2130 infile = scro
 oufile = ecpt
 
!     OPEN ECPT0, ECPT AND GPCT FILES
 
 2144 GO TO 3060
 
!     WRITE PIVOT GRID POINT ON ECPT
 
 2131 IF (ll-locsil >= nsil) GO TO 2179
 IF (iback <= 0) GO TO 21311
 CALL bckrec (gptt)
 
!     RESET /TA1ETT/ VARIABLES
 
 iback  = 0
 oldeid = 0
 oldel  = 0
 eorflg =.false.
 endid  =.true.
 CALL READ (*3201,*3202,gptt,iset,1,0,flag)
 IF (iset == tempid) GO TO 21311
 WRITE (outpt,3084) sfm,iset,tempid
 CALL mesage (-61,0,0)
 21311 npvt = z(ll)
 CALL WRITE (ecpt,npvt,1,0)
 IF (z(ll+1)-z(ll) == 1) npvt = -npvt
 i = locgpc
 
!     READ AN ECT LINE FROM ECPT0. SET POINTERS AS A FUNCTION OF ELEM
!     TYPE.  IF ELEMENT IS BAR, PROCESS ORIENTATION VECTOR.  AXIS AND
!     THE STRESS AXIS DEFINITION BASED ON GRID POINTS MA AND SA.
 
 2132 CALL READ (*3201,*2138,infile,buf(1),1,0,flag)
 IF (buf(1) < 0.0) THEN
   GO TO  3207
 ELSE IF (buf(1) == 0.0) THEN
   GO TO  2143
 END IF
 2133 CALL READ (*3201,*3202,infile,buf(2),buf(1),0,flag)
 ik = buf(2)/two24
 ii = ((ik-1)/2)*incr + 1
 lx = elem(ii+12) + 1
 m  = elem(ii+ 8)
 jscalr =  elem(ii+10)
 mm = lx + elem(ii+ 9) - 1
 lq = 4
 IF (m == 0) lq = 3
 NAME  = elem(ii   )
 jtemp = elem(ii+13)
 eltype= elem(ii+ 2)
 ntemp = 1
 IF (jtemp == 4) ntemp = elem(ii+14) - 1
 IF (eltype == quadts) GO TO 3083
 IF (eltype == triats) GO TO 30841
 IF (NAME   ==   cbar) GO TO 3080
 
!     SAVE INTERNAL GRID NOS AND CONVERT TO SIL NOS.
 
 2141 GO TO 3030
 
!     IF ONE   PASS, WRITE ECT       SECTION  OF ECPT LINE.
!     IF TWO PASSES, WRITE ECT + EPT SECTIONS OF ECPT LINE.
 
 2134 id = buf(3)
 nx = buf(1) + 2 - lq
 buf(1) = elem(ii+2)
 buf(2) = buf(2) - ik*two24
 elemid = buf(2)
 CALL WRITE (ecpt,buf(1),2,0)
 CALL WRITE (ecpt,buf(lq),nx,0)
 IF (ipass == 2) GO TO 2137
 
!     IF PROPERTY DATA IS DEFINED, LOOK UP AND WRITE EPT SECTION OF ECPT
 
 IF (m == 0) GO TO 2137
 ASSIGN 2137 TO ret
 GO TO 3040
 
!     IF ELEMENT IS NOT A SCALAR ELEMENT,
!     WRITE BGPDT AND ELEMENT TEMPERATURE SECTIONS OF ECPT LINE.
 
 2137 IF (jscalr /= 0) GO TO 2132
 GO TO 3090
 
!     CLOSE ECPT RECORD. WRITE GPCT RECORD.
 
 2138 CALL WRITE (ecpt,0,0,1)
 GO TO 3070
 
!     HERE IF NO ELEMENTS CONNECTED TO PIVOT.
 
 2143 CALL WRITE (ecpt,0,0,1)
 IF (nogpct /= 0) CALL WRITE (gpct,npvt,1,1)
 ll = ll + 1
 CALL fwdrec (*3202,infile)
 GO TO 2131
 
!     HERE IF ECPT CONSTRUCTION IS TWO PASSES.
!     PASS ECPT0 LINE BY LINE FOR EACH ECT ENTRY, ATTACH PROPERTY DATA
!     IF DEFINED
 
 2150 CALL OPEN (*3200,scro,z(buf1),rdrew)
 CALL OPEN (*3200,scri,z(buf2),wrtrew)
 oufile = scri
 
!     READ AN ECT LINE FROM ECT0. SET POINTERS AS FUNCTION OF ELEM TYPE.
 
 2152 CALL READ (*2159,*2156,scro,buf,1,0,flag)
 IF (buf(1) < 0.0) THEN
   GO TO  3207
 ELSE IF (buf(1) == 0.0) THEN
   GO TO  2158
 END IF
 2153 CALL READ (*3201,*3203,scro,buf(2),buf(1),0,flag)
 ik = buf(2)/two24
 ii = ((ik-1)/2)*incr + 1
 m  = elem(ii+8)
 nx = buf(1) + 1
 
!     IF PROPERTY DATA IS DEFINED FOR ELEMENT, WRITE ECT DATA ON SCRI,
!     THEN LOOK UP AND WRITE EPT DATA ON SCRI.
 
 IF (m == 0) GO TO 2155
 id = buf(3)
 buf(1) = buf(1) + m - 1
 CALL WRITE (scri,buf(1),nx,0)
 ASSIGN 2152 TO ret
 GO TO 3040
 
!     PROPERTY DATA NOT DEFINED. WRITE ECT LINE ON SCRI.
 
 2155 CALL WRITE (scri,buf,nx,0)
 GO TO 2152
 
!     CLOSE RECORD. RETURN FOR ANOTHER PIVOT.
 
 2156 CALL WRITE (scri,0,0,1)
 GO TO 2152
 
!     ALL PIVOTS COMPLETE. CLOSE FILES.
 
 2159 CALL CLOSE (scro,clsrew)
 CALL CLOSE (scri,clsrew)
 GO TO 2160
 
!     HERE IF NO ELEMENTS CONNECTED TO PIVOT.
 
 2158 CALL WRITE  (scri,0,1,1)
 CALL fwdrec (*3201,scro)
 GO TO 2152
 
!     READ THE BGPDT, SIL AND, IF TEMPERATURE PROBLEM,
!     THE GPTT INTO CORE.
 
 2160 l = 1
 ASSIGN 2170 TO ret
 GO TO 3050
 
!     SET POINTERS AND BRANCH TO COMMON CODE TO ASSEMBLE ECPT.
 
 2170 infile = scri
 oufile = ecpt
 ipass  = 2
 GO TO 2144
 
!     CLOSE FILES, WRITE TRAILERS AND EXIT.
 
 2179 CALL CLOSE (infile,clsrew)
 CALL CLOSE (gptt,clsrew)
 CALL CLOSE (ecpt,clsrew)
 buf(1) = ect
 CALL rdtrl (buf(1))
 buf(3) = 0
 k  = 8192
 k1 = andf(buf(5),k)
 IF (k1 /= k) GO TO 2180
 buf(3) = 1
 irigd  = 1
 2180 CONTINUE
 buf(1) = ecpt
 DO  i = 2,7
   buf(i) = 7
 END DO
 CALL wrttrl (buf)
 IF (nogpct == 0) RETURN
 CALL CLOSE (gpct,clsrew)
 buf(1) = gpct
 CALL wrttrl (buf)
 RETURN
 
 
!     INTERNAL BINARY SEARCH ROUTINE
 
 3000 klo = 1
 3001 k   = (klo+khi+1)/2
 3008 kx  = (k-1)*m + locx
 IF (id-z(kx) < 0) THEN
   GO TO  3002
 ELSE IF (id-z(kx) == 0) THEN
   GO TO  3009
 ELSE
   GO TO  3003
 END IF
 3002 khi = k
 GO TO 3004
 3003 klo = k
 3004 IF (khi-klo-1 < 0) THEN
   GO TO 30091
 ELSE IF (khi-klo-1 == 0) THEN
   GO TO  3005
 ELSE
   GO TO  3001
 END IF
 3005 IF (k == klo) GO TO 3006
 k   = klo
 GO TO 3007
 3006 k   = khi
 3007 klo = khi
 GO TO 3008
 3009 GO TO ret1, (3041)
 30091 GO TO ret2, (3205)
 
 
!     INTERNAL ROUTINE TO SAVE GRID PTS IN AN ECT LINE
!     AND CONVERT GRID PT NOS IN ECT LINE TO SIL VALUES
 
 3030 DO  l = lx,mm
   gpsav(l) = 0
   IF (buf(l) == 0) CYCLE
   gpsav(l) = buf(l)
   k  = gpsav(l) + locsil - 1
   buf(l) = z(k)
   ix = 0
   IF (z(k+1)-z(k) == 1) ix = 1
   z(i) = 2*z(k) + ix
   i  = i + 1
 END DO
 IF (i >= buf3) CALL mesage (-8,0,nam)
 GO TO 2134
 
 
!     INTERNAL ROUTINE TO ATTACH EPT DATA
 
 3040 locx = z(ik)
 IF (locx == 0) GO TO 3206
 khi  = z(ik+1)
 ASSIGN 3041 TO ret1
 ASSIGN 3205 TO ret2
 GO TO 3000
 3041 CALL WRITE (oufile,z(kx+1),m-1,0)
 GO TO ret, (2137,2152)
 
!     INTERNAL ROUTINE TO READ THE BGPDT, SIL AND GPTT INTO CORE
 
 3050 nbgp   = 0
 locbgp = l
 IF (kscalr == 0) GO TO 3059
 CALL OPEN (*3200,bgpdt,z(buf1),rdrew)
 CALL fwdrec (*3201,bgpdt)
 nbgp  = 4*nsil
 CALL READ (*3201,*3202,bgpdt,z(locbgp),nbgp,1,flag)
 CALL CLOSE (bgpdt,clsrew)
 3059 l = l + nbgp
 CALL OPEN (*3200,sil,z(buf1),rdrew)
 CALL fwdrec (*3201,sil)
 locsil = locbgp + nbgp
 CALL READ (*3201,*3203,sil,z(locsil),nsil,1,flag)
 CALL CLOSE (sil,clsrew)
 nx = locsil + nsil
 z(nx)  = luset + 1
 loctmp = nx + 1
 ntmp   = loctmp - 1
 record =.false.
 itemp  = tempid
 iback  = 0
 IF (tempid == 0) GO TO 3058
 FILE   = gptt
 CALL OPEN (*3200,gptt,z(buf4),rdrew)
 CALL READ (*3201,*3051,gptt,z(loctmp),buf3-loctmp,0,nid)
 CALL mesage (-8,0,nam)
 3051 itmpid = loctmp + 2
 ntmpid = loctmp + nid - 3
 DO  ijk = itmpid,ntmpid,3
   IF (tempid == z(ijk)) GO TO 3053
 END DO
 GO TO 3210
 3053 idftmp = z(ijk+1)
 IF (idftmp . NE. -1) deftmp = zz(ijk+1)
 n = z(ijk+2)
 IF (n == 0) GO TO 3058
 record =.true.
 n = n - 1
 IF (n == 0) GO TO 3055
 DO  ijk = 1,n
   CALL fwdrec (*3201,gptt)
 END DO
 
!     READ SET ID AND VERIFY FOR CORRECTNESS
 
 3055 CALL READ (*3201,*3202,gptt,iset,1,0,flag)
 IF (iset == tempid) GO TO 3061
 WRITE  (outpt,3084) sfm,iset,tempid
 3084 FORMAT (a25,' 4021, TA1B HAS PICKED UP TEMPERATURE SET',i9,  &
     ' AND NOT THE REQUESTED SET',i9,1H.)
 CALL mesage (-61,0,nam)
 
!     INITIALIZE /TA1ETT/ VARIABLES
 
 3061 oldeid = 0
 oldel  = 0
 eorflg =.false.
 endid  =.true.
 3058 GO TO ret, (2130,2170)
 
 
!     INTERNAL ROUTINE TO OPEN SCRATCH, ECPT AND GPCT FILES
 
 3060 CALL OPEN (*3200,infile,z(buf1),rdrew)
 CALL OPEN (*3200,ecpt,z(buf2),wrtrew)
 CALL fname (ecpt,buf)
 CALL WRITE (ecpt,buf,2,1)
 nogpct = 0
 CALL OPEN (*3062,gpct,z(buf3),wrtrew)
 nogpct = 1
 CALL fname (gpct,buf)
 CALL WRITE (gpct,buf,2,1)
 3062 ll = locsil
 locgpc = ntmp + 1
 GO TO 2131
 
!     INTERNAL ROUTINE TO SORT AND WRITE THE GPCT
 
 3070 IF (nogpct == 0) GO TO 3073
 n = i - locgpc
 CALL sort (0,0,1,1,z(locgpc),n)
 z(i) = 0
 j  = locgpc
 ii = locgpc
 3071 IF (z(ii) == z(ii+1)) GO TO 3072
 nx = z(ii)/2
 lx = z(ii) - 2*nx
 IF (lx /= 0) nx = -nx
 z(j) = nx
 j  = j  + 1
 3072 ii = ii + 1
 IF (ii < i) GO TO 3071
 n  = j - locgpc
 CALL WRITE (gpct,blk,2,0)
 CALL WRITE (gpct,z(locgpc),n,1)
 3073 ll = ll + 1
 GO TO 2131
 
!     FOR BAR ELEMENTS, STORE COORDINATES AND
!     COORDINATE SYSTEM ID FOR ORIENTATION VECTOR.
 
 3080 kx = 4*(buf(4)-1) + locbgp
 IF (buf(9) == 1) GO TO 3082
 buf(9) = buf(6)
 IF (buf(9) == 0) GO TO 3082
 k = 4*(buf(9)-1) + locbgp
 bufr(6) = zz(k+1) - zz(kx+1)
 bufr(7) = zz(k+2) - zz(kx+2)
 bufr(8) = zz(k+3) - zz(kx+3)
 buf(9)  = 0
 GO TO 2141
 3082 buf(9)  = z(kx)
 GO TO 2141
 
!     FOR QUADTS AND TRIATS ELEMENTS, STORE COORDINATES FOR MATERIAL
!     AND STRESS AXIS DEFINITION
 
 3083 is1 = 12
 GO TO 3085
 30841 is1 = 10
 3085 is2 = is1 + 9
 DO  ist = is1,is2,3
   igp = buf(ist)
   IF (igp == 0) CYCLE
   k = 4*(igp-1) + locbgp
   bufr(ist  ) = zz(k+1)
   bufr(ist+1) = zz(k+2)
   bufr(ist+2) = zz(k+3)
 END DO
 GO TO 2141
 
!     CODE TO WRITE BGPDT AND ELEMENT TEMPERATURE SECTIONS OF ECAT LINE.
 
 3090 DO  l = lx,mm
   IF (gpsav(l) == 0) GO TO 3094
   k = locbgp + 4*(gpsav(l)-1)
   CALL WRITE (ecpt,z(k),4,0)
   CYCLE
   3094 CALL WRITE (ecpt,zeros,4,0)
 END DO
 CALL ta1etd (elemid,tgrid,ntemp)
 IF (eltype == bar) tgrid(1) = (tgrid(1)+tgrid(2))/2.0
 CALL WRITE (ecpt,tgrid,ntemp,0)
 GO TO 2132
 
!     FATAL ERROR MESAGES
 
 3200 j = -1
 GO TO 3220
 3201 j = -2
 GO TO 3220
 3203 CONTINUE
 3202 j = -3
 GO TO 3220
 3205 buf(1) = elemid
 buf(2) = id
 n = 10
 GO TO 3219
 3206 buf(1) = elem(ii  )
 buf(2) = elem(ii+1)
 n = 11
 GO TO 3219
 3207 buf(1) = 0
 buf(2) = 0
 n = 14
 GO TO 3219
 3209 buf(1) = 0
 buf(2) = 0
 n = 13
 GO TO 3219
 3210 buf(1) = tempid
 buf(2) = 0
 n = 44
 3219 CALL mesage (-30,n,buf)
 3220 CALL mesage (j,FILE,nam)
 RETURN
END SUBROUTINE ta1b
