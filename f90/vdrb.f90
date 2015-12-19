SUBROUTINE vdrb (infil,outfl,ireqq)
     
!     VDRB PROCESSES VECTORS IN THE ANALYSIS OR MODAL SET. IN
!     ACCORDANCE WITH OUTPUT REQUESTS IN THE CASE CONTROL DATA BLOCK,
!     THESE VECTORS ARE FORMATTED FOR INPUT TO OFP WHERE ACTUAL OUTPUT
!     WILL OCCUR.
 
 
 INTEGER, INTENT(IN)                      :: infil
 INTEGER, INTENT(IN)                      :: outfl
 INTEGER, INTENT(IN)                      :: ireqq
 EXTERNAL        andf
 INTEGER :: app   ,FORM  ,sort2 ,output,z     ,sysbuf,date  ,  &
     time  ,ud    ,ue    ,two   ,qtype2,cei   ,frq   ,  &
     trn   , modal ,DIRECT,casecc,eqdyn ,usetd ,  &
      oeigs ,pp    ,buf   ,buf1  ,buf2  ,buf3  ,  &
     FILE  ,flag  ,sild  ,code  ,gptype,andf  ,branch,  &
     setno ,fsetno,word  ,ret   ,retx  ,FORMAT,eof   ,  &
     vdrcom,sdr2  ,xset0 ,xsetno,dest  ,axif  ,vdrreq, oharms
 DIMENSION       mcb(7)    ,buf(50)      ,bufr(50)     ,masks(6) ,  &
     zz(1)     ,cei(2)       ,frq(2)       ,trn(2)   ,  &
     modal(2)  ,DIRECT(2)    ,nam(2)       ,vdrcom(1)
 COMMON /condas/ consts(5)
 COMMON /BLANK / app(2),FORM(2),sort2,output,sdr2  ,imode
 COMMON /vdrcom/ vdrcom,idisp ,ivel  ,iacc  ,ispcf ,iloads,istr  ,  &
     ielf  ,iadisp,iavel ,iaacc ,ipnl  ,ittl  ,ilsym ,  &
     ifrout,idload,casecc,eqdyn ,usetd ,infile,oeigs ,  &
     pp    ,xycdb ,pnl   ,outfle,opnl1 ,scr1  ,scr2  ,  &
     buf1  ,buf2  ,buf3  ,nam   ,buf   ,masks ,cei   ,  &
     frq   ,trn   ,DIRECT,xset0 ,vdrreq,modal
 COMMON /zzzzzz/ z(1)
 COMMON /system/ sysbuf,xx(13),date(3),time,dum19(19),axif
 COMMON /names / rd    ,rdrew ,wrt   ,wrtrew,clsrew
 COMMON /bitpos/ um    ,uo    ,ur    ,usg   ,usb   ,ul    ,ua    ,  &
     uf    ,us    ,un    ,ug    ,ue    ,up    ,une   , ufe   ,ud
 COMMON /two   / two(32)
 COMMON /unpakx/ qtype2,i2    ,j2    ,incr2
 EQUIVALENCE     (consts(2),twopi)   ,(consts(3),raddeg)  ,  &
     (buf(1),bufr(1))    ,(z(1),zz(1))
 DATA    igpf  , iese,ireig   / 167   , 170, 4HREIG  /
 
!     PERFORM GENERAL INITIALIZATION.
 
 m8    = -8
 mskud = two(ud)
 mskue = two(ue)
 ilist = 1
 i2    = 1
 incr2 = 1
 ireq  = ireqq
 IF (FORM(1) /= modal(1) .AND. FORM(1) /= DIRECT(1)) GO TO 1432
 
!     READ TRAILER ON USETD. SET NO. OF EXTRA POINTS.
!     READ TRAILER ON INFIL. SET PARAMETERS.
!     IF MODAL PROBLEM, NO. OF MODES = NO. OF ROWS IN VECTOR - NO. XTRA
!     PTS.
 
 mcb(1) = usetd
 FILE   = usetd
 CALL rdtrl (mcb)
 IF (mcb(1) /= usetd) GO TO 2001
 nbrep  = mcb(3)
 mcb(1) = infil
 CALL rdtrl (mcb)
 IF (mcb(1) /= infil) GO TO 1431
 nvects = mcb(2)
 nrows  = mcb(3)
 IF (FORM(1) == modal(1)) nbrmod = imode + nrows - nbrep - 1
 IF (mcb(5) > 2) GO TO 1022
 IF (app(1) == frq(1)) GO TO 1022
 
!     REAL VECTOR.
 
 ktype  = 1
 qtype2 = 1
 nwds   = 8
 ktypex = 0
 GO TO 1030
 
!     COMPLEX VECTOR.
 
 1022 ktype  = 2
 qtype2 = 3
 nwds   = 14
 ktypex = 1000
 
!     IF DIRECT PROBLEM OR MODAL PROBLEM WITH EXTRA POINTS,
!     READ 2ND TABLE OF EQDYN INTO CORE. THEN READ USETD INTO CORE.
 
 1030 IF (FORM(1) == modal(1) .AND. nbrep == 0) GO TO 1050
 FILE = eqdyn
 CALL gopen (eqdyn,z(buf1),0)
 CALL fwdrec (*2002,eqdyn)
 CALL READ (*2002,*1031,eqdyn,z,buf1,1,neqd)
 CALL mesage (m8,0,nam)
 1031 CALL CLOSE (eqdyn,clsrew)
 iusetd = neqd + 1
 ncore  = buf1 - iusetd
 FILE   = usetd
 CALL gopen (usetd,z(buf1),0)
 CALL READ (*2002,*1032,usetd,z(iusetd),ncore,1,flag)
 CALL mesage (m8,0,nam)
 1032 CALL CLOSE (usetd,clsrrw)
 ilist  = iusetd
 neqdyn = neqd - 1
 kn     = neqd/2
 
!     BRANCH ON PROBLEM TYPE.
 
 IF (FORM(1) == modal(1)) GO TO 1049
 
!     DIRECT - PROCESS EACH ENTRY IN EQDYN. IF POINT IS NOT IN ANALYSIS
!              SET, REPLACE SILD NO. WITH ZERO. OTHERWISE, REPLACE SILD
!              NO. WITH POSITION IN ANALYSIS SET (I.E. ROW INDEX IN
!              VECTOR) AND CODE INDICATING WHICH COMPONENTS OF POINT ARE
!              IN ANALYSIS SET.
 
 DO  i = 1,neqdyn,2
   sild   = z(i+1)/10
   gptype = z(i+1) - 10*sild
   nusetd = iusetd + sild - 1
   k = 0
   m = 1
   IF (gptype == 1) m = 6
   j = nusetd
   DO  l = 1,m
     IF (andf(z(j),mskud) /= 0) k = k + masks(l)
     j = j + 1
   END DO
   IF (k == 0) GO TO 1043
   l = 1
   m = nusetd - 1
   IF (m < iusetd) GO TO 1045
   DO  j = iusetd,m
     IF (andf(z(j),mskud) /= 0) l = l + 1
   END DO
   1045 z(i+1) = gptype + k + 256*l
   CYCLE
   1043 z(i+1) = 0
 END DO
 GO TO 1050
 
!     MODAL - PROCESS EACH ENTRY IN EQDYN. IF POINT IS NOT AN EXTRA
!             POINT, REPLACE SILD NO. WITH ZERO. OTHERWISE, REPLACE SILD
!             NO. WITH POSITION IN MODAL SET (I.E. ROW INDEX IN VECTOR).
 
 1049 DO  i = 1,neqdyn,2
   sild   = z(i+1)/10
   gptype = z(i+1) - 10*sild
   IF (gptype /= 3) GO TO 1047
   nusetd = iusetd + sild - 1
   IF (andf(z(nusetd),mskue) == 0) GO TO 1047
   k = nbrmod - imode + 1
   DO  j = iusetd,nusetd
     IF (andf(z(j),mskue) /= 0) k = k + 1
   END DO
   z(i+1) = 10*k + 3
   CYCLE
   1047 z(i+1) = 0
 END DO
 
!     SET PARAMETER FOR APPROACH. THEN OPEN CASE CONTROL,
!     SKIP HEADER RECORD AND BRANCH ON APPROACH.
 
 1050 branch = 0
 IF (app(1) == cei(1)) branch = 1
 IF (app(1) == frq(1)) branch = 2
 IF (app(1) == trn(1)) branch = 3
 IF (app(1) == ireig ) branch = 4
 IF (branch == 0) GO TO 1432
 CALL gopen (casecc,z(buf1),0)
 SELECT CASE ( branch )
   CASE (    1)
     GO TO 1060
   CASE (    2)
     GO TO 1070
   CASE (    3)
     GO TO 1070
   CASE (    4)
     GO TO 1060
 END SELECT
 
!     COMPLEX EIGENVALUES - READ LIST OF MODE NOS. AND VALUES INTO CORE.
 
 1060 FILE = oeigs
 CALL gopen (oeigs,z(buf2),0)
 CALL fwdrec (*2002,oeigs)
 i = ilist
 m = 8 - ktype
 1061 CALL READ (*2002,*1062,oeigs,buf,m,0,flag)
 z(i  ) = buf(1)
 z(i+1) = buf(3)
 z(i+2) = buf(4)
 i = i + 3
 GO TO 1061
 1062 CALL CLOSE (oeigs,clsrew)
 nlist = i - 3
 icc   = i
 GO TO 1100
 
!     FREQUENCY OR TRANSIENT RESPONSE - READ LIST INTO CORE.
 
 1070 FILE = pp
 CALL OPEN (*2001,pp,z(buf2),rdrew)
 i  = ilist
 m  = 3
 ix = 1
 IF (app(1) == frq(1)) ix = 2
 1071 CALL READ (*2002,*1072,pp,buf,m,0,flag)
 z(i  ) = buf(m)
 z(i+1) = 0
 i = i + ix
 m = 1
 GO TO 1071
 1072 CALL CLOSE (pp,clsrew)
 nlist = i - ix
 icc   = i
 
!     OPEN OUTPUT FILE. WROTE HEADER RECORD.
 
 1100 FILE = outfl
 CALL OPEN (*1431,outfl,z(buf2),wrtrew)
 mcb(1) = outfl
 CALL fname (outfl,buf)
 DO  i = 1,3
   buf(i+2) = date(i)
 END DO
 buf(6) = time
 buf(7) = 1
 CALL WRITE (outfl,buf,7,1)
 
!     OPEN INPUT FILE. SKIP HEADER RECORD.
 
 FILE = infil
 CALL OPEN (*1430,infil,z(buf3),rdrew)
 CALL fwdrec (*2002,infil)
 
!     SET PARAMETERS TO KEEP CASE CONTROL AND VECTORS IN SYNCH.
 
 eof    = 0
 jcount = 0
 kcount = 1
 jlist  = ilist
 kfrq   = 0
 kwds   = 0
 incore = 0
 
!     READ A RECORD IN CASE CONTROL.
 
 1130 CALL READ (*1400,*1131,casecc,z(icc+1),buf3-icc,1,ncc)
 CALL mesage (m8,0,nam)
 1131 ivec  = icc + ncc + 1
 ireqx = icc + idisp
 IF (z(ireqx) /= 0) sdr2 = 1
 ireqx = icc + ivel
 IF (z(ireqx) /= 0) sdr2 = 1
 ireqx = icc + iacc
 IF (z(ireqx) /= 0) sdr2 = 1
 ireqx = icc + ispcf
 IF (z(ireqx) /= 0) sdr2 = 1
 ireqx = icc + iloads
 IF (z(ireqx) /= 0) sdr2 = 1
 ireqx = icc + istr
 IF (z(ireqx) /= 0) sdr2 = 1
 ireqx = icc + ielf
 IF (z(ireqx) /= 0) sdr2 = 1
 ireqx = icc + igpf
 IF (z(ireqx) /= 0) sdr2 = 1
 ireqx = icc + iese
 IF (z(ireqx) /= 0) sdr2 = 1
 
!     SET OUTPUT HARMONICS REQUEST WHICH IS USED IF FLUID ELEMENTS
!     ARE IN PROBLEM.
 
 oharms = z(icc+137)
 IF (oharms < 0 .AND. axif /= 0) oharms = axif
 
!     IN THE ABOVE IF OHARMS = -1  THEN ALL IS IMPLIED. IF OHARMS = 0
!     THEN NONE IS IMPLIED AND IF OHARMS IS POSITIVE THEN THAT VALUE
!     MINUS ONE IS IMPLIED.
 
 IF (axif   == 0) GO TO 1140
 IF (oharms == 0) GO TO 1140
 oharms =   oharms - 1
 oharms = 2*oharms + 3
 
!     DETERMINE IF OUTPUT REQUEST IS PRESENT. IF NOT, TEST FOR RECORD
!     SKIP ON INFIL, THEN GO TO END OF REQUEST. IF SO, SET POINTERS
!     TO SET DEFINING REQUEST.
 
 1140 ireqx = icc +ireq
 setno = z(ireqx  )
 dest  = z(ireqx+1)
 xsetno = -1
 IF (setno < 0.0) THEN
   GO TO  1150
 ELSE IF (setno == 0.0) THEN
   GO TO  1141
 ELSE
   GO TO  1143
 END IF
 1141 IF (app(1) /= frq(1)) GO TO 1142
 IF (kcount /=      1) GO TO 1350
 GO TO 1150
 1142 CALL fwdrec (*2002,infil)
 jcount = jcount + 1
 GO TO 1311
 1143 ix = icc + ilsym
 isetno = ix + z(ix) + 1
 1144 iset = isetno + 2
 nset = z(isetno+1) + iset - 1
 IF (z(isetno) == setno) GO TO 1145
 isetno = nset + 1
 IF (isetno < ivec) GO TO 1144
 GO TO 1150
 
!     IF REQUIRED, LOCATE PRINT/PUNCH SUBSET.
 
 1145 IF (setno < xset0) GO TO 1150
 xsetno = dest/10
 dest   = dest - 10*xsetno
 IF (xsetno == 0) GO TO 1150
 ixsetn = ix + z(ix) + 1
 1146 ixset  = ixsetn + 2
 nxset  = z(ixsetn+1) + ixset - 1
 IF (z(ixsetn) == xsetno) GO TO 1150
 ixsetn = nxset + 1
 IF (ixsetn < ivec) GO TO 1146
 xsetno = -1
 setno  = -1
 
!     UNPACK VECTOR INTO CORE (UNLESS VECTOR IS ALREADY IN CORE).
 
 1150 IF (incore /= 0) GO TO 1160
 ivecn = ivec + ktype*nrows - 1
 IF (ivecn >= buf3) CALL mesage (m8,0,nam)
 j2 = nrows
 CALL unpack (*1151,infil,z(ivec))
 GO TO 1153
 1151 DO  i = ivec,ivecn
   zz(i)  = 0.
 END DO
 1153 jcount = jcount + 1
 
!     TEST FOR CONTINUATION.
 
 1160 IF (app(1) == frq(1) .AND. setno == 0) GO TO 1350
 
!     PREPARE TO WRITE ID RECORD ON OUTPUT FILE.
 
 SELECT CASE ( branch )
   CASE (    1)
     GO TO 1190
   CASE (    2)
     GO TO 1200
   CASE (    3)
     GO TO 1220
   CASE (    4)
     GO TO 1190
 END SELECT
 
!     COMPLEX EIGENVALUES.
 
 1190 buf(2) = 1014
 buf(5) = z(jlist  )
 buf(6) = z(jlist+1)
 buf(7) = z(jlist+2)
 buf(8) = 0
 GO TO 1250
 
!     FREQUENCY RESPONSE.
 
 1200 ix = icc + idload
 buf(8) = z(ix)
 buf(6) = 0
 buf(7) = 0
 IF (kfrq /= 0) GO TO 1207
 
!     FIRST TIME FOR THIS LOAD VECTOR ONLY - MATCH LIST OF USER
!     REQUESTED FREQS WITH ACTUAL FREQS. MARK FOR OUTPUT EACH ACTUAL
!     FREQ WHICH IS CLOSEST TO USER REQUEST.
 
 kfrq = 1
 ix   = icc + ifrout
 fsetno = z(ix)
 IF (fsetno <= 0) GO TO 1202
 ix = icc + ilsym
 isetnf = ix + z(ix) + 1
 1201 isetf  = isetnf + 2
 nsetf  = z(isetnf+1) + isetf - 1
 IF (z(isetnf) == fsetno) GO TO 1204
 isetnf = nsetf + 1
 IF (isetnf < ivec) GO TO 1201
 fsetno = -1
 1202 DO  j = ilist,nlist,2
   z(j+1) = 1
 END DO
 GO TO 1207
 1204 DO  i = isetf,nsetf
   k    = 0
   diff = 1.e+25
   bufr(1) = zz(i)
   DO  j = ilist,nlist,2
     IF (z(j+1) /= 0) CYCLE
     diff1 = ABS(zz(j) - bufr(1))
     IF (diff1 >= diff) CYCLE
     diff = diff1
     k = j
   END DO
   IF (k /= 0) z(k+1) = 1
 END DO
 
!     DETERMINE IF CURRENT FREQ IS MARKED FOR OUTPUT.
 
 1207 IF (z(jlist+1) == 0) GO TO 1350
 buf(5) = z(jlist)
 buf(2) = kcount + 1014
 GO TO 1250
 
!     TRANSIENT RESPONSE.
 
 1220 buf(5) = z(jlist)
 buf(2) = kcount + 14
 IF (ireq == ipnl) buf(2) = 12
 ix     = icc + idload
 buf(8) = z(ix)
 buf(6) = 0
 buf(7) = 0
 
!     WRITE ID RECORD ON OUTPUT FILE.
 
 1250 ix = branch + 3
 IF (app(1) == cei(1)) ix = 9
 IF (app(1) == ireig ) ix = 2
 buf(1) = dest + 10*ix
 buf(3) = 0
 buf(4) = z(icc+1)
 IF (z(ireqx+2) < 0) sort2 = +1
 FORMAT  = IABS(z(ireqx+2))
 buf(9)  = FORMAT
 buf(10) = nwds
 CALL WRITE (outfl,buf,50,0)
 ix = icc + ittl
 CALL WRITE (outfl,z(ix),96,1)
 output = 1
 IF (z(ireqx+2) < 0) sort2 = 1
 
!     BUILD DATA RECORD ON OUTPUT FILE.
 
 IF (FORM(1) == modal(1)) GO TO 1270
 IF (setno /= -1) GO TO 1263
 
!     DIRECT PROBLEM SET .EQ. -ALL- - OUTPUT POINTS IN ANALYSIS SET
 
 kx = 1
 ASSIGN 1262 TO retx
 1261 word = z(kx+1)
 IF (word == 0) GO TO retx, (1262,1265,1268)
 j      = word/256
 buf(2) = andf(word,3)
 code   = word - 256*j - buf(2)
 buf(1) = z(kx)
 IF (buf(2) == 1) GO TO 1300
 GO TO 1290
 1262 kx = kx + 2
 IF (kx <= neqdyn) GO TO 1261
 GO TO 1310
 
!     DIRECT PROBLEM WITH SET .NE. -ALL- OUTPUT POINTS IN REQUESTED SET
!                                        WHICH ARE ALSO IN ANALYSIS SET.
 
 1263 jharm = 0
 1267 i = iset
 ASSIGN 1261 TO ret
 1264 buf(1) = z(i)
 IF (i   == nset) GO TO 1266
 IF (z(i+1) > 0) GO TO 1266
 n = -z(i+1)
 i = i + 1
 ASSIGN 1265 TO retx
 GO TO 3000
 1265 buf(1) = buf(1) + 1
 IF (buf(1) <= n) GO TO 3000
 GO TO 1268
 1266 ASSIGN 1268 TO retx
 GO TO 3000
 1268 i = i + 1
 IF (i <= nset) GO TO 1264
 IF (axif == 0) GO TO 1310
 jharm = jharm + 1
 IF (jharm <= oharms) GO TO 1267
 GO TO 1310
 
!     MODAL PROBLEM WITH SET .EQ. -ALL- OUTPUT ALL MODAL POINTS. THEN
!                                       IF EXTRA POINTS, OUTPUT THEM.
 
 1270 IF (setno /= -1) GO TO 1275
 buf(1) = imode
 buf(2) = 4
 j = 1
 ASSIGN 1271 TO retx
 GO TO 1290
 1271 buf(1) = buf(1) + 1
 j = buf(1) - imode + 1
 IF (buf(1) <= nbrmod) GO TO 1290
 IF (nbrep == 0) GO TO 1310
 kx = 1
 ASSIGN 1273 TO retx
 buf(2) = 3
 1272 j = z(kx+1)/10
 gptype = z(kx+1) - 10*j
 buf(1) = z(kx)
 IF (gptype == 3) GO TO 1290
 1273 kx = kx + 2
 IF (kx <= neqdyn) GO TO 1272
 GO TO 1310
 
!     MODAL PROBLEM WITH SET .NE. -ALL- ASSUME NUMBERS IN REQUESTED SET
!                                       WHICH ARE .LE. NO. OF MODES ARE
!                                       MODAL COORDINATES AND ANY OTHERS
!                                       ARE EXTRA POINTS.
 
 1275 jharm = 0
 1274 i = iset
 1276 buf(1) = z(i)
 IF (i   == nset) GO TO 1281
 IF (z(i+1) > 0) GO TO 1281
 n = -z(i+1)
 buf(2) = 4
 i = i + 1
 ASSIGN 1278 TO retx
 1277 IF (buf(1) < imode .OR. buf(1) > nbrmod) GO TO 1279
 j = buf(1) - imode + 1
 GO TO 1290
 1278 buf(1) = buf(1) + 1
 IF (buf(1) <= n) GO TO 1277
 GO TO 1284
 1279 IF (nbrep == 0) GO TO 1284
 ASSIGN 1280 TO ret
 buf(2) = 3
 GO TO 3000
 1280 j = z(kx+1)/10
 gptype = z(kx+1) - 10*j
 IF (gptype == 3) GO TO 1290
 GO TO 1278
 1281 ASSIGN 1284 TO retx
 IF (buf(1) < imode .OR. buf(1) > nbrmod) GO TO 1282
 ASSIGN 1284 TO retx
 j = buf(1) - imode + 1
 buf(2) = 4
 GO TO 1290
 1282 IF (nbrep == 0) GO TO 1284
 ASSIGN 1283 TO ret
 GO TO 3000
 1283 j = z(kx+1)/10
 buf(2) = z(kx+1) - 10*j
 IF (buf(2) == 3) GO TO 1290
 1284 i = i + 1
 IF (i <= nset) GO TO 1276
 IF (axif == 0) GO TO 1310
 jharm = jharm + 1
 IF (jharm <= oharms) GO TO 1274
 GO TO 1310
 
!     SCALAR, EXTRA OR MODAL POINT.
 
 1290 j = ivec + ktype*(j-1)
 bufr(3) = zz(j)
 DO  k = 4,nwds
   buf(k) = 0
 END DO
 IF (ktype == 1) GO TO 1309
 
!     COMPLEX SCALAR, EXTRA OR MODAL POINT.
 
 bufr(9) = zz(j+1)
 IF (FORMAT /= 3) GO TO 1309
 redner = SQRT(bufr(3)**2 + bufr(9)**2)
 IF (redner == 0.0) THEN
   GO TO  1309
 END IF
 12921 bufr(9) = ATAN2(bufr(9),bufr(3))*raddeg
 IF (bufr(9) < -0.00005) bufr(9) = bufr(9) + 360.0
 bufr(3) = redner
 GO TO 1309
 
!     GRID POINT.
 
 1300 DO  k = 3,nwds
   buf(k) = 1
 END DO
 j = ivec + ktype*(j-1)
 IF (ktype == 2) GO TO 1303
 DO  k = 1,6
   IF (andf(code,masks(k)) == 0) CYCLE
   bufr(k+2) = zz(j)
   j = j + 1
 END DO
 GO TO 1309
 
!     COMPLEX GRID POINT.
 
 1303 DO  k = 1,6
   IF (andf(code,masks(k)) == 0) CYCLE
   bufr(k+2) = zz(j  )
   bufr(k+8) = zz(j+1)
   j = j + 2
   IF (FORMAT /= 3) CYCLE
   redner = SQRT(bufr(k+2)**2 + bufr(k+8)**2)
   IF (redner == 0.0) THEN
     GO TO  1305
   END IF
   13031 bufr(k+8) = ATAN2(bufr(k+8),bufr(k+2))*raddeg
   IF (bufr(k+8) < -0.00005) bufr(k+8)= bufr(k+8) + 360.0
   bufr(k+2) = redner
   1305 CONTINUE
 END DO
 
!     DETERMINE DESTINATION FOR ENTRY.
 
 
!     IF A FLUID PROBLEM THEN A CHECK IS NOW MADE TO SEE IF THIS
!     HARMONIC IS TO BE OUTPUT
 
 1309 IF (axif == 0.0) THEN
   GO TO  1314
 END IF
 1315 IF (buf(1) < 500000) GO TO 1314
 itemp = buf(1) - MOD(buf(1),500000)
 itemp = itemp/500000
 IF (itemp >= oharms) GO TO 1310
 1314 id = buf(1)
 buf(1) = 10*id + dest
 IF (xsetno < 0.0) THEN
   GO TO  1304
 ELSE IF (xsetno == 0.0) THEN
   GO TO  1306
 ELSE
   GO TO  1307
 END IF
 1306 buf(1) = 10*id
 GO TO 1304
 1307 ix = ixset
 1313 IF (ix  == nxset) GO TO 1308
 IF (z(ix+1) > 0) GO TO 1308
 IF (id >= z(ix) .AND. id <= -z(ix+1)) GO TO 1304
 ix = ix + 2
 GO TO 1312
 1308 IF (id == z(ix)) GO TO 1304
 ix = ix + 1
 1312 IF (ix <= nxset) GO TO 1313
 GO TO 1306
 
!     WRITE ENTRY ON OUTPUT FILE.
 
 1304 CALL WRITE (outfl,buf,nwds,0)
 kwds = kwds + nwds
 buf(1) = id
 GO TO retx, (1262,1265,1268,1271,1273,1278,1284)
 
!     CONCLUDE PROCESSING OF THIS VECTOR.
 
 1310 CALL WRITE (outfl,0,0,1)
 1311 SELECT CASE ( branch )
   CASE (    1)
     GO TO 1340
   CASE (    2)
     GO TO 1350
   CASE (    3)
     GO TO 1360
   CASE (    4)
     GO TO 1340
 END SELECT
 
!     COMPLEX EIGENVALUES.
 
 1340 jlist = jlist + 3
 1341 IF (jcount >= nvects) GO TO 1410
 IF (eof == 0) GO TO 1130
 GO TO 1140
 
!     FREQUENCY RESPONSE.
 
 1350 IF (kcount == 3) GO TO 1356
 n = ivecn - 1
 omega = twopi*zz(jlist)
 DO  i = ivec,n,2
   bufr(1) = -omega*zz(i+1)
   zz(i+1) =  omega*zz(i  )
   zz(i  ) =  bufr(1)
 END DO
 IF (kcount == 2) GO TO 1352
 ireq = iavel
 GO TO 1353
 1352 ireq = iaacc
 1353 kcount = kcount + 1
 incore = 1
 GO TO 1140
 1356 kcount = 1
 incore = 0
 ireq   = iadisp
 jlist  = jlist + 2
 IF (jlist <= nlist .AND. jcount < nvects) GO TO 1140
 kfrq  = 0
 jlist = ilist
 DO  i = ilist,nlist,2
   z(i+1) = 0
 END DO
 IF (jcount < nvects) GO TO 1130
 GO TO 1410
 
!     TRANSIENT RESPONSE.
 
 1360 IF (ireq == ipnl) GO TO 1364
 IF (kcount-2 < 0) THEN
   GO TO  1361
 ELSE IF (kcount-2 == 0) THEN
   GO TO  1362
 ELSE
   GO TO  1363
 END IF
 1361 ireq   = iavel
 kcount = 2
 GO TO 1140
 1362 ireq   = iaacc
 kcount = 3
 GO TO 1140
 1363 ireq   = iadisp
 kcount = 1
 1364 jlist  = jlist + 1
 IF (jlist <= nlist .AND. jcount < nvects) GO TO 1140
 GO TO 1410
 
!     HERE WHEN EOF ENCOUNTERED ON CASE CONTROL.
 
 1400 eof = 1
 SELECT CASE ( branch )
   CASE (    1)
     GO TO 1341
   CASE (    2)
     GO TO 1410
   CASE (    3)
     GO TO 1410
   CASE (    4)
     GO TO 1341
 END SELECT
 
!     CONCLUDE PROCESSING.
 
 1410 CALL CLOSE (casecc,clsrew)
 CALL CLOSE (infil, clsrew)
 CALL CLOSE (outfl, clsrew)
 mcb(1) = outfl
 mcb(2) = kwds/65536
 mcb(3) = kwds - 65536*mcb(2)
 mcb(4) = 0
 mcb(5) = 0
 mcb(6) = 0
 mcb(7) = 0
 CALL wrttrl (mcb)
 RETURN
 
!     HERE IF ABNORMAL CONDITION.
 
 1430 CALL CLOSE (outfl,clsrew)
 1431 CALL mesage (30,78,0)
 1432 RETURN
 
!     FATAL FILE ERRORS
 
 2001 n = -1
 GO TO 2005
 2002 n = -2
 2005 CALL mesage (n,FILE,nam)
 RETURN
 
!     BINARY SEARCH ROUTINE
 
 3000 klo = 1
 khi = kn
 IF (axif == 0.0) THEN
   GO TO  3001
 END IF
 3011 buf(1) = jharm*500000 + buf(1)
 3001 k  = (klo+khi+1)/2
 3002 kx = 2*k - 1
 IF (buf(1)-z(kx) < 0.0) THEN
   GO TO  3003
 ELSE IF (buf(1)-z(kx) == 0.0) THEN
   GO TO  3009
 ELSE
   GO TO  3004
 END IF
 3003 khi = k
 GO TO 3005
 3004 klo = k
 3005 IF (khi-klo-1 < 0) THEN
   GO TO  3010
 ELSE IF (khi-klo-1 == 0) THEN
   GO TO  3006
 ELSE
   GO TO  3001
 END IF
 3006 IF (k == klo) GO TO 3007
 k = klo
 GO TO 3008
 3007 k = khi
 3008 klo = khi
 GO TO 3002
 3009 GO TO ret,  (1261,1280,1283)
 3010 GO TO retx, (1262,1265,1268,1273,1278,1284)
END SUBROUTINE vdrb
