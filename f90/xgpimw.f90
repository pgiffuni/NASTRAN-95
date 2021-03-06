SUBROUTINE xgpimw (msgno,i,j,a)
     
!     XGPIMW IS USED TO WRITE ALL NON-DIAGNOSTIC MESSAGES
!     GENERATED BY THE XGPI MODULE.
 
 
 INTEGER, INTENT(IN)                      :: msgno
 INTEGER, INTENT(IN)                      :: i
 INTEGER, INTENT(OUT)                     :: j
 INTEGER, INTENT(IN)                      :: a(6)
 EXTERNAL        rshift,andf
 DIMENSION       hdg1(18),hdg2(7),hdg3(2),hdg4(22),hdg5(26), line(12)
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 INTEGER            :: op
 
 COMMON /xmssg / ufm,uwm,uim
 COMMON /system/ zsys(90),lpch,ldict
 COMMON /output/ pghdg(1)
 COMMON /xgpic / icold,islsh,iequl,nblank,nxequi,  &
                 ndiag,nsol,ndmap,nestm1,nestm2,nexit,  &
                 nbegin,nend,njump,ncond,nrept,ntime,nsave,noutpt,  &
                 nchkpt,npurge,nequiv,ncpw,nbpc,nwpc,  &
                 maskhi,masklo,isgnon,nosgn,iallon,masks(1)
 COMMON /xgpi6 / medtp,fnmtp,cnmtp,medpnt,lmed,iplus,diag14,diag17
 COMMON /xgpid / icst,iunst,imst,ihapp,idsapp,idmapp,  &
     isave,itpflg,iapnd,idum(5),ieqflg
 COMMON /xgpi5 / iapp,start,alter(2),sol,subset,iflag,iestim,  &
     icftop,icfpnt,lctlfl,ictlfl(1)
 COMMON /moddmp/ iflg(6),namopt(26)
 COMMON /resdic/ irdict, iropen
 EQUIVALENCE     (zsys( 2),op    ), (zsys( 3),nogo  ),  &
     (zsys( 9),nlpp  ), (zsys(12),nlines),  &
     (zsys(77),bkdata), (zsys(26),cppgct), (zsys(19),iecho ), (zsys(24),icfiat)
 DATA    hdg1  / 17, 4H    ,4H COS,4HMIC ,4H/ NA,4HSTRA,4HN DM,  &
             4HAP C,4HOMPI,4HLER ,4H- SO,4HURCE,4H LIS, 4HTING,4*4H     /
 DATA    hdg2  / 6, 6*4H         /
 DATA    hdg3  / 1, 4H           /
 DATA    hdg4  / 21, 4HINTE,4HRPRE,4HTED ,4HFROM,4H THE,4H OSC,4HAR. ,  &
                     4HNEGA,4HTIVE,4H DMA,4HP IN,4HDICA,4HTES ,4HA NO,  &
                     4HN EX,4HECUT,4HABLE,4H INS,4HTRUC,4HTION,4H    /
 DATA    hdg5  / 25, 4H* * ,4H*  D,4H M A,4H P  ,4H C R,4H O S,4H S -,  &
                 4H R E,4H F E,4H R E,4H N C,4H E  ,4H * *,4H *  ,  &
                 7*4H    ,4HMODU,4HLE  ,4H  NA,4HMES  /
 DATA    kdlh   / 0     /
 DATA    istar  / 1H*   /
 DATA    ipage  / 0     /
 DATA    iblnk  / 4H    /
 
 
 diagi4 = diag14
 IF (diag14 == 10) diagi4 = 0
 IF (diagi4 /= 0 .OR. diag17 /= 0) GO TO 100
 IF (iecho == -2 .AND. msgno /= 9) RETURN
 GO TO 110
 100  IF (iecho == -2 .AND. (msgno <= 4 .OR. msgno == 8) .AND.  &
     msgno /= 9) RETURN
 110  IF (msgno <= 0 .OR. msgno > 13) GO TO 2240
 SELECT CASE ( msgno )
   CASE (    1)
     GO TO 1860
   CASE (    2)
     GO TO 1960
   CASE (    3)
     GO TO 2000
   CASE (    4)
     GO TO 2020
   CASE (    5)
     GO TO 2040
   CASE (    6)
     GO TO 2060
   CASE (    7)
     GO TO 2080
   CASE (    8)
     GO TO 2100
   CASE (    9)
     GO TO 2160
   CASE (   10)
     GO TO 2080
   CASE (   11)
     GO TO 1985
   CASE (   12)
     GO TO 1970
   CASE (   13)
     GO TO 2400
 END SELECT
 
!     MESSAGE 1 - INITIALIZE PAGE HEADING
!     ===================================
 
 1860 i1 = -i
 IF (i1 >= 0) GO TO 1940
 i1 = hdg1(1)
 i2 = i1 + 1
 DO  m = 1,i1
   pghdg(m+96) = hdg1(m+1)
 END DO
 DO  m = i2,32
   pghdg(m+96) = nblank
 END DO
 IF (i == 1) GO TO 1888
 IF (i == 2) GO TO 1884
 i1 = hdg5(1)
 DO  m = 1,i1
   pghdg(m+128) = hdg5(m+1)
 END DO
 GO TO 1886
 1884 i1 = hdg4(1)
 DO  m = 1,i1
   pghdg(m+128) = hdg4(m+1)
 END DO
 1886 i1 = i1 + 32
 GO TO 1940
 1888 i1 = hdg2(1)
 i2 = i1 + 1
 DO  m = 1,i1
   pghdg(m+128) = hdg2(m+1)
 END DO
 DO  m = i2,32
   pghdg(m+128) = nblank
 END DO
 i1 = hdg3(1)
 i2 = i1 + 1
 DO  m = 1,i1
   pghdg(m+160) = hdg3(m+1)
 END DO
 DO  m = i2,32
   pghdg(m+160) = nblank
 END DO
 IF (diagi4 == 0 .AND. start == icst) GO TO 3000
 ipage = 1
 CALL page
 GO TO 3000
 
!     BLANK OUT HEADING
 
 1940 i2 = i1 + 1
 DO  m = i2,96
   pghdg(m+ 96) = nblank
 END DO
 nlines = nlpp
 GO TO 3000
 
!     MESSAGE 2  (XGPI)
!     =================
 
 1960 IF (start == iunst) GO TO 1980
 CALL page2 (-6)
 WRITE  (op,1965)
 1965 FORMAT (1H0,/,'0  + INDICATES DMAP INSTRUCTIONS THAT ARE PROCESS',  &
     'ED ONLY AT DMAP COMPILATION TIME.')
 WRITE  (op,1967)
 1967 FORMAT ('0  * INDICATES DMAP INSTRUCTIONS THAT ARE FLAGGED FOR ',  &
     'EXECUTION IN THIS MODIFIED RESTART.')
 GO TO 3000
 
!     MESSAGE 12
!     ==========
 
 1970 CALL page2 (-5)
 WRITE  (op,1975) uim
 1975 FORMAT (a29,' 4147', /5X,'NOTE THAT ADDITIONAL DMAP INSTRUCTIONS',  &
     ' (NOT INDICATED BY AN * IN THE DMAP SOURCE LISTING)',/5X,  &
     'NEED TO BE FLAGGED FOR EXECUTION IN ORDER TO GENERATE ',  &
     'CERTAIN REQUIRED DATA BLOCKS.', /5X,  &
     'SUCH INSTRUCTIONS AND THE ASSOCIATED DATA BLOCKS ARE ', 'IDENTIFIED BELOW.')
 GO TO 3000
 
 1980 CALL page2 (-6)
 WRITE  (op,1965)
 WRITE  (op,1982)
 1982 FORMAT ('0  * INDICATES DMAP INSTRUCTIONS THAT ARE FLAGGED FOR ',  &
     'EXECUTION IN THIS UNMODIFIED RESTART.')
 GO TO 3000
 
!     MESSAGE 11
!     ==========
 
 1985 CALL page2 (-6)
 WRITE  (op,1988) uim,i
 1988 FORMAT (a29,' 4148', /5X,'NOTE THAT ADDITIONAL DMAP INSTRUCTIONS',  &
     ' (NOT INDICATED BY AN * IN THE DMAP SOURCE LISTING)',/5X,  &
     'NEED TO BE FLAGGED FOR EXECUTION SINCE THIS UNMODIFIED ',  &
     'RESTART INVOLVES DMAP LOOPING AND', /5X,  &
     'THE REENTRY POINT IS WITHIN A DMAP LOOP.  SUCH INSTRUCT',  &
     'IONS ARE IDENTIFIED BELOW.', /5X,  &
     'THE EXECUTION WILL, HOWEVER, RESUME AT THE LAST REENTRY',  &
     ' POINT (DMAP INSTRUCTION NO.',i5,2H).,/)
 GO TO 3000
 
!     MESSAGE 3  (KYXFLD)
!     ===================
 
 2000 nlines = nlines + 2
 IF (nlines >= nlpp) CALL page1
 WRITE  (op,2010) i,j
 2010 FORMAT ('0TO GENERATE DATA BLOCK ',2A4,' - TURN ON THE EXECUTE ',  &
     'FLAG FOR THE FOLLOWING DMAP INSTRUCTIONS',/)
 GO TO 3000
 
!     MESSAGE 4  (KYXFLD)
!     ===================
 
 2020 nlines = nlines + 1
 IF (nlines >= nlpp) CALL page1
 l = andf(a(6),maskhi)
 WRITE  (op,2030) l,a(4),a(5)
 2030 FORMAT (1X,i4,2X,2A4)
 GO TO 3000
 
!     MESSAGE 5 - WRITE DMAP INSTRUCTION, FIRST LINE   (XSCNDM)
!     =========================================================
 
 2040 IF (diagi4 == 0) GO TO 2050
 IF (kdlh   == 0) GO TO 2255
 2042 CALL page2 (-2)
 WRITE  (op,2045) j,(a(m),m=1,i)
 2045 FORMAT (/1X,i7,2X,20A4)
 2050 IF (diag17 /= 0) WRITE (lpch,2055) (a(m),m=1,i)
 2055 FORMAT (20A4)
 GO TO 3000
 
!     MESSAGE 6 - WRITE DMAP INSTRUCTION, CONTINUATION LINE   (XSCNDM)
!     ================================================================
 
 2060 IF (diagi4 == 0) GO TO 2070
 nlines = nlines + 1
 IF (nlines >= nlpp) CALL page
 WRITE  (op,2065) (a(m),m=1,i)
 2065 FORMAT (10X,20A4)
 2070 IF (diag17 /= 0) WRITE (lpch,2075) (a(m),m=1,i)
 2075 FORMAT (20A4)
 GO TO 3000
 
!     MESSAGE 7   (XLNKHD)
!     MESSAGE 10  (XLNKHD AND XOSGEN)
!     ===============================
 
 2080 CONTINUE
 IF (diagi4 == 0) GO TO 3000
 iwrite = istar
 IF (msgno == 10) iwrite = iplus
 WRITE  (op,2090) iplus,iwrite
 2090 FORMAT (a1,2X,a1)
 GO TO 3000
 
!     MESSAGE 8  (KYGPI)
!     ==================
 
 2100 CALL page1
 nlines = nlines + 4
 WRITE  (op,2110)
 2110 FORMAT ('0THE FOLLOWING FILES FROM THE OLD PROBLEM TAPE WERE ',  &
     'USED TO INITIATE RESTART', //4X, 'FILE NAME  REEL NO.  FILE NO.',/)
 DO  m = i,j,3
   ireel = andf(rshift(a(m+2),16),31)
   ifile = andf(a(m+2),maskhi)
   CALL page2 (-1)
   IF (ifile == 0) GO TO 2130
   WRITE  (op,2120) a(m),a(m+1),ireel,ifile
   2120 FORMAT (5X,2A4,2X,i8,2X,i8)
   CYCLE
   2130 WRITE  (op,2140) a(m),a(m+1)
   2140 FORMAT (5X,2A4,2X,8H(purged))
 END DO
 GO TO 3000
 
!     MESSAGE 9  (KYGPI)
!     ==================
 
 2160 IF (bkdata < 0) CALL page2 (-2)
 IF (bkdata == -2) GO TO 2171
 IF (iflg(1) == 1) GO TO 2165
 WRITE  (op,2161)
 2161 FORMAT (1H0,10X,'**COMPILATION COMPLETE - JOB TERMINATING WITH ',  &
     'NOGO STATUS**')
 GO TO 2175
 2165 IF (bkdata >= 0) GO TO 2175
 WRITE  (op,2170)
 2170 FORMAT (1H0,10X,'**NO ERRORS FOUND - EXECUTE NASTRAN PROGRAM**')
 GO TO 2175
 2171 WRITE  (op,2172)
 2172 FORMAT (1H0,9X,'**NO ERRORS FOUND - NASTRAN EXECUTION TERMINATED',  &
     ' BY USER REQUEST**')
 GO TO 3000
 
!     DUMP FIAT IF SENSE SWITCH 2 IS ON
 
 2175 CALL sswtch (2,l)
 IF (l == 0) GO TO 2210
 CALL page1
 nlines = nlines + 4
 WRITE  (op,2180) a(1),a(2),a(3)
2180 FORMAT (1H ,/5X,22HFIAT at END of preface,3I15, /,1H ,/5X,'EQUIV',  &
'  APPEND    LTU  TAPE  UNIT  FILE NAME',31X,'---TRAILER---')
l1 = a(3)*icfiat - 2
DO  l = 4,l1,icfiat
  iequiv = 0
  iappnd = 0
  itape  = 0
  IF (rshift(andf(a(l),ieqflg),1) > 0) iequiv = 1
  IF (andf(a(l),iapnd ) > 0) iappnd = 1
  IF (andf(a(l),itpflg) > 0) itape  = 1
  ltu = andf(rshift(a(l),16),16383)
  iunit = andf(a(l),16383)
  ia1 = a(l+1)
  ia2 = a(l+2)
  IF (ia1 /= 0) GO TO 2185
  ia1 = iblnk
  ia2 = iblnk
  2185 m1  = l + 3
  m2  = l + 5
  CALL page2 (-1)
  WRITE (op,2190) iequiv,iappnd,ltu,itape,iunit,ia1,ia2
  IF (icfiat ==  8) WRITE (op,2191) (a(m),m=m1,m2)
  IF (icfiat == 11) WRITE (op,2192) (a(m),m=m1,m2),(a(m+m2),m=3,5)
  2190 FORMAT (7X,i1,7X,i1,3X,i6,4X,i1,1X,i6,2X,2A4)
  2191 FORMAT (1H+,48X,3I20)
  2192 FORMAT (1H+,48X,6I10)
END DO
2210 IF (j == 0) GO TO 3000
IF (iropen == 1) GO TO 2215
!WKBD OPEN (UNIT=4, FILE='FORTDIC.ZAP', STATUS='UNKNOWN')
iropen = 1
2215 WRITE  (irdict,2220) i
2220 FORMAT (9X,'1,   XVPS    ,   FLAGS = 0,   REEL =  1,   FILE =',i7)
CALL sswtch (9,diag09)
IF (diag09 == 1) GO TO 3000
CALL page1
nlines = nlines + 3
WRITE  (op,2230)
2230 FORMAT (9X,'CONTINUATION OF CHECKPOINT DICTIONARY', /,1H )
WRITE  (op,2220) i
GO TO 3000
2240 nlines = 2 + nlines
IF (nlines >= nlpp) CALL page1
WRITE  (op,2250) msgno
2250 FORMAT (//,' NO MESSAGE AVAILABLE FOR MESSAGE NO. ',i4)
GO TO 3000

!     PROCESS DMAP COMPILER OPTIONS SUMMARY

2255 IF (kdlh /= 0) GO TO 3000
IF (ipage == 0) CALL page
nlines   = nlines + 4
line( 1) = namopt( 1)
line( 2) = namopt( 2)
line( 3) = namopt( 5)
line( 4) = iflg  ( 2)
line( 5) = namopt( 9)
line( 6) = namopt(10)
line( 7) = namopt(13)
line( 8) = namopt(14)
line( 9) = namopt(17)
line(10) = namopt(18)
line(11) = namopt(21)
line(12) = namopt(22)
IF (iflg(1) <= 0) line(1) = namopt(3)
IF (iflg(3) <= 0) GO TO 2300
line( 5) = namopt( 7)
line( 6) = namopt( 8)
2300 IF (iflg(4) <= 0) GO TO 2310
line( 7) = namopt(11)
line( 8) = namopt(12)
2310 IF (iflg(5) <= 0) GO TO 2330
line( 9) = namopt(15)
line(10) = namopt(16)
2330 IF (iflg(6) <= 0) GO TO 2345
line(11) = namopt(19)
line(12) = namopt(20)
2345 WRITE  (op,2350) (line(kkdj),kkdj=1,12)
2350 FORMAT ('0  OPTIONS IN EFFECT ',2A4,a3,1H=,i1,3X,8A4,/3X,17(1H-), /)
kdlh = 1
GO TO 2042

!     MESSAGE 13   (XGPI)
!     ===================

2400 nscr = 315
l    = 20
CALL OPEN (*3000,nscr,a(l),0)
CALL page1
CALL page3 (4)
WRITE  (op,2410) uim
2410 FORMAT (a29,' - DUE TO ERROR(S), POSSIBLY ORIGINATED FROM USER''S'  &
    ,      ' ALTER PACKAGE, THE UNMODIFIED RIGID FORMAT LISTING', /5X,  &
    'IS PRINTED FOR CROSS REFERENCE',/)
l = 0
2420 CALL READ (*2440,*2440,nscr,a(1),18,0,m)
IF (a(1) == iblnk) GO TO 2430
CALL page3 (2)
l = l + 1
WRITE (op,2045) l,(a(m),m=1,18)
GO TO 2420
2430 nlines = nlines + 1
IF (nlines >= nlpp) CALL page
WRITE (op,2065) (a(m),m=1,18)
GO TO 2420
2440 CALL CLOSE (nscr,1)

3000 RETURN
END SUBROUTINE xgpimw
