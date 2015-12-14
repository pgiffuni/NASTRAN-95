SUBROUTINE exo2
     
!     EXO2 PERFORMS EXTERNAL FORMAT SOFOUT OPERATIONS
 
 EXTERNAL rshift   ,andf
 LOGICAL :: univac
 INTEGER :: dry      ,cor(1)   ,uname    ,TYPE     ,UNIT     ,  &
     itms(50) ,sysbuf   ,a        ,eol      ,eor      ,  &
     ditsiz   ,z        ,all      ,q4       ,t3       ,  &
     matric   ,tables   ,phase3   ,whole(2) ,subr(2)  ,  &
     BLANK    ,sof      ,xxxx     ,srd      ,prc      ,  &
     swrt     ,eog      ,eoi      ,sp       ,bar      ,  &
     scr1     ,buf1     ,buf2     ,buf3     ,eltype   ,  &
     buf4     ,rc       ,hdr(7)   ,typout   ,bdit     ,  &
     bmdi     ,rshift   ,andf     ,offset
 INTEGER :: eqss     ,bgss     ,cstm     ,lods     ,loap     ,  &
     plts     ,soln     ,lams
 DOUBLE PRECISION :: dz(1)    ,da
 CHARACTER (LEN=1) ::
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm swm*27
 COMMON  /xmssg /   ufm      ,uwm      ,uim      ,sfm      , swm
 COMMON  /BLANK /   dry      ,x1(3)    ,uname(2) ,x2(6)    ,  &
     TYPE(2)  ,names(10),UNIT     ,univac   , lbuf     ,iadd
 COMMON  /system/   sysbuf   ,nout     ,x3(6)    ,nlpp     ,  &
     x4(2)    ,line     ,x6(9)    ,mach
 COMMON  /names /   rd       ,rdrew    ,wrt      ,wrtrew   , rew
 COMMON  /TYPE  /   prc(2)   ,nword(4)
 COMMON  /zntpkx/   a(4)     ,irow     ,eol      ,eor
 COMMON  /sof   /   x5(3)    ,ditsiz
 COMMON  /itemdt/   nitem    ,items(7,1)
 COMMON  /zzzzzz/   z(1)
 EQUIVALENCE        (cor(1)  ,z(1))
 EQUIVALENCE        (z(1)    ,dz(1))   ,(a(1)    ,da)
 DATA     all      ,matric   ,tables   ,phase3   ,whole          /  &
     4HALL    ,4HMATR   ,4HTABL   ,4HPHAS   ,4HWHOL  ,4HESOF/,  &
     subr               ,BLANK    ,sof      ,xxxx           /  &
     4HEXO2   ,4H       ,4H       ,4HSOF    ,4HXXXX         /,  &
     eqss     ,bgss     ,cstm     ,lods     ,loap           /  &
     4HEQSS   ,4HBGSS   ,4HCSTM   ,4HLODS   ,4HLOAP         /,  &
     plts     ,soln     ,lams     ,q4       ,t3      ,bar   /  &
     4HPLTS   ,4HSOLN   ,4HLAMS   ,2HQ4     ,2HT3    ,2HBR  /,  &
     srd      ,swrt     ,more     ,eog      ,eoi     ,sp    /  &
     1        ,2        ,1        ,2        ,3       ,1     /,  &
     jh       ,scr1     ,bdit     ,bmdi                     /  &
     1        ,301      ,4HDIT    ,4HMDI                    /
 
!     INITIALIZE
 
 IF (nitem > 50) CALL errmkn (23,10)
 ncore = korsz(z)
 i     = ncore - lbuf
 IF (mach == 4) i = i - lbuf
 ncore = i - 1
 irw   = iadd
 iadd  = i
 CALL exfort (3,UNIT,0,0,irw,0,0)
 buf1  = ncore - sysbuf + 1
 buf2  = buf1  - sysbuf - 1
 buf3  = buf2  - sysbuf
 buf4  = buf3  - sysbuf
 ncore = buf4  - 1
 IF (buf4 <= 0) GO TO 9008
 CALL sofopn (z(buf1),z(buf2),z(buf3))
 
!     CONSTRUCT ARRAY OF NAMES OF ITEMS TO BE COPIED
 
 IF (TYPE(1) /= all) GO TO 10
 nitems = nitem
 DO  i = 1,nitem
   itms(i) = items(1,i)
 END DO
 GO TO 70
 10 IF (TYPE(1) /= tables) GO TO 20
 nitems = 0
 DO  i = 1,nitem
   IF (items(2,i) > 0) CYCLE
   nitems = nitems + 1
   itms(nitems) = items(1,i)
 END DO
 GO TO 70
 20 IF (TYPE(1) /= matric) GO TO 50
 nitems = 0
 DO  i = 1,nitem
   IF (items(2,i) <= 0) CYCLE
   nitems = nitems + 1
   itms(nitems) = items(1,i)
 END DO
 GO TO 70
 50 IF (TYPE(1) /= phase3) GO TO 60
 nitems = 0
 DO  i = 1,nitem
   IF (andf(items(7,i),8) == 0) CYCLE
   nitems = nitems + 1
   itms(nitems) = items(1,i)
 END DO
 GO TO 70
 60 nitems  = 2
 itms(1) = TYPE(1)
 itms(2) = TYPE(2)
 IF (itms(2) == BLANK) nitems = 1
 
!     PUT NAMES OF ALL SUBSTRUCTURES TO BE COPIED AT TOP OF OPEN CORE
 
 70 nss = 0
 IF (names(1) == whole(1) .AND. names(2) == whole(2)) GO TO 90
 DO  i = 1,9,2
   IF (names(i) == xxxx) CYCLE
   nss = nss + 1
   IF (2*nss > ncore) GO TO 9008
   z(2*nss-1) = names(i  )
   z(2*nss  ) = names(i+1)
 END DO
 GO TO 110
 90 n = ditsiz/2
 DO  i = 1,n
   CALL fdit (i,j)
   IF (cor(j) == BLANK) CYCLE
   nss = nss + 1
   IF (2*nss > ncore) GO TO 9008
   z(2*nss-1) = cor(j)
   z(2*nss  ) = cor(j+1)
 END DO
 110 icore  = 2*nss + 3
 lcore  = ncore - icore + 1
 idpcor = icore/2 + 1
 CALL page
 
!     WRITE OUT DIT AND MDI CONTROL WORDS
 
 n = ditsiz/2
 IF (6*n > lcore) GO TO 9008
 hdr(1) = bdit
 hdr(2) = BLANK
 hdr(3) = BLANK
 hdr(4) = 2
 hdr(5) = ditsiz
 hdr(6) = sp
 hdr(7) = eog
 CALL exfort (swrt,UNIT,jh,hdr,7,sp,0)
 DO  i = 1,n
   CALL fdit (i,j)
   z(icore+2*i-2) = cor(j)
   z(icore+2*i-1) = cor(j+1)
 END DO
 CALL exfort (swrt,UNIT,2,z(icore),ditsiz,sp,0)
 hdr(1) = bmdi
 hdr(4) = 10
 hdr(5) = 6*n
 hdr(7) = eoi
 CALL exfort (swrt,UNIT,jh,hdr,7,sp,0)
 k = icore
 DO  i = 1,n
   CALL fmdi (i,j)
   z(k  ) = rshift(cor(j+1),20)
   z(k+1) = andf(rshift(cor(j+1),10),1023)
   z(k+2) = andf(cor(j+1),1023)
   z(k+3) = andf(rshift(cor(j+2),20),1023)
   z(k+4) = andf(rshift(cor(j+2),10),1023)
   z(k+5) = andf(cor(j+2),1023)
   k = k + 6
 END DO
 CALL exfort (swrt,UNIT,10,z(icore),6*n,sp,0)
 
!     LOOP OVER ALL SUBSTRUCTURES AND ITEMS, COPYING EACH ONE TO THE
!     EXTERNAL FILE
 
 DO  iss = 1,nss
   hdr(1) = z(2*iss-1)
   hdr(2) = z(2*iss)
   DO  item = 1,nitems
     hdr(3) = itms(item)
     itm = ittype (itms(item))
     IF (itm == 1) GO TO 800
     CALL sfetch (hdr,hdr(3),srd,rc)
     SELECT CASE ( rc )
       CASE (    1)
         GO TO 140
       CASE (    2)
         GO TO 120
       CASE (    3)
         GO TO 990
       CASE (    4)
         GO TO 120
       CASE (    5)
         GO TO 120
     END SELECT
     120 line = line + 2
     IF (line > nlpp) CALL page
     IF (rc > 3) GO TO 130
     WRITE (nout,6340) uwm,(hdr(i),i=1,3)
     CYCLE
     130 CALL smsg (rc-2,hdr(3),hdr)
     CYCLE
     140 CALL suread (z(icore),lcore,nwds,rc)
     IF (rc /= 2) GO TO 9008
     
     IF (itms(item) == eqss) GO TO 200
     IF (itms(item) == bgss) GO TO 300
     IF (itms(item) == cstm) GO TO 400
     IF (itms(item) == lods) GO TO 500
     IF (itms(item) == loap) GO TO 500
     IF (itms(item) == plts) GO TO 600
     IF (itms(item) == soln) GO TO 700
     IF (itms(item) == lams) GO TO 700
     GO TO 1100
     
!     EQSS
     
!     GROUP 0
     
     200 n  = nwds
     ns = z(icore+2)
     IF (ns > 13) n = 30
     hdr(4) = 3
     hdr(5) = n
     hdr(6) = sp
     hdr(7) = eog
     IF (n < nwds) hdr(7) = more
     CALL exfort (swrt,UNIT,jh,hdr,7,sp,0)
     CALL exfort (swrt,UNIT,3,z(icore),n,sp,0)
     IF (n == nwds) GO TO 210
     hdr(4) = 2
     hdr(5) = nwds - n
     hdr(7) = eog
     CALL exfort (swrt,UNIT,jh,hdr,7,sp,0)
     CALL exfort (swrt,UNIT,2,z(icore+n),nwds-n,sp,0)
     
!     GROUPS 1 TO NS + 1
     
     210 hdr(4) = 10
     ns = ns + 1
     DO  j = 1,ns
       CALL suread (z(icore),lcore,nwds,rc)
       IF (rc /= 2) GO TO 9008
       hdr(5) = nwds
       IF (j == ns) hdr(7) = eoi
       CALL exfort (swrt,UNIT,jh,hdr,7,sp,0)
       CALL exfort (swrt,UNIT,10,z(icore),nwds,sp,0)
     END DO
     GO TO 900
     
!     BGSS
     
!     GROUP 0
     
     300 hdr(4) = 3
     hdr(5) = 3
     hdr(6) = sp
     hdr(7) = eog
     CALL exfort (swrt,UNIT,jh,hdr,7,sp,0)
     CALL exfort (swrt,UNIT,3,z(icore),3,sp,0)
     
!     GROUP 1
     
     CALL suread (z(icore),lcore,nwds,rc)
     IF (rc /= 2) GO TO 9008
     hdr(4) = 6
     hdr(5) = nwds
     hdr(7) = eoi
     CALL exfort (swrt,UNIT,jh,hdr,7,sp,0)
     CALL exfort (swrt,UNIT,6,z(icore),nwds,sp,0)
     GO TO 900
     
!     CSTM
     
!     GROUP 0
     
     400 hdr(4) = 3
     hdr(5) = 2
     hdr(6) = sp
     hdr(7) = eog
     CALL exfort (swrt,UNIT,jh,hdr,7,sp,0)
     CALL exfort (swrt,UNIT,3,z(icore),2,sp,0)
     
!     GROUP 1
     
     IF (icore+13 > ncore) GO TO 9008
     420 CALL suread (z(icore),14,nwds,rc)
     IF (rc == 2) GO TO 430
     hdr(4) = 8
     hdr(5) = 4
     hdr(7) = more
     CALL exfort (swrt,UNIT,jh,hdr,7,sp,0)
     CALL exfort (swrt,UNIT,8,z(icore),4,sp,0)
     hdr(4) = 9
     hdr(5) = 10
     CALL exfort (swrt,UNIT,jh,hdr,7,sp,0)
     CALL exfort (swrt,UNIT,9,z(icore+4),10,sp,0)
     GO TO 420
     430 hdr(5) = 0
     hdr(7) = eog
     CALL exfort (swrt,UNIT,jh,hdr,7,sp,0)
     hdr(4) = 0
     hdr(7) = eoi
     CALL exfort (swrt,UNIT,jh,hdr,7,sp,0)
     GO TO 900
     
!     LODS AND LOAP
     
!     GROUP 0
     
     500 n  = nwds
     ns = z(icore+3)
     IF (ns > 13) n = 30
     hdr(4) = 3
     hdr(5) = n
     hdr(6) = sp
     hdr(7) = eog
     IF (n < nwds) hdr(7) = more
     CALL exfort (swrt,UNIT,jh,hdr,7,sp,0)
     CALL exfort (swrt,UNIT,3,z(icore),n,sp,0)
     IF (n == nwds) GO TO 510
     hdr(4) = 2
     hdr(5) = nwds -  n
     hdr(7) = eog
     CALL exfort (swrt,UNIT,jh,hdr,7,sp,0)
     CALL exfort (swrt,UNIT,2,z(icore+n),nwds-n,sp,0)
     
!     GROUP 1 TO NS
     
     510 hdr(4) = 10
     DO  j = 1,ns
       CALL suread (z(icore),lcore,nwds,rc)
       IF (rc /= 2) GO TO 9008
       hdr(5) = nwds
       IF (j == ns) hdr(7) = eoi
       CALL exfort (swrt,UNIT,jh,hdr,7,sp,0)
       CALL exfort (swrt,UNIT,10,z(icore),nwds,sp,0)
     END DO
     GO TO 900
     
!     PLTS
     
!     GROUP 0
     
     600 n  = nwds
     ns = z(icore+2)
     hdr(6) = sp
     hdr(4) = 3
     hdr(5) = 3
     hdr(7) = more
     CALL exfort (swrt,UNIT,jh,hdr,7,sp,0)
     CALL exfort (swrt,UNIT,3,z(icore),3,sp,0)
     DO  j = 1,ns
       hdr(4) = 13
       hdr(5) = 4
       CALL exfort (swrt,UNIT,jh,hdr,7,sp,0)
       CALL exfort (swrt,UNIT,13,z(icore+14*j-11),4,sp,0)
       hdr(4) = 9
       hdr(5) = 10
       IF (j == ns) hdr(7) = eog
       CALL exfort (swrt,UNIT,jh,hdr,7,sp,0)
       CALL exfort (swrt,UNIT,9,z(icore+14*j-7),10,sp,0)
     END DO
     
!     GROUP 1 -- BGPDT
     
     CALL suread (z(icore),lcore,nwds,rc)
     IF (rc == 3) GO TO 680
     IF (rc /= 2) GO TO 9008
     hdr(4) = 6
     hdr(5) = nwds
     CALL exfort (swrt,UNIT,jh,hdr,7,sp,0)
     CALL exfort (swrt,UNIT,6,z(icore),nwds,sp,0)
     
!     GROUP 2 -- EQEXIN
     
     CALL suread (z(icore),lcore,nwds,rc)
     IF (rc /= 2) GO TO 9008
     hdr(4) = 10
     hdr(5) = nwds
     CALL exfort (swrt,UNIT,jh,hdr,7,sp,0)
     CALL exfort (swrt,UNIT,10,z(icore),nwds,sp,0)
     
!     GROUP 3 -- GPSETS
     
     CALL suread (z(icore),lcore,nwds,rc)
     IF (rc /= 2) GO TO 9008
     hdr(4) = 10
     hdr(5) = nwds
     CALL exfort (swrt,UNIT,jh,hdr,7,sp,0)
     CALL exfort (swrt,UNIT,10,z(icore),nwds,sp,0)
     
!     GROUP 4 -- ELSETS
     
!     OUTPUT CHANGES MADE BY G.CHAN/UNISYS   4/91
     
!     IN 90 AND EARLIER VERSIONS, ONLY ONE ELEMENT PLOT SYMBOL WORD WAS
!     WRITTEN OUT USING FORMAT 2, AND ON NEXT ELSETS DATA LINE, FORMAT
!     10 WAS USED FOR ALL ELEMENTS. NO OFFSET DATA WAS PROVIDED FOR THE
!     BAR, QUAD4 AND TRIA3 ELEMENTS. THE NO. OF GRID POINT PER ELEMENT,
!     NGPEL, WAS THE FIRST WORD ON THE ELSETS DATA LINE.  (LINE=RECORD)
!     ALSO, THE 90 AND EARLIER VERSIONS DID NOT COUNT PROPERTY ID, PID,
!     ON THE ELSETS DATA LINE. THUS THE TOTAL NO. OF WORDS MAY BE IN
!     ERROR AND MAY CAUSE EXTRA ZEROS TO APPEAR AT THE END OF THE LINE.
     
!     IN 91 VERSION, ELEMENT PLOT SYMBOL LINE HAS 2 WORDS, SYMBOL AND
!     NGPEL, AND FORMAT 25 IS USED. ON NEXT ELSETS DATA LINE, FORMAT 10
!     IS USED FOR ALL ELEMENTS WITH NO OFFSETS. FORMAT 26 IS USED FOR
!     THE BAR WHICH HAS 6 OFFSET VALUES, AND FORMATS 27 AND 28 ARE USED
!     FOR TRIA3 AND QUAD4 WHICH HAVE 1 OFFSET VALUE EACH. NOTE THAT
!     NGPEL HAS BEEN MOVED, AND IS NO LONGER THE FIRST WORD ON THE
!     ELSETS DATA LINE.
     
     hdr(7) = more
     
!     READ PLOT SYMBOL, AND NO. OF GRID POINTS PER ELEMENT
!     SET UP NO. OF OFFSET DATA FOR BAR, QUAD4 AND TRIA3
     
     640 CALL suread (z(icore),2,nwds,rc)
     IF (rc >= 2) GO TO 670
     hdr(4) = 25
     hdr(5) = 2
     ngpel  = z(icore+1)
     eltype = z(icore  )
     offset = 0
     IF (eltype == bar) offset = 6
     IF (eltype == q4 .OR. eltype == t3) offset = 1
     CALL exfort (swrt,UNIT,jh,hdr,7,sp,0)
     CALL exfort (swrt,UNIT,25,z(icore),2,sp,0)
     
!     READ ELEMENT ID NUMBER, PROPERTY ID, GRID POINT CONNECTION INDICES
!     AND OFFSETS IF THEY EXIST
!     (ERROR IN 90 AND EARLIER VERSIONS, PROPERTY ID WAS LEFT OUT, AND
!     THEREFORE DATA COUNT PER ELEMENT WAS INCORRECT)
     
     n = icore - ngpel - 2 - offset
     650 n = n     + ngpel + 2 + offset
     IF (n > ncore) GO TO 9008
     CALL suread (z(n),1,nwds,rc)
     IF (z(n) == 0) GO TO 655
     IF (n+ngpel+2+offset > ncore) GO TO 9008
     CALL suread (z(n+1),ngpel+1,nwds,rc)
     IF (offset /= 0) CALL suread (z(n+ngpel+2),offset,nwds,rc)
     SELECT CASE ( rc )
       CASE (    1)
         GO TO 650
       CASE (    2)
         GO TO 6100
       CASE (    3)
         GO TO 6100
     END SELECT
     
!     ALL ELEMENTS OF ONE TYPE READ INTO CORE, NOW COPY OUT
     
     655 hdr(5) = n - icore + 1
     IF (offset-1 < 0.0) THEN
       GO TO   660
     ELSE IF (offset-1 == 0.0) THEN
       GO TO   661
     ELSE
       GO TO   663
     END IF
!               REGULAR  QUAD4  BAR
!               ELEMENT  TRIA3
     
     660 hdr(4) = 10
     GO TO 665
     661 hdr(4) = 27
     IF (eltype == q4) hdr(4) = 28
     GO TO 665
     663 hdr(4) = 26
     665 CALL exfort (swrt,UNIT,jh,hdr,7,sp,0)
     CALL exfort (swrt,UNIT,hdr(4),z(icore),hdr(5),sp,0)
     GO TO 640
     
!     WRITE END-OF-ITEM FOR PLTS
     
     670 hdr(5) = 0
     hdr(7) = eog
     CALL exfort (swrt,UNIT,jh,hdr,7,sp,0)
     680 hdr(4) = 0
     hdr(5) = 0
     hdr(7) = eoi
     CALL exfort (swrt,UNIT,jh,hdr,7,sp,0)
     GO TO 900
     
!     SOLN AND LAMS
     
     700 irfno = z(icore+2)
     IF (irfno == 1) GO TO 715
     IF (irfno == 2) GO TO 715
     IF (irfno == 3) GO TO 730
     IF (irfno == 8) GO TO 750
     IF (irfno == 9) GO TO 750
     line = line + 2
     IF (line > nlpp) CALL page
     WRITE (nout,6358) swm,irfno,hdr(1),hdr(2)
     GO TO 900
     
!     GROUP 0 -- STATICS
     
     715 n  = nwds
     ns = z(icore+3)
     IF (ns > 6) n = 23
     ns = z(icore+4)
     hdr(4) = 16
     hdr(5) = n
     hdr(6) = sp
     hdr(7) = eog
     IF (n < nwds) hdr(7) = more
     CALL exfort (swrt,UNIT,jh,hdr,7,sp,0)
     CALL exfort (swrt,UNIT,16,z(icore),n,sp,0)
     IF (n == nwds) GO TO 710
     hdr(4) = 17
     hdr(5) = nwds - n
     hdr(7) = eog
     CALL exfort (swrt,UNIT,jh,hdr,7,sp,0)
     CALL exfort (swrt,UNIT,17,z(icore+n),nwds-n,sp,0)
     
!     GROUPS 1 TO NS (ONE PER SUBCASE) -- STATICS
     
     710 DO  j = 1,ns
       CALL suread (z(icore),lcore,nwds,rc)
       IF (rc /= 2) GO TO 9008
       n = nwds
       IF (z(icore) > 5) n = 11
       hdr(4) = 18
       hdr(5) = n
       hdr(7) = eog
       IF (j == ns) hdr(7) = eoi
       IF (n < nwds) hdr(7) = more
       CALL exfort (swrt,UNIT,jh,hdr,7,sp,0)
       CALL exfort (swrt,UNIT,18,z(icore),n,sp,0)
       IF (n == nwds) CYCLE
       hdr(4) = 19
       hdr(5) = nwds - n
       hdr(7) = eog
       IF (j == ns) hdr(7) = eoi
       CALL exfort (swrt,UNIT,jh,hdr,7,sp,0)
       CALL exfort (swrt,UNIT,19,z(icore+n),nwds-n,sp,0)
     END DO
     GO TO 900
     
!     GROUP 0 -- NORMAL MODES (REAL OR COMPLEX)
     
     730 ns = z(icore+3)
     hdr(4) = 3
     hdr(5) = 4
     hdr(6) = sp
     hdr(7) = eog
     IF (ns <= 0) hdr(7) = eoi
     CALL exfort (swrt,UNIT,jh,hdr,7,sp,0)
     CALL exfort (swrt,UNIT,3,z(icore),4,sp,0)
     IF (ns <= 0) GO TO 900
     
!     GROUP 1 -- NORMAL MODES
     
     CALL suread (z(icore),lcore,nwds,rc)
     IF (rc /= 2) GO TO 9008
     hdr(4) = 20
     hdr(5) = nwds
     hdr(7) = eoi
     IF (itms(item) == lams) hdr(7) = eog
     CALL exfort (swrt,UNIT,jh,hdr,7,sp,0)
     CALL exfort (swrt,UNIT,20,z(icore),nwds,sp,0)
     IF (itms(item) /= lams) GO TO 900
     
!     GROUP 2 -- NORMAL MODES (LAMS ITEM ONLY)
     
     CALL suread (z(icore),lcore,nwds,rc)
     IF (rc /= 2) GO TO 9008
     hdr(4) = 10
     hdr(5) = nwds
     hdr(7) = eoi
     CALL exfort (swrt,UNIT,jh,hdr,7,sp,0)
     CALL exfort (swrt,UNIT,10,z(icore),nwds,sp,0)
     GO TO 900
     
!     GROUP 0 -- DYNAMICS
     
     750 ns    = z(icore+3)
     nwds0 = 3*ns + 5
     n     = nwds0
     IF (ns > 6) n = 23
     ns    = z(icore+4) + 1
     IF (z(icore+nwds0) == 0) ns = 1
     hdr(4) = 16
     hdr(5) = n
     hdr(6) = sp
     hdr(7) = more
     CALL exfort (swrt,UNIT,jh,hdr,7,sp,0)
     CALL exfort (swrt,UNIT,16,z(icore),n,sp,0)
     IF (n == nwds0) GO TO 760
     hdr(4) = 17
     hdr(5) = nwds0 - n
     CALL exfort (swrt,UNIT,jh,hdr,7,sp,0)
     CALL exfort (swrt,UNIT,17,z(icore+n),nwds0-n,sp,0)
     760 hdr(4) = 10
     hdr(5) = nwds - nwds0
     hdr(7) = eog
     CALL exfort (swrt,UNIT,jh,hdr,7,sp,0)
     CALL exfort (swrt,UNIT,10,z(icore+nwds0),nwds-nwds0,sp,0)
     
!     GROUP 1 TO NS+1 -- DYNAMICS
     
     DO  j = 1,ns
       CALL suread (z(icore),lcore,nwds,rc)
       IF (rc /= 2) GO TO 9008
       hdr(4) = 9
       hdr(5) = nwds
       hdr(7) = eog
       IF (j == ns) hdr(7) = eoi
       CALL exfort (swrt,UNIT,jh,hdr,7,sp,0)
       CALL exfort (swrt,UNIT,9,z(icore),nwds,sp,0)
     END DO
     GO TO 900
     
!     UNKNOWN TABLE ITME
     
     1100 line = line + 2
     IF (line > nlpp) CALL page
     WRITE (nout,6360) swm,itms(item)
     CYCLE
     
!     MATRICES
     
!     ON CDC MACHINE (NOT ANY 64-BIT MACHINE), FORCE ALL MATRIX DATA TO
!     BE DOUBLE PRECISION SO THE EXTRA DIGITS WONT BE LOST GOING TO
!     OTHER MACHINES
     
!     GROUP 0 -- MATRIX TRAILER
     
     800 CALL softrl (hdr,hdr(3),z(icore-1))
     rc = z(icore-1)
     SELECT CASE ( rc )
       CASE (    1)
         GO TO 805
       CASE (    2)
         GO TO 120
       CASE (    3)
         GO TO 990
       CASE (    4)
         GO TO 120
       CASE (    5)
         GO TO 120
     END SELECT
     805 typout = z(icore+3)
     IF (mach == 4 .AND. prc(typout) == 1) typout = typout + 1
     z(icore+3) = typout
     ncol   = z(icore)
     hdr(4) = 10
     hdr(5) = 6
     hdr(6) = sp
     hdr(7) = eog
     CALL exfort (swrt,UNIT,jh,hdr,7,sp,0)
     CALL exfort (swrt,UNIT,10,z(icore),6,sp,0)
     
!     MOVE MATRIX TO SCR2
     
     CALL mtrxi (scr1,hdr,hdr(3),z(buf4),rc)
     CALL gopen (scr1,z(buf4),rdrew)
     
!     COPY MATRIX OUT ONE COLUMN AT A TIME, NON-ZEROES ONLY.
     
!                        ROW NO.  +
!                        VALUE     +
!                        ROW NO.    +
!                        VALUE       I  FORMAT OF ONE MATRIX
!                          .         I  COLUMN ON THE EXTERNAL
!                          .         I  FILE.
!                          .        +
!                        -1        +
!                        0.0      +
     
     hdr(4) = 20 + typout
     hdr(6) = typout
     iprc   = prc(typout)
     n  = nword(typout) + iprc
     n2 = nword(typout) + 1
     DO  j = 1,ncol
       nwds = 0
       k  = icore
       CALL intpk (*820,scr1,0,typout,0)
       810 CALL zntpki
       z(k     ) = irow
       z(k+iprc) = a(1)
       IF (typout == 1) GO TO 815
       z(k+iprc+1) = a(2)
       IF (typout <= 3) GO TO 815
       z(k+4) = a(3)
       z(k+5) = a(4)
       815 nwds = nwds + n2
       k = k + n
       IF (k+n > ncore) GO TO 9008
       IF (eol ==     0) GO TO 810
       820 z(k) = -1
       z(k+iprc  ) = 0
       z(k+iprc+1) = 0
       z(k+4) = 0
       z(k+5) = 0
       nwds   = nwds + n2
       hdr(5) = nwds
       IF (j == ncol) hdr(7) = eoi
       CALL exfort (swrt,UNIT,jh,hdr,7,sp,0)
       CALL exfort (swrt,UNIT,20+typout,z(icore),nwds,typout,dz(idpcor))
     END DO
     CALL CLOSE (scr1,rew)
     
!     WRITE USER MESSAGE FOR SUCCESSFUL COPY
     
     900 line = line + 1
     IF (line > nlpp) CALL page
     WRITE (nout,6357) uim,hdr(1),hdr(2),hdr(3),sof,uname
   END DO
 END DO
 
!     NORMAL MODULE COMPLETION.  WRITE LOGICAL EOF
 
 CALL exfort (4,UNIT,0,0,1,0,0)
 CALL sofcls
 RETURN
 
!     ABNORMAL MODULE COMPLETION
 
 6100 CALL smsg (rc+4,itms(item),hdr)
 GO TO 9100
 9008 CALL mesage (8,0,subr)
 9100 dry = -2
 CALL sofcls
 RETURN
 
!     MESSAGE TEXT
 
 6340 FORMAT (a25,' 6340, SUBSTRUCTURE ',2A4,' ITEM ',a4, /5X,  &
     ' PSEUDO-EXISTS ONLY AND CANNOT BE COPIED OUT BY EXIO.')
 6357 FORMAT (a29,' 6357, SUBSTRUCTURE ',2A4,' ITEM ',a4,  &
     ' SUCCESSFULLY COPIED FROM ',a4,' TO ',2A4)
 6358 FORMAT (a27,' 6358, ILLEGAL RIGID FORMAT NUMBER ',i5,  &
     ' IN SOLN ITEM FOR SUBSTRUCTURE ',2A4,1H.,  &
     /34X,'THE ITEM WILL NOT BE COPIED.')
 6360 FORMAT (a27,' 6360, SOFOUT (EXTERNAL) ENCOUNTERS A UNSUPPORTED ',  &
     'TABLE ITEM ',a4, /35X,'THE ITEM WILL NOT BE COPIED.')
END SUBROUTINE exo2
