SUBROUTINE exi2
     
!     EXI2 PERFORMS EXTERNAL FORMAT SOFIN OPERATIONS
 
 EXTERNAL lshift
 LOGICAL :: usrmsg
 INTEGER :: dry      ,uname    ,UNIT     ,sysbuf   ,a        ,  &
     z        ,sof      ,prc      ,q4       ,t3       ,  &
     srd      ,swrt     ,eoi      ,sp       ,bar      ,  &
     scr1     ,subr(2)  ,buf1     ,buf2     ,buf3     ,  &
     buf4     ,hdr(7)   ,rc       ,mcb(7)   ,prec     ,  &
     dit      ,NAME(2)  ,eog      ,plts     ,offset
 REAL :: zero(6)
 DOUBLE PRECISION :: dz(1)    ,da
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON  /xmssg /   ufm      ,uwm      ,uim
 COMMON  /machin/   mach
 COMMON  /BLANK /   dry      ,x1(3)    ,uname(2) ,x2(18)   ,  &
     UNIT     ,univac   ,lbuf     ,iadd
 COMMON  /system/   sysbuf   ,nout     ,x3(6)    ,nlpp     , x4(2)    ,line
 COMMON  /zblpkx/   a(4)     ,irow
 COMMON  /names /   rd       ,rdrew    ,wrt      ,wrtrew   , rew
 COMMON  /TYPE  /   prc(2)   ,nword(4)
 COMMON  /zzzzzz/   z(1)
 EQUIVALENCE       (z(1),dz(1))        ,(a(1),da)
 DATA     sof      ,srd      ,swrt     ,eoi      ,sp       /  &
     4HSOF    ,1        ,2        ,3        ,1        /
 DATA     leof     ,jh       ,scr1     ,subr               /  &
     4H$eof   ,1        ,301      ,4HEXI2   ,4H       /
 DATA     dit      ,mdi      ,eog      ,zero               /  &
     4HDIT    ,4HMDI    ,2        ,6*0.0              /
 DATA     q4       ,t3       ,bar      ,plts               /  &
     2HQ4     ,2HT3     ,2HBR     ,4HPLTS             /
 
!     INITIALIZE
 
 ncore = korsz(z)
 i     = ncore - lbuf
 IF (mach == 12) i = i - lbuf
 ncore = i - 1
 irw   = iadd
 iadd  = i
 CALL exfort (3,UNIT,0,0,irw,0,0)
 buf1  = ncore - sysbuf + 1
 buf2  = buf1  - sysbuf - 1
 buf3  = buf2  - sysbuf
 buf4  = buf3  - sysbuf
 ncore = buf4  - 1
 nos   = 0
 idm   = 1
 usrmsg=.true.
 lcore = ncore
 IF (ncore <= 0) GO TO 9008
 CALL sofopn (z(buf1),z(buf2),z(buf3))
 CALL page
 
!     READ THE HEADER OF THE NEXT ITEM AND FETCH THE ITEM ON THE SOF
 
 10 CALL exfort (srd,UNIT,jh,hdr,7,sp,0)
 20 NAME(1) = hdr(1)
 NAME(2) = hdr(2)
 item    = hdr(3)
 itest   = hdr(7)
 IF (itest == eoi) itest = eog
 IF (hdr(1) == dit .OR. hdr(1) == mdi) GO TO 200
 IF (hdr(3) == -1 .OR. hdr(1) == leof) GO TO 300
 itm = ittype(hdr(3))
 IF (itm == 1) GO TO 100
 rc = 3
 CALL sfetch (hdr(1),hdr(3),swrt,rc)
 IF (rc == 3) GO TO 60
 line = line + 2
 IF (line > nlpp) CALL page
 SELECT CASE ( rc )
   CASE (    1)
     GO TO 30
   CASE (    2)
     GO TO 30
   CASE (    3)
     GO TO 60
   CASE (    4)
     GO TO 40
   CASE (    5)
     GO TO 50
 END SELECT
 30 WRITE (nout,6346) uwm,hdr(1),hdr(2),hdr(3)
 usrmsg = .false.
 GO TO 60
 40 CALL exlvl (nos,z(idm),hdr,z,lcore)
 rc = 3
 CALL sfetch (hdr(1),hdr(3),swrt,rc)
 IF (rc == 3) GO TO 60
 50 CALL smsg (rc-2,hdr(3),hdr)
 usrmsg = .false.
 60 CONTINUE
 
!     TABLES
 
 
!     ELSETS TABLE CORRECTION BY G.CHAN/UNISYS   4/91
 
!     IN 91 VERSION, ELEMENT PLOT SYMBOL LINE HAS 2 WORDS, SYMBOL AND
!     NO. OF GRID POINT PER ELEMENT, NGPEL, WRITTEN OUT BY EXO2 USING
!     FORMAT 25. THE ELSETS DATA LINE COMING UP NEXT USE FORMAT 10 FOR
!     ELEMENTS WITH NO OFFSETS, FORMAT 26 FOR BAR WHICH HAS 6 OFFSET
!     VALUES, AND FORMATS 27 AND 28 FOR TRIA3 AND QUAD4 WHICH HAS 1
!     OFFSET VALUE EACH.
!     IN 90 AND EARLIER VERSIONS, ONLY ONE ELEMENT PLOT SYMBOL WORD WAS
!     WRITTEN OUT, AND ON ELSETS DATA LINE COMING UP NEXT, FORMAT 10
!     WAS USED FOR ALL ELEMENTS. NO OFFSET DATA FOR THE BAR, QUAD4 AND
!     TRIA3 ELEMENTS. NGPEL WAS THE FIRST WORD ON THE ELSETS DATA LINE.
!     ALSO, THE 90 AND EARLIER VERSIONS DID NOT COUNT PROPERTY ID, PID,
!     ON THE ELSETS DATA LINE. THUS THE TOTAL NO. OF WORDS MAY BE IN
!     ERROR AND MAY CAUSE EXTRA ZEROS AT THE END OF THE DATA LINE.
 
!     THEREFORE, IF THE 90 OR EARLIER EXTERNAL SOF FILE WAS USED, WE
!     NEED TO ADD THE OFFSETS (1 OR 6 FLOATING POINTS ZEROS) TO THE BAR,
!     QUAD4 AND TRIA3 ELEMENTS FOR THE ELSETS TABLE.
!     (AS OF 4/91, THESE CHANGES HAVE NOT BEEN TESTED)
 
 offset = 0
 70 nwds = hdr(5)
 IF (nwds > lcore) GO TO 9008
 CALL exfort (srd,UNIT,hdr(4),z,nwds,sp,0)
 IF (offset == 0) GO TO 80
 j = 1
 CALL suwrt (z(1),1,j)
 np2 = z(1) + 2
 DO  k = 2,nwds,np2
   IF (z(k) == 0) EXIT
   CALL suwrt (z(k),np2,j)
   CALL suwrt (zero,offset,j)
 END DO
 75 z(1) = 0
 nwds = 1
 80 CALL suwrt (z,nwds,itest)
 IF (hdr(7) == eoi) GO TO 90
 CALL exfort (srd,UNIT,jh,hdr,7,sp,0)
 IF (hdr(1) /= NAME(1) .OR. hdr(2) /= NAME(2)) GO TO 160
 IF (hdr(3) /= item) GO TO 160
 itest = hdr(7)
 IF (itest == eoi) itest = eog
 IF (item /= plts .OR. hdr(5) /= 1 .OR. hdr(4) /= 10) GO TO 85
 offset = 0
 IF (z(1) == bar) offset = 6
 IF (z(1) == q4 .OR. z(1) == t3) offset = 1
 85 IF (hdr(4) > 0) GO TO 70
 90 itest = eoi
 CALL suwrt (0,0,itest)
 GO TO 140
 
!     MATRICES
 
 
!     READ TRAILER
 
 100 CALL softrl (hdr(1),hdr(3),mcb(1))
 rc = mcb(1)
 IF (rc == 3) GO TO 108
 line = line + 2
 IF (line > nlpp) CALL page
 SELECT CASE ( rc )
   CASE (    1)
     GO TO 102
   CASE (    2)
     GO TO 102
   CASE (    3)
     GO TO 108
   CASE (    4)
     GO TO 104
   CASE (    5)
     GO TO 106
 END SELECT
 102 WRITE (nout,6346) uwm,hdr(1),hdr(2),hdr(3)
 usrmsg = .false.
 GO TO 108
 104 CALL exlvl (nos,z(idm),hdr,z,lcore)
 GO TO 108
 106 CALL smsg (3,hdr(3),hdr)
 usrmsg = .false.
 108 CALL exfort (srd,UNIT,hdr(4),mcb(2),6,sp,0)
 ncol   = mcb(2)
 prec   = mcb(5)
 mcb(1) = scr1
 mcb(2) = 0
 mcb(6) = 0
 mcb(7) = 0
 IF (usrmsg) CALL gopen (scr1,z(buf4),wrtrew)
 
!     READ MATRIX ONE COLUMN AT A TIME AND PACK ON SCR2
 
 DO  j = 1,ncol
   CALL exfort (srd,UNIT,jh,hdr,7,sp,0)
   IF (hdr(1) /= NAME(1) .OR. hdr(2) /= NAME(2)) GO TO 160
   IF (hdr(3) /= item) GO TO 160
   nwds = hdr(5)
   IF (nwds*1.4 > ncore) GO TO 9008
   CALL exfort (srd,UNIT,hdr(4),z,nwds,prec,dz)
   IF (.NOT. usrmsg) CYCLE
   CALL bldpk (prec,prec,scr1,0,0)
   iprc = prc(prec)
   n    = nword(prec) + iprc
   k    = 1
   110 IF (z(k) < 0) GO TO 120
   irow = z(k)
   a(1) = z(k+iprc)
   IF (prec == 1) GO TO 115
   a(2) = z(k+iprc+1)
   IF (prec <= 3) GO TO 115
   a(3) = z(k+4)
   a(4) = z(k+5)
   115 CALL zblpki
   k = k + n
   GO TO 110
   120 CALL bldpkn (scr1,0,mcb)
 END DO
 IF (.NOT.usrmsg) GO TO 150
 CALL wrttrl (mcb)
 CALL CLOSE  (scr1,rew)
 CALL mtrxo  (scr1,hdr,hdr(3),0,rc)
 
!     WRITE USER MESSAGE
 
 140 IF (.NOT.usrmsg) GO TO 150
 line = line + 1
 IF (line > nlpp) CALL page
 WRITE (nout,6357) uim,hdr(1),hdr(2),hdr(3),uname,sof
 150 usrmsg = .true.
 GO TO 10
 
!     NO EOI FOR ITEM AND A NEW ITEM WAS READ
 
 160 line = line + 2
 IF (line > nlpp) CALL page
 WRITE (nout,6363) uwm,NAME(1),NAME(2),item,uname
 IF (itm == 0) CALL DELETE (NAME,item,rc)
 IF (itm == 1) CALL CLOSE (scr1,rew)
 usrmsg = .true.
 GO TO 20
 
!     READ DIT AND MDI
 
 200 nos   = hdr(5)/2
 lcore = ncore - hdr(5)*4
 idm   = lcore + 1
 IF (6*nos > lcore) GO TO 9008
 CALL exfort (srd,UNIT,hdr(4),z,hdr(5),sp,0)
 DO  i = 1,nos
   z(idm+4*i-4) = z(2*i-1)
   z(idm+4*i-3) = z(2*i  )
 END DO
 CALL exfort (srd,UNIT,jh,hdr,7,sp,0)
 CALL exfort (srd,UNIT,hdr(4),z,hdr(5),sp,0)
 DO  i = 1,nos
   j = idm + 4*i - 2
   k = 6*i - 6
   z(j  ) = lshift(z(k+1),20) + lshift(z(k+2),10) + z(k+3)
   z(j+1) = lshift(z(k+4),20) + lshift(z(k+5),10) + z(k+6)
 END DO
 GO TO 10
 
!     NORMAL MODULE COMPLETION
 
 300 CALL sofcls
 RETURN
 
!     ABNORMAL MODULE COMPLETION
 
 9008 CALL mesage (8,0,subr)
 dry = -2
 CALL sofcls
 RETURN
 
!     MESSAGE TEXTS
 
 6346 FORMAT (a25,' 6346, SUBSTRUCTURE ',2A4,' ITEM ',a4,  &
     ' NOT COPIED.  IT ALREADY EXISTS ON THE SOF.')
 6357 FORMAT (a29,' 6357, SUBSTRUCTURE ',2A4,' ITEM ',a4,  &
     ' SUCCESSFULLY COPIED FROM ',2A4,' TO ',a4)
 6363 FORMAT (a25,' 6363, INCOMPLETE DATA FOR SUBSTRUCTURE ',2A4,  &
     ' ITEM ',a4,' ON ',2A4,'. THE ITEM WILL NOT BE COPIED.')
END SUBROUTINE exi2
