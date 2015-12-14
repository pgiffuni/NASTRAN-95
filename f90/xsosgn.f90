SUBROUTINE xsosgn
     
!     THIS SUBROUTINE SCANS THE OSCAR TAPE AND GENERATES THE SOS + MD
 
!     LAST REVISED BY G.CHAN/UNISYS TO REMOVE THE VAX AND NOT-VAX
!     LOGICS, AND TO SYNCHRONIZE THE SCRATH FILE NAMES AS SET FORTH BY
!     THE XSEMX ROUTINES.   2/1990
 
 IMPLICIT INTEGER (a-z)
!     LOGICAL         DEC
 EXTERNAL        andf,orf,lshift,rshift
 DIMENSION       block1(93),str(30),nsosgn(2),fequ(1),fntu(1),  &
     fon(1),ford(1),minp(1),mlsn(1),mout(1),mscr(1),  &
     sal(1),sdbn(1),sntu(1),sord(1),BLOCK(100), numbr(10)
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim,sfm
 COMMON /system/ ibufsz,outtap
 COMMON /xfiat / fiat(1),fmxlg,fculg,FILE(1),fdbn(2),fmat(1)
 COMMON /xfist / fist
 COMMON /xdpl  / dpd(1),dmxlg,dculg,ddbn(2),dfnu(1)
 COMMON /zzzzzz/ buf1(1)
 COMMON /xsfa1 / md(401),sos(1501),comm(20),xf1at(1),fpun(1),  &
     fcum(1),fcus(1),fknd(1)
 COMMON /isosgn/ entn5,entn6,k,j,str
 EQUIVALENCE     (dpd(1),dnaf),(fiat(1),funlg),(FILE(1),fequ(1)),  &
     (FILE(1),ford(1)),(BLOCK(8),block1(1))
 EQUIVALENCE     (md(1),mlgn),(md(2),mlsn(1)),(md(3),minp(1)),  &
     (md(4),mout(1)),(md(5),mscr(1)), (sos(1),slgn),(sos(2),sdbn(1)),  &
     (sos(4),sal(1),sntu(1),sord(1)),  &
     (comm(1),almsk),(comm(2),apndmk),(comm(3),cursno),  &
     (comm(4),entn1),(comm(5),entn2 ),(comm (6),entn3),  &
     (comm(7),entn4),(comm(8),flag  ),(comm (9),fnx  ),  &
     (comm(10),lmsk),(comm(11),lxmsk),(comm(13),rmsk ),  &
     (comm(14),rxmsk),(comm(15),s  ),(comm(16),scornt),  &
     (comm(17),tapmsk),(comm(19),zap), (xf1at(1),fntu(1),fon(1))
 DATA    jump  / 4HJUMP/, rept /4HREPT/,  cond/4HCOND/
 DATA    oscar / 4HPOOL/, scrn1,scrn2 / 4HSCRA,4HTCH0/
 DATA    nsosgn/ 4HXSOS , 2HGN /
 DATA    numbr / 1H1,1H2,1H3,1H4,1H5,1H6,1H7,1H8,1H9,1H0  /
 
 iflag = 0
 CALL OPEN (*500,oscar,buf1,2)
 CALL bckrec (oscar)
 CALL READ (*400,*600,oscar,BLOCK,7,0,flag)
 IF (BLOCK(2) /= cursno) GO TO 900
 GO TO 103
 
!     READ OSCAR FORMAT HEADER + 1
 
 100 IF (j > 1400 .OR. k > 390) GO TO 410
 CALL READ (*400,*600,oscar,BLOCK,7,0,flag)
 103 BLOCK(3) = andf(rmsk,BLOCK(3))
 IF (BLOCK(6) >= 0) GO TO 108
 IF (BLOCK(3) <= 2) GO TO 110
 IF (BLOCK(3) /= 3) GO TO 108
 l = rshift(andf(lxmsk,BLOCK(7)),16) - BLOCK(2)
 IF (BLOCK(4) /= jump) GO TO 106
 IF (l <= 1) GO TO 107
 DO  i = 1,l
   CALL fwdrec (*400,oscar)
 END DO
 GO TO 100
 106 IF (BLOCK(4) /= rept .AND. BLOCK(4) /= cond) GO TO 108
 107 IF (l < 0) iflag = -1
 108 CALL fwdrec (*400,oscar)
 GO TO 100
 
!     INPUT FILES
 
 110 minp(k) = BLOCK(7)
 IF (BLOCK(7) == 0) GO TO 300
 nwds= BLOCK(7)*entn5
 ASSIGN 150 TO isw
 
!     FILES READER
 
 130 CALL READ (*400,*600,oscar,block1,nwds+1,0,flag)
 blkcnt = 0
 DO  i = 1,nwds,entn5
   IF (block1(i) ==  0) GO TO 140
   sos(j  ) = block1(i  )
   sos(j+1) = block1(i+1)
   sos(j+2) = block1(i+2)
   j = j+3
   IF (j > 1500) GO TO 460
   CYCLE
   140 blkcnt = blkcnt + 1
 END DO
 GO TO isw, (150,170)
 
 150 minp(k) = minp(k) - blkcnt
 IF (BLOCK(3) == 2) GO TO 310
 
!     OUTPUT FILES
 
 mout(k) = block1(nwds+1)
 155 IF (mout(k) == 0) GO TO 320
 nwds = mout(k)*entn6
 ASSIGN 170 TO isw
 GO TO 130
 
 170 mout(k) = mout(k) - blkcnt
 175 CALL fwdrec (*400,oscar)
 
!     SCRATCH FILES
 
 mscr(k) = block1(nwds+1)
 IF (mscr(k) == 0) GO TO 230
 l = mscr(k)
 scrn3  = scrn2
 lll = 1
 ll  = 0
 DO  i = 1,l
   ll  = ll + 1
   IF (ll == 10) scrn3 = khrfn1(scrn3,3,numbr(lll),1)
   sos(j  ) = scrn1
   sos(j+1) = khrfn1(scrn3,4,numbr(ll),1)
   IF (ll /= 10) GO TO 200
   ll  = 0
   lll = lll + 1
   200 IF (str(i) == 0) GO TO 210
   n1= str(i)
   sos(n1) = orf(lmsk,BLOCK(2))
   210 str(i)  = j + 2
   sos(j+2)= scornt + i
   j = j + 3
   IF (j > 1500) GO TO 460
 END DO
 
 230 mlsn(k) = BLOCK(2)
 IF (iflag == 0)  GO TO 240
 mlsn(k) = orf(s,mlsn(k))
 240 IF (minp(k)+mout(k)+mscr(k) == 0) GO TO 100
 k= k + entn3
 IF (k > 400) GO TO 460
 GO TO 100
 
!     ZERO INPUT FILES
 
 300 CALL READ (*400,*600,oscar,BLOCK(7),1,0,flag)
 IF (BLOCK(3) == 2) GO TO 310
 mout(k) = BLOCK(7)
 GO TO 155
 
!     TYPE O FORMAT - NO OUTPUTS
 
 310 mout(k) = 0
 GO TO 175
 
!     ZERO OUTPUT FILES
 
 320 CALL READ (*400,*600,oscar,block1(nwds+1),1,0,flag)
 GO TO 175
 
 400 CALL skpfil (oscar,-1)
 410 CALL CLOSE  (oscar, 2)
 slgn = (j-1)/entn2
 mlgn = (k-1)/entn3
 RETURN
 
!     SYSTEM FATAL MESSAGES
 
 460 WRITE  (outtap,461) sfm
 461 FORMAT (a25,' 1011, MD OR SOS TABLE OVERFLOW')
 GO TO  1000
 500 WRITE  (outtap,501) sfm
 501 FORMAT (a25,' 1012, POOL COULD NOT BE OPENED')
 GO TO  1000
 600 WRITE  (outtap,601) sfm
 601 FORMAT (a25,' 1013, ILLEGAL EOR ON POOL')
 GO TO  1000
 900 WRITE  (outtap,901) sfm,BLOCK(2),cursno
 901 FORMAT (a25,' 1014, POOL FILE MIS-POSITIONED ',2I7)
 1000 CALL mesage (-37,0,nsosgn)
 RETURN
END SUBROUTINE xsosgn
