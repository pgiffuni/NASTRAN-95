SUBROUTINE xyplot
     
!     XYPLOT IS AN OUTPUT MODULE
 
!     INFORMATION SUPPLIED BY XYTRAN THROUGH DATA BLOCK XYPLOT
!     IS INTERPRETED AND OUTPUT TO EITHER PLT1(BCD TAPE FILE) OR
!     PLT2(BINARY TAPE FILE) FOR PLOTTING ON AN OFF-LINE PLOTTER.
 
 
 EXTERNAL        lshift,rshift
 INTEGER :: expo,isym(2),ix(1),lttn(10),lttp(10),d4,  &
     outape,rshift,sysbuf,xyplt
 REAL :: nums,tltv(22),x(1),y(1),xy(2),chrscl,cscale
 CHARACTER (LEN=27) :: swm
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim,sfm,swm
 COMMON /machin/ mach,ihalf
 COMMON /system/ ksystm(65)
!ZZ   COMMON /ZZXYPL/ Z(1)
 COMMON /zzzzzz/ z(20000)
 COMMON /xxparm/ ipltbf,icmra,ifskp,pnam1,pnam2,iptdn,npens,  &
     papszx,papszy,ptyp1,ptyp2,jpsz(8),pc(8,2), spare,ya(115)
 COMMON /pltdat/ model,ipltnr,xwmin,ylow,axmax,yup,xwmax,ywmax,  &
     xedge,yedge,xa(9),chrscl, xymax(2),cntpi,cch,ccv,all,mnp,apo(2),itp,ltape
 COMMON /xyplin/ idsb,nfrm,ncrv,idpe,ncom,idmj,itbf,  &
     nwfr,iskp,d1  ,xmin,xmax,ymin,ymax, xtic,xdtc,xltc,nxdg,ixpr,nxtt,ixvs,  &
     ixdt,ytic,ydtc,yltc,nydg,iypr,nytt, iyvs,iydt,ittc,ibtc,iltc,irtc,logx,  &
     logy,ixax,xint,iyax,yint,icrv,d2(2),  &
     ipens,ipenn,skp5(5),titl(32),sbtl(32),  &
     clbl(32),cvtl(32),xatl(32),yatl(32),  &
     ixgd,iygd,d3(37),cscale,ipsz,nplt,xpap, ypap,ncmr,d4(13)
 EQUIVALENCE     (ksystm( 1),sysbuf), (ksystm( 2),outape),  &
     (ksystm( 9),nlpp  ), (ksystm(12),nlines),  &
     (z(1),x(1),ix(1),xy(1)), (xy(2),y(1))
 DATA    lpltmd, lcmr, xlpap, ylpap / -1, -1, -1.0, -1.0   /
 DATA    xyplt / 101   /
 DATA    nrwd  , irdrw ,iclsrw / 300   , 0     ,1      /
 DATA    iplus , ie, lep, lem  / 1H+, 1HE, 4H1E+ , 4H1E-   /
 DATA    lttn  / 8,  8, 5, 4, 3, 2, 2, 1, 1, 1 / ,  &
     lttp  / 15,15,10, 6, 3, 1, 1, 7, 7, 7 / ,  &
     tltv  / 3.,6.,2.,5.,8.,2.,4.,6.,8.,2.,3.,5.,7.,9.,2.,3., 4.,5.,6.,7.,8.,9. /
 
 
!     DEFINITION OF COMMON BLOCK /PLTDAT/ CONTENTS
 
!     MODEL  - MODEL NUMBER OF THE CURRENT PLOTTER.
!     IPLTNR - NUMBER OF CURRENT PLOTTER IN USE
!     XWMIN  - MINIMUM X VALUE OF PLOTTING REGION IN PLOTTER COUNTS
!     YLOW   - MIN. Y VALUE OF PLOT. REGION(AFTER TITLES)
!              IN PLOTTER COUNTS
!     AXMAX  - MAX. X VALUE OF PLOT. REGION(LESS MARGIN)
!              IN PLOTTER COUNTS
!     YUP    - MAX. Y VALUE OF PLOT. REGION(LESS MARGIN)
!              IN PLOTTER COUNTS
!     XWMAX  - ACTUAL MAXIMUM REGION SIZE IN X DIRECTION
!              IN PLOTTER COUNTS
!     YWMAX  - ACTUAL MAXIMUM REGION SIZE IN Y DIRECTION
!              IN PLOTTER COUNTS
!     XEDGE  - MARGIN OF X EDGE IN PLOTTER COUNTS (TABLE PLOTTERS ONLY)
!     YEDGE  - MARGIN OF Y EDGE IN PLOTTER COUNTS (TABLE PLOTTERS ONLY)
!     XA     - SPARES
 
!     THE FOLLOWING SYMBOLIC VALUES PERTAIN TO THE CURRENT PLOTTER.
!     AND ARE SET WHEN STPLOT OR PLTSET IS CALLED.
 
!     XYMAX - X AND Y FRAME LIMITS IN PLOTTER COUNTS.
!     CNTPI - PLOTTER COUNTS PER INCH.
!     CCH   - HORIZONTAL PLOTTER COUNTS PER SINGLE CHARACTER
!     CCV   - VERTICAL PLOTTER COUNTS PER SINGLE CHARACTER
!     ALL   - MAXIMUM LINE LENGTH DRAWN WITH SINGLE COMMAND
!             (PLOTTER COUNT)
!     MNP   - MAXIMUM NUMBER OF PENS
!     APOX  - ACTUAL PLOTTER X ORIGIN IN PLOTTER COUNTS
!     APOY  - ACTUAL PLOTTER Y ORIGIN IN PLOTTER COUNTS
!             NOTE - INCREMENTAL PLOTTERS USE AS CURRENT PEN POSITION.
!     ITP   - PLOTTER TYPE.
!     LTAPE - GINO NAME OF THE PLOT TAPE.
 
!     DEFINITION OF I.D. RECORD CONTENTS OF INPUT DATA FILE /XYPLIN/
 
!     IDSB - SUBCASE I.D.               NFRM - FRAME NUMBER
!     NCRV - CURVE NUMBER               IDPE - POINT OR ELEMENT I.D.
!     NCOM - COMPONENT NUMBER           IDMJ - VECTOR NUMBER
!     ITBF - BOTTOM TOP FULL FRAME IND. NWFR - NEW AXIS AND LABEL IND.
!     ISKP - FRAME SKIP NUMBER          D1   - SPARE
!     XMIN - MINIMUM X DATA FOR CURVE   XMAX - MAXIMUM X DATA FOR CURVE
!     YMIN - MINIMUM Y DATA FOR CURVE   YMAX - MAXIMUM Y DATA FOR CURVE
!     XTIC - FIRST X TICK VALUE         XDTC - VALUE BETWEEN X TICKS
!     XLTC - HIGHEST X-VALUE ON FRAME.  NXDG - MAX. DIGITS FOR X-TICKS
!     IXPR - 10 POWER ON PRINTED X TICK NXTT - TOTAL NUMBER OF X TICKS
!     IXVS - X TICKS BETWEEN LABELS     IXDT - DELTA PRINT VALUE X TICKS
!     YTIC - FIRST Y TICK VALUE         YDTC - VALUE BETWEEN Y TICKS
!     YLTC - HIGHEST Y-VALUE ON FRAME.  NYDG - MAX. DIGITS FOR Y-TICKS
!     IYPR - 10 POWER ON PRINTED Y TICK NYTT - TOTAL NUMBER OF Y TICKS
!     IYVS - Y TICKS BETWEEN LABELS     IYDT - DELTA PRINT VALUE Y TICKS
!     ITTC - TICKS W/WO VALUES - TOP    IBTC - TICKS W/WO VALUES - BOTTM
!     ILTC - TICKS W/WO VALUES - LEFT   IRTC - TICKS W/WO VALUES - RIGHT
!     LOGX - LINEAR/LOG - X DIRECTION   LOGY - LINEAR/LOG - Y DIRECTION
!     IXAX - X AXIS/NO AXIS INDICATOR   XINT - X AXIS  Y INTERCEPT
!     IYAX - Y AXIS/NO AXIS INDICATOR   YINT - Y AXIS  X INTERCEPT
!     ICRV - POINT/LINE PLOT INDICATOR  D2   - SPARES
!     TITL - PLOT TITLE                 SBTL - PLOT SUBTITLE
!     CLBL - PLOT LABEL                 CVTL - PLOT CURVE TITLE
!     XATL - X AXIS TITLE               YATL - Y AXIS TITLE
!     IXGD - X GRID LINES               IYGD - Y GRID LINES
!     D3   - SPARES                     IPNR - PEN COLOR
!     IPSZ - PEN SIZE                   NPLT - TYPE OF PLOTTER
!     XPAP - PAPER SIZE(IN.) X DIR.     YPAP - PAPER SIZE(IN.) Y DIR.
!     NCMR - CAMERA NR. FOR SC-4020     D4   - XYTRAN INTERNAL FLAGS
 
 
!     SET IOPN=0 (PLOT TAPE CLOSED) AND NERR=0 (NUMBER OF ID RECORDS
!     WITH WRONG WORD COUNT).  WHEN NERR=5, XYPLOT ASSUMES BAD INPUT
!     FILE AND ABANDONS OPERATION.
 
 mb1  = korsz(z) - sysbuf
 ipchg= 0
 iopn = 0
 CALL OPEN (*920,xyplt,z(mb1),irdrw)
 99 CALL fwdrec (*960,xyplt)
 nerr = 0
 
!     READ I.D. RECORD ON INPUT DATA FILE
 
 100 CALL READ (*960,*120,xyplt,idsb,nrwd+1,1,nact)
 110 nerr = nerr + 1
 IF (nerr >= 5) GO TO 940
 GO TO 100
 120 IF (nact /= nrwd) GO TO 110
 
!     SKIP DATA IF IT WAS FOR THE PAPERPLOTER ONLY
 
 IF (d4(2) <= 0) GO TO 99
 IF (nwfr  /= 0) GO TO 270
 
!     READ DATA PAIRS FROM INPUT DATA FILE FOR CURVE TO BE PLOTTED
 
 130 CALL READ (*960,*250,xyplt,z,mb3,0,nact)
 
!     SET IFIN TO SHOW MORE DATA REMAINING TO BE READ FROM RECORD.
!     SET L AS INDEX TO LAST LEGITIMATE X VALUE OF DATA PAIRS IN CORE.
 
 ifin = 0
 l    = mb3 - 1
 140 IF (ix(l) /= 1) GO TO 150
 l = l - 2
 IF (l <= 0) GO TO 240
 
!     CONVERT DATA POINTS TO PLOTTER COUNTS AND PLOT SYMBOL AT EACH
!     LEGITIMATE POINT WHEN REQUIRED.
 
 150 IF (icrv /= 0) CALL symbol (0,0,0,-1)
 
 isym(1) = IABS(icrv) + ncrv - 1
 isym(2) = 0
 
 DO  i = 1,l,2
   IF (ix(i) ==    1) CYCLE
   IF (x(i) > xmaxs) GO TO 180
   IF (x(i) < xmins) GO TO 180
   IF (logxs <=    0) GO TO 160
   x(i) = ALOG10(x(i))
   160 x(i) = xdr*x(i)+xc
   IF (y(i) > ymaxs) GO TO 180
   IF (y(i) < ymins) GO TO 180
   IF (logys <=    0) GO TO 170
   y(i) = ALOG10(y(i))
   170 y(i) = ydr*y(i) + yc
   IF (icrv /= 0) CALL symbol (x(i),y(i),isym,0)
   CYCLE
   180 ix(i  ) = 1
   ix(i+1) = 1
 END DO
 IF (icrv /= 0) CALL symbol (0,0,0,1)
 
!     PLOT LINES BETWEEN LEGITIMATE POINTS WHEN REQUIRED
 
 IF (icrv < 0 .AND. ipenn > 0) icrv = -icrv
 IF (icrv < 0) GO TO 240
 CALL line (0,0,0,0,0,-1)
 oldx = x(1)
 oldy = y(1)
 IF (ipchg == 1) GO TO 193
 icpen = ipsz
 IF (ipens == 0) GO TO 192
 icpen = ipens
 ipchg = 1
 GO TO 192
 193 IF (icpen == ipenn) icpen = ipens - 1
 icpen = icpen + 1
 192 CONTINUE
 DO  i = 1,l,2
   IF (ix(i) == 1) GO TO 220
   t1 = oldx - x(i)
   t2 = oldy - y(i)
   IF (t1 == 0.0) THEN
     GO TO   200
   ELSE
     GO TO   210
   END IF
   200 IF (t2 == 0.0) THEN
     GO TO   230
   END IF
   210 CALL line (oldx,oldy,x(i),y(i),icpen,0)
   oldx = x(i)
   oldy = y(i)
   CYCLE
   220 oldx = x(i+2)
   oldy = y(i+2)
 END DO
 CALL line (0,0,0,0,0,1)
 240 IF (ifin == 0) THEN
   GO TO   130
 ELSE
   GO TO   100
 END IF
 
!     ALL DATA PAIRS IN CORE, SET IFIN TO SHOW NO MORE DATA REMAINS
!     FOR PRESENT CURVE.  IF ODD NUMBER OF DATA VALUES OUTPUT WARNING
!     MESSAGE AND CONTINUE.  SET L AS INDEX TO LAST X VALUE OF DATA
!     PAIRS.
 
 250 ifin = 1
 IF (nact == (nact/2)*2) GO TO 260
 nact = nact - 1
 WRITE (outape,990) uwm,nfrm,ncrv
 nlines = nlines + 2
 IF (nlines >= nlpp) CALL page
 260 l = nact - 1
 IF (l > 0) THEN
   GO TO   140
 ELSE
   GO TO   100
 END IF
 
!     NEW AXIS, LABELS, ETC. ARE NEEDED.
 
!     NASTRAN PLOTTING SOFTWARE INITIALIZATION.
 
 270 IF (itbf >= 0 .AND. iopn /= 0) CALL stplot (-1)
 ipltnr = rshift(nplt,ihalf)
 model  = nplt - lshift(ipltnr,ihalf) - 100
 IF (ncmr > 0) icmra=ncmr
 ifskp  = iskp
 cscale = chrscl
 IF (cscale < 1.) cscale = 1.0
 IF (xpap   > 0.) papszx = xpap
 IF (ypap   > 0.) papszy = ypap
 DO  i = 1,npens
   jpsz(i) = ipsz
 END DO
 IF (itbf >= 0) GO TO 284
 
!     LOWER HALF MAY NOT CHANGE FRAME OR PLOTTER OR CALL PLTSET
 
!     IF (NCMR .NE. LCMR  ) GO TO 925
 IF (xpap /= xlpap ) GO TO 925
 IF (ypap /= ylpap ) GO TO 925
 IF (nplt /= lpltmd) GO TO 925
 GO TO 286
 
 284 CALL pltset
 lcmr  = ncmr
 xlpap = xpap
 ylpap = ypap
 lpltmd= nplt
 mb2   = mb1 - ipltbf
 mb3   = 2*((mb2-1)/2)
 
!     SET VALUES FOR FULL FRAME PLOTTING
 
 286 ywmin= 0.
 ylow = 4.*ccv
 yxtr = (ywmax+ylow)/2.
 
!     START A NEW PLOT IF NECESSARY.
 
 IF (itbf < 0) GO TO 290
 CALL sopen (*930,ltape,z(mb2),ipltbf)
 iopn = 1
 CALL stplot (nfrm)
 290 CALL PRINT (0,0,0,0,0,-1)
 IF (itbf < 0) THEN
   GO TO   300
 ELSE IF (itbf == 0) THEN
   GO TO   320
 ELSE
   GO TO   310
 END IF
 
!     MODIFY VALUE FOR LOWER HALF FRAME PLOTTING
 
 300 yup = yxtr
 GO TO 330
 
!     MODIFY VALUE FOR UPPER HALF FRAME PLOTTING
 
 310 ylow = yxtr
 
!     SAVE YLOW AND EXPAND REGION SIZE FOR PRINTING OF TITLES.  RESTORE
!     YLOW AFTER PRINTING THE FOUR CURVE TITLES AT BOTTOM OF FRAME.
 
 320 xprm = xwmin
 yprm = ywmin
 y1t  = ylow
 ylow = ywmin
 CALL PRINT (xprm,yprm,1,clbl(1),32,0)
 yprm = yprm + ccv
 CALL PRINT (xprm,yprm,1,sbtl(1),32,0)
 yprm = yprm + ccv
 CALL PRINT (xprm,yprm,1,titl(1),32,0)
 yprm = yprm + ccv
 CALL PRINT (xprm,yprm,1,cvtl(1),32,0)
 ylow = y1t
 
!     OUTPUT X AND Y AXES TITLES
 
 330 yprm = ylow
 xprm = xwmin + 8.*cch
 CALL PRINT (xprm,yprm,1,xatl(1),32,0)
 yprm = yup - 2*ccv
 xprm = xwmin
 CALL PRINT (xprm,yprm,2,yatl(1),32,0)
 CALL tipe (0,0,0,0,0,1)
 
!     MEANING OF SYMBOLS USED
!     XDR,XC,YDR,YC - FACTORS TO CONVERT ENGINEERING UNITS TO PLOTTER
!                     COUNTS IN X AND Y DIRECTIONS.
!     CONVERSION IS - PLOTTER COUNTS = ENG. UNITS * XDR  +  XC
 
!     JTC,J1T,J2T,J3T,J4T,J5T - TEMPORARY INTEGER VALUES
!     T1,T2,T3,T4,X1T,Y1T     - TEMPORARY REAL VALUES
 
!     TEST XMAX,XMIN,YMAX, AND YMIN FOR COMPATIBILITY
 
 n  = 0
 340 dx = xltc - xtic
 dy = yltc - ytic
 IF (dx > 0.0 .AND. dy > 0.0) GO TO 440
 IF (n /=  0) GO TO 350
 IF (dx <= 0.0) xltc = xtic + xdtc*FLOAT(nxtt+1)
 IF (dy <= 0.0) yltc = ytic + ydtc*FLOAT(nytt+1)
 n = 1
 GO TO 430
 350 n = 2
 IF (dx > 0.0) GO TO 360
 xltc = xtic + 10.0
 xdtc = 2.0
 nxtt = 0
 360 IF (dy > 0.0) GO TO 430
 yltc = ytic + 10.0
 ydtc = 2.0
 nytt = 4
 
!     PRINT WARNING (N=NO. OF PASSES TO CORRECT)
 
 430 WRITE (outape,1010) uwm,n,nfrm
 nlines = nlines + 2
 IF (nlines >= nlpp) CALL page
 IF (n == 1) GO TO 340
 
!     SAVE XMAX, XMIN, YMAX, YMIN, LOGX AND LOGY FOR USE IF NEXT
!     I.D. RECORD IS NOT A NEW FRAME
 
 440 logxs = logx
 logys = logy
 xmins = xtic
 xmaxs = xltc
 ymins = ytic
 ymaxs = yltc
 
!     CALCULATE CONVERSION FACTORS
 
 xpl = xwmax - 7.*cch
 xps = xwmin + 8.*cch
 ypl = yup   - 2.*ccv
 yps = ylow  + 2.*ccv
 
!     PUT FRAME AT X AND Y MAXIMUM AND MINIMUM LIMITS
 
 IF (ixgd == 0 .AND. iygd == 0) GO TO 450
 CALL axis (0,0,0,0,0,-1)
 CALL axis (xps,yps,xps,ypl,ipsz,0)
 CALL axis (xps,ypl,xpl,ypl,ipsz,0)
 CALL axis (xpl,ypl,xpl,yps,ipsz,0)
 CALL axis (xpl,yps,xps,yps,ipsz,0)
 CALL axis (0,0,0,0,0,+1)
 450 IF (logx <= 0) GO TO 460
 xtic = ALOG10(xtic)
 xltc = ALOG10(xltc)
 dx   = xltc - xtic
 IF (iyax == 1) yint = ALOG10(yint)
 460 IF (logy <= 0) GO TO 470
 ytic = ALOG10(ytic)
 yltc = ALOG10(yltc)
 dy   = yltc - ytic
 IF (ixax == 1) xint = ALOG10(xint)
 470 xdr = (xpl-xps)/dx
 xc  = (xps*xltc-xpl*xtic)/dx
 ydr = (ypl-yps)/dy
 yc  = (yps*yltc-ypl*ytic)/dy
 
!     PREPARE TO CREATE + LABEL ANY REQUESTED TIC MARKS IN THE
!     X-DIRECTION.
 
 IF (ittc == 0 .AND. ixax /= 1 .AND. ibtc == 0 .AND. iygd == 0) GO TO 575
 ndg = 0
 IF (logx > 0) GO TO 480
 dtc = xdtc
 IF (dtc > 0. .AND. nxtt > 0) GO TO 477
 ntt = 0
 GO TO 485
 477 ntt = nxtt
 xts = xtic*xdr + xc
 IF (ittc <= 0 .AND. ibtc <= 0) GO TO 485
 ndg  = MIN0(nxdg+1,6)
 expo = ndg + ixpr - 2
 nums = xtic/10.**expo
 dl   = dtc/10.**expo
 lstep= MAX0(ixvs+1,1)
 GO TO 485
 480 ntt = logx + 1
 xts = xtic*xdr + xc
 dtc = 1.
 ndg = 4
 IF (logx > 10) GO TO 485
 ill  = lttp(logx)
 nitk = lttn(logx)
 
 485 DO  k = 1,3
   label = 1
!WKBR 9/93 LOG = XTIC - 1.0 + SIGN(0.1,XTIC)
   LOG = xtic - 1.0 + SIGN(0.1,xtic) - 1
   SELECT CASE ( k )
     CASE (    1)
       GO TO 490
     CASE (    2)
       GO TO 495
     CASE (    3)
       GO TO 500
   END SELECT
   
!     TICS + LABELS AT THE TOP.
   
   490 itc= ittc
   yt = ypl
   yl = yt + ccv
   GO TO 505
   
!     TICS ALONG THE X-AXIS.
   
   495 itc = 0
   IF (ixax == 1) itc = -1
   IF (itc  == 0) GO TO 505
   yt = xint*ydr + yc
   CALL axis (0,0,0,0,0,-1)
   CALL axis (xpl,yt,xps,yt,ipsz,0)
   GO TO 505
   
!     TICS + LABELS AT THE BOTTOM.
   
   500 itc= ibtc
   yt = yps
   yl = yt - ccv
   
   505 IF (itc == 0 .OR. ntt <= 0) CYCLE
   CALL tipe (0,0,0,0,0,-1)
   DO  j = 1,ntt
     r = xts + dtc*xdr*FLOAT(j-1)
     CALL tipe (r,yt,1,iplus,1,0)
     IF (logx > 0) GO TO 530
     IF (itc < 0 .OR. label /= j) CYCLE
     
!     LABEL THIS LINEAR TIC MARK.
     
     ifield = ndg
     rnum   = nums + dl*FLOAT(j-1)
     IF (rnum < 0.0) THEN
       GO TO   510
     ELSE IF (rnum == 0.0) THEN
       GO TO   525
     ELSE
       GO TO   515
     END IF
     510 ifield = ifield + 1
     515 t = ABS(rnum)
     IF (t >= 1.e-4) GO TO 525
     IF (t >= 5.e-5) GO TO 520
     rnum = 0.
     GO TO 525
     520 rnum = SIGN(1.e-4,rnum)
     525 CALL typflt (r,yl,1,rnum,ifield,0)
     label = label + lstep
     IF (label <= ntt) CYCLE
     r = r + FLOAT(ifield)*cch
     CALL tipe (r,yl,1,ie,1,0)
     CALL typint(r+cch,yl,1,expo,0,0)
     CYCLE
     
!     LABEL THIS LOGARITHMIC CYCLE TIC MARK.
     
     530 LOG = LOG + 1
     IF (itc < 0) GO TO 535
     i = lep
     IF (LOG < 0) i = lem
     CALL PRINT (r-cch,yl,1,i,1,0)
     CALL typint (r+2.*cch,yl,1,IABS(LOG),0,0)
     535 IF (logx > 10 .OR. j == ntt) CYCLE
     
!     CREATE + LABEL THE LOGARITHMIC INTRACYCLE TIC MARKS WITHIN THIS
!     CYCLE.
     
     DO  i = 1,nitk
       l = ill + i - 1
       t = xdr*(ALOG10(tltv(l))+FLOAT(LOG)) + xc
       CALL tipe (t,yt,1,iplus,1,0)
       IF (itc < 0) CYCLE
       l = tltv(l) + .01
       CALL typint (t,yl,1,l,1,0)
     END DO
     
   END DO
   CALL tipe (0,0,0,0,0,+1)
 END DO
 IF (iygd == 0 .OR. ntt <= 0) GO TO 575
 
!     DRAW THE Y-DIRECTION GRID NETWORK.
 
 CALL axis (0,0,0,0,0,-1)
!WKBR 9/93 LOG = XTIC - 1.0 + SIGN(0.1,XTIC)
 LOG = xtic - 1.0 + SIGN(0.1,xtic) - 1
 k = 1
 DO  j = 1,ntt
   k = -k
   r = xts + dtc*xdr*FLOAT(ntt-j)
   IF (k > 0) CALL axis (r,ypl,r,yps,ipsz,0)
   IF (k < 0) CALL axis (r,yps,r,ypl,ipsz,0)
   IF (logx <= 0 .OR. logx > 10 .OR. j == ntt) CYCLE
   
!     DRAW THE Y-DIRECTION GRID LINES WITHIN THIS LOGARITHMIC CYCLE.
   
   LOG = LOG + 1
   DO  i = 1,nitk
     l = ill + nitk - i
     t = xdr*(ALOG10(tltv(l))+FLOAT(LOG)) + xc
     k = -k
     IF (k > 0) CALL axis (t,ypl,t,yps,ipsz,0)
     IF (k < 0) CALL axis (t,yps,t,ypl,ipsz,0)
   END DO
   
 END DO
 CALL axis (0,0,0,0,0,+1)
 
!     PREPARE TO CREATE + LABEL ANY REQUESTED TIC MARKS IN THE
!     Y-DIRECTION.
 
 575 IF (iltc == 0 .AND. iyax /= 1 .AND. irtc == 0 .AND. ixgd == 0) GO TO 130
 ndg = 0
 IF (logy > 0) GO TO 580
 dtc = ydtc
 IF (dtc > 0. .AND. nytt > 0) GO TO 577
 ntt = 0
 GO TO 585
 577 ntt = nytt
 yts = ytic*ydr + yc
 IF (iltc <= 0 .AND. irtc <= 0) GO TO 585
 ndg  = MIN0(nydg+1,6)
 expo = ndg + iypr - 2
 nums = ytic/10.**expo
 dl   = dtc/10.**expo
 lstep= MAX0(iyvs+1,1)
 GO TO 585
 580 ntt = logy + 1
 yts = ytic*ydr + yc
 dtc = 1.
 ndg = 4
 IF (logy > 10) GO TO 585
 ill = lttp(logy)
 nitk= lttn(logy)
 
 585 DO  k = 1,3
   label = 1
   LOG = ytic - 1.0 + SIGN(0.1,ytic)
   SELECT CASE ( k )
     CASE (    1)
       GO TO 590
     CASE (    2)
       GO TO 595
     CASE (    3)
       GO TO 600
   END SELECT
   
!     TICS + LABELS ON THE LEFT SIDE.
   
   590 itc= iltc
   xt = xps
   xl = xt - cch*FLOAT(ndg+1)
   GO TO 605
   
!     TICS ALONG THE Y-AXIS.
   
   595 itc = 0
   IF (iyax == 1) itc = -1
   IF (itc  == 0) GO TO 605
   xt = yint*xdr + xc
   CALL axis (0,0,0,0,0,-1)
   CALL axis (xt,ypl,xt,yps,ipsz,0)
   GO TO 605
   
!     TICS + LABELS ON THE RIGHT SIDE.
   
   600 itc= irtc
   xt = xpl
   xl = xt + cch
   
   605 IF (itc == 0 .OR. ntt <= 0) CYCLE
   CALL tipe (0,0,0,0,0,-1)
   DO  j = 1,ntt
     s = yts + dtc*ydr*FLOAT(j-1)
     CALL tipe (xt,s,1,iplus,1,0)
     IF (logy > 0) GO TO 630
     IF (itc < 0 .OR. label /= j) CYCLE
     
!     LABEL THIS LINEAR TIC MARK.
     
     ifield = ndg
     rnum   = nums + dl*FLOAT(j-1)
     IF (rnum < 0.0) THEN
       GO TO   610
     ELSE IF (rnum == 0.0) THEN
       GO TO   625
     ELSE
       GO TO   615
     END IF
     610 ifield = ifield + 1
     615 t = ABS(rnum)
     IF (t >= 1.e-4) GO TO 625
     IF (t >= 1.e-5) GO TO 620
     rnum = 0.
     GO TO 625
     620 rnum = SIGN(1.e-4,rnum)
     625 CALL typflt (xl,s,1,rnum,-ifield,0)
     label  = label + lstep
     ylabel = s
     CYCLE
     
!     LABEL THIS LOGARITHMIC CYCLE TIC MARK.
     
     630 LOG = LOG + 1
     IF (itc < 0) GO TO 635
     i = lep
     IF (LOG < 0) i = lem
     CALL PRINT (xl,s,1,i,1,0)
     CALL typint (xl+3.*cch,s,1,IABS(LOG),0,0)
     635 IF (logy > 10 .OR. j == ntt) CYCLE
     
!     CREATE + LABEL THE LOGARITHMIC INTRACYCLE TIC MARKS WITHIN THIS
!     CYCLE.
     
     DO  i = 1,nitk
       l = ill + i - 1
       t = ydr*(ALOG10(tltv(l))+FLOAT(LOG)) + yc
       CALL tipe (xt,t,1,iplus,1,0)
       IF (itc < 0) CYCLE
       l = tltv(l) + .01
       CALL typint (xl,t,1,l,1,0)
     END DO
     
   END DO
   IF (itc < 0 .OR. logy > 0) GO TO 650
   CALL tipe (xl,ylabel-ccv,1,ie,1,0)
   CALL typint (xl+cch,ylabel-ccv,1,expo,0,0)
   650 CALL tipe (0,0,0,0,0,+1)
 END DO
 IF (ixgd == 0 .OR. ntt <= 0) GO TO 130
 
!     DRAW THE X-DIRECTION GRID NETWORK.
 
 CALL axis (0,0,0,0,0,-1)
 LOG = ytic - 1.0 + SIGN(0.1,ytic)
 k = 1
 DO  j = 1,ntt
   k = -k
   s = yts + dtc*ydr*FLOAT(ntt-j)
   IF (k > 0) CALL axis (xps,s,xpl,s,ipsz,0)
   IF (k < 0) CALL axis (xpl,s,xps,s,ipsz,0)
   IF (logy <= 0 .OR. logy > 10 .OR. j == ntt) CYCLE
   
!     DRAW THE X-DIRECTION GRID LINES WITHIN THIS LOGARITHMIC CYCLE...
   
   LOG = LOG + 1
   DO  i = 1,nitk
     l = ill + nitk - i
     t = ydr*(ALOG10(tltv(l))+FLOAT(LOG)) + yc
     k = -k
     IF (k > 0) CALL axis (xps,t,xpl,t,ipsz,0)
     IF (k < 0) CALL axis (xpl,t,xps,t,ipsz,0)
   END DO
   
 END DO
 CALL axis (0,0,0,0,0,+1)
 GO TO 130
 
!     OUTPUT WARNING NESSAGES, CLOSE INPUT FILE AND PLOT TAPE AND RETURN
 
 920 RETURN
 925 WRITE (outape,1020) swm
 GO TO 950
 930 WRITE (outape,1000) uwm,ltape
 GO TO 950
 940 WRITE (outape,980) uwm
 950 nlines = nlines + 2
 IF (nlines >= nlpp) CALL page
 960 CALL CLOSE (xyplt,iclsrw)
 IF (iopn /= 0) CALL stplot (-1)
 RETURN
 
 980 FORMAT (a25,' 992, XYPLOT INPUT DATA FILE ID. RECORDS TOO SHORT.',  &
     '  XYPLOT ABANDONED.')
 990 FORMAT (a25,' 993, XYPLOT FOUND ODD NR. OF VALUES FOR DATA PAIRS',  &
     ' IN FRAME',i5,', CURVE NR.',i5,'.  LAST VALUE IGNORED.')
 1000 FORMAT (a25,' 994, XYPLOT OUTPUT FILE NAME ',a4,' NOT FOUND.',  &
     '  XYPLOT ABANDONED.')
 1010 FORMAT (a25,' 997, NR.',i4,'.  FRAME NR.',i5,' INPUT DATA ',  &
     'INCOMPATIBLE.  ASSUMPTIONS MAY PRODUCE INVALID PLOT.')
 1020 FORMAT (a27,' 998, XYPLOT PLOTTER OR FRAME MAY NOT CHANGE FOR ',  &
     'LOWER FRAME.  XYPLOT ABANDONED.')
END SUBROUTINE xyplot
