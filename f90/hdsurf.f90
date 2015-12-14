SUBROUTINE hdsurf (gplst,x,u,pen,deform,nmax,maxsf,iz,ib,pedge,  &
        iopcor)
     
!     THIS ROUTINE PREPARES THE ELEMENT SURFACES FOR HIDDEN LINE PLOT
!     IT ALSO GENERATES THE SHRINK PLOT IF SHRINK ALONE IS REQUESTED.
!     IF SHRINK AND HIDDEN ARE REQUESTED, THIS ROUTINE WILL PREPARE THE
!     SHRUNK SURFACES FOR HDPLOT.
 
!     REVISED  10/1990 BY G.CHAN/UNISYS
!     (1) HIDDEN PLOT WITH SOLID ELEMENTS BUGS
!     (2) HIDDEN AND SHRINK TOGETHER
!     (3) SKIP ANY OFFSET DATA IN ELSET FILE IF THEY ARE PRESENT
 
 
 INTEGER, INTENT(IN OUT)                  :: gplst(1)
 REAL, INTENT(IN)                         :: x(3,1)
 REAL, INTENT(IN)                         :: u(2,1)
 INTEGER, INTENT(IN OUT)                  :: pen
 INTEGER, INTENT(IN OUT)                  :: deform
 INTEGER, INTENT(OUT)                     :: nmax
 INTEGER, INTENT(OUT)                     :: maxsf
 INTEGER, INTENT(OUT)                     :: iz(14,1)
 INTEGER, INTENT(IN OUT)                  :: ib
 INTEGER, INTENT(IN)                      :: pedge
 INTEGER, INTENT(IN)                      :: iopcor
 LOGICAL :: shrink,hidden
 INTEGER :: etyp,g,NAME(2),gp,elid,  &
     elset, m1(16),ldx(9),FILE,solid,temp(27), offset
 
 COMMON /BLANK / ngp,skp11(11),elset,skp22(7),merr,idum(3),nscr1, nscr2,nscr3
 COMMON /system/ ibuf,iout
 COMMON /pltscr/ nnn,g(3)
 COMMON /hdrec / nofsur,ns,elid,lid,npers,p(3,13)
 
!     DIMENSIONS      TEMP, IZ, AND P ARE TEMP(2*N+1), IZ(N+1,1), AND
!                     P(3,N) WHERE N=LETSZ2=MAX OF LETSZ(2,I), I=1,9
 
 DIMENSION       let1(5),let2(4,4),let3(5,5),let4(5,6),let5(9,6),  &
     let6(13,6),let7(5),let8(7),let9(9),let(229), letsz(3,9)
 EQUIVALENCE     (let(  1),let1(  1)), (let(  6),let2(1,1)),  &
     (let( 22),let3(1,1)), (let( 47),let4(1,1)),  &
     (let( 77),let5(1,1)), (let(131),let6(1,1)),  &
     (let(209),let7(  1)), (let(214),let8(  1)), (let(221),let9(  1))
 
 DATA    NAME  / 4HHDSU, 4HRF   / ,  &
     nm1,m1/ 16,4H(33X, 4H,13H, 4HELEM, 4HENT , 4HTYPE, 4H a5,,  &
     4H4HWI, 4HTHI8, 4H,24H, 4H gri, 4HDS s, 4HKIPP, 4HED i,  &
     4HN li, 4HNEL., 4H)    /
 
!     SPECIAL ELEMENT CONNECTION PATTERNS
 
 DATA    ldx   / 2HD1,2HD2,2HD3,2HD4,2HD5,2HD6,2HD7,2HD8,2HD9  /
 DATA    ktet  / 2HTE /,   kweg  / 2HWG /,    khx1  / 2HH1 /,  &
     khx2  / 2HH2 /,   kix1  / 2HXL /,    kix2  / 2HXQ /,  &
     kix3  / 2HXC /,   kae   / 2HAE /,    ktrim6/ 2HT6 /,  &
     ktrplt/ 2HP6 /,   ktrshl/ 2HSL /,    kis2d8/ 2HD8 /,  &
     kfhex1/ 2HFA /,   kfhex2/ 2HFB /,    kfteta/ 2HFT /,  &
     kfwedg/ 2HFW /,   kbar  / 2HBR /,    kt3   / 2HT3 /, kq4   / 2HQ4 /
!    7        KELBOW/ 2HEB /
 
!     1   -   LINE,TRIANGLE,QUAD    5   -   IHEXA2
!     2   -   TETRA                 6   -   IHEXA3
!     3   -   WEDGE                 7   -   AERO
!     4   -   HEXA                  8   -   TRIM6 AND TRPLT1 AND TRSHL
 
 DATA    letsz2/ 13 /
 DATA    letsz / 1,      5,     1,  &
     4,      4,     6, 5,      5,    22,  &
     6,      5,    47, 6,      9,    77,  &
     6,     13,   131, 1,      5,   209,  &
     1,      7,   214, 1,      9,   221/
!         NELSRF,   NPTS,    IS
 DATA    let1  / 1,  2,  3,  4,  5/
 DATA    let2  / 1,  2,  3,  1,  &
     1,  2,  4,  1, 2,  3,  4,  2,  &
     1,  3,  4,  1/
 DATA    let3  / 1,  2,  3,  1,  0,  &
     4,  5,  6,  4,  0, 1,  3,  6,  4,  1,  &
     1,  2,  5,  4,  1, 2,  3,  6,  5,  2/
 DATA    let4  / 1,  2,  3,  4,  1,  &
     5,  6,  7,  8,  5, 3,  4,  8,  7,  3,  &
     1,  2,  6,  5,  1, 2,  3,  7,  6,  2,  &
     1,  4,  8,  5,  1/
 DATA    let5  / 1,  2,  3,  4,  5,  6,  7,  8,  1,  &
     13, 14, 15, 16, 17, 18, 19, 20, 13, 3, 10, 15, 16, 17, 11,  5,  4,  3,  &
     5, 11, 17, 18, 19, 12,  7,  6,  5, 7, 12, 19, 20, 13,  9,  1,  8,  7,  &
     1,  2,  3, 10, 15, 14, 13,  9,  1/
 DATA    let6  / 1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12,  1,  &
     21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 21,  &
     4,  5,  6,  7, 15, 19, 27, 26, 25, 24, 18, 14,  4,  &
     7,  8,  9, 10, 16, 20, 30, 29, 28, 27, 19, 15,  7,  &
     10, 11, 12,  1, 13, 17, 21, 32, 31, 30, 20, 16, 10,  &
     1,  2,  3,  9, 14, 18, 24, 23, 22, 21, 17, 13,  1/
 DATA    let7  / 1,  2,  3,  4,  1/
 DATA    let8  / 1,  2,  3,  4,  5,  6,  1/
 DATA    let9  / 1,  5,  2,  6,  3,  7,  4,  8,  1/
 
!     PEDGE FLAG = 2 OR 200    - HIDDEN LINE PLOT
!                = 10 THRU 100 - SHRINK PLOT.
!                = 100         - FILL, NOT USED HERE
!                = .GT.    200 - SHRINK AND HIDDEN LINE PLOT
!     E.G. PEDGE = 270 INDICATES HIDDEN LINE PLOT WITH EACH ELEMENT
!                  SHRUNK TO 70 PERCENT OF FULL SIZE.
 
 ipedge = MOD(pedge,200)
 nwds   = 0
 nmax   = 0
 ls     = 0
 lsmax  = iopcor/14
 IF (pedge > 200) lsmax = 0
 nofsur = 0
 shk    = 1.0
 shrink = .false.
 IF (pedge < 10) GO TO 10
 shrink = .true.
 shk = 1. - FLOAT(ipedge)/100.
 CALL line (0.,0.,0.,0.,1.,-1)
 10 hidden = .false.
 IF (pedge /= 2 .AND. pedge < 200) GO TO 20
 hidden = .true.
 CALL gopen (nscr2,gplst(ib),1)
 nwds = 3*letsz2 + 5
 
 20 CALL READ (*310,*190,elset,etyp,1,0,i)
 CALL fread (elset,i,1,0)
 ngpel  = IABS(i)
 ngpelx = ngpel
 solid  = 0
 
 offset = 0
 IF (etyp == kbar) offset = 6
 IF (etyp == kt3 .OR. etyp == kq4) offset = 1
 itype = 1
 IF (etyp == ktet    ) itype = 2
 IF (etyp == kfteta  ) itype = 2
 IF (etyp == kweg    ) itype = 3
 IF (etyp == kfwedg  ) itype = 3
 IF (etyp == khx1   .OR. etyp == khx2 .OR. etyp == kix1 .OR.  &
     etyp == kfhex1 .OR. etyp == kfhex2) itype = 4
 IF (etyp == kix2    ) itype = 5
 IF (etyp == kis2d8  ) itype = 9
 IF (etyp == kix3    ) itype = 6
 IF (etyp == kae     ) itype = 7
 IF (etyp == ktrim6 .OR. etyp == ktrplt .OR. etyp == ktrshl) itype = 8
 
 IF (itype /= 1) GO TO 40
 
!     SIMPLE ELEMENT
 
 IF (ngpel > 2 .AND. i > 0) ngpelx = ngpel + 1
 IF (ngpel > 4) GO TO 130
 npts = ngpelx
 GO TO 50
 
!     COMPLEX ELEMENT
 
 40 CONTINUE
 IF (itype >= 2 .AND. itype <= 6) solid = 1
 npts = letsz(2,itype)
 50 IF (npts-1 > nmax) nmax = npts - 1
 
!     READ THE ELEMENT DATA
 
 60 CALL fread (elset,elid,1,0)
 IF (elid <= 0) GO TO 20
 CALL fread (elset,lid,1,0)
 CALL fread (elset,g,ngpel,0)
 IF (offset /= 0) CALL fread (elset,0,-offset,0)
 IF (ngpel /= ngpelx) g(ngpelx) = g(1)
 IF (hidden .AND. .NOT.shrink) GO TO 80
 xc = 0.
 yc = 0.
 zc = 0.
 DO  i = 1,ngpel
   gp = g(i)
   gp = IABS(gplst(gp))
   xc = xc + x(2,gp)
   yc = yc + x(3,gp)
   zc = zc + x(1,gp)
 END DO
 xc = xc/ngpel
 yc = yc/ngpel
 zc = zc/ngpel
 
 80 nelsrf = letsz(1,itype)
 is = letsz(3,itype)
 
 DO  ns = 1,nelsrf
   nn = 0
   mm = (ns-1)*npts + is - 1
   npers = npts
   DO  i = 1,npts
     m = mm + i
     n = let(m)
     IF (n /= 0) GO TO 85
     82 npers = npers - 1
     CYCLE
     85 gp = g(n)
     IF (gp == 0) GO TO 82
     nn = nn + 1
     gp = IABS(gplst(gp))
     p(3,nn) = x(1,gp)
     IF (deform /= 0) GO TO 90
     p(1,nn) = x(2,gp)
     p(2,nn) = x(3,gp)
     GO TO 100
     90 p(1,nn) = u(1,gp)
     p(2,nn) = u(2,gp)
     100 CONTINUE
     IF (.NOT.shrink) CYCLE
     IF (     hidden) GO TO 105
     IF (nn == 1) CYCLE
     x1 = p(1,nn-1) - (p(1,nn-1)-xc)*shk
     y1 = p(2,nn-1) - (p(2,nn-1)-yc)*shk
     x2 = p(1,nn  ) - (p(1,nn  )-xc)*shk
     y2 = p(2,nn  ) - (p(2,nn  )-yc)*shk
     ipen = pen
     IF (ipedge == 100 .AND. pen > 31 .AND. i == npers) pen = 0
     IF (shrink) CALL line (x1,y1,x2,y2,pen,0)
     IF (pen == 0) pen = ipen
     CYCLE
     105 p(3,nn) = x(1,gp) - (x(1,gp)-zc)*shk
     p(1,nn) = x(2,gp) - (x(2,gp)-xc)*shk
     p(2,nn) = x(3,gp) - (x(3,gp)-yc)*shk
     IF (deform == 0) CYCLE
     p(1,nn) = u(1,gp) - (x(2,gp)-xc)*shk
     p(2,nn) = u(2,gp) - (x(3,gp)-yc)*shk
   END DO
   IF (shrink .AND. .NOT.hidden) CYCLE
   CALL WRITE (nscr2,nofsur,nwds,0)
   nofsur = nofsur + 1
   IF (solid == 0 .OR. .NOT.hidden) CYCLE
   
!     SAVE SOLID SURFACE DATA IN IZ SPACE FOR SECOND PROCESSING, HIDDEN
!     PLOT ONLY. SAVE AS MANY AS OPEN CORE CAN HOLD
   
   IF (ls >= lsmax) CYCLE
   ls = ls + 1
   nps1 = npers - 1
   DO  i = 1,nps1
     m  = mm + i
     n  = let(m)
     gp = g(n)
     temp(i     ) = gp
     temp(i+nps1) = gp
   END DO
   m  = 1
   MIN= temp(1)
   DO  i = 2,nps1
     IF (temp(i) >= MIN) CYCLE
     m  = i
     MIN= temp(i)
   END DO
   IF (m == 1) m = m + nps1
   n = + 1
   IF (temp(m-1) < temp(m+1)) n = -1
   IF (n == -1  .AND. m < nps1) m = m + nps1
   k = nps1 + 2
   DO  i = 3,k
     iz(i,ls) = temp(m)
     m = m + n
   END DO
   iz(1,ls) = nofsur
   iz(2,ls) = nps1
   
 END DO
 GO TO 60
 
!     CHECK FOR PDUM ELEMENTS BEFORE  EJECTING
 
 130 DO  i = 1,9
   IF (etyp == ldx(i)) GO TO 160
 END DO
 
!     ILLEGAL ELEMENT, NO CORE FOR 1 ELEMENT
 
 140 g(1) = 2
 g(2) = etyp
 g(3) = ngpel
 CALL wrtprt (merr,g,m1,nm1)
 
!     READ TO THE END OF THIS ELEMENT
 
 150 CALL READ (*180,*20,elset,elid,1,0,m)
 IF (elid <= 0) GO TO 20
 j = 1 + ngpel + offset
 CALL fread (elset,0,-j,0)
 GO TO 150
 160 WRITE  (iout,170) i
 170 FORMAT ('0*** MISSING PDUM',i1,' SUBROUTINE/HDSURF')
 GO TO 140
 180 CALL mesage (-8,elset,NAME)
 
 190 CONTINUE
 maxsf = nofsur
 CALL bckrec (elset)
 IF (shrink) CALL line (0.,0.,0.,0.,1.,+1)
 IF (.NOT.hidden) GO TO 300
 CALL WRITE (nscr2,0,0,1)
 IF (ls < 60) GO TO 280
 
!     REPROCESS NSCR2 TO REMOVE DUPLICATE SURFACES (INTERIOR-INTERFACES)
!     AND SAVE REDUCED DATA IN NSCR1.
!     INTERCHANGE NSCR1 AND NSCR2 INDICES
 
 j = (letsz2+1)*ls
 CALL sort2k (0,0,letsz2+1,3,iz,j)
 m = 0
 nps1 = 0
 loop240:  DO  i = 1,ls
   nps2 = iz(2,i) + 2
   IF (nps2 == nps1) GO TO 200
   nps1 = nps2
   CYCLE loop240
   200 im1  = i - 1
   DO  j = 3,nps1
     IF (iz(j,i) /= iz(j,im1)) CYCLE loop240
   END DO
   IF (m == 0) GO TO 220
   IF (iz(m,1) == iz(1,im1)) GO TO 230
   220 m = m + 1
   iz(m,1) = iz(1,im1)
   230 m = m + 1
   iz(m,1) = iz(1,i)
 END DO loop240
 
 IF (m < 20) GO TO 280
 CALL sort (0,0,1,1,iz,m)
 iz(m+1,1) = 999999999
 FILE = nscr1
 CALL gopen (nscr1,gplst(ib+ibuf),1)
 FILE = nscr2
 CALL CLOSE (nscr2,1)
 CALL gopen (nscr2,gplst(ib),0)
 n = 1
 DO  i = 1,maxsf
   CALL READ (*320,*330,nscr2,nofsur,nwds,0,j)
   IF (i-iz(n,1) < 0) THEN
     GO TO   250
   ELSE
     GO TO   260
   END IF
   250 CALL WRITE (nscr1,nofsur,nwds,0)
   CYCLE
   260 n = n + 1
 END DO
 
 CALL CLOSE (nscr2,1)
 j = nscr2
 nscr2 = nscr1
 nscr1 = j
 maxsf = maxsf - m
 CALL WRITE (nscr2,0,0,1)
 280 CALL CLOSE (nscr2,1)
 300 RETURN
 
 310 j = -1
 FILE = elset
 GO TO 340
 320 j = -2
 GO TO 340
 330 j = -3
 340 CALL mesage (j,FILE,NAME)
 GO TO 190
END SUBROUTINE hdsurf
